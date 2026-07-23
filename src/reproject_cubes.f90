! reproject_cubes -- standalone pre-rm-synthesis tool (geometry-matching
! project, planning/MULTI_BAND_TOMOGRAPHY_PLAN.md): reprojects a set of
! FITS cubes -- two or more, not tied to the multi-band-tomography "band"
! concept specifically -- onto a single reference cube's grid using
! Starlink AST for WCS handling and resampling, so the existing
! rm_synthesis exact-match ingestion (T1) can consume genuinely misaligned
! bands unchanged.
!
! Current stage: N-input footprint-mode computation (intersection/union/
! reference), built on top of two previously-verified building blocks:
! extract_sky_mapping (automatic sky-axis detection, correct for any axis
! order/adjacency -- see its own comment for why the earlier consecutive-
! axis-only search was wrong) and astMapBox-based footprint bounds
! (verified against a genuine partial-overlap case, not just full
! overlap). Cross-file alignment composes each file's own pixel->sky
! Mapping with an astConvert between the two SkyFrame *objects* (not
! whole FrameSets -- astConvert's domain search does not recurse into a
! CmpFrame's internal components, confirmed against the actual SUN/211
! manual, so it cannot align two compound "STOKES-SKY-SPECTRUM" current
! frames directly).
!
! Footprint policy: zero overlap between any input and the running output
! grid is always a hard failure, regardless of mode (loudly refuse before
! compute, matching the rest of this project's philosophy) -- this is
! different from *partial* overlap, which is the normal, expected case
! for combining bands with different sky coverage and is not rejected.
!
! Also demonstrates the actual regridding step: resample_first_plane
! reads one input's first plane and resamples it onto the final output
! grid via astResampleR. Caught a real bug on the first attempt: the
! documented 20-argument signature (this, ndim_in, lbnd_in, ubnd_in, in,
! in_var, interp, finterp, params, flags, tol, maxpix, badval, ndim_out,
! lbnd_out, ubnd_out, lbnd, ubnd, out, out_var) was missing "params"
! entirely, silently shifting every later argument by one position (a
! REAL array landing where an INTEGER "flags" scalar was expected, etc)
! -- segfaulted rather than erroring cleanly, so this needed bisecting
! with diagnostic prints rather than reading a status code. Verified
! after the fix: resampling the axis-swapped fixture onto the reference
! grid reproduces the exact known ground-truth value at a known pixel
! (read independently via Python beforehand), and the union-mode test's
! NaN count for uncovered output pixels matches the expected uncovered
! area exactly (16x32 = 512).
!
! Usage: reproject_cubes <intersection|union|reference> <reference_file> <input_file> [input_file ...]
program reproject_cubes
   implicit none
   ! AST_PAR (the vendor Fortran constants file, /usr/include/AST_PAR) is
   ! fixed-form Fortran 77 (`*`-column comments) and cannot be `include`d
   ! into a free-form .f90 file directly (gfortran misparses its comments
   ! as code). Only the handful of symbols actually used are declared
   ! directly instead, matching AST_PAR's own declared types exactly
   ! (checked against /usr/include/AST_PAR).
   external :: ast_null
   integer, parameter :: ast__null = 0
   integer, parameter :: ast__szchr = 200
   integer, parameter :: ast__base = 0
   integer, parameter :: ast__current = -1
   integer, external :: ast_fitschan, ast_read, ast_geti
   integer, external :: ast_getmapping, ast_simplify, ast_getframe
   integer, external :: ast_pickaxes, ast_cmpmap, ast_convert
   integer, external :: ast_resampler
   integer, parameter :: ast__linear = 5
   logical, external :: ast_isaframeset, ast_isaskyframe
   character(len=ast__szchr), external :: ast_getc

   integer, parameter :: max_axes = 10
   integer, parameter :: max_inputs = 50

   character(len=16) :: mode
   character(len=512) :: reffile
   character(len=512) :: infiles(max_inputs)
   integer :: n_inputs, iarg, i

   integer :: wcs_ref, skymap_ref, skyframe_ref
   integer :: naxes_ref(max_axes), pixaxes_ref(2)
   integer :: wcs_in, skymap_in, skyframe_in
   integer :: naxes_in(max_axes), pixaxes_in(2)
   integer :: map_in2ref
   double precision :: lbnd_out(2), ubnd_out(2)
   double precision :: this_lbnd(2), this_ubnd(2)
   integer :: status

   if (command_argument_count() < 3) then
      write(*,*) 'Usage: reproject_cubes <intersection|union|reference>',&
      &' <reference_file> <input_file> [input_file ...]'
      stop 1
   endif
   call get_command_argument(1, mode)
   call get_command_argument(2, reffile)
   n_inputs = command_argument_count() - 2
   if (n_inputs.gt.max_inputs) then
      write(*,*) 'ERROR: too many input files (max ', max_inputs, ')'
      stop 1
   endif
   do i = 1, n_inputs
      call get_command_argument(2+i, infiles(i))
   enddo

   if (trim(mode).ne.'intersection' .and. trim(mode).ne.'union' .and.&
   &trim(mode).ne.'reference') then
      write(*,*) 'ERROR: mode must be intersection, union, or reference'
      stop 1
   endif

   status = 0
   call ast_begin(status)

   call load_wcs(reffile, wcs_ref, naxes_ref, status)
   call extract_sky_mapping(wcs_ref, skymap_ref, skyframe_ref, pixaxes_ref, status)
   if (status.ne.0) then
      write(*,*) 'ERROR: failed to load the reference file''s WCS'
      stop 1
   endif

   ! Output grid starts as the reference's own full extent; intersection
   ! shrinks it, union grows it, reference mode leaves it untouched (the
   ! loop below is skipped entirely for reference mode).
   lbnd_out(1) = 1.0d0
   lbnd_out(2) = 1.0d0
   ubnd_out(1) = real(naxes_ref(pixaxes_ref(1)), kind=8)
   ubnd_out(2) = real(naxes_ref(pixaxes_ref(2)), kind=8)
   write(*,'(A,A,A,F0.0,A,F0.0,A,F0.0,A,F0.0,A)') 'Reference (', trim(reffile),&
   &') own extent: [', lbnd_out(1), ',', ubnd_out(1), '] x [',&
   &lbnd_out(2), ',', ubnd_out(2), ']'

   if (trim(mode).ne.'reference') then
      do i = 1, n_inputs
         call load_wcs(infiles(i), wcs_in, naxes_in, status)
         call extract_sky_mapping(wcs_in, skymap_in, skyframe_in, pixaxes_in, status)
         if (status.ne.0) then
            write(*,*) 'ERROR: failed to load input file: ', trim(infiles(i))
            stop 1
         endif

         call compose_pix2pix(skymap_in, skyframe_in, skymap_ref, skyframe_ref,&
         &map_in2ref, status)
         if (status.ne.0) then
            write(*,*) 'ERROR: failed to align input file to the reference: ',&
            &trim(infiles(i))
            stop 1
         endif

         call footprint_bounds(map_in2ref, naxes_in, pixaxes_in,&
         &this_lbnd, this_ubnd, status)
         write(*,'(A,A,A,F0.2,A,F0.2,A,F0.2,A,F0.2,A)') '  ', trim(infiles(i)),&
         &' footprint in reference space: [', this_lbnd(1), ',', this_ubnd(1),&
         &'] x [', this_lbnd(2), ',', this_ubnd(2), ']'

         ! Zero overlap with the running output grid is always a hard
         ! failure, regardless of mode -- a band sharing no sky at all
         ! with the rest is not legitimate partial coverage (that implies
         ! at least some shared sky), it is almost certainly the wrong
         ! file. Checked here, per-input, before folding it into the
         ! running bound, so the diagnostic names the specific offending
         ! file rather than a generic "empty result" at the end.
         if (this_ubnd(1).lt.lbnd_out(1) .or. this_lbnd(1).gt.ubnd_out(1) .or.&
         &this_ubnd(2).lt.lbnd_out(2) .or. this_lbnd(2).gt.ubnd_out(2)) then
            write(*,*) 'ERROR: zero sky overlap between the reference and: ',&
            &trim(infiles(i))
            write(*,*) 'Quitting now...'
            stop 1
         endif

         if (trim(mode).eq.'intersection') then
            lbnd_out = max(lbnd_out, this_lbnd)
            ubnd_out = min(ubnd_out, this_ubnd)
         else ! union
            lbnd_out = min(lbnd_out, this_lbnd)
            ubnd_out = max(ubnd_out, this_ubnd)
         endif

         call ast_annul(map_in2ref, status)
         call ast_annul(skymap_in, status)
         call ast_annul(skyframe_in, status)
         call ast_annul(wcs_in, status)
      enddo
   endif

   ! Round to the integer pixel grid: intersection wants the largest
   ! integer range fully CONTAINED within the real-valued bound (ceiling
   ! the lower edge, floor the upper edge); union wants the smallest
   ! integer range that fully CONTAINS it (floor/ceiling the other way).
   ! Reference mode's bound is already exactly integer (a file's own
   ! NAXIS), so this is a no-op for it either way.
   if (trim(mode).eq.'intersection') then
      lbnd_out = ceiling(lbnd_out)
      ubnd_out = floor(ubnd_out)
   else
      lbnd_out = floor(lbnd_out)
      ubnd_out = ceiling(ubnd_out)
   endif

   if (lbnd_out(1).gt.ubnd_out(1) .or. lbnd_out(2).gt.ubnd_out(2)) then
      write(*,*) 'ERROR: computed output grid is empty (', trim(mode), ' mode)'
      stop 1
   endif

   ! --- Resample the first input's first plane onto the final output
   ! grid, via astResampleR ---
   ! Proof-of-concept for the actual regridding step: re-derives the
   ! first input's own pixel->reference Mapping (cheap relative to the
   ! resample itself; the per-channel production loop will keep Mappings
   ! alive across a single pass instead of recomputing, once this
   ! correctness case is established). Only the first non-sky-axis plane
   ! (e.g. channel 1) is handled here -- looping over every plane is the
   ! next increment, not this one.
   call load_wcs(infiles(1), wcs_in, naxes_in, status)
   call extract_sky_mapping(wcs_in, skymap_in, skyframe_in, pixaxes_in, status)
   call compose_pix2pix(skymap_in, skyframe_in, skymap_ref, skyframe_ref,&
   &map_in2ref, status)
   if (status.ne.0) then
      write(*,*) 'ERROR: failed to re-derive the first input''s Mapping',&
      &' for resampling'
      stop 1
   endif
   call resample_first_plane(infiles(1), map_in2ref, naxes_in, pixaxes_in,&
   &lbnd_out, ubnd_out, status)
   call ast_annul(map_in2ref, status)
   call ast_annul(skymap_in, status)
   call ast_annul(skyframe_in, status)
   call ast_annul(wcs_in, status)

   write(*,'(A,A,A,F0.0,A,F0.0,A,F0.0,A,F0.0,A)') 'Final output grid (',&
   &trim(mode), ' mode): [', lbnd_out(1), ',', ubnd_out(1), '] x [',&
   &lbnd_out(2), ',', ubnd_out(2), ']'

   call ast_annul(skymap_ref, status)
   call ast_annul(skyframe_ref, status)
   call ast_annul(wcs_ref, status)
   call ast_end(status)

   if (status.ne.0) then
      write(*,*) 'ERROR: AST reported an error, final status=', status
      stop 1
   endif

   write(*,*) 'OK: footprint-mode output grid computed successfully.'

contains

   subroutine resample_first_plane(filename, map_in2ref, naxes_in,&
   &pixaxes_in, lbnd_out_d, ubnd_out_d, status)
      !! Read filename's first plane (every non-sky axis fixed at index 1)
      !! and resample it onto the output grid [lbnd_out_d,ubnd_out_d]
      !! (reference pixel space) via astResampleR, using map_in2ref (the
      !! forward input->reference pixel Mapping -- astResampleR uses its
      !! INVERSE internally to look up each output pixel's input value,
      !! per its own documented convention). Bad/uncovered output pixels
      !! are flagged with IEEE NaN, matching this codebase's existing
      !! NaN-based masking convention throughout rm_synthesis.
      !!
      !! Reading fpixels(k)=lpixels(k)=1 for every axis except the 2 sky
      !! ones (set to that axis's own full 1..NAXIS extent) naturally
      !! produces a flat array with the sky axes as its only two
      !! non-degenerate dimensions, in pixaxes_in's own order (first
      !! fastest) -- exactly the "Fortran array indexing" astResampleR
      !! expects -- regardless of which raw file-axis numbers they are,
      !! since a degenerate (extent-1) axis contributes no stride.
      use, intrinsic :: ieee_arithmetic
      character(len=*), intent(in) :: filename
      integer, intent(in) :: map_in2ref
      integer, intent(in) :: naxes_in(:), pixaxes_in(2)
      double precision, intent(in) :: lbnd_out_d(2), ubnd_out_d(2)
      integer, intent(inout) :: status

      integer :: unit, blocksize, fitsstat, group, naxis, k
      integer :: fpixels(max_axes), lpixels(max_axes), incs(max_axes)
      integer :: nx_in, ny_in, nx_out, ny_out
      real, allocatable :: data_in(:,:), data_out(:,:)
      logical :: anyflg
      integer :: lbnd_in(2), ubnd_in(2), lbnd_o(2), ubnd_o(2)
      integer :: nbad
      real :: badval
      double precision :: params_dummy(1)

      if (status.ne.0) return

      nx_in = naxes_in(pixaxes_in(1))
      ny_in = naxes_in(pixaxes_in(2))
      allocate(data_in(nx_in, ny_in))

      naxis = 0
      do k = 1, max_axes
         if (naxes_in(k).gt.0) naxis = k
      enddo
      fpixels(1:naxis) = 1
      lpixels(1:naxis) = 1
      incs(1:naxis) = 1
      lpixels(pixaxes_in(1)) = nx_in
      lpixels(pixaxes_in(2)) = ny_in

      fitsstat = 0
      blocksize = 1
      unit = 31
      group = 1
      badval = ieee_value(badval, ieee_quiet_nan)
      call FTOPEN(unit, trim(filename), 0, blocksize, fitsstat)
      call FTGSVE(unit, group, naxis, naxes_in(1:naxis),&
      &fpixels(1:naxis), lpixels(1:naxis), incs(1:naxis),&
      &badval, data_in, anyflg, fitsstat)
      call FTCLOS(unit, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to read data plane from ', trim(filename)
         call printerror(fitsstat)
         status = -1
         return
      endif

      nx_out = nint(ubnd_out_d(1) - lbnd_out_d(1)) + 1
      ny_out = nint(ubnd_out_d(2) - lbnd_out_d(2)) + 1
      allocate(data_out(nx_out, ny_out))

      lbnd_in(1) = 1
      lbnd_in(2) = 1
      ubnd_in(1) = nx_in
      ubnd_in(2) = ny_in
      lbnd_o(1) = nint(lbnd_out_d(1))
      lbnd_o(2) = nint(lbnd_out_d(2))
      ubnd_o(1) = nint(ubnd_out_d(1))
      ubnd_o(2) = nint(ubnd_out_d(2))
      params_dummy(1) = 0.0d0

      ! Full 20-argument signature (this, ndim_in, lbnd_in, ubnd_in, in,
      ! in_var, interp, finterp, params, flags, tol, maxpix, badval,
      ! ndim_out, lbnd_out, ubnd_out, lbnd, ubnd, out, out_var) -- an
      ! earlier attempt omitted "params" entirely, silently shifting every
      ! later argument by one position (a REAL array landing where an
      ! INTEGER "flags" scalar was expected, etc), which segfaulted.
      ! in_var/out_var are unused here (no variance requested) but AST's
      ! Fortran binding needs a validly-typed/sized array regardless of a
      ! true C NULL, so data_in/data_out are reused as harmless
      ! placeholders rather than risk an invalid "null" convention guess.
      nbad = ast_resampler(map_in2ref, 2, lbnd_in, ubnd_in, data_in, data_in,&
      &ast__linear, ast_null, params_dummy, 0, 0.0d0, 100, badval,&
      &2, lbnd_o, ubnd_o, lbnd_o, ubnd_o, data_out, data_out, status)

      write(*,'(A,A,A,I0,A,I0,A)') 'Resampled ', trim(filename), ': ',&
      &nx_out, ' x ', ny_out, ' output plane'
      write(*,'(A,I0,A)') '  bad (uncovered) output pixels: ', nbad

      ! Diagnostic: print the resampled value at reference pixel (5,9),
      ! for comparison against a known ground-truth value read directly
      ! from the original data with a Python script beforehand.
      if (5.ge.lbnd_o(1) .and. 5.le.ubnd_o(1) .and.&
      &9.ge.lbnd_o(2) .and. 9.le.ubnd_o(2)) then
         write(*,'(A,F0.9)') '  value at reference pixel (5,9): ',&
         &data_out(5-lbnd_o(1)+1, 9-lbnd_o(2)+1)
      endif

      deallocate(data_in)
      deallocate(data_out)
   end subroutine resample_first_plane

   subroutine compose_pix2pix(skymap_from, skyframe_from, skymap_to,&
   &skyframe_to, map_out, status)
      !! Compose pixel_from -> sky_from -> sky_to -> pixel_to, aligning
      !! sky_from/sky_to via a Frame-to-Frame astConvert between the two
      !! SkyFrame objects (handles any axis-order/equinox/system
      !! difference between them; see extract_sky_mapping's own comment
      !! for why this cannot be skipped -- composing the two pixel->sky
      !! Mappings directly, without this alignment step, silently
      !! produces wrong results whenever the two files' SkyFrames present
      !! their axes in a different order).
      integer, intent(in) :: skymap_from, skyframe_from
      integer, intent(in) :: skymap_to, skyframe_to
      integer, intent(out) :: map_out
      integer, intent(inout) :: status

      integer :: sky2sky

      map_out = ast__null
      if (status.ne.0) return

      sky2sky = ast_convert(skyframe_from, skyframe_to, ' ', status)
      if (status.ne.0 .or. sky2sky.eq.ast__null) then
         write(*,*) 'ERROR: failed to align two SkyFrames, status=', status
         status = -1
         return
      endif

      call ast_invert(skymap_to, status)
      map_out = ast_cmpmap(skymap_from, sky2sky, .true., ' ', status)
      map_out = ast_cmpmap(map_out, skymap_to, .true., ' ', status)
      call ast_invert(skymap_to, status)
      if (status.ne.0 .or. map_out.eq.ast__null) then
         write(*,*) 'ERROR: failed to compose the pixel->pixel Mapping,',&
         &' status=', status
         status = -1
      endif

      call ast_annul(sky2sky, status)
   end subroutine compose_pix2pix

   subroutine footprint_bounds(map_from_to, naxes_from, pixaxes_from,&
   &lbnd, ubnd, status)
      !! "from" file's full pixel extent (its own NAXIS on the 2 axes its
      !! sky Mapping depends on), expressed in "to" pixel space, via
      !! astMapBox (the true enclosing bound of each output coordinate --
      !! not a naive 4-corner check, which can underestimate the true
      !! extent for a non-axis-aligned Mapping). map_from_to must be the
      !! forward "from"->"to" pixel Mapping (as returned by
      !! compose_pix2pix with "from" as its first, "to" as its second
      !! argument pair).
      integer, intent(in) :: map_from_to
      integer, intent(in) :: naxes_from(:), pixaxes_from(2)
      double precision, intent(out) :: lbnd(2), ubnd(2)
      integer, intent(inout) :: status

      double precision :: lbnd_in(2), ubnd_in(2), xl(2), xu(2)

      lbnd = 0.0d0
      ubnd = 0.0d0
      if (status.ne.0) return

      lbnd_in(1) = 1.0d0
      lbnd_in(2) = 1.0d0
      ubnd_in(1) = real(naxes_from(pixaxes_from(1)), kind=8)
      ubnd_in(2) = real(naxes_from(pixaxes_from(2)), kind=8)

      call ast_mapbox(map_from_to, lbnd_in, ubnd_in, .true., 1,&
      &lbnd(1), ubnd(1), xl, xu, status)
      call ast_mapbox(map_from_to, lbnd_in, ubnd_in, .true., 2,&
      &lbnd(2), ubnd(2), xl, xu, status)
   end subroutine footprint_bounds

   subroutine extract_sky_mapping(wcs, skymap, skyframe, pixel_axes, status)
      !! Extract the pixel-grid -> sky (RA/Dec) Mapping from a WCS
      !! FrameSet, with the sky axes' positions in the (possibly compound
      !! Stokes+Sky+Spectrum) current frame detected automatically -- no
      !! assumption about which positions they occupy. Every axis PAIR
      !! (not just consecutive ones -- an earlier version assumed a
      !! SkyFrame's 2 axes must be consecutive within a CmpFrame, since
      !! decompose splits a CmpFrame into two axis-contiguous components;
      !! that assumption turned out to be wrong, confirmed by direct
      !! test: a file with RA/Dec on non-adjacent pixel axes 2 and 4 still
      !! has a genuine SkyFrame recoverable via ast_pickaxes(2,4), so
      !! whatever internal structure connects them is not simply "the two
      !! components of one decompose split") is probed with ast_pickaxes
      !! + ast_isaskyframe -- a genuine AST class check, not a guess -- to
      !! find which pair it is. Once known, astMapSplit selects a
      !! Mapping's INPUT axes, but the sky axes are on the OUTPUT side
      !! (current frame) of the pixel->compound Mapping -- so invert first
      !! (making sky axes selectable as inputs), simplify (helps AST
      !! recognise separability), split, then invert the result back to a
      !! forward pixel-subset -> sky Mapping. Also returns the isolated
      !! SkyFrame object itself (not just the Mapping) -- a SkyFrame's own
      !! axis order is NOT a fixed RA-then-Dec convention (it reflects
      !! whichever axis the header declared as longitude vs latitude
      !! first, confirmed by direct comparison: a file with CTYPE1=DEC
      !! presents (Dec,RA), not (RA,Dec)), so the caller needs the actual
      !! Frame object to align two files' sky axes correctly via
      !! astConvert rather than assuming a shared order.
      integer, intent(in) :: wcs
      integer, intent(out) :: skymap, skyframe
      integer, intent(out) :: pixel_axes(2)
      integer, intent(inout) :: status

      integer :: curframe, nout, i, j
      integer :: probe_axes(2), probe_frame, probe_map
      integer :: sky_axes_in(2), out_axes(4)
      integer :: fullmap, simplemap
      logical :: found_sky

      skymap = ast__null
      skyframe = ast__null
      pixel_axes = 0
      if (status.ne.0) return

      nout = ast_geti(wcs, 'Nout', status)
      curframe = ast_getframe(wcs, ast__current, status)
      found_sky = .false.
      outer: do i = 1, nout - 1
         do j = i + 1, nout
            probe_axes(1) = i
            probe_axes(2) = j
            probe_frame = ast_pickaxes(curframe, 2, probe_axes, probe_map, status)
            if (ast_isaskyframe(probe_frame, status)) then
               sky_axes_in = probe_axes
               found_sky = .true.
               skyframe = probe_frame
               call ast_annul(probe_map, status)
               exit outer
            endif
            call ast_annul(probe_frame, status)
            call ast_annul(probe_map, status)
         enddo
      enddo outer
      call ast_annul(curframe, status)

      if (.not. found_sky) then
         write(*,*) 'ERROR: no axis-pair SkyFrame found in the WCS'
         status = -1
         return
      endif

      fullmap = ast_getmapping(wcs, ast__base, ast__current, status)
      call ast_invert(fullmap, status)
      simplemap = ast_simplify(fullmap, status)
      call ast_mapsplit(simplemap, 2, sky_axes_in, out_axes, skymap, status)
      if (status.ne.0 .or. skymap.eq.ast__null) then
         write(*,*) 'ERROR: ast_mapsplit failed to isolate the sky Mapping,',&
         &' status=', status
         status = -1
         call ast_annul(fullmap, status)
         call ast_annul(simplemap, status)
         return
      endif
      call ast_invert(skymap, status)
      pixel_axes = out_axes(1:2)

      call ast_annul(fullmap, status)
      call ast_annul(simplemap, status)
   end subroutine extract_sky_mapping

   subroutine load_wcs(filename, wcs, naxes, status)
      !! Read filename's FITS header via CFITSIO, load it into an AST
      !! FitsChan, and return the WCS FrameSet recovered from it plus this
      !! file's own per-axis pixel-grid extent (NAXISn), needed later to
      !! bound each axis's own footprint for astMapBox.
      character(len=*), intent(in) :: filename
      integer, intent(out) :: wcs
      integer, intent(out) :: naxes(:)
      integer, intent(inout) :: status

      integer :: unit, blocksize, fitsstat, nkeys, nmore, i, fitschan
      character(len=80) :: card

      wcs = ast__null
      naxes = 0
      if (status.ne.0) return

      fitsstat = 0
      blocksize = 1
      unit = 11
      call FTOPEN(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to open FITS file: ', trim(filename)
         call printerror(fitsstat)
         status = -1
         return
      endif

      call FTGHSP(unit, nkeys, nmore, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: FTGHSP failed for ', trim(filename)
         call printerror(fitsstat)
         call FTCLOS(unit, fitsstat)
         status = -1
         return
      endif

      call FTGISZ(unit, size(naxes), naxes, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: FTGISZ failed for ', trim(filename)
         call printerror(fitsstat)
         call FTCLOS(unit, fitsstat)
         status = -1
         return
      endif

      fitschan = ast_fitschan(ast_null, ast_null, ' ', status)
      if (status.ne.0) then
         call FTCLOS(unit, fitsstat)
         return
      endif

      do i = 1, nkeys
         fitsstat = 0
         call FTGREC(unit, i, card, fitsstat)
         if (fitsstat.ne.0) then
            write(*,*) 'ERROR: FTGREC failed at card ', i, ' for ', trim(filename)
            call printerror(fitsstat)
            call FTCLOS(unit, fitsstat)
            status = -1
            return
         endif
         call ast_putfits(fitschan, card, .false., status)
      enddo
      ! Rewind the channel's internal read pointer before ast_read
      ! consumes it (Card is a 1-based cursor into the card list;
      ! ast_putfits above leaves it sitting past the last card written).
      call ast_seti(fitschan, 'Card', 1, status)

      call FTCLOS(unit, fitsstat)

      wcs = ast_read(fitschan, status)
      if (status.ne.0 .or. wcs.eq.ast__null) then
         write(*,*) 'ERROR: ast_read failed to recover a WCS FrameSet for ',&
         &trim(filename)
         status = -1
         return
      endif
      if (.not. ast_isaframeset(wcs, status)) then
         write(*,*) 'ERROR: object read from FitsChan is not a FrameSet for ',&
         &trim(filename)
         status = -1
         return
      endif
      call ast_annul(fitschan, status)
   end subroutine load_wcs

end program reproject_cubes
