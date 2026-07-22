! reproject_cubes -- standalone pre-rm-synthesis tool (geometry-matching
! project, planning/MULTI_BAND_TOMOGRAPHY_PLAN.md): reprojects a set of
! FITS cubes -- two or more, not tied to the multi-band-tomography "band"
! concept specifically -- onto a single reference cube's grid using
! Starlink AST for WCS handling and resampling, so the existing
! rm_synthesis exact-match ingestion (T1) can consume genuinely misaligned
! bands unchanged.
!
! Current stage: cross-file pixel(A)->pixel(B) Mapping, composed from each
! file's own pixel->sky Mapping (extract_sky_mapping, proven in the
! previous commit: automatic axis detection via ast_isaskyframe + a clean
! astMapSplit extraction) rather than via astConvert. astConvert was tried
! first (the textbook cross-WCS-alignment primitive) but, per the actual
! SUN/211 manual (installed via libstarlink-ast-doc), its domain search
! only examines each FrameSet's own top-level registered frames (base +
! current) by their own Domain attribute -- it does not recurse into a
! CmpFrame's internal components, so a domainlist of 'SKY' can never match
! a compound "STOKES-SKY-SPECTRUM" current frame. Composing pixel_A->sky
! (this file's own Mapping) with sky->pixel_B (the other file's own
! Mapping, inverted) sidesteps that limitation entirely.
!
! Usage: reproject_cubes <fits_file_a> <fits_file_b>
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
   logical, external :: ast_isaframeset, ast_isaskyframe
   character(len=ast__szchr), external :: ast_getc

   character(len=512) :: infile_a, infile_b
   integer :: wcs_a, wcs_b, skymap_a, skymap_b, pix2pix
   integer :: skyframe_a, skyframe_b, sky2sky
   integer :: status
   double precision :: xin(1), yin(1), xout(1), yout(1)

   if (command_argument_count() < 2) then
      write(*,*) 'Usage: reproject_cubes <fits_file_a> <fits_file_b>'
      stop 1
   endif
   call get_command_argument(1, infile_a)
   call get_command_argument(2, infile_b)

   status = 0
   call ast_begin(status)

   call load_wcs(infile_a, wcs_a, status)
   call load_wcs(infile_b, wcs_b, status)
   if (status.ne.0) then
      write(*,*) 'ERROR: failed to load one or both WCS FrameSets'
      stop 1
   endif

   call extract_sky_mapping(wcs_a, skymap_a, skyframe_a, status)
   call extract_sky_mapping(wcs_b, skymap_b, skyframe_b, status)
   if (status.ne.0) then
      write(*,*) 'ERROR: failed to extract one or both pixel->sky Mappings'
      stop 1
   endif
   write(*,*) 'Both files'' pixel -> sky Mappings extracted successfully.'

   ! --- Align the two SkyFrames themselves ---
   ! A SkyFrame's own axis order is NOT a fixed RA-then-Dec convention --
   ! it reflects whichever axis the header declared as longitude vs
   ! latitude first (confirmed by direct comparison: file B's CTYPE1=DEC
   ! makes its SkyFrame present (Dec,RA), not (RA,Dec) like file A's).
   ! astConvert between two whole FrameSets failed earlier because domain
   ! search only checks each FrameSet's own top-level registered frames --
   ! but skyframe_a/skyframe_b are genuine Frame objects (not FrameSets),
   ! so Frame-to-Frame astConvert applies here directly and correctly
   ! resolves any axis-order/equinox/system difference via AST's own
   ! SkyFrame alignment logic, rather than this code assuming an order.
   sky2sky = ast_convert(skyframe_a, skyframe_b, ' ', status)
   if (status.ne.0 .or. sky2sky.eq.ast__null) then
      write(*,*) 'ERROR: failed to align the two SkyFrames, status=', status
      stop 1
   endif

   ! --- Compose pixel_A -> sky_A -> sky_B -> pixel_B ---
   call ast_invert(skymap_b, status)
   pix2pix = ast_cmpmap(skymap_a, sky2sky, .true., ' ', status)
   pix2pix = ast_cmpmap(pix2pix, skymap_b, .true., ' ', status)
   call ast_invert(skymap_b, status)
   if (status.ne.0 .or. pix2pix.eq.ast__null) then
      write(*,*) 'ERROR: failed to compose the pixel(A)->pixel(B) Mapping,',&
      &' status=', status
      stop 1
   endif

   write(*,*) 'pixel(A) -> pixel(B) Mapping composed:'
   write(*,'(A,I0)') '  Nin  (pixel axes in A): ', ast_geti(pix2pix, 'Nin', status)
   write(*,'(A,I0)') '  Nout (pixel axes in B): ', ast_geti(pix2pix, 'Nout', status)

   xin(1) = 5.0d0
   yin(1) = 9.0d0
   call ast_tran2(pix2pix, 1, xin, yin, .true., xout, yout, status)
   write(*,'(A,F0.3,A,F0.3,A)') '  A pixel (5.0,9.0) -> B pixel (',&
   &xout(1), ' , ', yout(1), ')'

   call ast_annul(pix2pix, status)
   call ast_annul(sky2sky, status)
   call ast_annul(skymap_a, status)
   call ast_annul(skymap_b, status)
   call ast_annul(skyframe_a, status)
   call ast_annul(skyframe_b, status)
   call ast_annul(wcs_a, status)
   call ast_annul(wcs_b, status)
   call ast_end(status)

   if (status.ne.0) then
      write(*,*) 'ERROR: AST reported an error, final status=', status
      stop 1
   endif

   write(*,*) 'OK: cross-file pixel(A)->pixel(B) Mapping verified.'

contains

   subroutine extract_sky_mapping(wcs, skymap, skyframe, status)
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
      integer, intent(inout) :: status

      integer :: curframe, nout, i, j
      integer :: probe_axes(2), probe_frame, probe_map
      integer :: sky_axes_in(2), out_axes(4)
      integer :: fullmap, simplemap
      logical :: found_sky

      skymap = ast__null
      skyframe = ast__null
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

      call ast_annul(fullmap, status)
      call ast_annul(simplemap, status)
   end subroutine extract_sky_mapping

   subroutine load_wcs(filename, wcs, status)
      !! Read filename's FITS header via CFITSIO, load it into an AST
      !! FitsChan, and return the WCS FrameSet recovered from it.
      character(len=*), intent(in) :: filename
      integer, intent(out) :: wcs
      integer, intent(inout) :: status

      integer :: unit, blocksize, fitsstat, nkeys, nmore, i, fitschan
      character(len=80) :: card

      wcs = ast__null
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
