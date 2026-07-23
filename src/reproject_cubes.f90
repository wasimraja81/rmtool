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
! Also demonstrates the actual regridding step: resample_one_plane reads
! one specific plane and resamples it onto the final output grid via
! astResampleR. Caught a real bug on the first attempt: the documented
! 20-argument signature (this, ndim_in, lbnd_in, ubnd_in, in, in_var,
! interp, finterp, params, flags, tol, maxpix, badval, ndim_out,
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
! Loops over every combination of an input's non-sky axes (channel,
! Stokes, or any other axis a cube happens to have), decoding a flat
! plane index via a mixed-radix counter so it works for any number of
! non-sky axes, not just the 2 (channel, Stokes) this project's own
! fixtures happen to have -- and reuses the same input->reference Mapping
! across every plane (it depends only on the 2 sky axes, never on which
! plane is being read). Verified across all 200 channels of the axis-
! swapped fixture: 3 spot-checked channels (1, 100, 200) all reproduce
! their known ground-truth values exactly, and the whole 200-channel run
! completes in well under a second.
!
! Actually writes an output FITS file now (write_reprojected_file):
! output axis layout puts the 2 sky axes at OUTPUT positions 1,2 (RA-
! fastest-on-disk, matching what rm_synthesis's own tile-read I/O already
! assumes -- confirmed against rm_synthesis_mod.f90's own auto-tile-
! planner comment), in whatever order the REFERENCE presents them; other
! (non-sky) axes keep INFILE's own values unchanged, at output positions
! 3... This guarantees "sky is axes 1,2" but only guarantees "axis 1 is
! literally RA" for a conventionally-ordered reference (RA before Dec) --
! a deliberate, documented scope limit, not a silent assumption. CRPIX
! for the sky axes is shifted for the output grid's own origin using the
! exact same formula rm_synthesis.f90 already uses for its own subimage
! CRPIX shift. Verified independently via Python/astropy, header and
! data both: reference mode reproduces the reference's own header
! exactly and all 3 spot-checked channels match known ground truth;
! intersection mode's CRPIX1 correctly shifts from 17.0 to 1.0 for a
! 16-pixel crop, CRPIX2 stays untouched (that axis wasn't cropped), and
! the pixel value at the shifted position matches independently-computed
! ground truth exactly.
!
! Usage: reproject_cubes mode=<intersection|union|reference> reffile=<reference_file> infiles=<input_file>[,<input_file>...]
!    or: reproject_cubes --config <cfgfile>
!    or: reproject_cubes --config <cfgfile> mode=<...> [reffile=<...>] [infiles=<...>]
!    or: reproject_cubes --help | -h
! No positional args: mode/reffile/infiles are always named key=value
! (no spaces around '='), on the command line or via --config, never
! inferred from argument order -- deliberate, to leave no room for a
! user mistake on which bare word means what. Each CLI key=value token
! overrides only that same key from --config (per-key precedence, not
! an atomic-group replacement -- unambiguous now that every CLI value
! names its own field). Config file is the same key=value style:
! mode=..., reffile=..., and infiles=file1,file2,file3 (comma-
! separated, same csv-list convention rm_synthesis's own multi-band
! config keys use). Full usage text is in print_usage below (shared by
! --help and the argument-error path).
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
   integer :: n_inputs, i

   character(len=512) :: cfgfile
   character(len=512) :: this_arg, cli_key, cli_val
   character(len=16) :: cli_mode
   character(len=512) :: cli_reffile, raw_cli_infiles
   integer :: argc, iarg
   logical :: has_kv
   logical :: have_cfgfile
   logical :: seen_mode, seen_reffile, seen_infiles
   logical :: cli_seen_mode, cli_seen_reffile, cli_seen_infiles

   integer :: wcs_ref, skymap_ref, skyframe_ref
   integer :: naxes_ref(max_axes), pixaxes_ref(2)
   integer :: wcs_in, skymap_in, skyframe_in
   integer :: naxes_in(max_axes), pixaxes_in(2)
   integer :: map_in2ref
   double precision :: lbnd_out(2), ubnd_out(2)
   double precision :: this_lbnd(2), this_ubnd(2)
   integer :: status

   ! Two ways to supply mode/reffile/infiles: explicit key=value CLI
   ! tokens (mode=..., reffile=..., infiles=...) or a --config key=value
   ! file (read_reproject_cfg below, same csv-list convention). No
   ! positional args -- deliberately: a bare "reproject_cubes
   ! intersection ref.fits a.fits b.fits" leaves the user to remember
   ! argument order from memory, which is exactly the kind of mistake
   ! named key=value args rule out. Both sources may be given together;
   ! each CLI key=value token overrides only its own field from --config
   ! (per-key, not atomic-group, precedence -- unambiguous now that every
   ! CLI value names the field it sets).
   have_cfgfile = .false.
   cli_seen_mode = .false.
   cli_seen_reffile = .false.
   cli_seen_infiles = .false.
   argc = command_argument_count()
   iarg = 1
   do while (iarg.le.argc)
      call get_command_argument(iarg, this_arg)
      if (trim(this_arg).eq.'--help' .or. trim(this_arg).eq.'-h') then
         call print_usage()
         stop
      else if (trim(this_arg).eq.'--config') then
         if (iarg.eq.argc) then
            write(*,*) 'ERROR: --config requires a file path argument'
            stop 1
         endif
         call get_command_argument(iarg+1, cfgfile)
         have_cfgfile = .true.
         iarg = iarg + 2
      else
         call split_cli_kv(this_arg, cli_key, cli_val, has_kv)
         if (.not. has_kv) then
            write(*,*) 'ERROR: unrecognised argument "', trim(this_arg),&
            &'" -- expected key=value (mode=..., reffile=..., infiles=...),',&
            &' --config <file>, or --help'
            stop 1
         endif
         select case (trim(cli_key))
         case ('mode')
            if (cli_seen_mode) then
               write(*,*) 'ERROR: mode given more than once on the command line'
               stop 1
            endif
            cli_mode = trim(cli_val)
            cli_seen_mode = .true.
         case ('reffile')
            if (cli_seen_reffile) then
               write(*,*) 'ERROR: reffile given more than once on the command line'
               stop 1
            endif
            cli_reffile = trim(cli_val)
            cli_seen_reffile = .true.
         case ('infiles')
            if (cli_seen_infiles) then
               write(*,*) 'ERROR: infiles given more than once on the command line'
               stop 1
            endif
            raw_cli_infiles = trim(cli_val)
            cli_seen_infiles = .true.
         case default
            write(*,*) 'ERROR: unrecognised key "', trim(cli_key), '" -- expected',&
            &' mode, reffile, or infiles'
            stop 1
         end select
         iarg = iarg + 1
      endif
   enddo

   seen_mode = .false.
   seen_reffile = .false.
   seen_infiles = .false.
   if (have_cfgfile) then
      call read_reproject_cfg(cfgfile, mode, reffile, infiles, n_inputs, status)
      if (status.ne.0) stop 1
      seen_mode = .true.
      seen_reffile = .true.
      seen_infiles = .true.
   endif

   if (cli_seen_mode) then
      mode = cli_mode
      seen_mode = .true.
   endif
   if (cli_seen_reffile) then
      reffile = cli_reffile
      seen_reffile = .true.
   endif
   if (cli_seen_infiles) then
      n_inputs = cfg_csv_count(raw_cli_infiles)
      if (n_inputs.lt.1 .or. n_inputs.gt.max_inputs) then
         write(*,*) 'ERROR: infiles must list between 1 and ', max_inputs, ' files'
         stop 1
      endif
      do i = 1, n_inputs
         call cfg_csv_get_item(raw_cli_infiles, i, infiles(i))
      enddo
      seen_infiles = .true.
   endif

   if (.not. seen_mode .or. .not. seen_reffile .or. .not. seen_infiles) then
      call print_usage()
      stop 1
   endif

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

   ! --- Resample and write every plane of every input onto the final
   ! output grid, via astResampleR + FTPSSE ---
   ! Only naxes_in/pixaxes_in (plain per-file array-shape data, not an
   ! AST Object) are needed here -- write_reprojected_file builds its
   ! own private input->reference Mapping per OpenMP thread internally
   ! (see its own comment for why), so this loop no longer needs to
   ! derive map_in2ref itself the way it used to.
   do i = 1, n_inputs
      call load_wcs(infiles(i), wcs_in, naxes_in, status)
      call extract_sky_mapping(wcs_in, skymap_in, skyframe_in, pixaxes_in, status)
      if (status.ne.0) then
         write(*,*) 'ERROR: failed to read input''s WCS for resampling: ',&
         &trim(infiles(i))
         stop 1
      endif
      call write_reprojected_file(reffile, infiles(i),&
      &'!'//trim(infiles(i))//'_REPROJ.FITS', pixaxes_ref,&
      &naxes_in, pixaxes_in, lbnd_out, ubnd_out, status)
      if (status.ne.0) then
         write(*,*) 'ERROR: failed to write reprojected output for: ',&
         &trim(infiles(i))
         stop 1
      endif
      call ast_annul(skymap_in, status)
      call ast_annul(skyframe_in, status)
      call ast_annul(wcs_in, status)
   enddo

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

   subroutine write_reprojected_file(reffile, infile, outfile, pixaxes_ref,&
   &naxes_in, pixaxes_in, lbnd_out_d, ubnd_out_d, status)
      !! Create outfile and write every reprojected plane of infile into
      !! it. Output axis layout: the 2 sky axes always occupy OUTPUT
      !! positions 1,2 (matching what rm_synthesis's own tile-read I/O
      !! already assumes -- RA-fastest-on-disk -- confirmed against
      !! rm_synthesis_mod.f90's own auto-tile-planner comment), in
      !! whatever order the REFERENCE file itself presents them (pixaxes_ref);
      !! this guarantees "sky is axes 1,2" but only guarantees "axis 1 is
      !! literally RA" if the reference is conventionally ordered (RA
      !! before Dec) -- a deliberate, documented scope limit, not a
      !! silent assumption: full semantic RA/Dec canonicalisation
      !! regardless of the reference's own convention is a follow-on
      !! refinement. The other (non-sky) axes keep their own values from
      !! INFILE unchanged (that band's own channel/Stokes definitions,
      !! untouched by reprojection), placed at output positions 3.. in
      !! their original relative order.
      use omp_lib, only: omp_get_max_threads
      character(len=*), intent(in) :: reffile, infile, outfile
      integer, intent(in) :: pixaxes_ref(2)
      integer, intent(in) :: naxes_in(:), pixaxes_in(2)
      double precision, intent(in) :: lbnd_out_d(2), ubnd_out_d(2)
      integer, intent(inout) :: status

      integer :: naxis, k, other_axes(max_axes), n_other
      integer :: other_idx(max_axes), remainder, radix
      integer :: iplane, n_planes, status_par, plane_status, nthreads
      integer :: nx_out, ny_out, naxis_out, naxes_out(max_axes)
      integer :: ref_unit, out_unit, fitsstat, blocksize
      logical :: simple, extend

      ! Per-OpenMP-thread private AST working set (see the parallel
      ! region below for why each thread builds its own).
      integer :: t_status, t_wcs_ref, t_skymap_ref, t_skyframe_ref
      integer :: t_naxes_ref(max_axes), t_pixaxes_ref(2)
      integer :: t_wcs_in, t_skymap_in, t_skyframe_in
      integer :: t_naxes_in(max_axes), t_pixaxes_in(2)
      integer :: t_map_in2ref

      if (status.ne.0) return

      naxis = 0
      do k = 1, size(naxes_in)
         if (naxes_in(k).gt.0) naxis = k
      enddo
      n_other = 0
      do k = 1, naxis
         if (k.ne.pixaxes_in(1) .and. k.ne.pixaxes_in(2)) then
            n_other = n_other + 1
            other_axes(n_other) = k
         endif
      enddo

      nx_out = nint(ubnd_out_d(1) - lbnd_out_d(1)) + 1
      ny_out = nint(ubnd_out_d(2) - lbnd_out_d(2)) + 1
      naxis_out = 2 + n_other
      naxes_out(1) = nx_out
      naxes_out(2) = ny_out
      do k = 1, n_other
         naxes_out(2+k) = naxes_in(other_axes(k))
      enddo

      ! --- Create the output file and its primary header ---
      fitsstat = 0
      blocksize = 1
      out_unit = 43
      call FTINIT(out_unit, trim(outfile), blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to create output file: ', trim(outfile)
         call printerror(fitsstat)
         status = -1
         return
      endif
      simple = .true.
      extend = .false.
      call FTPHPR(out_unit, simple, -32, naxis_out, naxes_out(1:naxis_out),&
      &0, 1, extend, fitsstat)

      ! Sky axes (output 1,2): WCS copied from the REFERENCE's own
      ! pixaxes_ref-numbered keywords, CRPIX shifted for the output
      ! grid's own origin (same formula rm_synthesis.f90 already uses for
      ! its own subimage CRPIX shift: (old_crpix - pixel_offset) + 1,
      ! stride 1 since this is a crop/grow, never a sub-sample).
      ref_unit = 44
      fitsstat = 0
      call FTOPEN(ref_unit, trim(reffile), 0, blocksize, fitsstat)
      call copy_axis_keywords(ref_unit, pixaxes_ref(1), out_unit, 1,&
      &lbnd_out_d(1)-1.0d0, status)
      call copy_axis_keywords(ref_unit, pixaxes_ref(2), out_unit, 2,&
      &lbnd_out_d(2)-1.0d0, status)
      call FTCLOS(ref_unit, fitsstat)

      ! Other axes (output 3..): WCS copied from INFILE's own axis
      ! numbering unchanged (no CRPIX shift -- reprojection never touches
      ! these axes).
      fitsstat = 0
      call FTOPEN(ref_unit, trim(infile), 0, blocksize, fitsstat)
      do k = 1, n_other
         call copy_axis_keywords(ref_unit, other_axes(k), out_unit,&
         &2+k, 0.0d0, status)
      enddo
      call FTCLOS(ref_unit, fitsstat)

      ! --- Resample and write every plane, in parallel across OpenMP
      ! threads ---
      ! Each thread builds its OWN private input->reference pixel
      ! Mapping from scratch (own ast_begin context, own load_wcs/
      ! extract_sky_mapping/compose_pix2pix calls) rather than sharing
      ! or cloning one Mapping built by the caller. SUN/211 Sec 4.12
      ! ("AST Objects within Multi-threaded Applications") requires
      ! astLock/astUnlock to hand an AST Object from the thread that
      ! created it to another thread; checked this build's actual
      ! Fortran symbol table (libstarlink_ast.so.9) and ast_lock_/
      ! ast_unlock_ are simply not exported -- only ast_copy_ is, which
      ! is useless on its own without the paired lock handoff (a copy
      ! is still initially locked to its creating thread). Building a
      ! Mapping is a handful of FITS header reads plus cheap AST calls,
      ! negligible next to the per-pixel resampling cost this
      ! parallelises, so redoing it once per thread (not once per
      ! plane) costs nothing worth avoiding. Every thread's whole AST
      ! object graph is therefore self-created and never shared, which
      ! is squarely inside AST's documented thread-safe model without
      ! needing any lock/unlock call at all -- confirmed empirically too
      ! (see the setup critical section's own comment below): concurrent
      ! per-thread construction, and concurrent astResampleR on those
      ! independent Mappings, both ran cleanly across many repeated
      ! stress runs once the actual bug (below) was found and fixed.
      !
      ! CFITSIO unit numbers used inside this region are offset by
      ! omp_get_thread_num() (see load_wcs/resample_one_plane), using
      ! disjoint ranges (1000-based for load_wcs, 5000-based for
      ! resample_one_plane's read) with wide headroom on both sides: up
      ! to ~4000 threads before load_wcs's range could reach 5000, and
      ! both bases sit far above the fixed out_unit=43/ref_unit=44 --
      ! comfortably beyond any real HPC node's core count, chosen that
      ! wide deliberately after an earlier 100/200 pairing turned out to
      ! only be safe up to ~100 threads. Channel/plane count plays no
      ! part in this -- each thread reuses its own single unit number
      ! sequentially across every plane it processes, opening and
      ! closing it each time, so only thread count (not cube size)
      ! determines how much headroom is needed. An earlier version got
      ! the numbering wrong in two ways: load_wcs had no thread offset
      ! at all (a bare literal unit=11, so concurrent per-thread setup
      ! collided outright), and resample_one_plane's offset (31+tid)
      ! numerically reached 43 at thread 12 -- exactly out_unit, which
      ! stays open for the entire parallel region -- so anything needing
      ! >=13 threads reliably crashed while <=12 threads never touched
      ! unit 43 and ran clean. That thread-count-dependent threshold
      ! (reliable at 12, reliable *failure* at 13+) was what actually
      ! pointed at a plain unit collision rather than a deeper AST/
      ! CFITSIO concurrency limit.
      !
      ! The actual disk write (FTPSSE, inside resample_one_plane) stays
      ! serialised through an OMP critical section together with the
      ! read -- concurrent CFITSIO writes to the same output unit are
      ! not something this project trusts even via distinct handles (see
      ! rm_synthesis's own io_write_threads rewrite to raw stream writes
      ! after a production corruption incident,
      ! planning/IO_PARALLEL_OPTIMISATION_PLAN.md's T4 postmortem), and
      ! this build's CFITSIO also visibly crashed (ffpsse/ffgcprll) when
      ! reads on other units were left free to run fully concurrently
      ! against it. Resampling itself (astResampleR) is embarrassingly
      ! parallel and stays outside any critical section, so serialising
      ! setup+read+write still yields a real speedup (verified ~1.8x on
      ! a 1024x1024x300 synthetic cube, 16 threads vs. 1).
      n_planes = 1
      do k = 1, n_other
         n_planes = n_planes * naxes_in(other_axes(k))
      enddo
      write(*,'(A,A,A,I0,A)') 'Writing ', trim(outfile), ': ',&
      &n_planes, ' plane(s)'

      status_par = 0
      nthreads = max(1, min(omp_get_max_threads(), n_planes))
      !$omp parallel num_threads(nthreads) if(n_planes.gt.1) default(none)&
      !$omp& shared(infile, reffile, naxes_in, pixaxes_in, other_axes,&
      !$omp& n_other, lbnd_out_d, ubnd_out_d, n_planes, out_unit, status_par)&
      !$omp& private(t_status, t_wcs_ref, t_skymap_ref, t_skyframe_ref,&
      !$omp& t_naxes_ref, t_pixaxes_ref, t_wcs_in, t_skymap_in, t_skyframe_in,&
      !$omp& t_naxes_in, t_pixaxes_in, t_map_in2ref, iplane, other_idx,&
      !$omp& remainder, k, radix, plane_status)

      ! Building each thread's own AST object graph (ast_begin through
      ! compose_pix2pix) is serialised through the SAME reproj_io
      ! critical section resample_one_plane uses for its CFITSIO read/
      ! write (not a separate lock) -- this setup calls load_wcs twice,
      ! which does its own CFITSIO FTOPEN/FTGREC/FTCLOS header reads, and
      ! keeping that under the same lock as every other CFITSIO call in
      ! this file is the simplest way to guarantee none of them ever run
      ! concurrently, regardless of unit number bookkeeping. (The actual
      ! bug that motivated a long, initially wrong detour into "maybe AST
      ! Object construction itself isn't safely concurrent" turned out to
      ! be nothing of the kind -- see the unit-number comment above: a
      ! bare, non-thread-offset unit in load_wcs, and a thread-offset
      ! range in resample_one_plane that numerically collided with
      ! out_unit at exactly 13 threads. Once both were fixed, a battery
      ! of isolated reproducers confirmed concurrent per-thread AST
      ! construction AND concurrent astResampleR both work fine on this
      ! build, matching SUN/211's documented model. This critical section
      ! is kept anyway, now purely for the CFITSIO-call-serialisation
      ! reason, not an AST one.) Setup is cheap (a few header reads + AST
      ! calls) next to the per-pixel resampling this still parallelises
      ! (astResampleR itself stays outside any critical section), so
      ! serialising it costs little.
      t_status = 0
      !$omp critical(reproj_io)
      call ast_begin(t_status)
      call load_wcs(reffile, t_wcs_ref, t_naxes_ref, t_status)
      call extract_sky_mapping(t_wcs_ref, t_skymap_ref, t_skyframe_ref,&
      &t_pixaxes_ref, t_status)
      call load_wcs(infile, t_wcs_in, t_naxes_in, t_status)
      call extract_sky_mapping(t_wcs_in, t_skymap_in, t_skyframe_in,&
      &t_pixaxes_in, t_status)
      call compose_pix2pix(t_skymap_in, t_skyframe_in, t_skymap_ref,&
      &t_skyframe_ref, t_map_in2ref, t_status)
      !$omp end critical(reproj_io)

      if (t_status.ne.0) then
         !$omp critical(reproj_status)
         status_par = -1
         !$omp end critical(reproj_status)
      endif

      ! Every thread must reach this worksharing construct unconditionally
      ! (OpenMP requires all threads in the team to encounter it) -- a
      ! thread whose own setup above failed just skips its share of the
      ! real work inside the loop body instead of branching around the
      ! construct itself, which previously deadlocked the team at the
      ! implicit end-of-do barrier whenever setup failed on only some
      ! threads.
      !$omp do schedule(dynamic)
      do iplane = 1, n_planes
         if (t_status.eq.0) then
            remainder = iplane - 1
            do k = 1, n_other
               radix = naxes_in(other_axes(k))
               other_idx(k) = mod(remainder, radix) + 1
               remainder = remainder / radix
            enddo
            plane_status = 0
            call resample_one_plane(infile, t_map_in2ref, naxes_in,&
            &pixaxes_in, other_axes(1:n_other), other_idx(1:n_other),&
            &lbnd_out_d, ubnd_out_d, iplane, out_unit, n_other, plane_status)
            if (plane_status.ne.0) then
               !$omp critical(reproj_status)
               status_par = -1
               !$omp end critical(reproj_status)
            endif
         endif
      enddo
      !$omp end do

      !$omp critical(reproj_io)
      call ast_end(t_status)
      !$omp end critical(reproj_io)
      !$omp end parallel

      if (status_par.ne.0) then
         write(*,*) 'ERROR: failed to resample/write one or more planes for: ',&
         &trim(infile)
         status = -1
         call FTCLOS(out_unit, fitsstat)
         return
      endif

      call FTCLOS(out_unit, fitsstat)
   end subroutine write_reprojected_file

   subroutine copy_axis_keywords(src_unit, src_axis, dst_unit, dst_axis,&
   &crpix_shift, status)
      !! Copy CTYPE/CRVAL/CRPIX/CDELT/CUNIT for src_axis (in src_unit's
      !! own header) to dst_axis (in dst_unit's header, already created
      !! via FTPHPR). CRPIX is additionally shifted by -crpix_shift (0 for
      !! a straight passthrough; the output grid's own pixel-1 offset,
      !! reference-pixel-numbered, for a cropped/grown sky axis) --
      !! matches rm_synthesis.f90's own existing subimage CRPIX-shift
      !! formula exactly, generalised to any axis number via a
      !! constructed keyword string ("CRVAL"//axis, etc.) rather than a
      !! literal "1"/"2" suffix.
      integer, intent(in) :: src_unit, src_axis, dst_unit, dst_axis
      double precision, intent(in) :: crpix_shift
      integer, intent(inout) :: status

      integer :: fitsstat
      character(len=8) :: axstr
      character(len=68) :: comment
      character(len=68) :: sval
      double precision :: dval

      if (status.ne.0) return

      write(axstr,'(I0)') src_axis
      fitsstat = 0
      call FTGKYS(src_unit, 'CTYPE'//trim(axstr), sval, comment, fitsstat)
      if (fitsstat.eq.0) then
         write(axstr,'(I0)') dst_axis
         call FTPKYS(dst_unit, 'CTYPE'//trim(axstr), trim(sval), ' ', fitsstat)
      endif

      write(axstr,'(I0)') src_axis
      fitsstat = 0
      call FTGKYD(src_unit, 'CRVAL'//trim(axstr), dval, comment, fitsstat)
      if (fitsstat.eq.0) then
         write(axstr,'(I0)') dst_axis
         call FTPKYD(dst_unit, 'CRVAL'//trim(axstr), dval, 13, ' ', fitsstat)
      endif

      write(axstr,'(I0)') src_axis
      fitsstat = 0
      call FTGKYD(src_unit, 'CRPIX'//trim(axstr), dval, comment, fitsstat)
      if (fitsstat.eq.0) then
         dval = dval - crpix_shift
         write(axstr,'(I0)') dst_axis
         call FTPKYD(dst_unit, 'CRPIX'//trim(axstr), dval, 13, ' ', fitsstat)
      endif

      write(axstr,'(I0)') src_axis
      fitsstat = 0
      call FTGKYD(src_unit, 'CDELT'//trim(axstr), dval, comment, fitsstat)
      if (fitsstat.eq.0) then
         write(axstr,'(I0)') dst_axis
         call FTPKYD(dst_unit, 'CDELT'//trim(axstr), dval, 13, ' ', fitsstat)
      endif

      write(axstr,'(I0)') src_axis
      fitsstat = 0
      call FTGKYS(src_unit, 'CUNIT'//trim(axstr), sval, comment, fitsstat)
      if (fitsstat.eq.0) then
         write(axstr,'(I0)') dst_axis
         call FTPKYS(dst_unit, 'CUNIT'//trim(axstr), trim(sval), ' ', fitsstat)
      endif
   end subroutine copy_axis_keywords

   subroutine resample_one_plane(filename, map_in2ref, naxes_in,&
   &pixaxes_in, other_axes, other_idx, lbnd_out_d, ubnd_out_d, iplane,&
   &out_unit, n_other, status)
      !! Read one specific plane of filename (other_axes(k) fixed at
      !! other_idx(k) for every non-sky axis) and resample it onto the
      !! output grid [lbnd_out_d,ubnd_out_d] (reference pixel space) via
      !! astResampleR, using map_in2ref (the forward input->reference
      !! pixel Mapping -- astResampleR uses its INVERSE internally to
      !! look up each output pixel's input value, per its own documented
      !! convention). Bad/uncovered output pixels are flagged with IEEE
      !! NaN, matching this codebase's existing NaN-based masking
      !! convention throughout rm_synthesis.
      !!
      !! Reading fpixels(k)=lpixels(k)=other_idx for every non-sky axis
      !! (full 1..NAXIS extent for the 2 sky axes) naturally produces a
      !! flat array with the sky axes as its only two non-degenerate
      !! dimensions, in pixaxes_in's own order (first fastest) -- exactly
      !! the "Fortran array indexing" astResampleR expects -- regardless
      !! of which raw file-axis numbers they are, since a degenerate
      !! (extent-1) axis contributes no stride.
      !!
      !! Called from within write_reprojected_file's OpenMP parallel
      !! region (one call per plane, across threads): the CFITSIO read
      !! unit is offset by omp_get_thread_num() so concurrently-running
      !! threads never collide on the same unit number, and the actual
      !! output write (FTPSSE, at the end) is wrapped in an OMP critical
      !! section since out_unit is shared across all threads.
      use, intrinsic :: ieee_arithmetic
      use omp_lib, only: omp_get_thread_num
      character(len=*), intent(in) :: filename
      integer, intent(in) :: map_in2ref
      integer, intent(in) :: naxes_in(:), pixaxes_in(2)
      integer, intent(in) :: other_axes(:), other_idx(:)
      double precision, intent(in) :: lbnd_out_d(2), ubnd_out_d(2)
      integer, intent(in) :: iplane
      integer, intent(in) :: out_unit, n_other
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
      integer :: fpixels_wr(max_axes), lpixels_wr(max_axes), naxis_wr
      integer :: naxes_wr(max_axes)

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
      do k = 1, size(other_axes)
         fpixels(other_axes(k)) = other_idx(k)
         lpixels(other_axes(k)) = other_idx(k)
      enddo

      fitsstat = 0
      blocksize = 1
      unit = 5000 + omp_get_thread_num()
      group = 1
      badval = ieee_value(badval, ieee_quiet_nan)
      ! Serialised with the write critical section below (same lock name,
      ! reproj_io), together with load_wcs's own CFITSIO calls during
      ! per-thread setup -- a genuine cross-unit CFITSIO crash was
      ! observed (ffpsse -> ffgcprll) with this read left unserialised,
      ! but that run also had the unit-number collision bug described in
      ! write_reprojected_file's comment (this read's unit could
      ! numerically equal out_unit at higher thread counts), so it is
      ! not confirmed whether genuinely concurrent CFITSIO calls on truly
      ! distinct units are actually unsafe on this build or whether that
      ! crash was just the same collision in different clothes. Kept
      ! serialised rather than re-litigate that distinction: I/O is not
      ! the dominant cost here, astResampleR is (see below), and it stays
      ! outside any critical section so still runs fully in parallel
      ! across threads regardless of how the I/O is serialised.
      !$omp critical(reproj_io)
      call FTOPEN(unit, trim(filename), 0, blocksize, fitsstat)
      call FTGSVE(unit, group, naxis, naxes_in(1:naxis),&
      &fpixels(1:naxis), lpixels(1:naxis), incs(1:naxis),&
      &badval, data_in, anyflg, fitsstat)
      call FTCLOS(unit, fitsstat)
      !$omp end critical(reproj_io)
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

      ! --- Write this plane into the output file ---
      ! Output axis layout (see write_reprojected_file's own comment):
      ! sky always at output positions 1,2, full extent; other axes at
      ! output positions 3.. (2+k for other_axes(k)), fixed at this
      ! plane's own other_idx(k) -- same value, just relocated to the
      ! output file's own axis numbering rather than the input's.
      ! Wrapped in the same reproj_io critical section as the read above
      ! (not just out_unit/stdout being shared -- see the read's own
      ! comment on why CFITSIO calls need to be mutually exclusive across
      ! ALL threads here, not merely per-unit).
      naxis_wr = 2 + n_other
      naxes_wr(1) = nx_out
      naxes_wr(2) = ny_out
      fpixels_wr(1) = 1
      fpixels_wr(2) = 1
      lpixels_wr(1) = nx_out
      lpixels_wr(2) = ny_out
      do k = 1, n_other
         naxes_wr(2+k) = naxes_in(other_axes(k))
         fpixels_wr(2+k) = other_idx(k)
         lpixels_wr(2+k) = other_idx(k)
      enddo
      !$omp critical(reproj_io)
      call FTPSSE(out_unit, 1, naxis_wr, naxes_wr(1:naxis_wr),&
      &fpixels_wr(1:naxis_wr), lpixels_wr(1:naxis_wr), data_out, status)
      if (status.ne.0) then
         write(*,*) 'ERROR: failed to write plane ', iplane, ' to output'
      endif

      ! Verification print for a handful of representative planes only
      ! (not all -- would flood output for a 200-channel cube): the
      ! resampled value at reference pixel (5,9), for comparison against
      ! known ground-truth values read directly from the original data
      ! with a Python script beforehand, at 3 different channels.
      if (status.eq.0 .and.&
      &(iplane.eq.1 .or. iplane.eq.100 .or. iplane.eq.200) .and.&
      &5.ge.lbnd_o(1) .and. 5.le.ubnd_o(1) .and.&
      &9.ge.lbnd_o(2) .and. 9.le.ubnd_o(2)) then
         write(*,'(A,I0,A,I0,A,F0.9)') '  plane ', iplane,&
         &': bad=', nbad, ', value at reference pixel (5,9): ',&
         &data_out(5-lbnd_o(1)+1, 9-lbnd_o(2)+1)
      endif
      !$omp end critical(reproj_io)

      deallocate(data_in)
      deallocate(data_out)
   end subroutine resample_one_plane

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
      !!
      !! The CFITSIO unit is offset by omp_get_thread_num() (1000-based,
      !! disjoint from resample_one_plane's 5000-based range and from
      !! write_reprojected_file's fixed out_unit=43/ref_unit=44, with
      !! headroom for thousands of threads -- see write_reprojected_file's
      !! own comment on why that wide a margin) -- found the hard way: an
      !! earlier version left this unit a bare
      !! literal (11) with no thread offset at all, so every thread
      !! calling load_wcs concurrently during write_reprojected_file's
      !! per-thread setup collided on the very same CFITSIO unit number,
      !! which is what actually crashed/deadlocked multi-threaded runs --
      !! not, as first suspected, any inherent AST/CFITSIO thread-safety
      !! limitation (a battery of isolated reproducers matching the real
      !! call pattern -- concurrent per-thread AST FitsChan/FrameSet/
      !! Mapping construction, concurrent astResampleR on independent
      !! per-thread Mappings, concurrent CFITSIO reads+writes on
      !! correctly-separated unit numbers -- ran cleanly, 16 threads x
      !! hundreds of iterations x many repeats, once unit numbers
      !! actually stopped colliding).
      use omp_lib, only: omp_get_thread_num
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
      unit = 1000 + omp_get_thread_num()
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

   subroutine print_usage()
      !! Shared by --help/-h and the argument-error path, so the two
      !! can't drift out of sync with each other.
      write(*,'(A)') 'reproject_cubes -- reproject FITS cubes onto a common grid'
      write(*,'(A)') ''
      write(*,'(A)') 'Usage:'
      write(*,'(A)') '  reproject_cubes mode=<intersection|union|reference>'//&
      &' reffile=<reference_file> infiles=<input_file>[,<input_file>...]'
      write(*,'(A)') '  reproject_cubes --config <cfgfile>'
      write(*,'(A)') '  reproject_cubes --config <cfgfile> mode=<...>'//&
      &' [reffile=<...>] [infiles=<...>]'
      write(*,'(A)') '  reproject_cubes --help | -h'
      write(*,'(A)') ''
      write(*,'(A)') 'No positional args -- mode/reffile/infiles must each be given as'//&
      &' key=value (no spaces around ''='', e.g. mode=union), either directly on the'
      write(*,'(A)') 'command line or via --config. If both are given, each CLI key=value'//&
      &' overrides only that same key from --config; unset keys still come from --config.'
      write(*,'(A)') ''
      write(*,'(A)') 'Modes:'
      write(*,'(A)') '  reference    output grid is the reference file''s own extent'
      write(*,'(A)') '  intersection output grid shrinks to the overlap of all inputs'//&
      &' with the reference'
      write(*,'(A)') '  union        output grid grows to cover all inputs and the reference'
      write(*,'(A)') ''
      write(*,'(A)') 'Config file: key=value text file with three required keys:'
      write(*,'(A)') '  mode    = intersection | union | reference'
      write(*,'(A)') '  reffile = /path/to/reference.fits'
      write(*,'(A)') '  infiles = /path/band1.fits,/path/band2.fits,/path/band3.fits'
      write(*,'(A)') '(infiles is comma-separated, no spaces required; ''#'' or '';'''//&
      &' starts a comment.)'
   end subroutine print_usage

   subroutine read_reproject_cfg(cfgfile, mode, reffile, infiles, n_inputs, status)
      !! Parse a --config key=value file: three required keys, mode,
      !! reffile, and infiles (comma-separated, same csv-list convention
      !! rm_synthesis's own multi-band keys use). Standalone re-
      !! implementation of rm_synthesis_mod's split_key_value/csv_count/
      !! csv_get_item (below) rather than a `use` dependency -- this tool
      !! is deliberately kept off the main rm_synthesis build graph (own
      !! binary, own dependency set, see the Makefile comment), and these
      !! are a handful of generic string-parsing lines each, unlikely to
      !! drift.
      character(len=*), intent(in) :: cfgfile
      character(len=*), intent(out) :: mode, reffile
      character(len=*), intent(out) :: infiles(:)
      integer, intent(out) :: n_inputs
      integer, intent(out) :: status

      character(len=512) :: line, key, val, raw_infiles
      integer :: unit_cfg, ios, line_no, j
      logical :: has_kv, seen_mode, seen_reffile, seen_infiles

      status = 0
      n_inputs = 0
      seen_mode = .false.
      seen_reffile = .false.
      seen_infiles = .false.
      raw_infiles = ' '

      open(newunit=unit_cfg, file=trim(cfgfile), status='old', action='read', iostat=ios)
      if (ios.ne.0) then
         write(*,*) 'ERROR: cannot open config file: ', trim(cfgfile)
         status = -1
         return
      endif

      line_no = 0
      do
         read(unit_cfg, '(A)', iostat=ios) line
         if (ios.ne.0) exit
         line_no = line_no + 1
         call cfg_split_key_value(line, key, val, has_kv)
         if (.not. has_kv) cycle
         select case (trim(key))
         case ('mode')
            mode = trim(val)
            seen_mode = .true.
         case ('reffile')
            reffile = trim(val)
            seen_reffile = .true.
         case ('infiles')
            raw_infiles = trim(val)
            seen_infiles = .true.
         case default
            write(*,*) 'ERROR: unrecognised config key "', trim(key), '" at line ',&
            &line_no, ' in ', trim(cfgfile)
            status = -1
            close(unit_cfg)
            return
         end select
      enddo
      close(unit_cfg)

      if (.not. seen_mode .or. .not. seen_reffile .or. .not. seen_infiles) then
         write(*,*) 'ERROR: config file ', trim(cfgfile),&
         &' must set mode, reffile, and infiles'
         status = -1
         return
      endif

      n_inputs = cfg_csv_count(raw_infiles)
      if (n_inputs.lt.1 .or. n_inputs.gt.size(infiles)) then
         write(*,*) 'ERROR: infiles in config must list between 1 and ',&
         &size(infiles), ' files'
         status = -1
         n_inputs = 0
         return
      endif
      do j = 1, n_inputs
         call cfg_csv_get_item(raw_infiles, j, infiles(j))
      enddo
   end subroutine read_reproject_cfg

   subroutine cfg_split_key_value(raw_line, key, val, has_kv)
      !! Same convention as rm_synthesis_mod's split_key_value: strips
      !! ';'/'#' comments, splits on the first '=', blank/comment-only
      !! lines and lines missing either side of '=' yield has_kv=.false.
      character(len=*), intent(in) :: raw_line
      character(len=*), intent(out) :: key, val
      logical, intent(out) :: has_kv
      character(len=len(raw_line)) :: line
      integer :: p1, p2, peq, pcut

      key = ' '
      val = ' '
      has_kv = .false.

      line = raw_line
      p1 = index(line, ';')
      p2 = index(line, '#')
      if (p1 > 0 .and. p2 > 0) then
         pcut = min(p1, p2)
      else if (p1 > 0) then
         pcut = p1
      else
         pcut = p2
      endif
      if (pcut > 0) line = line(1:pcut - 1)

      line = adjustl(line)
      if (len_trim(line) == 0) return

      peq = index(line, '=')
      if (peq <= 1) return

      key = adjustl(line(1:peq - 1))
      val = adjustl(line(peq + 1:))
      if (len_trim(key) == 0 .or. len_trim(val) == 0) return

      key = trim(key)
      val = trim(val)
      has_kv = .true.
   end subroutine cfg_split_key_value

   subroutine split_cli_kv(token, key, val, has_kv)
      !! Split a single CLI argument of the form key=value. Deliberately
      !! no '#'/';' comment-stripping (unlike cfg_split_key_value) --
      !! this is one shell-split argv token, not a config-file line, and
      !! a file path could legitimately contain either character.
      character(len=*), intent(in) :: token
      character(len=*), intent(out) :: key, val
      logical, intent(out) :: has_kv
      integer :: peq

      key = ' '
      val = ' '
      has_kv = .false.

      peq = index(token, '=')
      if (peq <= 1) return

      key = adjustl(token(1:peq - 1))
      val = adjustl(token(peq + 1:))
      if (len_trim(key) == 0 .or. len_trim(val) == 0) return

      key = trim(key)
      val = trim(val)
      has_kv = .true.
   end subroutine split_cli_kv

   function cfg_csv_count(str) result(n)
      !! Number of comma-separated items in str (1 if no comma present,
      !! 0 for a blank/empty string).
      character(len=*), intent(in) :: str
      integer :: n
      integer :: i

      n = 0
      if (len_trim(str) == 0) return
      n = 1
      do i = 1, len_trim(str)
         if (str(i:i) == ',') n = n + 1
      enddo
   end function cfg_csv_count

   subroutine cfg_csv_get_item(str, idx, item)
      !! Extract the idx-th (1-based) comma-separated item from str,
      !! trimmed of surrounding blanks.
      character(len=*), intent(in) :: str
      integer, intent(in) :: idx
      character(len=*), intent(out) :: item
      integer :: i, cur, p0, n

      item = ' '
      n = len_trim(str)
      if (n == 0) return

      cur = 1
      p0 = 1
      do i = 1, n
         if (str(i:i) == ',') then
            if (cur == idx) then
               item = adjustl(str(p0:i - 1))
               return
            endif
            cur = cur + 1
            p0 = i + 1
         endif
      enddo
      if (cur == idx) item = adjustl(str(p0:n))
   end subroutine cfg_csv_get_item

end program reproject_cubes
