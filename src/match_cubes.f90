! match_cubes -- consolidates reproject_cubes (src/reproject_cubes.f90) and
! convolve_cubes (src/convolve_cubes.f90) into one tool that can run either
! stage alone, or both chained THROUGH MEMORY with no intermediate FITS file
! written to disk. Built for real multi-band data at real scale (200GB+
! cubes): reproject_cubes and convolve_cubes run back-to-back already
! produce a correct pipeline, but the intermediate _REPROJ.FITS file is
! written in full and immediately read back in full for no reason other
! than being two separate programs -- doubling disk I/O and disk space for
! an artifact nobody actually wants.
!
! The two existing standalone tools are NOT touched by this file (a
! deliberate choice, confirmed with the user): reproject_cubes.f90 and
! convolve_cubes.f90 remain fully independent, already-tested tools for
! anyone who wants just one stage without this one. This file therefore
! duplicates (adapts, not `use`s) the subroutines it needs from both,
! rather than extracting a shared module -- a real, accepted maintenance
! cost in exchange for zero regression risk to two already-shipped tools.
! rm_synthesis itself is out of scope here -- feeding it directly, without
! any intermediate file at all, is a separate, harder design question for
! later.
!
! Order matters, and not just for tidiness: convolving to the common
! target beam BEFORE reprojection low-pass-filters the image before
! astResampleR's linear interpolation ever touches it, so resampling
! operates on smooth, well-sampled data rather than a band's own native
! (possibly only marginally Nyquist-sampled) sharp PSF -- avoiding
! interpolation/aliasing error that convolving afterward cannot undo,
! since the error is already baked into the resampled pixel values by
! then. This is the same reasoning behind an anti-alias filter before
! downsampling in ordinary signal processing. It also usually costs less:
! in `union` footprint mode the reprojected output grid is larger than any
! input's own native grid, so convolving first does the expensive FFT work
! on the smaller native footprint rather than the inflated union one.
! Confirmed with the user: default chain order is convolve-then-reproject;
! reproject-then-convolve is also correct (not wrong), just not the
! default, and remains selectable.
!
! Axis-scope handling is deliberately asymmetric by stage. `stages=
! reproject` alone keeps reproject_cubes' own fully general N-dimensional
! "other axes" handling (any number of non-sky axes, e.g. Stokes AND
! frequency both varying) -- no new restriction versus today's standalone
! tool. `stages=convolve` or `stages=both` adopt convolve_cubes' own
! existing restriction instead (exactly 2 sky axes + 1 FREQ axis, every
! other axis degenerate) -- not a new limitation introduced here, but the
! scope gaussft_mod's own per-channel convolution already has today (it
! only knows how to convolve a plane indexed by frequency).
!
! Usage: match_cubes stages=reproject|convolve|both
!    [order=convolve_reproject|reproject_convolve]
!    infiles=<file1>[,<file2>...]
!    [footprint_mode=intersection|union|reference] [reffile=<reference_file>]
!    [beamfiles=<spec1>[,<spec2>...]] [badchan_file=<file>]
!    [target_bmaj=<arcsec> target_bmin=<arcsec> target_bpa=<deg>]
!    [max_common_bmaj=<arcsec>] [mem_frac_ram=<fraction>] [outsuffix=<suffix>]
!    [npts=<n>] [khachiyan_tol=<tol>]
!    or: match_cubes --config <cfgfile>
!    or: match_cubes --help | -h
! Full usage text in print_usage below (shared by --help and the
! argument-error path, same convention as reproject_cubes.f90/
! convolve_cubes.f90).
program match_cubes
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none
   ! AST_PAR (the vendor Fortran constants file, /usr/include/AST_PAR) is
   ! fixed-form Fortran 77 and cannot be `include`d into a free-form .f90
   ! file directly -- same issue documented in reproject_cubes.f90's own
   ! comment. Only the handful of symbols actually used are declared
   ! directly instead, matching AST_PAR's own declared types.
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
   integer, parameter :: max_channels = 20000

   character(len=16) :: stages
   character(len=32) :: order
   logical :: do_reproject, do_convolve, convolve_first

   character(len=512) :: infiles(max_inputs), beamfiles(max_inputs)
   integer :: n_inputs
   character(len=64) :: outsuffix
   logical :: seen_outsuffix

   character(len=16) :: footprint_mode
   character(len=512) :: reffile
   logical :: seen_footprint_mode, seen_reffile

   character(len=512) :: badchan_file
   logical :: have_badchan_file
   logical :: have_target
   real(dp) :: target_bmaj, target_bmin, target_bpa
   real(dp) :: max_common_bmaj
   logical :: have_max_common_bmaj
   real :: mem_frac_ram
   integer :: npts
   real(dp) :: khachiyan_tol

   integer :: i, status

   ! --- Convolve-stage per-file bookkeeping (used whenever do_convolve) ---
   integer :: naxis_f(max_inputs), sky1_f(max_inputs), sky2_f(max_inputs)
   integer :: freq_axis_f(max_inputs), naxes_f(max_inputs, max_axes)
   real(dp) :: cdelt1_f(max_inputs), cdelt2_f(max_inputs)
   integer :: nfreq_f(max_inputs)
   real(dp), allocatable :: bmaj_f(:,:), bmin_f(:,:), bpa_f(:,:)
   logical, allocatable :: isbad_f(:,:)
   real(dp), allocatable :: pool_bmaj(:), pool_bmin(:), pool_bpa(:)
   integer :: n_pool
   real(dp) :: common_bmaj, common_bmin, common_bpa
   integer :: badchan_list(max_channels), n_badchan

   ! --- Reproject-stage bookkeeping (used whenever do_reproject) ---
   integer :: wcs_ref, skymap_ref, skyframe_ref
   integer :: naxes_ref(max_axes), pixaxes_ref(2)
   integer :: wcs_in, skymap_in, skyframe_in
   integer :: naxes_in(max_axes), pixaxes_in(2)
   integer :: map_in2ref
   double precision :: lbnd_out(2), ubnd_out(2)
   double precision :: this_lbnd(2), this_ubnd(2)
   integer :: ast_status
   integer :: nx_out_common, ny_out_common

   call parse_args(status)
   if (status.ne.0) stop 1

   do_reproject = (trim(stages).eq.'reproject' .or. trim(stages).eq.'both')
   do_convolve = (trim(stages).eq.'convolve' .or. trim(stages).eq.'both')
   convolve_first = (trim(order).eq.'convolve_reproject')

   ! === Pre-scan phase 1: reproject footprint (order-independent of the
   ! convolve pre-scan below; only needs WCS/geometry, not beam metadata) ===
   nx_out_common = 0
   ny_out_common = 0
   if (do_reproject) then
      ast_status = 0
      call ast_begin(ast_status)

      call load_wcs(reffile, wcs_ref, naxes_ref, ast_status)
      call extract_sky_mapping(wcs_ref, skymap_ref, skyframe_ref, pixaxes_ref, ast_status)
      if (ast_status.ne.0) then
         write(*,*) 'ERROR: failed to load the reference file''s WCS'
         stop 1
      endif

      lbnd_out(1) = 1.0d0
      lbnd_out(2) = 1.0d0
      ubnd_out(1) = real(naxes_ref(pixaxes_ref(1)), kind=8)
      ubnd_out(2) = real(naxes_ref(pixaxes_ref(2)), kind=8)
      write(*,'(A,A,A,F0.0,A,F0.0,A,F0.0,A,F0.0,A)') 'Reference (', trim(reffile),&
      &') own extent: [', lbnd_out(1), ',', ubnd_out(1), '] x [',&
      &lbnd_out(2), ',', ubnd_out(2), ']'

      if (trim(footprint_mode).ne.'reference') then
         do i = 1, n_inputs
            call load_wcs(infiles(i), wcs_in, naxes_in, ast_status)
            call extract_sky_mapping(wcs_in, skymap_in, skyframe_in, pixaxes_in, ast_status)
            if (ast_status.ne.0) then
               write(*,*) 'ERROR: failed to load input file: ', trim(infiles(i))
               stop 1
            endif

            call compose_pix2pix(skymap_in, skyframe_in, skymap_ref, skyframe_ref,&
            &map_in2ref, ast_status)
            if (ast_status.ne.0) then
               write(*,*) 'ERROR: failed to align input file to the reference: ',&
               &trim(infiles(i))
               stop 1
            endif

            call footprint_bounds(map_in2ref, naxes_in, pixaxes_in,&
            &this_lbnd, this_ubnd, ast_status)
            write(*,'(A,A,A,F0.2,A,F0.2,A,F0.2,A,F0.2,A)') '  ', trim(infiles(i)),&
            &' footprint in reference space: [', this_lbnd(1), ',', this_ubnd(1),&
            &'] x [', this_lbnd(2), ',', this_ubnd(2), ']'

            if (this_ubnd(1).lt.lbnd_out(1) .or. this_lbnd(1).gt.ubnd_out(1) .or.&
            &this_ubnd(2).lt.lbnd_out(2) .or. this_lbnd(2).gt.ubnd_out(2)) then
               write(*,*) 'ERROR: zero sky overlap between the reference and: ',&
               &trim(infiles(i))
               write(*,*) 'Quitting now...'
               stop 1
            endif

            if (trim(footprint_mode).eq.'intersection') then
               lbnd_out = max(lbnd_out, this_lbnd)
               ubnd_out = min(ubnd_out, this_ubnd)
            else ! union
               lbnd_out = min(lbnd_out, this_lbnd)
               ubnd_out = max(ubnd_out, this_ubnd)
            endif

            call ast_annul(map_in2ref, ast_status)
            call ast_annul(skymap_in, ast_status)
            call ast_annul(skyframe_in, ast_status)
            call ast_annul(wcs_in, ast_status)
         enddo
      endif

      if (trim(footprint_mode).eq.'intersection') then
         lbnd_out = ceiling(lbnd_out)
         ubnd_out = floor(ubnd_out)
      else
         lbnd_out = floor(lbnd_out)
         ubnd_out = ceiling(ubnd_out)
      endif

      if (lbnd_out(1).gt.ubnd_out(1) .or. lbnd_out(2).gt.ubnd_out(2)) then
         write(*,*) 'ERROR: computed output grid is empty (', trim(footprint_mode), ' mode)'
         stop 1
      endif

      nx_out_common = nint(ubnd_out(1) - lbnd_out(1)) + 1
      ny_out_common = nint(ubnd_out(2) - lbnd_out(2)) + 1
      write(*,'(A,A,A,F0.0,A,F0.0,A,F0.0,A,F0.0,A)') 'Final output grid (',&
      &trim(footprint_mode), ' mode): [', lbnd_out(1), ',', ubnd_out(1), '] x [',&
      &lbnd_out(2), ',', ubnd_out(2), ']'
   endif

   ! === Pre-scan phase 2: convolve beam metadata (independent of phase 1) ===
   if (do_convolve) then
      allocate(bmaj_f(max_inputs, max_channels), bmin_f(max_inputs, max_channels))
      allocate(bpa_f(max_inputs, max_channels), isbad_f(max_inputs, max_channels))

      n_badchan = 0
      if (have_badchan_file) then
         call read_badchan_file(badchan_file, badchan_list, n_badchan, status)
         if (status.ne.0) stop 1
      endif

      do i = 1, n_inputs
         call read_axis_info(infiles(i), naxis_f(i), sky1_f(i), sky2_f(i),&
         &freq_axis_f(i), naxes_f(i,:), cdelt1_f(i), cdelt2_f(i), status)
         if (status.ne.0) then
            write(*,*) 'ERROR: failed to read axis info for: ', trim(infiles(i))
            stop 1
         endif
         nfreq_f(i) = naxes_f(i, freq_axis_f(i))
         if (nfreq_f(i).gt.max_channels) then
            write(*,*) 'ERROR: ', trim(infiles(i)), ' has ', nfreq_f(i),&
            &' channels, exceeding this program''s max_channels=', max_channels
            stop 1
         endif

         call read_beams(infiles(i), trim(beamfiles(i)), nfreq_f(i),&
         &bmaj_f(i,1:nfreq_f(i)), bmin_f(i,1:nfreq_f(i)), bpa_f(i,1:nfreq_f(i)),&
         &isbad_f(i,1:nfreq_f(i)), status)
         if (status.ne.0) then
            write(*,*) 'ERROR: failed to read per-channel beams for: ', trim(infiles(i))
            stop 1
         endif

         call apply_badchan_list(badchan_list, n_badchan, nfreq_f(i), isbad_f(i,1:nfreq_f(i)))

         write(*,'(A,A,A,I0,A,I0,A)') 'Read ', trim(infiles(i)), ': ', nfreq_f(i),&
         &' channels, ', count(isbad_f(i,1:nfreq_f(i))), ' flagged bad'
      enddo

      if (have_target) then
         common_bmaj = target_bmaj
         common_bmin = target_bmin
         common_bpa = target_bpa
         write(*,'(A,F0.4,A,F0.4,A,F0.4)') 'Using explicit target beam: BMAJ=',&
         &common_bmaj, ' BMIN=', common_bmin, ' PA=', common_bpa
      else
         n_pool = 0
         do i = 1, n_inputs
            n_pool = n_pool + count(.not. isbad_f(i,1:nfreq_f(i)))
         enddo
         if (n_pool.lt.1) then
            write(*,*) 'ERROR: no good (non-bad) channels across any input file'
            stop 1
         endif
         allocate(pool_bmaj(n_pool), pool_bmin(n_pool), pool_bpa(n_pool))
         call pool_good_beams(n_inputs, nfreq_f, bmaj_f, bmin_f, bpa_f, isbad_f,&
         &max_inputs, max_channels, n_pool, pool_bmaj, pool_bmin, pool_bpa)

         call find_common_beam_wrap(n_pool, pool_bmaj, pool_bmin, pool_bpa,&
         &npts, khachiyan_tol, common_bmaj, common_bmin, common_bpa, status)
         if (status.ne.0) then
            write(*,*) 'ERROR: could not find a common beam deconvolvable from',&
            &' every good input channel'
            stop 1
         endif
         write(*,'(A,I0,A,F0.4,A,F0.4,A,F0.4)') 'Derived common beam from ',&
         &n_pool, ' good channels: BMAJ=', common_bmaj, ' BMIN=', common_bmin,&
         &' PA=', common_bpa
         deallocate(pool_bmaj, pool_bmin, pool_bpa)

         if (have_max_common_bmaj .and. common_bmaj.gt.max_common_bmaj) then
            write(*,'(A,F0.4,A,F0.4,A)') 'ERROR: derived common beam BMAJ=',&
            &common_bmaj, ' arcsec exceeds max_common_bmaj=', max_common_bmaj,&
            &' arcsec -- refusing to proceed. Investigate why the required'//&
            &' common resolution is this coarse (e.g. an outlier per-channel'//&
            &' beam not already flagged bad), or raise max_common_bmaj if this'//&
            &' resolution is genuinely intended.'
            stop 1
         endif
      endif
   endif

   ! === Per-file processing ===
   do i = 1, n_inputs
      if (.not. do_convolve) then
         ! stages=reproject alone: fully general N-dimensional axis
         ! handling, unrestricted -- exactly today's standalone
         ! reproject_cubes behaviour.
         call load_wcs(infiles(i), wcs_in, naxes_in, ast_status)
         call extract_sky_mapping(wcs_in, skymap_in, skyframe_in, pixaxes_in, ast_status)
         if (ast_status.ne.0) then
            write(*,*) 'ERROR: failed to read input''s WCS for resampling: ',&
            &trim(infiles(i))
            stop 1
         endif
         call process_one_file_general(reffile, infiles(i),&
         &'!'//trim(infiles(i))//trim(outsuffix), pixaxes_ref,&
         &naxes_in, pixaxes_in, lbnd_out, ubnd_out, mem_frac_ram, status)
         if (status.ne.0) then
            write(*,*) 'ERROR: failed to write reprojected output for: ',&
            &trim(infiles(i))
            stop 1
         endif
         call ast_annul(skymap_in, ast_status)
         call ast_annul(skyframe_in, ast_status)
         call ast_annul(wcs_in, ast_status)
      else
         ! stages=convolve or stages=both: restricted 2-sky+1-freq axis
         ! handling (convolve_cubes' own existing scope).
         call process_one_file_restricted(infiles(i),&
         &trim(infiles(i))//trim(outsuffix), do_reproject, convolve_first,&
         &naxis_f(i), sky1_f(i), sky2_f(i), freq_axis_f(i), naxes_f(i,:),&
         &cdelt1_f(i), cdelt2_f(i), nfreq_f(i), bmaj_f(i,1:nfreq_f(i)),&
         &bmin_f(i,1:nfreq_f(i)), bpa_f(i,1:nfreq_f(i)), isbad_f(i,1:nfreq_f(i)),&
         &common_bmaj, common_bmin, common_bpa, reffile, pixaxes_ref,&
         &nx_out_common, ny_out_common, lbnd_out, ubnd_out, mem_frac_ram, status)
         if (status.ne.0) then
            write(*,*) 'ERROR: failed to write output for: ', trim(infiles(i))
            stop 1
         endif
      endif
      write(*,*) 'OK: wrote ', trim(infiles(i))//trim(outsuffix)
   enddo

   if (do_reproject) then
      call ast_annul(skymap_ref, ast_status)
      call ast_annul(skyframe_ref, ast_status)
      call ast_annul(wcs_ref, ast_status)
      call ast_end(ast_status)
   endif

   if (allocated(bmaj_f)) deallocate(bmaj_f, bmin_f, bpa_f, isbad_f)
   write(*,*) 'OK: all inputs processed.'

contains

   !===========================================================
   ! CLI / config parsing (adapted from convolve_cubes.f90, extended
   ! with reproject-stage keys)
   !===========================================================

   subroutine parse_args(status)
      integer, intent(out) :: status
      character(len=512) :: this_arg, cli_key, cli_val, cfgfile
      character(len=512) :: raw_infiles, raw_beamfiles
      integer :: argc, iarg
      logical :: has_kv, have_cfgfile, seen_infiles, seen_stages

      status = 0
      n_inputs = 0
      outsuffix = ' '
      seen_outsuffix = .false.
      stages = ' '
      seen_stages = .false.
      order = 'convolve_reproject'
      footprint_mode = ' '
      seen_footprint_mode = .false.
      reffile = ' '
      seen_reffile = .false.
      have_badchan_file = .false.
      badchan_file = ' '
      have_target = .false.
      target_bmaj = 0.0d0
      target_bmin = 0.0d0
      target_bpa = 0.0d0
      have_max_common_bmaj = .false.
      max_common_bmaj = 0.0d0
      mem_frac_ram = 0.25
      npts = 2000
      khachiyan_tol = 1.0d-5
      raw_infiles = ' '
      raw_beamfiles = ' '
      have_cfgfile = .false.
      seen_infiles = .false.

      argc = command_argument_count()
      if (argc.eq.0) then
         call print_usage()
         status = -1
         return
      endif

      iarg = 1
      do while (iarg.le.argc)
         call get_command_argument(iarg, this_arg)
         if (trim(this_arg).eq.'--help' .or. trim(this_arg).eq.'-h') then
            call print_usage()
            status = -1
            return
         else if (trim(this_arg).eq.'--config') then
            if (iarg.eq.argc) then
               write(*,*) 'ERROR: --config requires a file path argument'
               status = -1
               return
            endif
            call get_command_argument(iarg+1, cfgfile)
            have_cfgfile = .true.
            iarg = iarg + 2
         else
            call split_cli_kv(this_arg, cli_key, cli_val, has_kv)
            if (.not. has_kv) then
               write(*,*) 'ERROR: unrecognised argument "', trim(this_arg),&
               &'" -- expected key=value, --config <file>, or --help'
               status = -1
               return
            endif
            call apply_kv(trim(cli_key), trim(cli_val), raw_infiles,&
            &raw_beamfiles, seen_infiles, seen_stages, status)
            if (status.ne.0) return
            iarg = iarg + 1
         endif
      enddo

      if (have_cfgfile) call read_cfg_file(cfgfile, raw_infiles, raw_beamfiles,&
      &seen_infiles, seen_stages, status)
      if (status.ne.0) return

      if (.not. seen_infiles .or. .not. seen_stages) then
         call print_usage()
         status = -1
         return
      endif

      if (trim(stages).ne.'reproject' .and. trim(stages).ne.'convolve'&
      &.and. trim(stages).ne.'both') then
         write(*,*) 'ERROR: stages must be reproject, convolve, or both'
         status = -1
         return
      endif
      if (trim(order).ne.'convolve_reproject' .and. trim(order).ne.'reproject_convolve') then
         write(*,*) 'ERROR: order must be convolve_reproject or reproject_convolve'
         status = -1
         return
      endif
      if (trim(order).ne.'convolve_reproject' .and. trim(stages).ne.'both') then
         write(*,*) 'NOTE: order= only matters when stages=both -- ignored'
      endif

      if (.not. seen_outsuffix) then
         if (trim(stages).eq.'reproject') then
            outsuffix = '_REPROJ.FITS'
         else if (trim(stages).eq.'convolve') then
            outsuffix = '_CONV.FITS'
         else
            outsuffix = '_MATCHED.FITS'
         endif
      endif

      n_inputs = cfg_csv_count(raw_infiles)
      if (n_inputs.lt.1 .or. n_inputs.gt.max_inputs) then
         write(*,*) 'ERROR: infiles must list between 1 and ', max_inputs, ' files'
         status = -1
         return
      endif
      do i = 1, n_inputs
         call cfg_csv_get_item(raw_infiles, i, infiles(i))
         beamfiles(i) = 'auto'
      enddo
      if (len_trim(raw_beamfiles).gt.0) then
         if (cfg_csv_count(raw_beamfiles).ne.n_inputs) then
            write(*,*) 'ERROR: beamfiles must list exactly ', n_inputs,&
            &' entries (one per infile; use ''auto'' for a file''s own',&
            &' BEAMS table)'
            status = -1
            return
         endif
         do i = 1, n_inputs
            call cfg_csv_get_item(raw_beamfiles, i, beamfiles(i))
         enddo
      endif

      if ((trim(stages).eq.'reproject' .or. trim(stages).eq.'both')) then
         if (.not. seen_footprint_mode) then
            write(*,*) 'ERROR: footprint_mode is required when stages includes reproject'
            status = -1
            return
         endif
         if (.not. seen_reffile) then
            write(*,*) 'ERROR: reffile is required when stages includes reproject'
            status = -1
            return
         endif
         if (trim(footprint_mode).ne.'intersection' .and.&
         &trim(footprint_mode).ne.'union' .and. trim(footprint_mode).ne.'reference') then
            write(*,*) 'ERROR: footprint_mode must be intersection, union, or reference'
            status = -1
            return
         endif
      endif

      if (mem_frac_ram.le.0.0 .or. mem_frac_ram.gt.0.95) then
         write(*,*) 'ERROR: mem_frac_ram must be > 0 and <= 0.95, got ', mem_frac_ram
         status = -1
         return
      endif
      if (npts.lt.12) then
         write(*,*) 'ERROR: npts must be at least 12, got ', npts
         status = -1
         return
      endif
   end subroutine parse_args

   subroutine apply_kv(key, val, raw_infiles, raw_beamfiles, seen_infiles,&
   &seen_stages, status)
      character(len=*), intent(in) :: key, val
      character(len=*), intent(inout) :: raw_infiles, raw_beamfiles
      logical, intent(inout) :: seen_infiles, seen_stages
      integer, intent(out) :: status
      integer :: ios

      status = 0
      select case (key)
      case ('stages')
         stages = val
         seen_stages = .true.
      case ('order')
         order = val
      case ('infiles')
         raw_infiles = val
         seen_infiles = .true.
      case ('beamfiles')
         raw_beamfiles = val
      case ('outsuffix')
         outsuffix = val
         seen_outsuffix = .true.
      case ('footprint_mode')
         footprint_mode = val
         seen_footprint_mode = .true.
      case ('reffile')
         reffile = val
         seen_reffile = .true.
      case ('badchan_file')
         badchan_file = val
         have_badchan_file = .true.
      case ('target_bmaj')
         read(val, *, iostat=ios) target_bmaj
         if (ios.ne.0) then
            write(*,*) 'ERROR: target_bmaj must be a number'
            status = -1
            return
         endif
         have_target = .true.
      case ('target_bmin')
         read(val, *, iostat=ios) target_bmin
         if (ios.ne.0) then
            write(*,*) 'ERROR: target_bmin must be a number'
            status = -1
            return
         endif
         have_target = .true.
      case ('target_bpa')
         read(val, *, iostat=ios) target_bpa
         if (ios.ne.0) then
            write(*,*) 'ERROR: target_bpa must be a number'
            status = -1
            return
         endif
         have_target = .true.
      case ('max_common_bmaj')
         read(val, *, iostat=ios) max_common_bmaj
         if (ios.ne.0) then
            write(*,*) 'ERROR: max_common_bmaj must be a number'
            status = -1
            return
         endif
         have_max_common_bmaj = .true.
      case ('mem_frac_ram')
         read(val, *, iostat=ios) mem_frac_ram
         if (ios.ne.0) then
            write(*,*) 'ERROR: mem_frac_ram must be a number'
            status = -1
            return
         endif
      case ('npts')
         read(val, *, iostat=ios) npts
         if (ios.ne.0) then
            write(*,*) 'ERROR: npts must be an integer'
            status = -1
            return
         endif
      case ('khachiyan_tol')
         read(val, *, iostat=ios) khachiyan_tol
         if (ios.ne.0) then
            write(*,*) 'ERROR: khachiyan_tol must be a number'
            status = -1
            return
         endif
      case default
         write(*,*) 'ERROR: unrecognised key "', key, '"'
         status = -1
         return
      end select
   end subroutine apply_kv

   subroutine read_cfg_file(cfgfile, raw_infiles, raw_beamfiles, seen_infiles,&
   &seen_stages, status)
      character(len=*), intent(in) :: cfgfile
      character(len=*), intent(inout) :: raw_infiles, raw_beamfiles
      logical, intent(inout) :: seen_infiles, seen_stages
      integer, intent(out) :: status
      character(len=512) :: line, key, val
      integer :: unit_cfg, ios, line_no
      logical :: has_kv

      status = 0
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
         call apply_kv(trim(key), trim(val), raw_infiles, raw_beamfiles,&
         &seen_infiles, seen_stages, status)
         if (status.ne.0) then
            write(*,*) '  (at line ', line_no, ' in ', trim(cfgfile), ')'
            close(unit_cfg)
            return
         endif
      enddo
      close(unit_cfg)
   end subroutine read_cfg_file

   subroutine print_usage()
      write(*,'(A)') 'match_cubes -- reproject and/or convolve FITS cubes,'//&
      &' chained through memory when both are requested'
      write(*,'(A)') ''
      write(*,'(A)') 'Usage:'
      write(*,'(A)') '  match_cubes stages=reproject|convolve|both'
      write(*,'(A)') '    [order=convolve_reproject|reproject_convolve]'
      write(*,'(A)') '    infiles=<file1>[,<file2>...]'
      write(*,'(A)') '    [footprint_mode=intersection|union|reference] [reffile=<file>]'
      write(*,'(A)') '    [beamfiles=<spec1>[,<spec2>...]] [badchan_file=<file>]'
      write(*,'(A)') '    [target_bmaj=<arcsec> target_bmin=<arcsec> target_bpa=<deg>]'
      write(*,'(A)') '    [max_common_bmaj=<arcsec>] [mem_frac_ram=<fraction>]'
      write(*,'(A)') '    [outsuffix=<suffix>] [npts=<n>] [khachiyan_tol=<tol>]'
      write(*,'(A)') '  match_cubes --config <cfgfile>'
      write(*,'(A)') '  match_cubes --help | -h'
      write(*,'(A)') ''
      write(*,'(A)') 'stages: reproject (align sky grids only, fully general axis'//&
      &' handling -- identical scope to the standalone reproject_cubes tool),'
      write(*,'(A)') '  convolve (common-resolution convolution only, identical scope'//&
      &' to the standalone convolve_cubes tool), or both -- chained through'
      write(*,'(A)') '  memory, no intermediate FITS file written. stages=convolve/both'//&
      &' require exactly 2 sky axes + 1 FREQ axis, every other axis degenerate'
      write(*,'(A)') '  (convolve_cubes''/gaussft_mod''s own existing scope).'
      write(*,'(A)') ''
      write(*,'(A)') 'order (default convolve_reproject, only meaningful for'//&
      &' stages=both): convolving before resampling avoids interpolation error'
      write(*,'(A)') '  on data that may only marginally sample its own native beam,'//&
      &' and is usually cheaper too (FFT work on the smaller native grid rather'
      write(*,'(A)') '  than a possibly-larger reprojected one). reproject_convolve is'//&
      &' not wrong, just not the default.'
      write(*,'(A)') ''
      write(*,'(A)') 'footprint_mode/reffile: required when stages includes reproject,'//&
      &' same semantics as the standalone reproject_cubes tool -- see its own --help.'
      write(*,'(A)') ''
      write(*,'(A)') 'beamfiles/badchan_file/target_bmaj/target_bmin/target_bpa/'//&
      &'max_common_bmaj/npts/khachiyan_tol: used when stages includes convolve,'
      write(*,'(A)') '  same semantics as the standalone convolve_cubes tool -- see its'//&
      &' own --help, and cfg/example_beamLog.txt/.csv for the ASCII beam format.'
      write(*,'(A)') ''
      write(*,'(A)') 'mem_frac_ram (default 0.25): fraction of total system RAM'//&
      &' budgeted for one read/process/write block of planes at a time.'
      write(*,'(A)') ''
      write(*,'(A)') 'outsuffix: appended to each infile''s own path for its output'//&
      &' filename. Default depends on stages: _REPROJ.FITS, _CONV.FITS, or'
      write(*,'(A)') '  _MATCHED.FITS for stages=both.'
   end subroutine print_usage

   subroutine cfg_split_key_value(raw_line, key, val, has_kv)
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
      character(len=*), intent(in) :: str
      integer :: n, ii

      n = 0
      if (len_trim(str) == 0) return
      n = 1
      do ii = 1, len_trim(str)
         if (str(ii:ii) == ',') n = n + 1
      enddo
   end function cfg_csv_count

   subroutine cfg_csv_get_item(str, idx, item)
      character(len=*), intent(in) :: str
      integer, intent(in) :: idx
      character(len=*), intent(out) :: item
      integer :: ii, cur, p0, n

      item = ' '
      n = len_trim(str)
      if (n == 0) return
      cur = 1
      p0 = 1
      do ii = 1, n
         if (str(ii:ii) == ',') then
            if (cur == idx) then
               item = adjustl(str(p0:ii - 1))
               return
            endif
            cur = cur + 1
            p0 = ii + 1
         endif
      enddo
      if (cur == idx) item = adjustl(str(p0:n))
   end subroutine cfg_csv_get_item

   !===========================================================
   ! Convolve-stage helpers (adapted verbatim from convolve_cubes.f90)
   !===========================================================

   subroutine read_badchan_file(filename, list, n, status)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: list(:)
      integer, intent(out) :: n
      integer, intent(out) :: status
      integer :: unit_bc, ios

      status = 0
      n = 0
      open(newunit=unit_bc, file=trim(filename), status='old', iostat=ios)
      if (ios.ne.0) then
         write(*,*) 'ERROR: cannot open badchan_file: ', trim(filename)
         status = -1
         return
      endif
      do
         if (n.ge.size(list)) then
            write(*,*) 'ERROR: too many entries in badchan_file (max ', size(list), ')'
            status = -1
            close(unit_bc)
            return
         endif
         n = n + 1
         read(unit_bc, *, iostat=ios) list(n)
         if (ios.ne.0) then
            n = n - 1
            exit
         endif
      enddo
      close(unit_bc)
   end subroutine read_badchan_file

   subroutine apply_badchan_list(list, n, nfreq, isbad)
      integer, intent(in) :: list(:), n, nfreq
      logical, intent(inout) :: isbad(nfreq)
      integer :: ii

      do ii = 1, n
         if (list(ii).ge.1 .and. list(ii).le.nfreq) isbad(list(ii)) = .true.
      enddo
   end subroutine apply_badchan_list

   subroutine pool_good_beams(n_inputs_l, nfreq, bmaj, bmin, bpa, isbad,&
   &dim1, dim2, n_pool_l, pool_bmaj_l, pool_bmin_l, pool_bpa_l)
      integer, intent(in) :: n_inputs_l, dim1, dim2, nfreq(dim1)
      real(dp), intent(in) :: bmaj(dim1,dim2), bmin(dim1,dim2), bpa(dim1,dim2)
      logical, intent(in) :: isbad(dim1,dim2)
      integer, intent(in) :: n_pool_l
      real(dp), intent(out) :: pool_bmaj_l(n_pool_l), pool_bmin_l(n_pool_l), pool_bpa_l(n_pool_l)
      integer :: ii, jj, k

      k = 0
      do ii = 1, n_inputs_l
         do jj = 1, nfreq(ii)
            if (.not. isbad(ii,jj)) then
               k = k + 1
               pool_bmaj_l(k) = bmaj(ii,jj)
               pool_bmin_l(k) = bmin(ii,jj)
               pool_bpa_l(k) = bpa(ii,jj)
            endif
         enddo
      enddo
   end subroutine pool_good_beams

   subroutine find_common_beam_wrap(n, bmaj, bmin, bpa, npts_in, tol_in,&
   &out_bmaj, out_bmin, out_bpa, status)
      use commonbeam_mod, only: find_common_beam
      integer, intent(in) :: n, npts_in
      real(dp), intent(in) :: bmaj(n), bmin(n), bpa(n), tol_in
      real(dp), intent(out) :: out_bmaj, out_bmin, out_bpa
      integer, intent(out) :: status

      call find_common_beam(n, bmaj, bmin, bpa, npts_in, tol_in,&
      &out_bmaj, out_bmin, out_bpa, status)
   end subroutine find_common_beam_wrap

   subroutine read_axis_info(filename, naxis, sky1, sky2, freq_axis, naxes,&
   &cdelt1, cdelt2, status)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: naxis, sky1, sky2, freq_axis
      integer, intent(out) :: naxes(max_axes)
      real(dp), intent(out) :: cdelt1, cdelt2
      integer, intent(out) :: status

      integer :: unit, blocksize, fitsstat, k
      character(len=68) :: ctype, comment
      character(len=8) :: axstr

      status = 0
      fitsstat = 0
      unit = 200
      call FTOPEN(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(filename)
         status = -1
         return
      endif

      call FTGKYJ(unit, 'NAXIS', naxis, comment, fitsstat)
      if (fitsstat.ne.0 .or. naxis.lt.2 .or. naxis.gt.max_axes) then
         write(*,*) 'ERROR: bad or missing NAXIS in: ', trim(filename)
         status = -1
         call FTCLOS(unit, fitsstat)
         return
      endif

      sky1 = 0
      sky2 = 0
      freq_axis = 0
      naxes = 0
      do k = 1, naxis
         write(axstr,'(I0)') k
         fitsstat = 0
         call FTGKYJ(unit, 'NAXIS'//trim(axstr), naxes(k), comment, fitsstat)
         if (fitsstat.ne.0) then
            write(*,*) 'ERROR: missing NAXIS', k, ' in: ', trim(filename)
            status = -1
            call FTCLOS(unit, fitsstat)
            return
         endif
         fitsstat = 0
         call FTGKYS(unit, 'CTYPE'//trim(axstr), ctype, comment, fitsstat)
         if (fitsstat.ne.0) cycle
         ctype = adjustl(ctype)
         if (ctype(1:2).eq.'RA') then
            if (sky1.eq.0) then
               sky1 = k
            else
               write(*,*) 'ERROR: more than one RA-like axis in: ', trim(filename)
               status = -1
               call FTCLOS(unit, fitsstat)
               return
            endif
         else if (ctype(1:3).eq.'DEC') then
            if (sky2.eq.0) then
               sky2 = k
            else
               write(*,*) 'ERROR: more than one DEC-like axis in: ', trim(filename)
               status = -1
               call FTCLOS(unit, fitsstat)
               return
            endif
         else if (ctype(1:4).eq.'FREQ') then
            freq_axis = k
         endif
      enddo

      if (sky1.eq.0 .or. sky2.eq.0) then
         write(*,*) 'ERROR: could not identify RA/DEC sky axes in: ', trim(filename)
         status = -1
         call FTCLOS(unit, fitsstat)
         return
      endif
      if (freq_axis.eq.0) then
         write(*,*) 'ERROR: could not identify a FREQ axis in: ', trim(filename)
         status = -1
         call FTCLOS(unit, fitsstat)
         return
      endif

      do k = 1, naxis
         if (k.ne.sky1 .and. k.ne.sky2 .and. k.ne.freq_axis) then
            if (naxes(k).gt.1) then
               write(*,*) 'ERROR: ', trim(filename), ' axis ', k,&
               &' has extent ', naxes(k), ' > 1 -- only 2 sky axes plus one'//&
               &' FREQ axis are supported; run separate slices (e.g. per'//&
               &' Stokes) as separate infiles'
               status = -1
               call FTCLOS(unit, fitsstat)
               return
            endif
         endif
      enddo

      write(axstr,'(I0)') sky1
      fitsstat = 0
      call FTGKYD(unit, 'CDELT'//trim(axstr), cdelt1, comment, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: missing CDELT for the RA axis in: ', trim(filename)
         status = -1
         call FTCLOS(unit, fitsstat)
         return
      endif
      write(axstr,'(I0)') sky2
      fitsstat = 0
      call FTGKYD(unit, 'CDELT'//trim(axstr), cdelt2, comment, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: missing CDELT for the DEC axis in: ', trim(filename)
         status = -1
         call FTCLOS(unit, fitsstat)
         return
      endif

      call check_no_rotation(unit, sky1, sky2, filename, status)
      if (status.ne.0) then
         call FTCLOS(unit, fitsstat)
         return
      endif

      call FTCLOS(unit, fitsstat)
   end subroutine read_axis_info

   subroutine check_no_rotation(unit, sky1, sky2, filename, status)
      integer, intent(in) :: unit, sky1, sky2
      character(len=*), intent(in) :: filename
      integer, intent(out) :: status
      integer :: fitsstat
      character(len=68) :: comment
      character(len=8) :: a1, a2
      real(dp) :: dval

      status = 0
      write(a1,'(I0)') sky1
      fitsstat = 0
      call FTGKYD(unit, 'CROTA'//trim(a1), dval, comment, fitsstat)
      if (fitsstat.eq.0 .and. dval.ne.0.0d0) then
         write(*,*) 'ERROR: ', trim(filename), ' has a nonzero CROTA',&
         &trim(a1), ' -- rotated sky grids are not supported by this'//&
         &' program''s sky-to-pixel BPA conversion; reproject onto an'//&
         &' axis-aligned grid first (see reproject_cubes)'
         status = -1
         return
      endif
      write(a2,'(I0)') sky2
      fitsstat = 0
      call FTGKYD(unit, 'CROTA'//trim(a2), dval, comment, fitsstat)
      if (fitsstat.eq.0 .and. dval.ne.0.0d0) then
         write(*,*) 'ERROR: ', trim(filename), ' has a nonzero CROTA',&
         &trim(a2), ' -- rotated sky grids are not supported'
         status = -1
         return
      endif

      call check_offdiag(unit, 'PC', sky1, sky2, filename, status)
      if (status.ne.0) return
      call check_offdiag(unit, 'CD', sky1, sky2, filename, status)
   end subroutine check_no_rotation

   subroutine check_offdiag(unit, prefix, sky1, sky2, filename, status)
      integer, intent(in) :: unit, sky1, sky2
      character(len=*), intent(in) :: prefix, filename
      integer, intent(out) :: status
      integer :: fitsstat
      character(len=68) :: comment
      character(len=16) :: k12, k21
      character(len=8) :: a1, a2
      real(dp) :: dval

      status = 0
      write(a1,'(I0)') sky1
      write(a2,'(I0)') sky2
      k12 = trim(prefix)//trim(a1)//'_'//trim(a2)
      fitsstat = 0
      call FTGKYD(unit, trim(k12), dval, comment, fitsstat)
      if (fitsstat.eq.0 .and. dval.ne.0.0d0) then
         write(*,*) 'ERROR: ', trim(filename), ' has a nonzero ', trim(k12),&
         &' -- rotated/sheared sky grids are not supported'
         status = -1
         return
      endif
      k21 = trim(prefix)//trim(a2)//'_'//trim(a1)
      fitsstat = 0
      call FTGKYD(unit, trim(k21), dval, comment, fitsstat)
      if (fitsstat.eq.0 .and. dval.ne.0.0d0) then
         write(*,*) 'ERROR: ', trim(filename), ' has a nonzero ', trim(k21),&
         &' -- rotated/sheared sky grids are not supported'
         status = -1
         return
      endif
   end subroutine check_offdiag

   subroutine read_beams(filename, beamspec, nfreq, bmaj, bmin, bpa, isbad, status)
      character(len=*), intent(in) :: filename, beamspec
      integer, intent(in) :: nfreq
      real(dp), intent(out) :: bmaj(nfreq), bmin(nfreq), bpa(nfreq)
      logical, intent(out) :: isbad(nfreq)
      integer, intent(out) :: status

      if (trim(beamspec).eq.'auto' .or. len_trim(beamspec).eq.0) then
         call read_beams_table(filename, nfreq, bmaj, bmin, bpa, isbad, status)
      else
         call read_beams_ascii(beamspec, nfreq, bmaj, bmin, bpa, isbad, status)
      endif
   end subroutine read_beams

   subroutine read_beams_table(filename, nfreq, bmaj, bmin, bpa, isbad, status)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nfreq
      real(dp), intent(out) :: bmaj(nfreq), bmin(nfreq), bpa(nfreq)
      logical, intent(out) :: isbad(nfreq)
      integer, intent(out) :: status

      integer :: unit, blocksize, fitsstat, nrows, k
      integer :: col_bmaj, col_bmin, col_bpa, col_chan
      real, allocatable :: rb_bmaj(:), rb_bmin(:), rb_bpa(:)
      integer, allocatable :: rb_chan(:)
      logical :: anyflag
      integer :: chan1

      status = 0
      fitsstat = 0
      unit = 201
      call FTOPEN(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(filename)
         status = -1
         return
      endif

      fitsstat = 0
      call FTMNHD(unit, -1, 'BEAMS', 0, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: no BEAMS binary table extension found in: ',&
         &trim(filename), ' -- pass an ASCII beamfile instead (see --help)'
         status = -1
         call FTCLOS(unit, fitsstat)
         return
      endif

      fitsstat = 0
      call FTGNRW(unit, nrows, fitsstat)
      if (fitsstat.ne.0 .or. nrows.lt.1) then
         write(*,*) 'ERROR: could not read row count of BEAMS table in: ', trim(filename)
         status = -1
         call FTCLOS(unit, fitsstat)
         return
      endif

      fitsstat = 0
      call FTGCNO(unit, .false., 'BMAJ', col_bmaj, fitsstat)
      call FTGCNO(unit, .false., 'BMIN', col_bmin, fitsstat)
      call FTGCNO(unit, .false., 'BPA', col_bpa, fitsstat)
      call FTGCNO(unit, .false., 'CHAN', col_chan, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: BEAMS table in ', trim(filename),&
         &' missing one of BMAJ/BMIN/BPA/CHAN columns'
         status = -1
         call FTCLOS(unit, fitsstat)
         return
      endif

      allocate(rb_bmaj(nrows), rb_bmin(nrows), rb_bpa(nrows), rb_chan(nrows))
      fitsstat = 0
      call FTGCVE(unit, col_bmaj, 1, 1, nrows, 0.0, rb_bmaj, anyflag, fitsstat)
      call FTGCVE(unit, col_bmin, 1, 1, nrows, 0.0, rb_bmin, anyflag, fitsstat)
      call FTGCVE(unit, col_bpa, 1, 1, nrows, 0.0, rb_bpa, anyflag, fitsstat)
      call FTGCVJ(unit, col_chan, 1, 1, nrows, 0, rb_chan, anyflag, fitsstat)
      call FTCLOS(unit, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed reading BEAMS table columns in: ', trim(filename)
         status = -1
         deallocate(rb_bmaj, rb_bmin, rb_bpa, rb_chan)
         return
      endif

      isbad = .true.
      bmaj = 0.0d0
      bmin = 0.0d0
      bpa = 0.0d0
      do k = 1, nrows
         chan1 = rb_chan(k) + 1
         if (chan1.lt.1 .or. chan1.gt.nfreq) cycle
         bmaj(chan1) = real(rb_bmaj(k), kind=8)
         bmin(chan1) = real(rb_bmin(k), kind=8)
         bpa(chan1) = real(rb_bpa(k), kind=8)
         isbad(chan1) = (bmaj(chan1).lt.1.0d-6 .or. bmin(chan1).lt.1.0d-6)
      enddo
      deallocate(rb_bmaj, rb_bmin, rb_bpa, rb_chan)
   end subroutine read_beams_table

   subroutine read_beams_ascii(beamfile, nfreq, bmaj, bmin, bpa, isbad, status)
      character(len=*), intent(in) :: beamfile
      integer, intent(in) :: nfreq
      real(dp), intent(out) :: bmaj(nfreq), bmin(nfreq), bpa(nfreq)
      logical, intent(out) :: isbad(nfreq)
      integer, intent(out) :: status

      integer :: unit_bf, ios
      character(len=512) :: line
      integer :: ich
      real(dp) :: a, b, p

      status = 0
      isbad = .true.
      bmaj = 0.0d0
      bmin = 0.0d0
      bpa = 0.0d0

      open(newunit=unit_bf, file=trim(beamfile), status='old', iostat=ios)
      if (ios.ne.0) then
         write(*,*) 'ERROR: cannot open ASCII beam file: ', trim(beamfile)
         status = -1
         return
      endif
      do
         read(unit_bf, '(A)', iostat=ios) line
         if (ios.ne.0) exit
         line = adjustl(line)
         if (len_trim(line).eq.0) cycle
         if (line(1:1).eq.'#') cycle
         read(line, *, iostat=ios) ich, a, b, p
         if (ios.ne.0) then
            write(*,*) 'ERROR: malformed line in ', trim(beamfile), ': ', trim(line)
            status = -1
            close(unit_bf)
            return
         endif
         if (ich.lt.1 .or. ich.gt.nfreq) then
            write(*,*) 'ERROR: channel ', ich, ' in ', trim(beamfile),&
            &' is out of range 1..', nfreq
            status = -1
            close(unit_bf)
            return
         endif
         bmaj(ich) = a
         bmin(ich) = b
         bpa(ich) = p
         isbad(ich) = (a.lt.1.0d-6 .or. b.lt.1.0d-6)
      enddo
      close(unit_bf)
   end subroutine read_beams_ascii

   subroutine sky_to_pixel_bpa(bpa_sky_deg, cdelt1, cdelt2, bpa_pixel_deg)
      real(dp), intent(in) :: bpa_sky_deg, cdelt1, cdelt2
      real(dp), intent(out) :: bpa_pixel_deg
      real(dp), parameter :: pi = 3.14159265358979323846d0
      real(dp) :: theta, s1, s2

      theta = bpa_sky_deg*pi/180.0d0
      s1 = sign(1.0d0, cdelt1)
      s2 = sign(1.0d0, cdelt2)
      bpa_pixel_deg = atan2(s2*cos(theta), s1*sin(theta))*180.0d0/pi
   end subroutine sky_to_pixel_bpa

   !===========================================================
   ! Reproject-stage AST helpers (adapted verbatim from reproject_cubes.f90)
   !===========================================================

   subroutine compose_pix2pix(skymap_from, skyframe_from, skymap_to,&
   &skyframe_to, map_out, status)
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

   !===========================================================
   ! Header propagation (adapted from reproject_cubes.f90; the BMAJ/BMIN/
   ! BPA/CASAMBM exclusion additionally applied whenever convolve is
   ! active, matching convolve_cubes.f90's own exclusion list)
   !===========================================================

   subroutine copy_axis_keywords(src_unit, src_axis, dst_unit, dst_axis,&
   &crpix_shift, status)
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

      write(axstr,'(I0)') src_axis
      fitsstat = 0
      call FTGKYD(src_unit, 'CROTA'//trim(axstr), dval, comment, fitsstat)
      if (fitsstat.eq.0) then
         write(axstr,'(I0)') dst_axis
         call FTPKYD(dst_unit, 'CROTA'//trim(axstr), dval, 13, ' ', fitsstat)
      endif
   end subroutine copy_axis_keywords

   subroutine copy_sky_rotation_matrix(src_unit, src_axis1, src_axis2,&
   &dst_unit, status)
      integer, intent(in) :: src_unit, src_axis1, src_axis2, dst_unit
      integer, intent(inout) :: status

      if (status.ne.0) return

      call copy_one_matrix_entry(src_unit, 'PC', src_axis1, src_axis1,&
      &dst_unit, 1, 1)
      call copy_one_matrix_entry(src_unit, 'PC', src_axis1, src_axis2,&
      &dst_unit, 1, 2)
      call copy_one_matrix_entry(src_unit, 'PC', src_axis2, src_axis1,&
      &dst_unit, 2, 1)
      call copy_one_matrix_entry(src_unit, 'PC', src_axis2, src_axis2,&
      &dst_unit, 2, 2)
      call copy_one_matrix_entry(src_unit, 'CD', src_axis1, src_axis1,&
      &dst_unit, 1, 1)
      call copy_one_matrix_entry(src_unit, 'CD', src_axis1, src_axis2,&
      &dst_unit, 1, 2)
      call copy_one_matrix_entry(src_unit, 'CD', src_axis2, src_axis1,&
      &dst_unit, 2, 1)
      call copy_one_matrix_entry(src_unit, 'CD', src_axis2, src_axis2,&
      &dst_unit, 2, 2)
   end subroutine copy_sky_rotation_matrix

   subroutine copy_one_matrix_entry(su, prefix, sa, sb, du, da, db)
      integer, intent(in) :: su, sa, sb, du, da, db
      character(len=*), intent(in) :: prefix

      character(len=16) :: srckey, dstkey
      character(len=68) :: comment
      double precision :: dval
      integer :: fitsstat
      character(len=4) :: si, sj

      write(si,'(I0)') sa
      write(sj,'(I0)') sb
      srckey = trim(prefix)//trim(si)//'_'//trim(sj)
      fitsstat = 0
      call FTGKYD(su, trim(srckey), dval, comment, fitsstat)
      if (fitsstat.eq.0) then
         write(si,'(I0)') da
         write(sj,'(I0)') db
         dstkey = trim(prefix)//trim(si)//'_'//trim(sj)
         call FTPKYD(du, trim(dstkey), dval, 13, ' ', fitsstat)
      endif
   end subroutine copy_one_matrix_entry

   subroutine copy_wcs_system_keywords(src_unit, dst_unit, status)
      integer, intent(in) :: src_unit, dst_unit
      integer, intent(inout) :: status

      integer :: fitsstat
      character(len=68) :: comment, sval
      double precision :: dval

      if (status.ne.0) return

      fitsstat = 0
      call FTGKYD(src_unit, 'EQUINOX', dval, comment, fitsstat)
      if (fitsstat.eq.0) call FTPKYD(dst_unit, 'EQUINOX', dval, 13, ' ', fitsstat)

      fitsstat = 0
      call FTGKYS(src_unit, 'RADESYS', sval, comment, fitsstat)
      if (fitsstat.eq.0) call FTPKYS(dst_unit, 'RADESYS', trim(sval), ' ', fitsstat)
   end subroutine copy_wcs_system_keywords

   subroutine copy_generic_header_match(src_unit, dst_unit, exclude_axis_wcs,&
   &exclude_beam, status)
      !! Generalization of reproject_cubes' copy_generic_header (excludes
      !! structural keywords + axis-indexed WCS, keeping BMAJ/BMIN/BPA/
      !! CASAMBM passthrough) and convolve_cubes' copy_generic_header_convolve
      !! (excludes structural + BMAJ/BMIN/BPA/CASAMBM, keeps axis WCS
      !! passthrough): here both exclusions are independently selectable,
      !! since match_cubes can be in either situation depending on which
      !! stage(s) are active for a given output file. exclude_axis_wcs
      !! should be true whenever reproject is active (its own explicit
      !! copy_axis_keywords/copy_sky_rotation_matrix/
      !! copy_wcs_system_keywords calls handle those keywords instead);
      !! exclude_beam should be true whenever convolve is active (its own
      !! explicit BMAJ/BMIN/BPA/CASAMBM writing handles those instead).
      integer, intent(in) :: src_unit, dst_unit
      logical, intent(in) :: exclude_axis_wcs, exclude_beam
      integer, intent(inout) :: status

      integer :: nkeys, nmore, i, fitsstat
      character(len=80) :: card
      character(len=8) :: key

      if (status.ne.0) return

      fitsstat = 0
      call FTGHSP(src_unit, nkeys, nmore, fitsstat)
      do i = 1, nkeys
         fitsstat = 0
         call FTGREC(src_unit, i, card, fitsstat)
         if (fitsstat.ne.0) cycle
         key = adjustl(card(1:8))
         if (skip_generic_header_key(key, exclude_axis_wcs, exclude_beam)) cycle
         fitsstat = 0
         call FTPREC(dst_unit, card, fitsstat)
      enddo
   end subroutine copy_generic_header_match

   logical function skip_generic_header_key(key, exclude_axis_wcs, exclude_beam)
      character(len=8), intent(in) :: key
      logical, intent(in) :: exclude_axis_wcs, exclude_beam

      skip_generic_header_key = .true.
      select case (trim(key))
      case ('SIMPLE', 'BITPIX', 'NAXIS', 'EXTEND', 'PCOUNT', 'GCOUNT', 'END')
         return
      end select
      if (is_indexed_keyword(key, 'NAXIS')) return
      if (exclude_axis_wcs) then
         if (is_indexed_keyword(key, 'CTYPE')) return
         if (is_indexed_keyword(key, 'CRVAL')) return
         if (is_indexed_keyword(key, 'CRPIX')) return
         if (is_indexed_keyword(key, 'CDELT')) return
         if (is_indexed_keyword(key, 'CUNIT')) return
         if (is_indexed_keyword(key, 'CROTA')) return
         if (is_indexed_keyword(key, 'PC')) return
         if (is_indexed_keyword(key, 'CD')) return
         select case (trim(key))
         case ('EQUINOX', 'RADESYS', 'WCSAXES', 'LONPOLE', 'LATPOLE')
            return
         end select
      endif
      if (exclude_beam) then
         select case (trim(key))
         case ('BMAJ', 'BMIN', 'BPA', 'CASAMBM')
            return
         end select
      endif
      skip_generic_header_key = .false.
   end function skip_generic_header_key

   logical function is_indexed_keyword(key, prefix)
      character(len=8), intent(in) :: key
      character(len=*), intent(in) :: prefix

      character(len=8) :: rest
      integer :: plen, klen, i
      logical :: seen_underscore

      is_indexed_keyword = .false.
      plen = len_trim(prefix)
      klen = len_trim(key)
      if (klen.le.plen) return
      if (key(1:plen).ne.prefix(1:plen)) return

      rest = key(plen+1:klen)
      seen_underscore = .false.
      do i = 1, len_trim(rest)
         if (rest(i:i).eq.'_') then
            if (seen_underscore) return
            seen_underscore = .true.
         else if (rest(i:i).lt.'0' .or. rest(i:i).gt.'9') then
            return
         endif
      enddo
      if (len_trim(rest).eq.0) return
      if (rest(1:1).eq.'_' .or. rest(len_trim(rest):len_trim(rest)).eq.'_') return
      is_indexed_keyword = .true.
   end function is_indexed_keyword

   !===========================================================
   ! stages=reproject alone: fully general N-dimensional axis handling,
   ! adapted verbatim from reproject_cubes.f90's write_reprojected_file/
   ! read_one_block/write_one_block.
   !===========================================================

   subroutine process_one_file_general(reffile_l, infile, outfile, pixaxes_ref_l,&
   &naxes_in_l, pixaxes_in_l, lbnd_out_d, ubnd_out_d, mem_frac_ram_l, status)
      use, intrinsic :: ieee_arithmetic
      use omp_lib, only: omp_get_max_threads
      character(len=*), intent(in) :: reffile_l, infile, outfile
      integer, intent(in) :: pixaxes_ref_l(2)
      integer, intent(in) :: naxes_in_l(:), pixaxes_in_l(2)
      double precision, intent(in) :: lbnd_out_d(2), ubnd_out_d(2)
      real, intent(in) :: mem_frac_ram_l
      integer, intent(inout) :: status

      integer :: naxis, k, other_axes(max_axes), n_other
      integer :: other_idx(max_axes), remainder, radix
      integer :: n_planes, status_par, nthreads
      integer :: nx_out, ny_out, naxis_out, naxes_out(max_axes)
      integer :: nx_in, ny_in
      integer :: ref_unit, out_unit, fitsstat, blocksize
      logical :: simple, extend
      integer :: beams_unit, beams_status, casambm_status, hdutype_dum
      logical :: casambm_val
      character(len=80) :: comment_dum

      integer(kind=8) :: mem_total_kb, bytes_per_plane, mem_safe_bytes
      integer(kind=8) :: block_planes64
      integer :: block_planes, n_groups, igroup, axis1_extent
      integer :: chan_start, chan_len, local_iplane
      real, allocatable :: block_data_in(:,:,:), block_data_out(:,:,:)

      integer :: t_status, t_wcs_ref, t_skymap_ref, t_skyframe_ref
      integer :: t_naxes_ref(max_axes), t_pixaxes_ref(2)
      integer :: t_wcs_in, t_skymap_in, t_skyframe_in
      integer :: t_naxes_in(max_axes), t_pixaxes_in(2)
      integer :: t_map_in2ref

      integer :: lbnd_in(2), ubnd_in(2), lbnd_o(2), ubnd_o(2), nbad
      real :: badval
      double precision :: params_dummy(1)

      if (status.ne.0) return

      naxis = 0
      do k = 1, size(naxes_in_l)
         if (naxes_in_l(k).gt.0) naxis = k
      enddo
      n_other = 0
      do k = 1, naxis
         if (k.ne.pixaxes_in_l(1) .and. k.ne.pixaxes_in_l(2)) then
            n_other = n_other + 1
            other_axes(n_other) = k
         endif
      enddo

      nx_in = naxes_in_l(pixaxes_in_l(1))
      ny_in = naxes_in_l(pixaxes_in_l(2))
      nx_out = nint(ubnd_out_d(1) - lbnd_out_d(1)) + 1
      ny_out = nint(ubnd_out_d(2) - lbnd_out_d(2)) + 1
      naxis_out = 2 + n_other
      naxes_out(1) = nx_out
      naxes_out(2) = ny_out
      do k = 1, n_other
         naxes_out(2+k) = naxes_in_l(other_axes(k))
      enddo

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

      ref_unit = 44
      fitsstat = 0
      call FTOPEN(ref_unit, trim(reffile_l), 0, blocksize, fitsstat)
      call copy_axis_keywords(ref_unit, pixaxes_ref_l(1), out_unit, 1,&
      &lbnd_out_d(1)-1.0d0, status)
      call copy_axis_keywords(ref_unit, pixaxes_ref_l(2), out_unit, 2,&
      &lbnd_out_d(2)-1.0d0, status)
      call copy_sky_rotation_matrix(ref_unit, pixaxes_ref_l(1),&
      &pixaxes_ref_l(2), out_unit, status)
      call copy_wcs_system_keywords(ref_unit, out_unit, status)
      call FTCLOS(ref_unit, fitsstat)

      fitsstat = 0
      call FTOPEN(ref_unit, trim(infile), 0, blocksize, fitsstat)
      do k = 1, n_other
         call copy_axis_keywords(ref_unit, other_axes(k), out_unit,&
         &2+k, 0.0d0, status)
      enddo
      call copy_generic_header_match(ref_unit, out_unit, .true., .false., status)
      call FTPHIS(out_unit, 'match_cubes: reprojected from '//&
      &trim(infile)//' onto the grid of '//trim(reffile_l), fitsstat)
      call FTCLOS(ref_unit, fitsstat)

      ! CASAMBM/BEAMS: reproject-alone (stages=reproject) never touches
      ! the beam itself -- see reproject_cubes.f90's own identical block,
      ! which this one is a verbatim port of. copy_generic_header_match
      ! above already copied the scalar CASAMBM keyword verbatim as a
      ! raw header card (only PRIMARY-header cards, though); this
      ! attaches the actual BEAMS extension HDU it refers to. Own
      ! dedicated unit (45), disjoint from every other unit number this
      ! file uses (43/44/45/200/201/210/211/212/220/1000+thread/5000/
      ! 5100).
      casambm_status = 0
      call ftgkyl(out_unit, 'CASAMBM', casambm_val, comment_dum, casambm_status)
      if (casambm_status.eq.0 .and. casambm_val) then
         beams_unit = 45
         beams_status = 0
         call ftopen(beams_unit, trim(infile), 0, blocksize, beams_status)
         call ftmnhd(beams_unit, -1, 'BEAMS', 0, beams_status)
         if (beams_status.eq.0) then
            status = 0
            call ftcopy(beams_unit, out_unit, 0, status)
            call ftmahd(out_unit, 1, hdutype_dum, status)
         else
            write(*,*) 'WARNING: CASAMBM=T but no BEAMS extension found in: ',&
            &trim(infile), ' -- output keeps CASAMBM=T with no BEAMS table.'
         endif
         beams_status = 0
         call ftclos(beams_unit, beams_status)
      endif

      n_planes = 1
      do k = 1, n_other
         n_planes = n_planes * naxes_in_l(other_axes(k))
      enddo
      call get_mem_total_kb(mem_total_kb)
      bytes_per_plane = int(4,8) * (int(nx_in,8)*int(ny_in,8) +&
      &int(nx_out,8)*int(ny_out,8))
      mem_safe_bytes = int(real(mem_frac_ram_l,8) * real(mem_total_kb,8) *&
      &1024.0d0, 8)
      block_planes64 = max(1_8, mem_safe_bytes / bytes_per_plane)
      block_planes = int(min(block_planes64, int(n_planes,8)))
      axis1_extent = 1
      if (n_other.ge.1) then
         axis1_extent = naxes_in_l(other_axes(1))
         block_planes = min(block_planes, axis1_extent)
      endif
      if (block_planes.lt.1) block_planes = 1

      write(*,'(A,A,A,I0,A,I0,A)') 'Writing ', trim(outfile), ': ',&
      &n_planes, ' plane(s), in blocks of up to ', block_planes, ' plane(s)'

      if (block_planes.lt.omp_get_max_threads()) then
         write(*,'(A,I0,A,I0,A)') 'WARNING: mem_frac_ram limits blocks to ',&
         &block_planes, ' plane(s), below the ', omp_get_max_threads(),&
         &' threads available -- parallelism (not just I/O) is reduced;'//&
         &' raise mem_frac_ram for full speedup if memory allows.'
      endif

      allocate(block_data_in(nx_in, ny_in, block_planes))
      allocate(block_data_out(nx_out, ny_out, block_planes))

      n_groups = 1
      do k = 2, n_other
         n_groups = n_groups * naxes_in_l(other_axes(k))
      enddo

      status_par = 0
      nthreads = max(1, min(omp_get_max_threads(), block_planes))
      !$omp parallel num_threads(nthreads) default(none)&
      !$omp& shared(infile, reffile_l, naxes_in_l, pixaxes_in_l, other_axes,&
      !$omp& n_other, lbnd_out_d, ubnd_out_d, out_unit, status_par,&
      !$omp& nx_in, ny_in, nx_out, ny_out, block_planes, block_data_in,&
      !$omp& block_data_out, n_groups, axis1_extent)&
      !$omp& private(t_status, t_wcs_ref, t_skymap_ref, t_skyframe_ref,&
      !$omp& t_naxes_ref, t_pixaxes_ref, t_wcs_in, t_skymap_in, t_skyframe_in,&
      !$omp& t_naxes_in, t_pixaxes_in, t_map_in2ref, other_idx, remainder,&
      !$omp& radix, k, igroup, chan_start, chan_len, local_iplane, nbad,&
      !$omp& lbnd_in, ubnd_in, lbnd_o, ubnd_o, badval, params_dummy)

      t_status = 0
      call ast_begin(t_status)
      call load_wcs(reffile_l, t_wcs_ref, t_naxes_ref, t_status)
      call extract_sky_mapping(t_wcs_ref, t_skymap_ref, t_skyframe_ref,&
      &t_pixaxes_ref, t_status)
      call load_wcs(infile, t_wcs_in, t_naxes_in, t_status)
      call extract_sky_mapping(t_wcs_in, t_skymap_in, t_skyframe_in,&
      &t_pixaxes_in, t_status)
      call compose_pix2pix(t_skymap_in, t_skyframe_in, t_skymap_ref,&
      &t_skyframe_ref, t_map_in2ref, t_status)
      if (t_status.ne.0) then
         !$omp atomic write
         status_par = -1
      endif

      do igroup = 1, n_groups
         remainder = igroup - 1
         do k = 2, n_other
            radix = naxes_in_l(other_axes(k))
            other_idx(k) = mod(remainder, radix) + 1
            remainder = remainder / radix
         enddo

         chan_start = 1
         do while (chan_start.le.axis1_extent)
            chan_len = min(block_planes, axis1_extent - chan_start + 1)

            !$omp single
            call read_one_block(infile, naxes_in_l, pixaxes_in_l, other_axes,&
            &other_idx, n_other, chan_start, chan_len, nx_in, ny_in,&
            &block_data_in(:,:,1:chan_len), status_par)
            !$omp end single

            !$omp do schedule(static)
            do local_iplane = 1, chan_len
               if (t_status.eq.0 .and. status_par.eq.0) then
                  lbnd_in(1) = 1
                  lbnd_in(2) = 1
                  ubnd_in(1) = nx_in
                  ubnd_in(2) = ny_in
                  lbnd_o(1) = nint(lbnd_out_d(1))
                  lbnd_o(2) = nint(lbnd_out_d(2))
                  ubnd_o(1) = nint(ubnd_out_d(1))
                  ubnd_o(2) = nint(ubnd_out_d(2))
                  badval = ieee_value(badval, ieee_quiet_nan)
                  params_dummy(1) = 0.0d0
                  nbad = ast_resampler(t_map_in2ref, 2, lbnd_in, ubnd_in,&
                  &block_data_in(:,:,local_iplane),&
                  &block_data_in(:,:,local_iplane),&
                  &ast__linear, ast_null, params_dummy, 0, 0.0d0, 100, badval,&
                  &2, lbnd_o, ubnd_o, lbnd_o, ubnd_o,&
                  &block_data_out(:,:,local_iplane),&
                  &block_data_out(:,:,local_iplane), t_status)
                  if (t_status.ne.0) then
                     !$omp atomic write
                     status_par = -1
                  endif
               endif
            enddo
            !$omp end do

            !$omp single
            call write_one_block(out_unit, naxes_in_l, other_axes, other_idx,&
            &n_other, chan_start, chan_len, nx_out, ny_out,&
            &block_data_out(:,:,1:chan_len), status_par)
            !$omp end single

            chan_start = chan_start + chan_len
         enddo
      enddo

      call ast_end(t_status)
      !$omp end parallel

      deallocate(block_data_in)
      deallocate(block_data_out)

      if (status_par.ne.0) then
         write(*,*) 'ERROR: failed to resample/write one or more planes for: ',&
         &trim(infile)
         status = -1
         call FTCLOS(out_unit, fitsstat)
         return
      endif

      call FTCLOS(out_unit, fitsstat)
   end subroutine process_one_file_general

   subroutine read_one_block(filename, naxes_in_l, pixaxes_in_l, other_axes,&
   &other_idx, n_other, chan_start, chan_len, nx_in, ny_in, block_data_in,&
   &status)
      use, intrinsic :: ieee_arithmetic
      character(len=*), intent(in) :: filename
      integer, intent(in) :: naxes_in_l(:), pixaxes_in_l(2)
      integer, intent(in) :: other_axes(:), other_idx(:), n_other
      integer, intent(in) :: chan_start, chan_len, nx_in, ny_in
      real, intent(out) :: block_data_in(:,:,:)
      integer, intent(inout) :: status

      integer :: unit, blocksize, fitsstat, group, naxis, k
      integer :: fpixels(max_axes), lpixels(max_axes), incs(max_axes)
      logical :: anyflg
      real :: badval
      integer :: ax_sky1, ax_sky2, ax_block
      integer :: rank_sky1, rank_sky2, rank_block
      integer :: dims(3), idxvec(3), i, j, c
      logical :: natural_order
      real, allocatable :: natural_buf(:,:,:)

      if (status.ne.0) return

      naxis = 0
      do k = 1, max_axes
         if (naxes_in_l(k).gt.0) naxis = k
      enddo
      fpixels(1:naxis) = 1
      lpixels(1:naxis) = 1
      incs(1:naxis) = 1
      lpixels(pixaxes_in_l(1)) = nx_in
      lpixels(pixaxes_in_l(2)) = ny_in
      if (n_other.ge.1) then
         fpixels(other_axes(1)) = chan_start
         lpixels(other_axes(1)) = chan_start + chan_len - 1
      endif
      do k = 2, n_other
         fpixels(other_axes(k)) = other_idx(k)
         lpixels(other_axes(k)) = other_idx(k)
      enddo

      ax_sky1 = pixaxes_in_l(1)
      ax_sky2 = pixaxes_in_l(2)
      ax_block = 0
      if (n_other.ge.1) ax_block = other_axes(1)
      rank_sky1 = 1
      if (ax_sky2.lt.ax_sky1) rank_sky1 = rank_sky1 + 1
      if (n_other.ge.1 .and. ax_block.lt.ax_sky1) rank_sky1 = rank_sky1 + 1
      rank_sky2 = 1
      if (ax_sky1.lt.ax_sky2) rank_sky2 = rank_sky2 + 1
      if (n_other.ge.1 .and. ax_block.lt.ax_sky2) rank_sky2 = rank_sky2 + 1
      rank_block = 3
      if (n_other.ge.1) then
         rank_block = 1
         if (ax_sky1.lt.ax_block) rank_block = rank_block + 1
         if (ax_sky2.lt.ax_block) rank_block = rank_block + 1
      endif
      natural_order = (rank_sky1.eq.1 .and. rank_sky2.eq.2 .and. rank_block.eq.3)

      fitsstat = 0
      blocksize = 1
      unit = 5000
      group = 1
      badval = ieee_value(badval, ieee_quiet_nan)
      call FTOPEN(unit, trim(filename), 0, blocksize, fitsstat)
      if (natural_order) then
         call FTGSVE(unit, group, naxis, naxes_in_l(1:naxis),&
         &fpixels(1:naxis), lpixels(1:naxis), incs(1:naxis),&
         &badval, block_data_in, anyflg, fitsstat)
      else
         dims(rank_sky1) = nx_in
         dims(rank_sky2) = ny_in
         dims(rank_block) = chan_len
         allocate(natural_buf(dims(1), dims(2), dims(3)))
         call FTGSVE(unit, group, naxis, naxes_in_l(1:naxis),&
         &fpixels(1:naxis), lpixels(1:naxis), incs(1:naxis),&
         &badval, natural_buf, anyflg, fitsstat)
         if (fitsstat.eq.0) then
            do c = 1, chan_len
               do j = 1, ny_in
                  do i = 1, nx_in
                     idxvec(rank_sky1) = i
                     idxvec(rank_sky2) = j
                     idxvec(rank_block) = c
                     block_data_in(i,j,c) = natural_buf(idxvec(1), idxvec(2), idxvec(3))
                  enddo
               enddo
            enddo
         endif
         deallocate(natural_buf)
      endif
      call FTCLOS(unit, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to read block (planes ', chan_start,&
         &'-', chan_start+chan_len-1, ') from ', trim(filename)
         call printerror(fitsstat)
         status = -1
      endif
   end subroutine read_one_block

   subroutine write_one_block(out_unit, naxes_in_l, other_axes, other_idx,&
   &n_other, chan_start, chan_len, nx_out, ny_out, block_data_out, status)
      integer, intent(in) :: out_unit, n_other
      integer, intent(in) :: naxes_in_l(:), other_axes(:), other_idx(:)
      integer, intent(in) :: chan_start, chan_len, nx_out, ny_out
      real, intent(in) :: block_data_out(:,:,:)
      integer, intent(inout) :: status

      integer :: fpixels_wr(max_axes), lpixels_wr(max_axes)
      integer :: naxes_wr(max_axes), naxis_wr, k, fitsstat

      if (status.ne.0) return

      naxis_wr = 2 + n_other
      naxes_wr(1) = nx_out
      naxes_wr(2) = ny_out
      fpixels_wr(1) = 1
      fpixels_wr(2) = 1
      lpixels_wr(1) = nx_out
      lpixels_wr(2) = ny_out
      if (n_other.ge.1) then
         naxes_wr(3) = naxes_in_l(other_axes(1))
         fpixels_wr(3) = chan_start
         lpixels_wr(3) = chan_start + chan_len - 1
      endif
      do k = 2, n_other
         naxes_wr(2+k) = naxes_in_l(other_axes(k))
         fpixels_wr(2+k) = other_idx(k)
         lpixels_wr(2+k) = other_idx(k)
      enddo

      fitsstat = 0
      call FTPSSE(out_unit, 1, naxis_wr, naxes_wr(1:naxis_wr),&
      &fpixels_wr(1:naxis_wr), lpixels_wr(1:naxis_wr), block_data_out, status)
      if (status.ne.0) then
         write(*,*) 'ERROR: failed to write block (planes ', chan_start,&
         &'-', chan_start+chan_len-1, ') to output'
      endif
   end subroutine write_one_block

   !===========================================================
   ! stages=convolve or stages=both: restricted 2-sky+1-freq axis handling.
   ! The core new logic -- combines gaussft_mod convolution and AST
   ! resampling in ONE block-processing loop, in either order, with no
   ! intermediate FITS file when both are active.
   !===========================================================

   subroutine process_one_file_restricted(infile, outfile, do_reproject_l,&
   &convolve_first_l, naxis, sky1, sky2, freq_axis, naxes, cdelt1, cdelt2,&
   &nfreq, bmaj_in, bmin_in, bpa_in, isbad, tgt_bmaj, tgt_bmin, tgt_bpa,&
   &reffile_l, pixaxes_ref_l, nx_out_in, ny_out_in, lbnd_out_d, ubnd_out_d,&
   &mem_frac_ram_l, status)
      use, intrinsic :: ieee_arithmetic
      use omp_lib, only: omp_get_max_threads
      use gaussft_mod, only: plan_convolution, convolve_to_beam, destroy_convolution_plan
      character(len=*), intent(in) :: infile, outfile, reffile_l
      logical, intent(in) :: do_reproject_l, convolve_first_l
      integer, intent(in) :: naxis, sky1, sky2, freq_axis, naxes(max_axes)
      integer, intent(in) :: pixaxes_ref_l(2)
      real(dp), intent(in) :: cdelt1, cdelt2
      integer, intent(in) :: nfreq
      real(dp), intent(in) :: bmaj_in(nfreq), bmin_in(nfreq), bpa_in(nfreq)
      logical, intent(in) :: isbad(nfreq)
      real(dp), intent(in) :: tgt_bmaj, tgt_bmin, tgt_bpa
      integer, intent(in) :: nx_out_in, ny_out_in
      double precision, intent(in) :: lbnd_out_d(2), ubnd_out_d(2)
      real, intent(in) :: mem_frac_ram_l
      integer, intent(out) :: status

      integer :: nx_in, ny_in, nx_out, ny_out
      integer :: in_unit, ref_unit, out_unit, fitsstat, blocksize
      logical :: simple, extend
      integer :: naxes_out(max_axes)
      integer(kind=8) :: mem_total_kb, bytes_per_plane, mem_safe_bytes, block_planes64
      integer :: block_planes, chan_start, chan_len, local_iplane, nthreads
      real, allocatable :: block_in(:,:,:), block_out(:,:,:)
      real(dp) :: bpa_in_pixel(nfreq), tgt_bpa_pixel_native, tgt_bpa_pixel_out
      real(dp) :: ref_cdelt1, ref_cdelt2
      integer(kind=8) :: plan_fwd, plan_bwd
      integer :: status_par, ich, k
      real(dp) :: dx_deg, dy_deg, ref_dx_deg, ref_dy_deg
      real(dp) :: nanval
      real(dp) :: plane_native(naxes(sky1), naxes(sky2))
      real(dp) :: plane_out_arr(nx_out_in, ny_out_in)
      ! astResampleR (ast_resampler) is a REAL*4 Fortran interface (matching
      ! reproject_cubes.f90's own block_data_in/block_data_out, both plain
      ! `real`, never real(dp)) -- these two single-precision scratch
      ! arrays exist purely to hand it genuine REAL*4 storage in the
      ! convolve_reproject order, where the plane being resampled comes
      ! from a real(dp) convolution result, not directly from a block_in/
      ! block_out slice the way reproject_convolve order can use as-is.
      real :: plane_native_sp(naxes(sky1), naxes(sky2))
      real :: plane_out_sp(nx_out_in, ny_out_in)

      integer :: t_status, t_wcs_ref, t_skymap_ref, t_skyframe_ref
      integer :: t_naxes_ref(max_axes), t_pixaxes_ref(2)
      integer :: t_wcs_in, t_skymap_in, t_skyframe_in
      integer :: t_naxes_in2(max_axes), t_pixaxes_in(2)
      integer :: t_map_in2ref
      integer :: lbnd_in(2), ubnd_in(2), lbnd_o(2), ubnd_o(2), nbad
      real :: badval_sp
      double precision :: params_dummy(1)

      status = 0
      nx_in = naxes(sky1)
      ny_in = naxes(sky2)
      dx_deg = abs(cdelt1)
      dy_deg = abs(cdelt2)

      if (do_reproject_l) then
         nx_out = nx_out_in
         ny_out = ny_out_in
      else
         nx_out = nx_in
         ny_out = ny_in
      endif

      do ich = 1, nfreq
         call sky_to_pixel_bpa(bpa_in(ich), cdelt1, cdelt2, bpa_in_pixel(ich))
      enddo
      ! Native-grid pixel-frame target BPA (order=convolve_reproject: this
      ! file's own CDELT); output-grid pixel-frame target BPA
      ! (order=reproject_convolve: the reference's CDELT, since every file
      ! shares the reference grid once reprojected) -- only the one that's
      ! actually needed for the selected order gets computed for real, the
      ! other stays at a harmless default when do_reproject_l is false.
      call sky_to_pixel_bpa(tgt_bpa, cdelt1, cdelt2, tgt_bpa_pixel_native)
      tgt_bpa_pixel_out = tgt_bpa_pixel_native
      ref_dx_deg = dx_deg
      ref_dy_deg = dy_deg
      if (do_reproject_l) then
         call read_ref_cdelt(reffile_l, pixaxes_ref_l, ref_cdelt1, ref_cdelt2, status)
         if (status.ne.0) return
         ref_dx_deg = abs(ref_cdelt1)
         ref_dy_deg = abs(ref_cdelt2)
         call sky_to_pixel_bpa(tgt_bpa, ref_cdelt1, ref_cdelt2, tgt_bpa_pixel_out)
      endif

      fitsstat = 0
      in_unit = 210
      call FTOPEN(in_unit, trim(infile), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot reopen input for output header: ', trim(infile)
         status = -1
         return
      endif

      out_unit = 211
      fitsstat = 0
      call FTINIT(out_unit, '!'//trim(outfile), blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot create output file: ', trim(outfile)
         status = -1
         call FTCLOS(in_unit, fitsstat)
         return
      endif

      naxes_out(1:naxis) = naxes(1:naxis)
      naxes_out(sky1) = nx_out
      naxes_out(sky2) = ny_out
      simple = .true.
      extend = .false.
      call FTPHPR(out_unit, simple, -32, naxis, naxes_out(1:naxis), 0, 1, extend, fitsstat)

      if (do_reproject_l) then
         ref_unit = 212
         fitsstat = 0
         call FTOPEN(ref_unit, trim(reffile_l), 0, blocksize, fitsstat)
         call copy_axis_keywords(ref_unit, sky1, out_unit, sky1,&
         &lbnd_out_d(1)-1.0d0, status)
         call copy_axis_keywords(ref_unit, sky2, out_unit, sky2,&
         &lbnd_out_d(2)-1.0d0, status)
         call copy_sky_rotation_matrix(ref_unit, sky1, sky2, out_unit, status)
         call copy_wcs_system_keywords(ref_unit, out_unit, status)
         call FTCLOS(ref_unit, fitsstat)
         call copy_axis_keywords(in_unit, freq_axis, out_unit, freq_axis, 0.0d0, status)
         ! Any OTHER axis (degenerate by this scope's own requirement --
         ! see read_axis_info -- e.g. a size-1 STOKES axis) also needs its
         ! own WCS keywords copied through explicitly: exclude_axis_wcs=
         ! true below skips ALL axis-indexed keywords in the generic copy,
         ! not just sky1/sky2/freq_axis, so without this loop a degenerate
         ! axis's CTYPE/CRVAL/CRPIX/CDELT would simply be lost. Caught by
         ! this tool's own chaining-equivalence verification (a diff
         ! against the two-step disk-based reference, which does carry
         ! this through via reproject_cubes' own general "other axes"
         ! handling).
         do k = 1, naxis
            if (k.ne.sky1 .and. k.ne.sky2 .and. k.ne.freq_axis) then
               call copy_axis_keywords(in_unit, k, out_unit, k, 0.0d0, status)
            endif
         enddo
      endif

      call copy_generic_header_match(in_unit, out_unit, do_reproject_l, .true., status)
      call FTPKYD(out_unit, 'BMAJ', tgt_bmaj/3600.0d0, 13,&
      &'common-resolution major axis FWHM (deg)', fitsstat)
      call FTPKYD(out_unit, 'BMIN', tgt_bmin/3600.0d0, 13,&
      &'common-resolution minor axis FWHM (deg)', fitsstat)
      call FTPKYD(out_unit, 'BPA', tgt_bpa, 13,&
      &'common-resolution position angle (deg)', fitsstat)
      if (do_reproject_l) then
         call FTPHIS(out_unit, 'match_cubes: convolved and reprojected from '//&
         &trim(infile), fitsstat)
      else
         call FTPHIS(out_unit, 'match_cubes: convolved from '//trim(infile)//&
         &' to a common resolution', fitsstat)
      endif
      call FTCLOS(in_unit, fitsstat)

      call get_mem_total_kb(mem_total_kb)
      bytes_per_plane = int(4,8) * (int(nx_in,8)*int(ny_in,8) +&
      &int(nx_out,8)*int(ny_out,8))
      mem_safe_bytes = int(real(mem_frac_ram_l,8) * real(mem_total_kb,8) * 1024.0d0, 8)
      block_planes64 = max(1_8, mem_safe_bytes / bytes_per_plane)
      block_planes = int(min(block_planes64, int(nfreq,8)))
      if (block_planes.lt.1) block_planes = 1

      write(*,'(A,A,A,I0,A,I0,A)') 'Writing ', trim(outfile), ': ', nfreq,&
      &' plane(s), in blocks of up to ', block_planes, ' plane(s)'
      if (block_planes.lt.omp_get_max_threads()) then
         write(*,'(A,I0,A,I0,A)') 'WARNING: mem_frac_ram limits blocks to ',&
         &block_planes, ' plane(s), below the ', omp_get_max_threads(),&
         &' threads available -- parallelism is reduced; raise mem_frac_ram',&
         &' for full speedup if memory allows.'
      endif

      allocate(block_in(nx_in, ny_in, block_planes))
      allocate(block_out(nx_out, ny_out, block_planes))

      ! FFTW's plan is sized for one specific (nx,ny) and MUST NOT be
      ! executed against arrays of a different size (silent heap
      ! corruption, not a clean error -- caught the hard way via a
      ! munmap_chunk() crash before this fix). Convolution happens on the
      ! native grid for order=convolve_reproject (or always, for
      ! stages=convolve alone, where nx_out/ny_out were already set equal
      ! to nx_in/ny_in above) -- but on the OUTPUT grid for
      ! order=reproject_convolve, which can be a genuinely different size
      ! (footprint_mode=intersection/union crops or grows it). Plan
      ! whichever grid convolution will actually run on.
      if (.not. do_reproject_l .or. convolve_first_l) then
         call plan_convolution(nx_in, ny_in, plan_fwd, plan_bwd)
      else
         call plan_convolution(nx_out, ny_out, plan_fwd, plan_bwd)
      endif

      status_par = 0
      nthreads = max(1, min(omp_get_max_threads(), block_planes))
      chan_start = 1
      do while (chan_start.le.nfreq)
         chan_len = min(block_planes, nfreq - chan_start + 1)

         call read_freq_block(infile, naxis, sky1, sky2, freq_axis,&
         &naxes, chan_start, chan_len, nx_in, ny_in, block_in(:,:,1:chan_len), status_par)
         if (status_par.ne.0) exit

         !$omp parallel num_threads(nthreads) default(none)&
         !$omp& shared(chan_len, nx_in, ny_in, nx_out, ny_out, block_in,&
         !$omp& block_out, isbad, chan_start, plan_fwd, plan_bwd, dx_deg,&
         !$omp& dy_deg, ref_dx_deg, ref_dy_deg, bmaj_in, bmin_in,&
         !$omp& bpa_in_pixel, tgt_bmaj, tgt_bmin, tgt_bpa_pixel_native,&
         !$omp& tgt_bpa_pixel_out, status_par, do_reproject_l, convolve_first_l,&
         !$omp& reffile_l, infile, lbnd_out_d, ubnd_out_d)&
         !$omp& private(local_iplane, ich, k, nanval, plane_native,&
         !$omp& plane_out_arr, plane_native_sp, plane_out_sp, t_status,&
         !$omp& t_wcs_ref, t_skymap_ref, t_skyframe_ref, t_naxes_ref,&
         !$omp& t_pixaxes_ref, t_wcs_in, t_skymap_in, t_skyframe_in,&
         !$omp& t_naxes_in2, t_pixaxes_in, t_map_in2ref, lbnd_in, ubnd_in,&
         !$omp& lbnd_o, ubnd_o, badval_sp, params_dummy, nbad)

         t_status = 0
         if (do_reproject_l) then
            ! Own private AST Mapping per thread -- required, see
            ! process_one_file_general's own comment on why (this Fortran
            ! AST binding has no lock/unlock for cross-thread sharing).
            call ast_begin(t_status)
            call load_wcs(reffile_l, t_wcs_ref, t_naxes_ref, t_status)
            call extract_sky_mapping(t_wcs_ref, t_skymap_ref, t_skyframe_ref,&
            &t_pixaxes_ref, t_status)
            call load_wcs(infile, t_wcs_in, t_naxes_in2, t_status)
            call extract_sky_mapping(t_wcs_in, t_skymap_in, t_skyframe_in,&
            &t_pixaxes_in, t_status)
            call compose_pix2pix(t_skymap_in, t_skyframe_in, t_skymap_ref,&
            &t_skyframe_ref, t_map_in2ref, t_status)
            if (t_status.ne.0) then
               !$omp atomic write
               status_par = -1
            endif
         endif

         !$omp do schedule(dynamic)
         do local_iplane = 1, chan_len
            ich = chan_start + local_iplane - 1
            if (t_status.ne.0 .or. status_par.ne.0) cycle
            if (isbad(ich)) then
               nanval = ieee_value(1.0_dp, ieee_quiet_nan)
               block_out(:,:,local_iplane) = real(nanval)
               cycle
            endif

            if (.not. do_reproject_l) then
               ! stages=convolve alone: no resampling, output grid ==
               ! native grid, identical to today's standalone
               ! convolve_cubes.
               plane_native = real(block_in(:,:,local_iplane), dp)
               call convolve_to_beam(plan_fwd, plan_bwd, plane_native, nx_in, ny_in,&
               &dx_deg, dy_deg, bmaj_in(ich)/3600.0_dp, bmin_in(ich)/3600.0_dp,&
               &bpa_in_pixel(ich), tgt_bmaj/3600.0_dp, tgt_bmin/3600.0_dp,&
               &tgt_bpa_pixel_native, plane_native, k)
               block_out(:,:,local_iplane) = real(plane_native)
               if (k.ne.0) then
                  !$omp atomic write
                  status_par = -1
               endif
               cycle
            endif

            lbnd_in(1) = 1
            lbnd_in(2) = 1
            ubnd_in(1) = nx_in
            ubnd_in(2) = ny_in
            lbnd_o(1) = nint(lbnd_out_d(1))
            lbnd_o(2) = nint(lbnd_out_d(2))
            ubnd_o(1) = nint(ubnd_out_d(1))
            ubnd_o(2) = nint(ubnd_out_d(2))
            params_dummy(1) = 0.0d0

            if (convolve_first_l) then
               ! order=convolve_reproject (default): convolve on the
               ! native grid using this file's OWN dx/dy, then resample
               ! the (now low-pass-filtered) native-grid plane onto the
               ! output grid.
               plane_native = real(block_in(:,:,local_iplane), dp)
               call convolve_to_beam(plan_fwd, plan_bwd, plane_native, nx_in, ny_in,&
               &dx_deg, dy_deg, bmaj_in(ich)/3600.0_dp, bmin_in(ich)/3600.0_dp,&
               &bpa_in_pixel(ich), tgt_bmaj/3600.0_dp, tgt_bmin/3600.0_dp,&
               &tgt_bpa_pixel_native, plane_native, k)
               if (k.ne.0) then
                  !$omp atomic write
                  status_par = -1
                  cycle
               endif
               ! ast_resampler needs genuine REAL*4 storage, not a real(dp)
               ! array or a throwaway real(plane_native) conversion
               ! expression (which is a temporary, not a reference it can
               ! read/write through) -- convert into the dedicated
               ! single-precision scratch arrays declared above.
               plane_native_sp = real(plane_native)
               badval_sp = ieee_value(badval_sp, ieee_quiet_nan)
               nbad = ast_resampler(t_map_in2ref, 2, lbnd_in, ubnd_in,&
               &plane_native_sp, plane_native_sp,&
               &ast__linear, ast_null, params_dummy, 0, 0.0d0, 100, badval_sp,&
               &2, lbnd_o, ubnd_o, lbnd_o, ubnd_o,&
               &plane_out_sp, plane_out_sp, t_status)
               if (t_status.ne.0) then
                  !$omp atomic write
                  status_par = -1
                  cycle
               endif
               block_out(:,:,local_iplane) = plane_out_sp
            else
               ! order=reproject_convolve: resample first onto the output
               ! grid, then convolve there using the OUTPUT grid's own
               ! (reference) dx/dy.
               badval_sp = ieee_value(badval_sp, ieee_quiet_nan)
               nbad = ast_resampler(t_map_in2ref, 2, lbnd_in, ubnd_in,&
               &block_in(:,:,local_iplane), block_in(:,:,local_iplane),&
               &ast__linear, ast_null, params_dummy, 0, 0.0d0, 100, badval_sp,&
               &2, lbnd_o, ubnd_o, lbnd_o, ubnd_o,&
               &block_out(:,:,local_iplane), block_out(:,:,local_iplane), t_status)
               if (t_status.ne.0) then
                  !$omp atomic write
                  status_par = -1
                  cycle
               endif
               plane_out_arr = real(block_out(:,:,local_iplane), dp)
               call convolve_to_beam(plan_fwd, plan_bwd, plane_out_arr, nx_out, ny_out,&
               &ref_dx_deg, ref_dy_deg, bmaj_in(ich)/3600.0_dp, bmin_in(ich)/3600.0_dp,&
               &bpa_in_pixel(ich), tgt_bmaj/3600.0_dp, tgt_bmin/3600.0_dp,&
               &tgt_bpa_pixel_out, plane_out_arr, k)
               if (k.ne.0) then
                  !$omp atomic write
                  status_par = -1
                  cycle
               endif
               block_out(:,:,local_iplane) = real(plane_out_arr)
            endif
         enddo
         !$omp end do

         if (do_reproject_l) call ast_end(t_status)
         !$omp end parallel
         if (status_par.ne.0) exit

         call write_freq_block(out_unit, naxis, sky1, sky2, freq_axis,&
         &naxes_out, chan_start, chan_len, nx_out, ny_out, block_out(:,:,1:chan_len), status_par)
         if (status_par.ne.0) exit

         chan_start = chan_start + chan_len
      enddo

      call destroy_convolution_plan(plan_fwd, plan_bwd)
      deallocate(block_in, block_out)

      if (status_par.ne.0) then
         write(*,*) 'ERROR: failed to process/write one or more planes for: ', trim(infile)
         status = -1
         call FTCLOS(out_unit, fitsstat)
         return
      endif

      ! CASAMBM/BEAMS: always attached whenever convolve is active
      ! (stages=convolve or both, alone or chained) -- see
      ! write_beams_table_match's own comment for why, and
      ! convolve_cubes.f90's identical write_beams_table this is a
      ! verbatim port of.
      fitsstat = 0
      call FTPKYL(out_unit, 'CASAMBM', .true.,&
      &'Multiple beams per plane (see BEAMS ext)', fitsstat)
      call write_beams_table_match(out_unit, nfreq, isbad, tgt_bmaj,&
      &tgt_bmin, tgt_bpa, status_par)
      if (status_par.ne.0) then
         status = -1
         call FTCLOS(out_unit, fitsstat)
         return
      endif

      call FTCLOS(out_unit, fitsstat)
   end subroutine process_one_file_restricted

   subroutine write_beams_table_match(unit, nfreq, isbad, tgt_bmaj,&
   &tgt_bmin, tgt_bpa, status)
      !! Verbatim port of convolve_cubes.f90's own write_beams_table --
      !! see that subroutine's comment for the full rationale (CASA-style
      !! 5-column BEAMS layout confirmed against a real ASKAP cube, one
      !! row per channel: common target beam for a good channel, the
      !! same tiny(1.0) degenerate sentinel CASA itself uses for a bad
      !! one). Applies whenever this file's own convolve stage is active
      !! (stages=convolve or both), regardless of whether reproject also
      !! ran, since reprojection never changes which channels are
      !! good/bad or what beam they end up at.
      integer, intent(inout) :: unit
      integer, intent(in) :: nfreq
      logical, intent(in) :: isbad(nfreq)
      real(dp), intent(in) :: tgt_bmaj, tgt_bmin, tgt_bpa
      integer, intent(out) :: status

      character(len=8) :: ttype(5), tform(5), tunit_(5)
      real, allocatable :: col_bmaj(:), col_bmin(:), col_bpa(:)
      integer, allocatable :: col_chan(:), col_pol(:)
      integer :: fitsstat, ich, colnum

      ttype = (/'BMAJ    ', 'BMIN    ', 'BPA     ', 'CHAN    ', 'POL     '/)
      tform = (/'1E      ', '1E      ', '1E      ', '1J      ', '1J      '/)
      tunit_ = (/'arcsec  ', 'arcsec  ', 'deg     ', '        ', '        '/)

      allocate(col_bmaj(nfreq), col_bmin(nfreq), col_bpa(nfreq),&
      &col_chan(nfreq), col_pol(nfreq))
      do ich = 1, nfreq
         col_chan(ich) = ich - 1
         col_pol(ich) = 0
         if (isbad(ich)) then
            col_bmaj(ich) = tiny(1.0)
            col_bmin(ich) = tiny(1.0)
            col_bpa(ich) = 0.0
         else
            col_bmaj(ich) = real(tgt_bmaj)
            col_bmin(ich) = real(tgt_bmin)
            col_bpa(ich) = real(tgt_bpa)
         endif
      enddo

      fitsstat = 0
      call FTIBIN(unit, nfreq, 5, ttype, tform, tunit_, 'BEAMS', 0, fitsstat)
      call FTGCNO(unit, .false., 'BMAJ', colnum, fitsstat)
      call FTPCLE(unit, colnum, 1, 1, nfreq, col_bmaj, fitsstat)
      call FTGCNO(unit, .false., 'BMIN', colnum, fitsstat)
      call FTPCLE(unit, colnum, 1, 1, nfreq, col_bmin, fitsstat)
      call FTGCNO(unit, .false., 'BPA', colnum, fitsstat)
      call FTPCLE(unit, colnum, 1, 1, nfreq, col_bpa, fitsstat)
      call FTGCNO(unit, .false., 'CHAN', colnum, fitsstat)
      call FTPCLJ(unit, colnum, 1, 1, nfreq, col_chan, fitsstat)
      call FTGCNO(unit, .false., 'POL', colnum, fitsstat)
      call FTPCLJ(unit, colnum, 1, 1, nfreq, col_pol, fitsstat)
      deallocate(col_bmaj, col_bmin, col_bpa, col_chan, col_pol)

      status = 0
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to write BEAMS binary table extension'
         status = -1
      endif
   end subroutine write_beams_table_match

   subroutine read_ref_cdelt(reffile_l, pixaxes_ref_l, ref_cdelt1, ref_cdelt2, status)
      !! Only the reference file's own CDELT for its 2 sky axes -- needed
      !! for order=reproject_convolve's post-resample convolution step,
      !! which must use the OUTPUT grid's own pixel scale, not the
      !! per-file native one. Deliberately NOT read_axis_info: that
      !! subroutine (a) requires a FREQ axis to be present, wrong for a
      !! reference file that may be a plain 2D image with no spectral
      !! axis at all, and (b) re-derives which axis is RA/DEC from CTYPE,
      !! redundant here since the caller already has pixaxes_ref_l -- the
      !! AST-derived sky axis numbers for this exact file, computed once
      !! during the reproject pre-scan (extract_sky_mapping) -- so this
      !! reads CDELT directly for those two axis numbers instead of
      !! re-detecting them a second, independent way.
      character(len=*), intent(in) :: reffile_l
      integer, intent(in) :: pixaxes_ref_l(2)
      real(dp), intent(out) :: ref_cdelt1, ref_cdelt2
      integer, intent(out) :: status
      integer :: unit, blocksize, fitsstat
      character(len=68) :: comment
      character(len=8) :: axstr

      status = 0
      fitsstat = 0
      unit = 220
      call FTOPEN(unit, trim(reffile_l), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open reference file: ', trim(reffile_l)
         status = -1
         return
      endif
      write(axstr,'(I0)') pixaxes_ref_l(1)
      fitsstat = 0
      call FTGKYD(unit, 'CDELT'//trim(axstr), ref_cdelt1, comment, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: missing CDELT for the reference''s RA axis in: ',&
         &trim(reffile_l)
         status = -1
         call FTCLOS(unit, fitsstat)
         return
      endif
      write(axstr,'(I0)') pixaxes_ref_l(2)
      fitsstat = 0
      call FTGKYD(unit, 'CDELT'//trim(axstr), ref_cdelt2, comment, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: missing CDELT for the reference''s DEC axis in: ',&
         &trim(reffile_l)
         status = -1
         call FTCLOS(unit, fitsstat)
         return
      endif
      call FTCLOS(unit, fitsstat)
   end subroutine read_ref_cdelt

   subroutine read_freq_block(filename, naxis, sky1, sky2, freq_axis,&
   &naxes, chan_start, chan_len, nx, ny, block_data, status)
      use, intrinsic :: ieee_arithmetic
      character(len=*), intent(in) :: filename
      integer, intent(in) :: naxis, sky1, sky2, freq_axis, naxes(max_axes)
      integer, intent(in) :: chan_start, chan_len, nx, ny
      real, intent(out) :: block_data(:,:,:)
      integer, intent(inout) :: status

      integer :: unit, blocksize, fitsstat, group
      integer :: fpixels(max_axes), lpixels(max_axes), incs(max_axes)
      logical :: anyflg
      real :: badval
      integer :: rank_sky1, rank_sky2, rank_freq
      integer :: dims(3), idxvec(3), i, j, c
      logical :: natural_order
      real, allocatable :: natural_buf(:,:,:)

      if (status.ne.0) return

      fpixels(1:naxis) = 1
      lpixels(1:naxis) = 1
      incs(1:naxis) = 1
      lpixels(sky1) = nx
      lpixels(sky2) = ny
      fpixels(freq_axis) = chan_start
      lpixels(freq_axis) = chan_start + chan_len - 1

      rank_sky1 = 1
      if (sky2.lt.sky1) rank_sky1 = rank_sky1 + 1
      if (freq_axis.lt.sky1) rank_sky1 = rank_sky1 + 1
      rank_sky2 = 1
      if (sky1.lt.sky2) rank_sky2 = rank_sky2 + 1
      if (freq_axis.lt.sky2) rank_sky2 = rank_sky2 + 1
      rank_freq = 1
      if (sky1.lt.freq_axis) rank_freq = rank_freq + 1
      if (sky2.lt.freq_axis) rank_freq = rank_freq + 1
      natural_order = (rank_sky1.eq.1 .and. rank_sky2.eq.2 .and. rank_freq.eq.3)

      fitsstat = 0
      blocksize = 1
      unit = 5100
      group = 1
      badval = ieee_value(badval, ieee_quiet_nan)
      call FTOPEN(unit, trim(filename), 0, blocksize, fitsstat)
      if (natural_order) then
         call FTGSVE(unit, group, naxis, naxes(1:naxis), fpixels(1:naxis),&
         &lpixels(1:naxis), incs(1:naxis), badval, block_data, anyflg, fitsstat)
      else
         dims(rank_sky1) = nx
         dims(rank_sky2) = ny
         dims(rank_freq) = chan_len
         allocate(natural_buf(dims(1), dims(2), dims(3)))
         call FTGSVE(unit, group, naxis, naxes(1:naxis), fpixels(1:naxis),&
         &lpixels(1:naxis), incs(1:naxis), badval, natural_buf, anyflg, fitsstat)
         if (fitsstat.eq.0) then
            do c = 1, chan_len
               do j = 1, ny
                  do i = 1, nx
                     idxvec(rank_sky1) = i
                     idxvec(rank_sky2) = j
                     idxvec(rank_freq) = c
                     block_data(i,j,c) = natural_buf(idxvec(1), idxvec(2), idxvec(3))
                  enddo
               enddo
            enddo
         endif
         deallocate(natural_buf)
      endif
      call FTCLOS(unit, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to read block (planes ', chan_start, '-',&
         &chan_start+chan_len-1, ') from ', trim(filename)
         status = -1
      endif
   end subroutine read_freq_block

   subroutine write_freq_block(out_unit, naxis, sky1, sky2, freq_axis,&
   &naxes_out_l, chan_start, chan_len, nx, ny, block_data, status)
      integer, intent(in) :: out_unit, naxis, sky1, sky2, freq_axis, naxes_out_l(max_axes)
      integer, intent(in) :: chan_start, chan_len, nx, ny
      real, intent(in) :: block_data(:,:,:)
      integer, intent(inout) :: status

      integer :: fpixels_wr(max_axes), lpixels_wr(max_axes), fitsstat

      if (status.ne.0) return

      fpixels_wr(1:naxis) = 1
      lpixels_wr(1:naxis) = 1
      lpixels_wr(sky1) = nx
      lpixels_wr(sky2) = ny
      fpixels_wr(freq_axis) = chan_start
      lpixels_wr(freq_axis) = chan_start + chan_len - 1

      fitsstat = 0
      call FTPSSE(out_unit, 1, naxis, naxes_out_l(1:naxis), fpixels_wr(1:naxis),&
      &lpixels_wr(1:naxis), block_data, status)
      if (status.ne.0) then
         write(*,*) 'ERROR: failed to write block (planes ', chan_start, '-',&
         &chan_start+chan_len-1, ') to output'
      endif
   end subroutine write_freq_block

   subroutine get_mem_total_kb(mem_total_kb)
      integer(kind=8), intent(out) :: mem_total_kb
      integer :: mem_unit, ios_mem
      character(len=128) :: mem_line
      integer(kind=8) :: mem_kb_tmp

      mem_total_kb = 0_8
      open(newunit=mem_unit, file='/proc/meminfo', status='old', iostat=ios_mem)
      if (ios_mem.eq.0) then
         do
            read(mem_unit, '(A)', iostat=ios_mem) mem_line
            if (ios_mem.ne.0) exit
            if (index(mem_line, 'MemTotal:').eq.1) then
               read(mem_line(index(mem_line,':')+1:), *, iostat=ios_mem) mem_kb_tmp
               if (ios_mem.eq.0) mem_total_kb = mem_kb_tmp
               exit
            endif
         enddo
         close(mem_unit)
      endif
      if (mem_total_kb.le.0_8) mem_total_kb = 4194304_8
   end subroutine get_mem_total_kb

end program match_cubes
