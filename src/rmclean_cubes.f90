! rmclean_cubes -- standalone RM-CLEAN tool, driving rmclean_mod (src/
! rmclean.f90) against a real dirty AMP/PHA cube pair rm_synthesis itself
! wrote (<outfile>.AMP.RMCUBE.FITS/.PHA.RMCUBE.FITS), plus its own
! .MASK.CUBE.FITS and CHANFREQ binary table (rm_synthesis.f90:4762-4826).
! planning/RMCLEAN_INTEGRATION_PLAN.md T2. rmclean_mod itself is pure
! computation (no FITS I/O, no config) -- this program is exactly the
! "main program" its own header comment anticipates, mirroring
! convolve_cubes.f90/reproject_cubes.f90's own split.
!
! Two axis-3 sizes matter here and must not be confused: ampfile/phafile's
! own NAXIS3 is nrm (the RM/Faraday-depth axis rm_synthesis already
! synthesized, CTYPE3='FDEP'); maskfile's own NAXIS3 is nchan (the INPUT
! frequency-channel axis, CTYPE3='FREQ', one 0/1 flag per pixel per
! channel saying whether that channel contributed to that pixel's own
! dirty spectrum), with per-channel true frequency/band given by its
! CHANFREQ table (columns CHAN/FREQ/BAND, 0-indexed CHAN/BAND, FREQ in
! Hz). This program CLEANs the (nrm-length) dirty spectrum at each pixel;
! nchan only matters for building that pixel's own RMSF (a function of
! which channels were actually valid there).
!
! lsq_ref_native: the phase reference the dirty AMP/PHA cube was ACTUALLY
! built at, read from its own LSQREF header keyword (rm_synthesis now
! writes this -- rm_synthesis.f90's own cfg%lsq_ref_mode/lsq_ref_fixed_value,
! recorded so a reader never has to assume a value). Falls back to 0.0
! (with a printed warning, not a silent assumption) for cubes written
! before this keyword existed -- rm_synthesis's own historical convention
! was unconditionally lsq_ref=0 (verified directly at rm_synthesis_mod.
! f90's own extract_general_setup: phi_tmp=omega*t(kk) used raw L_sq, no
! subtraction), so 0.0 is the CORRECT default for those older cubes, not
! an arbitrary guess.
!
! lsq_ref_compute -- the reference this program's OWN RMSF table/CLEAN
! computation actually uses -- IS a free, independent config choice
! (mode+value, mirroring rmclean_mod's own get_lsq_ref_compute exactly),
! defaulting to 'native' (= lsq_ref_native, no extra work). Choosing a
! DIFFERENT value is fully supported: derotate_to_lsq_ref is an exact,
! lossless phase rotation of the ALREADY-SAMPLED dirty spectrum (verified
! |P|-invariant to 1e-4 relative, tests/test_rmclean_lsqref_flex.f90), so
! re-expressing it at any other reference before building the table costs
! nothing in accuracy. It also costs nothing in COMPUTE beyond one cheap
! O(nrm) derotation per pixel: once the RM grid (CDELT3/nrm) already
! exists, build_rmsf_offset_table/clean_complex/restore_clean's own cost
! depends only on nrm/niter/table_oversample, none of which depend on
! lsq_ref_compute -- the grid-size WIN a favourable reference buys
! (get_drm's own bound, minimized by e.g. mode=mid) only applies UPSTREAM,
! when rm_synthesis itself is deciding what CDELT3 to use (its own
! lsq_ref_mode option, added alongside this one). So there is no reason
! to forbid the choice here, only to apply it correctly: build the RMSF
! table at whatever lsq_ref_compute ends up in force, matching the
! (possibly-derotated) dirty spectrum exactly (build_rmsf_offset_table's
! own correctness requirement: R(delta) is only a pure function of offset
! for a FIXED, SHARED reference between the dirty map and its own
! matched-filter template -- a mismatch is a real phase distortion, not
! just a rotation, per compute_dirty_rmbeam_direct's own doc comment).
!
! Gate 0 (T3b, planning/RMCLEAN_INTEGRATION_PLAN.md): validates the
! cube's existing RM grid against the RM-RESOLUTION requirement --
! |CDELT3| <= fwhm/min_samples_per_fwhm, with the fwhm from
! compute_rmsf_fwhm_multiband, an lsq_ref-INDEPENDENT quantity set by the
! lambda^2 span alone. This is a deliberate T3b replacement of the
! original get_drm-based (max_offset/carrier-driven, lsq_ref-dependent)
! gate: since T3 replaced clean_complex's own peak refinement with the
! model-based matched filter (refine_peak_matched_filter's own doc
! comment in src/rmclean.f90), the stored grid no longer needs to
! resolve the fast Re/Im carrier at all -- it only needs enough samples
! per RMSF fwhm for the discrete peak (and its anchor neighbours) to
! land reliably inside the right resolution element (confirmed
! empirically down to 2 samples/fwhm, tests/test_matched_filter_refine.
! f90). Unlike reproject_cubes/convolve_cubes, this tool cannot resample
! the RM axis -- it is fixed by whatever CDELT3 rm_synthesis already
! wrote -- so a too-coarse grid is refused loudly (no silent fallback --
! re-gridding would mean re-synthesizing, out of scope for a CLEAN-only
! tool).
!
! Usage: rmclean_cubes ampfile=<f> phafile=<f> maskfile=<f> outfile=<base>
!    [niter=<n>] [gain=<g>] threshold=<t>
!    [lsq_ref_report_mode=intrinsic|centroid|min|max|mid|fixed]
!    [lsq_ref_report_value=<v>] [min_samples_per_fwhm=<f>]
!    [refine_nsigma=<f>] [table_oversample=<n>] [restore_fwhm=<v>]
!    or: rmclean_cubes --config <cfgfile>
!    or: rmclean_cubes --help | -h
! Full usage text in print_usage below (shared by --help and the
! argument-error path, same convention as reproject_cubes.f90/
! convolve_cubes.f90).
!
! Stage A (this version): whole-cube-in-memory, serial, one rmsf_table_t
! built fresh per pixel from that pixel's own valid-channel subset --
! correct, but without the memory-budgeted block I/O, OpenMP
! parallelization, or mask-pattern table cache planning/
! RMCLEAN_INTEGRATION_PLAN.md's T2 design calls for for production-scale
! cubes (those are Stage B, added on top of this once Stage A's
! correctness is verified against tests/thesis_scenario_rmclean.f90's own
! tolerances -- see that ticket's own verification plan).
program rmclean_cubes
   use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
   use rmclean_mod
   implicit none

   character(len=512) :: ampfile, phafile, maskfile, outfile
   integer :: niter
   real(sp) :: gain, threshold
   logical :: have_threshold
   real(sp) :: min_samples_per_fwhm, refine_nsigma
   integer :: table_oversample
   logical :: have_restore_fwhm
   real(sp) :: restore_fwhm_override
   integer :: lsq_ref_report_mode_sel
   logical :: have_lsq_ref_report_value
   real(sp) :: lsq_ref_report_value
   character(len=16) :: lsq_ref_compute_mode
   logical :: have_lsq_ref_compute_value
   real(sp) :: lsq_ref_compute_value

   integer :: status
   integer :: nx, ny, nrm, nchan
   real(dp) :: cdelt3_amp, crval3_amp
   real(sp), allocatable :: l_sq(:)
   integer, allocatable :: band_id(:)
   integer, allocatable :: band_offset(:), band_nz(:)
   integer :: n_bands

   real(sp) :: lsq_ref_native, lsq_ref_compute, lsq_ref_report
   real(sp) :: drm_required, fwhm_rm, fwhm_data
   logical :: have_lsqref_keyword
   real(sp), allocatable :: rm_samp(:)

   real(sp), allocatable :: re_cube(:,:,:), im_cube(:,:,:)
   integer(kind=1), allocatable :: mask_cube(:,:,:)

   real(sp), allocatable :: clean_re_cube(:,:,:), clean_im_cube(:,:,:)
   real(sp), allocatable :: resid_re_cube(:,:,:), resid_im_cube(:,:,:)
   real(sp), allocatable :: restored_re_cube(:,:,:), restored_im_cube(:,:,:)

   integer :: ix, iy, k
   integer(kind=8) :: restore_plan_fwd, restore_plan_bwd
   integer :: n_pixels_done
   integer :: mask_pattern_cache_max

   ! --- Mask-pattern -> rmsf_table_t cache (planning/RMCLEAN_INTEGRATION_
   ! PLAN.md decision 10) -- built once, serially, by build_mask_pattern_
   ! cache below, BEFORE the parallel per-pixel CLEAN loop starts, so
   ! every thread's own lookup (cache_lookup_readonly) is read-only and
   ! race-free. See table_cache_entry_t's own comment.
   type :: table_cache_entry_t
      integer(kind=8) :: hash = 0_8
      integer(kind=1), allocatable :: pattern(:)
      type(rmsf_table_t) :: table
   end type table_cache_entry_t
   type(table_cache_entry_t), allocatable :: cache_entries(:)
   integer, allocatable :: cache_buckets(:)
   integer :: n_cache_entries, n_cache_buckets

   call parse_args(status)
   if (status.ne.0) stop 1

   call read_chanfreq(maskfile, l_sq, band_id, nchan, status)
   if (status.ne.0) stop 1
   call derive_band_segmentation(band_id, nchan, band_offset, band_nz, n_bands, status)
   if (status.ne.0) stop 1

   call read_amp_pha_geometry(ampfile, phafile, nx, ny, nrm, cdelt3_amp,&
   &crval3_amp, status)
   if (status.ne.0) stop 1

   call read_lsqref_keyword(ampfile, lsq_ref_native, have_lsqref_keyword)
   if (have_lsqref_keyword) then
      write(*,'(A,F0.6)') 'Read LSQREF from '//trim(ampfile)//': ',&
      &lsq_ref_native
   else
      write(*,'(A)') 'WARNING: no LSQREF keyword in '//trim(ampfile)//&
      &' -- assuming lsq_ref_native=0.0 (this project''s historical'//&
      &' convention for cubes written before this keyword existed). If'//&
      &' this cube was actually built at a different reference, this'//&
      &' program''s own Gate 0 check and RMSF table will be WRONG.'
      lsq_ref_native = 0.0_sp
   endif

   ! The DATA's own RM resolution -- needed by Gate 0 just below and as
   ! the default restoring-beam width further down.
   call compute_rmsf_fwhm_multiband(l_sq, nchan, band_offset, band_nz,&
   &n_bands, fwhm_data)

   ! --- Gate 0 (T3b): validate the EXISTING RM grid against the
   ! RM-RESOLUTION requirement, do not resample it ---
   ! |CDELT3| <= fwhm/min_samples_per_fwhm, fwhm from
   ! compute_rmsf_fwhm_multiband -- an lsq_ref-INDEPENDENT quantity set
   ! by the lambda^2 span alone (this file's own top comment: since T3's
   ! model-based peak refinement, the stored grid no longer needs to
   ! resolve the fast, lsq_ref-dependent Re/Im carrier at all; it only
   ! needs the discrete peak and its anchor neighbours to land reliably
   ! within the right resolution element). The DATA's own fwhm is used
   ! here (never the restore_fwhm override below, which only shapes the
   ! restoring beam, not what the data can actually resolve).
   drm_required = fwhm_data/max(min_samples_per_fwhm, 1.0e-3_sp)
   if (abs(real(cdelt3_amp, sp)).gt.drm_required) then
      write(*,'(A)') 'FATAL: Gate 0 failed -- the existing RM grid in '//&
      &trim(ampfile)//' does not resolve the RMSF.'
      write(*,'(A,F0.6,A,F0.6,A,F0.3,A,F0.6,A)') 'Existing |CDELT3| = ',&
      &abs(real(cdelt3_amp, sp)), ' rad/m^2, but fwhm=', fwhm_data,&
      &' rad/m^2 at min_samples_per_fwhm=', min_samples_per_fwhm,&
      &' requires <= ', drm_required, ' rad/m^2.'
      write(*,'(A)') 'Re-run rm_synthesis with a finer cdelt3 (smaller'//&
      &' RM spacing), or lower min_samples_per_fwhm= here (floor 1) if'//&
      &' the coarser grid is genuinely acceptable -- this program does'//&
      &' not resample the RM axis itself.'
      stop 1
   endif
   write(*,'(A,F0.6,A,F0.6,A)') 'Gate 0 OK: existing |CDELT3|=',&
   &abs(real(cdelt3_amp, sp)), ' rad/m^2 <= required ', drm_required,&
   &' rad/m^2 (RM-resolution criterion).'

   allocate(rm_samp(nrm))
   do k = 1, nrm
      rm_samp(k) = real(crval3_amp, sp) + real(k-1, sp)*real(cdelt3_amp, sp)
   enddo

   ! lsq_ref_compute: this program's OWN free choice of reference for the
   ! RMSF table/CLEAN computation (this file's own top comment) --
   ! 'native' (default) means "use lsq_ref_native, no derotation needed".
   if (trim(lsq_ref_compute_mode).eq.'native') then
      lsq_ref_compute = lsq_ref_native
   else
      call get_lsq_ref_compute(l_sq, nchan,&
      &lsq_ref_compute_mode_to_int(lsq_ref_compute_mode), lsq_ref_compute,&
      &fixed_value=merge(lsq_ref_compute_value, 0.0_sp,&
      &have_lsq_ref_compute_value))
   endif
   if (lsq_ref_compute.ne.lsq_ref_native) then
      write(*,'(A,F0.6,A,F0.6)') 'lsq_ref_compute=', lsq_ref_compute,&
      &' differs from lsq_ref_native=', lsq_ref_native,&
      &' -- each pixel''s dirty spectrum will be derotated before CLEANing.'
   endif

   call get_lsq_ref_report(l_sq, nchan, lsq_ref_report_mode_sel,&
   &lsq_ref_report, fixed_value=merge(lsq_ref_report_value, 0.0_sp,&
   &have_lsq_ref_report_value))

   if (have_restore_fwhm) then
      fwhm_rm = restore_fwhm_override
   else
      fwhm_rm = fwhm_data
   endif
   write(*,'(A,F0.6,A)') 'Restoring beam FWHM = ', fwhm_rm, ' rad/m^2'

   call read_cube(ampfile, phafile, nx, ny, nrm, re_cube, im_cube, status)
   if (status.ne.0) stop 1
   call read_mask_cube(maskfile, nx, ny, nchan, mask_cube, status)
   if (status.ne.0) stop 1

   allocate(clean_re_cube(nx,ny,nrm), clean_im_cube(nx,ny,nrm))
   allocate(resid_re_cube(nx,ny,nrm), resid_im_cube(nx,ny,nrm))
   allocate(restored_re_cube(nx,ny,nrm), restored_im_cube(nx,ny,nrm))

   call plan_fourier_interp(nrm, nrm, restore_plan_fwd, restore_plan_bwd)

   call build_mask_pattern_cache()

   n_pixels_done = 0
   !$omp parallel do collapse(2) schedule(dynamic) default(shared) private(ix,iy)
   do iy = 1, ny
      do ix = 1, nx
         call clean_one_pixel(ix, iy)
      enddo
   enddo
   !$omp end parallel do
   write(*,'(A,I0,A)') 'CLEANed ', n_pixels_done, ' pixels.'

   call destroy_fourier_interp_plan(restore_plan_fwd, restore_plan_bwd)

   call write_output_cube(ampfile, trim(outfile)//'.CLEAN.AMP.RMCUBE.FITS',&
   &trim(outfile)//'.CLEAN.PHA.RMCUBE.FITS', nx, ny, nrm, clean_re_cube,&
   &clean_im_cube, status)
   if (status.ne.0) stop 1
   call write_output_cube(ampfile, trim(outfile)//'.RESID.AMP.RMCUBE.FITS',&
   &trim(outfile)//'.RESID.PHA.RMCUBE.FITS', nx, ny, nrm, resid_re_cube,&
   &resid_im_cube, status)
   if (status.ne.0) stop 1
   call write_output_cube(ampfile, trim(outfile)//'.RESTORED.AMP.RMCUBE.FITS',&
   &trim(outfile)//'.RESTORED.PHA.RMCUBE.FITS', nx, ny, nrm,&
   &restored_re_cube, restored_im_cube, status)
   if (status.ne.0) stop 1

   write(*,'(A)') 'OK: rmclean_cubes complete.'

contains

   subroutine parse_args(status)
      integer, intent(out) :: status
      character(len=512) :: this_arg, cli_key, cli_val, cfgfile
      integer :: argc, iarg
      logical :: has_kv, have_cfgfile
      logical :: seen_ampfile, seen_phafile, seen_maskfile, seen_outfile

      status = 0
      ampfile = ' '
      phafile = ' '
      maskfile = ' '
      outfile = ' '
      niter = 500
      gain = 0.1_sp
      have_threshold = .false.
      threshold = 0.0_sp
      min_samples_per_fwhm = 2.0_sp
      refine_nsigma = 3.0_sp
      table_oversample = 20
      have_restore_fwhm = .false.
      restore_fwhm_override = 0.0_sp
      lsq_ref_report_mode_sel = lsq_ref_report_intrinsic
      have_lsq_ref_report_value = .false.
      lsq_ref_report_value = 0.0_sp
      lsq_ref_compute_mode = 'native'
      have_lsq_ref_compute_value = .false.
      lsq_ref_compute_value = 0.0_sp
      mask_pattern_cache_max = 4096
      have_cfgfile = .false.
      seen_ampfile = .false.
      seen_phafile = .false.
      seen_maskfile = .false.
      seen_outfile = .false.

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
            call apply_kv(trim(cli_key), trim(cli_val), seen_ampfile,&
            &seen_phafile, seen_maskfile, seen_outfile, status)
            if (status.ne.0) return
            iarg = iarg + 1
         endif
      enddo

      if (have_cfgfile) then
         call read_cfg_file(cfgfile, seen_ampfile, seen_phafile,&
         &seen_maskfile, seen_outfile, status)
         if (status.ne.0) return
      endif

      if (.not. seen_ampfile .or. .not. seen_phafile .or.&
      &.not. seen_maskfile .or. .not. seen_outfile) then
         write(*,*) 'ERROR: ampfile=, phafile=, maskfile=, and outfile='//&
         &' are all required'
         call print_usage()
         status = -1
         return
      endif
      if (.not. have_threshold) then
         write(*,*) 'ERROR: threshold= is required (CLEAN stopping flux,'//&
         &' same units as the dirty AMP cube)'
         status = -1
         return
      endif
   end subroutine parse_args

   subroutine apply_kv(key, val, seen_ampfile, seen_phafile, seen_maskfile,&
   &seen_outfile, status)
      character(len=*), intent(in) :: key, val
      logical, intent(inout) :: seen_ampfile, seen_phafile, seen_maskfile
      logical, intent(inout) :: seen_outfile
      integer, intent(out) :: status
      integer :: ios

      status = 0
      select case (key)
      case ('ampfile')
         ampfile = val
         seen_ampfile = .true.
      case ('phafile')
         phafile = val
         seen_phafile = .true.
      case ('maskfile')
         maskfile = val
         seen_maskfile = .true.
      case ('outfile')
         outfile = val
         seen_outfile = .true.
      case ('niter')
         read(val, *, iostat=ios) niter
         if (ios.ne.0) then
            write(*,*) 'ERROR: niter must be an integer'
            status = -1
            return
         endif
      case ('gain')
         read(val, *, iostat=ios) gain
         if (ios.ne.0) then
            write(*,*) 'ERROR: gain must be a number'
            status = -1
            return
         endif
      case ('threshold')
         read(val, *, iostat=ios) threshold
         if (ios.ne.0) then
            write(*,*) 'ERROR: threshold must be a number'
            status = -1
            return
         endif
         have_threshold = .true.
      case ('min_samples_per_fwhm')
         read(val, *, iostat=ios) min_samples_per_fwhm
         if (ios.ne.0 .or. min_samples_per_fwhm.lt.1.0_sp) then
            write(*,*) 'ERROR: min_samples_per_fwhm must be a number >= 1'
            status = -1
            return
         endif
      case ('refine_nsigma')
         read(val, *, iostat=ios) refine_nsigma
         if (ios.ne.0 .or. refine_nsigma.le.0.0_sp) then
            write(*,*) 'ERROR: refine_nsigma must be a positive number'
            status = -1
            return
         endif
      case ('table_oversample')
         read(val, *, iostat=ios) table_oversample
         if (ios.ne.0) then
            write(*,*) 'ERROR: table_oversample must be an integer'
            status = -1
            return
         endif
      case ('restore_fwhm')
         read(val, *, iostat=ios) restore_fwhm_override
         if (ios.ne.0) then
            write(*,*) 'ERROR: restore_fwhm must be a number'
            status = -1
            return
         endif
         have_restore_fwhm = .true.
      case ('lsq_ref_report_mode')
         select case (trim(val))
         case ('intrinsic')
            lsq_ref_report_mode_sel = lsq_ref_report_intrinsic
         case ('centroid')
            lsq_ref_report_mode_sel = lsq_ref_report_centroid
         case ('min')
            lsq_ref_report_mode_sel = lsq_ref_report_min
         case ('max')
            lsq_ref_report_mode_sel = lsq_ref_report_max
         case ('mid')
            lsq_ref_report_mode_sel = lsq_ref_report_mid
         case ('fixed')
            lsq_ref_report_mode_sel = lsq_ref_report_fixed
         case default
            write(*,*) 'ERROR: lsq_ref_report_mode must be one of'//&
            &' intrinsic|centroid|min|max|mid|fixed'
            status = -1
            return
         end select
      case ('lsq_ref_report_value')
         read(val, *, iostat=ios) lsq_ref_report_value
         if (ios.ne.0) then
            write(*,*) 'ERROR: lsq_ref_report_value must be a number'
            status = -1
            return
         endif
         have_lsq_ref_report_value = .true.
      case ('lsq_ref_compute_mode')
         select case (trim(val))
         case ('native', 'zero', 'mid', 'centroid', 'min', 'max', 'fixed')
            lsq_ref_compute_mode = trim(val)
         case default
            write(*,*) 'ERROR: lsq_ref_compute_mode must be one of'//&
            &' native|zero|mid|centroid|min|max|fixed'
            status = -1
            return
         end select
      case ('lsq_ref_compute_value')
         read(val, *, iostat=ios) lsq_ref_compute_value
         if (ios.ne.0) then
            write(*,*) 'ERROR: lsq_ref_compute_value must be a number'
            status = -1
            return
         endif
         have_lsq_ref_compute_value = .true.
      case ('mask_pattern_cache_max')
         read(val, *, iostat=ios) mask_pattern_cache_max
         if (ios.ne.0 .or. mask_pattern_cache_max.lt.0) then
            write(*,*) 'ERROR: mask_pattern_cache_max must be a'//&
            &' non-negative integer'
            status = -1
            return
         endif
      case default
         write(*,*) 'ERROR: unrecognised key "', key, '"'
         status = -1
         return
      end select
   end subroutine apply_kv

   subroutine read_cfg_file(cfgfile, seen_ampfile, seen_phafile,&
   &seen_maskfile, seen_outfile, status)
      character(len=*), intent(in) :: cfgfile
      logical, intent(inout) :: seen_ampfile, seen_phafile, seen_maskfile
      logical, intent(inout) :: seen_outfile
      integer, intent(out) :: status
      integer :: cunit, ios, line_no
      character(len=1024) :: line
      character(len=512) :: key, val
      logical :: has_kv

      status = 0
      open(newunit=cunit, file=trim(cfgfile), status='old', action='read',&
      &iostat=ios)
      if (ios.ne.0) then
         write(*,*) 'ERROR: cannot open --config file: ', trim(cfgfile)
         status = -1
         return
      endif
      line_no = 0
      do
         read(cunit, '(A)', iostat=ios) line
         if (ios.ne.0) exit
         line_no = line_no + 1
         call cfg_split_key_value(line, key, val, has_kv)
         if (.not. has_kv) cycle
         call apply_kv(trim(key), trim(val), seen_ampfile, seen_phafile,&
         &seen_maskfile, seen_outfile, status)
         if (status.ne.0) then
            write(*,*) '  (at ', trim(cfgfile), ' line ', line_no, ')'
            close(cunit)
            return
         endif
      enddo
      close(cunit)
   end subroutine read_cfg_file

   subroutine print_usage()
      write(*,'(A)') 'Usage: rmclean_cubes ampfile=<f> phafile=<f>'//&
      &' maskfile=<f> outfile=<base> threshold=<t>'
      write(*,'(A)') '    [niter=<n>] [gain=<g>] [min_samples_per_fwhm=<f>]'//&
      &' [refine_nsigma=<f>] [table_oversample=<n>] [restore_fwhm=<v>]'
      write(*,'(A)') '    [lsq_ref_report_mode=intrinsic|centroid|min|max'//&
      &'|mid|fixed] [lsq_ref_report_value=<v>]'
      write(*,'(A)') '    [lsq_ref_compute_mode=native|zero|centroid|min'//&
      &'|max|mid|fixed] [lsq_ref_compute_value=<v>]'
      write(*,'(A)') '    [mask_pattern_cache_max=<n>]'
      write(*,'(A)') '   or: rmclean_cubes --config <cfgfile>'
      write(*,'(A)') '   or: rmclean_cubes --help | -h'
      write(*,'(A)') ''
      write(*,'(A)') 'ampfile/phafile: the dirty AMP.RMCUBE.FITS/'//&
      &'PHA.RMCUBE.FITS pair rm_synthesis itself wrote.'
      write(*,'(A)') 'maskfile: the matching MASK.CUBE.FITS (with its own'//&
      &' CHANFREQ binary table).'
      write(*,'(A)') 'outfile: base name for the 6 output cubes'//&
      &' (<outfile>.CLEAN/.RESID/.RESTORED.AMP/PHA.RMCUBE.FITS).'
      write(*,'(A)') 'threshold: CLEAN stopping flux, same units as the'//&
      &' dirty AMP cube (required, no default -- see clean_complex''s own'//&
      &' doc comment in src/rmclean.f90).'
      write(*,'(A)') 'niter (default 500), gain (default 0.1): Hogbom'//&
      &' CLEAN parameters, passed straight to clean_complex.'
      write(*,'(A)') 'min_samples_per_fwhm (default 2, floor 1): Gate 0''s'//&
      &' own RM-resolution criterion -- the cube''s |CDELT3| must be <='//&
      &' fwhm/min_samples_per_fwhm (fwhm from compute_rmsf_fwhm_'//&
      &'multiband, lsq_ref-independent). Since T3''s model-based peak'//&
      &' refinement, resolution-level sampling is ALL the stored grid'//&
      &' needs -- the old carrier-driven (get_drm/max_offset) gate is'//&
      &' retired (planning ticket T3b).'
      write(*,'(A)') 'refine_nsigma (default 3): escalation threshold for'//&
      &' the tiered peak refinement (T3c) -- the cheap fixed-location'//&
      &' fit is accepted when its leftover misfit is within'//&
      &' refine_nsigma x the per-iteration data-driven noise estimate;'//&
      &' beyond that the full local search runs for that iteration.'
      write(*,'(A)') 'table_oversample (default 20): rmsf_table_t''s own'//&
      &' interpolation-table fineness (build_rmsf_offset_table''s own'//&
      &' oversample argument).'
      write(*,'(A)') 'restore_fwhm (optional): override the restoring'//&
      &' beam FWHM (rad/m^2); default derived from'//&
      &' compute_rmsf_fwhm_multiband.'
      write(*,'(A)') 'lsq_ref_report_mode (default intrinsic, i.e.'//&
      &' lsq_ref_report=0.0): where to report the derotated chi0/'//&
      &' restored phase -- a safe, independent post-processing choice.'
      write(*,'(A)') 'lsq_ref_compute_mode (default native, i.e. whatever'//&
      &' the cube''s own LSQREF header says): the reference this'//&
      &' program''s own RMSF table/CLEAN computation uses -- an exact,'//&
      &' free choice (derotate_to_lsq_ref), see this file''s own top'//&
      &' comment for why choosing a non-native value costs nothing but'//&
      &' also saves nothing once the RM grid already exists.'
      write(*,'(A)') 'mask_pattern_cache_max (default 4096): pixels'//&
      &' sharing the same valid-channel mask pattern share one'//&
      &' rmsf_table_t, built once during a serial pre-scan; past this'//&
      &' many DISTINCT patterns, additional patterns fall back to a'//&
      &' one-off table per pixel (safety valve, not a correctness'//&
      &' issue -- just loses the reuse benefit).'
   end subroutine print_usage

   subroutine cfg_split_key_value(raw_line, key, val, has_kv)
      character(len=*), intent(in) :: raw_line
      character(len=*), intent(out) :: key, val
      logical, intent(out) :: has_kv
      character(len=1024) :: line
      integer :: eqpos

      has_kv = .false.
      key = ' '
      val = ' '
      line = adjustl(raw_line)
      if (len_trim(line).eq.0) return
      if (line(1:1).eq.'#') return
      eqpos = index(line, '=')
      if (eqpos.lt.2) return
      key = adjustl(line(1:eqpos-1))
      val = adjustl(line(eqpos+1:))
      has_kv = .true.
   end subroutine cfg_split_key_value

   subroutine split_cli_kv(token, key, val, has_kv)
      character(len=*), intent(in) :: token
      character(len=*), intent(out) :: key, val
      logical, intent(out) :: has_kv
      integer :: eqpos

      has_kv = .false.
      key = ' '
      val = ' '
      eqpos = index(token, '=')
      if (eqpos.lt.2) return
      key = token(1:eqpos-1)
      val = token(eqpos+1:)
      has_kv = .true.
   end subroutine split_cli_kv

   subroutine read_chanfreq(filename, l_sq_out, band_id_out, nchan_out, status)
      !! Read maskfile's own CHANFREQ binary table (rm_synthesis.f90:
      !! 4762-4826: columns CHAN[1J]/FREQ[1D]/BAND[1J], row i lines up
      !! with MASK.CUBE.FITS's own axis-3 plane i) and invert FREQ (Hz)
      !! back to L_sq (m^2) exactly as rm_synthesis derived it forwards
      !! (freq_Hz = c_velocity(Mm/s)*1e6/sqrt(L_sq) there --> L_sq =
      !! (c_velocity*1e6/freq_Hz)**2, using the SAME c_velocity constant
      !! rmclean_mod itself already carries, src/rmclean.f90's own
      !! duplicated-not-imported convention).
      character(len=*), intent(in) :: filename
      real(sp), allocatable, intent(out) :: l_sq_out(:)
      integer, allocatable, intent(out) :: band_id_out(:)
      integer, intent(out) :: nchan_out
      integer, intent(out) :: status

      integer :: unit, blocksize, fitsstat, hdutype
      integer :: col_chan, col_freq, col_band
      logical :: anyflag
      integer, allocatable :: chan_col(:), band_col(:)
      real(dp), allocatable :: freq_col(:)
      real(dp), parameter :: c_velocity_dp = 299.792458_dp
      character(len=68) :: comment

      status = 0
      fitsstat = 0
      unit = 210
      call FTOPEN(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(filename)
         status = -1
         return
      endif

      call FTMNHD(unit, -1, 'CHANFREQ', 0, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: no CHANFREQ binary table in: ', trim(filename)
         call FTCLOS(unit, fitsstat)
         status = -1
         return
      endif
      call FTGKYJ(unit, 'NAXIS2', nchan_out, comment, fitsstat)
      if (fitsstat.ne.0 .or. nchan_out.lt.1) then
         write(*,*) 'ERROR: bad or missing NAXIS2 on CHANFREQ table in: ',&
         &trim(filename)
         call FTCLOS(unit, fitsstat)
         status = -1
         return
      endif

      allocate(chan_col(nchan_out), band_col(nchan_out), freq_col(nchan_out))
      call FTGCNO(unit, .false., 'CHAN', col_chan, fitsstat)
      call FTGCNO(unit, .false., 'FREQ', col_freq, fitsstat)
      call FTGCNO(unit, .false., 'BAND', col_band, fitsstat)
      call FTGCVJ(unit, col_chan, 1, 1, nchan_out, 0, chan_col, anyflag, fitsstat)
      call FTGCVD(unit, col_freq, 1, 1, nchan_out, 0.0d0, freq_col, anyflag, fitsstat)
      call FTGCVJ(unit, col_band, 1, 1, nchan_out, 0, band_col, anyflag, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to read CHANFREQ columns from: ', trim(filename)
         deallocate(chan_col, band_col, freq_col)
         call FTCLOS(unit, fitsstat)
         status = -1
         return
      endif
      call FTCLOS(unit, fitsstat)

      allocate(l_sq_out(nchan_out), band_id_out(nchan_out))
      l_sq_out = real((c_velocity_dp*1.0d6/freq_col)**2, sp)
      band_id_out = band_col
      deallocate(chan_col, band_col, freq_col)
   end subroutine read_chanfreq

   subroutine derive_band_segmentation(band_id, nchan, band_offset_out,&
   &band_nz_out, n_bands_out, status)
      !! rm_synthesis's own band_offset/band_nz convention (rm_synthesis.
      !! f90's CHANFREQ-writing comment): band k occupies a CONTIGUOUS
      !! run of channel indices, reference band (BAND=0) first at
      !! offset 0, every other band appended after in increasing BAND
      !! order. Verified here (not assumed) by checking BAND is
      !! non-decreasing along the channel axis -- refuses loudly if not,
      !! rather than silently mis-segmenting compute_rmsf_fwhm_
      !! multiband's own per-band span calculation.
      integer, intent(in) :: band_id(:), nchan
      integer, allocatable, intent(out) :: band_offset_out(:), band_nz_out(:)
      integer, intent(out) :: n_bands_out
      integer, intent(out) :: status
      integer :: k

      status = 0
      do k = 2, nchan
         if (band_id(k).lt.band_id(k-1)) then
            write(*,*) 'ERROR: CHANFREQ BAND column is not non-decreasing'//&
            &' -- expected rm_synthesis''s own contiguous-per-band'//&
            &' channel ordering'
            status = -1
            return
         endif
      enddo
      n_bands_out = band_id(nchan) + 1
      allocate(band_offset_out(n_bands_out), band_nz_out(n_bands_out))
      band_nz_out = 0
      do k = 1, nchan
         band_nz_out(band_id(k)+1) = band_nz_out(band_id(k)+1) + 1
      enddo
      band_offset_out(1) = 0
      do k = 2, n_bands_out
         band_offset_out(k) = band_offset_out(k-1) + band_nz_out(k-1)
      enddo
   end subroutine derive_band_segmentation

   subroutine read_amp_pha_geometry(ampfile, phafile, nx_out, ny_out,&
   &nrm_out, cdelt3_out, crval3_out, status)
      character(len=*), intent(in) :: ampfile, phafile
      integer, intent(out) :: nx_out, ny_out, nrm_out
      real(dp), intent(out) :: cdelt3_out, crval3_out
      integer, intent(out) :: status
      integer :: unit, blocksize, fitsstat, naxis
      integer :: nx2, ny2, nrm2
      real(dp) :: cdelt3_2, crval3_2
      character(len=68) :: comment

      status = 0
      fitsstat = 0
      unit = 211
      call FTOPEN(unit, trim(ampfile), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(ampfile)
         status = -1
         return
      endif
      call FTGKYJ(unit, 'NAXIS', naxis, comment, fitsstat)
      if (fitsstat.ne.0 .or. naxis.ne.3) then
         write(*,*) 'ERROR: expected NAXIS=3 in: ', trim(ampfile)
         call FTCLOS(unit, fitsstat)
         status = -1
         return
      endif
      call FTGKYJ(unit, 'NAXIS1', nx_out, comment, fitsstat)
      call FTGKYJ(unit, 'NAXIS2', ny_out, comment, fitsstat)
      call FTGKYJ(unit, 'NAXIS3', nrm_out, comment, fitsstat)
      call FTGKYD(unit, 'CDELT3', cdelt3_out, comment, fitsstat)
      call FTGKYD(unit, 'CRVAL3', crval3_out, comment, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: missing NAXIS1/2/3 or CDELT3/CRVAL3 in: ',&
         &trim(ampfile)
         call FTCLOS(unit, fitsstat)
         status = -1
         return
      endif
      call FTCLOS(unit, fitsstat)

      fitsstat = 0
      unit = 212
      call FTOPEN(unit, trim(phafile), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(phafile)
         status = -1
         return
      endif
      call FTGKYJ(unit, 'NAXIS1', nx2, comment, fitsstat)
      call FTGKYJ(unit, 'NAXIS2', ny2, comment, fitsstat)
      call FTGKYJ(unit, 'NAXIS3', nrm2, comment, fitsstat)
      call FTGKYD(unit, 'CDELT3', cdelt3_2, comment, fitsstat)
      call FTGKYD(unit, 'CRVAL3', crval3_2, comment, fitsstat)
      call FTCLOS(unit, fitsstat)
      if (fitsstat.ne.0 .or. nx2.ne.nx_out .or. ny2.ne.ny_out .or.&
      &nrm2.ne.nrm_out .or. cdelt3_2.ne.cdelt3_out .or.&
      &crval3_2.ne.crval3_out) then
         write(*,*) 'ERROR: ', trim(phafile), ' geometry does not match ',&
         &trim(ampfile)
         status = -1
         return
      endif
   end subroutine read_amp_pha_geometry

   subroutine read_lsqref_keyword(filename, lsq_ref_out, found)
      !! Reads the LSQREF keyword rm_synthesis now writes onto its own
      !! AMP/PHA cubes (rm_synthesis.f90's own cfg%lsq_ref_mode/
      !! lsq_ref_fixed_value, via compute_lsq_ref). found=.false. (with
      !! lsq_ref_out left undefined) signals the keyword is absent --
      !! the caller, not this subroutine, decides the correct fallback
      !! (0.0, this project's historical convention -- see this file's
      !! own top comment), so the caller can print an explicit warning
      !! rather than this routine silently assuming a value.
      character(len=*), intent(in) :: filename
      real(sp), intent(out) :: lsq_ref_out
      logical, intent(out) :: found
      integer :: unit, blocksize, fitsstat
      real(dp) :: lsq_ref_dp
      character(len=68) :: comment

      found = .false.
      fitsstat = 0
      unit = 219
      call FTOPEN(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) return
      fitsstat = 0
      call FTGKYD(unit, 'LSQREF', lsq_ref_dp, comment, fitsstat)
      if (fitsstat.eq.0) then
         found = .true.
         lsq_ref_out = real(lsq_ref_dp, sp)
      endif
      fitsstat = 0
      call FTCLOS(unit, fitsstat)
   end subroutine read_lsqref_keyword

   function lsq_ref_compute_mode_to_int(mode_str) result(mode_int)
      !! Maps this program's own lsq_ref_compute_mode string (its 'native'
      !! value has no rmclean_mod counterpart -- handled by the caller
      !! before this function is ever invoked) onto rmclean_mod's own
      !! get_lsq_ref_compute mode constants.
      character(len=*), intent(in) :: mode_str
      integer :: mode_int
      select case (trim(mode_str))
      case ('zero')
         mode_int = lsq_ref_compute_intrinsic
      case ('mid')
         mode_int = lsq_ref_compute_mid
      case ('centroid')
         mode_int = lsq_ref_compute_centroid
      case ('min')
         mode_int = lsq_ref_compute_min
      case ('max')
         mode_int = lsq_ref_compute_max
      case ('fixed')
         mode_int = lsq_ref_compute_fixed
      case default
         write(*,*) 'FATAL: lsq_ref_compute_mode_to_int: unrecognized mode: ',&
         &trim(mode_str)
         stop 1
      end select
   end function lsq_ref_compute_mode_to_int

   subroutine read_cube(ampfile, phafile, nx_in, ny_in, nrm_in, re_out,&
   &im_out, status)
      !! Whole-cube read (Stage A -- see this file's own top comment):
      !! amp/pha -> re/im, matching rm_synthesis's own p_tile_arr=sqrt(re^2
      !! +im^2)/phi_tile_arr=atan2(im,re) forward convention exactly
      !! (rm_synthesis_mod.f90:1462-1465, its own output_mode=0/
      !! ap_angle_mode=0 branch -- the default this program assumes,
      !! since ampfile/phafile are always the AMP/PHA pair by this
      !! program's own required config keys): re=amp*cos(pha),
      !! im=amp*sin(pha).
      character(len=*), intent(in) :: ampfile, phafile
      integer, intent(in) :: nx_in, ny_in, nrm_in
      real(sp), allocatable, intent(out) :: re_out(:,:,:), im_out(:,:,:)
      integer, intent(out) :: status
      integer :: unit, blocksize, fitsstat
      logical :: anyflag
      real(sp), allocatable :: amp_cube(:,:,:), pha_cube(:,:,:)

      status = 0
      allocate(amp_cube(nx_in,ny_in,nrm_in), pha_cube(nx_in,ny_in,nrm_in))

      fitsstat = 0
      unit = 213
      call FTOPEN(unit, trim(ampfile), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(ampfile)
         status = -1
         return
      endif
      call FTGPVE(unit, 1, 1, nx_in*ny_in*nrm_in, 0.0_sp, amp_cube, anyflag, fitsstat)
      call FTCLOS(unit, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to read data from: ', trim(ampfile)
         status = -1
         return
      endif

      fitsstat = 0
      unit = 214
      call FTOPEN(unit, trim(phafile), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(phafile)
         status = -1
         return
      endif
      call FTGPVE(unit, 1, 1, nx_in*ny_in*nrm_in, 0.0_sp, pha_cube, anyflag, fitsstat)
      call FTCLOS(unit, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to read data from: ', trim(phafile)
         status = -1
         return
      endif

      allocate(re_out(nx_in,ny_in,nrm_in), im_out(nx_in,ny_in,nrm_in))
      re_out = amp_cube*cos(pha_cube)
      im_out = amp_cube*sin(pha_cube)
      deallocate(amp_cube, pha_cube)
   end subroutine read_cube

   subroutine read_mask_cube(filename, nx_in, ny_in, nchan_in, mask_out, status)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nx_in, ny_in, nchan_in
      integer(kind=1), allocatable, intent(out) :: mask_out(:,:,:)
      integer, intent(out) :: status
      integer :: unit, blocksize, fitsstat
      logical :: anyflag

      status = 0
      allocate(mask_out(nx_in,ny_in,nchan_in))
      fitsstat = 0
      unit = 215
      call FTOPEN(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(filename)
         status = -1
         return
      endif
      call FTGPVB(unit, 1, 1, nx_in*ny_in*nchan_in, 0_1, mask_out, anyflag, fitsstat)
      call FTCLOS(unit, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to read data from: ', trim(filename)
         status = -1
         return
      endif
   end subroutine read_mask_cube

   subroutine write_output_cube(template_file, amp_outname, pha_outname,&
   &nx_in, ny_in, nrm_in, re_in, im_in, status)
      !! re/im -> amp/pha (inverse of read_cube's own forward convention),
      !! header copied verbatim from template_file (always ampfile: same
      !! WCS on every axis, since this program never resamples anything --
      !! Gate 0 validates the existing RM grid rather than changing it).
      character(len=*), intent(in) :: template_file, amp_outname, pha_outname
      integer, intent(in) :: nx_in, ny_in, nrm_in
      real(sp), intent(in) :: re_in(nx_in,ny_in,nrm_in), im_in(nx_in,ny_in,nrm_in)
      integer, intent(out) :: status
      integer :: src_unit, amp_unit, pha_unit, fitsstat, blocksize
      integer :: naxes_out(3)
      logical :: simple, extend
      real(sp), allocatable :: amp_cube(:,:,:), pha_cube(:,:,:)

      status = 0
      allocate(amp_cube(nx_in,ny_in,nrm_in), pha_cube(nx_in,ny_in,nrm_in))
      amp_cube = sqrt(re_in**2 + im_in**2)
      pha_cube = atan2(im_in, re_in)

      naxes_out(1) = nx_in
      naxes_out(2) = ny_in
      naxes_out(3) = nrm_in
      simple = .true.
      extend = .false.

      fitsstat = 0
      src_unit = 216
      call FTOPEN(src_unit, trim(template_file), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(template_file)
         status = -1
         deallocate(amp_cube, pha_cube)
         return
      endif

      ! FTINIT fails (nonzero fitsstat) if amp_outname already exists on
      ! disk -- checked and bailed out on IMMEDIATELY, before any further
      ! CFITSIO call on amp_unit: every call after a failed FTINIT
      ! operates on a unit CFITSIO never actually set up, which is not a
      ! clean no-op but undefined behaviour (confirmed directly: this
      ! previously crashed with a SIGSEGV inside CFITSIO when an output
      ! file from an earlier run was left on disk -- tests/run_tests.sh's
      ! own section 28 must therefore also clean up its own rmc_* outputs
      ! before each run, which it now does).
      fitsstat = 0
      amp_unit = 217
      call FTINIT(amp_unit, trim(amp_outname), blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot create (already exists?): ', trim(amp_outname)
         status = -1
         call FTCLOS(src_unit, fitsstat)
         deallocate(amp_cube, pha_cube)
         return
      endif
      call FTPHPR(amp_unit, simple, -32, 3, naxes_out, 0, 1, extend, fitsstat)
      call copy_generic_header_rmclean(src_unit, amp_unit, status)
      call FTPPRE(amp_unit, 1, 1, nx_in*ny_in*nrm_in, amp_cube, fitsstat)
      call FTCLOS(amp_unit, fitsstat)
      if (fitsstat.ne.0 .or. status.ne.0) then
         write(*,*) 'ERROR: failed to write: ', trim(amp_outname)
         status = -1
         call FTCLOS(src_unit, fitsstat)
         deallocate(amp_cube, pha_cube)
         return
      endif

      fitsstat = 0
      pha_unit = 218
      call FTINIT(pha_unit, trim(pha_outname), blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot create (already exists?): ', trim(pha_outname)
         status = -1
         call FTCLOS(src_unit, fitsstat)
         deallocate(amp_cube, pha_cube)
         return
      endif
      call FTPHPR(pha_unit, simple, -32, 3, naxes_out, 0, 1, extend, fitsstat)
      call copy_generic_header_rmclean(src_unit, pha_unit, status)
      call FTPPRE(pha_unit, 1, 1, nx_in*ny_in*nrm_in, pha_cube, fitsstat)
      call FTCLOS(pha_unit, fitsstat)
      if (fitsstat.ne.0 .or. status.ne.0) then
         write(*,*) 'ERROR: failed to write: ', trim(pha_outname)
         status = -1
      endif

      call FTCLOS(src_unit, fitsstat)
      deallocate(amp_cube, pha_cube)
   end subroutine write_output_cube

   subroutine copy_generic_header_rmclean(src_unit, dst_unit, status)
      !! Verbatim adaptation of convolve_cubes.f90's own
      !! copy_generic_header_convolve: input and output share IDENTICAL
      !! axis layout (this program resamples nothing, Gate 0 validates
      !! rather than changes the RM grid), so every header card copies
      !! through except the structural keywords FTPHPR already wrote.
      integer, intent(in) :: src_unit, dst_unit
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
         select case (trim(key))
         case ('SIMPLE', 'BITPIX', 'NAXIS', 'EXTEND', 'PCOUNT', 'GCOUNT', 'END')
            cycle
         end select
         if (is_naxis_keyword(key)) cycle
         fitsstat = 0
         call FTPREC(dst_unit, card, fitsstat)
      enddo
   end subroutine copy_generic_header_rmclean

   logical function is_naxis_keyword(key)
      character(len=8), intent(in) :: key
      character(len=8) :: k
      integer :: klen
      k = adjustl(key)
      klen = len_trim(k)
      is_naxis_keyword = .false.
      if (klen.ge.6 .and. klen.le.7 .and. k(1:5).eq.'NAXIS') then
         is_naxis_keyword = .true.
      endif
   end function is_naxis_keyword

   function fnv1a_hash(bytes, n) result(h)
      !! Standard 64-bit FNV-1a over a byte pattern (here: one pixel's
      !! own valid-channel mask row, mask_cube(ix,iy,:)) -- used as the
      !! bucket key for table_cache_entry_t below. Collisions are
      !! possible (any hash is), so every lookup below ALSO does a full
      !! byte-for-byte pattern compare on a hash match -- never trusts
      !! the hash alone.
      integer(kind=1), intent(in) :: bytes(n)
      integer, intent(in) :: n
      integer(kind=8) :: h
      integer(kind=8), parameter :: fnv_offset = -3750763034362895579_8
      integer(kind=8), parameter :: fnv_prime = 1099511628211_8
      integer :: i
      h = fnv_offset
      do i = 1, n
         h = ieor(h, int(bytes(i), 8))
         h = h * fnv_prime
      end do
   end function fnv1a_hash

   subroutine cache_lookup_readonly(pattern, entry_idx)
      !! Open-addressing (linear probing) lookup into cache_buckets/
      !! cache_entries -- PURE READ access, safe to call concurrently
      !! from every OMP thread in the main CLEAN loop, since
      !! build_mask_pattern_cache below has already fully populated the
      !! cache serially before that parallel region ever starts (no
      !! insertion ever happens in here). entry_idx=0 means "not
      !! cached" -- either this exact pattern never occurred during the
      !! pre-scan (impossible, since the pre-scan visits the same
      !! mask_cube this is called against) or the cache was already at
      !! mask_pattern_cache_max when this pattern was first seen --
      !! either way, the caller's own correct response is to build a
      !! one-off throwaway table (see clean_one_pixel).
      integer(kind=1), intent(in) :: pattern(:)
      integer, intent(out) :: entry_idx
      integer(kind=8) :: h, probe
      integer :: tries, slot

      h = fnv1a_hash(pattern, size(pattern))
      probe = modulo(h, int(n_cache_buckets, 8))
      entry_idx = 0
      do tries = 1, n_cache_buckets
         slot = int(probe)
         if (cache_buckets(slot).eq.0) return
         if (cache_entries(cache_buckets(slot))%hash.eq.h) then
            if (size(cache_entries(cache_buckets(slot))%pattern).eq.size(pattern)) then
               if (all(cache_entries(cache_buckets(slot))%pattern.eq.pattern)) then
                  entry_idx = cache_buckets(slot)
                  return
               endif
            endif
         endif
         probe = modulo(probe+1, int(n_cache_buckets, 8))
      enddo
   end subroutine cache_lookup_readonly

   subroutine build_mask_pattern_cache()
      !! Serial pre-scan (planning/RMCLEAN_INTEGRATION_PLAN.md decision
      !! 10), called ONCE before the parallel per-pixel CLEAN loop:
      !! canonicalize each pixel's own valid-channel bit pattern
      !! (mask_cube(ix,iy,:) itself, already a 0/1 byte array -- no
      !! further packing needed), hash it, and build ONE rmsf_table_t
      !! per DISTINCT pattern, up to mask_pattern_cache_max. This is
      !! the ONLY place cache_entries/cache_buckets/n_cache_entries are
      !! ever written -- by construction single-threaded, so no lock is
      !! needed here either.
      integer :: ix_l, iy_l, entry_idx
      integer(kind=8) :: h, probe
      integer :: tries, slot
      integer :: nvalid_l
      integer, allocatable :: valid_idx_l(:)
      real(sp), allocatable :: l_sq_valid_l(:)
      integer :: n_distinct_patterns, n_overflow_pixels

      n_cache_buckets = max(16, 4*mask_pattern_cache_max)
      allocate(cache_buckets(0:n_cache_buckets-1))
      cache_buckets = 0
      allocate(cache_entries(max(1, mask_pattern_cache_max)))
      n_cache_entries = 0
      allocate(valid_idx_l(nchan), l_sq_valid_l(nchan))
      n_distinct_patterns = 0
      n_overflow_pixels = 0

      do iy_l = 1, ny
         do ix_l = 1, nx
            if (all(mask_cube(ix_l,iy_l,:).eq.0_1)) cycle

            h = fnv1a_hash(mask_cube(ix_l,iy_l,:), nchan)
            probe = modulo(h, int(n_cache_buckets, 8))
            entry_idx = -1
            do tries = 1, n_cache_buckets
               slot = int(probe)
               if (cache_buckets(slot).eq.0) exit
               if (cache_entries(cache_buckets(slot))%hash.eq.h) then
                  if (all(cache_entries(cache_buckets(slot))%pattern.eq.&
                  &mask_cube(ix_l,iy_l,:))) then
                     entry_idx = cache_buckets(slot)
                     exit
                  endif
               endif
               probe = modulo(probe+1, int(n_cache_buckets, 8))
            enddo
            if (entry_idx.ne.-1) cycle ! already cached (0 means bucket array full, handled below)

            if (n_cache_entries.ge.mask_pattern_cache_max) then
               n_overflow_pixels = n_overflow_pixels + 1
               cycle ! safety valve: clean_one_pixel builds a throwaway table
            endif

            n_cache_entries = n_cache_entries + 1
            n_distinct_patterns = n_distinct_patterns + 1
            cache_buckets(slot) = n_cache_entries
            cache_entries(n_cache_entries)%hash = h
            allocate(cache_entries(n_cache_entries)%pattern(nchan))
            cache_entries(n_cache_entries)%pattern = mask_cube(ix_l,iy_l,:)

            nvalid_l = 0
            do k = 1, nchan
               if (mask_cube(ix_l,iy_l,k).ne.0) then
                  nvalid_l = nvalid_l + 1
                  valid_idx_l(nvalid_l) = k
               endif
            enddo
            l_sq_valid_l(1:nvalid_l) = l_sq(valid_idx_l(1:nvalid_l))
            call build_rmsf_offset_table(l_sq_valid_l(1:nvalid_l), nvalid_l,&
            &lsq_ref_compute, rm_samp(nrm)-rm_samp(1), real(cdelt3_amp, sp),&
            &table_oversample, cache_entries(n_cache_entries)%table)
         enddo
      enddo

      deallocate(valid_idx_l, l_sq_valid_l)
      write(*,'(A,I0,A,I0,A)') 'Mask-pattern cache: ', n_distinct_patterns,&
      &' distinct pattern(s) cached'
      if (n_overflow_pixels.gt.0) then
         write(*,'(A,I0,A)') '  (', n_overflow_pixels, ' pixel(s) past'//&
         &' mask_pattern_cache_max fall back to a one-off table each)'
      endif
   end subroutine build_mask_pattern_cache

   subroutine clean_one_pixel(ix_p, iy_p)
      !! One pixel's full CLEAN+restore, called from the main program's
      !! own `!$omp parallel do` over (ix,iy). Every array declared here
      !! is a genuine LOCAL (automatic per-call) variable -- Fortran
      !! gives each concurrent call its own independent storage for
      !! these, with no `save` and no module-level state written here,
      !! so this subroutine is thread-safe purely by construction, no
      !! explicit `private()` bookkeeping needed for them. Everything
      !! this subroutine reads via host association (re_cube, mask_cube,
      !! rm_samp, l_sq, cache_entries/cache_buckets, restore_plan_fwd/
      !! bwd, ...) is READ-ONLY here; every array it WRITES (the 6
      !! output cubes) is written only at (ix_p,iy_p,:), a disjoint
      !! location per call -- no race either way.
      integer, intent(in) :: ix_p, iy_p
      integer :: nvalid_p, n_iter_used_p, k_p, entry_idx_p
      integer, allocatable :: valid_idx_p(:)
      real(sp), allocatable :: l_sq_valid_p(:)
      real(sp) :: dirty_re_p(nrm), dirty_im_p(nrm)
      real(sp) :: comp_re_p(nrm), comp_im_p(nrm)
      real(sp) :: resid_re_p(nrm), resid_im_p(nrm)
      real(sp) :: comp_rm_refined_p(nrm)
      real(sp) :: restored_re_p(nrm), restored_im_p(nrm)
      type(rmsf_table_t) :: throwaway_table
      logical :: used_throwaway

      nvalid_p = 0
      allocate(valid_idx_p(nchan), l_sq_valid_p(nchan))
      do k_p = 1, nchan
         if (mask_cube(ix_p,iy_p,k_p).ne.0) then
            nvalid_p = nvalid_p + 1
            valid_idx_p(nvalid_p) = k_p
         endif
      enddo

      if (nvalid_p.lt.1) then
         ! No valid channel ever contributed at this pixel -- the dirty
         ! spectrum itself is already all-NaN (rm_synthesis's own
         ! no-valid-data convention), nothing to CLEAN. Pass NaN
         ! straight through to every output, matching rm_synthesis's own
         ! bad-pixel policy rather than inventing a table for an empty
         ! channel set.
         clean_re_cube(ix_p,iy_p,:) = re_cube(ix_p,iy_p,:)
         clean_im_cube(ix_p,iy_p,:) = im_cube(ix_p,iy_p,:)
         resid_re_cube(ix_p,iy_p,:) = re_cube(ix_p,iy_p,:)
         resid_im_cube(ix_p,iy_p,:) = im_cube(ix_p,iy_p,:)
         restored_re_cube(ix_p,iy_p,:) = re_cube(ix_p,iy_p,:)
         restored_im_cube(ix_p,iy_p,:) = im_cube(ix_p,iy_p,:)
         deallocate(valid_idx_p, l_sq_valid_p)
         return
      endif

      ! l_sq_valid_p is filled unconditionally now (not just for the
      ! throwaway-table branch): clean_complex's own T3 refinement step
      ! (refine_peak_matched_filter) needs the valid-channel l_sq list
      ! regardless of which table (cached or throwaway) was used to
      ! build the SUBTRACTION beam.
      l_sq_valid_p(1:nvalid_p) = l_sq(valid_idx_p(1:nvalid_p))

      call cache_lookup_readonly(mask_cube(ix_p,iy_p,:), entry_idx_p)
      used_throwaway = (entry_idx_p.eq.0)
      if (used_throwaway) then
         call build_rmsf_offset_table(l_sq_valid_p(1:nvalid_p), nvalid_p,&
         &lsq_ref_compute, rm_samp(nrm)-rm_samp(1), real(cdelt3_amp, sp),&
         &table_oversample, throwaway_table)
      endif

      dirty_re_p = re_cube(ix_p,iy_p,:)
      dirty_im_p = im_cube(ix_p,iy_p,:)

      if (lsq_ref_compute.ne.lsq_ref_native) then
         ! Exact, lossless phase rotation of the ALREADY-SAMPLED dirty
         ! spectrum from the reference it was actually built at
         ! (lsq_ref_native) to this program's own chosen compute
         ! reference -- see this file's own top comment. Reusing
         ! derotate_to_lsq_ref's own (from, to) argument order: its own
         ! parameter names are lsq_ref_compute/lsq_ref_report, but the
         ! routine is fully general (re-express a spectrum built at ANY
         ! reference at any other), not tied to those two particular use
         ! cases.
         call derotate_to_lsq_ref(rm_samp, nrm, lsq_ref_native,&
         &lsq_ref_compute, dirty_re_p, dirty_im_p, dirty_re_p, dirty_im_p)
      endif

      if (used_throwaway) then
         call clean_complex(l_sq_valid_p(1:nvalid_p), nvalid_p,&
         &lsq_ref_compute, rm_samp, nrm, dirty_re_p, dirty_im_p,&
         &throwaway_table, niter, gain, threshold, comp_re_p, comp_im_p,&
         &resid_re_p, resid_im_p, n_iter_used_p, comp_rm_refined_p,&
         &nsigma_refine=refine_nsigma)
         call destroy_rmsf_offset_table(throwaway_table)
      else
         call clean_complex(l_sq_valid_p(1:nvalid_p), nvalid_p,&
         &lsq_ref_compute, rm_samp, nrm, dirty_re_p, dirty_im_p,&
         &cache_entries(entry_idx_p)%table, niter, gain, threshold,&
         &comp_re_p, comp_im_p, resid_re_p, resid_im_p, n_iter_used_p,&
         &comp_rm_refined_p, nsigma_refine=refine_nsigma)
      endif

      call restore_clean(rm_samp, nrm, comp_re_p, comp_im_p, resid_re_p,&
      &resid_im_p, fwhm_rm, restore_plan_fwd, restore_plan_bwd,&
      &restored_re_p, restored_im_p)

      if (lsq_ref_report.ne.lsq_ref_compute) then
         call derotate_to_lsq_ref(rm_samp, nrm, lsq_ref_compute,&
         &lsq_ref_report, comp_re_p, comp_im_p, comp_re_p, comp_im_p)
         call derotate_to_lsq_ref(rm_samp, nrm, lsq_ref_compute,&
         &lsq_ref_report, resid_re_p, resid_im_p, resid_re_p, resid_im_p)
         call derotate_to_lsq_ref(rm_samp, nrm, lsq_ref_compute,&
         &lsq_ref_report, restored_re_p, restored_im_p, restored_re_p,&
         &restored_im_p)
      endif

      clean_re_cube(ix_p,iy_p,:) = comp_re_p
      clean_im_cube(ix_p,iy_p,:) = comp_im_p
      resid_re_cube(ix_p,iy_p,:) = resid_re_p
      resid_im_cube(ix_p,iy_p,:) = resid_im_p
      restored_re_cube(ix_p,iy_p,:) = restored_re_p
      restored_im_cube(ix_p,iy_p,:) = restored_im_p

      deallocate(valid_idx_p, l_sq_valid_p)
      !$omp atomic
      n_pixels_done = n_pixels_done + 1
   end subroutine clean_one_pixel

end program rmclean_cubes
