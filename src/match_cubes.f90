! match_cubes -- consolidates reproject_cubes (src/reproject_cubes.f90) and
! convolve_cubes (src/convolve_cubes.f90) into one tool that can run either
! stage alone, or both chained THROUGH MEMORY with no intermediate FITS file
! written to disk. Built for real multi-band data at real scale (200GB+
! cubes): reproject_cubes and convolve_cubes run back-to-back already
! produce a correct pipeline, but the intermediate _REPROJ.FITS file is
! written in full and immediately read back in full for no reason other
! than being two separate programs -- doubling disk I/O and disk space for
! an artifact nobody wants.
!
! The two existing standalone tools are NOT touched by this file (a
! deliberate choice, confirmed with the user): reproject_cubes.f90 and
! convolve_cubes.f90 remain independent, already-tested tools for
! anyone who wants just one stage without this one. This file therefore
! duplicates (adapts, not `use`s) the subroutines it needs from both,
! rather than extracting a shared module -- an accepted maintenance
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
!    [beamfiles=<spec1>[,<spec2>...]] [badchan_file=<file1>[,<file2>...]]
!    [target_bmaj=<arcsec> target_bmin=<arcsec> target_bpa=<deg>]
!    [max_common_bmaj=<arcsec>] [mem_frac_ram=<fraction>] [outsuffix=<suffix>]
!    [npts=<n>] [khachiyan_tol=<tol>] [manifest=<path>]
!    or: match_cubes --config <cfgfile>
!    or: match_cubes --help | -h
! Full usage text in print_usage below (shared by --help and the
! argument-error path, same convention as reproject_cubes.f90/
! convolve_cubes.f90).
program match_cubes
   use, intrinsic :: iso_fortran_env, only: dp => real64, int8, int32, int64
   use, intrinsic :: iso_c_binding, only: c_int, c_long, c_ptr, c_funptr,&
   &c_null_ptr, c_funloc, c_loc, c_f_pointer
   use logging_mod
   use fitsio_unit_mod
   ! Per-channel beam readers, shared with match_cubes/convolve_cubes and
   ! the reprojection beam-sampling check (one copy, src/beam_io_mod.f90).
   use beam_io_mod, only: read_beams
   use wcs_match_mod
   use thread_safety_mod, only: check_and_clamp_omp_threads
   ! ast_interp_mod: shared with reproject_cubes.f90 (docs/dev/
   ! REPROJECT_INTERP_OPTIONS_PLAN.md) -- the one deliberate exception
   ! to this file's own top-of-file "duplicates, does not `use`, the
   ! standalone tools' own logic" policy: ast_interp_mod has no FITS-
   ! reading side effects and touches neither existing tool's own
   ! tested code, and a second independent copy of this resolution
   ! logic is exactly the kind of drift that produced AST's own
   ! SincSinc/SincCos cutoff-constant mismatch (REPROJECT_INTERP_
   ! OPTIONS_PLAN.md section 3).
   use ast_interp_mod, only: resolve_interp, ast__ukern1, interp_exact
   use exact_overlap_mod, only: exact_weights_t, exact_init, exact_build_rows,&
   &exact_apply_rows, exact_free, exact_weights_bytes, exact_estimate_bytes, mib_str
   use target_grid_mod, only: resolve_target_grid, write_target_wcs_cards,&
   &check_output_grid_validity,&
   &build_target_skymap, check_overlap_and_size_grid, check_disk_space,&
   &output_cell_size_arcsec, check_beam_sampling, common_beam_sampling,&
   &fits_reference_cell_arcsec, cell_equal_rtol
   implicit none

   ! --- Logging & timing (planning-doc ticket) -- see convolve_cubes.f90
   ! /logging_mod.f90's own header comments. Stage names: startup,
   ! block_read, block_convolve/block_resample, block_write,
   ! block_write_join, finalize.
   character(len=16) :: log_level
   logical :: timing_enabled
   character(len=272) :: log_output_file

   interface
      function c_pthread_create(thread, attr, start_routine, arg)&
      &bind(C, name="pthread_create") result(rc)
         import :: c_int, c_long, c_ptr, c_funptr
         integer(c_long) :: thread
         type(c_ptr), value :: attr
         type(c_funptr), value :: start_routine
         type(c_ptr), value :: arg
         integer(c_int) :: rc
      end function c_pthread_create

      function c_pthread_join(thread, retval)&
      &bind(C, name="pthread_join") result(rc)
         import :: c_int, c_long, c_ptr
         integer(c_long), value :: thread
         type(c_ptr), value :: retval
         integer(c_int) :: rc
      end function c_pthread_join
   end interface
   ! AST_PAR (the vendor Fortran constants file, /usr/include/AST_PAR) is
   ! fixed-form Fortran 77 and cannot be `include`d into a free-form .f90
   ! file directly -- same issue documented in reproject_cubes.f90's own
   ! comment. Only the handful of symbols used are declared
   ! directly instead, matching AST_PAR's own declared types.
   external :: ast_null, lanczos_kernel
   integer, parameter :: ast__null = 0
   integer, parameter :: ast__szchr = 200
   integer, parameter :: ast__base = 0
   integer, parameter :: ast__current = -1
   integer, external :: ast_fitschan, ast_read, ast_geti
   integer, external :: ast_getmapping, ast_simplify, ast_getframe
   integer, external :: ast_pickaxes, ast_cmpmap, ast_convert
   integer, external :: ast_resampler
   logical, external :: ast_isaframeset, ast_isaskyframe
   character(len=ast__szchr), external :: ast_getc

   integer, parameter :: max_axes = 10
   integer, parameter :: max_inputs = 50
   integer, parameter :: max_channels = 20000
   ! Hard ceiling on elements-per-block -- see convolve_cubes.f90's own
   ! identical constant/comment (CFITSIO's Fortran wrapper takes a
   ! default-INTEGER element count; this clamp keeps nx*ny*block_planes
   ! safely under 2^31-1 for ANY image size/mem_frac_ram combination).
   integer(kind=8), parameter :: max_elements_per_block = 2000000000_8

   ! --- io_overlap: background-thread block write (planning-doc ticket)
   ! --- same design as convolve_cubes.f90's own (see write_convolved_
   ! file's own comment there for the full single-writer-at-a-time
   ! rationale); block_write_job_t is self-contained for the same reason
   ! (process_one_file_restricted's own out_unit/naxes/nx_out/ny_out are
   ! LOCAL to that subroutine, called once per input file, not
   ! program-level state).
   type :: block_write_job_t
      ! file_path/datastart (not out_unit): raw Fortran stream I/O at a
      ! computed byte offset (T19, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.
      ! md), not a live CFITSIO handle -- see process_one_file_restricted's
      ! own comment for why. "Restricted" scope guarantees every axis
      ! other than sky1/sky2/freq_axis is degenerate (size 1), so nx*ny*
      ! chan_len is always the exact contiguous byte range for this
      ! block; naxis/sky1/sky2/freq_axis/naxes are no longer needed here.
      character(len=1024) :: file_path = ' '
      integer(kind=8) :: datastart = 0_8
      integer :: chan_start = 0, chan_len = 0, nx = 0, ny = 0
      integer :: nwriters_eff = 1
      real, pointer :: data(:,:,:) => null()
   end type block_write_job_t
   type(block_write_job_t), target, save :: write_job
   integer(c_long) :: write_thread_id = 0
   logical :: write_pending = .false.
   logical :: write_failed = .false.

   ! Second, separately-shaped job type for process_one_file_general (the
   ! stages=reproject-only path) -- same field layout as reproject_cubes.
   ! f90's own block_write_job_t (general multi-non-sky-axis other_axes/
   ! other_idx, not the single freq_axis the convolve path above uses).
   ! Own write_thread_id/write_pending/write_failed are safe to share
   ! with the convolve path's above: one match_cubes run's own stages=
   ! setting is global, so a single run only ever exercises ONE of these
   ! two paths, never both concurrently.
   ! file_path/datastart (not out_unit): raw Fortran stream I/O at a
   ! computed byte offset (T19, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md)
   ! -- see process_one_file_general's own comment for why, same
   ! rationale/design as reproject_cubes.f90's own block_write_job_t.
   type :: block_write_job_general_t
      character(len=1024) :: file_path = ' '
      integer(kind=8) :: datastart = 0_8
      integer :: naxes_in(max_axes) = 0, other_axes(max_axes) = 0
      integer :: other_idx(max_axes) = 0, n_other = 0
      integer :: chan_start = 0, chan_len = 0, nx = 0, ny = 0
      integer :: nwriters_eff = 1
      real, pointer :: data(:,:,:) => null()
   end type block_write_job_general_t
   type(block_write_job_general_t), target, save :: write_job_general

   character(len=16) :: stages
   character(len=32) :: order
   logical :: do_reproject, do_convolve, convolve_first

   character(len=512) :: infiles(max_inputs), beamfiles(max_inputs)
   integer :: n_inputs
   character(len=64) :: outsuffix
   logical :: seen_outsuffix
   character(len=512) :: manifest_path
   logical :: have_manifest

   character(len=16) :: footprint_mode
   character(len=512) :: reffile
   logical :: seen_footprint_mode, seen_reffile

   ! badchan_file: one entry per infile, same comma-list convention as
   ! beamfiles -- per-band (T16, docs/dev/MULTI_BAND_TOMOGRAPHY_
   ! PLAN.md), unlike the single shared list this used to be.
   character(len=512) :: badchan_file(max_inputs)
   logical :: have_target
   real(dp) :: target_bmaj, target_bmin, target_bpa
   real(dp) :: max_common_bmaj
   logical :: have_max_common_bmaj
   real :: mem_frac_ram
   logical :: io_overlap
   ! nwriters (default 1, T21 docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md):
   ! same key/clamp formula as rm_synthesis/rmclean_cubes' own nwriters
   ! (T7/T12) and convolve_cubes' own copy of this option -- see
   ! do_block_write's own comment for how it composes with io_overlap.
   integer :: nwriters
   ! dry_run (default n, T24 docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md):
   ! same key name as rm_synthesis' own dry_run -- checks the target
   ! output disk's rotational status and writes a suggested
   ! match_cubes_dryrun.cfg (io_overlap/nwriters) instead of processing
   ! any real file, mirroring rm_synthesis' own tile_autotune.cfg
   ! convention. See run_dry_run_check's own comment for the full
   ! rationale.
   logical :: dry_run
   ! min_free_pct_after_write (default 10, UR-15 docs/dev/BRANCH_STATUS_
   ! update-reproject.md): percent of the output disk's total size that
   ! must still be free once every output is written, checked before
   ! any file is created when stages includes reproject.
   double precision :: min_free_pct_after_write
   integer :: npts
   real(dp) :: khachiyan_tol
   ! interp/interp_params (docs/dev/REPROJECT_INTERP_OPTIONS_PLAN.md):
   ! same keys/semantics as reproject_cubes.f90's own, resolved once
   ! via ast_interp_mod's resolve_interp (parse_args, below) -- the
   ! same AST scheme constant and params(:) apply to every plane of
   ! every input this run processes, through all 3 of this file's own
   ! ast_resampler call sites.
   character(len=16) :: interp
   character(len=256) :: interp_params
   integer :: interp_scheme
   double precision :: interp_params_arr(2)

   ! target_system/target_projection/target_center (docs/dev/
   ! REPROJECT_TARGET_FRAME_PLAN.md): same keys/semantics as
   ! reproject_cubes.f90's own, resolved once in parse_args (right
   ! after resolve_interp, below) via the module shared with it, so the
   ! two tools' target-frame handling cannot diverge. Only meaningful
   ! when reffile is given (stages includes reproject) -- parse_args
   ! errors clearly if target_system/target_projection are given
   ! without reffile, rather than silently ignoring them.
   character(len=16) :: target_system, target_projection, target_center
   character(len=32) :: target_lon, target_lat
   character(len=32) :: target_cdelt_lon, target_cdelt_lat
   character(len=16) :: target_hpx_nside
   integer :: status_target_grid
   ! used_target_grid/target_center_*_deg/target_cdelt*_mag/
   ! target_final_crpix (docs/dev/BRANCH_STATUS_update-reproject.md
   ! UR-7): PLAIN VALUES resolve_target_grid returns when reffile is
   ! absent -- never live AST object handles, since no AST object may
   ! cross an OpenMP thread boundary (confirmed by direct test; see
   ! target_grid_mod.f90's own header comment). Every per-thread setup
   ! block below that needs the actual Mapping calls build_target_
   ! skymap itself, from these, rather than sharing one built here.
   logical :: used_target_grid
   ! The system and projection the no-reffile grid was built in, as
   ! resolve_target_grid returns them -- target_system/
   ! target_projection as given can be blank (projection alone
   ! inherits the infiles' system; an omitted projection means sin),
   ! so every later call that rebuilds the grid takes these instead.
   character(len=16) :: grid_system, grid_projection
   double precision :: target_center_lon_deg, target_center_lat_deg
   double precision :: target_cdelt1_mag, target_cdelt2_mag
   double precision :: target_final_crpix(2)

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
   double precision :: lbnd_out(2), ubnd_out(2)

   ! --- UR-10 (docs/dev/BRANCH_STATUS_update-reproject.md): invalid-
   ! region info/warning, once per run on the final output grid ---
   integer :: ur10_skymap, ur10_skyframe, ur10_pixaxes(2), ur10_status
   ! UR-21: the output grid's cell size on the sky, this run's grid keys
   ! quoted back in the advice, and the common beam's sampling verdict.
   double precision :: out_cell1_arcsec, out_cell2_arcsec, suggested_beam_arcsec
   double precision :: in_cell1_arcsec, in_cell2_arcsec, finest_in_cell_arcsec
   logical :: grid_coarser_than_inputs
   integer :: cell_status, common_verdict
   character(len=1024) :: grid_keys_for_advice
   double precision :: ur10_invalid_frac
   integer :: ast_status
   integer :: nx_out_common, ny_out_common

   ! --- UR-15 (docs/dev/BRANCH_STATUS_update-reproject.md): unified
   ! overlap check + disk-space pre-check, both before any file write ---
   integer :: n_planes_per_infile(max_inputs)
   integer(kind=8) :: required_bytes, avail_bytes, total_bytes
   integer :: diskcheck_status

   ! --- Skip-if-already-matched (planning-doc ticket) ---
   logical :: needs_processing(max_inputs)
   integer :: manifest_unit
   logical :: manifest_exists, out_exists, matches_geom, matches_beam
   integer :: fitsstat_skip, blocksize_skip, ref_skip_unit, cand_skip_unit
   integer :: cand_ax1, cand_ax2, ref_ax1, ref_ax2

   call check_and_clamp_omp_threads()

   call parse_args(status)
   if (status.ne.0) stop 1

   call init_logging(log_level, timing_enabled, log_output_file, status)
   if (status.ne.0) then
      write(*,*) 'ERROR: cannot open log_output_file: ', trim(log_output_file)
      stop 1
   endif
   call log_message('info', 'startup', 'match_cubes run started')

   do_reproject = (trim(stages).eq.'reproject' .or. trim(stages).eq.'both')
   do_convolve = (trim(stages).eq.'convolve' .or. trim(stages).eq.'both')
   convolve_first = (trim(order).eq.'convolve_reproject')

   ! === Pre-scan phase 1: reproject footprint (order-independent of the
   ! convolve pre-scan below; only needs WCS/geometry, not beam metadata).
   ! Moved above dry_run=y's own early stop (UR-15, docs/dev/
   ! BRANCH_STATUS_update-reproject.md): reporting required-vs-available
   ! disk space in dry_run needs the grid size, which needs this pre-scan
   ! to have already run -- stages=convolve alone has no reprojection
   ! grid at all, so dry_run skips straight to the rotational-disk
   ! suggestion in that case, exactly as before. ===
   nx_out_common = 0
   ny_out_common = 0
   diskcheck_status = 0
   avail_bytes = -1_8
   if (do_reproject) then
      ast_status = 0
      call ast_begin(ast_status)

      if (.not. used_target_grid) then
         call load_wcs(reffile, wcs_ref, naxes_ref, ast_status)
         call extract_sky_mapping(wcs_ref, skymap_ref, skyframe_ref, pixaxes_ref, ast_status)
         if (ast_status.ne.0) then
            write(*,*) 'ERROR: failed to load the reference file''s WCS'
            stop 1
         endif
      endif

      lbnd_out(1) = 1.0d0
      lbnd_out(2) = 1.0d0
      ubnd_out(1) = real(naxes_ref(pixaxes_ref(1)), kind=8)
      ubnd_out(2) = real(naxes_ref(pixaxes_ref(2)), kind=8)
      if (used_target_grid) then
         write(*,'(A,F0.0,A,F0.0,A,F0.0,A,F0.0,A)') 'Target grid (no reffile)'//&
         &' own extent: [', lbnd_out(1), ',', ubnd_out(1), '] x [',&
         &lbnd_out(2), ',', ubnd_out(2), ']'
      else
         write(*,'(A,A,A,F0.0,A,F0.0,A,F0.0,A,F0.0,A)') 'Reference (', trim(reffile),&
         &') own extent: [', lbnd_out(1), ',', ubnd_out(1), '] x [',&
         &lbnd_out(2), ',', ubnd_out(2), ']'
      endif

      ! UR-15: unified overlap sanity-check across every mode, including
      ! the footprint_mode=reference check that never existed before --
      ! see check_overlap_and_size_grid's own header (target_grid_mod.f90)
      ! for the collect-every-offender-before-reporting behaviour and the
      ! daisy-chain-preserving fold. used_target_grid already ran this
      ! internally (resolve_target_grid_noref) and populated n_planes_
      ! per_infile; only run it again here when a real reffile was given.
      if (.not. used_target_grid) then
         call check_overlap_and_size_grid(skymap_ref, skyframe_ref, infiles,&
         &n_inputs, footprint_mode, .true., lbnd_out, ubnd_out,&
         &n_planes_per_infile, ast_status)
         if (ast_status.ne.0) stop 1
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

      ! === UR-15: disk-space pre-check, before any file is written. ===
      call check_disk_space(infiles(1), n_inputs, n_planes_per_infile,&
      &nx_out_common, ny_out_common, min_free_pct_after_write, required_bytes,&
      &avail_bytes, total_bytes, diskcheck_status)
   endif

   if (dry_run) then
      call run_dry_run_check(infiles(1), 'match_cubes')
      if (do_reproject .and. avail_bytes.ge.0_8) then
         write(*,'(A,F0.2,A,F0.2,A,F0.2,A,F0.2,A)') 'Disk space: output needs ',&
         &dble(required_bytes) / (1024.0d0**3), ' GB; ',&
         &dble(avail_bytes) / (1024.0d0**3), ' GB free now; ',&
         &dble(avail_bytes - required_bytes) / (1024.0d0**3),&
         &' GB would remain; minimum to keep free is ',&
         &(min_free_pct_after_write / 100.0d0) * dble(total_bytes) / (1024.0d0**3), ' GB.'
         if (diskcheck_status.ne.0) then
            write(*,'(A)') 'A real run would refuse (min_free_pct_after_write).'
         endif
      endif
      stop
   endif

   if (diskcheck_status.ne.0) stop 1

   if (do_reproject) then
      ! === UR-10: invalid-region info/warning on the final output
      ! grid, once per run -- see reproject_cubes.f90's own identical
      ! block for the full rationale. ===
      ur10_status = 0
      if (used_target_grid) then
         call build_target_skymap(grid_system, grid_projection, target_center_lon_deg,&
         &target_center_lat_deg, target_cdelt1_mag, target_cdelt2_mag,&
         &target_final_crpix, ur10_skymap, ur10_skyframe, ur10_pixaxes,&
         &ur10_status)
         call check_output_grid_validity(ur10_skymap, 0.0d0, 0.0d0,&
         &nx_out_common, ny_out_common, ur10_invalid_frac, ur10_status)
         cell_status = 0
         call output_cell_size_arcsec(ur10_skymap, ur10_skyframe,&
         &target_final_crpix(1), target_final_crpix(2),&
         &out_cell1_arcsec, out_cell2_arcsec, cell_status)
         call ast_annul(ur10_skymap, ur10_status)
         call ast_annul(ur10_skyframe, ur10_status)
      else
         call check_output_grid_validity(skymap_ref, lbnd_out(1)-1.0d0,&
         &lbnd_out(2)-1.0d0, nx_out_common, ny_out_common,&
         &ur10_invalid_frac, ur10_status)
         cell_status = 0
         call fits_reference_cell_arcsec(reffile, out_cell1_arcsec, out_cell2_arcsec,&
         &cell_status)
      endif
      if (ur10_status.eq.0) then
         call log_message('info', 'reproject', 'output grid invalid-mapping'//&
         &' fraction (no valid sky position, e.g. a projection facet seam):')
         write(*,'(A,F0.2,A)') '  ', ur10_invalid_frac*100.0d0, '% of pixels'
         if (ur10_invalid_frac.gt.0.20d0) then
            write(*,'(A,F0.2,A)') 'WARNING: ', ur10_invalid_frac*100.0d0,&
            &'% of the output grid has no valid sky mapping -- this is a'//&
            &' projection-domain artifact (e.g. a facet seam), not missing'//&
            &' input data coverage.'
         endif
      endif
   endif

   ! === UR-21: refuse an output grid that cannot sample the beam it
   ! resamples. stages=reproject, and stages=both with
   ! order=reproject_convolve, resample at each infile's own beam, so
   ! that beam is checked here; order=convolve_reproject resamples at
   ! the common beam, checked once that is known, below. ===
   if (do_reproject) then
      if (cell_status.ne.0) stop 1
      if (used_target_grid) then
         grid_keys_for_advice = 'target_system='//trim(grid_system)//&
         &' target_projection='//trim(grid_projection)
         if (len_trim(target_center).ne.0) grid_keys_for_advice =&
         &trim(grid_keys_for_advice)//' target_center='//trim(target_center)
         if (len_trim(target_lon).ne.0) grid_keys_for_advice =&
         &trim(grid_keys_for_advice)//' target_lon='//trim(target_lon)//' target_lat='//trim(target_lat)
         if (len_trim(target_cdelt_lon).ne.0) grid_keys_for_advice =&
         &trim(grid_keys_for_advice)//' target_cdelt_lon='//trim(target_cdelt_lon)//&
         &' target_cdelt_lat='//trim(target_cdelt_lat)
         if (len_trim(target_hpx_nside).ne.0) grid_keys_for_advice =&
         &trim(grid_keys_for_advice)//' target_hpx_nside='//trim(target_hpx_nside)
      else
         grid_keys_for_advice = 'footprint_mode='//trim(footprint_mode)//' reffile='//trim(reffile)
      endif
      if (.not. (do_convolve .and. convolve_first)) then
         cell_status = 0
         call check_beam_sampling(infiles, n_inputs, out_cell1_arcsec,&
         &out_cell2_arcsec, grid_keys_for_advice, cell_status)
         if (cell_status.ne.0) stop 1
      endif
   endif

   ! === UR-17: cross-file beam advisory, for stages=reproject only.
   ! This stage changes the pixel grid, not the angular resolution, so
   ! infiles with different beams produce outputs with different beams.
   ! stages=both/convolve bring every input to one common beam
   ! themselves, so the advisory would be wrong there. ===
   if (do_reproject .and. .not. do_convolve) then
      call advise_beam_consistency(infiles, n_inputs)
   endif

   ! === Pre-scan phase 2: convolve beam metadata (independent of phase 1) ===
   if (do_convolve) then
      allocate(bmaj_f(max_inputs, max_channels), bmin_f(max_inputs, max_channels))
      allocate(bpa_f(max_inputs, max_channels), isbad_f(max_inputs, max_channels))

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

         n_badchan = 0
         if (len_trim(badchan_file(i)).gt.0 .and. trim(badchan_file(i)).ne.'none') then
            call read_badchan_file(badchan_file(i), badchan_list, n_badchan, status)
            if (status.ne.0) stop 1
         endif
         call apply_badchan_list(badchan_list, n_badchan, nfreq_f(i), isbad_f(i,1:nfreq_f(i)))

         write(*,'(A,A,A,I0,A,I0,A)') 'Read ', trim(infiles(i)), ': ', nfreq_f(i),&
         &' channels, ', count(isbad_f(i,1:nfreq_f(i))), ' flagged bad'
      enddo

      ! T30: advisory only, never changes reffile= itself -- see
      ! advise_reference_band's own header comment and
      ! docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md T30.
      if (do_reproject) then
         call advise_reference_band(infiles, n_inputs, reffile, nfreq_f,&
         &sky1_f, sky2_f, cdelt1_f, cdelt2_f)
      endif

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

   ! === UR-21: with order=convolve_reproject the resampling runs at the
   ! common beam, so that beam must span at least 2.0 output cells (2.5
   ! to avoid the barely-Nyquist warning). It is refused, not enlarged:
   ! on a grid whose cells differ per axis there is no unique PA-
   ! preserving way to widen an elliptical beam, so the choice is the
   ! user's. ===
   if (do_reproject .and. do_convolve .and. convolve_first) then
      ! As for the per-infile check: the rule applies only when the
      ! output cells are larger than an infile's own, since a grid no
      ! coarser than the convolved data's cannot alias it.
      finest_in_cell_arcsec = huge(1.0d0)
      do i = 1, n_inputs
         cell_status = 0
         call fits_reference_cell_arcsec(infiles(i), in_cell1_arcsec, in_cell2_arcsec, cell_status)
         if (cell_status.ne.0) stop 1
         finest_in_cell_arcsec = min(finest_in_cell_arcsec, in_cell1_arcsec, in_cell2_arcsec)
      enddo
      grid_coarser_than_inputs = max(out_cell1_arcsec, out_cell2_arcsec) .gt.&
      &finest_in_cell_arcsec*(1.0d0 + cell_equal_rtol)
      common_verdict = common_beam_sampling(common_bmin,&
      &max(out_cell1_arcsec, out_cell2_arcsec))
      if (.not. grid_coarser_than_inputs .and. common_verdict.eq.2) then
         write(*,'(A,F0.4,A,F0.2,A)') 'WARNING: the convolution beam BMIN=', common_bmin,&
         &' arcsec spans only ', common_bmin/max(out_cell1_arcsec, out_cell2_arcsec),&
         &' output cells; the output grid is no coarser than the inputs'', so'//&
         &' reprojection adds no aliasing, but interpolated values of'//&
         &' undersampled structure are inaccurate.'
         common_verdict = 0
      endif
      if (common_verdict.eq.2) then
         suggested_beam_arcsec = max(2.5d0*max(out_cell1_arcsec, out_cell2_arcsec), common_bmaj)
         suggested_beam_arcsec = ceiling(suggested_beam_arcsec*100.0d0)/100.0d0
         if (have_target) then
            write(*,'(A,F0.4,A,F0.2,A)') 'ERROR: target_bmin=', common_bmin,&
            &' arcsec spans only ', common_bmin/max(out_cell1_arcsec, out_cell2_arcsec),&
            &' output cells; at least 2.0 are needed.'
         else
            write(*,'(A,F0.4,A,F0.2,A)') 'ERROR: the derived common beam (BMIN=',&
            &common_bmin, ' arcsec) spans only ',&
            &common_bmin/max(out_cell1_arcsec, out_cell2_arcsec),&
            &' output cells; at least 2.0 are needed.'
         endif
         write(*,'(A,F0.4,A,F0.4,A,F0.4,A)') '  Output cells: ', out_cell1_arcsec, ' x ',&
         &out_cell2_arcsec, ' arcsec, so BMIN >= ',&
         &2.0d0*max(out_cell1_arcsec, out_cell2_arcsec), ' arcsec is required.'
         write(*,'(A,F0.2,A,F0.2,A)') '  Give an explicit target beam, e.g. target_bmaj=',&
         &suggested_beam_arcsec, ' target_bmin=', suggested_beam_arcsec, ' target_bpa=0'
         write(*,'(A)') '  (2.5 output cells, or the common beam''s BMAJ if larger,'//&
         &' so every channel can still be convolved to it), or use coarser-'//&
         &'beam data or a finer output grid.'
         write(*,'(A)') 'Quitting now...'
         stop 1
      else if (common_verdict.eq.1) then
         write(*,'(A,F0.4,A,F0.2,A)') 'WARNING: the convolution beam BMIN=', common_bmin,&
         &' arcsec spans ', common_bmin/max(out_cell1_arcsec, out_cell2_arcsec),&
         &' output cells -- barely Nyquist (2.5 or more is recommended).'
      endif
   endif

   ! === Skip-if-already-matched pre-flight (planning-doc ticket) ===
   ! Runs AFTER both pre-scans (target grid/beam already known) but
   ! BEFORE any file is processed -- every safety check and skip
   ! decision for the WHOLE batch happens up front, so a bad run (a
   ! stale output already on disk, a stale manifest already on disk)
   ! fails fast rather than partway through a multi-hour job. Standing
   ! rule: a file this program is about to write that already exists on
   ! disk is always refused, never silently reused or overwritten --
   ! this applies regardless of what THIS run's own skip decision would
   ! have been, since ambiguous pre-existing state is refused outright,
   ! never interpreted.
   if (have_manifest) then
      inquire(file=trim(manifest_path), exist=manifest_exists)
      if (manifest_exists) then
         write(*,*) 'ERROR: manifest already exists, refusing to overwrite: ',&
         &trim(manifest_path)
         write(*,*) 'Remove it first if you intend to regenerate it.'
         stop 1
      endif
   endif

   do i = 1, n_inputs
      inquire(file=trim(strip_fits_ext(infiles(i)))//trim(outsuffix), exist=out_exists)
      if (out_exists) then
         write(*,*) 'ERROR: output path already exists, refusing to proceed'//&
         &' (stale output from a previous run? remove it first): ',&
         &trim(strip_fits_ext(infiles(i)))//trim(outsuffix)
         stop 1
      endif
   enddo

   ! needs_processing(i): does file i need EITHER stage it was asked to
   ! run? Compared against the ALREADY-COMPUTED shared target grid/beam
   ! from the pre-scans above -- independently per file, no cross-file
   ! coupling in the decision itself (see planning doc's own "worked
   ! through" section). Axis numbering matches whichever code path will
   ! process this file: sky1_f(i)/sky2_f(i) (from read_axis_info,
   ! assumed identical on the reference too -- the SAME assumption
   ! process_one_file_restricted's own copy_axis_keywords calls already
   ! make) when do_convolve, since that's what process_one_file_restricted
   ! uses; AST-derived pixaxes_in/pixaxes_ref (independently extracted per
   ! file) only for the reproject-ALONE path, matching
   ! process_one_file_general's own convention.
   ! Skip-if-already-matched needs a real reffile file to compare
   ! against (sky_wcs_matches_target reads its CTYPE/CRVAL/etc via
   ! CFITSIO) -- not available when used_target_grid (no reffile, a
   ! freshly-built in-memory grid). Skipped in that case, same as
   ! reproject_cubes.f90's own identical case: a pure performance
   ! optimisation, not a correctness requirement, and the freshly-built
   ! grid is derived from these exact infiles' own footprints, so an
   ! exact pre-existing match would be a rare coincidence anyway.
   if (do_reproject .and. .not. used_target_grid) then
      fitsstat_skip = 0
      call safe_ftopen(ref_skip_unit, trim(reffile), 0, blocksize_skip, fitsstat_skip)
      if (fitsstat_skip.ne.0) then
         write(*,*) 'ERROR: cannot reopen reference file for the geometry check: ',&
         &trim(reffile)
         stop 1
      endif
   endif

   do i = 1, n_inputs
      matches_geom = .true.
      matches_beam = .true.

      if (do_reproject .and. used_target_grid) then
         matches_geom = .false.
      else if (do_reproject) then
         fitsstat_skip = 0
         call safe_ftopen(cand_skip_unit, trim(infiles(i)), 0, blocksize_skip, fitsstat_skip)
         if (fitsstat_skip.ne.0) then
            write(*,*) 'ERROR: cannot reopen input for the geometry check: ',&
            &trim(infiles(i))
            stop 1
         endif
         if (do_convolve) then
            cand_ax1 = sky1_f(i)
            cand_ax2 = sky2_f(i)
            ref_ax1 = sky1_f(i)
            ref_ax2 = sky2_f(i)
            matches_geom = sky_wcs_matches_target(cand_skip_unit, cand_ax1, cand_ax2,&
            &naxes_f(i,sky1_f(i)), naxes_f(i,sky2_f(i)), ref_skip_unit, ref_ax1, ref_ax2,&
            &lbnd_out(1)-1.0d0, lbnd_out(2)-1.0d0, nx_out_common, ny_out_common)
         else
            call load_wcs(infiles(i), wcs_in, naxes_in, ast_status)
            call extract_sky_mapping(wcs_in, skymap_in, skyframe_in, pixaxes_in, ast_status)
            if (ast_status.ne.0) then
               write(*,*) 'ERROR: failed to read input''s WCS for the geometry check: ',&
               &trim(infiles(i))
               stop 1
            endif
            matches_geom = sky_wcs_matches_target(cand_skip_unit, pixaxes_in(1), pixaxes_in(2),&
            &naxes_in(pixaxes_in(1)), naxes_in(pixaxes_in(2)), ref_skip_unit,&
            &pixaxes_ref(1), pixaxes_ref(2), lbnd_out(1)-1.0d0, lbnd_out(2)-1.0d0,&
            &nx_out_common, ny_out_common)
            call ast_annul(skymap_in, ast_status)
            call ast_annul(skyframe_in, ast_status)
            call ast_annul(wcs_in, ast_status)
         endif
         call safe_ftclos(cand_skip_unit, fitsstat_skip)
      endif

      if (do_convolve) then
         matches_beam = beam_matches_target(nfreq_f(i), bmaj_f(i,1:nfreq_f(i)),&
         &bmin_f(i,1:nfreq_f(i)), bpa_f(i,1:nfreq_f(i)), isbad_f(i,1:nfreq_f(i)),&
         &common_bmaj, common_bmin, common_bpa)
      endif

      needs_processing(i) = .not. (matches_geom .and. matches_beam)
      if (.not. needs_processing(i)) then
         write(*,'(A,A,A)') 'SKIP: ', trim(infiles(i)),&
         &' already matches target geometry/beam -- no output written, use it directly'
      endif
   enddo

   if (do_reproject .and. .not. used_target_grid) then
      fitsstat_skip = 0
      call safe_ftclos(ref_skip_unit, fitsstat_skip)
   endif

   ! === Per-file processing ===
   do i = 1, n_inputs
      if (.not. needs_processing(i)) cycle
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
         ! Plain filename (NOT '!'-prefixed) -- see process_one_file_
         ! restricted's own FTINIT comment below for the full story;
         ! this call site had the identical clobber bug, just baked into
         ! the filename string here rather than inside the FTINIT call
         ! itself (process_one_file_general's own FTINIT was already
         ! correct -- the bug was in what the caller passed it).
         call process_one_file_general(reffile, infiles(i),&
         &trim(strip_fits_ext(infiles(i)))//trim(outsuffix), pixaxes_ref,&
         &naxes_in, pixaxes_in, lbnd_out, ubnd_out, mem_frac_ram,&
         &io_overlap, nwriters, interp_scheme, interp_params_arr,&
         &used_target_grid, grid_system, grid_projection, target_center_lon_deg,&
         &target_center_lat_deg, target_cdelt1_mag, target_cdelt2_mag,&
         &target_final_crpix, status)
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
         &trim(strip_fits_ext(infiles(i)))//trim(outsuffix), do_reproject, convolve_first,&
         &naxis_f(i), sky1_f(i), sky2_f(i), freq_axis_f(i), naxes_f(i,:),&
         &cdelt1_f(i), cdelt2_f(i), nfreq_f(i), bmaj_f(i,1:nfreq_f(i)),&
         &bmin_f(i,1:nfreq_f(i)), bpa_f(i,1:nfreq_f(i)), isbad_f(i,1:nfreq_f(i)),&
         &common_bmaj, common_bmin, common_bpa, reffile, pixaxes_ref,&
         &nx_out_common, ny_out_common, lbnd_out, ubnd_out, mem_frac_ram,&
         &io_overlap, nwriters, interp_scheme, interp_params_arr,&
         &used_target_grid, grid_system, grid_projection, target_center_lon_deg,&
         &target_center_lat_deg, target_cdelt1_mag, target_cdelt2_mag,&
         &target_final_crpix, status)
         if (status.ne.0) then
            write(*,*) 'ERROR: failed to write output for: ', trim(infiles(i))
            stop 1
         endif
      endif
      write(*,*) 'OK: wrote ', trim(strip_fits_ext(infiles(i)))//trim(outsuffix)
   enddo

   ! Manifest: one line per input, machine-readable record of whether it
   ! was skipped or processed and the effective path to use downstream --
   ! never inferred from filesystem state by anything consuming this
   ! file (planning-doc ticket's own explicit rationale: a stray/stale
   ! output left by an unrelated earlier run must never be silently
   ! mistaken for this run's own result).
   if (have_manifest) then
      open(newunit=manifest_unit, file=trim(manifest_path), status='new',&
      &action='write', iostat=fitsstat_skip)
      if (fitsstat_skip.ne.0) then
         write(*,*) 'ERROR: failed to create manifest: ', trim(manifest_path)
         stop 1
      endif
      do i = 1, n_inputs
         if (needs_processing(i)) then
            write(manifest_unit,'(A,A,A,A,A)') trim(infiles(i)), char(9),&
            &'PROCESSED', char(9), trim(strip_fits_ext(infiles(i)))//trim(outsuffix)
         else
            write(manifest_unit,'(A,A,A,A,A)') trim(infiles(i)), char(9),&
            &'SKIPPED', char(9), trim(infiles(i))
         endif
      enddo
      close(manifest_unit)
      write(*,*) 'OK: wrote manifest: ', trim(manifest_path)
   endif

   if (do_reproject) then
      if (.not. used_target_grid) then
         call ast_annul(skymap_ref, ast_status)
         call ast_annul(skyframe_ref, ast_status)
         call ast_annul(wcs_ref, ast_status)
      endif
      call ast_end(ast_status)
   endif

   if (allocated(bmaj_f)) deallocate(bmaj_f, bmin_f, bpa_f, isbad_f)
   write(*,*) 'OK: all inputs processed.'
   call timer_report_summary()
   call log_message('info', 'finalize', 'match_cubes run completed')

contains

   !===========================================================
   ! CLI / config parsing (adapted from convolve_cubes.f90, extended
   ! with reproject-stage keys)
   !===========================================================

   subroutine parse_args(status)
      integer, intent(out) :: status
      ! this_arg/cli_val and the raw_* CSV-list staging variables below
      ! hold a WHOLE comma-separated infiles=/beamfiles=/badchan_file=
      ! argument, not one path -- up to max_inputs (50) entries, each up
      ! to 512 chars (infiles(:)'s own per-entry length). 16384 gives
      ! generous headroom over the worst case (~25,600 chars); 512 here
      ! silently truncated a real run's own infiles= argument (found via
      ! scripts/run_pipeline.sh's symlink-redirection scheme, which
      ! lengthens every path enough to cross 512 for just 4 real files).
      character(len=16384) :: this_arg, cli_val
      character(len=512) :: cli_key, cfgfile
      character(len=16384) :: raw_infiles, raw_beamfiles, raw_badchan_file
      integer :: argc, iarg
      logical :: has_kv, have_cfgfile, seen_infiles, seen_stages

      status = 0
      n_inputs = 0
      outsuffix = ' '
      seen_outsuffix = .false.
      manifest_path = ' '
      have_manifest = .false.
      stages = ' '
      seen_stages = .false.
      order = 'convolve_reproject'
      footprint_mode = ' '
      seen_footprint_mode = .false.
      reffile = ' '
      seen_reffile = .false.
      badchan_file = ' '
      have_target = .false.
      target_bmaj = 0.0d0
      target_bmin = 0.0d0
      target_bpa = 0.0d0
      have_max_common_bmaj = .false.
      max_common_bmaj = 0.0d0
      mem_frac_ram = 0.25
      io_overlap = .false.
      nwriters = 1
      dry_run = .false.
      min_free_pct_after_write = 10.0d0
      log_level = 'info'
      timing_enabled = .false.
      log_output_file = ' '
      npts = 2000
      khachiyan_tol = 1.0d-5
      ! Default changed from linear to lanczos3, 2026-10-13 -- see
      ! reproject_cubes.f90's own identical comment (UR-12 benchmark).
      interp = 'lanczos3'
      interp_params = ' '
      target_system = ' '
      target_projection = ' '
      target_center = ' '
      target_lon = ' '
      target_lat = ' '
      target_cdelt_lon = ' '
      target_cdelt_lat = ' '
      target_hpx_nside = ' '
      raw_infiles = ' '
      raw_beamfiles = ' '
      raw_badchan_file = ' '
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
            &raw_beamfiles, raw_badchan_file, seen_infiles, seen_stages, status)
            if (status.ne.0) return
            iarg = iarg + 1
         endif
      enddo

      if (have_cfgfile) call read_cfg_file(cfgfile, raw_infiles, raw_beamfiles,&
      &raw_badchan_file, seen_infiles, seen_stages, status)
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
      if (len_trim(raw_badchan_file).gt.0) then
         if (cfg_csv_count(raw_badchan_file).ne.n_inputs) then
            write(*,*) 'ERROR: badchan_file must list exactly ', n_inputs,&
            &' entries (one per infile; give an entry the literal value',&
            &' ''none'' -- e.g. ''none,file2.txt'' -- for a file with no',&
            &' manual bad-channel list)'
            status = -1
            return
         endif
         do i = 1, n_inputs
            call cfg_csv_get_item(raw_badchan_file, i, badchan_file(i))
         enddo
         ! T34 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): once badchan_file=
         ! is given at all, every position must be explicit -- a real path
         ! or the literal value 'none' -- matching rm_synthesis's own T33
         ! convention. Omitting badchan_file= entirely is untouched (this
         ! whole block is skipped, raw_badchan_file stays blank); this only
         ! rejects a blank POSITION inside an otherwise-given list, which
         ! was previously silently accepted as "nothing to list."
         do i = 1, n_inputs
            if (len_trim(badchan_file(i)).eq.0) then
               write(*,*) 'ERROR: badchan_file entry ', i, ' is blank --',&
               &' use the literal value ''none'' for a file with nothing',&
               &' to list'
               status = -1
               return
            endif
         enddo
      endif

      ! footprint_mode and reffile (UR-7: reffile= is now optional) must
      ! be given together -- reffile= fully determines the output grid,
      ! full stop -- or both omitted, in which case target_system=/
      ! target_projection= build one instead (resolve_target_grid,
      ! below, enforces the actual mutual-exclusivity with every
      ! target_* key). Same pairing reproject_cubes.f90 enforces between
      ! its own mode=/reffile=.
      if ((trim(stages).eq.'reproject' .or. trim(stages).eq.'both')) then
         if (seen_footprint_mode .neqv. seen_reffile) then
            write(*,*) 'ERROR: footprint_mode and reffile must be given together,'//&
            &' or both omitted (then target_system=/target_projection= build the'//&
            &' output grid instead)'
            status = -1
            return
         endif
         if (seen_footprint_mode .and. trim(footprint_mode).ne.'intersection' .and.&
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
      if (min_free_pct_after_write.lt.0.0d0 .or. min_free_pct_after_write.gt.100.0d0) then
         write(*,*) 'ERROR: min_free_pct_after_write must be between 0 and 100 (percent), got ',&
         &min_free_pct_after_write
         status = -1
         return
      endif
      if (npts.lt.12) then
         write(*,*) 'ERROR: npts must be at least 12, got ', npts
         status = -1
         return
      endif
      call resolve_interp(interp, interp_params, interp_scheme,&
      &interp_params_arr, status)
      if (status.ne.0) return

      ! target_system/target_projection/target_center/target_lon/
      ! target_lat/target_cdelt_lon/target_cdelt_lat (docs/dev/
      ! BRANCH_STATUS_update-reproject.md UR-7): resolved once here,
      ! before reffile is first used (process_one_file_general/
      ! process_one_file_restricted below), via the module shared with
      ! reproject_cubes.f90 so the two tools' target-frame handling
      ! cannot diverge. Only meaningful when stages includes reproject
      ! -- there is nothing for these to apply to otherwise (stages=
      ! convolve alone). When reffile is absent, used_target_grid comes
      ! back true and naxes_ref/pixaxes_ref (host-associated with the
      ! main program, like every other reproject-stage variable this
      ! subroutine already sets) are populated exactly as load_wcs(
      ! reffile,...)+extract_sky_mapping would have, plus the PLAIN
      ! VALUES (target_center_lon_deg/target_center_lat_deg/
      ! target_cdelt1_mag/target_cdelt2_mag/target_final_crpix) needed
      ! to rebuild the actual Mapping/SkyFrame later -- see this
      ! program's own declaration comment for why plain values, not AST
      ! object handles. footprint_mode is forced to 'reference' in that
      ! case: there is no reffile for union/intersection to be relative
      ! to, since the built grid is already exactly the infiles' own
      ! combined footprint.
      used_target_grid = .false.
      if ((trim(stages).eq.'reproject' .or. trim(stages).eq.'both')) then
         status_target_grid = 0
         call resolve_target_grid(seen_reffile, target_system, target_projection,&
         &target_center, target_lon, target_lat, target_cdelt_lon, target_cdelt_lat,&
         &target_hpx_nside, infiles, n_inputs, used_target_grid, grid_system,&
         &grid_projection, target_center_lon_deg,&
         &target_center_lat_deg, target_cdelt1_mag, target_cdelt2_mag,&
         &target_final_crpix, pixaxes_ref, naxes_ref, n_planes_per_infile,&
         &status_target_grid)
         if (status_target_grid.ne.0) then
            status = -1
            return
         endif
         if (used_target_grid) footprint_mode = 'reference'
      else if (len_trim(target_system).ne.0 .or. len_trim(target_projection).ne.0&
      &.or. len_trim(target_center).ne.0 .or. len_trim(target_lon).ne.0 .or.&
      &len_trim(target_lat).ne.0 .or. len_trim(target_cdelt_lon).ne.0 .or.&
      &len_trim(target_cdelt_lat).ne.0 .or. len_trim(target_hpx_nside).ne.0) then
         write(*,*) 'ERROR: target_system/target_projection/target_center/'//&
         &'target_lon/target_lat/target_cdelt_lon/target_cdelt_lat/'//&
         &'target_hpx_nside require stages to include reproject'
         status = -1
         return
      endif
   end subroutine parse_args

   subroutine apply_kv(key, val, raw_infiles, raw_beamfiles, raw_badchan_file,&
   &seen_infiles, seen_stages, status)
      character(len=*), intent(in) :: key, val
      character(len=*), intent(inout) :: raw_infiles, raw_beamfiles, raw_badchan_file
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
      case ('manifest')
         manifest_path = val
         have_manifest = .true.
      case ('footprint_mode')
         footprint_mode = val
         seen_footprint_mode = .true.
      case ('reffile')
         reffile = val
         seen_reffile = .true.
      case ('badchan_file')
         raw_badchan_file = val
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
      case ('io_overlap')
         io_overlap = flag_from_value_match(val)
      case ('nwriters')
         read(val, *, iostat=ios) nwriters
         if (ios.ne.0 .or. nwriters.lt.1) then
            write(*,*) 'ERROR: nwriters must be an integer >= 1'
            status = -1
            return
         endif
      case ('dry_run')
         dry_run = flag_from_value_match(val)
      case ('min_free_pct_after_write')
         read(val, *, iostat=ios) min_free_pct_after_write
         if (ios.ne.0) then
            write(*,*) 'ERROR: min_free_pct_after_write must be a number'
            status = -1
            return
         endif
      case ('log_level')
         log_level = trim(val)
      case ('timing_enabled')
         timing_enabled = flag_from_value_logging(val)
      case ('log_output_file')
         log_output_file = trim(val)
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
      case ('interp')
         interp = trim(val)
      case ('interp_params')
         interp_params = trim(val)
      case ('target_system')
         target_system = trim(val)
      case ('target_projection')
         target_projection = trim(val)
      case ('target_center')
         target_center = trim(val)
      case ('target_lon')
         target_lon = trim(val)
      case ('target_lat')
         target_lat = trim(val)
      case ('target_cdelt_lon')
         target_cdelt_lon = trim(val)
      case ('target_cdelt_lat')
         target_cdelt_lat = trim(val)
      case ('target_hpx_nside')
         target_hpx_nside = trim(val)
      case default
         write(*,*) 'ERROR: unrecognised key "', key, '"'
         status = -1
         return
      end select
   end subroutine apply_kv

   subroutine read_cfg_file(cfgfile, raw_infiles, raw_beamfiles, raw_badchan_file,&
   &seen_infiles, seen_stages, status)
      character(len=*), intent(in) :: cfgfile
      character(len=*), intent(inout) :: raw_infiles, raw_beamfiles, raw_badchan_file
      logical, intent(inout) :: seen_infiles, seen_stages
      integer, intent(out) :: status
      character(len=16384) :: line, val
      character(len=512) :: key
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
         &raw_badchan_file, seen_infiles, seen_stages, status)
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
      write(*,'(A)') '    [beamfiles=<spec1>[,<spec2>...]] [badchan_file=<file1>[,<file2>...]]'
      write(*,'(A)') '    [target_bmaj=<arcsec> target_bmin=<arcsec> target_bpa=<deg>]'
      write(*,'(A)') '    [max_common_bmaj=<arcsec>] [mem_frac_ram=<fraction>]'
      write(*,'(A)') '    [interp=<scheme>] [interp_params=<p1>[,<p2>]]'
      write(*,'(A)') '  -- when stages includes reproject, footprint_mode=/reffile= must'//&
      &' be given together, or both omitted and target_system=/'
      write(*,'(A)') '     target_projection= build the output grid instead:'
      write(*,'(A)') '    [target_system=<system>] [target_projection=<proj>]'//&
      &' [target_center=<center>] [target_lon=<lon> target_lat=<lat>]'
      write(*,'(A)') '    [target_cdelt_lon=<deg> target_cdelt_lat=<deg>]'//&
      &' [target_hpx_nside=<n>]'
      write(*,'(A)') '    [outsuffix=<suffix>] [npts=<n>] [khachiyan_tol=<tol>]'
      write(*,'(A)') '    [manifest=<path>] [io_overlap=y|n] [nwriters=<n>] [dry_run=y|n]'
      write(*,'(A)') '    [min_free_pct_after_write=<percent>]'
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
      write(*,'(A)') 'footprint_mode/reffile: when stages includes reproject, must be'//&
      &' given together -- reffile already fully determines the output grid,'
      write(*,'(A)') '  full stop; combining it with ANY target_* key below is a hard'//&
      &' error, not a silent ignore -- or both omitted, in which case'
      write(*,'(A)') '  target_system=/target_projection= below build a fresh output'//&
      &' grid instead, fit to the infiles'' own combined footprint. Same'
      write(*,'(A)') '  semantics as the standalone reproject_cubes tool -- see its own'//&
      &' --help.'
      write(*,'(A)') ''
      write(*,'(A)') 'interp/interp_params (default lanczos3 since 2026-10-13, used'//&
      &' whenever stages includes reproject): same semantics as the standalone'
      write(*,'(A)') '  reproject_cubes tool -- see its own --help, including the'//&
      &' measured kernel-width defaults for gauss/sincgauss/blockave (UR-19),'
      write(*,'(A)') '  which are reported when used, and interp=exact, the'//&
      &' flux-conserving overlap method, in every stages/order combination.'
      write(*,'(A)') ''
      write(*,'(A)') 'target_system/target_projection/target_center/target_lon/target_lat/'//&
      &'target_cdelt_lon/target_cdelt_lat/target_hpx_nside: only meaningful when'
      write(*,'(A)') '  stages includes reproject AND reffile is absent (see'//&
      &' footprint_mode/reffile above) -- same semantics as the standalone'
      write(*,'(A)') '  reproject_cubes tool, see its own --help (incl.'//&
      &' target_center=user_defined''s target_lon=/target_lat= notation,'
      write(*,'(A)') '  target_cdelt_lon=/target_cdelt_lat= defaulting to the smallest'//&
      &' |CDELT| across the infiles, target_projection=tan|sin|car|ait|hpx (sin'
      write(*,'(A)') '  when omitted), and target_hpx_nside=<power of 2>, only valid'//&
      &' with target_projection=hpx, a HEALPix-aligned grid whose centre, pixel'
      write(*,'(A)') '  scale and reference pixel all follow from Nside alone --'//&
      &' target_center=/target_lon=/target_lat=/target_cdelt_lon=/'
      write(*,'(A)') '  target_cdelt_lat= are then errors).'
      write(*,'(A)') ''
      write(*,'(A)') 'Beam sampling: onto output cells larger than an infile''s own,'//&
      &' its beam minor axis (BMIN, header or narrowest BEAMS-table channel) must'
      write(*,'(A)') '  span at least 2.0 output cells, or the run is refused before'//&
      &' anything is written; 2.0 to 2.5 cells gives a barely-Nyquist WARNING.'
      write(*,'(A)') '  An infile with no beam information is refused onto coarser'//&
      &' cells. No beam is ever widened automatically: convolve first'
      write(*,'(A)') '  (stages=both with order=convolve_reproject, which applies the same'//&
      &' rule to the convolution beam); the refusal prints what to use.'
      write(*,'(A)') ''
      write(*,'(A)') 'RM synthesis advisory: for multi-band data feeding RM synthesis,'//&
      &' prefer reffile= -- footprint_mode=union/intersection only work as'
      write(*,'(A)') '  intended relative to a reference file. The no-reffile'//&
      &' target_system= path DOES also produce one shared grid for every'
      write(*,'(A)') '  band, but only within a SINGLE match_cubes invocation -- pass'//&
      &' every band together in one infiles= list, never run match_cubes'
      write(*,'(A)') '  separately per band with no reffile (each separate run would'//&
      &' independently resolve its own centroid/footprint-union and could'
      write(*,'(A)') '  land on a different grid).'
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
      write(*,'(A)') ''
      write(*,'(A)') 'Skip-if-already-matched: before processing, each input is checked'//&
      &' against the already-computed target grid/beam (whichever is relevant'
      write(*,'(A)') '  for the requested stages). A file that already matches -- e.g.'//&
      &' it IS the reference, or already went through an identical prior run --'
      write(*,'(A)') '  is left untouched: no output is written for it, and its own'//&
      &' original path is what downstream tools should use directly. Tight'
      write(*,'(A)') '  tolerances only ("already processed identically", not "close'//&
      &' enough") -- a false match would silently misalign downstream RM'
      write(*,'(A)') '  synthesis, so anything short of a near-exact match is processed'//&
      &' as normal. A pre-existing output file at the path this run would'
      write(*,'(A)') '  write to (regardless of whether this run''s own decision was'//&
      &' skip or process) always aborts the whole run before anything is'
      write(*,'(A)') '  touched -- never silently reused or overwritten.'
      write(*,'(A)') ''
      write(*,'(A)') 'manifest=<path> (optional): write a machine-readable record --'//&
      &' one line per input, tab-separated "<infile> SKIPPED|PROCESSED'
      write(*,'(A)') '  <effective_path>" -- so a caller (e.g. scripts/run_pipeline.sh)'//&
      &' can chain the right path per file without guessing from the'
      write(*,'(A)') '  filesystem. Aborts if the manifest path already exists, same'//&
      &' "never silently overwrite" rule as everywhere else.'
      write(*,'(A)') ''
      write(*,'(A)') 'io_overlap (default n): y -- write each block on a background thread,'//&
      &' overlapped with the NEXT block''s own read+process, instead of blocking'
      write(*,'(A)') '  on the write before starting it. Only one background write is'//&
      &' ever in flight at a time.'
      write(*,'(A)') ''
      write(*,'(A)') 'nwriters (default 1): split one block''s own planes into this many'//&
      &' disjoint concurrent writers -- same key/clamp as rm_synthesis'' and'
      write(*,'(A)') '  rmclean_cubes'' own nwriters (max(1, min(nwriters, OMP thread count))).'//&
      &' Orthogonal to io_overlap: io_overlap decides WHEN a block''s write runs;'
      write(*,'(A)') '  nwriters decides how many concurrent writers do it once dispatched.'
      write(*,'(A)') ''
      write(*,'(A)') 'dry_run (default n): y -- check the target output disk''s own'//&
      &' rotational status and write a suggested match_cubes_dryrun.cfg'
      write(*,'(A)') '  (io_overlap/nwriters), instead of processing any real file.'//&
      &' Advisory only -- does not change any clamp, just suggests a'
      write(*,'(A)') '  starting value. Touches no data. When stages includes reproject,'//&
      &' also reports the disk-space check below.'
      write(*,'(A)') ''
      write(*,'(A)') 'min_free_pct_after_write (default 10): percent [0,100] of the output'//&
      &' disk''s TOTAL size that must still be free once every output file'
      write(*,'(A)') '  is written. Checked before any file is created when stages includes'//&
      &' reproject; the run refuses, printing required/free/remaining/minimum,'
      write(*,'(A)') '  if it would leave less.'
      write(*,'(A)') ''
      write(*,'(A)') 'Optional keys (CLI or config, logging/timing):'
      write(*,'(A)') '  log_level       = error|warn|info|debug (default info)'
      write(*,'(A)') '  timing_enabled  = y|n -- print a stage timing summary (default n)'
      write(*,'(A)') '  log_output_file = path -- append log/timing output to this file'//&
      &' instead of stdout (default empty = stdout)'
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

   logical function beam_matches_target(nfreq, bmaj, bmin, bpa, isbad,&
   &common_bmaj, common_bmin, common_bpa) result(matches)
      !! Verbatim duplicate of convolve_cubes.f90's own beam_matches_target
      !! (this file's own established convention: adapt/duplicate small
      !! shared logic rather than factor out a module -- see this file's
      !! own top-of-file comment). True only if EVERY good channel of
      !! this one input already has (bmaj,bmin,bpa) equal to the common
      !! target beam within a tight absolute tolerance -- "already-
      !! processed-identically", not "close enough to convolve
      !! negligibly".
      integer, intent(in) :: nfreq
      real(dp), intent(in) :: bmaj(nfreq), bmin(nfreq), bpa(nfreq)
      logical, intent(in) :: isbad(nfreq)
      real(dp), intent(in) :: common_bmaj, common_bmin, common_bpa
      real(dp), parameter :: tol_beam = 1.0d-6
      integer :: k

      matches = .true.
      do k = 1, nfreq
         if (isbad(k)) cycle
         if (abs(bmaj(k)-common_bmaj).gt.tol_beam .or.&
         &abs(bmin(k)-common_bmin).gt.tol_beam .or.&
         &abs(bpa(k)-common_bpa).gt.tol_beam) then
            matches = .false.
            return
         endif
      enddo
   end function beam_matches_target

   function flag_from_value_match(val) result(flag)
      !! Same convention as rm_synthesis_mod.f90's own flag_from_value:
      !! first non-blank character '1'/'y'/'Y'/'t'/'T' -> true, anything
      !! else (including blank) -> false.
      character(len=*), intent(in) :: val
      logical :: flag
      character(len=64) :: t
      integer :: i

      t = adjustl(val)
      do i = 1, len_trim(t)
         t(i:i) = achar(iachar(t(i:i)) + merge(32, 0, t(i:i).ge.'A'.and.t(i:i).le.'Z'))
      enddo
      flag = .false.
      if (len_trim(t).eq.0) return
      if (t(1:1).eq.'1' .or. t(1:1).eq.'y' .or. t(1:1).eq.'t') flag = .true.
   end function flag_from_value_match

   function strip_fits_ext(filename) result(base)
      !! Output-name helper: verbatim port of convolve_cubes.f90's own --
      !! see there for the full rationale (avoids the double-extension
      !! "name.fits_CONV.FITS" this used to produce). Strips whatever
      !! the trailing extension is, not just ".fits".
      character(len=*), intent(in) :: filename
      character(len=len(filename)) :: base
      integer :: n, i, slash, dot

      base = filename
      n = len_trim(filename)
      if (n.lt.1) return

      slash = 0
      do i = n, 1, -1
         if (filename(i:i).eq.'/') then
            slash = i
            exit
         endif
      enddo

      dot = 0
      do i = n, slash+1, -1
         if (filename(i:i).eq.'.') then
            dot = i
            exit
         endif
      enddo

      if (dot.gt.slash+1) base = filename(1:dot-1)
   end function strip_fits_ext

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
      call safe_ftopen(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(filename)
         status = -1
         call free_fits_unit(unit)
         return
      endif

      call FTGKYJ(unit, 'NAXIS', naxis, comment, fitsstat)
      if (fitsstat.ne.0 .or. naxis.lt.2 .or. naxis.gt.max_axes) then
         write(*,*) 'ERROR: bad or missing NAXIS in: ', trim(filename)
         status = -1
         call safe_ftclos(unit, fitsstat)
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
            call safe_ftclos(unit, fitsstat)
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
               call safe_ftclos(unit, fitsstat)
               return
            endif
         else if (ctype(1:3).eq.'DEC') then
            if (sky2.eq.0) then
               sky2 = k
            else
               write(*,*) 'ERROR: more than one DEC-like axis in: ', trim(filename)
               status = -1
               call safe_ftclos(unit, fitsstat)
               return
            endif
         else if (ctype(1:4).eq.'FREQ') then
            freq_axis = k
         endif
      enddo

      if (sky1.eq.0 .or. sky2.eq.0) then
         write(*,*) 'ERROR: could not identify RA/DEC sky axes in: ', trim(filename)
         status = -1
         call safe_ftclos(unit, fitsstat)
         return
      endif
      if (freq_axis.eq.0) then
         write(*,*) 'ERROR: could not identify a FREQ axis in: ', trim(filename)
         status = -1
         call safe_ftclos(unit, fitsstat)
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
               call safe_ftclos(unit, fitsstat)
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
         call safe_ftclos(unit, fitsstat)
         return
      endif
      write(axstr,'(I0)') sky2
      fitsstat = 0
      call FTGKYD(unit, 'CDELT'//trim(axstr), cdelt2, comment, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: missing CDELT for the DEC axis in: ', trim(filename)
         status = -1
         call safe_ftclos(unit, fitsstat)
         return
      endif

      call check_no_rotation(unit, sky1, sky2, filename, status)
      if (status.ne.0) then
         call safe_ftclos(unit, fitsstat)
         return
      endif

      call safe_ftclos(unit, fitsstat)
   end subroutine read_axis_info

   subroutine advise_reference_band(infiles, n_inputs, reffile, nfreq_f,&
   &sky1_f, sky2_f, cdelt1_f, cdelt2_f)
      !! T30: advisory only -- never changes reffile itself. Groups
      !! infiles by identical sky grid (CRVAL1/CRVAL2/CDELT1/CDELT2 --
      !! same pointing AND same pixel scale, confirmed directly on the
      !! real WALLABY+EMU data to be exactly what makes reproject cheap:
      !! WALLABY-Q and WALLABY-U share identical CRVAL/CRPIX/NAXIS, EMU-Q
      !! and EMU-U likewise share identical CRVAL/CRPIX/NAXIS with each
      !! other but NOT with WALLABY -- reprojecting within a group is a
      !! near-identity transform (~0.44s/plane measured), across groups
      !! is genuine resampling (~22.75s/plane measured, ~52x). Sums each
      !! group's total plane count and recommends reffile come from
      !! whichever group has the most, since that group's planes get the
      !! cheap treatment for free while every other group's planes pay
      !! the expensive one regardless of which file within it is chosen.
      character(len=*), intent(in) :: infiles(:)
      integer, intent(in) :: n_inputs
      character(len=*), intent(in) :: reffile
      integer, intent(in) :: nfreq_f(:)
      integer, intent(in) :: sky1_f(:), sky2_f(:)
      real(dp), intent(in) :: cdelt1_f(:), cdelt2_f(:)

      real(dp) :: crval1_f(max_inputs), crval2_f(max_inputs)
      integer :: group_id(max_inputs), group_planes(max_inputs)
      integer :: group_first(max_inputs)
      integer :: n_groups, i, j, this_group, ref_group, best_group
      integer :: unit, blocksize, fitsstat, status
      character(len=68) :: comment
      character(len=8) :: axstr
      real(dp), parameter :: tol = 1.0d-9

      do i = 1, n_inputs
         fitsstat = 0
         call safe_ftopen(unit, trim(infiles(i)), 0, blocksize, fitsstat)
         if (fitsstat.ne.0) then
            call free_fits_unit(unit)
            return ! advisory only -- silently skip on any read trouble
         endif
         write(axstr,'(I0)') sky1_f(i)
         call FTGKYD(unit, 'CRVAL'//trim(axstr), crval1_f(i), comment, fitsstat)
         write(axstr,'(I0)') sky2_f(i)
         call FTGKYD(unit, 'CRVAL'//trim(axstr), crval2_f(i), comment, fitsstat)
         call safe_ftclos(unit, fitsstat)
         if (fitsstat.ne.0) return
      enddo

      ! Group by matching (CRVAL1, CRVAL2, CDELT1, CDELT2) within tol.
      n_groups = 0
      do i = 1, n_inputs
         this_group = 0
         do j = 1, n_groups
            if (abs(crval1_f(i)-crval1_f(group_first(j))).lt.tol .and.&
            &abs(crval2_f(i)-crval2_f(group_first(j))).lt.tol .and.&
            &abs(cdelt1_f(i)-cdelt1_f(group_first(j))).lt.tol .and.&
            &abs(cdelt2_f(i)-cdelt2_f(group_first(j))).lt.tol) then
               this_group = j
               exit
            endif
         enddo
         if (this_group.eq.0) then
            n_groups = n_groups + 1
            this_group = n_groups
            group_first(this_group) = i
            group_planes(this_group) = 0
         endif
         group_id(i) = this_group
         group_planes(this_group) = group_planes(this_group) + nfreq_f(i)
      enddo

      if (n_groups.lt.2) return ! nothing to recommend -- single grid

      ref_group = 0
      do i = 1, n_inputs
         if (trim(infiles(i)).eq.trim(reffile)) then
            ref_group = group_id(i)
            exit
         endif
      enddo
      if (ref_group.eq.0) return ! reffile not one of infiles -- skip

      best_group = 1
      do j = 2, n_groups
         if (group_planes(j).gt.group_planes(best_group)) best_group = j
      enddo

      if (best_group.ne.ref_group) then
         write(*,'(A)') 'ADVISORY: reffile= is not from the band with'//&
         &' the most total planes.'
         write(*,'(A,I0,A,I0,A,A,A)') '  Current reference band totals ',&
         &group_planes(ref_group), ' planes; ', group_planes(best_group),&
         &' planes are available from another band (e.g. ',&
         &trim(infiles(group_first(best_group))), ')'
         write(*,'(A)') '  Reprojecting within a band (same pointing/'//&
         &'pixel scale) is a near-identity transform; reprojecting'//&
         &' across bands is genuine resampling and much more'//&
         &' expensive per plane -- using the largest-plane-count band'//&
         &' as reffile= minimizes total reproject cost. Not changed'//&
         &' automatically -- reffile= also determines the output grid,'//&
         &' a choice this program leaves to you.'
      endif
   end subroutine advise_reference_band

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

   subroutine build_exact_weights_file(reffile_l, infile, used_target_grid_l,&
   &target_system_arg, target_projection_arg, target_center_lon_arg,&
   &target_center_lat_arg, target_cdelt1_arg, target_cdelt2_arg,&
   &target_final_crpix_arg, lbnd_out_d, ubnd_out_d, nx_in, ny_in, wts, status)
      !! interp=exact: builds the overlap weights of one input file onto
      !! the output grid, once, before any block is read -- every
      !! reprojection path of this program then applies them per plane.
      !! Same per-thread AST setup as the paths' own regions (each thread
      !! creates its own objects inside critical (ast_setup); see
      !! process_one_file_general); bands of output rows are shared among
      !! the threads, each writing only its own rows of the table.
      use omp_lib, only: omp_get_max_threads
      character(len=*), intent(in) :: reffile_l, infile
      logical, intent(in) :: used_target_grid_l
      character(len=*), intent(in) :: target_system_arg, target_projection_arg
      double precision, intent(in) :: target_center_lon_arg, target_center_lat_arg
      double precision, intent(in) :: target_cdelt1_arg, target_cdelt2_arg
      double precision, intent(in) :: target_final_crpix_arg(2)
      double precision, intent(in) :: lbnd_out_d(2), ubnd_out_d(2)
      integer, intent(in) :: nx_in, ny_in
      type(exact_weights_t), intent(inout) :: wts
      integer, intent(inout) :: status

      integer, parameter :: band_rows = 8
      integer :: lbnd_o(2), ubnd_o(2), ny_out, n_bands, iband, status_par
      integer :: t_status, t_wcs_ref, t_skymap_ref, t_skyframe_ref
      integer :: t_naxes_ref(max_axes), t_pixaxes_ref(2)
      integer :: t_wcs_in, t_skymap_in, t_skyframe_in
      integer :: t_naxes_in(max_axes), t_pixaxes_in(2)
      integer :: t_map_in2ref, t_sky2sky, t_map_in2sky, t_lonaxis

      if (status.ne.0) return
      lbnd_o(1) = nint(lbnd_out_d(1)); lbnd_o(2) = nint(lbnd_out_d(2))
      ubnd_o(1) = nint(ubnd_out_d(1)); ubnd_o(2) = nint(ubnd_out_d(2))
      ny_out = ubnd_o(2) - lbnd_o(2) + 1
      n_bands = (ny_out + band_rows - 1) / band_rows
      call exact_init(wts, nx_in, ny_in, lbnd_o, ubnd_o)
      status_par = 0
      !$omp parallel num_threads(max(1, omp_get_max_threads())) default(none)&
      !$omp& shared(reffile_l, infile, used_target_grid_l, target_system_arg,&
      !$omp& target_projection_arg, target_center_lon_arg, target_center_lat_arg,&
      !$omp& target_cdelt1_arg, target_cdelt2_arg, target_final_crpix_arg,&
      !$omp& wts, n_bands, ny_out, status_par)&
      !$omp& private(t_status, t_wcs_ref, t_skymap_ref, t_skyframe_ref,&
      !$omp& t_naxes_ref, t_pixaxes_ref, t_wcs_in, t_skymap_in, t_skyframe_in,&
      !$omp& t_naxes_in, t_pixaxes_in, t_map_in2ref, t_sky2sky, t_map_in2sky,&
      !$omp& t_lonaxis, iband)
      t_status = 0
      !$omp critical (ast_setup)
      call ast_begin(t_status)
      if (used_target_grid_l) then
         call build_target_skymap(target_system_arg, target_projection_arg,&
         &target_center_lon_arg, target_center_lat_arg, target_cdelt1_arg,&
         &target_cdelt2_arg, target_final_crpix_arg, t_skymap_ref, t_skyframe_ref,&
         &t_pixaxes_ref, t_status)
      else
         call load_wcs(reffile_l, t_wcs_ref, t_naxes_ref, t_status)
         call extract_sky_mapping(t_wcs_ref, t_skymap_ref, t_skyframe_ref,&
         &t_pixaxes_ref, t_status)
      endif
      call load_wcs(infile, t_wcs_in, t_naxes_in, t_status)
      call extract_sky_mapping(t_wcs_in, t_skymap_in, t_skyframe_in,&
      &t_pixaxes_in, t_status)
      call compose_pix2pix(t_skymap_in, t_skyframe_in, t_skymap_ref,&
      &t_skyframe_ref, t_map_in2ref, t_status)
      ! input pixel -> output sky, without passing through output pixel
      ! coordinates (see exact_overlap_mod's exact_build_rows)
      t_sky2sky = ast_convert(t_skyframe_in, t_skyframe_ref, ' ', t_status)
      t_map_in2sky = ast_cmpmap(t_skymap_in, t_sky2sky, .true., ' ', t_status)
      t_lonaxis = ast_geti(t_skyframe_ref, 'LonAxis', t_status)
      !$omp end critical (ast_setup)
      if (t_status.ne.0) then
         !$omp atomic write
         status_par = -1
      endif
      !$omp do schedule(dynamic, 1)
      do iband = 1, n_bands
         if (t_status.eq.0 .and. status_par.eq.0) then
            call exact_build_rows(wts, t_map_in2ref, t_skymap_ref, t_map_in2sky,&
            &t_lonaxis, (iband - 1) * band_rows + 1, min(iband * band_rows, ny_out),&
            &t_status)
            if (t_status.ne.0) then
               !$omp atomic write
               status_par = -1
            endif
         endif
      enddo
      !$omp end do
      call ast_end(t_status)
      !$omp end parallel
      if (status_par.ne.0) then
         write(*,*) 'ERROR: failed to build the interp=exact overlap weights for: ',&
         &trim(infile)
         status = -1
         return
      endif
      write(*,'(A)') 'interp=exact: overlap weights built, '//&
      &trim(mib_str(exact_weights_bytes(wts)))//' MiB'
   end subroutine build_exact_weights_file

   subroutine exact_budget(infile, reffile_l, used_target_grid_l, target_cdelt1_arg,&
   &target_cdelt2_arg, nx_out, ny_out, mem_safe_bytes, status)
      !! interp=exact: the weight table is held for the whole file, so its
      !! estimated size comes out of the block budget (mem_safe_bytes).
      !! Output and input pixel sizes are measured at their reference
      !! pixels, as the beam check does.
      character(len=*), intent(in) :: infile, reffile_l
      logical, intent(in) :: used_target_grid_l
      double precision, intent(in) :: target_cdelt1_arg, target_cdelt2_arg
      integer, intent(in) :: nx_out, ny_out
      integer(kind=8), intent(inout) :: mem_safe_bytes
      integer, intent(inout) :: status
      double precision :: in_cell1, in_cell2, out_cell1, out_cell2
      integer(kind=8) :: est
      integer :: cell_status
      if (status.ne.0) return
      cell_status = 0
      call fits_reference_cell_arcsec(infile, in_cell1, in_cell2, cell_status)
      if (used_target_grid_l) then
         out_cell1 = target_cdelt1_arg * 3600.0d0; out_cell2 = target_cdelt2_arg * 3600.0d0
      else
         call fits_reference_cell_arcsec(reffile_l, out_cell1, out_cell2, cell_status)
      endif
      if (cell_status.ne.0) then
         status = -1
         return
      endif
      est = exact_estimate_bytes(nx_out, ny_out, max(out_cell1, out_cell2) / min(in_cell1, in_cell2))
      write(*,'(A)') 'interp=exact: overlap weights estimated at '//&
      &trim(mib_str(est))//' MiB, of a mem_frac_ram budget of '//&
      &trim(mib_str(mem_safe_bytes))//' MiB'
      if (est.ge.mem_safe_bytes) then
         write(*,*) 'ERROR: interp=exact needs more memory for its overlap'//&
         &' weights than mem_frac_ram allows; raise mem_frac_ram or use a'//&
         &' smaller output grid.'
         status = -1
         return
      endif
      mem_safe_bytes = mem_safe_bytes - est
   end subroutine exact_budget

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
      !! Called concurrently by every OpenMP thread during the reproject
      !! path's per-thread setup -- the CFITSIO unit comes from
      !! fitsio_unit_mod's safe_ftopen (FTGIOU+FTOPEN inside one critical
      !! section, same as reproject_cubes.f90's own load_wcs; see there
      !! for why FTOPEN/FTCLOS must share the lock with FTGIOU/FTFIOU,
      !! not just the unit-number bookkeeping alone).
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
      call safe_ftopen(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to open FITS file: ', trim(filename)
         call printerror(fitsstat)
         status = -1
         call free_fits_unit(unit)
         return
      endif

      call FTGHSP(unit, nkeys, nmore, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: FTGHSP failed for ', trim(filename)
         call printerror(fitsstat)
         call safe_ftclos(unit, fitsstat)
         status = -1
         return
      endif

      call FTGISZ(unit, size(naxes), naxes, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: FTGISZ failed for ', trim(filename)
         call printerror(fitsstat)
         call safe_ftclos(unit, fitsstat)
         status = -1
         return
      endif

      fitschan = ast_fitschan(ast_null, ast_null, ' ', status)
      if (status.ne.0) then
         call safe_ftclos(unit, fitsstat)
         return
      endif

      do i = 1, nkeys
         fitsstat = 0
         call FTGREC(unit, i, card, fitsstat)
         if (fitsstat.ne.0) then
            write(*,*) 'ERROR: FTGREC failed at card ', i, ' for ', trim(filename)
            call printerror(fitsstat)
            call safe_ftclos(unit, fitsstat)
            status = -1
            return
         endif
         call ast_putfits(fitschan, card, .false., status)
      enddo
      call ast_seti(fitschan, 'Card', 1, status)

      call safe_ftclos(unit, fitsstat)
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

   ! sky_wcs_matches_target: now shared (T36, docs/dev/
   ! MULTI_BAND_TOMOGRAPHY_PLAN.md), see src/wcs_match_mod.f90. Formerly
   ! this file's own verbatim duplicate of reproject_cubes.f90's own
   ! copy of the same function.

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
   &naxes_in_l, pixaxes_in_l, lbnd_out_d, ubnd_out_d, mem_frac_ram_l,&
   &io_overlap_l, nwriters_l, interp_scheme_l, interp_params_arr_l,&
   &used_target_grid_l, target_system_arg, target_projection_arg,&
   &target_center_lon_arg,&
   &target_center_lat_arg, target_cdelt1_arg, target_cdelt2_arg,&
   &target_final_crpix_arg, status)
      use, intrinsic :: ieee_arithmetic
      use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_wtime
      character(len=*), intent(in) :: reffile_l, infile, outfile
      integer, intent(in) :: pixaxes_ref_l(2)
      integer, intent(in) :: naxes_in_l(:), pixaxes_in_l(2)
      double precision, intent(in) :: lbnd_out_d(2), ubnd_out_d(2)
      real, intent(in) :: mem_frac_ram_l
      logical, intent(in) :: io_overlap_l
      integer, intent(in) :: nwriters_l
      integer, intent(in) :: interp_scheme_l
      double precision, intent(in) :: interp_params_arr_l(2)
      ! used_target_grid_l/target_*_arg (docs/dev/BRANCH_STATUS_update-
      ! reproject.md UR-7): PLAIN VALUES, not AST object handles -- see
      ! this program's own top-level declaration comment and
      ! target_grid_mod.f90's own header comment for why no AST object
      ! may cross the OpenMP thread boundary below.
      logical, intent(in) :: used_target_grid_l
      character(len=*), intent(in) :: target_system_arg, target_projection_arg
      double precision, intent(in) :: target_center_lon_arg, target_center_lat_arg
      double precision, intent(in) :: target_cdelt1_arg, target_cdelt2_arg
      double precision, intent(in) :: target_final_crpix_arg(2)
      integer, intent(inout) :: status

      ! interp_finterp_l: Tier 2 (AST__UKERN1, REPROJECT_INTERP_
      ! OPTIONS_PLAN.md section 6d) -- deliberately a LOCAL variable,
      ! resolved from interp_scheme_l just before the parallel region
      ! below, never received as a dummy argument. Passing a
      ! procedure(),pointer as THIS subroutine's own dummy argument,
      ! shared into its own !$omp parallel region, segfaulted inside
      ! AST when this was first wired into reproject_cubes.f90 --
      ! isolated to a gfortran code-generation issue specific to a
      ! procedure-pointer dummy argument of a contains-ed procedure
      ! shared into that procedure's own parallel region, not the
      ! finterp mechanism itself (verified correct independently --
      ! see section 6b/6c/6d and scratch/build_check/
      ! test_ukern1_internal_proc_bug.f90).
      procedure(), pointer :: interp_finterp_l

      integer :: cur_slot
      logical :: write_dispatched_ok
      integer :: naxis, k, other_axes(max_axes), n_other
      integer :: other_idx(max_axes), remainder, radix
      integer :: n_planes, status_par, nthreads, nwriters_eff
      integer :: nx_out, ny_out, naxis_out, naxes_out(max_axes)
      integer :: nx_in, ny_in
      integer :: ref_unit, out_unit, fitsstat, blocksize
      integer(kind=8) :: datastart, headstart_dum, dataend_dum
      logical :: simple, extend
      integer :: beams_unit, beams_status, casambm_status, hdutype_dum
      logical :: casambm_val
      character(len=80) :: comment_dum

      integer(kind=8) :: mem_total_kb, bytes_per_plane, mem_safe_bytes
      integer(kind=8) :: block_planes64
      integer :: block_planes, n_groups, igroup, axis1_extent
      integer :: chan_start, chan_len, local_iplane
      real, allocatable, target :: block_data_in(:,:,:), block_data_out(:,:,:,:)

      integer :: t_status, t_wcs_ref, t_skymap_ref, t_skyframe_ref
      integer :: t_naxes_ref(max_axes), t_pixaxes_ref(2)
      integer :: t_wcs_in, t_skymap_in, t_skyframe_in
      integer :: t_naxes_in(max_axes), t_pixaxes_in(2)
      integer :: t_map_in2ref

      integer :: lbnd_in(2), ubnd_in(2), lbnd_o(2), ubnd_o(2), nbad
      ! T30 item 3/7: spatial/tile reproject parallelism -- each thread
      ! resamples its own row-range of the SAME plane (AST's native
      ! REGION-vs-FULL output-bounds support), rather than each thread
      ! taking a whole separate plane. See docs/dev/
      ! MULTI_BAND_TOMOGRAPHY_PLAN.md T30 item 7 for the full
      ! measurement/derivation (verified against the real cross-band
      ! WCS transform, not just a synthetic one).
      integer :: lbnd_reg(2), ubnd_reg(2), row0, row1, rows_per, rem
      real :: badval
      real(dp) :: t_stage
      ! interp=exact: overlap weights built once for this file (before the
      ! region below) and shared read-only while rows are applied
      logical :: is_exact
      type(exact_weights_t) :: exact_wts
      ! Per-thread swim-lane instrumentation (planning-doc ticket) -- see
      ! convolve_cubes.f90's own write_convolved_file for the full
      ! rationale. iblock is PRIVATE (not shared): every thread in this
      ! persistent parallel region independently, redundantly computes
      ! the identical chan_start/chan_len progression each do-while
      ! iteration (the same pattern this subroutine already uses), so
      ! iblock is incremented the same way -- no atomic/shared needed.
      integer :: iblock, tid_local
      real(dp) :: t_thread_start, t_thread_elapsed
      character(len=160) :: thread_msg

      if (status.ne.0) return
      call timer_reset_file_stages()
      call log_message('info', 'reproject', 'starting: '//trim(infile))

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
      call safe_ftinit(out_unit, trim(outfile), blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to create output file: ', trim(outfile)
         call printerror(fitsstat)
         status = -1
         call free_fits_unit(out_unit)
         return
      endif
      simple = .true.
      extend = .false.
      call FTPHPR(out_unit, simple, -32, naxis_out, naxes_out(1:naxis_out),&
      &0, 1, extend, fitsstat)

      ! used_target_grid_l (UR-7, no reffile given): write_target_wcs_
      ! cards rebuilds the final Mapping itself from the plain values
      ! and derives CTYPE/CRVAL/CRPIX/CDELT/PC/EQUINOX/RADESYS from it --
      ! see reproject_cubes.f90's own identical branch for the full
      ! rationale.
      if (used_target_grid_l) then
         call write_target_wcs_cards(target_system_arg, target_projection_arg,&
         &target_center_lon_arg,&
         &target_center_lat_arg, target_cdelt1_arg, target_cdelt2_arg,&
         &target_final_crpix_arg, 1, 2, out_unit, status)
      else
         fitsstat = 0
         call safe_ftopen(ref_unit, trim(reffile_l), 0, blocksize, fitsstat)
         call copy_axis_keywords(ref_unit, pixaxes_ref_l(1), out_unit, 1,&
         &lbnd_out_d(1)-1.0d0, status)
         call copy_axis_keywords(ref_unit, pixaxes_ref_l(2), out_unit, 2,&
         &lbnd_out_d(2)-1.0d0, status)
         call copy_sky_rotation_matrix(ref_unit, pixaxes_ref_l(1),&
         &pixaxes_ref_l(2), out_unit, status)
         call copy_wcs_system_keywords(ref_unit, out_unit, status)
         call safe_ftclos(ref_unit, fitsstat)
      endif
      fitsstat = 0
      call safe_ftopen(ref_unit, trim(infile), 0, blocksize, fitsstat)
      do k = 1, n_other
         call copy_axis_keywords(ref_unit, other_axes(k), out_unit,&
         &2+k, 0.0d0, status)
      enddo
      call copy_generic_header_match(ref_unit, out_unit, .true., .false., status)
      if (used_target_grid_l) then
         call FTPHIS(out_unit, 'match_cubes: reprojected from '//&
         &trim(infile)//' onto a '//trim(target_system_arg)//'/'//&
         &trim(target_projection_arg)//' target grid (no reffile)', fitsstat)
      else
         call FTPHIS(out_unit, 'match_cubes: reprojected from '//&
         &trim(infile)//' onto the grid of '//trim(reffile_l), fitsstat)
      endif
      call safe_ftclos(ref_unit, fitsstat)
      ! CASAMBM/BEAMS: reproject-alone (stages=reproject) never touches
      ! the beam itself -- see reproject_cubes.f90's own identical block,
      ! which this one is a verbatim port of. copy_generic_header_match
      ! above already copied the scalar CASAMBM keyword verbatim as a
      ! raw header card (only PRIMARY-header cards, though); this
      ! attaches the actual BEAMS extension HDU it refers to.
      casambm_status = 0
      call ftgkyl(out_unit, 'CASAMBM', casambm_val, comment_dum, casambm_status)
      if (casambm_status.eq.0 .and. casambm_val) then
         beams_status = 0
         call safe_ftopen(beams_unit, trim(infile), 0, blocksize, beams_status)
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
         call safe_ftclos(beams_unit, beams_status)
      endif

      ! T19 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): CFITSIO is not
      ! guaranteed thread-safe even across different file units/handles
      ! without a specific reentrant build -- confirmed the hard way, a
      ! real, timing-dependent hang on real WALLABY+EMU data. Fetch this
      ! file's own pixel-data byte offset ONCE, right here (all header
      ! writes, including the BEAMS extension copy just above, are
      ! already done), then close out_unit immediately -- CFITSIO's job
      ! for this file is done. Every pixel write from here on goes
      ! through write_one_block_raw_general (plain Fortran stream I/O,
      ! computed byte offsets), same design as reproject_cubes.f90's own
      ! write_one_block_raw and this file's own process_one_file_
      ! restricted/write_freq_block_raw -- CFITSIO is never touched
      ! concurrently because it is not touched AT ALL during the write
      ! loop.
      datastart = 0_8
      headstart_dum = 0_8
      dataend_dum = 0_8
      fitsstat = 0
      call ftghad(out_unit, headstart_dum, datastart, dataend_dum, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to get data-start offset (FTGHAD) for: ', trim(outfile)
         call printerror(fitsstat)
         status = -1
         call safe_ftclos(out_unit, fitsstat)
         return
      endif
      call safe_ftclos(out_unit, fitsstat)

      n_planes = 1
      do k = 1, n_other
         n_planes = n_planes * naxes_in_l(other_axes(k))
      enddo
      call get_mem_total_kb(mem_total_kb)
      ! block_data_out's own term is x2 (not x1) -- double-buffered for
      ! io_overlap (see below), budgeted whether or not it's on.
      bytes_per_plane = int(4,8) * (int(nx_in,8)*int(ny_in,8) +&
      &2_8*int(nx_out,8)*int(ny_out,8))
      mem_safe_bytes = int(real(mem_frac_ram_l,8) * real(mem_total_kb,8) *&
      &1024.0d0, 8)
      is_exact = interp_scheme_l.eq.interp_exact
      if (is_exact) then
         call exact_budget(infile, reffile_l, used_target_grid_l, target_cdelt1_arg,&
         &target_cdelt2_arg, nx_out, ny_out, mem_safe_bytes, status)
         if (status.ne.0) return
      endif
      block_planes64 = max(1_8, mem_safe_bytes / bytes_per_plane)
      block_planes64 = min(block_planes64, max_elements_per_block /&
      &max(1_8, max(int(nx_in,8)*int(ny_in,8), int(nx_out,8)*int(ny_out,8))))
      block_planes64 = max(1_8, block_planes64)
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
      allocate(block_data_out(nx_out, ny_out, block_planes, 0:1))

      n_groups = 1
      do k = 2, n_other
         n_groups = n_groups * naxes_in_l(other_axes(k))
      enddo

      cur_slot = 0
      write_pending = .false.
      write_failed = .false.
      status_par = 0
      nthreads = max(1, min(omp_get_max_threads(), block_planes))
      nwriters_eff = max(1, min(nwriters_l, omp_get_max_threads()))
      if (interp_scheme_l.eq.ast__ukern1) then
         interp_finterp_l => lanczos_kernel
      else
         interp_finterp_l => ast_null
      endif
      if (is_exact) then
         call timer_start(t_stage)
         call build_exact_weights_file(reffile_l, infile, used_target_grid_l,&
         &target_system_arg, target_projection_arg, target_center_lon_arg,&
         &target_center_lat_arg, target_cdelt1_arg, target_cdelt2_arg,&
         &target_final_crpix_arg, lbnd_out_d, ubnd_out_d, nx_in, ny_in, exact_wts, status)
         call timer_stop('exact_weights', t_stage)
         if (status.ne.0) return
      endif
      !$omp parallel num_threads(nthreads) default(none)&
      !$omp& shared(infile, reffile_l, outfile, naxes_in_l, pixaxes_in_l,&
      !$omp& other_axes, n_other, lbnd_out_d, ubnd_out_d, datastart, status_par,&
      !$omp& nx_in, ny_in, nx_out, ny_out, block_planes, block_data_in,&
      !$omp& block_data_out, n_groups, axis1_extent, cur_slot, io_overlap_l,&
      !$omp& write_dispatched_ok, write_pending, write_thread_id,&
      !$omp& write_failed, write_job_general, t_stage, nwriters_eff, nthreads,&
      !$omp& interp_scheme_l, interp_params_arr_l, interp_finterp_l,&
      !$omp& used_target_grid_l, target_system_arg, target_projection_arg,&
      !$omp& target_center_lon_arg,&
      !$omp& target_center_lat_arg, target_cdelt1_arg, target_cdelt2_arg,&
      !$omp& target_final_crpix_arg, is_exact, exact_wts)&
      !$omp& private(t_status, t_wcs_ref, t_skymap_ref, t_skyframe_ref,&
      !$omp& t_naxes_ref, t_pixaxes_ref, t_wcs_in, t_skymap_in, t_skyframe_in,&
      !$omp& t_naxes_in, t_pixaxes_in, t_map_in2ref, other_idx, remainder,&
      !$omp& radix, k, igroup, chan_start, chan_len, local_iplane, nbad,&
      !$omp& lbnd_in, ubnd_in, lbnd_o, ubnd_o, lbnd_reg, ubnd_reg, row0, row1,&
      !$omp& rows_per, rem, badval,&
      !$omp& iblock, tid_local, t_thread_start, t_thread_elapsed, thread_msg)

      ! T20 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): each thread still
      ! builds its OWN private Mapping (AST enforces per-thread object
      ! ownership at runtime, and this Fortran binding exports neither
      ! astLock nor astUnlock to hand one off -- see reproject_cubes.
      ! f90's own write_reprojected_file for the full account of why a
      ! build-once-and-copy design was tried and reverted). The `!$omp
      ! critical` below only serializes the CREATION step (load_wcs/
      ! extract_sky_mapping/compose_pix2pix) -- the confirmed
      ! (via gdb) cause of the real WALLABY+EMU deadlock was 6 threads
      ! calling ast_read_ simultaneously, contending on
      ! libstarlink_ast.so.9's own internal object-creation bookkeeping,
      ! not object ownership. Each thread still ends up with its own,
      ! properly-owned object; the actual per-plane resampling (ast_
      ! resampler, below) stays parallel and untouched -- it was
      ! never the problem.
      iblock = 0
      t_status = 0
      !$omp critical (ast_setup)
      call ast_begin(t_status)
      if (used_target_grid_l) then
         call build_target_skymap(target_system_arg, target_projection_arg,&
         &target_center_lon_arg,&
         &target_center_lat_arg, target_cdelt1_arg, target_cdelt2_arg,&
         &target_final_crpix_arg, t_skymap_ref, t_skyframe_ref,&
         &t_pixaxes_ref, t_status)
      else
         call load_wcs(reffile_l, t_wcs_ref, t_naxes_ref, t_status)
         call extract_sky_mapping(t_wcs_ref, t_skymap_ref, t_skyframe_ref,&
         &t_pixaxes_ref, t_status)
      endif
      call load_wcs(infile, t_wcs_in, t_naxes_in, t_status)
      call extract_sky_mapping(t_wcs_in, t_skymap_in, t_skyframe_in,&
      &t_pixaxes_in, t_status)
      call compose_pix2pix(t_skymap_in, t_skyframe_in, t_skymap_ref,&
      &t_skyframe_ref, t_map_in2ref, t_status)
      !$omp end critical (ast_setup)
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
            call timer_start(t_stage)
            call read_one_block(infile, naxes_in_l, pixaxes_in_l, other_axes,&
            &other_idx, n_other, chan_start, chan_len, nx_in, ny_in,&
            &block_data_in(:,:,1:chan_len), status_par)
            call timer_stop('block_read', t_stage)
            call timer_start(t_stage)
            !$omp end single

            iblock = iblock + 1
            tid_local = omp_get_thread_num()
            ! T30 item 3/7: spatial/tile parallelism -- every thread
            ! processes every plane (no worksharing over local_iplane
            ! any more), each working its own row-range of the OUTPUT
            ! footprint via AST's native REGION-vs-FULL output-bounds
            ! support. Verified no worse than plane-parallel when
            ! chan_len is a multiple of nthreads, win otherwise --
            ! see docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md T30 item 7.
            ! Explicit barrier at the end of each plane iteration:
            ! block_data_out(:,:,local_iplane,cur_slot) must be
            ! written by every thread before the loop moves on.
            rows_per = ny_out / nthreads
            rem = mod(ny_out, nthreads)
            if (tid_local < rem) then
               row0 = tid_local*(rows_per+1) + 1
               row1 = row0 + rows_per
            else
               row0 = rem*(rows_per+1) + (tid_local-rem)*rows_per + 1
               row1 = row0 + rows_per - 1
            endif
            ! row1 < row0 here means nthreads > ny_out -- more threads
            ! than output rows to hand out (confirmed live: 36 OMP
            ! threads against this project's own 32-row test fixture,
            ! ny_out=32, gives every thread from 32 onward row0=33 >
            ! row1=32, an empty region AST's own resampler correctly
            ! rejects rather than silently accepting). This thread
            ! simply has no rows assigned; skip its resample work for
            ! every plane rather than feeding astResampleR a malformed
            ! region. Loop/barrier structure below is unchanged so
            ! every thread still reaches the same barriers the same
            ! number of times.
            do local_iplane = 1, chan_len
               t_thread_start = omp_get_wtime()
               if (t_status.eq.0 .and. status_par.eq.0 .and. row1.ge.row0) then
                  lbnd_in(1) = 1
                  lbnd_in(2) = 1
                  ubnd_in(1) = nx_in
                  ubnd_in(2) = ny_in
                  lbnd_o(1) = nint(lbnd_out_d(1))
                  lbnd_o(2) = nint(lbnd_out_d(2))
                  ubnd_o(1) = nint(ubnd_out_d(1))
                  ubnd_o(2) = nint(ubnd_out_d(2))
                  lbnd_reg(1) = lbnd_o(1)
                  lbnd_reg(2) = lbnd_o(2) + row0 - 1
                  ubnd_reg(1) = ubnd_o(1)
                  ubnd_reg(2) = lbnd_o(2) + row1 - 1
                  badval = ieee_value(badval, ieee_quiet_nan)
                  if (is_exact) then
                     call exact_apply_rows(exact_wts, block_data_in(:,:,local_iplane),&
                     &block_data_out(:,:,local_iplane,cur_slot), row0, row1)
                  else
                  nbad = ast_resampler(t_map_in2ref, 2, lbnd_in, ubnd_in,&
                  &block_data_in(:,:,local_iplane),&
                  &block_data_in(:,:,local_iplane),&
                  &interp_scheme_l, interp_finterp_l, interp_params_arr_l, 0, 0.0d0, 100, badval,&
                  &2, lbnd_o, ubnd_o, lbnd_reg, ubnd_reg,&
                  &block_data_out(:,:,local_iplane,cur_slot),&
                  &block_data_out(:,:,local_iplane,cur_slot), t_status)
                  endif
                  if (t_status.ne.0) then
                     !$omp atomic write
                     status_par = -1
                  endif
               endif
               !$omp barrier
               t_thread_elapsed = (omp_get_wtime() - t_thread_start) * 1000.0_dp
               ! T30 item 7 follow-up (2026-08-09): per-PLANE, not
               ! per-block, team timing -- every thread in the team
               ! works the SAME plane together (row-split), so every
               ! thread logs its OWN line, all sharing the identical
               ! team-elapsed dur_ms for that one plane -- keeps a
               ! swim-lane plot honest (nthreads simultaneous lanes for
               ! this one plane, not 1 -- see scripts/
               ! plot_tile_async_swimlane.py's THREAD_CPU_RE comment).
               ! Matches convolve's own per-plane 'block=...plane=...
               ! nthreads=...' shape exactly; replaces the old per-block
               ! start/done pair (redundant with this).
               write(thread_msg,'(A,I0,A,I0,A,I0,A,I0,A,F10.3)') 'thread_timing stage=resample event=done tid=',&
               &tid_local,' block=',iblock,' plane=',chan_start+local_iplane-1,&
               &' nthreads=',nthreads,' dur_ms=',t_thread_elapsed
               call log_message('debug','tile_thread',trim(thread_msg))
            enddo

            !$omp single
            call timer_stop('block_resample', t_stage)

            call timer_start(t_stage)
            if (write_pending) then
               call block_write_join(write_thread_id)
               write_pending = .false.
               if (write_failed) status_par = -1
            endif
            call timer_stop('block_write_join', t_stage)

            call timer_start(t_stage)
            if (status_par.eq.0) then
               write_job_general%file_path = outfile
               write_job_general%datastart = datastart
               write_job_general%naxes_in(1:max_axes) = naxes_in_l(1:max_axes)
               write_job_general%other_axes(1:max_axes) = other_axes(1:max_axes)
               write_job_general%other_idx(1:max_axes) = other_idx(1:max_axes)
               write_job_general%n_other = n_other
               write_job_general%chan_start = chan_start
               write_job_general%chan_len = chan_len
               write_job_general%nwriters_eff = max(1, min(nwriters_eff, chan_len))
               write_job_general%nx = nx_out
               write_job_general%ny = ny_out
               write_job_general%data => block_data_out(:,:,1:chan_len,cur_slot)
               if (io_overlap_l) then
                  call block_write_dispatch_async_general(write_job_general,&
                  &write_thread_id, write_dispatched_ok)
                  write_pending = write_dispatched_ok
               else
                  call do_block_write_general(write_job_general)
                  if (write_failed) status_par = -1
               endif
            endif
            call timer_stop('block_write', t_stage)
            cur_slot = 1 - cur_slot
            !$omp end single

            chan_start = chan_start + chan_len
         enddo
      enddo

      call ast_end(t_status)
      !$omp end parallel

      call timer_start(t_stage)
      if (write_pending) then
         call block_write_join(write_thread_id)
         write_pending = .false.
         if (write_failed) status_par = -1
      endif
      call timer_stop('block_write_join', t_stage)

      deallocate(block_data_in)
      deallocate(block_data_out)
      if (is_exact) call exact_free(exact_wts)

      if (status_par.ne.0) then
         write(*,*) 'ERROR: failed to resample/write one or more planes for: ',&
         &trim(infile)
         status = -1
         return
      endif

      call log_message('info', 'reproject', 'finished: '//trim(infile))
      call timer_report_file_summary(infile)
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
      group = 1
      badval = ieee_value(badval, ieee_quiet_nan)
      call safe_ftopen(unit, trim(filename), 0, blocksize, fitsstat)
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
      call safe_ftclos(unit, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to read block (planes ', chan_start,&
         &'-', chan_start+chan_len-1, ') from ', trim(filename)
         call printerror(fitsstat)
         status = -1
      endif
   end subroutine read_one_block

   subroutine write_one_block_raw_general(file_path, datastart, naxes_in_l,&
   &other_axes, other_idx, n_other, chan_start, chan_len, nx_out, ny_out,&
   &block_data_out, status)
      !! Raw Fortran stream I/O at a computed byte offset -- CFITSIO is
      !! never touched here at all (T19, docs/dev/MULTI_BAND_
      !! TOMOGRAPHY_PLAN.md). Verbatim port of reproject_cubes.f90's own
      !! write_one_block_raw -- see there for the full byte-offset/
      !! contiguity rationale (this file's general-axis output layout,
      !! sky at output axes 1,2 full extent, other_axes(1) at output
      !! axis 3 spanning chan_start:chan_start+chan_len-1, every slower
      !! other axis fixed at other_idx(2:n_other), is identical). Reuses
      !! this file's own host_is_big_endian_mc/swap_bytes_r4_inplace_mc
      !! (process_one_file_restricted's own helpers, further down) --
      !! not duplicated here.
      character(len=*), intent(in) :: file_path
      integer(kind=8), intent(in) :: datastart
      integer, intent(in) :: naxes_in_l(:), other_axes(:), other_idx(:)
      integer, intent(in) :: n_other, chan_start, chan_len, nx_out, ny_out
      real, intent(in) :: block_data_out(:,:,:)
      integer, intent(inout) :: status

      integer :: u, ios, ip, k
      logical :: need_swap
      integer(kind=8) :: plane_elems, stride, fixed_offset, byte_pos
      real, allocatable :: plane_buf(:,:)

      if (status.ne.0) return

      need_swap = .not. host_is_big_endian_mc()
      plane_elems = int(nx_out,8) * int(ny_out,8)

      stride = plane_elems
      if (n_other.ge.1) stride = stride * int(naxes_in_l(other_axes(1)),8)
      fixed_offset = 0_8
      do k = 2, n_other
         fixed_offset = fixed_offset + int(other_idx(k)-1,8) * stride
         stride = stride * int(naxes_in_l(other_axes(k)),8)
      enddo

      open(newunit=u, file=trim(file_path), access='stream',&
      &form='unformatted', status='old', action='write', iostat=ios)
      if (ios.ne.0) then
         write(*,*) 'ERROR: write_one_block_raw_general: failed to open ',&
         &trim(file_path)
         status = -1
         return
      endif

      allocate(plane_buf(nx_out,ny_out))
      do ip = 1, chan_len
         plane_buf = block_data_out(:,:,ip)
         if (need_swap) call swap_bytes_r4_inplace_mc(plane_buf, plane_elems)
         byte_pos = datastart +&
         &(int(chan_start-1+ip-1,8)*plane_elems + fixed_offset)*4_8 + 1_8
         write(u, pos=byte_pos, iostat=ios) plane_buf
         if (ios.ne.0) then
            write(*,*) 'ERROR: write_one_block_raw_general: write failed (plane ',&
            &chan_start+ip-1, ') for ', trim(file_path)
            status = -1
            exit
         endif
      enddo
      deallocate(plane_buf)
      close(u)
   end subroutine write_one_block_raw_general

   subroutine do_block_write_general(job)
      !! T4d-style: verbatim port of reproject_cubes.f90's own
      !! do_block_write -- see convolve_cubes.f90's write_convolved_file
      !! for the full single-writer-at-a-time rationale.
      type(block_write_job_general_t), intent(inout) :: job
      integer :: status_local, wk, wbase, wrem, wlen_k, woff_k, wchan_start_k
      logical :: any_failed

      any_failed = .false.
      if (job%nwriters_eff.le.1) then
         status_local = 0
         call write_one_block_raw_general(trim(job%file_path), job%datastart,&
         &job%naxes_in(1:max_axes), job%other_axes(1:max_axes),&
         &job%other_idx(1:max_axes), job%n_other, job%chan_start, job%chan_len,&
         &job%nx, job%ny, job%data, status_local)
         if (status_local.ne.0) any_failed = .true.
      else
         ! nwriters_eff>1 (T21, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md):
         ! same even-split scheme as do_block_write's own restricted-path
         ! copy -- splits along other_axes(1) (the axis varying within
         ! this block), every slower axis already fixed at other_idx(:)
         ! for the whole job, so each sub-chunk is still a plain,
         ! contiguous other_axes(1) sub-range.
         wbase = job%chan_len / job%nwriters_eff
         wrem = mod(job%chan_len, job%nwriters_eff)
         !$omp parallel do num_threads(job%nwriters_eff) default(none)&
         !$omp& shared(job, wbase, wrem, any_failed)&
         !$omp& private(wk, wlen_k, woff_k, wchan_start_k, status_local)
         do wk = 0, job%nwriters_eff-1
            wlen_k = wbase + merge(1, 0, wk.lt.wrem)
            woff_k = wk*wbase + min(wk, wrem)
            wchan_start_k = job%chan_start + woff_k
            status_local = 0
            if (wlen_k.gt.0) then
               call write_one_block_raw_general(trim(job%file_path),&
               &job%datastart, job%naxes_in(1:max_axes),&
               &job%other_axes(1:max_axes), job%other_idx(1:max_axes),&
               &job%n_other, wchan_start_k, wlen_k, job%nx, job%ny,&
               &job%data(:,:,woff_k+1:woff_k+wlen_k), status_local)
               if (status_local.ne.0) then
                  !$omp atomic write
                  any_failed = .true.
               endif
            endif
         enddo
         !$omp end parallel do
      endif
      if (any_failed) then
         write(*,*) 'ERROR: background write failed for planes ',&
         &job%chan_start, '-', job%chan_start+job%chan_len-1
         write_failed = .true.
      endif
   end subroutine do_block_write_general

   function block_write_thread_entry_general(arg) bind(C) result(res)
      type(c_ptr), value :: arg
      type(c_ptr) :: res
      type(block_write_job_general_t), pointer :: job

      call c_f_pointer(arg, job)
      call do_block_write_general(job)
      res = c_null_ptr
   end function block_write_thread_entry_general

   subroutine block_write_dispatch_async_general(job, thread_id, dispatched)
      type(block_write_job_general_t), intent(inout), target :: job
      integer(c_long), intent(out) :: thread_id
      logical, intent(out) :: dispatched
      integer(c_int) :: rc

      rc = c_pthread_create(thread_id, c_null_ptr,&
      &c_funloc(block_write_thread_entry_general), c_loc(job))
      if (rc.ne.0_c_int) then
         write(*,*) 'WARNING: pthread_create failed for async block write;'//&
         &' running inline'
         call do_block_write_general(job)
         dispatched = .false.
      else
         dispatched = .true.
      endif
   end subroutine block_write_dispatch_async_general

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
   &mem_frac_ram_l, io_overlap_l, nwriters_l, interp_scheme_l,&
   &interp_params_arr_l, used_target_grid_l, target_system_arg,&
   &target_projection_arg,&
   &target_center_lon_arg, target_center_lat_arg, target_cdelt1_arg,&
   &target_cdelt2_arg, target_final_crpix_arg, status)
      use, intrinsic :: ieee_arithmetic
      use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_wtime
      use gaussft_mod, only: plan_convolution, convolve_to_beam,&
      &destroy_convolution_plan, next_fast_fft_size
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
      logical, intent(in) :: io_overlap_l
      integer, intent(in) :: nwriters_l
      integer, intent(in) :: interp_scheme_l
      double precision, intent(in) :: interp_params_arr_l(2)
      ! used_target_grid_l/target_*_arg (docs/dev/BRANCH_STATUS_update-
      ! reproject.md UR-7): PLAIN VALUES, not AST object handles -- see
      ! process_one_file_general's own identical declaration comment.
      logical, intent(in) :: used_target_grid_l
      character(len=*), intent(in) :: target_system_arg, target_projection_arg
      double precision, intent(in) :: target_center_lon_arg, target_center_lat_arg
      double precision, intent(in) :: target_cdelt1_arg, target_cdelt2_arg
      double precision, intent(in) :: target_final_crpix_arg(2)
      integer, intent(out) :: status
      integer :: nwriters_eff
      ! interp_finterp_l: Tier 2 (AST__UKERN1, REPROJECT_INTERP_
      ! OPTIONS_PLAN.md section 6d) -- LOCAL variable, same rationale
      ! as process_one_file_general's own identical comment. Resolved
      ! once below, before the do-while block loop, reused by
      ! whichever of this subroutine's two !$omp parallel regions
      ! (convolve_reproject/reproject_convolve order) runs.
      procedure(), pointer :: interp_finterp_l

      integer :: nx_in, ny_in, nx_out, ny_out, nx_pad, ny_pad, conv_nx, conv_ny
      integer :: in_unit, ref_unit, out_unit, fitsstat, blocksize
      logical :: simple, extend
      integer :: naxes_out(max_axes)
      ! T19: this file's own pixel-data byte offset (FTGHAD, fetched
      ! once), for write_freq_block_raw -- see the comment right before
      ! its own use, below, for why out_unit is closed immediately after
      ! header setup instead of kept open for the whole block loop.
      integer(kind=8) :: datastart, headstart_dum, dataend_dum
      integer(kind=8) :: mem_total_kb, bytes_per_plane, mem_safe_bytes, block_planes64
      integer(kind=8) :: transient_bytes_per_plane, mem_safe_bytes_io
      integer :: block_planes, chan_start, chan_len, local_iplane, nthreads
      integer :: omp_threads_cap
      ! interp=exact: overlap weights built once for this file, before the
      ! block loop, and shared read-only by each block's reproject region
      logical :: is_exact
      type(exact_weights_t) :: exact_wts
      integer :: cur_slot
      logical :: write_dispatched_ok
      real, allocatable, target :: block_in(:,:,:), block_out(:,:,:,:)
      real(dp) :: bpa_in_pixel(nfreq), tgt_bpa_pixel_native, tgt_bpa_pixel_out
      real(dp) :: ref_cdelt1, ref_cdelt2
      integer(kind=8) :: plan_fwd, plan_bwd
      integer :: status_par, ich, k
      real(dp) :: dx_deg, dy_deg, ref_dx_deg, ref_dy_deg
      real(dp) :: nanval
      ! Per-thread swim-lane instrumentation (planning-doc ticket) -- see
      ! convolve_cubes.f90's own write_convolved_file for the full
      ! rationale. iblock lives in the SERIAL outer do-while scope here
      ! (this parallel region opens/closes fresh each block, unlike
      ! process_one_file_general's persistent one), so it is shared and
      ! read-only inside the region.
      integer :: iblock, tid_local
      real(dp) :: t_thread_start, t_thread_elapsed
      character(len=160) :: thread_msg
      ! T22 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): per-thread private
      ! start times for convolve-vs-reproject sub-stage timing, split
      ! out of the single lumped 'block_convolve' timer below (which
      ! wraps BOTH convolve_to_beam and ast_resampler for stages=both,
      ! giving no way to tell which one dominates a real run's own
      ! wall-clock time -- exactly what motivated this).
      real(dp) :: t_conv_local, t_resamp_local
      ! ALLOCATABLE, not automatic/stack -- see convolve_cubes.f90's own
      ! write_convolved_file's identical comment: for a real full-size
      ! image these are ~160MB (real(dp)) or ~80MB (real) EACH, and as
      ! automatic arrays private to each OMP thread they were silently
      ! stack-allocated, guaranteed to overflow a worker thread's default
      ! stack on real production-scale data (found via a ~46GB
      ! end-to-end verification run; invisible on this project's
      ! own tiny 32x32 test fixtures). Allocated on first use per thread
      ! below (nx_in/ny_in/nx_out/ny_out are fixed for this whole call, so
      ! one allocation per thread suffices).
      real(dp), allocatable :: plane_native(:,:), plane_out_arr(:,:)
      ! astResampleR (ast_resampler) is a REAL*4 Fortran interface (matching
      ! reproject_cubes.f90's own block_data_in/block_data_out, both plain
      ! `real`, never real(dp)) -- these two single-precision scratch
      ! arrays exist to hand it REAL*4 storage in the
      ! convolve_reproject order, where the plane being resampled comes
      ! from a real(dp) convolution result, not directly from a block_in/
      ! block_out slice the way reproject_convolve order can use as-is.
      real, allocatable :: plane_native_sp(:,:), plane_out_sp(:,:)

      integer :: t_status, t_wcs_ref, t_skymap_ref, t_skyframe_ref
      integer :: t_naxes_ref(max_axes), t_pixaxes_ref(2)
      integer :: t_wcs_in, t_skymap_in, t_skyframe_in
      integer :: t_naxes_in2(max_axes), t_pixaxes_in(2)
      integer :: t_map_in2ref
      integer :: lbnd_in(2), ubnd_in(2), lbnd_o(2), ubnd_o(2), nbad
      ! T30 item 3/7: spatial/tile reproject parallelism -- see
      ! process_one_file_general's own identical comment.
      integer :: lbnd_reg(2), ubnd_reg(2), row0, row1, rows_per, rem
      real :: badval_sp
      real(dp) :: t_stage

      status = 0
      call timer_reset_file_stages()
      call log_message('info', 'convolve', 'starting: '//trim(infile))
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
      ! needed for the selected order gets computed, the
      ! other stays at a harmless default when do_reproject_l is false.
      call sky_to_pixel_bpa(tgt_bpa, cdelt1, cdelt2, tgt_bpa_pixel_native)
      tgt_bpa_pixel_out = tgt_bpa_pixel_native
      ref_dx_deg = dx_deg
      ref_dy_deg = dy_deg
      if (do_reproject_l) then
         if (used_target_grid_l) then
            ! No reffile to reopen -- target_cdelt1_arg/target_cdelt2_arg
            ! are already the output grid's own pixel-scale MAGNITUDES
            ! (UR-7); apply the same sign convention build_trial_wcs
            ! (target_grid_mod.f90) itself used to build that grid
            ! (longitude axis negative, latitude axis positive), so
            ! sky_to_pixel_bpa below sees the same convention read_ref_
            ! cdelt's raw CDELT1/CDELT2 would have given it.
            ref_cdelt1 = -target_cdelt1_arg
            ref_cdelt2 = target_cdelt2_arg
         else
            call read_ref_cdelt(reffile_l, pixaxes_ref_l, ref_cdelt1, ref_cdelt2, status)
            if (status.ne.0) return
         endif
         ref_dx_deg = abs(ref_cdelt1)
         ref_dy_deg = abs(ref_cdelt2)
         call sky_to_pixel_bpa(tgt_bpa, ref_cdelt1, ref_cdelt2, tgt_bpa_pixel_out)
      endif

      fitsstat = 0
      call safe_ftopen(in_unit, trim(infile), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot reopen input for output header: ', trim(infile)
         status = -1
         call free_fits_unit(in_unit)
         return
      endif

      fitsstat = 0
      ! Plain filename (NOT '!'-prefixed): FTINIT fails if outfile already
      ! exists rather than silently clobbering it. Fixed bug, found while
      ! implementing the skip-if-already-matched feature (planning-doc
      ! ticket) and confirmed to be pervasive by auditing every FTINIT
      ! call site in this project for the same pattern (also found and
      ! fixed in process_one_file_general's own caller just above in
      ! this file, in reproject_cubes.f90's own equivalent call site,
      ! and in convolve_cubes.f90's own internal FTINIT): this call
      ! previously used the '!'-prefix CLOBBER convention, silently
      ! deleting and overwriting a pre-existing output with no warning --
      ! inconsistent with this project's own standing rule that a
      ! pre-existing output path is always refused, never silently
      ! reused or overwritten.
      call safe_ftinit(out_unit, trim(outfile), blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot create output file (already exists?): ',&
         &trim(outfile)
         status = -1
         call safe_ftclos(in_unit, fitsstat)
         call free_fits_unit(out_unit)
         return
      endif

      naxes_out(1:naxis) = naxes(1:naxis)
      naxes_out(sky1) = nx_out
      naxes_out(sky2) = ny_out
      simple = .true.
      extend = .false.
      call FTPHPR(out_unit, simple, -32, naxis, naxes_out(1:naxis), 0, 1, extend, fitsstat)

      if (do_reproject_l) then
         ! used_target_grid_l (UR-7, no reffile given): see process_one_
         ! file_general's own identical branch for the full rationale.
         ! dst_axis1/dst_axis2 are sky1/sky2 here (NOT always 1,2) --
         ! this subroutine's own "restricted" scope keeps each axis in
         ! its original input position, unlike process_one_file_
         ! general/reproject_cubes.f90's own always-1,2 convention.
         if (used_target_grid_l) then
            call write_target_wcs_cards(target_system_arg, target_projection_arg,&
         &target_center_lon_arg,&
            &target_center_lat_arg, target_cdelt1_arg, target_cdelt2_arg,&
            &target_final_crpix_arg, sky1, sky2, out_unit, status)
         else
            fitsstat = 0
            call safe_ftopen(ref_unit, trim(reffile_l), 0, blocksize, fitsstat)
            call copy_axis_keywords(ref_unit, sky1, out_unit, sky1,&
            &lbnd_out_d(1)-1.0d0, status)
            call copy_axis_keywords(ref_unit, sky2, out_unit, sky2,&
            &lbnd_out_d(2)-1.0d0, status)
            call copy_sky_rotation_matrix(ref_unit, sky1, sky2, out_unit, status)
            call copy_wcs_system_keywords(ref_unit, out_unit, status)
            call safe_ftclos(ref_unit, fitsstat)
         endif
         call copy_axis_keywords(in_unit, freq_axis, out_unit, freq_axis, 0.0d0, status)
         ! Any OTHER axis (degenerate by this scope's own requirement --
         ! see read_axis_info -- e.g. a size-1 STOKES axis) also needs its
         ! own WCS keywords copied through explicitly: exclude_axis_wcs=
         ! true below skips ALL axis-indexed keywords in the generic copy,
         ! not just sky1/sky2/freq_axis, so without this loop a degenerate
         ! axis's CTYPE/CRVAL/CRPIX/CDELT would be lost. Caught by
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
      call safe_ftclos(in_unit, fitsstat)

      ! T19 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): CFITSIO is not
      ! guaranteed thread-safe even across different file units/handles
      ! without a specific reentrant build -- confirmed the hard way, a
      ! real, timing-dependent hang on real WALLABY+EMU data: the main
      ! thread's own next-block read (read_freq_block, CFITSIO) blocked
      ! concurrently with the background write thread's own write
      ! (formerly FTPSSE, also CFITSIO), both on what is almost
      ! certainly an internal CFITSIO lock -- not an environment/
      ! systemd issue (confirmed absent even with no such wrapper in
      ! use at all). Fetch this file's own pixel-data byte offset ONCE,
      ! right here, then close out_unit immediately: CFITSIO's job for
      ! this file is done until the BEAMS table gets appended, long
      ! after every block write/join has completed (see the end of this
      ! subroutine) -- no concurrency risk there. Every pixel write
      ! between now and then goes through write_freq_block_raw (plain
      ! Fortran stream I/O, computed byte offsets), exactly mirroring
      ! rm_synthesis_mod.f90's own io_write_threads>1 design
      ! (write_rm_chunk_raw) -- CFITSIO is never touched concurrently
      ! because it is not touched AT ALL during the write loop.
      datastart = 0_8
      headstart_dum = 0_8
      dataend_dum = 0_8
      fitsstat = 0
      call ftghad(out_unit, headstart_dum, datastart, dataend_dum, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to get data-start offset (FTGHAD) for: ', trim(outfile)
         call printerror(fitsstat)
         status = -1
         call safe_ftclos(out_unit, fitsstat)
         return
      endif
      call safe_ftclos(out_unit, fitsstat)

      ! Convolution happens on the native grid for order=convolve_reproject
      ! (or always, for stages=convolve alone, where nx_out/ny_out were
      ! already set equal to nx_in/ny_in above) -- but on the OUTPUT grid
      ! for order=reproject_convolve, which can be a different
      ! size (footprint_mode=intersection/union crops or grows it).
      ! Determined here (moved ahead of the block/memory budget below,
      ! which needs nx_pad/ny_pad) rather than right before
      ! plan_convolution's own call further down.
      if (.not. do_reproject_l .or. convolve_first_l) then
         conv_nx = nx_in
         conv_ny = ny_in
      else
         conv_nx = nx_out
         conv_ny = ny_out
      endif
      nx_pad = next_fast_fft_size(conv_nx)
      ny_pad = next_fast_fft_size(conv_ny)

      call get_mem_total_kb(mem_total_kb)
      ! block_out's own term is x2 (not x1) -- double-buffered for
      ! io_overlap (see below), budgeted whether or not it's on.
      bytes_per_plane = int(4,8) * (int(nx_in,8)*int(ny_in,8) +&
      &2_8*int(nx_out,8)*int(ny_out,8))
      mem_safe_bytes = int(real(mem_frac_ram_l,8) * real(mem_total_kb,8) * 1024.0d0, 8)

      ! T26 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): the block I/O buffer
      ! (bytes_per_plane above) is NOT the only memory a plane in flight
      ! needs -- convolve_to_beam's own NaN-aware path (src/gaussft.f90)
      ! allocates several full-size temporaries of its own (a native-size
      ! NaN mask, two padded-size complex FFT buffers alive at once, two
      ! native-size real buffers), on top of this caller's own native-size
      ! plane_native array (real(dp)) held for the whole call. Confirmed
      ! empirically on the real WALLABY+EMU run: VmPeak hit 71.65GB against a
      ! 62.7GB host with block_planes=14 and OMP_NUM_THREADS=6, entirely
      ! from this gap (system swap climbed to 12GB).
      !
      ! 16*(conv_nx*conv_ny) + 16*(nx_pad*ny_pad) bytes, PER THREAD
      ! actively convolving a plane (thread-private locals, not shared).
      ! The two 16s are NOT the same "16" -- each is its own sum, per
      ! T27 (gaussft_mod's FFT buffers/intermediates now single
      ! precision internally):
      !   native-size term  = nan_mask(4, logical) + c_d-or-c_m(4,
      !                        real(sp)) + this caller's own
      !                        plane_native(8, real(dp) -- deliberately
      !                        NOT converted, so convolve_to_beam's
      !                        external interface stays unchanged) = 16
      !   padded-size term  = g_final(8, complex(sp)) + cimg-or-cmsk(8,
      !                        complex(sp)) = 16
      ! The padded term halved (32->16, both contributors
      ! converted to sp); the native term only partly halved (20->16),
      ! held back by plane_native's own real(dp) -- converting that too
      ! (a caller-side change, not just gaussft_mod's own internals)
      ! would drop the native term further, to 8. NOT computed from
      ! gaussft_mod's own type declarations, so this drifts stale if
      ! those kinds change again (as it did once already across T27,
      ! caught and re-derived here) -- re-derive by hand, don't assume.
      !
      ! T28 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): reserved ONCE now,
      ! not omp_get_max_threads() copies of it. The per-plane loop below no longer
      ! runs convolve_to_beam concurrently across planes -- it's called
      ! serially, once per plane, with each call internally
      ! multi-threaded via plan_convolution's own nthreads argument
      ! (sfftw_plan_with_nthreads) instead. So this is a per-PLANE cost,
      ! not a per-THREAD one -- T26's own decrementing search (find the
      ! largest safe thread count) no longer applies to THIS
      ! reservation; renamed accordingly (transient_bytes_per_plane, not
      ! _per_thread). omp_threads_cap below is a plain performance
      ! choice now (how many threads help each plane's own convolve,
      ! and separately, how many planes' worth of reproject work runs
      ! concurrently -- reproject KEPT its original parallel-across-
      ! planes strategy, see the per-plane loop below), not a
      ! memory-safety cap.
      transient_bytes_per_plane = 16_8*int(conv_nx,8)*int(conv_ny,8) +&
      &16_8*int(nx_pad,8)*int(ny_pad,8)
      omp_threads_cap = omp_get_max_threads()
      mem_safe_bytes_io = mem_safe_bytes - transient_bytes_per_plane
      if (mem_safe_bytes_io.le.0_8) then
         write(*,'(A,F0.2,A)') 'WARNING: convolve working memory alone (',&
         &real(transient_bytes_per_plane,8)/(1024.0d0**3),&
         &' GB/plane) exceeds the mem_frac_ram budget -- block shrunk to'//&
         &' 1 plane, but even that may not fit; raise mem_frac_ram or add'//&
         &' more RAM.'
         mem_safe_bytes_io = bytes_per_plane
      endif
      is_exact = do_reproject_l .and. interp_scheme_l.eq.interp_exact
      if (is_exact) then
         call exact_budget(infile, reffile_l, used_target_grid_l, target_cdelt1_arg,&
         &target_cdelt2_arg, nx_out, ny_out, mem_safe_bytes_io, status)
         if (status.ne.0) return
      endif
      block_planes64 = max(1_8, mem_safe_bytes_io / bytes_per_plane)
      block_planes64 = min(block_planes64, max_elements_per_block /&
      &max(1_8, max(int(nx_in,8)*int(ny_in,8), int(nx_out,8)*int(ny_out,8))))
      block_planes64 = max(1_8, block_planes64)
      block_planes = int(min(block_planes64, int(nfreq,8)))
      if (block_planes.lt.1) block_planes = 1

      write(*,'(A,A,A,I0,A,I0,A)') 'Writing ', trim(outfile), ': ', nfreq,&
      &' plane(s), in blocks of up to ', block_planes, ' plane(s)'

      allocate(block_in(nx_in, ny_in, block_planes))
      allocate(block_out(nx_out, ny_out, block_planes, 0:1))

      ! FFTW's plan is sized for one specific (nx,ny) and MUST NOT be
      ! executed against arrays of a different size (silent heap
      ! corruption, not a clean error -- caught the hard way via a
      ! munmap_chunk() crash before this fix). conv_nx/conv_ny (and the
      ! resulting nx_pad/ny_pad) were already determined above, ahead of
      ! the memory budget that needs them. omp_threads_cap (T28): this
      ! plan now executes each transform using this many threads
      ! INTERNALLY -- see plan_convolution's own comment.
      call plan_convolution(conv_nx, conv_ny, plan_fwd, plan_bwd, nx_pad, ny_pad, omp_threads_cap)
      if (nx_pad.ne.conv_nx .or. ny_pad.ne.conv_ny) then
         write(*,'(A,I0,A,I0,A,I0,A,I0,A)') 'Convolution FFT padded from ',&
         &conv_nx, 'x', conv_ny, ' to ', nx_pad, 'x', ny_pad,&
         &' (next 7-smooth size -- avoids a large-prime-factor slowdown).'
      endif

      cur_slot = 0
      write_pending = .false.
      write_failed = .false.
      status_par = 0
      ! nthreads here is for the REPROJECT pass specifically (still
      ! parallel-across-planes, T28 kept its original strategy -- see
      ! the per-plane loop below) -- capped at block_planes since more
      ! threads than planes-in-flight would just sit idle under the
      ! dynamic schedule, not a memory-safety concern (omp_threads_cap
      ! itself no longer is one, see above).
      nthreads = max(1, min(omp_threads_cap, block_planes))
      nwriters_eff = max(1, min(nwriters_l, omp_get_max_threads()))
      if (interp_scheme_l.eq.ast__ukern1) then
         interp_finterp_l => lanczos_kernel
      else
         interp_finterp_l => ast_null
      endif
      if (is_exact) then
         call timer_start(t_stage)
         call build_exact_weights_file(reffile_l, infile, used_target_grid_l,&
         &target_system_arg, target_projection_arg, target_center_lon_arg,&
         &target_center_lat_arg, target_cdelt1_arg, target_cdelt2_arg,&
         &target_final_crpix_arg, lbnd_out_d, ubnd_out_d, nx_in, ny_in, exact_wts, status)
         call timer_stop('exact_weights', t_stage)
         if (status.ne.0) return
      endif
      chan_start = 1
      iblock = 0
      do while (chan_start.le.nfreq)
         chan_len = min(block_planes, nfreq - chan_start + 1)
         iblock = iblock + 1

         call timer_start(t_stage)
         call read_freq_block(infile, naxis, sky1, sky2, freq_axis,&
         &naxes, chan_start, chan_len, nx_in, ny_in, block_in(:,:,1:chan_len), status_par)
         call timer_stop('block_read', t_stage)
         if (status_par.ne.0) exit

         call timer_start(t_stage)
         ! T28 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): convolve and
         ! reproject are now TWO SEPARATE passes over this block's own
         ! planes, not one combined per-plane loop -- convolve always
         ! runs SERIAL-over-planes/THREADED-within-each-plane (via
         ! plan_convolution's own nthreads), reproject always runs
         ! PARALLEL-over-planes exactly as before (T20's own AST
         ! setup fix, unchanged). Which one is pass 1 vs pass 2 follows
         ! convolve_first_l, same as the original single-pass ordering.
         ! Cross-checked before implementing (see T28's own ticket):
         ! this preserves reproject's existing cross-plane parallelism
         ! (risk found -- a naive single serial pass for BOTH
         ! steps would have made reproject-heavy files slower, since
         ! ast_resampler has no in-call threading of its own to give
         ! it the same treatment convolve_to_beam gets).
         if (.not. do_reproject_l) then
            ! stages=convolve alone: one serial/in-plane-threaded pass,
            ! straight into block_out -- no reproject step exists.
            do local_iplane = 1, chan_len
               if (status_par.ne.0) exit
               ich = chan_start + local_iplane - 1
               if (isbad(ich)) then
                  nanval = ieee_value(1.0_dp, ieee_quiet_nan)
                  block_out(:,:,local_iplane,cur_slot) = real(nanval)
                  cycle
               endif
               if (.not. allocated(plane_native)) allocate(plane_native(nx_in,ny_in))
               plane_native = real(block_in(:,:,local_iplane), dp)
               t_thread_start = omp_get_wtime()
               call timer_start(t_conv_local)
               call convolve_to_beam(plan_fwd, plan_bwd, plane_native, nx_in, ny_in, nx_pad, ny_pad,&
               &dx_deg, dy_deg, bmaj_in(ich)/3600.0_dp, bmin_in(ich)/3600.0_dp,&
               &bpa_in_pixel(ich), tgt_bmaj/3600.0_dp, tgt_bmin/3600.0_dp,&
               &tgt_bpa_pixel_native, plane_native, k)
               call timer_stop('convolve_compute', t_conv_local)
               t_thread_elapsed = (omp_get_wtime() - t_thread_start) * 1000.0_dp
               write(thread_msg,'(A,I0,A,I0,A,I0,A,F10.3)')&
               &'thread_timing stage=convolve event=done tid=0 block=',&
               &iblock,' plane=',ich,' nthreads=',omp_threads_cap,&
               &' dur_ms=',t_thread_elapsed
               call log_message('debug','tile_thread',trim(thread_msg))
               block_out(:,:,local_iplane,cur_slot) = real(plane_native)
               if (k.ne.0) status_par = -1
            enddo

         else if (convolve_first_l) then
            ! order=convolve_reproject (default): PASS 1 = convolve,
            ! serial/in-plane-threaded, overwriting block_in IN PLACE
            ! (native grid in, native grid convolved result out -- same
            ! shape, so the existing buffer is reused, no new
            ! allocation). Bad channels are left untouched here (still
            ! whatever read_freq_block put there) -- pass 2 below
            ! NaN-fills block_out for them directly and never reads
            ! this plane's own block_in content, so it does not matter
            ! that it was never convolved.
            do local_iplane = 1, chan_len
               if (status_par.ne.0) exit
               ich = chan_start + local_iplane - 1
               if (isbad(ich)) cycle
               if (.not. allocated(plane_native)) allocate(plane_native(nx_in,ny_in))
               plane_native = real(block_in(:,:,local_iplane), dp)
               t_thread_start = omp_get_wtime()
               call timer_start(t_conv_local)
               call convolve_to_beam(plan_fwd, plan_bwd, plane_native, nx_in, ny_in, nx_pad, ny_pad,&
               &dx_deg, dy_deg, bmaj_in(ich)/3600.0_dp, bmin_in(ich)/3600.0_dp,&
               &bpa_in_pixel(ich), tgt_bmaj/3600.0_dp, tgt_bmin/3600.0_dp,&
               &tgt_bpa_pixel_native, plane_native, k)
               call timer_stop('convolve_compute', t_conv_local)
               t_thread_elapsed = (omp_get_wtime() - t_thread_start) * 1000.0_dp
               write(thread_msg,'(A,I0,A,I0,A,I0,A,F10.3)')&
               &'thread_timing stage=convolve event=done tid=0 block=',&
               &iblock,' plane=',ich,' nthreads=',omp_threads_cap,&
               &' dur_ms=',t_thread_elapsed
               call log_message('debug','tile_thread',trim(thread_msg))
               block_in(:,:,local_iplane) = real(plane_native)
               if (k.ne.0) status_par = -1
            enddo

            ! PASS 2 = reproject, parallel-across-planes, EXACTLY the
            ! original strategy (T20's own AST-setup critical section
            ! unchanged) -- reads the now-convolved block_in.
            if (status_par.eq.0) then
               !$omp parallel num_threads(nthreads) default(none)&
               !$omp& shared(chan_len, nx_in, ny_in, nx_out, ny_out, cur_slot,&
               !$omp& block_in, block_out, isbad, chan_start, status_par,&
               !$omp& reffile_l, infile, lbnd_out_d, ubnd_out_d, iblock, nthreads,&
               !$omp& interp_scheme_l, interp_params_arr_l, interp_finterp_l,&
               !$omp& used_target_grid_l, target_system_arg, target_projection_arg,&
      !$omp& target_center_lon_arg,&
               !$omp& target_center_lat_arg, target_cdelt1_arg, target_cdelt2_arg,&
               !$omp& target_final_crpix_arg, is_exact, exact_wts)&
               !$omp& private(local_iplane, ich, nanval, plane_native_sp,&
               !$omp& plane_out_sp, t_status, t_wcs_ref, t_skymap_ref,&
               !$omp& t_skyframe_ref, t_naxes_ref, t_pixaxes_ref, t_wcs_in,&
               !$omp& t_skymap_in, t_skyframe_in, t_naxes_in2, t_pixaxes_in,&
               !$omp& t_map_in2ref, lbnd_in, ubnd_in, lbnd_o, ubnd_o,&
               !$omp& lbnd_reg, ubnd_reg, row0, row1, rows_per, rem,&
               !$omp& badval_sp, nbad, tid_local,&
               !$omp& t_thread_start, t_thread_elapsed, thread_msg,&
               !$omp& t_resamp_local)

               tid_local = omp_get_thread_num()

               ! Own private AST Mapping per thread -- required, see
               ! process_one_file_general's own comment on why (this
               ! Fortran AST binding has no lock/unlock for cross-thread
               ! sharing). T20: the creation step itself (not the
               ! resulting per-thread objects) is what deadlocks under
               ! concurrency.
               t_status = 0
               !$omp critical (ast_setup)
               call ast_begin(t_status)
               if (used_target_grid_l) then
                  call build_target_skymap(target_system_arg, target_projection_arg,&
         &target_center_lon_arg,&
                  &target_center_lat_arg, target_cdelt1_arg, target_cdelt2_arg,&
                  &target_final_crpix_arg, t_skymap_ref, t_skyframe_ref,&
                  &t_pixaxes_ref, t_status)
               else
                  call load_wcs(reffile_l, t_wcs_ref, t_naxes_ref, t_status)
                  call extract_sky_mapping(t_wcs_ref, t_skymap_ref, t_skyframe_ref,&
                  &t_pixaxes_ref, t_status)
               endif
               call load_wcs(infile, t_wcs_in, t_naxes_in2, t_status)
               call extract_sky_mapping(t_wcs_in, t_skymap_in, t_skyframe_in,&
               &t_pixaxes_in, t_status)
               call compose_pix2pix(t_skymap_in, t_skyframe_in, t_skymap_ref,&
               &t_skyframe_ref, t_map_in2ref, t_status)
               !$omp end critical (ast_setup)
               if (t_status.ne.0) then
                  !$omp atomic write
                  status_par = -1
               endif

               ! T30 item 3/7: spatial/tile parallelism -- see
               ! process_one_file_general's own identical comment. No
               ! worksharing over local_iplane any more -- every thread
               ! processes every plane, each working its own row-range.
               rows_per = ny_out / nthreads
               rem = mod(ny_out, nthreads)
               if (tid_local < rem) then
                  row0 = tid_local*(rows_per+1) + 1
                  row1 = row0 + rows_per
               else
                  row0 = rem*(rows_per+1) + (tid_local-rem)*rows_per + 1
                  row1 = row0 + rows_per - 1
               endif
               ! row1 < row0 here means nthreads > ny_out -- see the
               ! identical guard's own comment above process_one_file_
               ! general's row-split (this file). This thread has no
               ! rows assigned; the isbad NaN-fill below is already a
               ! harmless empty Fortran slice in that case, but the
               ! resample call is not (AST rejects an empty region), so
               ! skip this thread's work entirely rather than only the
               ! resample half.
               do local_iplane = 1, chan_len
                  ich = chan_start + local_iplane - 1
                  if (t_status.eq.0 .and. status_par.eq.0 .and. row1.ge.row0) then
                     if (isbad(ich)) then
                        nanval = ieee_value(1.0_dp, ieee_quiet_nan)
                        block_out(:,row0:row1,local_iplane,cur_slot) = real(nanval)
                     else
                        lbnd_in(1) = 1
                        lbnd_in(2) = 1
                        ubnd_in(1) = nx_in
                        ubnd_in(2) = ny_in
                        lbnd_o(1) = nint(lbnd_out_d(1))
                        lbnd_o(2) = nint(lbnd_out_d(2))
                        ubnd_o(1) = nint(ubnd_out_d(1))
                        ubnd_o(2) = nint(ubnd_out_d(2))
                        lbnd_reg(1) = lbnd_o(1)
                        lbnd_reg(2) = lbnd_o(2) + row0 - 1
                        ubnd_reg(1) = ubnd_o(1)
                        ubnd_reg(2) = lbnd_o(2) + row1 - 1
                        if (.not. allocated(plane_native_sp)) allocate(plane_native_sp(nx_in,ny_in))
                        if (.not. allocated(plane_out_sp)) allocate(plane_out_sp(nx_out,ny_out))
                        ! block_in already holds the CONVOLVED result
                        ! (pass 1, above) -- same real(4) kind as
                        ! plane_native_sp, no cast needed, unlike pass
                        ! 1's own real(dp) conversion.
                        plane_native_sp = block_in(:,:,local_iplane)
                        badval_sp = ieee_value(badval_sp, ieee_quiet_nan)
                        t_thread_start = omp_get_wtime()
                        call timer_start(t_resamp_local)
                        if (is_exact) then
                           call exact_apply_rows(exact_wts, plane_native_sp, plane_out_sp, row0, row1)
                        else
                        nbad = ast_resampler(t_map_in2ref, 2, lbnd_in, ubnd_in,&
                        &plane_native_sp, plane_native_sp,&
                        &interp_scheme_l, interp_finterp_l, interp_params_arr_l, 0, 0.0d0, 100, badval_sp,&
                        &2, lbnd_o, ubnd_o, lbnd_reg, ubnd_reg,&
                        &plane_out_sp, plane_out_sp, t_status)
                        endif
                        call timer_stop('reproject_compute', t_resamp_local)
                        if (t_status.ne.0) then
                           !$omp atomic write
                           status_par = -1
                        else
                           block_out(:,row0:row1,local_iplane,cur_slot) = plane_out_sp(:,row0:row1)
                        endif
                     endif
                  endif
                  !$omp barrier
                  ! T30 item 7 follow-up (2026-08-09): per-PLANE, not
                  ! per-block, team timing -- see process_one_file_
                  ! general's own identical comment. Skips isbad planes,
                  ! matching convolve's own per-plane logging convention
                  ! (no line for a plane it never worked).
                  if (t_status.eq.0 .and. status_par.eq.0 .and. .not.isbad(ich)) then
                     t_thread_elapsed = (omp_get_wtime() - t_thread_start) * 1000.0_dp
                     write(thread_msg,'(A,I0,A,I0,A,I0,A,I0,A,F10.3)') 'thread_timing stage=reproject event=done tid=',&
                     &tid_local,' block=',iblock,' plane=',ich,' nthreads=',nthreads,&
                     &' dur_ms=',t_thread_elapsed
                     call log_message('debug','tile_thread',trim(thread_msg))
                  endif
               enddo
               call ast_end(t_status)
               !$omp end parallel
            endif

         else
            ! order=reproject_convolve: PASS 1 = reproject,
            ! parallel-across-planes, EXACTLY the original strategy,
            ! writing into block_out (not yet convolved).
            !$omp parallel num_threads(nthreads) default(none)&
            !$omp& shared(chan_len, nx_in, ny_in, nx_out, ny_out, cur_slot,&
            !$omp& block_in, block_out, isbad, chan_start, status_par,&
            !$omp& reffile_l, infile, lbnd_out_d, ubnd_out_d, iblock, nthreads,&
            !$omp& interp_scheme_l, interp_params_arr_l, interp_finterp_l,&
            !$omp& used_target_grid_l, target_system_arg, target_projection_arg,&
      !$omp& target_center_lon_arg,&
            !$omp& target_center_lat_arg, target_cdelt1_arg, target_cdelt2_arg,&
            !$omp& target_final_crpix_arg, is_exact, exact_wts)&
            !$omp& private(local_iplane, ich, nanval, t_status,&
            !$omp& t_wcs_ref, t_skymap_ref, t_skyframe_ref, t_naxes_ref,&
            !$omp& t_pixaxes_ref, t_wcs_in, t_skymap_in, t_skyframe_in,&
            !$omp& t_naxes_in2, t_pixaxes_in, t_map_in2ref, lbnd_in, ubnd_in,&
            !$omp& lbnd_o, ubnd_o, lbnd_reg, ubnd_reg, row0, row1, rows_per, rem,&
            !$omp& badval_sp, nbad, tid_local,&
            !$omp& t_thread_start, t_thread_elapsed, thread_msg, t_resamp_local)

            tid_local = omp_get_thread_num()

            t_status = 0
            !$omp critical (ast_setup)
            call ast_begin(t_status)
            if (used_target_grid_l) then
               call build_target_skymap(target_system_arg, target_projection_arg,&
         &target_center_lon_arg,&
               &target_center_lat_arg, target_cdelt1_arg, target_cdelt2_arg,&
               &target_final_crpix_arg, t_skymap_ref, t_skyframe_ref,&
               &t_pixaxes_ref, t_status)
            else
               call load_wcs(reffile_l, t_wcs_ref, t_naxes_ref, t_status)
               call extract_sky_mapping(t_wcs_ref, t_skymap_ref, t_skyframe_ref,&
               &t_pixaxes_ref, t_status)
            endif
            call load_wcs(infile, t_wcs_in, t_naxes_in2, t_status)
            call extract_sky_mapping(t_wcs_in, t_skymap_in, t_skyframe_in,&
            &t_pixaxes_in, t_status)
            call compose_pix2pix(t_skymap_in, t_skyframe_in, t_skymap_ref,&
            &t_skyframe_ref, t_map_in2ref, t_status)
            !$omp end critical (ast_setup)
            if (t_status.ne.0) then
               !$omp atomic write
               status_par = -1
            endif

            ! T30 item 3/7: spatial/tile parallelism -- see
            ! process_one_file_general's own identical comment. No
            ! worksharing over local_iplane any more -- every thread
            ! processes every plane, each working its own row-range.
            rows_per = ny_out / nthreads
            rem = mod(ny_out, nthreads)
            if (tid_local < rem) then
               row0 = tid_local*(rows_per+1) + 1
               row1 = row0 + rows_per
            else
               row0 = rem*(rows_per+1) + (tid_local-rem)*rows_per + 1
               row1 = row0 + rows_per - 1
            endif
            ! row1 < row0 here means nthreads > ny_out -- see the
            ! identical guard's own comment above process_one_file_
            ! general's row-split (this file). Skip this thread's work
            ! entirely rather than feeding astResampleR an empty region.
            do local_iplane = 1, chan_len
               ich = chan_start + local_iplane - 1
               if (t_status.eq.0 .and. status_par.eq.0 .and. row1.ge.row0) then
                  if (isbad(ich)) then
                     nanval = ieee_value(1.0_dp, ieee_quiet_nan)
                     block_out(:,row0:row1,local_iplane,cur_slot) = real(nanval)
                  else
                     lbnd_in(1) = 1
                     lbnd_in(2) = 1
                     ubnd_in(1) = nx_in
                     ubnd_in(2) = ny_in
                     lbnd_o(1) = nint(lbnd_out_d(1))
                     lbnd_o(2) = nint(lbnd_out_d(2))
                     ubnd_o(1) = nint(ubnd_out_d(1))
                     ubnd_o(2) = nint(ubnd_out_d(2))
                     lbnd_reg(1) = lbnd_o(1)
                     lbnd_reg(2) = lbnd_o(2) + row0 - 1
                     ubnd_reg(1) = ubnd_o(1)
                     ubnd_reg(2) = lbnd_o(2) + row1 - 1
                     badval_sp = ieee_value(badval_sp, ieee_quiet_nan)
                     t_thread_start = omp_get_wtime()
                     call timer_start(t_resamp_local)
                     if (is_exact) then
                        call exact_apply_rows(exact_wts, block_in(:,:,local_iplane),&
                        &block_out(:,:,local_iplane,cur_slot), row0, row1)
                     else
                     nbad = ast_resampler(t_map_in2ref, 2, lbnd_in, ubnd_in,&
                     &block_in(:,:,local_iplane), block_in(:,:,local_iplane),&
                     &interp_scheme_l, interp_finterp_l, interp_params_arr_l, 0, 0.0d0, 100, badval_sp,&
                     &2, lbnd_o, ubnd_o, lbnd_reg, ubnd_reg,&
                     &block_out(:,:,local_iplane,cur_slot), block_out(:,:,local_iplane,cur_slot), t_status)
                     endif
                     call timer_stop('reproject_compute', t_resamp_local)
                     if (t_status.ne.0) then
                        !$omp atomic write
                        status_par = -1
                     endif
                  endif
               endif
               !$omp barrier
               ! T30 item 7 follow-up (2026-08-09): per-PLANE, not
               ! per-block, team timing -- see process_one_file_
               ! general's own identical comment. Skips isbad planes,
               ! matching convolve's own per-plane logging convention.
               if (t_status.eq.0 .and. status_par.eq.0 .and. .not.isbad(ich)) then
                  t_thread_elapsed = (omp_get_wtime() - t_thread_start) * 1000.0_dp
                  write(thread_msg,'(A,I0,A,I0,A,I0,A,I0,A,F10.3)') 'thread_timing stage=reproject event=done tid=',&
                  &tid_local,' block=',iblock,' plane=',ich,' nthreads=',nthreads,&
                  &' dur_ms=',t_thread_elapsed
                  call log_message('debug','tile_thread',trim(thread_msg))
               endif
            enddo
            call ast_end(t_status)
            !$omp end parallel

            ! PASS 2 = convolve, serial/in-plane-threaded, reading FROM
            ! and overwriting block_out IN PLACE (reprojected-not-yet-
            ! convolved in, convolved final result out -- same shape).
            if (status_par.eq.0) then
               do local_iplane = 1, chan_len
                  if (status_par.ne.0) exit
                  ich = chan_start + local_iplane - 1
                  if (isbad(ich)) cycle
                  if (.not. allocated(plane_out_arr)) allocate(plane_out_arr(nx_out,ny_out))
                  plane_out_arr = real(block_out(:,:,local_iplane,cur_slot), dp)
                  t_thread_start = omp_get_wtime()
                  call timer_start(t_conv_local)
                  call convolve_to_beam(plan_fwd, plan_bwd, plane_out_arr, nx_out, ny_out, nx_pad, ny_pad,&
                  &ref_dx_deg, ref_dy_deg, bmaj_in(ich)/3600.0_dp, bmin_in(ich)/3600.0_dp,&
                  &bpa_in_pixel(ich), tgt_bmaj/3600.0_dp, tgt_bmin/3600.0_dp,&
                  &tgt_bpa_pixel_out, plane_out_arr, k)
                  call timer_stop('convolve_compute', t_conv_local)
                  t_thread_elapsed = (omp_get_wtime() - t_thread_start) * 1000.0_dp
                  write(thread_msg,'(A,I0,A,I0,A,I0,A,F10.3)')&
                  &'thread_timing stage=convolve event=done tid=0 block=',&
                  &iblock,' plane=',ich,' nthreads=',omp_threads_cap,&
                  &' dur_ms=',t_thread_elapsed
                  call log_message('debug','tile_thread',trim(thread_msg))
                  block_out(:,:,local_iplane,cur_slot) = real(plane_out_arr)
                  if (k.ne.0) status_par = -1
               enddo
            endif
         endif
         call timer_stop('block_convolve', t_stage)
         if (status_par.ne.0) exit

         call timer_start(t_stage)
         if (write_pending) then
            call block_write_join(write_thread_id)
            write_pending = .false.
            if (write_failed) then
               status_par = -1
               exit
            endif
         endif
         call timer_stop('block_write_join', t_stage)
         write_job%file_path = outfile
         write_job%datastart = datastart
         write_job%chan_start = chan_start
         write_job%chan_len = chan_len
         write_job%nwriters_eff = max(1, min(nwriters_eff, chan_len))
         write_job%nx = nx_out
         write_job%ny = ny_out
         write_job%data => block_out(:,:,1:chan_len,cur_slot)
         call timer_start(t_stage)
         if (io_overlap_l) then
            call block_write_dispatch_async(write_job, write_thread_id,&
            &write_dispatched_ok)
            write_pending = write_dispatched_ok
         else
            call do_block_write(write_job)
            if (write_failed) then
               status_par = -1
               exit
            endif
         endif
         call timer_stop('block_write', t_stage)

         chan_start = chan_start + chan_len
         cur_slot = 1 - cur_slot
      enddo

      call timer_start(t_stage)
      if (write_pending) then
         call block_write_join(write_thread_id)
         write_pending = .false.
         if (write_failed) status_par = -1
      endif
      call timer_stop('block_write_join', t_stage)

      call destroy_convolution_plan(plan_fwd, plan_bwd)
      deallocate(block_in, block_out)
      if (is_exact) call exact_free(exact_wts)

      if (status_par.ne.0) then
         write(*,*) 'ERROR: failed to process/write one or more planes for: ', trim(infile)
         status = -1
         ! out_unit is NOT open here -- closed right after FTGHAD, above,
         ! before the block-write loop ever started.
         return
      endif

      ! CASAMBM/BEAMS: always attached whenever convolve is active
      ! (stages=convolve or both, alone or chained) -- see
      ! write_beams_table_match's own comment for why, and
      ! convolve_cubes.f90's identical write_beams_table this is a
      ! verbatim port of. Reopen out_unit now that every block write/
      ! join above has finished (block_write_join, right
      ! before this point, guarantees no writer thread is still active)
      ! -- CFITSIO writes this extra HDU, which raw stream writes have
      ! no way to add.
      fitsstat = 0
      call safe_ftopen(out_unit, trim(outfile), 1, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to reopen output for BEAMS table: ', trim(outfile)
         call printerror(fitsstat)
         status = -1
         return
      endif
      call FTPKYL(out_unit, 'CASAMBM', .true.,&
      &'Multiple beams per plane (see BEAMS ext)', fitsstat)
      call write_beams_table_match(out_unit, nfreq, isbad, tgt_bmaj,&
      &tgt_bmin, tgt_bpa, status_par)
      if (status_par.ne.0) then
         status = -1
         call safe_ftclos(out_unit, fitsstat)
         return
      endif

      call safe_ftclos(out_unit, fitsstat)
      call log_message('info', 'convolve', 'finished: '//trim(infile))
      call timer_report_file_summary(infile)
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
      call safe_ftopen(unit, trim(reffile_l), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open reference file: ', trim(reffile_l)
         status = -1
         call free_fits_unit(unit)
         return
      endif
      write(axstr,'(I0)') pixaxes_ref_l(1)
      fitsstat = 0
      call FTGKYD(unit, 'CDELT'//trim(axstr), ref_cdelt1, comment, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: missing CDELT for the reference''s RA axis in: ',&
         &trim(reffile_l)
         status = -1
         call safe_ftclos(unit, fitsstat)
         return
      endif
      write(axstr,'(I0)') pixaxes_ref_l(2)
      fitsstat = 0
      call FTGKYD(unit, 'CDELT'//trim(axstr), ref_cdelt2, comment, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: missing CDELT for the reference''s DEC axis in: ',&
         &trim(reffile_l)
         status = -1
         call safe_ftclos(unit, fitsstat)
         return
      endif
      call safe_ftclos(unit, fitsstat)
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
      group = 1
      badval = ieee_value(badval, ieee_quiet_nan)
      call safe_ftopen(unit, trim(filename), 0, blocksize, fitsstat)
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
      call safe_ftclos(unit, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to read block (planes ', chan_start, '-',&
         &chan_start+chan_len-1, ') from ', trim(filename)
         status = -1
      endif
   end subroutine read_freq_block

   logical function host_is_big_endian_mc() result(is_be)
      !! Verbatim port of rm_synthesis_mod.f90's own host_is_big_endian
      !! -- see there for the full rationale (FITS mandates big-endian;
      !! every realistic host here is little-endian, so this always
      !! ends up swapping in practice, but is checked at runtime rather
      !! than hard-coded so a big-endian host isn't silently
      !! double-swapped and corrupted).
      integer(int32) :: probe
      integer(int8) :: bytes(4)

      probe = 1_int32
      bytes = transfer(probe, bytes)
      is_be = (bytes(1) == 0_int8)
   end function host_is_big_endian_mc

   subroutine swap_bytes_r4_inplace_mc(buf, n)
      !! Verbatim port of rm_synthesis_mod.f90's own
      !! swap_bytes_r4_inplace -- reverses the 4 bytes of every real
      !! element of buf, in place. Self-inverse.
      integer(kind=int64), intent(in) :: n
      real, intent(inout) :: buf(n)
      integer(kind=int64) :: i
      integer(int8) :: b(4), t

      do i = 1_int64, n
         b = transfer(buf(i), b)
         t = b(1); b(1) = b(4); b(4) = t
         t = b(2); b(2) = b(3); b(3) = t
         buf(i) = transfer(b, buf(i))
      end do
   end subroutine swap_bytes_r4_inplace_mc

   subroutine write_freq_block_raw(file_path, datastart, nx, ny, chan_start,&
   &chan_len, block_data, status)
      !! Raw Fortran stream I/O at a computed byte offset -- CFITSIO is
      !! never touched here at all (see process_one_file_restricted's
      !! own comment for why: a real, timing-dependent hang, found via
      !! this exact code path on real WALLABY+EMU data, root-caused to
      !! two OS threads both inside CFITSIO concurrently -- the main
      !! thread's own next-block read and this write, on different
      !! handles, which CFITSIO does not guarantee is safe without a
      !! specific reentrant build). Verbatim design port of rm_synthesis
      !! _mod.f90's own write_rm_chunk_raw, simplified for this tool's
      !! own "restricted" scope: every axis other than sky1/sky2/
      !! freq_axis is degenerate (size 1) by construction (read_axis_
      !! info's own requirement), and there is no spatial sub-tiling at
      !! all in this tool (always the full nx*ny extent) -- so every
      !! channel plane is exactly nx*ny contiguous elements, planes
      !! back-to-back in channel order, i.e. always rm_synthesis_mod.
      !! f90's own "full-width" case, never its "partial-width" one.
      !!
      !! Only ever one writer in flight at a time by this whole file's
      !! own design (block_write_join always runs before the next
      !! dispatch) -- unlike rm_synthesis' own io_write_threads>1 (N
      !! concurrent writers), so no critical section is
      !! needed around the newunit= open here: there is no other thread
      !! that could be calling it at the same moment.
      character(len=*), intent(in) :: file_path
      integer(kind=8), intent(in) :: datastart
      integer, intent(in) :: nx, ny, chan_start, chan_len
      real, intent(in) :: block_data(:,:,:)
      integer, intent(inout) :: status

      integer :: u, ios, ip
      logical :: need_swap
      integer(kind=8) :: plane_stride_bytes, plane_elems, byte_pos
      real, allocatable :: plane_buf(:,:)

      if (status.ne.0) return

      need_swap = .not. host_is_big_endian_mc()
      plane_elems = int(nx,8) * int(ny,8)
      plane_stride_bytes = plane_elems * 4_8

      open(newunit=u, file=trim(file_path), access='stream',&
      &form='unformatted', status='old', action='write', iostat=ios)
      if (ios.ne.0) then
         write(*,*) 'ERROR: write_freq_block_raw: failed to open ', trim(file_path)
         status = -1
         return
      endif

      allocate(plane_buf(nx,ny))
      do ip = 1, chan_len
         plane_buf = block_data(:,:,ip)
         if (need_swap) call swap_bytes_r4_inplace_mc(plane_buf, plane_elems)
         byte_pos = datastart + int(chan_start-1+ip-1,8)*plane_stride_bytes + 1_8
         write(u, pos=byte_pos, iostat=ios) plane_buf
         if (ios.ne.0) then
            write(*,*) 'ERROR: write_freq_block_raw: write failed (plane ',&
            &chan_start+ip-1, ') for ', trim(file_path)
            status = -1
            exit
         endif
      enddo
      deallocate(plane_buf)
      close(u)
   end subroutine write_freq_block_raw

   subroutine do_block_write(job)
      !! T19 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): raw stream write,
      !! not CFITSIO -- see write_freq_block_raw's own comment.
      type(block_write_job_t), intent(inout) :: job
      integer :: status_local, wk, wbase, wrem, wlen_k, woff_k, wchan_start_k
      logical :: any_failed

      any_failed = .false.
      if (job%nwriters_eff.le.1) then
         status_local = 0
         call write_freq_block_raw(trim(job%file_path), job%datastart,&
         &job%nx, job%ny, job%chan_start, job%chan_len, job%data, status_local)
         if (status_local.ne.0) any_failed = .true.
      else
         ! nwriters_eff>1 (T21, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md):
         ! same even-split (base+remainder) scheme as convolve_cubes.f90's
         ! own do_block_write / rm_synthesis_mod.f90's own do_tile_write.
         wbase = job%chan_len / job%nwriters_eff
         wrem = mod(job%chan_len, job%nwriters_eff)
         !$omp parallel do num_threads(job%nwriters_eff) default(none)&
         !$omp& shared(job, wbase, wrem, any_failed)&
         !$omp& private(wk, wlen_k, woff_k, wchan_start_k, status_local)
         do wk = 0, job%nwriters_eff-1
            wlen_k = wbase + merge(1, 0, wk.lt.wrem)
            woff_k = wk*wbase + min(wk, wrem)
            wchan_start_k = job%chan_start + woff_k
            status_local = 0
            if (wlen_k.gt.0) then
               call write_freq_block_raw(trim(job%file_path), job%datastart,&
               &job%nx, job%ny, wchan_start_k, wlen_k,&
               &job%data(:,:,woff_k+1:woff_k+wlen_k), status_local)
               if (status_local.ne.0) then
                  !$omp atomic write
                  any_failed = .true.
               endif
            endif
         enddo
         !$omp end parallel do
      endif
      if (any_failed) then
         write(*,*) 'ERROR: background write failed for channels ',&
         &job%chan_start, '-', job%chan_start+job%chan_len-1
         write_failed = .true.
      endif
   end subroutine do_block_write

   function block_write_thread_entry(arg) bind(C) result(res)
      type(c_ptr), value :: arg
      type(c_ptr) :: res
      type(block_write_job_t), pointer :: job

      call c_f_pointer(arg, job)
      call do_block_write(job)
      res = c_null_ptr
   end function block_write_thread_entry

   subroutine block_write_dispatch_async(job, thread_id, dispatched)
      type(block_write_job_t), intent(inout), target :: job
      integer(c_long), intent(out) :: thread_id
      logical, intent(out) :: dispatched
      integer(c_int) :: rc

      rc = c_pthread_create(thread_id, c_null_ptr,&
      &c_funloc(block_write_thread_entry), c_loc(job))
      if (rc.ne.0_c_int) then
         write(*,*) 'WARNING: pthread_create failed for async block write;'//&
         &' running inline'
         call do_block_write(job)
         dispatched = .false.
      else
         dispatched = .true.
      endif
   end subroutine block_write_dispatch_async

   subroutine block_write_join(thread_id)
      integer(c_long), intent(in) :: thread_id
      integer(c_int) :: rc
      rc = c_pthread_join(thread_id, c_null_ptr)
   end subroutine block_write_join

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

   subroutine advise_beam_consistency(infiles, n_inputs)
      !! UR-17 (docs/dev/BRANCH_STATUS_update-reproject.md): own copy
      !! of reproject_cubes.f90's subroutine of the same name, per this
      !! project's adapt-don't-share convention for tool-local helpers
      !! (target_grid_mod.f90 is scoped to the output GRID, so a beam
      !! check does not belong there). Called for stages=reproject
      !! only -- see this program's own call site. Advisory
      !! only -- changes nothing, never fails the run. This tool
      !! resamples onto a common pixel grid and does not touch angular
      !! resolution (convolve_cubes' own job), so infiles with different
      !! beams keep them.
      !!
      !! Warned about here because rm_synthesis does not refuse such
      !! data: its own cross-band beam check (src/rm_synthesis.f90, the
      !! n_bands_t2>1 block, read directly) warns, keeps processing, and
      !! propagates only the reference band's BMAJ/BMIN/BPA to its
      !! outputs, treating a beam mismatch as a metadata concern rather
      !! than a numerical-correctness one. A mismatch is therefore
      !! survivable downstream and easy to miss.
      !!
      !! Compares each infile's nominal header BMAJ/BMIN/BPA against the
      !! first infile's (relative tolerance, not equality -- these are
      !! written from 32-bit values in practice), flags a presence
      !! mismatch (one has them, another does not), and flags CASAMBM=T,
      !! which means that file carries a per-channel BEAMS table and so
      !! has no single nominal beam. Silent for a single infile: there
      !! is no cross-file comparison to make.
      character(len=*), intent(in) :: infiles(:)
      integer, intent(in) :: n_inputs

      integer :: i, unit_b, blocksize_b, fitsstat_b
      character(len=68) :: comment
      double precision :: bmaj_i, bmin_i, bpa_i
      double precision :: bmaj_1, bmin_1, bpa_1
      logical :: have_i, have_1, casambm_i, any_casambm, mismatch
      logical :: casambm_present

      if (n_inputs.lt.2) return

      have_1 = .false.
      bmaj_1 = 0.0d0; bmin_1 = 0.0d0; bpa_1 = 0.0d0
      any_casambm = .false.
      mismatch = .false.

      do i = 1, n_inputs
         fitsstat_b = 0; blocksize_b = 1
         call safe_ftopen(unit_b, trim(infiles(i)), 0, blocksize_b, fitsstat_b)
         if (fitsstat_b.ne.0) then
            ! Unreadable here is not this check's business to report:
            ! every other stage opens the same file and will say so.
            call free_fits_unit(unit_b)
            cycle
         endif

         have_i = .true.
         fitsstat_b = 0
         call FTGKYD(unit_b, 'BMAJ', bmaj_i, comment, fitsstat_b)
         if (fitsstat_b.ne.0) have_i = .false.
         fitsstat_b = 0
         call FTGKYD(unit_b, 'BMIN', bmin_i, comment, fitsstat_b)
         if (fitsstat_b.ne.0) have_i = .false.
         fitsstat_b = 0
         call FTGKYD(unit_b, 'BPA', bpa_i, comment, fitsstat_b)
         if (fitsstat_b.ne.0) have_i = .false.

         fitsstat_b = 0
         casambm_i = .false.
         call FTGKYL(unit_b, 'CASAMBM', casambm_i, comment, fitsstat_b)
         casambm_present = (fitsstat_b.eq.0)
         fitsstat_b = 0
         call safe_ftclos(unit_b, fitsstat_b)

         if (casambm_present .and. casambm_i) then
            any_casambm = .true.
            write(*,'(A,A,A)') '  ', trim(infiles(i)), ': CASAMBM=T'//&
            &' (per-channel BEAMS table -- no single nominal beam)'
         endif

         if (i.eq.1) then
            have_1 = have_i
            if (have_i) then
               bmaj_1 = bmaj_i; bmin_1 = bmin_i; bpa_1 = bpa_i
            endif
            cycle
         endif

         if (have_i .neqv. have_1) then
            write(*,'(A,A,A)') '  ', trim(infiles(i)), ': BMAJ/BMIN/BPA'//&
            &' presence differs from the first infile'
            mismatch = .true.
         else if (have_i .and. have_1) then
            if (abs(bmaj_i-bmaj_1).gt.1.0d-6*max(abs(bmaj_1),abs(bmaj_i),1.0d-30) .or.&
            &abs(bmin_i-bmin_1).gt.1.0d-6*max(abs(bmin_1),abs(bmin_i),1.0d-30) .or.&
            &abs(bpa_i-bpa_1).gt.1.0d-6*max(abs(bpa_1),abs(bpa_i),1.0d-30)) then
               write(*,'(A,A,A,ES12.5,A,ES12.5,A,F0.3)') '  ', trim(infiles(i)),&
               &': beam BMAJ=', bmaj_i, ' BMIN=', bmin_i, ' BPA=', bpa_i
               mismatch = .true.
            endif
         endif
      enddo

      if (mismatch .or. any_casambm) then
         write(*,'(A)') 'WARNING: the infiles do not share one angular resolution'//&
         &' (see the lines above).'
         write(*,'(A)') '  This tool matches the pixel grid only -- it does not change'//&
         &' resolution, so the outputs will still differ.'
         write(*,'(A)') '  rm_synthesis will combine them, warning but not refusing,'//&
         &' and will record only one band''s beam.'
         write(*,'(A)') '  Run convolve_cubes (or match_cubes stages=both) across all'//&
         &' of them for one common beam.'
      endif
   end subroutine advise_beam_consistency

   subroutine run_dry_run_check(target_path, tool_name)
      !! dry_run=y (T24, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): planning-
      !! only pass, mirroring rm_synthesis' own dry_run/tile_autotune.cfg
      !! convention -- touches no real data, just checks the target output
      !! disk's own rotational status (via /sys/block/<dev>/queue/
      !! rotational, resolved from target_path's own mount point -- the
      !! same, direct mechanism this session used to confirm /data1 is a
      !! spinning disk and /home is NVMe, no benchmarking needed)
      !! and writes a suggested <tool_name>_dryrun.cfg with recommended
      !! io_overlap/nwriters settings. Advisory only: this project's own
      !! nwriters clamp formula was deliberately left tied to
      !! omp_get_max_threads(), not disk type (see T7/T12's own
      !! discussion) -- this dry-run pass suggests a starting VALUE within
      !! that clamp, it does not change the clamp itself.
      !!
      !! Shells out (execute_command_line) rather than reimplementing
      !! mount-point/block-device resolution in Fortran -- df/basename/
      !! sed are already exactly what this session used interactively to
      !! confirm /data1 vs /home, so this reuses that same, already-
      !! verified logic instead of a parallel, untested reimplementation.
      !! Captures the shell command's own stdout via a plain output
      !! redirect to a file in the CURRENT working directory (never /tmp
      !! or any other system directory, matching this project's own
      !! standing convention for its own scratch files, applied here to
      !! the tool's own runtime behaviour too since HPC compute nodes
      !! often restrict or omit /tmp) -- read back, then deleted
      !! immediately.
      character(len=*), intent(in) :: target_path, tool_name

      character(len=2048) :: shell_cmd
      character(len=600) :: capture_file
      character(len=600) :: target_dir
      character(len=16) :: rot_str
      integer :: capture_unit, ios_cap, exitstat_cmd, cfg_unit
      integer :: rot_val, slash_pos
      logical :: is_rotational, have_answer

      ! Capture file and the suggested <tool_name>_dryrun.cfg both land
      ! next to target_path itself (guaranteed writable -- that's where
      ! the real output would go), not a bare CWD-relative name -- the
      ! latter broke whenever the running process doesn't own the
      ! current working directory (e.g. a non-root container's
      ! WORKDIR).
      slash_pos = index(trim(target_path), '/', back=.true.)
      if (slash_pos.gt.0) then
         target_dir = target_path(1:slash_pos)
      else
         target_dir = ''
      endif
      write(capture_file, '(A,A,I0,A)') trim(target_dir), '.dryrun_rotational_', getpid_wrap(), '.tmp'

      write(shell_cmd, '(A,A,A,A,A)')&
      &'dev=$(df --output=source ''', trim(target_path), ''' 2>/dev/null | tail -1); ',&
      &'base=$(basename "$dev"); ',&
      &'if echo "$base" | grep -Eq ''^nvme[0-9]+n[0-9]+p[0-9]+$''; then parent="${base%p*}"; '//&
      &'elif echo "$base" | grep -Eq ''^nvme[0-9]+n[0-9]+$''; then parent="$base"; '//&
      &'else parent=$(echo "$base" | sed -E ''s/[0-9]+$//''); fi; '//&
      &'cat "/sys/block/$parent/queue/rotational" 2>/dev/null > '//trim(capture_file)

      call execute_command_line(trim(shell_cmd), wait=.true., exitstat=exitstat_cmd)

      have_answer = .false.
      rot_val = -1
      open(newunit=capture_unit, file=trim(capture_file), status='old',&
      &action='read', iostat=ios_cap)
      if (ios_cap.eq.0) then
         read(capture_unit, '(A)', iostat=ios_cap) rot_str
         close(capture_unit, status='delete')
         if (ios_cap.eq.0) then
            read(rot_str, *, iostat=ios_cap) rot_val
            if (ios_cap.eq.0 .and. (rot_val.eq.0 .or. rot_val.eq.1)) have_answer = .true.
         endif
      else
         close(capture_unit, status='delete', iostat=ios_cap)
      endif

      is_rotational = (rot_val.eq.1)

      write(*,'(A)') 'dry_run=y: no files will be processed.'
      write(*,'(A,A)') 'Target path checked: ', trim(target_path)

      open(newunit=cfg_unit, file=trim(target_dir)//trim(tool_name)//'_dryrun.cfg',&
      &status='replace', action='write', iostat=ios_cap)
      if (ios_cap.ne.0) then
         write(*,'(A)') 'WARNING: could not write suggested cfg file.'
         return
      endif
      write(cfg_unit, '(A)') '# Autogenerated by '//trim(tool_name)//' dry_run=y'
      write(cfg_unit, '(A)') '# Target path checked: '//trim(target_path)

      if (.not. have_answer) then
         ! Common inside a container: df reports the container's own
         ! overlay filesystem, which has no backing block device, so
         ! /sys/block lookup fails -- not a bug, an isolation boundary (see
         ! src/thread_safety_mod.f90's own header comment for the
         ! analogous CPU-topology case). Rather than silently suggesting
         ! nothing, give both suggestions and let the user pick the one
         ! that matches their own target storage, which they know and
         ! this process cannot see.
         write(*,'(A)') 'Could not determine the target disk''s rotational status'//&
         &' (df/sysfs check failed or path not found).'
         write(*,'(A)') 'If your target disk is spinning (HDD): io_overlap=n,'//&
         &' nwriters=1 -- concurrent reads/writes contend for the same physical head.'
         write(*,'(A)') 'If your target disk is SSD/NVMe (non-rotational):'//&
         &' io_overlap=y, nwriters=2 -- concurrent I/O is generally safe; raise'//&
         &' nwriters further if your own machine has more idle cores available.'
         write(*,'(A)') '  see docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md T24 for the full'//&
         &' reasoning (this is advisory, not enforced).'
         write(cfg_unit, '(A)') '# Could not automatically determine whether the target'
         write(cfg_unit, '(A)') '# disk is spinning (HDD) or solid-state (SSD/NVMe) --'
         write(cfg_unit, '(A)') '# uncomment ONE block below based on your own knowledge'
         write(cfg_unit, '(A)') '# of the target storage, and delete the other.'
         write(cfg_unit, '(A)') '#'
         write(cfg_unit, '(A)') '# If the target disk is spinning (HDD):'
         write(cfg_unit, '(A)') '# io_overlap=n'
         write(cfg_unit, '(A)') '# nwriters=1'
         write(cfg_unit, '(A)') '#'
         write(cfg_unit, '(A)') '# If the target disk is SSD/NVMe (non-rotational):'
         write(cfg_unit, '(A)') '# io_overlap=y'
         write(cfg_unit, '(A)') '# nwriters=2'
         close(cfg_unit)
         write(*,'(A)') 'Wrote both suggestions to '//trim(target_dir)//trim(tool_name)//&
         &'_dryrun.cfg -- uncomment the one that matches your target storage.'
         return
      endif

      if (is_rotational) then
         write(*,'(A)') 'Target disk: spinning (rotational=1).'
         write(*,'(A)') 'Suggestion: io_overlap=n, nwriters=1 -- on a single spinning'//&
         &' disk, concurrent reads/writes contend for the same physical head;'
         write(*,'(A)') '  see docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md T24 for the full'//&
         &' reasoning (this is advisory, not enforced).'
         write(cfg_unit, '(A)') '# Target disk: spinning (rotational=1)'
         write(cfg_unit, '(A)') 'io_overlap=n'
         write(cfg_unit, '(A)') 'nwriters=1'
      else
         write(*,'(A)') 'Target disk: non-rotational (SSD/NVMe, rotational=0).'
         write(*,'(A)') 'Suggestion: io_overlap=y, nwriters=2 -- concurrent I/O is'//&
         &' generally safe on SSD/NVMe; nwriters starts conservative (this machine''s'
         write(*,'(A)') '  own spare-core headroom outside OMP_NUM_THREADS), raise it'//&
         &' if your own machine has more idle cores available.'
         write(cfg_unit, '(A)') '# Target disk: non-rotational (SSD/NVMe, rotational=0)'
         write(cfg_unit, '(A)') 'io_overlap=y'
         write(cfg_unit, '(A)') 'nwriters=2'
      endif
      write(cfg_unit, '(A)') '# Copy these KEY=VALUE lines into your own cfg, or pass'
      write(cfg_unit, '(A)') '# them directly on the command line.'
      close(cfg_unit)
      write(*,'(A)') 'Wrote suggested settings to '//trim(target_dir)//trim(tool_name)//'_dryrun.cfg'
   end subroutine run_dry_run_check

   integer function getpid_wrap() result(pid)
      !! Wraps the getpid() C library call (via iso_c_binding) to
      !! make the dry_run capture-file name collision-proof if multiple
      !! dry_run invocations somehow overlap in the same directory --
      !! not needed for correctness of a single run, just cheap safety.
      use, intrinsic :: iso_c_binding, only: c_int
      interface
         function c_getpid() bind(C, name="getpid") result(r)
            import :: c_int
            integer(c_int) :: r
         end function c_getpid
      end interface
      pid = int(c_getpid())
   end function getpid_wrap

end program match_cubes
