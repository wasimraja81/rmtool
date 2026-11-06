! reproject_cubes -- standalone pre-rm-synthesis tool (geometry-matching
! project, planning/MULTI_BAND_TOMOGRAPHY_PLAN.md): reprojects a set of
! FITS cubes -- two or more, not tied to the multi-band-tomography "band"
! concept specifically -- onto a single reference cube's grid using
! Starlink AST for WCS handling and resampling, so the existing
! rm_synthesis exact-match ingestion (T1) can consume misaligned
! bands unchanged.
!
! Current stage: N-input footprint-mode computation (intersection/union/
! reference), built on top of two previously-verified building blocks:
! extract_sky_mapping (automatic sky-axis detection, correct for any axis
! order/adjacency -- see its own comment for why the earlier consecutive-
! axis-only search was wrong) and astMapBox-based footprint bounds
! (verified against a partial-overlap case, not just full
! overlap). Cross-file alignment composes each file's own pixel->sky
! Mapping with an astConvert between the two SkyFrame *objects* (not
! whole FrameSets -- astConvert's domain search does not recurse into a
! CmpFrame's internal components, confirmed against the SUN/211
! manual, so it cannot align two compound "STOKES-SKY-SPECTRUM" current
! frames directly).
!
! Footprint policy, revised 2026-10-14 (UR-15, docs/dev/BRANCH_STATUS_
! update-reproject.md): zero overlap between any infile and the
! output grid is a hard failure in every mode, including mode=
! reference -- which earlier had no check at all, and before that had
! a header comment claiming a universal refusal the code never ran
! for it. A band sharing no sky at all with the rest is
! never legitimate partial coverage (that implies at least some
! shared sky); it is almost always the wrong file, so letting it
! through as a silent all-NaN plane was reversed in favour of
! refusing it the same way union/intersection already do. The check
! (target_grid_mod.f90's own check_overlap_and_size_grid, shared by
! every mode and by the no-reffile path too) scans every infile
! before reporting, rather than stopping at the first offender, so a
! systematic problem -- e.g. a whole band's own Q/U pair -- can be
! told apart from a single mismatched file. Partial overlap (the
! normal, expected case for combining bands with different sky
! coverage) is never rejected, in any mode.
!
! Also demonstrates the regridding step: reads planes and
! resamples them onto the final output grid via astResampleR (originally
! one plane at a time, via a since-removed resample_one_plane -- see the
! OpenMP/block paragraph below for its replacement). Caught a real bug on
! the first attempt: the documented
! 20-argument signature (this, ndim_in, lbnd_in, ubnd_in, in, in_var,
! interp, finterp, params, flags, tol, maxpix, badval, ndim_out,
! lbnd_out, ubnd_out, lbnd, ubnd, out, out_var) was missing "params",
! silently shifting every later argument by one position (a
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
! Writes an output FITS file now (write_reprojected_file):
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
! OpenMP-parallelised across planes, and reads/resamples/writes in
! memory-budgeted BLOCKS of planes (mem_frac_ram, same concept as
! rm_synthesis's own tile planner -- see get_mem_total_kb) rather than
! one plane at a time. Each block goes through three strictly separated
! phases: read (one thread, OMP `single`), resample (every thread, in
! parallel -- each thread builds its own private AST Mapping, since this
! Fortran AST binding has no astLock_/astUnlock_ to hand one Mapping
! between threads), write (one thread, `single`). Blocks mean a whole
! block's CFITSIO I/O happens on a single thread BY CONSTRUCTION, so
! there is nothing to lock -- an earlier plane-at-a-time version
! serialised every CFITSIO call behind an OMP critical section instead,
! which this replaced (see write_reprojected_file's own comment). A
! bug surfaced getting here: FTGSVE fills its output array in ascending-
! axis-number order among the axes being read, and reading a
! whole block (not a single degenerate plane) exposed a case the old
! per-plane code never hit -- TEST_NONADJACENT.Q.FITSCUBE has FREQ on
! axis 1, *before* RA/DEC on axes 2 and 4, so a block read's fastest
! dimension is the block axis, not the sky axes; read_one_block now
! ranks the 3 relevant axes and only pays for a permute copy on
! non-conventional orderings like that one, not the common case.
! Verified byte-identical (header and data) against the pre-blocking,
! pre-OpenMP committed output across all fixtures, multiple block sizes
! (including a 3-plane block against 200 channels, exercising a
! non-exact remainder), and 25 repeated stress runs at default thread
! count with no failures.
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
   use, intrinsic :: iso_c_binding, only: c_int, c_long, c_ptr, c_funptr,&
   &c_null_ptr, c_funloc, c_loc, c_f_pointer
   use, intrinsic :: iso_fortran_env, only: int8, int32, int64
   use logging_mod
   use fitsio_unit_mod
   use wcs_match_mod
   use thread_safety_mod, only: check_and_clamp_omp_threads
   use ast_interp_mod, only: resolve_interp, ast__ukern1, interp_exact
   use exact_overlap_mod, only: exact_weights_t, exact_init, exact_build_rows,&
   &exact_apply, exact_free, exact_weights_bytes, exact_estimate_bytes, mib_str
   use target_grid_mod, only: resolve_target_grid, write_target_wcs_cards,&
   &build_target_skymap, check_output_grid_validity, check_overlap_and_size_grid,&
   &check_disk_space, output_cell_size_arcsec, check_beam_sampling,&
   &fits_reference_cell_arcsec
   implicit none

   ! --- Logging & timing (planning-doc ticket) -- see convolve_cubes.f90
   ! /logging_mod.f90's own header comments. Stage names used below:
   ! startup, block_read, block_resample, block_write, block_write_join,
   ! finalize.
   character(len=16) :: log_level
   logical :: timing_enabled
   character(len=272) :: log_output_file
   ! AST_PAR (the vendor Fortran constants file, /usr/include/AST_PAR) is
   ! fixed-form Fortran 77 (`*`-column comments) and cannot be `include`d
   ! into a free-form .f90 file directly (gfortran misparses its comments
   ! as code). Only the handful of symbols used are declared
   ! directly instead, matching AST_PAR's own declared types exactly
   ! (checked against /usr/include/AST_PAR).
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
   ! Hard ceiling on elements-per-block -- see convolve_cubes.f90's own
   ! identical constant/comment (CFITSIO's Fortran wrapper takes a
   ! default-INTEGER element count; this clamp keeps
   ! max(nx_in*ny_in,nx_out*ny_out)*block_planes safely under 2^31-1 for
   ! ANY image size/mem_frac_ram combination).
   integer(kind=8), parameter :: max_elements_per_block = 2000000000_8

   ! --- io_overlap: background-thread block write (planning-doc ticket)
   ! --- same design as convolve_cubes.f90's own (see write_convolved_
   ! file's own comment there for the full single-writer-at-a-time
   ! rationale). Dispatched from inside write_reprojected_file's own
   ! `!$omp single` write phase (see there): the single thread kicks off
   ! the pthread and returns immediately, so the `single` construct's own
   ! implicit barrier only waits for the (fast) dispatch call, not the
   ! full write -- the NEXT block's read+resample then proceeds
   ! concurrently with the write running on the background thread.
   type :: block_write_job_t
      ! file_path/datastart (not out_unit): raw Fortran stream I/O at a
      ! computed byte offset (T19, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.
      ! md), not a live CFITSIO handle -- see write_reprojected_file's
      ! own comment for why. naxes_in/other_axes/other_idx/n_other are
      ! still needed (unlike convolve_cubes.f90's simpler single-freq-
      ! axis case) since reproject's output axis layout is general.
      character(len=1024) :: file_path = ' '
      integer(kind=8) :: datastart = 0_8
      integer :: naxes_in(max_axes) = 0, other_axes(max_axes) = 0
      integer :: other_idx(max_axes) = 0, n_other = 0
      integer :: chan_start = 0, chan_len = 0, nx = 0, ny = 0
      integer :: nwriters_eff = 1
      real, pointer :: data(:,:,:) => null()
   end type block_write_job_t
   type(block_write_job_t), target, save :: write_job
   integer(c_long) :: write_thread_id = 0
   logical :: write_pending = .false.
   logical :: write_failed = .false.

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

   character(len=16) :: mode
   character(len=512) :: reffile
   character(len=512) :: infiles(max_inputs)
   integer :: n_inputs, i

   character(len=512) :: cfgfile
   ! this_arg/cli_val/raw_cli_infiles hold a WHOLE comma-separated
   ! infiles= argument, not one path -- up to max_inputs (50) entries,
   ! each up to 512 chars (infiles(:)'s own per-entry length). 16384
   ! gives generous headroom over the worst case (~25,600 chars); 512
   ! here silently truncated a real run's own infiles= argument (found
   ! via scripts/run_pipeline.sh's symlink-redirection scheme, which
   ! lengthens every path enough to cross 512 for just 4 real files).
   character(len=16384) :: this_arg, cli_val
   character(len=512) :: cli_key
   character(len=16) :: cli_mode
   character(len=512) :: cli_reffile
   character(len=16384) :: raw_cli_infiles
   integer :: argc, iarg
   logical :: has_kv
   logical :: have_cfgfile
   logical :: seen_mode, seen_reffile, seen_infiles
   logical :: cli_seen_mode, cli_seen_reffile, cli_seen_infiles

   ! mem_frac_ram: optional (unlike mode/reffile/infiles, has a default
   ! and is never required) -- fraction of total system RAM (see
   ! get_mem_total_kb) budgeted for one read/resample/write block of
   ! planes at a time, same concept and default as rm_synthesis's own
   ! cfg%mem_frac_ram (rm_synthesis_mod.f90's plan_tile).
   real :: mem_frac_ram, cli_mem_frac_ram
   logical :: cli_seen_mem_frac_ram
   integer :: ios_mfr
   ! io_overlap (default n): background-thread block write, overlapped
   ! with the NEXT block's read+resample -- same scheme/key name as rm_
   ! synthesis/rmclean_cubes/convolve_cubes' own io_overlap (planning-doc
   ! ticket, added after a ~46GB end-to-end run measured
   ! this class of block write as serial with compute on real
   ! storage). Same CLI-then-cfg-then-override merge pattern as
   ! mem_frac_ram just above.
   logical :: io_overlap, cli_io_overlap
   logical :: cli_seen_io_overlap
   ! nwriters (default 1, T21 docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md):
   ! same key/clamp formula as rm_synthesis/rmclean_cubes' own nwriters
   ! (T7/T12) and convolve_cubes'/match_cubes' own copy of this option.
   ! Same CLI-then-cfg-then-override merge pattern as mem_frac_ram/
   ! io_overlap above.
   integer :: nwriters, cli_nwriters
   logical :: cli_seen_nwriters
   integer :: ios_nw
   ! dry_run (default n, T24 docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md):
   ! same key name as rm_synthesis' own dry_run -- checks the target
   ! output disk's rotational status and writes a suggested
   ! reproject_cubes_dryrun.cfg instead of processing any real file.
   ! Same CLI-then-cfg-then-override merge pattern as the others above.
   logical :: dry_run, cli_dry_run
   logical :: cli_seen_dry_run
   ! min_free_pct_after_write (default 10, UR-15 docs/dev/BRANCH_STATUS_
   ! update-reproject.md): percent of the output disk's total size that
   ! must still be free once every output is written, checked before
   ! any file is created. Same CLI-then-cfg-then-override merge pattern
   ! as the others above.
   double precision :: min_free_pct_after_write, cli_min_free_pct_after_write
   logical :: cli_seen_min_free_pct_after_write
   integer :: ios_mfp
   ! Logging & timing (planning-doc ticket) -- log_level/timing_enabled/
   ! log_output_file themselves declared near the top of this program
   ! (alongside the `use logging_mod`); only their CLI-override
   ! counterparts are declared here, same CLI-then-cfg-then-override
   ! merge pattern as mem_frac_ram/io_overlap above.
   character(len=16) :: cli_log_level
   logical :: cli_timing_enabled
   character(len=272) :: cli_log_output_file
   logical :: cli_seen_log_level, cli_seen_timing_enabled
   logical :: cli_seen_log_output_file
   ! interp/interp_params (docs/dev/REPROJECT_INTERP_OPTIONS_PLAN.md):
   ! same CLI-then-cfg-then-override merge pattern as mem_frac_ram/
   ! io_overlap above. interp_scheme/interp_params_arr are resolve_
   ! interp's own output (ast_interp_mod), resolved once here, not once
   ! per OpenMP thread/plane -- the same AST scheme constant and
   ! params(:) apply to every plane of a given run.
   character(len=16) :: interp, cli_interp
   character(len=256) :: interp_params, cli_interp_params
   logical :: cli_seen_interp, cli_seen_interp_params
   integer :: interp_scheme, status_interp
   double precision :: interp_params_arr(2)

   ! target_system/target_projection/target_center/target_lon/
   ! target_lat/target_cdelt_lon/target_cdelt_lat (docs/dev/
   ! BRANCH_STATUS_update-reproject.md UR-7): same CLI-then-cfg-then-
   ! override merge pattern as interp/interp_params above. reffile= is
   ! now optional and, when given, is exclusive with every one of
   ! these keys (resolve_target_grid enforces this) -- reffile=
   ! already fully determines the output grid, full stop. When reffile
   ! is absent, resolve_target_grid builds skymap_ref/skyframe_ref/
   ! pixaxes_ref/naxes_ref itself, entirely in memory (used_target_
   ! frame=.true.); the existing load_wcs(reffile,...)+extract_sky_
   ! mapping call further below is then skipped, not replaced.
   character(len=16) :: target_system, cli_target_system
   character(len=16) :: target_projection, cli_target_projection
   character(len=16) :: target_center, cli_target_center
   character(len=32) :: target_lon, cli_target_lon
   character(len=32) :: target_lat, cli_target_lat
   character(len=32) :: target_cdelt_lon, cli_target_cdelt_lon
   character(len=32) :: target_cdelt_lat, cli_target_cdelt_lat
   ! target_hpx_nside (UR-23): same CLI-then-cfg-then-override merge
   ! pattern as every other target_* key above -- see target_grid_mod.
   ! f90 for why it is only valid with target_projection=hpx.
   character(len=16) :: target_hpx_nside, cli_target_hpx_nside
   logical :: cli_seen_target_system, cli_seen_target_projection
   logical :: cli_seen_target_center
   logical :: cli_seen_target_lon, cli_seen_target_lat
   logical :: cli_seen_target_cdelt_lon, cli_seen_target_cdelt_lat
   logical :: cli_seen_target_hpx_nside
   integer :: status_target_grid
   logical :: used_target_grid
   ! The system and projection the no-reffile grid was built in, as
   ! resolve_target_grid returns them: target_system/target_projection
   ! as given can be blank (projection alone inherits the infiles'
   ! system; an omitted projection means sin), so every later call that
   ! rebuilds the grid takes these instead.
   character(len=16) :: grid_system, grid_projection
   ! Plain values, not AST object handles (docs/dev/
   ! BRANCH_STATUS_update-reproject.md UR-7) -- AST enforces per-thread
   ! object ownership on EVERY operation, confirmed by direct
   ! multi-threaded test (scratch/build_check/
   ! test_ast_copy_threads.f90), so no AST object built by this
   ! program's own serial setup may ever be touched by write_
   ! reprojected_file's own OpenMP worker threads. Each thread rebuilds
   ! its own Mapping from these numbers via target_grid_mod.f90's own
   ! build_target_skymap instead.
   double precision :: target_center_lon_deg, target_center_lat_deg
   double precision :: target_cdelt1_mag, target_cdelt2_mag
   double precision :: target_final_crpix(2)

   integer :: wcs_ref, skymap_ref, skyframe_ref
   integer :: naxes_ref(max_axes), pixaxes_ref(2)
   integer :: wcs_in, skymap_in, skyframe_in
   integer :: naxes_in(max_axes), pixaxes_in(2)

   ! --- UR-10 (docs/dev/BRANCH_STATUS_update-reproject.md): invalid-
   ! region info/warning, once per run on the final output grid ---
   integer :: ur10_skymap, ur10_skyframe, ur10_pixaxes(2), ur10_status
   ! UR-21: the output grid's cell size on the sky, and the run's own
   ! grid keys as match_cubes would take them, for the beam check.
   double precision :: out_cell1_arcsec, out_cell2_arcsec
   integer :: cell_status
   character(len=1024) :: grid_keys_for_advice
   double precision :: ur10_invalid_frac

   ! --- UR-15 (docs/dev/BRANCH_STATUS_update-reproject.md): unified
   ! overlap check + disk-space pre-check, both before any file write ---
   integer :: n_planes_per_infile(max_inputs)
   integer(kind=8) :: required_bytes, avail_bytes, total_bytes
   integer :: diskcheck_status

   ! --- Skip-if-already-matched (planning-doc ticket, same logic as
   ! match_cubes.f90's own -- see there for the full rationale) ---
   logical :: needs_processing(max_inputs), out_exists
   integer :: nx_out_common, ny_out_common
   integer :: ref_skip_unit, cand_skip_unit, fitsstat_skip, blocksize_skip
   double precision :: lbnd_out(2), ubnd_out(2)
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
   call check_and_clamp_omp_threads()

   have_cfgfile = .false.
   cli_seen_mode = .false.
   cli_seen_reffile = .false.
   cli_seen_infiles = .false.
   cli_seen_mem_frac_ram = .false.
   mem_frac_ram = 0.25
   cli_seen_io_overlap = .false.
   io_overlap = .false.
   cli_seen_nwriters = .false.
   nwriters = 1
   cli_seen_dry_run = .false.
   dry_run = .false.
   cli_seen_min_free_pct_after_write = .false.
   min_free_pct_after_write = 10.0d0
   cli_seen_log_level = .false.
   log_level = 'info'
   cli_seen_timing_enabled = .false.
   timing_enabled = .false.
   cli_seen_log_output_file = .false.
   log_output_file = ' '
   cli_seen_interp = .false.
   ! Default changed from linear to lanczos3, 2026-10-13 (UR-12
   ! benchmark, docs/dev/BRANCH_STATUS_update-reproject.md): mean RMS
   ! error ~6x lower than linear on realistic (smooth/compact) sources,
   ! confirmed against analytic ground truth and cross-checked against
   ! astropy's own reproject package. No interp_params= needed --
   ! lanczos3 fixes its own support radius (see ast_interp_mod.f90).
   interp = 'lanczos3'
   cli_seen_interp_params = .false.
   interp_params = ' '
   cli_seen_target_system = .false.
   target_system = ' '
   cli_seen_target_projection = .false.
   target_projection = ' '
   cli_seen_target_center = .false.
   target_center = ' '
   cli_seen_target_lon = .false.
   target_lon = ' '
   cli_seen_target_lat = .false.
   target_lat = ' '
   cli_seen_target_cdelt_lon = .false.
   target_cdelt_lon = ' '
   cli_seen_target_cdelt_lat = .false.
   target_cdelt_lat = ' '
   cli_seen_target_hpx_nside = .false.
   target_hpx_nside = ' '
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
         case ('mem_frac_ram')
            if (cli_seen_mem_frac_ram) then
               write(*,*) 'ERROR: mem_frac_ram given more than once on the command line'
               stop 1
            endif
            read(cli_val, *, iostat=ios_mfr) cli_mem_frac_ram
            if (ios_mfr.ne.0) then
               write(*,*) 'ERROR: mem_frac_ram must be a number, got "',&
               &trim(cli_val), '"'
               stop 1
            endif
            cli_seen_mem_frac_ram = .true.
         case ('io_overlap')
            if (cli_seen_io_overlap) then
               write(*,*) 'ERROR: io_overlap given more than once on the command line'
               stop 1
            endif
            cli_io_overlap = flag_from_value_reproject(cli_val)
            cli_seen_io_overlap = .true.
         case ('nwriters')
            if (cli_seen_nwriters) then
               write(*,*) 'ERROR: nwriters given more than once on the command line'
               stop 1
            endif
            read(cli_val, *, iostat=ios_nw) cli_nwriters
            if (ios_nw.ne.0 .or. cli_nwriters.lt.1) then
               write(*,*) 'ERROR: nwriters must be an integer >= 1, got "',&
               &trim(cli_val), '"'
               stop 1
            endif
            cli_seen_nwriters = .true.
         case ('dry_run')
            if (cli_seen_dry_run) then
               write(*,*) 'ERROR: dry_run given more than once on the command line'
               stop 1
            endif
            cli_dry_run = flag_from_value_reproject(cli_val)
            cli_seen_dry_run = .true.
         case ('min_free_pct_after_write')
            if (cli_seen_min_free_pct_after_write) then
               write(*,*) 'ERROR: min_free_pct_after_write given more than once on the command line'
               stop 1
            endif
            read(cli_val, *, iostat=ios_mfp) cli_min_free_pct_after_write
            if (ios_mfp.ne.0) then
               write(*,*) 'ERROR: min_free_pct_after_write must be a number, got "',&
               &trim(cli_val), '"'
               stop 1
            endif
            cli_seen_min_free_pct_after_write = .true.
         case ('log_level')
            if (cli_seen_log_level) then
               write(*,*) 'ERROR: log_level given more than once on the command line'
               stop 1
            endif
            cli_log_level = trim(cli_val)
            cli_seen_log_level = .true.
         case ('timing_enabled')
            if (cli_seen_timing_enabled) then
               write(*,*) 'ERROR: timing_enabled given more than once on the command line'
               stop 1
            endif
            cli_timing_enabled = flag_from_value_logging(cli_val)
            cli_seen_timing_enabled = .true.
         case ('log_output_file')
            if (cli_seen_log_output_file) then
               write(*,*) 'ERROR: log_output_file given more than once on the command line'
               stop 1
            endif
            cli_log_output_file = trim(cli_val)
            cli_seen_log_output_file = .true.
         case ('interp')
            if (cli_seen_interp) then
               write(*,*) 'ERROR: interp given more than once on the command line'
               stop 1
            endif
            cli_interp = trim(cli_val)
            cli_seen_interp = .true.
         case ('interp_params')
            if (cli_seen_interp_params) then
               write(*,*) 'ERROR: interp_params given more than once on the command line'
               stop 1
            endif
            cli_interp_params = trim(cli_val)
            cli_seen_interp_params = .true.
         case ('target_system')
            if (cli_seen_target_system) then
               write(*,*) 'ERROR: target_system given more than once on the command line'
               stop 1
            endif
            cli_target_system = trim(cli_val)
            cli_seen_target_system = .true.
         case ('target_projection')
            if (cli_seen_target_projection) then
               write(*,*) 'ERROR: target_projection given more than once on the command line'
               stop 1
            endif
            cli_target_projection = trim(cli_val)
            cli_seen_target_projection = .true.
         case ('target_center')
            if (cli_seen_target_center) then
               write(*,*) 'ERROR: target_center given more than once on the command line'
               stop 1
            endif
            cli_target_center = trim(cli_val)
            cli_seen_target_center = .true.
         case ('target_lon')
            if (cli_seen_target_lon) then
               write(*,*) 'ERROR: target_lon given more than once on the command line'
               stop 1
            endif
            cli_target_lon = trim(cli_val)
            cli_seen_target_lon = .true.
         case ('target_lat')
            if (cli_seen_target_lat) then
               write(*,*) 'ERROR: target_lat given more than once on the command line'
               stop 1
            endif
            cli_target_lat = trim(cli_val)
            cli_seen_target_lat = .true.
         case ('target_cdelt_lon')
            if (cli_seen_target_cdelt_lon) then
               write(*,*) 'ERROR: target_cdelt_lon given more than once on the command line'
               stop 1
            endif
            cli_target_cdelt_lon = trim(cli_val)
            cli_seen_target_cdelt_lon = .true.
         case ('target_cdelt_lat')
            if (cli_seen_target_cdelt_lat) then
               write(*,*) 'ERROR: target_cdelt_lat given more than once on the command line'
               stop 1
            endif
            cli_target_cdelt_lat = trim(cli_val)
            cli_seen_target_cdelt_lat = .true.
         case ('target_hpx_nside')
            if (cli_seen_target_hpx_nside) then
               write(*,*) 'ERROR: target_hpx_nside given more than once on the command line'
               stop 1
            endif
            cli_target_hpx_nside = trim(cli_val)
            cli_seen_target_hpx_nside = .true.
         case default
            write(*,*) 'ERROR: unrecognised key "', trim(cli_key), '" -- expected',&
            &' mode, reffile, infiles, mem_frac_ram, io_overlap, nwriters,'//&
            &' dry_run, min_free_pct_after_write, log_level, timing_enabled,'//&
            &' log_output_file, interp,'//&
            &' interp_params, target_system, target_projection, target_center,'//&
            &' target_lon, target_lat, target_cdelt_lon, target_cdelt_lat, or'//&
            &' target_hpx_nside'
            stop 1
         end select
         iarg = iarg + 1
      endif
   enddo

   seen_mode = .false.
   seen_reffile = .false.
   seen_infiles = .false.
   if (have_cfgfile) then
      call read_reproject_cfg(cfgfile, mode, reffile, infiles, n_inputs,&
      &mem_frac_ram, io_overlap, nwriters, dry_run, min_free_pct_after_write,&
      &log_level, timing_enabled, log_output_file,&
      &interp, interp_params, target_system, target_projection, target_center,&
      &target_lon, target_lat, target_cdelt_lon, target_cdelt_lat, target_hpx_nside,&
      &seen_mode, seen_reffile, seen_infiles, status)
      if (status.ne.0) stop 1
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
   if (cli_seen_mem_frac_ram) mem_frac_ram = cli_mem_frac_ram
   if (cli_seen_io_overlap) io_overlap = cli_io_overlap
   if (cli_seen_nwriters) nwriters = cli_nwriters
   if (cli_seen_dry_run) dry_run = cli_dry_run
   if (cli_seen_min_free_pct_after_write) min_free_pct_after_write = cli_min_free_pct_after_write
   if (cli_seen_log_level) log_level = cli_log_level
   if (cli_seen_timing_enabled) timing_enabled = cli_timing_enabled
   if (cli_seen_log_output_file) log_output_file = cli_log_output_file
   if (cli_seen_interp) interp = cli_interp
   if (cli_seen_interp_params) interp_params = cli_interp_params
   if (cli_seen_target_system) target_system = cli_target_system
   if (cli_seen_target_projection) target_projection = cli_target_projection
   if (cli_seen_target_center) target_center = cli_target_center
   if (cli_seen_target_lon) target_lon = cli_target_lon
   if (cli_seen_target_lat) target_lat = cli_target_lat
   if (cli_seen_target_cdelt_lon) target_cdelt_lon = cli_target_cdelt_lon
   if (cli_seen_target_cdelt_lat) target_cdelt_lat = cli_target_cdelt_lat
   if (cli_seen_target_hpx_nside) target_hpx_nside = cli_target_hpx_nside

   ! infiles is always required. mode/reffile (UR-7: reffile= is now
   ! optional) must be given together -- reffile= fully determines the
   ! output grid, full stop -- or both omitted, in which case
   ! target_system=/target_projection= build one instead
   ! (resolve_target_grid, below, enforces the actual mutual-
   ! exclusivity with every target_* key).
   if (.not. seen_infiles) then
      call print_usage()
      stop 1
   endif
   if (seen_mode .neqv. seen_reffile) then
      write(*,*) 'ERROR: mode and reffile must be given together, or both'//&
      &' omitted (then target_system=/target_projection= build the output'//&
      &' grid instead)'
      stop 1
   endif

   if (mem_frac_ram.le.0.0 .or. mem_frac_ram.gt.0.95) then
      write(*,*) 'ERROR: mem_frac_ram must be > 0 and <= 0.95, got ', mem_frac_ram
      stop 1
   endif
   if (min_free_pct_after_write.lt.0.0d0 .or. min_free_pct_after_write.gt.100.0d0) then
      write(*,*) 'ERROR: min_free_pct_after_write must be between 0 and 100 (percent), got ',&
      &min_free_pct_after_write
      stop 1
   endif

   ! interp/interp_params (docs/dev/REPROJECT_INTERP_OPTIONS_PLAN.md):
   ! resolved once here (not once per OpenMP thread/plane below -- the
   ! same scheme/params apply to every plane of this run), via the
   ! module shared with match_cubes.f90's own 3 independent resampling
   ! call sites, so the name-to-AST-constant mapping cannot drift
   ! between the two files the way AST's own SincSinc/SincCos drifted
   ! from each other.
   status_interp = 0
   call resolve_interp(interp, interp_params, interp_scheme,&
   &interp_params_arr, status_interp)
   if (status_interp.ne.0) stop 1

   ! target_system/target_projection/target_center/target_lon/
   ! target_lat/target_cdelt_lon/target_cdelt_lat (docs/dev/
   ! BRANCH_STATUS_update-reproject.md UR-7): resolved once here. When
   ! reffile is absent, returns the PLAIN VALUES (not AST object
   ! handles -- see this program's own declaration comment above)
   ! describing the fresh grid, plus pixaxes_ref=(1,2)/naxes_ref --
   ! the load_wcs(reffile,...)+extract_sky_mapping call further below
   ! is then skipped outright. mode is forced to 'reference' in that
   ! case: there is no reffile for union/intersection to be relative
   ! to, since the built grid is already exactly the infiles' own
   ! combined footprint.
   status_target_grid = 0
   call resolve_target_grid(seen_reffile, target_system, target_projection,&
   &target_center, target_lon, target_lat, target_cdelt_lon, target_cdelt_lat,&
   &target_hpx_nside, infiles, n_inputs, used_target_grid, grid_system, grid_projection,&
   &target_center_lon_deg,&
   &target_center_lat_deg, target_cdelt1_mag, target_cdelt2_mag,&
   &target_final_crpix, pixaxes_ref, naxes_ref, n_planes_per_infile,&
   &status_target_grid)
   if (status_target_grid.ne.0) stop 1
   if (used_target_grid) mode = 'reference'

   if (trim(mode).ne.'intersection' .and. trim(mode).ne.'union' .and.&
   &trim(mode).ne.'reference') then
      write(*,*) 'ERROR: mode must be intersection, union, or reference'
      stop 1
   endif

   call init_logging(log_level, timing_enabled, log_output_file, status)
   if (status.ne.0) then
      write(*,*) 'ERROR: cannot open log_output_file: ', trim(log_output_file)
      stop 1
   endif
   call log_message('info', 'startup', 'reproject_cubes run started')

   ! dry_run=y (UR-15): the rotational-disk suggestion below needs no
   ! grid at all, but reporting required-vs-available disk space does --
   ! so the grid itself (and the unified overlap check that produces it)
   ! must run before dry_run's own early stop, not after. Moved up from
   ! this program's own earlier structure for exactly that reason.
   status = 0
   call ast_begin(status)

   if (.not. used_target_grid) then
      call load_wcs(reffile, wcs_ref, naxes_ref, status)
      call extract_sky_mapping(wcs_ref, skymap_ref, skyframe_ref, pixaxes_ref, status)
      if (status.ne.0) then
         write(*,*) 'ERROR: failed to load the reference file''s WCS'
         stop 1
      endif
   endif

   ! Output grid starts as the reference's own full extent; intersection
   ! shrinks it, union grows it, reference mode leaves it untouched (the
   ! unified check below, UR-15, never folds for reference mode).
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

   ! UR-15 (docs/dev/BRANCH_STATUS_update-reproject.md): unified overlap
   ! sanity-check across every mode, including the reference-mode check
   ! that never existed before -- see check_overlap_and_size_grid's own
   ! header for the collect-every-offender-before-reporting behaviour
   ! and the daisy-chain-preserving fold. used_target_grid already ran
   ! this internally (resolve_target_grid_noref) and populated n_planes_
   ! per_infile; only run it again here when a real reffile was given.
   if (.not. used_target_grid) then
      call check_overlap_and_size_grid(skymap_ref, skyframe_ref, infiles,&
      &n_inputs, mode, .true., lbnd_out, ubnd_out, n_planes_per_infile, status)
      if (status.ne.0) stop 1
   endif

   if (lbnd_out(1).gt.ubnd_out(1) .or. lbnd_out(2).gt.ubnd_out(2)) then
      write(*,*) 'ERROR: computed output grid is empty (', trim(mode), ' mode)'
      stop 1
   endif
   nx_out_common = nint(ubnd_out(1) - lbnd_out(1)) + 1
   ny_out_common = nint(ubnd_out(2) - lbnd_out(2)) + 1

   ! === UR-15: disk-space pre-check, before any file is written. Runs
   ! even during dry_run=y, so its own report can be included below. ===
   call check_disk_space(infiles(1), n_inputs, n_planes_per_infile,&
   &nx_out_common, ny_out_common, min_free_pct_after_write, required_bytes,&
   &avail_bytes, total_bytes, diskcheck_status)

   if (dry_run) then
      call run_dry_run_check(infiles(1), 'reproject_cubes')
      if (avail_bytes.ge.0_8) then
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

   ! === UR-10: invalid-region info/warning on the final output grid,
   ! once per run, before any plane is resampled -- a facet seam (e.g.
   ! -HPX) or any other projection's own undefined region can leave a
   ! part of the output with no valid sky mapping at all, for reasons
   ! having nothing to do with missing input data coverage. used_
   ! target_grid has no serial skymap_ref of its own at this point
   ! (every other caller builds its own per-thread copy) -- build one
   ! here, single-threaded, solely for this check, then annul it. ===
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
      &lbnd_out(2)-1.0d0, nx_out_common, ny_out_common, ur10_invalid_frac,&
      &ur10_status)
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

   ! === UR-21: refuse an output grid that cannot sample an infile's beam,
   ! before any output is written (see target_grid_mod's
   ! check_beam_sampling). The advice quotes this run's grid as
   ! match_cubes keys, so it can be rerun as printed. ===
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
      grid_keys_for_advice = 'footprint_mode='//trim(mode)//' reffile='//trim(reffile)
   endif
   cell_status = 0
   call check_beam_sampling(infiles, n_inputs, out_cell1_arcsec, out_cell2_arcsec,&
   &grid_keys_for_advice, cell_status)
   if (cell_status.ne.0) stop 1

   ! === UR-17: cross-file beam advisory. This tool changes the pixel
   ! grid, not the angular resolution, so infiles with different beams
   ! produce outputs with different beams. ===
   call advise_beam_consistency(infiles, n_inputs)

   ! === Skip-if-already-matched pre-flight (planning-doc ticket) ===
   ! Same scheme as match_cubes.f90's own (see there for the full
   ! rationale): every safety check and skip decision for the WHOLE
   ! batch happens up front, before any file is processed, so a bad run
   ! (a stale output already on disk) fails fast. A pre-existing output
   ! path is always refused, never silently reused or overwritten,
   ! regardless of what this run's own skip decision would have been.
   do i = 1, n_inputs
      inquire(file=trim(strip_fits_ext(infiles(i)))//'_REPROJ.FITS', exist=out_exists)
      if (out_exists) then
         write(*,*) 'ERROR: output path already exists, refusing to proceed'//&
         &' (stale output from a previous run? remove it first): ',&
         &trim(strip_fits_ext(infiles(i)))//'_REPROJ.FITS'
         stop 1
      endif
   enddo

   ! Skip-if-already-matched needs a real reffile file to compare
   ! against (sky_wcs_matches_target reads its CTYPE/CRVAL/etc via
   ! CFITSIO) -- not available when used_target_grid (no reffile, a
   ! freshly-built in-memory grid). Skipped in that case: a pure
   ! performance optimisation, not a correctness requirement, and the
   ! freshly-built grid is derived from these exact infiles' own
   ! footprints, so an exact pre-existing match would be a rare
   ! coincidence anyway.
   if (.not. used_target_grid) then
      fitsstat_skip = 0
      call safe_ftopen(ref_skip_unit, trim(reffile), 0, blocksize_skip, fitsstat_skip)
      if (fitsstat_skip.ne.0) then
         write(*,*) 'ERROR: cannot reopen reference file for the geometry check: ',&
         &trim(reffile)
         stop 1
      endif

      do i = 1, n_inputs
         call load_wcs(infiles(i), wcs_in, naxes_in, status)
         call extract_sky_mapping(wcs_in, skymap_in, skyframe_in, pixaxes_in, status)
         if (status.ne.0) then
            write(*,*) 'ERROR: failed to read input''s WCS for the geometry check: ',&
            &trim(infiles(i))
            stop 1
         endif
         fitsstat_skip = 0
         call safe_ftopen(cand_skip_unit, trim(infiles(i)), 0, blocksize_skip, fitsstat_skip)
         if (fitsstat_skip.ne.0) then
            write(*,*) 'ERROR: cannot reopen input for the geometry check: ',&
            &trim(infiles(i))
            stop 1
         endif
         needs_processing(i) = .not. sky_wcs_matches_target(cand_skip_unit,&
         &pixaxes_in(1), pixaxes_in(2), naxes_in(pixaxes_in(1)), naxes_in(pixaxes_in(2)),&
         &ref_skip_unit, pixaxes_ref(1), pixaxes_ref(2), lbnd_out(1)-1.0d0,&
         &lbnd_out(2)-1.0d0, nx_out_common, ny_out_common)
         call safe_ftclos(cand_skip_unit, fitsstat_skip)
         call ast_annul(skymap_in, status)
         call ast_annul(skyframe_in, status)
         call ast_annul(wcs_in, status)
         if (.not. needs_processing(i)) then
            write(*,'(A,A,A)') 'SKIP: ', trim(infiles(i)),&
            &' already matches the target grid -- no output written, use it directly'
         endif
      enddo
      fitsstat_skip = 0
      call safe_ftclos(ref_skip_unit, fitsstat_skip)
   else
      needs_processing(1:n_inputs) = .true.
   endif

   ! --- Resample and write every plane of every input onto the final
   ! output grid, via astResampleR + FTPSSE ---
   ! Only naxes_in/pixaxes_in (plain per-file array-shape data, not an
   ! AST Object) are needed here -- write_reprojected_file builds its
   ! own private input->reference Mapping per OpenMP thread internally
   ! (see its own comment for why), so this loop no longer needs to
   ! derive map_in2ref itself the way it used to.
   do i = 1, n_inputs
      if (.not. needs_processing(i)) cycle
      call load_wcs(infiles(i), wcs_in, naxes_in, status)
      call extract_sky_mapping(wcs_in, skymap_in, skyframe_in, pixaxes_in, status)
      if (status.ne.0) then
         write(*,*) 'ERROR: failed to read input''s WCS for resampling: ',&
         &trim(infiles(i))
         stop 1
      endif
      ! Plain filename (NOT '!'-prefixed): write_reprojected_file's own
      ! FTINIT fails if the output already exists rather than silently
      ! clobbering it. Fixed bug, found while implementing match_cubes'
      ! own skip-if-already-matched feature (planning-doc ticket) and
      ! auditing every FTINIT call site for the same pattern: this call
      ! previously baked the '!'-prefix CLOBBER convention into the
      ! output filename itself, silently deleting and overwriting a
      ! pre-existing output with no warning.
      call write_reprojected_file(reffile, infiles(i),&
      &trim(strip_fits_ext(infiles(i)))//'_REPROJ.FITS', pixaxes_ref,&
      &naxes_in, pixaxes_in, lbnd_out, ubnd_out, mem_frac_ram, io_overlap,&
      &nwriters, interp_scheme, interp_params_arr, used_target_grid,&
      &grid_system, grid_projection, target_center_lon_deg, target_center_lat_deg,&
      &target_cdelt1_mag, target_cdelt2_mag, target_final_crpix, status)
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

   if (.not. used_target_grid) then
      call ast_annul(skymap_ref, status)
      call ast_annul(skyframe_ref, status)
      call ast_annul(wcs_ref, status)
   endif
   call ast_end(status)

   if (status.ne.0) then
      write(*,*) 'ERROR: AST reported an error, final status=', status
      stop 1
   endif

   write(*,*) 'OK: footprint-mode output grid computed successfully.'
   call timer_report_summary()
   call log_message('info', 'finalize', 'reproject_cubes run completed')

contains

   subroutine write_reprojected_file(reffile, infile, outfile, pixaxes_ref,&
   &naxes_in, pixaxes_in, lbnd_out_d, ubnd_out_d, mem_frac_ram, io_overlap_l,&
   &nwriters, interp_scheme, interp_params_arr, used_target_grid,&
   &target_system_arg, target_projection_arg, target_center_lon_arg,&
   &target_center_lat_arg,&
   &target_cdelt1_arg, target_cdelt2_arg, target_final_crpix_arg, status)
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
      !! Planes are read, resampled, and written in BLOCKS, not one at a
      !! time: block size comes from mem_frac_ram (get_mem_total_kb
      !! below mirrors rm_synthesis's own /proc/meminfo MemTotal read
      !! exactly, same concept as its plan_tile). Each block goes through
      !! three strictly separated phases -- read (one thread, via OMP
      !! `single`), resample (every thread, in parallel, one
      !! astResampleR call per plane in the block), write (one thread,
      !! via `single`) -- relying on the implicit barrier every `single`/
      !! `do` construct has by default. This replaced an earlier
      !! plane-at-a-time version that serialised every CFITSIO call
      !! behind an OMP critical section: with blocks, a whole block's I/O
      !! happens on a single thread BY CONSTRUCTION, so there is no
      !! concurrent CFITSIO access to guard against in the first place --
      !! nothing to lock, rather than a lock relied on to make concurrent
      !! access safe. Batching also amortises CFITSIO's per-call
      !! (FTOPEN/FTGSVE/FTPSSE) overhead across many planes instead of
      !! paying it every single one.
      use, intrinsic :: ieee_arithmetic
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_wtime
      character(len=*), intent(in) :: reffile, infile, outfile
      integer, intent(in) :: pixaxes_ref(2)
      integer, intent(in) :: naxes_in(:), pixaxes_in(2)
      double precision, intent(in) :: lbnd_out_d(2), ubnd_out_d(2)
      real, intent(in) :: mem_frac_ram
      logical, intent(in) :: io_overlap_l
      integer, intent(in) :: nwriters
      integer, intent(in) :: interp_scheme
      double precision, intent(in) :: interp_params_arr(2)
      logical, intent(in) :: used_target_grid
      character(len=*), intent(in) :: target_system_arg, target_projection_arg
      double precision, intent(in) :: target_center_lon_arg, target_center_lat_arg
      double precision, intent(in) :: target_cdelt1_arg, target_cdelt2_arg
      double precision, intent(in) :: target_final_crpix_arg(2)
      integer, intent(inout) :: status

      ! interp_finterp: Tier 2 (AST__UKERN1, REPROJECT_INTERP_OPTIONS_
      ! PLAN.md section 6d) -- deliberately a LOCAL variable, resolved
      ! from the plain-integer interp_scheme dummy argument just below
      ! (outside the parallel region), never received as a dummy
      ! argument itself. Passing a procedure(),pointer as THIS
      ! subroutine's own dummy argument, shared into its own !$omp
      ! parallel region below, segfaulted inside AST -- isolated to a
      ! gfortran code-generation issue specific to that one
      ! combination (contains-ed procedure + procedure-pointer dummy
      ! argument + OMP sharing), not the finterp mechanism itself
      ! (verified correct independently -- see section 6b/6c/6d and
      ! scratch/build_check/test_ukern1_internal_proc_bug.f90).
      procedure(), pointer :: interp_finterp

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
      character(len=80) :: comment

      integer(kind=8) :: mem_total_kb, bytes_per_plane, mem_safe_bytes
      integer(kind=8) :: block_planes64
      ! interp=exact: the overlap weights, built once for this file and
      ! shared read-only by every thread while planes are applied.
      logical :: is_exact
      type(exact_weights_t) :: exact_wts
      integer(kind=8) :: exact_est_bytes
      double precision :: in_cell1, in_cell2, out_cell1, out_cell2
      integer :: cell_status, n_bands, iband, exact_lbnd(2), exact_ubnd(2)
      integer :: t_map_in2sky, t_lonaxis, t_sky2sky
      integer, parameter :: exact_band_rows = 8
      integer :: block_planes, n_groups, igroup, axis1_extent
      integer :: chan_start, chan_len, local_iplane
      real, allocatable, target :: block_data_in(:,:,:), block_data_out(:,:,:,:)

      ! Per-OpenMP-thread private AST working set (see the parallel
      ! region below, T20, for why each thread still builds its own,
      ! and why the creation step -- not the objects themselves -- is
      ! what needs the critical section).
      integer :: t_status, t_wcs_ref, t_skymap_ref, t_skyframe_ref
      integer :: t_naxes_ref(max_axes), t_pixaxes_ref(2)
      integer :: t_wcs_in, t_skymap_in, t_skyframe_in
      integer :: t_naxes_in(max_axes), t_pixaxes_in(2)
      integer :: t_map_in2ref

      ! Per-plane resample working set (private, one astResampleR call
      ! per plane inside a block's parallel do).
      integer :: lbnd_in(2), ubnd_in(2), lbnd_o(2), ubnd_o(2), nbad
      real :: badval
      real(dp) :: t_stage
      ! Per-thread swim-lane instrumentation (planning-doc ticket) -- see
      ! convolve_cubes.f90's own write_convolved_file for the full
      ! rationale, and match_cubes.f90's own process_one_file_general
      ! for why iblock is PRIVATE here (persistent parallel region, every
      ! thread redundantly computes the same chan_start/chan_len
      ! progression).
      integer :: iblock, tid_local
      real(dp) :: t_thread_start, t_thread_elapsed
      character(len=160) :: thread_msg

      if (status.ne.0) return
      call timer_reset_file_stages()
      call log_message('info', 'reproject', 'starting: '//trim(infile))

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

      nx_in = naxes_in(pixaxes_in(1))
      ny_in = naxes_in(pixaxes_in(2))
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

      ! Sky axes (output 1,2): WCS copied from the REFERENCE's own
      ! pixaxes_ref-numbered keywords, CRPIX shifted for the output
      ! grid's own origin (same formula rm_synthesis.f90 already uses for
      ! its own subimage CRPIX shift: (old_crpix - pixel_offset) + 1,
      ! stride 1 since this is a crop/grow, never a sub-sample). When
      ! used_target_grid (UR-7, no reffile given), there is no file to
      ! reopen -- target_*_arg carry the plain values (system, centre,
      ! cdelt, final CRPIX) that fully determine the final Mapping, and
      ! write_target_wcs_cards (target_grid_mod.f90) rebuilds that
      ! Mapping itself (via build_target_skymap, safe to call from this
      ! thread since it owns every AST object it creates) before
      ! deriving CTYPE/CRVAL/CRPIX/CDELT/PC/EQUINOX/RADESYS from it via
      ! AST's own FitsChan writer in one call, covering what
      ! copy_axis_keywords/copy_sky_rotation_matrix/
      ! copy_wcs_system_keywords do between them for the file-based case.
      if (used_target_grid) then
         call write_target_wcs_cards(target_system_arg, target_projection_arg,&
         &target_center_lon_arg,&
         &target_center_lat_arg, target_cdelt1_arg, target_cdelt2_arg,&
         &target_final_crpix_arg, 1, 2, out_unit, status)
      else
         fitsstat = 0
         call safe_ftopen(ref_unit, trim(reffile), 0, blocksize, fitsstat)
         call copy_axis_keywords(ref_unit, pixaxes_ref(1), out_unit, 1,&
         &lbnd_out_d(1)-1.0d0, status)
         call copy_axis_keywords(ref_unit, pixaxes_ref(2), out_unit, 2,&
         &lbnd_out_d(2)-1.0d0, status)
         ! PCi_j/CDi_j (the 2x2-matrix alternative to CROTA, which
         ! copy_axis_keywords just above already handles per-axis): also
         ! from the reference, same reasoning as CROTA -- the output grid's
         ! orientation is the reference's own, unaffected by cropping/
         ! growing its extent.
         call copy_sky_rotation_matrix(ref_unit, pixaxes_ref(1),&
         &pixaxes_ref(2), out_unit, status)
         ! EQUINOX/RADESYS describe the coordinate SYSTEM, not a specific
         ! axis or the data -- the output is expressed in the REFERENCE's
         ! sky frame (compose_pix2pix converts every input's sky coordinates
         ! into it), so these follow reffile, not infile, same reasoning as
         ! the sky axes themselves just above.
         call copy_wcs_system_keywords(ref_unit, out_unit, status)
         call safe_ftclos(ref_unit, fitsstat)
      endif

      ! Other axes (output 3..): WCS copied from INFILE's own axis
      ! numbering unchanged (no CRPIX shift -- reprojection never touches
      ! these axes).
      fitsstat = 0
      call safe_ftopen(ref_unit, trim(infile), 0, blocksize, fitsstat)
      do k = 1, n_other
         call copy_axis_keywords(ref_unit, other_axes(k), out_unit,&
         &2+k, 0.0d0, status)
      enddo
      ! Everything else non-structural and non-axis-indexed (BUNIT,
      ! BMAJ/BMIN/BPA, OBJECT, TELESCOP, DATE-OBS, RESTFRQ, HISTORY,
      ! COMMENT, ...) describes the DATA, unchanged by reprojection --
      ! follows infile, matching the other axes just above. Previously
      ! dropped (only per-axis WCS keywords were ever copied),
      ! silently losing units/beam/provenance from the output.
      call copy_generic_header(ref_unit, out_unit, status)
      if (used_target_grid) then
         call FTPHIS(out_unit, 'reproject_cubes: reprojected from '//&
         &trim(infile)//' onto a '//trim(target_system_arg)//'/'//&
         &trim(target_projection_arg)//' target grid (no reffile)', fitsstat)
      else
         call FTPHIS(out_unit, 'reproject_cubes: reprojected from '//&
         &trim(infile)//' onto the grid of '//trim(reffile), fitsstat)
      endif
      call safe_ftclos(ref_unit, fitsstat)

      ! CASAMBM/BEAMS: reprojection is a pure spatial resample of each
      ! plane -- it never touches the beam itself (BMAJ/BMIN/BPA are
      ! stored in sky degrees, not pixels), so a per-channel
      ! BEAMS table on infile is still exactly correct for the
      ! reprojected output and should follow it unchanged. The scalar
      ! CASAMBM keyword itself already rode along verbatim as a raw
      ! header card inside copy_generic_header just above; that copy
      ! only ever touches the PRIMARY header, though, and can't reach
      ! the separate BEAMS extension HDU the keyword refers to -- this
      ! attaches that extension explicitly, matching the ftcopy pattern
      ! rm_synthesis.f90 already uses for its own CASAMBM/BEAMS
      ! passthrough.
      casambm_status = 0
      call ftgkyl(out_unit, 'CASAMBM', casambm_val, comment, casambm_status)
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
         call safe_ftclos(beams_unit, beams_status)
      endif

      ! T19 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): CFITSIO is not
      ! guaranteed thread-safe even across different file units/handles
      ! without a specific reentrant build -- confirmed the hard way, a
      ! real, timing-dependent hang on real WALLABY+EMU data (found via
      ! match_cubes.f90's own identical design): the main thread's own
      ! next-block read blocked concurrently with the background write
      ! thread's own write, both on what is almost certainly an internal
      ! CFITSIO lock. Fetch this file's own pixel-data byte offset ONCE,
      ! right here (all header writes, including the BEAMS extension
      ! copy just above, are already done), then close out_unit
      ! immediately -- CFITSIO's job for this file is done. Every pixel
      ! write from here on goes through write_one_block_raw (plain
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

      ! --- Block size: mem_frac_ram fraction of total system RAM,
      ! divided by the bytes one plane's worth of input+output costs
      ! (same budgeting concept as rm_synthesis's plan_tile, applied to
      ! "planes in a block" instead of "output pixels in a tile"). Never
      ! larger than a single cycle of other_axes(1) (the fastest-varying
      ! non-sky axis) -- a block always covers one CONTIGUOUS range of
      ! that axis with every slower other axis held fixed (see the
      ! group/block loop below), so it can never usefully span past
      ! where that axis wraps.
      n_planes = 1
      do k = 1, n_other
         n_planes = n_planes * naxes_in(other_axes(k))
      enddo
      call get_mem_total_kb(mem_total_kb)
      ! block_data_out's own term is x2 (not x1) -- double-buffered for
      ! io_overlap (see below), budgeted whether or not it's on.
      bytes_per_plane = int(4,8) * (int(nx_in,8)*int(ny_in,8) +&
      &2_8*int(nx_out,8)*int(ny_out,8))
      mem_safe_bytes = int(real(mem_frac_ram,8) * real(mem_total_kb,8) *&
      &1024.0d0, 8)
      is_exact = interp_scheme.eq.interp_exact
      if (is_exact) then
         ! The weight table is held for the whole file, so it comes out of
         ! the same budget as the blocks. Estimated from the output size and
         ! the ratio of output to input pixel size (both measured at their
         ! reference pixels, as the beam check does).
         cell_status = 0
         call fits_reference_cell_arcsec(infile, in_cell1, in_cell2, cell_status)
         if (used_target_grid) then
            out_cell1 = target_cdelt1_arg * 3600.0d0; out_cell2 = target_cdelt2_arg * 3600.0d0
         else
            call fits_reference_cell_arcsec(reffile, out_cell1, out_cell2, cell_status)
         endif
         if (cell_status.ne.0) then
            status = -1
            return
         endif
         exact_est_bytes = exact_estimate_bytes(nx_out, ny_out,&
         &max(out_cell1, out_cell2) / min(in_cell1, in_cell2))
         write(*,'(A)') 'interp=exact: overlap weights estimated at '//&
         &trim(mib_str(exact_est_bytes))//' MiB, of a mem_frac_ram budget of '//&
         &trim(mib_str(mem_safe_bytes))//' MiB'
         if (exact_est_bytes.ge.mem_safe_bytes) then
            write(*,*) 'ERROR: interp=exact needs more memory for its overlap'//&
            &' weights than mem_frac_ram allows; raise mem_frac_ram or use a'//&
            &' smaller output grid.'
            status = -1
            return
         endif
         mem_safe_bytes = mem_safe_bytes - exact_est_bytes
      endif
      block_planes64 = max(1_8, mem_safe_bytes / bytes_per_plane)
      block_planes64 = min(block_planes64, max_elements_per_block /&
      &max(1_8, max(int(nx_in,8)*int(ny_in,8), int(nx_out,8)*int(ny_out,8))))
      block_planes64 = max(1_8, block_planes64)
      block_planes = int(min(block_planes64, int(n_planes,8)))
      axis1_extent = 1
      if (n_other.ge.1) then
         axis1_extent = naxes_in(other_axes(1))
         block_planes = min(block_planes, axis1_extent)
      endif
      if (block_planes.lt.1) block_planes = 1

      write(*,'(A,A,A,I0,A,I0,A)') 'Writing ', trim(outfile), ': ',&
      &n_planes, ' plane(s), in blocks of up to ', block_planes, ' plane(s)'

      ! Parallelism only happens WITHIN a block (its plane range, via the
      ! `!$omp do` below) -- blocks themselves are strictly sequential,
      ! never overlapped/pipelined, so block_planes is a hard ceiling on
      ! how many threads can ever do useful work at once, not just a
      ! memory knob. A too-small mem_frac_ram therefore doesn't just add
      ! CFITSIO call overhead -- it can silently throw away most of the
      ! OpenMP speedup: measured 12.0s (block_planes=1, forced to 1
      ! thread) vs. 5.2s (one block covering the whole cube, 16 threads)
      ! on the same 1024x1024x300 synthetic cube, over 2x slower from
      ! losing parallelism, not from smaller reads/writes. Warn rather
      ! than silently eat that cost or override the user's own budget.
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
         n_groups = n_groups * naxes_in(other_axes(k))
      enddo

      cur_slot = 0
      write_pending = .false.
      write_failed = .false.
      status_par = 0
      nthreads = max(1, min(omp_get_max_threads(), block_planes))
      ! interp=exact builds its weights over bands of output rows, not
      ! planes, so a single-plane file still uses every thread for that.
      if (is_exact) nthreads = max(1, omp_get_max_threads())
      exact_lbnd(1) = nint(lbnd_out_d(1)); exact_lbnd(2) = nint(lbnd_out_d(2))
      exact_ubnd(1) = nint(ubnd_out_d(1)); exact_ubnd(2) = nint(ubnd_out_d(2))
      n_bands = (ny_out + exact_band_rows - 1) / exact_band_rows
      nwriters_eff = max(1, min(nwriters, omp_get_max_threads()))
      if (interp_scheme.eq.ast__ukern1) then
         interp_finterp => lanczos_kernel
      else
         interp_finterp => ast_null
      endif
      !$omp parallel num_threads(nthreads) default(none)&
      !$omp& shared(infile, reffile, outfile, naxes_in, pixaxes_in, other_axes,&
      !$omp& n_other, lbnd_out_d, ubnd_out_d, datastart, status_par,&
      !$omp& nx_in, ny_in, nx_out, ny_out, block_planes, block_data_in,&
      !$omp& block_data_out, n_groups, axis1_extent, cur_slot,&
      !$omp& io_overlap_l, write_dispatched_ok, write_pending,&
      !$omp& write_thread_id, write_failed, write_job, t_stage, nwriters_eff,&
      !$omp& interp_scheme, interp_params_arr, interp_finterp, used_target_grid,&
      !$omp& target_system_arg, target_projection_arg,&
      !$omp& target_center_lon_arg, target_center_lat_arg,&
      !$omp& target_cdelt1_arg, target_cdelt2_arg, target_final_crpix_arg,&
      !$omp& is_exact, exact_wts, exact_lbnd, exact_ubnd, n_bands)&
      !$omp& private(t_status, t_wcs_ref, t_skymap_ref, t_skyframe_ref,&
      !$omp& t_naxes_ref, t_pixaxes_ref, t_wcs_in, t_skymap_in, t_skyframe_in,&
      !$omp& t_naxes_in, t_pixaxes_in, t_map_in2ref, other_idx, remainder,&
      !$omp& radix, k, igroup, chan_start, chan_len, local_iplane, nbad,&
      !$omp& lbnd_in, ubnd_in, lbnd_o, ubnd_o, badval,&
      !$omp& iblock, tid_local, t_thread_start, t_thread_elapsed, thread_msg,&
      !$omp& iband, t_map_in2sky, t_lonaxis, t_sky2sky)

      ! T20 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): each thread still
      ! builds its OWN private input->reference pixel Mapping from
      ! scratch (own ast_begin context, own load_wcs/extract_sky_
      ! mapping/compose_pix2pix calls) -- an earlier attempt to build
      ! ONE Mapping and hand out ast_copy'd clones to every thread was
      ! tried and reverted: AST enforces per-thread object OWNERSHIP at
      ! runtime regardless of how the object was obtained (confirmed
      ! directly -- ast_resampler raised "Invalid Object pointer given
      ! ... currently owned by another thread" for every thread other
      ! than the one that ran ast_copy itself). SUN/211 Sec 4.12 ("AST
      ! Objects within Multi-threaded Applications") documents the
      ! correct handoff as astUnlock (by the creating thread) then
      ! astLock (by the receiving thread) around any ast_copy -- but
      ! this Fortran binding exports neither (checked directly against
      ! libstarlink_ast.so.9's own symbol table: no ast_lock_/
      ! ast_unlock_/ast_thread_ at all), so that handoff isn't
      ! available here, and per-thread private objects are the only
      ! option this binding supports.
      !
      ! What DOES need protecting: the earlier, real WALLABY+EMU
      ! deadlock (found via a live gdb backtrace: all 6 threads parked
      ! on two of libstarlink_ast.so.9's own internal mutexes, all
      ! inside ast_read_) was never about object OWNERSHIP -- it was 6
      ! threads calling ast_read_ simultaneously, contending
      ! on the library's own internal object-creation bookkeeping
      ! (shared across all threads regardless of how independent the
      ! resulting objects are). The `!$omp critical` below serializes
      ! just that creation step -- each thread still ends up with its
      ! own, self-created, properly-owned object, exactly as before;
      ! only the brief window where 6 threads used to call ast_read_/
      ! compose_pix2pix at the exact same instant is now one-at-a-time.
      ! The actual per-plane resampling (ast_resampler, in the block
      ! loop below) is untouched and stays parallel -- it was
      ! never the problem (band 1 of the real run completed its own
      ! ~44GB output correctly before band 2's setup step deadlocked).
      iblock = 0
      t_status = 0
      !$omp critical (ast_setup)
      call ast_begin(t_status)
      if (used_target_grid) then
         ! No reffile to reopen -- build_target_skymap rebuilds the
         ! final Mapping from the plain values alone, entirely within
         ! this thread's own AST context, so it is as safe here as the
         ! file-based load_wcs/extract_sky_mapping pair below (see T20
         ! above: every AST object this creates is owned by THIS
         ! thread, never handed in from another one).
         call build_target_skymap(target_system_arg, target_projection_arg,&
         &target_center_lon_arg,&
         &target_center_lat_arg, target_cdelt1_arg, target_cdelt2_arg,&
         &target_final_crpix_arg, t_skymap_ref, t_skyframe_ref,&
         &t_pixaxes_ref, t_status)
      else
         call load_wcs(reffile, t_wcs_ref, t_naxes_ref, t_status)
         call extract_sky_mapping(t_wcs_ref, t_skymap_ref, t_skyframe_ref,&
         &t_pixaxes_ref, t_status)
      endif
      call load_wcs(infile, t_wcs_in, t_naxes_in, t_status)
      call extract_sky_mapping(t_wcs_in, t_skymap_in, t_skyframe_in,&
      &t_pixaxes_in, t_status)
      call compose_pix2pix(t_skymap_in, t_skyframe_in, t_skymap_ref,&
      &t_skyframe_ref, t_map_in2ref, t_status)
      if (is_exact .and. t_status.eq.0) then
         ! input pixel -> output sky, without passing through output pixel
         ! coordinates (see exact_build_rows)
         t_sky2sky = ast_convert(t_skyframe_in, t_skyframe_ref, ' ', t_status)
         t_map_in2sky = ast_cmpmap(t_skymap_in, t_sky2sky, .true., ' ', t_status)
         t_lonaxis = ast_geti(t_skyframe_ref, 'LonAxis', t_status)
      endif
      !$omp end critical (ast_setup)
      if (t_status.ne.0) then
         !$omp atomic write
         status_par = -1
      endif

      ! interp=exact: build the overlap weights once for this file. Bands
      ! of output rows are shared among the threads; each uses its own AST
      ! objects (above) and writes only its own rows of the shared table.
      if (is_exact) then
         !$omp single
         call exact_init(exact_wts, nx_in, ny_in, exact_lbnd, exact_ubnd)
         call timer_start(t_stage)
         !$omp end single
         !$omp do schedule(dynamic, 1)
         do iband = 1, n_bands
            if (t_status.eq.0 .and. status_par.eq.0) then
               call exact_build_rows(exact_wts, t_map_in2ref, t_skymap_ref,&
               &t_map_in2sky, t_lonaxis, (iband - 1) * exact_band_rows + 1,&
               &min(iband * exact_band_rows, ny_out), t_status)
               if (t_status.ne.0) then
                  !$omp atomic write
                  status_par = -1
               endif
            endif
         enddo
         !$omp end do
         !$omp single
         call timer_stop('exact_weights', t_stage)
         write(*,'(A)') 'interp=exact: overlap weights built, '//&
         &trim(mib_str(exact_weights_bytes(exact_wts)))//' MiB'
         !$omp end single
      endif

      ! Outer "group" loop: one iteration per combination of the SLOWER
      ! non-sky axes (other_axes(2:n_other)), decoded from igroup via the
      ! same mixed-radix scheme used elsewhere in this file. Every thread
      ! computes the same igroup/other_idx/chan_start/chan_len values
      ! redundantly (cheap, deterministic from shared inputs) so they all
      ! hit the same single/do/single sequence together, as OpenMP
      ! worksharing constructs require.
      do igroup = 1, n_groups
         remainder = igroup - 1
         do k = 2, n_other
            radix = naxes_in(other_axes(k))
            other_idx(k) = mod(remainder, radix) + 1
            remainder = remainder / radix
         enddo

         chan_start = 1
         do while (chan_start.le.axis1_extent)
            chan_len = min(block_planes, axis1_extent - chan_start + 1)

            !$omp single
            call timer_start(t_stage)
            call read_one_block(infile, naxes_in, pixaxes_in, other_axes,&
            &other_idx, n_other, chan_start, chan_len, nx_in, ny_in,&
            &block_data_in(:,:,1:chan_len), status_par)
            call timer_stop('block_read', t_stage)
            call timer_start(t_stage)
            !$omp end single

            iblock = iblock + 1
            tid_local = omp_get_thread_num()
            t_thread_start = omp_get_wtime()
            write(thread_msg,'(A,I0,A,I0,A,I0)') 'thread_timing stage=resample event=start tid=',&
            &tid_local,' block=',iblock,' unit_count=',chan_len
            call log_message('debug','tile_thread',trim(thread_msg))
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
                  ! Same 20-argument astResampleR signature/pitfall as
                  ! before (see git history for the "missing params"
                  ! segfault this avoids); in_var/out_var unused, block
                  ! arrays reused as harmless placeholders. interp_scheme/
                  ! interp_params_arr resolved once, outside this parallel
                  ! region (docs/dev/REPROJECT_INTERP_OPTIONS_PLAN.md) --
                  ! ast__linear/params_dummy hardcoded here previously.
                  ! interp_finterp resolved just above (outside this
                  ! parallel region) as this subroutine's own LOCAL
                  ! variable -- never received as a dummy argument (see
                  ! REPROJECT_INTERP_OPTIONS_PLAN.md section 6d: passing
                  ! it in as this subroutine's own dummy argument
                  ! segfaulted inside AST itself; a gfortran code-
                  ! generation issue specific to a procedure(),pointer
                  ! dummy argument of a contains-ed procedure shared
                  ! into that procedure's own parallel region).
                  if (is_exact) then
                     call exact_apply(exact_wts, block_data_in(:,:,local_iplane),&
                     &block_data_out(:,:,local_iplane,cur_slot))
                  else
                  nbad = ast_resampler(t_map_in2ref, 2, lbnd_in, ubnd_in,&
                  &block_data_in(:,:,local_iplane),&
                  &block_data_in(:,:,local_iplane),&
                  &interp_scheme, interp_finterp, interp_params_arr, 0, 0.0d0, 100, badval,&
                  &2, lbnd_o, ubnd_o, lbnd_o, ubnd_o,&
                  &block_data_out(:,:,local_iplane,cur_slot),&
                  &block_data_out(:,:,local_iplane,cur_slot), t_status)
                  endif
                  if (t_status.ne.0) then
                     !$omp atomic write
                     status_par = -1
                  endif
               endif
            enddo
            !$omp end do
            t_thread_elapsed = (omp_get_wtime() - t_thread_start) * 1000.0_dp
            write(thread_msg,'(A,I0,A,I0,A,I0,A,F10.3)') 'thread_timing stage=resample event=done tid=',&
            &tid_local,' block=',iblock,' unit_count=',chan_len,' dur_ms=',t_thread_elapsed
            call log_message('debug','tile_thread',trim(thread_msg))

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
               write_job%file_path = outfile
               write_job%datastart = datastart
               write_job%naxes_in(1:max_axes) = naxes_in(1:max_axes)
               write_job%other_axes(1:max_axes) = other_axes(1:max_axes)
               write_job%other_idx(1:max_axes) = other_idx(1:max_axes)
               write_job%n_other = n_other
               write_job%chan_start = chan_start
               write_job%chan_len = chan_len
               write_job%nwriters_eff = max(1, min(nwriters_eff, chan_len))
               write_job%nx = nx_out
               write_job%ny = ny_out
               write_job%data => block_data_out(:,:,1:chan_len,cur_slot)
               if (io_overlap_l) then
                  call block_write_dispatch_async(write_job, write_thread_id,&
                  &write_dispatched_ok)
                  write_pending = write_dispatched_ok
               else
                  call do_block_write(write_job)
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
   end subroutine write_reprojected_file

   subroutine copy_axis_keywords(src_unit, src_axis, dst_unit, dst_axis,&
   &crpix_shift, status)
      !! Copy CTYPE/CRVAL/CRPIX/CDELT/CUNIT/CROTA for src_axis (in
      !! src_unit's own header) to dst_axis (in dst_unit's header,
      !! already created via FTPHPR). CRPIX is additionally shifted by
      !! -crpix_shift (0 for a straight passthrough; the output grid's
      !! own pixel-1 offset, reference-pixel-numbered, for a cropped/
      !! grown sky axis) -- matches rm_synthesis.f90's own existing
      !! subimage CRPIX-shift formula exactly, generalised to any axis
      !! number via a constructed keyword string ("CRVAL"//axis, etc.)
      !! rather than a literal "1"/"2" suffix. CROTA needs no shift --
      !! cropping/growing the same grid changes CRPIX (a position within
      !! it), never its rotation. (PCi_j/CDi_j, the 2x2-matrix
      !! alternative to CROTA, are handled separately by
      !! copy_sky_rotation_matrix below -- they're indexed by an axis
      !! *pair*, not a single axis, so they don't fit this per-axis
      !! copy loop.)
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
      !! PCi_j/CDi_j: the 2x2-matrix alternative to CROTA (a file uses
      !! one convention or the other, essentially never both -- copying
      !! whichever is present, doing nothing if neither is, is
      !! the correct behaviour either way). Indexed by an axis PAIR
      !! (e.g. PC2_4), so this reads the 4 entries for the sky axis pair
      !! specifically (src_axis1,src_axis2 -- in the REFERENCE file's own
      !! numbering) and writes them at the output's fixed sky positions
      !! (1,2), the same axis-renumbering copy_axis_keywords already does
      !! for CTYPE/CRVAL/etc, just for a 2-axis-indexed keyword instead
      !! of a 1-axis-indexed one. This is intentionally NOT attempted for
      !! the generic (infile, non-sky-axis) header copy elsewhere
      !! (copy_generic_header skips PC/CD there) -- a rotation
      !! entry for a channel or Stokes axis would be unusual, and
      !! blindly relocating it under the input's own arbitrary axis
      !! numbers, unlike this deliberate sky-axis-pair copy, risks
      !! landing it on the wrong axis.
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
      !! Single-entry helper for copy_sky_rotation_matrix (a sibling
      !! internal subroutine, not nested inside it -- Fortran doesn't
      !! allow a second level of CONTAINS). Copies keyword
      !! "<prefix><sa>_<sb>" from unit su to "<prefix><da>_<db>" on unit
      !! du, only if present on su (absent is not an error -- most files
      !! use CROTA or nothing, not PC/CD, so 0-of-8 entries found here is
      !! the common case).
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
   ! this file's own copy, duplicated a second time in match_cubes.f90.

   subroutine copy_wcs_system_keywords(src_unit, dst_unit, status)
      !! EQUINOX/RADESYS describe the sky coordinate SYSTEM (not a
      !! specific axis), copied verbatim if present -- absent from
      !! src_unit is not an error (many FITS files omit them,
      !! relying on standard defaults; nothing to propagate then).
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

   subroutine copy_generic_header(src_unit, dst_unit, status)
      !! Copy every header card from src_unit (always called with this
      !! opened on INFILE, for the "other" -- non-sky -- axes) to
      !! dst_unit EXCEPT: structural keywords FTPHPR already wrote
      !! (SIMPLE/BITPIX/NAXIS/EXTEND/PCOUNT/GCOUNT/END), and any
      !! axis-indexed keyword (NAXISn/CTYPEn/CRVALn/CRPIXn/CDELTn/
      !! CUNITn/CROTAn/PCn_m/CDn_m, for ANY n/m). CTYPE/CRVAL/CRPIX/
      !! CDELT/CUNIT/CROTA are handled per-axis elsewhere
      !! (copy_axis_keywords), PCi_j/CDi_j elsewhere too
      !! (copy_sky_rotation_matrix), and EQUINOX/RADESYS elsewhere too
      !! (copy_wcs_system_keywords) -- but ALL THREE of those, unlike
      !! this generic copy, are only ever called against the REFERENCE
      !! file for the 2 sky axes specifically, with correct axis
      !! renumbering built in. A rotation/WCS-system keyword found here,
      !! on one of INFILE's own non-sky axes, is skipped rather than
      !! blindly relocated under this file's own arbitrary axis
      !! numbering (which the output does not share) -- landing it on
      !! the wrong axis would be worse than leaving it absent, and a
      !! rotation entry on a channel/Stokes axis would be
      !! unusual in practice anyway. Everything else -- BUNIT,
      !! BMAJ/BMIN/BPA, OBJECT, TELESCOP, INSTRUME, DATE-OBS, RESTFRQ,
      !! HISTORY, COMMENT, and anything this project has never heard of
      !! -- is copied verbatim via a raw 80-column header record
      !! (FTGREC/FTPREC), not decoded/re-encoded through a type-specific
      !! FTGKYx/FTPKYx pair, so formatting/precision/comments survive
      !! exactly. This is what fixes losing
      !! units/beam/provenance from the output (previously ONLY the
      !! axis-indexed WCS keywords were ever copied at all).
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
         if (skip_generic_header_key(key)) cycle
         fitsstat = 0
         call FTPREC(dst_unit, card, fitsstat)
      enddo
   end subroutine copy_generic_header

   logical function skip_generic_header_key(key)
      !! True for any keyword copy_generic_header must NOT copy verbatim
      !! -- see that subroutine's own comment for why each category is
      !! excluded.
      character(len=8), intent(in) :: key

      skip_generic_header_key = .true.
      select case (trim(key))
      case ('SIMPLE', 'BITPIX', 'NAXIS', 'EXTEND', 'PCOUNT', 'GCOUNT',&
      &'END', 'EQUINOX', 'RADESYS', 'WCSAXES', 'LONPOLE', 'LATPOLE')
         return
      end select
      if (is_indexed_keyword(key, 'NAXIS')) return
      if (is_indexed_keyword(key, 'CTYPE')) return
      if (is_indexed_keyword(key, 'CRVAL')) return
      if (is_indexed_keyword(key, 'CRPIX')) return
      if (is_indexed_keyword(key, 'CDELT')) return
      if (is_indexed_keyword(key, 'CUNIT')) return
      if (is_indexed_keyword(key, 'CROTA')) return
      if (is_indexed_keyword(key, 'PC')) return
      if (is_indexed_keyword(key, 'CD')) return
      skip_generic_header_key = .false.
   end function skip_generic_header_key

   logical function is_indexed_keyword(key, prefix)
      !! True if key is exactly prefix followed by one or two groups of
      !! digits separated by at most one underscore (e.g. prefix='CD':
      !! matches "CD1_2", "CD10_2"; prefix='NAXIS': matches "NAXIS1",
      !! "NAXIS12"; does not match "CDELT1" -- the "ELT1" tail after
      !! stripping "CD" is not all-digits-and-at-most-one-underscore).
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

   subroutine read_one_block(filename, naxes_in, pixaxes_in, other_axes,&
   &other_idx, n_other, chan_start, chan_len, nx_in, ny_in, block_data_in,&
   &status)
      !! Read chan_len consecutive planes of filename in one CFITSIO
      !! call: sky axes full extent, other_axes(1) spanning
      !! [chan_start, chan_start+chan_len-1], every slower other axis
      !! (other_axes(2:n_other)) fixed at other_idx(2:n_other) -- the
      !! write-side mirror of write_one_block below, same fpixels/lpixels
      !! construction as the old per-plane resample_one_plane just with a
      !! range instead of a single index on other_axes(1). Always called
      !! from a single OpenMP thread (write_reprojected_file's `!$omp
      !! single` region), so a fixed CFITSIO unit number is safe here --
      !! no other thread can be touching CFITSIO at the same moment, by
      !! construction, not by locking convention.
      !!
      !! FTGSVE fills its output array in ascending-axis-number order
      !! among the non-degenerate (extent>1) axes of THIS read: the 2 sky
      !! axes (always non-degenerate) and other_axes(1) (non-degenerate
      !! whenever chan_len>1). The single-plane version this replaced
      !! never had to care which of the 3 was numerically smallest --
      !! every non-sky axis was degenerate then (extent exactly 1), and a
      !! degenerate axis contributes no stride regardless of where it
      !! sits in the axis order. That stopped being true here: whenever
      !! other_axes(1) itself is numerically BEFORE one or both sky axes
      !! (e.g. TEST_NONADJACENT.Q.FITSCUBE: FREQ on axis 1, RA/DEC on
      !! axes 2 and 4 -- FREQ reads fastest, not slowest), the natural
      !! read order stops matching the caller's fixed
      !! (nx_in,ny_in,chan_len) block_data_in layout, and reading
      !! straight into it silently scrambles the data (caught by
      !! comparing block output against the pre-batching serial output
      !! byte-for-byte -- TEST_NONADJACENT differed in 102397 of 102400
      !! elements). The output side never has this problem -- its axis
      !! layout is always canonicalised to sky-first (see
      !! write_one_block) regardless of the input's own numbering.
      !!
      !! Fixed generally: read into a buffer shaped in the ACTUAL natural
      !! order (via pairwise-comparison ranking of the 3 axis numbers),
      !! then copy into the caller's fixed layout explicitly. Common case
      !! (other_axes(1) numerically after both sky axes, i.e. a
      !! conventionally-ordered cube with channels/Stokes after RA/Dec)
      !! reads directly into block_data_in with no extra buffer or copy
      !! -- the permute path only costs anything on the non-conventional
      !! axis orderings that need it.
      use, intrinsic :: ieee_arithmetic
      character(len=*), intent(in) :: filename
      integer, intent(in) :: naxes_in(:), pixaxes_in(2)
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
         if (naxes_in(k).gt.0) naxis = k
      enddo
      fpixels(1:naxis) = 1
      lpixels(1:naxis) = 1
      incs(1:naxis) = 1
      lpixels(pixaxes_in(1)) = nx_in
      lpixels(pixaxes_in(2)) = ny_in
      if (n_other.ge.1) then
         fpixels(other_axes(1)) = chan_start
         lpixels(other_axes(1)) = chan_start + chan_len - 1
      endif
      do k = 2, n_other
         fpixels(other_axes(k)) = other_idx(k)
         lpixels(other_axes(k)) = other_idx(k)
      enddo

      ! Rank the 3 axes that matter (sky1, sky2, block) by ascending raw
      ! axis number -- pairwise comparison, works for any 3 distinct
      ! integers. natural_order (ranks already 1,2,3 in that order) is
      ! the fast path; anything else needs the permute buffer.
      ax_sky1 = pixaxes_in(1)
      ax_sky2 = pixaxes_in(2)
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
         call FTGSVE(unit, group, naxis, naxes_in(1:naxis),&
         &fpixels(1:naxis), lpixels(1:naxis), incs(1:naxis),&
         &badval, block_data_in, anyflg, fitsstat)
      else
         dims(rank_sky1) = nx_in
         dims(rank_sky2) = ny_in
         dims(rank_block) = chan_len
         allocate(natural_buf(dims(1), dims(2), dims(3)))
         call FTGSVE(unit, group, naxis, naxes_in(1:naxis),&
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

   logical function host_is_big_endian_rp() result(is_be)
      !! Verbatim port of rm_synthesis_mod.f90's own host_is_big_endian
      !! (see convolve_cubes.f90's own host_is_big_endian_cv, same design).
      integer(int32) :: probe
      integer(int8) :: bytes(4)

      probe = 1_int32
      bytes = transfer(probe, bytes)
      is_be = (bytes(1) == 0_int8)
   end function host_is_big_endian_rp

   subroutine swap_bytes_r4_inplace_rp(buf, n)
      !! Verbatim port of rm_synthesis_mod.f90's own
      !! swap_bytes_r4_inplace.
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
   end subroutine swap_bytes_r4_inplace_rp

   subroutine write_one_block_raw(file_path, datastart, naxes_in,&
   &other_axes, other_idx, n_other, chan_start, chan_len, nx_out, ny_out,&
   &block_data_out, status)
      !! Raw Fortran stream I/O at a computed byte offset -- CFITSIO is
      !! never touched here at all (T19, docs/dev/MULTI_BAND_
      !! TOMOGRAPHY_PLAN.md -- see write_reprojected_file's own comment
      !! for the full rationale). Verbatim design port of rm_synthesis_
      !! mod.f90's own write_rm_chunk_raw, generalised for this file's
      !! own N-other-axes output layout. Output axis layout (same
      !! convention the old FTPSSE-based write_one_block used): sky
      !! always at output positions 1,2, full extent; other axes at
      !! output positions 3.. (2+k for other_axes(k)) -- other_axes(1)'s
      !! output slot spans the block's own chan_start:chan_start+
      !! chan_len-1 range, every slower other axis fixed at
      !! other_idx(2:n_other). Since axes 1,2 are always written at
      !! full extent and every axis slower than other_axes(1) is held
      !! fixed for the whole block, chan_len consecutive other_axes(1)
      !! values are always exactly nx_out*ny_out contiguous elements
      !! each, back-to-back in FITS's own axis-1-fastest layout --
      !! fixed_offset below is the constant contribution from the
      !! slower (other_idx(2:n_other)) axes, computed once per call;
      !! only the chan_start term varies block to block. Always called
      !! from a single OpenMP thread by construction (see
      !! read_one_block's own comment) -- no locking needed here.
      character(len=*), intent(in) :: file_path
      integer(kind=8), intent(in) :: datastart
      integer, intent(in) :: naxes_in(:), other_axes(:), other_idx(:)
      integer, intent(in) :: n_other, chan_start, chan_len, nx_out, ny_out
      real, intent(in) :: block_data_out(:,:,:)
      integer, intent(inout) :: status

      integer :: u, ios, ip, k
      logical :: need_swap
      integer(kind=8) :: plane_elems, stride, fixed_offset, byte_pos
      real, allocatable :: plane_buf(:,:)

      if (status.ne.0) return

      need_swap = .not. host_is_big_endian_rp()
      plane_elems = int(nx_out,8) * int(ny_out,8)

      stride = plane_elems
      if (n_other.ge.1) stride = stride * int(naxes_in(other_axes(1)),8)
      fixed_offset = 0_8
      do k = 2, n_other
         fixed_offset = fixed_offset + int(other_idx(k)-1,8) * stride
         stride = stride * int(naxes_in(other_axes(k)),8)
      enddo

      open(newunit=u, file=trim(file_path), access='stream',&
      &form='unformatted', status='old', action='write', iostat=ios)
      if (ios.ne.0) then
         write(*,*) 'ERROR: write_one_block_raw: failed to open ', trim(file_path)
         status = -1
         return
      endif

      allocate(plane_buf(nx_out,ny_out))
      do ip = 1, chan_len
         plane_buf = block_data_out(:,:,ip)
         if (need_swap) call swap_bytes_r4_inplace_rp(plane_buf, plane_elems)
         byte_pos = datastart +&
         &(int(chan_start-1+ip-1,8)*plane_elems + fixed_offset)*4_8 + 1_8
         write(u, pos=byte_pos, iostat=ios) plane_buf
         if (ios.ne.0) then
            write(*,*) 'ERROR: write_one_block_raw: write failed (plane ',&
            &chan_start+ip-1, ') for ', trim(file_path)
            status = -1
            exit
         endif
      enddo
      deallocate(plane_buf)
      close(u)
   end subroutine write_one_block_raw

   subroutine do_block_write(job)
      !! T4d-style: verbatim port of convolve_cubes.f90's own
      !! do_block_write -- see write_convolved_file's own comment there
      !! for the full single-writer-at-a-time rationale. Uses its own
      !! LOCAL status (not write_reprojected_file's own status_par --
      !! that variable is touched by the main thread's `!$omp do` region
      !! concurrently with a background write, so sharing it directly
      !! would race; write_failed is the safe, checked-after-join
      !! handoff instead).
      type(block_write_job_t), intent(inout) :: job
      integer :: status_local, wk, wbase, wrem, wlen_k, woff_k, wchan_start_k
      logical :: any_failed

      any_failed = .false.
      if (job%nwriters_eff.le.1) then
         status_local = 0
         call write_one_block_raw(trim(job%file_path), job%datastart,&
         &job%naxes_in(1:max_axes), job%other_axes(1:max_axes),&
         &job%other_idx(1:max_axes), job%n_other, job%chan_start, job%chan_len,&
         &job%nx, job%ny, job%data, status_local)
         if (status_local.ne.0) any_failed = .true.
      else
         ! nwriters_eff>1 (T21, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md):
         ! same even-split scheme as convolve_cubes.f90's/match_cubes.f90's
         ! own copies -- splits along other_axes(1) (the axis varying
         ! within this block), every slower axis already fixed for the
         ! whole job.
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
               call write_one_block_raw(trim(job%file_path), job%datastart,&
               &job%naxes_in(1:max_axes), job%other_axes(1:max_axes),&
               &job%other_idx(1:max_axes), job%n_other, wchan_start_k, wlen_k,&
               &job%nx, job%ny, job%data(:,:,woff_k+1:woff_k+wlen_k),&
               &status_local)
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
      !! has a SkyFrame recoverable via ast_pickaxes(2,4), so
      !! whatever internal structure connects them is not simply "the two
      !! components of one decompose split") is probed with ast_pickaxes
      !! + ast_isaskyframe -- an AST class check, not a guess -- to
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
      !! presents (Dec,RA), not (RA,Dec)), so the caller needs the
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
      !! The CFITSIO unit comes from fitsio_unit_mod's safe_ftopen --
      !! called concurrently by every OpenMP thread during
      !! write_reprojected_file's per-thread setup, so each thread needs
      !! a distinct unit with no collision risk. Found the hard
      !! way, twice:
      !! (1) an earlier version left this unit a bare literal (11) with
      !! no per-thread distinction at all, so every thread calling
      !! load_wcs concurrently collided on the same CFITSIO unit
      !! number, which is what crashed/deadlocked multi-threaded
      !! runs -- not, as first suspected, any inherent AST/CFITSIO
      !! thread-safety limitation (a battery of isolated reproducers
      !! matching the real call pattern -- concurrent per-thread AST
      !! FitsChan/FrameSet/Mapping construction, concurrent astResampleR
      !! on independent per-thread Mappings, concurrent CFITSIO
      !! reads+writes on correctly-separated unit numbers -- ran cleanly,
      !! 16 threads x hundreds of iterations x many repeats, once unit
      !! numbers stopped colliding). A later, hand-picked
      !! 1000+thread_num offset scheme fixed that.
      !! (2) replacing the hand-picked offset with plain FTGIOU-acquired
      !! units (critical-section-guarded, but with FTOPEN/FTCLOS
      !! themselves OUTSIDE that critical section) reintroduced the exact
      !! same class of crash via a different mechanism: FTGIOU/FTFIOU's
      !! own internal bookkeeping shares mutable global state with
      !! FTOPEN/FTCLOS's own open-file table, so two threads racing a
      !! critical-guarded FTGIOU against another thread's unguarded
      !! FTOPEN still corrupts that shared state. safe_ftopen/safe_ftclos
      !! (fitsio_unit_mod.f90) fix this by keeping FTGIOU+FTOPEN (and
      !! FTCLOS+FTFIOU) inside the SAME critical section -- verified with
      !! an isolated 8-thread reproducer that crashed reliably within a
      !! few hundred iterations without this, and ran clean with it.
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
      ! Rewind the channel's internal read pointer before ast_read
      ! consumes it (Card is a 1-based cursor into the card list;
      ! ast_putfits above leaves it sitting past the last card written).
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

   subroutine get_mem_total_kb(mem_total_kb)
      !! Total system RAM in kB, from /proc/meminfo's MemTotal -- same
      !! source and same reasoning as rm_synthesis.f90's own memory
      !! planner (search that file for "Tile planning for memory-
      !! efficient cube processing"): budgeting against TOTAL RAM rather
      !! than instantaneously-available RAM makes the chosen block size
      !! deterministic for a given cube/mem_frac_ram (reproducible across
      !! runs) instead of fluctuating with whatever else the machine is
      !! doing -- with the same caveat rm_synthesis documents: on a busy/
      !! shared node, a large mem_frac_ram can over-commit, since memory
      !! used by other jobs is not subtracted here. 4 GiB fallback if
      !! /proc/meminfo can't be read (e.g. non-Linux), matching
      !! rm_synthesis's own fallback constant exactly.
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

   subroutine print_usage()
      !! Shared by --help/-h and the argument-error path, so the two
      !! can't drift out of sync with each other.
      write(*,'(A)') 'reproject_cubes -- reproject FITS cubes onto a common grid'
      write(*,'(A)') ''
      write(*,'(A)') 'Usage:'
      write(*,'(A)') '  reproject_cubes mode=<intersection|union|reference>'//&
      &' reffile=<reference_file> infiles=<input_file>[,<input_file>...]'
      write(*,'(A)') '  reproject_cubes [target_system=<system>]'//&
      &' [target_projection=<proj>] infiles=<input_file>[,<input_file>...]'//&
      &' [target_center=<center>] [target_cdelt_lon=<deg> target_cdelt_lat=<deg>]'//&
      &' [target_hpx_nside=<n>]'
      write(*,'(A)') '  reproject_cubes --config <cfgfile>'
      write(*,'(A)') '  reproject_cubes --help | -h'
      write(*,'(A)') ''
      write(*,'(A)') 'No positional args -- every key is key=value (no spaces'//&
      &' around ''='', e.g. mode=union), either directly on the command line'
      write(*,'(A)') 'or via --config. If both are given, each CLI key=value'//&
      &' overrides only that same key from --config; unset keys still come'//&
      &' from --config.'
      write(*,'(A)') ''
      write(*,'(A)') 'infiles= is always required. Exactly one of these two forms'//&
      &' is required too: mode=+reffile= (together -- see Modes below), or'
      write(*,'(A)') 'target_system=/target_projection= with no reffile= at all'//&
      &' (see the target_system= section below) -- giving reffile= together'
      write(*,'(A)') 'with any target_* key is a hard error.'
      write(*,'(A)') ''
      write(*,'(A)') 'Modes:'
      write(*,'(A)') '  reference    output grid is the reference file''s own extent'
      write(*,'(A)') '  intersection output grid shrinks to the overlap of all inputs'//&
      &' with the reference'
      write(*,'(A)') '  union        output grid grows to cover all inputs and the reference'
      write(*,'(A)') ''
      write(*,'(A)') 'Config file: key=value text file. infiles= is always'//&
      &' required; mode=+reffile= together, or target_system=/target_projection='
      write(*,'(A)') 'with no reffile= at all (see above):'
      write(*,'(A)') '  mode    = intersection | union | reference'
      write(*,'(A)') '  reffile = /path/to/reference.fits'
      write(*,'(A)') '  infiles = /path/band1.fits,/path/band2.fits,/path/band3.fits'
      write(*,'(A)') '(infiles is comma-separated, no spaces required; ''#'' or '';'''//&
      &' starts a comment.)'
      write(*,'(A)') ''
      write(*,'(A)') 'Optional key (CLI or config, default 0.25):'
      write(*,'(A)') '  mem_frac_ram = fraction (0,0.95] of total system RAM budgeted'//&
      &' for one read/resample/write block of planes at a time -- same concept'
      write(*,'(A)') '  as rm_synthesis''s own mem_frac_ram. Smaller = more, smaller'//&
      &' blocks (less peak memory, more CFITSIO calls); larger = fewer, bigger blocks.'
      write(*,'(A)') '  Threads only parallelise WITHIN one block, never across blocks'//&
      &' (blocks are processed strictly one after another) -- too small a'
      write(*,'(A)') '  mem_frac_ram therefore also throws away most of the OpenMP'//&
      &' speedup, not just increases I/O calls (a printed WARNING flags this).'
      write(*,'(A)') ''
      write(*,'(A)') 'Optional key (CLI or config, default n):'
      write(*,'(A)') '  io_overlap = y|n -- write each block on a background thread,'//&
      &' overlapped with the NEXT block''s own read+resample, instead of blocking'
      write(*,'(A)') '  on the write before starting it. Only one background write is'//&
      &' ever in flight at a time (concurrent CFITSIO handle use from two'
      write(*,'(A)') '  threads is unsafe) -- still overlaps the write with the'//&
      &' following block''s read+compute, which is where the dead time is.'
      write(*,'(A)') ''
      write(*,'(A)') 'Optional key (CLI or config, default 1):'
      write(*,'(A)') '  nwriters = split one block''s own planes into this many disjoint'//&
      &' concurrent writers -- same key/clamp as rm_synthesis'' and'
      write(*,'(A)') '  rmclean_cubes'' own nwriters (max(1, min(nwriters, OMP thread'//&
      &' count))). Orthogonal to io_overlap: io_overlap decides WHEN a'
      write(*,'(A)') '  block''s write runs; nwriters decides how many concurrent'//&
      &' writers do it once dispatched.'
      write(*,'(A)') ''
      write(*,'(A)') 'Optional key (CLI or config, default n):'
      write(*,'(A)') '  dry_run = y|n -- check the target output disk''s own rotational'//&
      &' status and write a suggested reproject_cubes_dryrun.cfg'
      write(*,'(A)') '  (io_overlap/nwriters), instead of processing any real file.'//&
      &' Advisory only. Touches no data. Also reports the disk-space check below.'
      write(*,'(A)') ''
      write(*,'(A)') 'Optional key (CLI or config, default 10):'
      write(*,'(A)') '  min_free_pct_after_write = percent [0,100] of the output disk''s'//&
      &' TOTAL size that must still be free once every output file is'
      write(*,'(A)') '  written. Checked before any file is created; the run refuses,'//&
      &' printing required/free/remaining/minimum, if it would leave less.'
      write(*,'(A)') ''
      write(*,'(A)') 'Optional keys (CLI or config, logging/timing):'
      write(*,'(A)') '  log_level       = error|warn|info|debug (default info)'
      write(*,'(A)') '  timing_enabled  = y|n -- print a stage timing summary (default n)'
      write(*,'(A)') '  log_output_file = path -- append log/timing output to this file'//&
      &' instead of stdout (default empty = stdout)'
      write(*,'(A)') ''
      write(*,'(A)') 'Optional key (CLI or config, default lanczos3):'
      write(*,'(A)') '  interp = nearest|linear|sinc|sincsinc|sinccos|sincgauss|blockave|'//&
      &'gauss|somb|sombcos|lanczos2|lanczos3|lanczos|exact -- interpolation kernel'
      write(*,'(A)') '  for resampling. lanczos3 is the default since 2026-10-13'//&
      &' (docs/dev/BRANCH_STATUS_update-reproject.md UR-12: mean RMS error ~6x'
      write(*,'(A)') '  lower than linear on realistic sources, measured against'//&
      &' analytic ground truth and cross-checked against astropy''s own'
      write(*,'(A)') '  reproject package); linear is AST__LINEAR (this tool''s only'//&
      &' behaviour before interp= existed); the rest are AST''s own named'
      write(*,'(A)') '  schemes. gauss/sincgauss/blockave take a kernel width AST'//&
      &' cannot choose itself; omitting interp_params= uses a measured'
      write(*,'(A)') '  default (UR-19: gauss FWHM 0.8 px, sincgauss FWHM 8.0 px,'//&
      &' blockave half-width 1.0 px -- lowest worst-case error over beam'
      write(*,'(A)') '  samplings and output geometries, tests/'//&
      &'check_interp_width_defaults.py) and the value used is reported.'
      write(*,'(A)') '  The widths are in pixels and do not scale with the beam.'//&
      &' sincsinc carries a known AST bug (confirmed'
      write(*,'(A)') '  against AST''s own source: it truncates before its'//&
      &' own window function reaches zero, unlike sinccos); sinccos is'
      write(*,'(A)') '  the clean alternative with similar character.'
      write(*,'(A)') ''
      write(*,'(A)') 'Optional key (CLI or config, scheme-specific default):'
      write(*,'(A)') '  interp_params = comma-separated list of up to 2 numbers:'//&
      &' neighbour/support-radius count (<=0: let AST choose), then the'
      write(*,'(A)') '  scheme''s own shape value where it has one (sinc width for'//&
      &' sincsinc/sinccos/sombcos, Gaussian FWHM for sincgauss/gauss).'
      write(*,'(A)') '  sinc/somb/sincsinc/sinccos/sombcos use the first number as'//&
      &' their own neighbour/support-radius count; blockave uses it as its own'
      write(*,'(A)') '  box half-width (see the default-width note above). Unused'//&
      &' only for nearest/linear. See docs/dev/REPROJECT_INTERP_OPTIONS_PLAN.md'
      write(*,'(A)') '  section 5 for the full table.'
      write(*,'(A)') ''
      write(*,'(A)') 'interp also accepts lanczos2|lanczos3|lanczos -- Lanczos-a'//&
      &' interpolation (a corrected alternative to sincsinc/sinccos above:'
      write(*,'(A)') '  continuously tapers to exactly zero at its own support edge,'//&
      &' rather than cutting off partway like sincsinc does). lanczos2/lanczos3'
      write(*,'(A)') '  fix a=2/a=3 with support locked at the correct radius for'//&
      &' that a (not configurable -- interp_params is rejected for these two);'
      write(*,'(A)') '  lanczos is the general form for any other a.'
      write(*,'(A)') ''
      write(*,'(A)') 'interp also accepts exact -- flux-conserving reprojection, not a'//&
      &' kernel: each output pixel is the'
      write(*,'(A)') '  overlap-area-weighted mean of the input pixels it covers on the'//&
      &' sky (every pixel a spherical'
      write(*,'(A)') '  quadrilateral with great-circle edges, the method of astropy''s'//&
      &' reproject_exact). The sum of value'
      write(*,'(A)') '  x pixel solid angle -- a source''s flux in Jy/beam data -- is'//&
      &' carried over wherever output and'
      write(*,'(A)') '  input pixel sizes differ, which the kernels above do only'//&
      &' approximately. Averaging over each'
      write(*,'(A)') '  output pixel''s area adds that pixel to the effective beam.'//&
      &' interp_params is rejected. The overlap'
      write(*,'(A)') '  weights are built once per input file; their memory comes out of'//&
      &' the mem_frac_ram budget, and the'
      write(*,'(A)') '  run prints the estimate and the size built.'
      write(*,'(A)') ''
      write(*,'(A)') 'Optional key (CLI or config, lanczos only):'
      write(*,'(A)') '  interp_params = <support>,<a> -- required for interp=lanczos'//&
      &' (no default a exists); support<=0 uses the correct default, ceil(a).'
      write(*,'(A)') '  A smaller support is accepted as a deliberate speed/accuracy'//&
      &' trade-off (same kind of knob sincsinc/sinccos/sincgauss already allow),'
      write(*,'(A)') '  but is never auto-derived for lanczos the way it is for'//&
      &' those -- give it explicitly. Rejected for lanczos2/lanczos3'
      write(*,'(A)') '  (their support is fixed; use interp=lanczos for an'//&
      &' explicit a/support combination instead).'
      write(*,'(A)') ''
      write(*,'(A)') 'reffile= is optional. When given, it already fully'//&
      &' determines the output grid -- combining it with ANY of the'
      write(*,'(A)') '  target_system=/target_projection=/target_center=/'//&
      &'target_lon=/target_lat=/target_cdelt_lon=/target_cdelt_lat=/'//&
      &'target_hpx_nside= keys below'
      write(*,'(A)') '  is a hard error, not a silent ignore. mode= is'//&
      &' required together with reffile=, and has no meaning without one.'
      write(*,'(A)') ''
      write(*,'(A)') 'When reffile= is absent, these keys build a fresh output'//&
      &' grid from scratch, fit to the infiles'' own combined footprint'
      write(*,'(A)') '  (the same footprint-bounding machinery mode=union'//&
      &' already uses) -- at least one of target_system=/target_projection='
      write(*,'(A)') '  is then required (nothing requested otherwise). mode='//&
      &' is not used in this case.'
      write(*,'(A)') '  target_system = icrs|fk5|fk4|galactic -- when omitted,'//&
      &' the infiles'' own system is kept (they must all share one).'
      write(*,'(A)') '  target_projection = tan|sin|car|ait|hpx (any case) --'//&
      &' sin when omitted. tan/sin suit a pointed field; car/ait/hpx a'
      write(*,'(A)') '    wide or all-sky one. sin samples the sky more coarsely'//&
      &' away from the centre (by 1/cos of the angle from it) and cannot'
      write(*,'(A)') '    reach past 90 deg from it.'
      write(*,'(A)') '  target_center = infiles_centroid|user_defined --'//&
      &' defaults to infiles_centroid when omitted.'
      write(*,'(A)') '    infiles_centroid: the new field centre is the'//&
      &' centroid of every infile''s own pixel-grid-centre sky position'
      write(*,'(A)') '    (vector-averaged, not naive lon/lat averaging).'
      write(*,'(A)') '    user_defined: the new field centre is given directly'//&
      &' via target_lon=/target_lat= (required together, rejected otherwise).'
      write(*,'(A)') '  target_lon/target_lat = coordinate value, already'//&
      &' expressed in target_system -- accepts the same flexible notation as'
      write(*,'(A)') '    scripts/casda_fetch.py''s own --lon/--lat: plain decimal'//&
      &' degrees, colon-sexagesimal (12:34:56.7), or letter-suffixed'
      write(*,'(A)') '    sexagesimal (12h34m56.7s, -45d30m00s, 180.5d). A bare'//&
      &' colon-sexagesimal target_lon is read as hour angle for icrs/fk5/fk4'
      write(*,'(A)') '    (the standard RA convention), as degrees for galactic'//&
      &' (no hour-angle convention of its own); target_lat is always degrees.'
      write(*,'(A)') '  target_cdelt_lon/target_cdelt_lat = pixel scale (degrees'//&
      &' per pixel) for the new grid''s own 2 axes. Optional, and must be'
      write(*,'(A)') '    given together: the default is the smallest |CDELT| across'//&
      &' every infile, per axis, so no input is coarsened at the grid'
      write(*,'(A)') '    centre (the projection changes the sky per pixel'//&
      &' away from it; see target_projection). Where the'
      write(*,'(A)') '    infiles'' own scales differ, the chosen value and the'//&
      &' spread are both printed.'
      write(*,'(A)') '  target_hpx_nside = power of 2. Only valid with'//&
      &' target_projection=hpx -- a HEALPix-aligned grid''s centre, pixel'
      write(*,'(A)') '    scale and reference pixel all follow from Nside alone'//&
      &' (target_center=/target_lon=/target_lat=/target_cdelt_lon=/'
      write(*,'(A)') '    target_cdelt_lat= are then errors). Optional: default'//&
      &' Nside is the one closest to the finest infile pixel, coarser only'
      write(*,'(A)') '    if every infile has a beam spanning >= 2.5 of that'//&
      &' Nside''s cell, never coarser for a no-beam infile.'
      write(*,'(A)') ''
      write(*,'(A)') 'Beam sampling: onto output cells larger than an infile''s own,'//&
      &' its beam minor axis (BMIN, header or narrowest BEAMS-table channel) must'
      write(*,'(A)') '  span at least 2.0 output cells, or the run is refused before'//&
      &' anything is written; 2.0 to 2.5 cells gives a barely-Nyquist WARNING.'
      write(*,'(A)') '  An infile with no beam information is refused onto coarser'//&
      &' cells. No beam is ever widened automatically: convolve first'
      write(*,'(A)') '  (match_cubes stages=both, or convolve_cubes); the refusal'//&
      &' prints the commands to use.'
      write(*,'(A)') ''
      write(*,'(A)') 'Advisory: for rm_synthesis''s own multi-band use case'//&
      &' (via match_cubes), always supply reffile= -- mode=union/intersection'
      write(*,'(A)') '  need a stable reference grid, which the no-reffile'//&
      &' path does not provide (it has no mode= concept at all).'
   end subroutine print_usage

   subroutine read_reproject_cfg(cfgfile, mode, reffile, infiles, n_inputs,&
   &mem_frac_ram, io_overlap_l, nwriters_l, dry_run_l, min_free_pct_after_write_l,&
   &log_level_l, timing_enabled_l,&
   &log_output_file_l, interp_l, interp_params_l, target_system_l, target_projection_l,&
   &target_center_l, target_lon_l, target_lat_l, target_cdelt_lon_l,&
   &target_cdelt_lat_l, target_hpx_nside_l, cfg_seen_mode, cfg_seen_reffile,&
   &cfg_seen_infiles, status)
      !! Parse a --config key=value file: infiles always required;
      !! mode/reffile required together or both omitted (UR-7 --
      !! reffile= is now optional); comma-separated infiles, same
      !! csv-list convention rm_synthesis's own multi-band keys use;
      !! plus optional keys mem_frac_ram/io_overlap/min_free_pct_after_
      !! write/log_level/timing_enabled/log_output_file/interp/interp_params/
      !! target_system/target_projection/target_center/target_lon/
      !! target_lat/target_cdelt_lon/target_cdelt_lat (all intent
      !! (inout) -- left untouched, keeping whatever default the
      !! caller set, if the file doesn't mention them). cfg_seen_mode/
      !! cfg_seen_reffile/cfg_seen_infiles (intent(out)) let the caller
      !! correctly merge CLI-vs-cfg precedence for the now-optional
      !! mode/reffile pair, the same way it already tracks its own
      !! cli_seen_mode/cli_seen_reffile.
      !! Standalone re-implementation of rm_synthesis_mod's
      !! split_key_value/csv_count/csv_get_item (below) rather than a
      !! `use` dependency -- this tool is deliberately kept off the main
      !! rm_synthesis build graph (own binary, own dependency set, see
      !! the Makefile comment), and these are a handful of generic
      !! string-parsing lines each, unlikely to drift.
      character(len=*), intent(in) :: cfgfile
      character(len=*), intent(out) :: mode, reffile
      character(len=*), intent(out) :: infiles(:)
      integer, intent(out) :: n_inputs
      real, intent(inout) :: mem_frac_ram
      logical, intent(inout) :: io_overlap_l
      integer, intent(inout) :: nwriters_l
      logical, intent(inout) :: dry_run_l
      double precision, intent(inout) :: min_free_pct_after_write_l
      character(len=*), intent(inout) :: log_level_l
      logical, intent(inout) :: timing_enabled_l
      character(len=*), intent(inout) :: log_output_file_l
      character(len=*), intent(inout) :: interp_l, interp_params_l
      character(len=*), intent(inout) :: target_system_l, target_projection_l
      character(len=*), intent(inout) :: target_center_l
      character(len=*), intent(inout) :: target_lon_l, target_lat_l
      character(len=*), intent(inout) :: target_cdelt_lon_l, target_cdelt_lat_l
      character(len=*), intent(inout) :: target_hpx_nside_l
      logical, intent(out) :: cfg_seen_mode, cfg_seen_reffile, cfg_seen_infiles
      integer, intent(out) :: status

      character(len=16384) :: line, val, raw_infiles
      character(len=512) :: key
      integer :: unit_cfg, ios, line_no, j, ios_mfr, ios_nw, ios_mfp
      logical :: has_kv

      status = 0
      n_inputs = 0
      cfg_seen_mode = .false.
      cfg_seen_reffile = .false.
      cfg_seen_infiles = .false.
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
            cfg_seen_mode = .true.
         case ('reffile')
            reffile = trim(val)
            cfg_seen_reffile = .true.
         case ('infiles')
            raw_infiles = trim(val)
            cfg_seen_infiles = .true.
         case ('mem_frac_ram')
            read(val, *, iostat=ios_mfr) mem_frac_ram
            if (ios_mfr.ne.0) then
               write(*,*) 'ERROR: mem_frac_ram must be a number, at line ',&
               &line_no, ' in ', trim(cfgfile)
               status = -1
               close(unit_cfg)
               return
            endif
         case ('io_overlap')
            io_overlap_l = flag_from_value_reproject(val)
         case ('nwriters')
            read(val, *, iostat=ios_nw) nwriters_l
            if (ios_nw.ne.0 .or. nwriters_l.lt.1) then
               write(*,*) 'ERROR: nwriters must be an integer >= 1, at line ',&
               &line_no, ' in ', trim(cfgfile)
               status = -1
               close(unit_cfg)
               return
            endif
         case ('dry_run')
            dry_run_l = flag_from_value_reproject(val)
         case ('min_free_pct_after_write')
            read(val, *, iostat=ios_mfp) min_free_pct_after_write_l
            if (ios_mfp.ne.0) then
               write(*,*) 'ERROR: min_free_pct_after_write must be a number, at line ',&
               &line_no, ' in ', trim(cfgfile)
               status = -1
               close(unit_cfg)
               return
            endif
         case ('log_level')
            log_level_l = trim(val)
         case ('timing_enabled')
            timing_enabled_l = flag_from_value_logging(val)
         case ('log_output_file')
            log_output_file_l = trim(val)
         case ('interp')
            interp_l = trim(val)
         case ('interp_params')
            interp_params_l = trim(val)
         case ('target_system')
            target_system_l = trim(val)
         case ('target_projection')
            target_projection_l = trim(val)
         case ('target_center')
            target_center_l = trim(val)
         case ('target_lon')
            target_lon_l = trim(val)
         case ('target_lat')
            target_lat_l = trim(val)
         case ('target_cdelt_lon')
            target_cdelt_lon_l = trim(val)
         case ('target_cdelt_lat')
            target_cdelt_lat_l = trim(val)
         case ('target_hpx_nside')
            target_hpx_nside_l = trim(val)
         case default
            write(*,*) 'ERROR: unrecognised config key "', trim(key), '" at line ',&
            &line_no, ' in ', trim(cfgfile)
            status = -1
            close(unit_cfg)
            return
         end select
      enddo
      close(unit_cfg)

      if (.not. cfg_seen_infiles) then
         write(*,*) 'ERROR: config file ', trim(cfgfile), ' must set infiles'
         status = -1
         return
      endif
      if (cfg_seen_mode .neqv. cfg_seen_reffile) then
         write(*,*) 'ERROR: config file ', trim(cfgfile),&
         &' sets only one of mode/reffile -- give both (reffile= fully'//&
         &' determines the output grid) or neither (target_system=/'//&
         &'target_projection= build one instead)'
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

   function flag_from_value_reproject(val) result(flag)
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
   end function flag_from_value_reproject

   function strip_fits_ext(filename) result(base)
      !! Output-name helper: verbatim port of convolve_cubes.f90's own --
      !! see there for the full rationale (avoids the double-extension
      !! "name.fits_REPROJ.FITS" this used to produce). Strips whatever
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

   subroutine advise_beam_consistency(infiles, n_inputs)
      !! UR-17 (docs/dev/BRANCH_STATUS_update-reproject.md): advisory
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
      !! rotational, resolved from target_path's own mount point) and
      !! writes a suggested <tool_name>_dryrun.cfg with recommended
      !! io_overlap/nwriters settings. Advisory only: this project's own
      !! nwriters clamp formula was deliberately left tied to
      !! omp_get_max_threads(), not disk type (see T7/T12's own
      !! discussion) -- this dry-run pass suggests a starting VALUE within
      !! that clamp, it does not change the clamp itself. See match_cubes.
      !! f90's own copy of this subroutine for the full rationale (adapt,
      !! don't share, per this project's own module convention).
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

end program reproject_cubes
