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
! defaulting to 'mid' (matching rmclean_mod's own get_lsq_ref_compute
! doc comment: "the RECOMMENDED DEFAULT... pick it unless you have a
! specific reason not to" -- this program's own default previously said
! 'native' here, a real drift from that recommendation, not a
! deliberate departure; fixed together with this comment). Choosing a
! DIFFERENT value is fully supported: derotate_to_lsq_ref is an exact,
! lossless phase rotation of the ALREADY-SAMPLED dirty spectrum (verified
! |P|-invariant to 1e-4 relative, tests/test_rmclean_lsqref_flex.f90), so
! re-expressing it at any other reference before building the table costs
! nothing in accuracy.
!
! COMPUTE cost, corrected: build_rmsf_offset_table/clean_complex/
! restore_clean's own FAST-PATH cost (the common case: refine_peak_
! matched_filter's fixed-location closed-form fit, no search) depends
! only on nrm/niter/table_oversample, independent of lsq_ref_compute --
! matching this comment's own earlier claim, which was correct as far as
! it went. But refine_peak_matched_filter's ESCALATION path (full
! search, whenever the fast fit's own leftover misfit exceeds threshold)
! is NOT independent of lsq_ref_compute: its own m_search/
! cycles_in_window are set by max_offset=max_k|l_sq(k)-lsq_ref_compute|
! (src/rmclean.f90's own refine_peak_matched_filter doc comment), the
! exact same quantity get_drm minimizes at mode=mid -- so a favourable
! lsq_ref_compute reduces the PER-ESCALATION search cost here too, not
! only the upstream grid-size choice at rm_synthesis time. Previously
! documented as "costs nothing in COMPUTE... none of which depend on
! lsq_ref_compute", which is only true for the fast path, not the
! escalation search -- corrected here since real (non-synthetic,
! possibly noisy) data escalates often enough for this to matter, unlike
! this project's own clean synthetic test scenarios. The grid SIZE
! (nrm/CDELT3) itself remains fixed by whatever rm_synthesis already
! wrote (its own lsq_ref_mode option, added alongside this one) --
! lsq_ref_compute chosen here cannot change that, only the per-pixel
! escalation cost against it.
!
! Apply it correctly regardless of mode: build the RMSF table at
! whatever lsq_ref_compute ends up in force, matching the (possibly-
! derotated) dirty spectrum exactly (build_rmsf_offset_table's own
! correctness requirement: R(delta) is only a pure function of offset
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
!    [niter=<n>] [gain=<g>] [abs_flux_floor=<v>] [auto_nsigma=<n>]
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
   use, intrinsic :: iso_c_binding, only: c_int, c_long, c_ptr, c_funptr,&
   &c_null_ptr, c_funloc, c_loc, c_f_pointer
   use rmclean_mod
   use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_wtime,&
   &omp_set_num_threads
   use logging_mod
   use fitsio_unit_mod
   use rmclean_io_mod
   use rmclean_cache_mod, only: fnv1a_hash, linear_scan_extreme,&
   &pattern_registry_t, registry_init, registry_lookup_or_insert,&
   &registry_lookup, registry_advance, registry_next_occurrence
   implicit none

   ! --- Logging & timing (planning-doc ticket) -- see convolve_cubes.f90
   ! /logging_mod.f90's own header comments. Stage names: startup,
   ! tile_read, tile_compute, tile_write_join, tile_write, finalize.
   character(len=16) :: log_level
   logical :: timing_enabled
   character(len=272) :: log_output_file

   character(len=512) :: ampfile, phafile, maskfile, outfile
   integer :: niter
   real(sp) :: gain
   ! --- CLEAN stopping criteria (planning/RMCLEAN_INTEGRATION_PLAN.md
   ! T8 -- three independent, unambiguously-named conditions; whichever
   ! fires first stops CLEAN for that pixel; niter (above) is always the
   ! hard backstop, the other two are each independently optional and
   ! MAY be combined (not mutually exclusive -- a real design change
   ! from the old threshold=/threshold_snr= pair, which required
   ! exactly one). See clean_complex's own top comment (src/rmclean.f90)
   ! for the full story on why these were split out of one overloaded,
   ! always-multiplicative `thresh` parameter that silently broke the
   ! absolute-flux mode.
   ! abs_flux_floor=<v>: stop the instant a pixel's own peak amplitude
   ! drops to/below this literal fixed flux value -- accepts a bare
   ! number (native AMP-cube units) or a Jy/mJy/uJy suffix, converted
   ! against the AMP cube's own BUNIT once it's read (resolve_abs_flux_
   ! floor, called after read_amp_pha_geometry).
   logical :: have_abs_flux_floor
   real(sp) :: abs_flux_floor, abs_flux_floor_raw_value, abs_flux_floor_unit_scale
   ! auto_nsigma=<n>: stop when a pixel's own peak amplitude drops
   ! to/below n x that SAME pixel's own noise sigma. Sigma is estimated
   ! once per pixel from its own FULL dirty spectrum's own interquartile
   ! range (clean_complex/estimate_iqr_sigma, src/rmclean.f90 -- see its
   ! own comment for why the full spectrum, not a lowest-percentile
   ! subset, turned out to be the more accurate estimator) -- not a
   ! whole-cube pre-scan, and not recomputed per iteration, unlike the
   ! old threshold_snr's noise_floor + thresh*rms_val combination.
   logical :: have_auto_nsigma
   real(sp) :: auto_nsigma_mult
   real(sp) :: min_samples_per_fwhm, refine_nsigma
   ! trace_ix/trace_iy (default 0, disabled): 1-indexed GLOBAL pixel to
   ! log a full per-iteration CLEAN trend for (peak_val/rms_val/stop
   ! reason) -- a debugging/tuning aid (planning/RMCLEAN_INTEGRATION_
   ! PLAN.md, T6-adjacent), not meant for whole-cube production runs.
   integer :: trace_ix, trace_iy
   ! log_every (default 50): throttles trace_ix/trace_iy's own per-
   ! iteration log lines to every log_every-th iteration -- iteration 1
   ! (the pre-CLEAN dirty-spectrum state, see clean_complex's own
   ! trace_flux_val comment) and the final iteration are ALWAYS logged
   ! regardless of log_every, so the trend's start/end points are never
   ! missed even if log_every doesn't divide n_iter_used.
   integer :: log_every
   integer :: table_oversample
   logical :: have_restore_fwhm
   real(sp) :: restore_fwhm_override
   integer :: lsq_ref_report_mode_sel
   logical :: have_lsq_ref_report_value
   real(sp) :: lsq_ref_report_value
   character(len=16) :: lsq_ref_compute_mode
   logical :: have_lsq_ref_compute_value
   real(sp) :: lsq_ref_compute_value

   ! --- T4a: tile geometry (planning/RMCLEAN_INTEGRATION_PLAN.md T4) ---
   ! Same scheme, same key names/defaults, as rm_synthesis's own plan_tile
   ! (rm_synthesis_mod.f90) -- see plan_rmclean_tile's own comment for why
   ! (RA-strips-first, mem_frac_ram-budgeted, safety-shrink). tile_ra/
   ! tile_dec here double as BOTH the config-requested override (0 = auto,
   ! read by parse_args) AND, after plan_rmclean_tile runs, the actual
   ! planned tile size every subsequent tile in the main loop allocates
   ! against (the max size; the last tile in each direction may use less).
   integer :: tile_ra, tile_dec
   logical :: tile_auto
   real(sp) :: mem_frac_ram

   ! --- T17 (docs/dev/RMCLEAN_INTEGRATION_PLAN.md): coverage-fraction
   ! pixel skip -- pixels whose own valid-channel fraction (summed
   ! across ALL bands) falls below this threshold are never CLEANed at
   ! all (output NaN across every RM-bin), never registered in Pass
   ! 0's own pattern registry, and never inserted into the runtime
   ! mask-pattern cache -- a sparse-coverage pixel's own RMSF is
   ! already coarser/noisier than a well-covered one regardless of how
   ! much compute is spent on it. Same "_frac" fraction-value naming
   ! convention as mem_frac_ram; default 0.0 (off -- every pixel with
   ! at least one valid channel is CLEANed, today's original
   ! behaviour), matching abs_flux_floor's own "0.0 = inert" precedent
   ! elsewhere in this file.
   real(sp) :: min_valid_chan_frac

   ! --- T4b: io_read_threads (planning/RMCLEAN_INTEGRATION_PLAN.md T4) ---
   ! Same scheme/key name/clamping convention as rm_synthesis's own
   ! io_read_threads: each thread opens its OWN readonly CFITSIO handle
   ! (safe -- readonly opens are exempt from the same-file handle-aliasing
   ! hazard read-write handles hit, see nwriters/T4c's own
   ! comment) and reads its own disjoint RM-slice of the tile.
   integer :: io_read_threads, io_read_threads_eff

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

   ! T10b: mask_tile is now tiled exactly like re_tile/im_tile below
   ! (replacing the old whole-cube-resident mask_cube + its own separate
   ! read_mask_cube, T4's original deliberate divergence) -- read fresh
   ! per tile via read_mask_tile, and the mask-pattern cache is now built
   ! INCREMENTALLY, one tile's own (1:tx,1:ty,:) subrange at a time
   ! (update_mask_pattern_cache_for_tile), rather than one upfront global
   ! pre-scan requiring the whole mask resident in memory. This also
   ! makes the mask's own memory footprint part of the SAME
   ! mem_frac_ram-budgeted tile-sizing formula as everything else
   ! (plan_rmclean_tile's own bytes_per_tile_pixel) -- a big machine
   ! naturally gets a tile large enough to recover the old whole-cube-
   ! resident behaviour for free (up to tx=nx,ty=ny, if it fits the
   ! budget); a small machine naturally gets smaller tiles instead of
   ! failing outright.
   integer(kind=1), allocatable :: mask_tile(:,:,:)

   ! T4a: per-tile input buffers, allocated ONCE at the planned
   ! (tile_ra,tile_dec,nrm) max size and reused across every tile -- a
   ! tile at the image's own right/bottom edge just uses a (tx,ty,:)
   ! subrange of this same storage (tx<=tile_ra, ty<=tile_dec), same
   ! convention reproject_cubes.f90's own block_data_in/out use for a
   ! short final block. Never double-buffered (T4d's own comment on the
   ! write_job/_buf arrays below): once a tile's compute step has consumed
   ! re_tile/im_tile, nothing later in that tile's own pipeline needs
   ! them again, so the next tile's read can safely reuse the same
   ! storage immediately -- unlike the OUTPUT arrays below, which a
   ! still-in-flight background write may still be reading.
   real(sp), allocatable :: re_tile(:,:,:), im_tile(:,:,:)

   ! T4d: output-tile storage is double-buffered (io_overlap) -- 2
   ! physical copies (trailing dimension, 1-based slot 1/2) of each of
   ! the 6 output arrays, with clean_re_tile etc as POINTERs re-targeted
   ! at slot cur_slot+1 every tile (same _s0/_s1 ping-pong concept as
   ! rm_synthesis.f90's own p_tile_arr/phi_tile_arr, generalized here to
   ! a trailing array dimension instead of named pairs, and to 6 outputs
   ! instead of rm_synthesis's 2). With io_overlap=n, cur_slot is always
   ! 0 and this costs one extra (unused) physical copy of each array --
   ! negligible next to a whole tile's own float-cube footprint, and
   ! avoids a separate no-overlap code path entirely.
   real(sp), allocatable, target :: clean_re_buf(:,:,:,:), clean_im_buf(:,:,:,:)
   real(sp), allocatable, target :: resid_re_buf(:,:,:,:), resid_im_buf(:,:,:,:)
   real(sp), allocatable, target :: restored_re_buf(:,:,:,:), restored_im_buf(:,:,:,:)
   real(sp), pointer :: clean_re_tile(:,:,:) => null(), clean_im_tile(:,:,:) => null()
   real(sp), pointer :: resid_re_tile(:,:,:) => null(), resid_im_tile(:,:,:) => null()
   real(sp), pointer :: restored_re_tile(:,:,:) => null(), restored_im_tile(:,:,:) => null()

   integer :: ix_tile_beg, iy_tile_beg, tx, ty

   ! --- T4d: io_overlap (planning/RMCLEAN_INTEGRATION_PLAN.md T4) ---
   ! Same scheme as rm_synthesis's own io_overlap: a raw POSIX thread (not
   ! an OpenMP task -- see tile_write_job_t's own comment for why) runs
   ! one tile's write while the NEXT tile's read/compute proceeds using
   ! the other buffer slot.
   logical :: io_overlap
   integer :: tile_seq, cur_slot
   integer :: n_blocks_total
   integer(c_long) :: write_thread_id(0:1)
   logical :: write_pending(0:1) = .false.
   logical :: write_dispatched_ok

   ! One pointer per output array a job needs to write, wrapped in a
   ! derived type since Fortran has no array-of-pointers primitive.
   type :: tile_ptr_t
      real(sp), pointer :: p(:,:,:) => null()
   end type tile_ptr_t

   ! Asynchronous tile-write job -- see rm_synthesis_mod.f90's own
   ! tile_write_job_t for the full "why a pthread, not an omp task"
   ! rationale (identical here: keeping read/compute's own `!$omp
   ! parallel do` regions completely undisturbed by the write's own
   ! lifetime, which must be able to outlive them). re_ptr/im_ptr(k),
   ! k=1..3, address (clean,resid,restored) in that fixed order --
   ! do_tile_write pairs each with out_unit/out_path/out_datastart
   ! (idx_clean_amp,idx_clean_pha)/(idx_resid_amp,idx_resid_pha)/
   ! (idx_restored_amp,idx_restored_pha) by the same fixed order.
   type :: tile_write_job_t
      integer :: ix_tile_beg = 0, iy_tile_beg = 0, tx = 0, ty = 0
      type(tile_ptr_t) :: re_ptr(3), im_ptr(3)
   end type tile_write_job_t
   ! target+save: must remain valid (untouched, undeallocated) from
   ! dispatch until tile_write_join returns for that slot -- see
   ! tile_write_dispatch_async's own comment.
   type(tile_write_job_t), target, save :: write_job(0:1)

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

   ! --- T4c: output-cube bookkeeping, generalized to n_outputs=6 (rather
   ! than rm_synthesis's own 2 named AMP/PHA fields) since rmclean_cubes
   ! writes CLEAN/RESID/RESTORED x AMP/PHA. Index convention below.
   integer, parameter :: n_outputs = 6
   integer, parameter :: idx_clean_amp = 1, idx_clean_pha = 2
   integer, parameter :: idx_resid_amp = 3, idx_resid_pha = 4
   integer, parameter :: idx_restored_amp = 5, idx_restored_pha = 6
   integer :: out_unit(n_outputs)
   character(len=600) :: out_path(n_outputs)
   integer(kind=8) :: out_datastart(n_outputs)
   logical :: out_is_open(n_outputs)
   integer :: nwriters, nwriters_eff

   integer :: ix, iy, k
   integer(kind=8) :: restore_plan_fwd, restore_plan_bwd
   integer :: n_pixels_done
   ! Aggregate CLEAN stopping-reason tally (T8, planning/
   ! RMCLEAN_INTEGRATION_PLAN.md): every CLEANed pixel (nvalid_p>=1 in
   ! clean_one_pixel) exits clean_complex with one of the 3 named
   ! stop_reason values ('niter'/'abs_flux'/'auto_nsigma') -- tallied
   ! here (atomic increments, same pattern as n_pixels_done) and
   ! reported once at run end, so "why did CLEAN stop" is answered at
   ! the population level without per-pixel log spam. n_iter_used_sum/
   ! min/max feed the same summary's mean/range.
   integer(kind=8) :: n_stopped_niter, n_stopped_abs_flux, n_stopped_auto_nsigma
   integer(kind=8) :: n_iter_used_sum
   integer :: n_iter_used_min, n_iter_used_max
   ! total_cleaned_flux_sum: sum, over every CLEANed pixel, of that
   ! pixel's own sum(|comp|) across all nrm bins -- how much flux CLEAN
   ! actually extracted into components, run-wide. peak_residual_*: each
   ! pixel's own final max(|resid|) (the peak left behind once CLEAN
   ! stopped), aggregated mean/min/max across pixels -- the population-
   ! level "how clean is clean" signal the user asked for alongside the
   ! stop-reason tally above.
   real(dp) :: total_cleaned_flux_sum
   real(dp) :: peak_residual_sum
   real(sp) :: peak_residual_max, peak_residual_min
   integer :: mask_pattern_cache_max
   ! T14 (docs/dev/RMCLEAN_INTEGRATION_PLAN.md): which policy decides
   ! what to evict once mask_pattern_cache_max is reached. 'belady'
   ! (default) is the offline-optimal policy, informed by a one-time
   ! whole-mask-cube pre-scan (run_pattern_prescan); 'hitcount' is a
   ! plain least-frequently-used fallback needing no pre-scan at all --
   ! see either policy's own call site for the full story.
   character(len=16) :: cache_eviction_policy
   ! T14: the whole-mask-cube pre-scan's own output -- every distinct
   ! pattern's future occurrence timeline, built once by
   ! run_pattern_prescan before real compute starts, consulted (never
   ! rebuilt) by the Belady eviction decision during the real run.
   ! Allocated/populated only when cache_eviction_policy='belady'; left
   ! entirely untouched (registry%n_entries stays 0) under 'hitcount',
   ! which needs no registry at all.
   type(pattern_registry_t) :: pattern_registry
   real(dp) :: t_stage
   ! Per-thread swim-lane instrumentation (planning-doc ticket) -- see
   ! convolve_cubes.f90's own write_convolved_file for the full
   ! rationale. tile_seq (already an existing running tile counter,
   ! incremented once per tile before this parallel region) doubles as
   ! the block index here -- no separate counter needed.
   integer :: tid_local
   real(dp) :: t_thread_start, t_thread_elapsed
   character(len=160) :: thread_msg

   ! --- Mask-pattern -> rmsf_table_t cache (planning/RMCLEAN_INTEGRATION_
   ! PLAN.md decision 10) -- T10b: built INCREMENTALLY now, one tile's
   ! own mask subrange at a time (update_mask_pattern_cache_for_tile),
   ! called serially right after that tile's own mask_tile is read and
   ! BEFORE that tile's own parallel per-pixel CLEAN loop starts, so
   ! every thread's own lookup (cache_lookup_readonly) is still read-only
   ! and race-free within that tile -- insertion for a LATER tile happens
   ! strictly after all of THIS tile's own lookups have completed
   ! (sequential tile loop), so the same "insertion never races a
   ! lookup" safety argument as the old one-upfront-pass design still
   ! holds, just repeated once per tile instead of once globally. See
   ! table_cache_entry_t's own comment.
   type :: table_cache_entry_t
      integer(kind=8) :: hash = 0_8
      ! T14 increment 9: allocated (holds this slot's own pattern
      ! bytes) under cache_eviction_policy='hitcount' only -- the only
      ! policy under which no registry exists to hold them instead.
      ! Under 'belady', this stays UNALLOCATED: the same bytes already
      ! live in pattern_registry%entries(registry_id)%pattern (Pass 0
      ! put them there before Pass 1 ever runs), so storing a second
      ! copy here would be pure duplication -- registry_id below
      ! reaches them instead, via cached_pattern_matches, the one place
      ! either lookup site compares a cached entry's own pattern.
      integer(kind=1), allocatable :: pattern(:)
      type(rmsf_table_t) :: table
      ! T14 (docs/dev/RMCLEAN_INTEGRATION_PLAN.md): cache_eviction_
      ! policy='hitcount' only -- counts lookups that hit this entry
      ! (starts at 1 on insertion, since that first use IS a hit).
      ! Least-frequently-used entry is the eviction candidate once the
      ! cache is full. Unused under policy='belady' (see that policy's
      ! own eviction logic instead).
      integer(kind=8) :: hit_count = 0_8
      ! T14 increment 7: cache_eviction_policy='belady' only -- which
      ! pattern_registry entry this cache row corresponds to, set once
      ! at insertion (registry_lookup, never re-looked-up on a hit).
      ! Lets a cache HIT advance that pattern's own registry timeline
      ! (registry_advance) without paying for a second hash lookup into
      ! a SEPARATE hash table on every hit -- only insertions/overflow
      ! (the rare case once the cache is warm) need a fresh
      ! registry_lookup. Unused (stays 0) under policy='hitcount'.
      integer :: registry_id = 0
   end type table_cache_entry_t
   type(table_cache_entry_t), allocatable :: cache_entries(:)
   integer, allocatable :: cache_buckets(:)
   integer :: n_cache_entries, n_cache_buckets
   ! Running totals across ALL tiles (T10b), printed once after the
   ! whole tile loop completes -- the old design printed this once,
   ! right after its own single global pre-scan; there is no longer one
   ! single pre-scan moment to print it at.
   integer :: n_distinct_patterns_total, n_overflow_pixels_total

   call parse_args(status)
   if (status.ne.0) stop 1

   call init_logging(log_level, timing_enabled, log_output_file, status)
   if (status.ne.0) then
      write(*,*) 'ERROR: cannot open log_output_file: ', trim(log_output_file)
      stop 1
   endif
   call log_message('info', 'startup', 'rmclean_cubes run started')
   write(*,'(A)') 'cache_eviction_policy='//trim(cache_eviction_policy)

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

   call resolve_abs_flux_floor(status)
   if (status.ne.0) stop 1
   if (have_auto_nsigma) then
      write(*,'(A,F0.4,A)') 'auto_nsigma: enabled, multiplier=',&
      &auto_nsigma_mult, ' (per-pixel noise sigma from each pixel''s'//&
      &' own full dirty spectrum, IQR-based -- see estimate_iqr_sigma,'//&
      &' src/rmclean.f90).'
   endif

   allocate(rm_samp(nrm))
   do k = 1, nrm
      rm_samp(k) = real(crval3_amp, sp) + real(k-1, sp)*real(cdelt3_amp, sp)
   enddo

   ! lsq_ref_compute: this program's OWN free choice of reference for the
   ! RMSF table/CLEAN computation (this file's own top comment) --
   ! defaults to 'mid' (recommended, minimizes the per-escalation search
   ! cost); 'native' means "use lsq_ref_native, no derotation needed",
   ! still available as an explicit opt-out.
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

   call plan_fourier_interp(nrm, nrm, restore_plan_fwd, restore_plan_bwd)

   ! T10b: plan_rmclean_tile now runs BEFORE the mask cache is
   ! initialised (it no longer needs the mask read first -- there is no
   ! more whole-cube read_mask_cube call at all), since its own
   ! bytes_per_tile_pixel budget now includes the mask's own per-pixel
   ! footprint too, and mask_tile below is allocated at the tile size it
   ! decides.
   call plan_rmclean_tile()

   ! T4b: same clamp convention as rm_synthesis's own io_read_threads --
   ! never more threads than OMP has, never more than there are RM planes
   ! to split across them (a tile's own nrm, not the outer tile count).
   ! Moved earlier than T4b's own original spot (still right after
   ! plan_rmclean_tile, just ahead of nwriters_eff now) since T14's own
   ! Pass-0 pre-scan needs io_read_threads_eff for its own read_mask_tile
   ! calls, and runs before nwriters_eff (an OUTPUT-side concern,
   ! irrelevant to Pass 0) is otherwise needed.
   io_read_threads_eff = max(1, min(io_read_threads, omp_get_max_threads()))
   io_read_threads_eff = min(io_read_threads_eff, nrm)
   if (io_read_threads.gt.io_read_threads_eff) then
      write(*,'(A,I0,A,I0)') 'WARNING: io_read_threads=', io_read_threads,&
      &' clamped to ', io_read_threads_eff
   endif

   ! T14: mask_tile itself moved earlier too (was allocated further
   ! down, after init_mask_pattern_cache) so Pass 0 and the real Pass-1
   ! tile loop share ONE allocation instead of doubling peak mask-tile
   ! memory.
   allocate(mask_tile(tile_ra,tile_dec,nchan))

   if (cache_eviction_policy.eq.'belady') call run_pattern_prescan()

   call init_mask_pattern_cache()

   ! T4c: same clamp convention as io_read_threads_eff above.
   nwriters_eff = max(1, min(nwriters, omp_get_max_threads()))
   nwriters_eff = min(nwriters_eff, nrm)
   if (nwriters.gt.nwriters_eff) then
      write(*,'(A,I0,A,I0)') 'WARNING: nwriters=', nwriters,&
      &' clamped to ', nwriters_eff
   endif

   allocate(re_tile(tile_ra,tile_dec,nrm), im_tile(tile_ra,tile_dec,nrm))
   allocate(clean_re_buf(tile_ra,tile_dec,nrm,2), clean_im_buf(tile_ra,tile_dec,nrm,2))
   allocate(resid_re_buf(tile_ra,tile_dec,nrm,2), resid_im_buf(tile_ra,tile_dec,nrm,2))
   allocate(restored_re_buf(tile_ra,tile_dec,nrm,2), restored_im_buf(tile_ra,tile_dec,nrm,2))

   ! Each out_unit(idx) is assigned its own genuinely distinct CFITSIO
   ! unit inside open_output_cube itself, via fitsio_unit_mod's
   ! safe_ftinit (FTGIOU+FTINIT) -- no caller-chosen numbering needed.
   call open_output_cube(ampfile, trim(outfile)//'.CLEAN.AMP.RMCUBE.FITS',&
   &nx, ny, nrm, idx_clean_amp, status)
   if (status.ne.0) stop 1
   call open_output_cube(ampfile, trim(outfile)//'.CLEAN.PHA.RMCUBE.FITS',&
   &nx, ny, nrm, idx_clean_pha, status)
   if (status.ne.0) stop 1
   call open_output_cube(ampfile, trim(outfile)//'.RESID.AMP.RMCUBE.FITS',&
   &nx, ny, nrm, idx_resid_amp, status)
   if (status.ne.0) stop 1
   call open_output_cube(ampfile, trim(outfile)//'.RESID.PHA.RMCUBE.FITS',&
   &nx, ny, nrm, idx_resid_pha, status)
   if (status.ne.0) stop 1
   call open_output_cube(ampfile, trim(outfile)//'.RESTORED.AMP.RMCUBE.FITS',&
   &nx, ny, nrm, idx_restored_amp, status)
   if (status.ne.0) stop 1
   call open_output_cube(ampfile, trim(outfile)//'.RESTORED.PHA.RMCUBE.FITS',&
   &nx, ny, nrm, idx_restored_pha, status)
   if (status.ne.0) stop 1

   ! T4a/T4d: sequential tiles, each a strict read (single-threaded or T4b
   ! io_read_threads-parallel) / compute (OpenMP-parallel across the
   ! tile's own pixels) / write (inline, or -- io_overlap=y -- on a
   ! background thread while the NEXT tile's read/compute proceeds)
   ! sequence -- see plan_rmclean_tile's own comment and this ticket's
   ! own top-of-file/planning-doc rationale for why this shape (spatial
   ! tile, full RM depth) rather than reproject_cubes.f90's own
   ! depth-plane blocking.
   n_pixels_done = 0
   n_stopped_niter = 0_8
   n_stopped_abs_flux = 0_8
   n_stopped_auto_nsigma = 0_8
   n_iter_used_sum = 0_8
   n_iter_used_min = huge(n_iter_used_min)
   n_iter_used_max = 0
   total_cleaned_flux_sum = 0.0_dp
   peak_residual_sum = 0.0_dp
   peak_residual_max = 0.0_sp
   peak_residual_min = huge(peak_residual_min)
   tile_seq = 0
   iy_tile_beg = 1
   do while (iy_tile_beg.le.ny)
      ty = min(tile_dec, ny-iy_tile_beg+1)
      ix_tile_beg = 1
      do while (ix_tile_beg.le.nx)
         tx = min(tile_ra, nx-ix_tile_beg+1)

         ! T4d buffer-reuse join: slot cur_slot is about to be (re)written
         ! by this tile's own compute step; if the write from two tiles
         ! ago (the last user of this same slot) is still in flight, join
         ! it now -- the one synchronisation point that makes the double
         ! buffering safe. With io_overlap=n, write_pending is never set,
         ! so this is a no-op (see tile_write_job_t's own comment).
         cur_slot = mod(tile_seq, 2)
         tile_seq = tile_seq + 1
         call timer_start(t_stage)
         if (io_overlap .and. write_pending(cur_slot)) then
            call tile_write_join(write_thread_id(cur_slot))
            write_pending(cur_slot) = .false.
         endif
         call timer_stop('tile_write_join', t_stage)
         clean_re_tile => clean_re_buf(:,:,:,cur_slot+1)
         clean_im_tile => clean_im_buf(:,:,:,cur_slot+1)
         resid_re_tile => resid_re_buf(:,:,:,cur_slot+1)
         resid_im_tile => resid_im_buf(:,:,:,cur_slot+1)
         restored_re_tile => restored_re_buf(:,:,:,cur_slot+1)
         restored_im_tile => restored_im_buf(:,:,:,cur_slot+1)

         call timer_start(t_stage)
         call read_amp_pha_tile(ampfile, phafile, ix_tile_beg, iy_tile_beg,&
         &tx, ty, re_tile(1:tx,1:ty,:), im_tile(1:tx,1:ty,:), status)
         if (status.ne.0) stop 1
         ! T10b: this tile's own mask subregion, read via the same
         ! tiled/io_read_threads-aware mechanism as AMP/PHA above,
         ! replacing the old once-upfront whole-cube read_mask_cube.
         call read_mask_tile(maskfile, nx, ny, nchan, ix_tile_beg,&
         &iy_tile_beg, tx, ty, io_read_threads_eff,&
         &mask_tile(1:tx,1:ty,:), status)
         call timer_stop('tile_read', t_stage)
         if (status.ne.0) stop 1

         ! T10b: incremental mask-pattern-cache pre-scan for THIS tile
         ! only, serial, before the parallel CLEAN loop below starts --
         ! replaces the old one-upfront-global build_mask_pattern_cache
         ! call (see update_mask_pattern_cache_for_tile's own comment).
         call update_mask_pattern_cache_for_tile(tx, ty)

         call timer_start(t_stage)
         !$omp parallel default(shared)&
         !$omp& private(ix, iy, tid_local, t_thread_start, t_thread_elapsed, thread_msg)
         tid_local = omp_get_thread_num()
         t_thread_start = omp_get_wtime()
         write(thread_msg,'(A,I0,A,I0,A,I0)') 'thread_timing stage=clean event=start tid=',&
         &tid_local,' block=',tile_seq,' unit_count=',tx*ty
         call log_message('debug','tile_thread',trim(thread_msg))
         !$omp do collapse(2) schedule(dynamic)
         do iy = 1, ty
            do ix = 1, tx
               call clean_one_pixel(ix, iy, ix_tile_beg, iy_tile_beg)
            enddo
         enddo
         !$omp end do
         t_thread_elapsed = (omp_get_wtime() - t_thread_start) * 1000.0_dp
         ! F10.3 only holds up to 6 integer digits (max 999999.999 ms ~=
         ! 16.7 min) before Fortran prints '*' fill instead of the value
         ! -- too narrow once a single block's own compute stage runs
         ! longer than that (found on a real full-cube run: one block
         ! took 55min). F0.3 (auto-width) has no overflow ceiling and no
         ! fixed-width padding, matching this file's own convention for
         ! other floating-point log fields (F0.6, F0.4).
         write(thread_msg,'(A,I0,A,I0,A,I0,A,F0.3)') 'thread_timing stage=clean event=done tid=',&
         &tid_local,' block=',tile_seq,' unit_count=',tx*ty,' dur_ms=',t_thread_elapsed
         call log_message('debug','tile_thread',trim(thread_msg))
         !$omp end parallel
         call timer_stop('tile_compute', t_stage)

         write(*,'(A,I0,A,I0,A)') 'Block ', tile_seq, ' of ', n_blocks_total,&
         &' processed.'

         write_job(cur_slot)%ix_tile_beg = ix_tile_beg
         write_job(cur_slot)%iy_tile_beg = iy_tile_beg
         write_job(cur_slot)%tx = tx
         write_job(cur_slot)%ty = ty
         write_job(cur_slot)%re_ptr(1)%p => clean_re_tile(1:tx,1:ty,:)
         write_job(cur_slot)%im_ptr(1)%p => clean_im_tile(1:tx,1:ty,:)
         write_job(cur_slot)%re_ptr(2)%p => resid_re_tile(1:tx,1:ty,:)
         write_job(cur_slot)%im_ptr(2)%p => resid_im_tile(1:tx,1:ty,:)
         write_job(cur_slot)%re_ptr(3)%p => restored_re_tile(1:tx,1:ty,:)
         write_job(cur_slot)%im_ptr(3)%p => restored_im_tile(1:tx,1:ty,:)

         call timer_start(t_stage)
         if (io_overlap) then
            ! T4d handle-safety join: before any NEW write is dispatched,
            ! whichever write is currently outstanding (either slot) is
            ! joined first, unconditionally -- all 6 output files' CFITSIO/
            ! raw-write handles are shared across slots regardless of
            ! which slot is being written, so the per-slot join above only
            ! guards buffer reuse (2 tiles apart); it does NOT guarantee
            ! the *previous* tile's write (a different slot) has finished.
            ! Ported directly from rm_synthesis.f90's own identical rule
            ! (see its own comment for the SIGSEGV this prevents).
            if (write_pending(0)) then
               call tile_write_join(write_thread_id(0))
               write_pending(0) = .false.
            endif
            if (write_pending(1)) then
               call tile_write_join(write_thread_id(1))
               write_pending(1) = .false.
            endif
            call timer_stop('tile_write_join', t_stage)
            call timer_start(t_stage)
            call tile_write_dispatch_async(write_job(cur_slot),&
            &write_thread_id(cur_slot), write_dispatched_ok)
            write_pending(cur_slot) = write_dispatched_ok
         else
            call do_tile_write(write_job(cur_slot))
         endif
         call timer_stop('tile_write', t_stage)

         ix_tile_beg = ix_tile_beg + tx
      enddo
      iy_tile_beg = iy_tile_beg + ty
   enddo
   ! T10b: printed once here now, accumulated across every tile's own
   ! update_mask_pattern_cache_for_tile call -- the old design printed
   ! this once, right after its own single global pre-scan; there is no
   ! longer one single pre-scan moment to print it at.
   write(*,'(A,I0,A,I0,A)') 'Mask-pattern cache: ', n_distinct_patterns_total,&
   &' distinct pattern(s) cached'
   if (n_overflow_pixels_total.gt.0) then
      write(*,'(A,I0,A)') '  (', n_overflow_pixels_total, ' pixel(s) past'//&
      &' mask_pattern_cache_max fall back to a one-off table each)'
   endif
   ! T14 increment 9: direct, printed evidence for the actual point of
   ! this increment -- how many resident cache slots still carry their
   ! own local pattern-bytes copy. Expected 0 under belady (every
   ! slot's own pattern lives in the registry instead) and equal to
   ! n_cache_entries under hitcount (no registry exists to delegate to,
   ! so every slot must keep its own copy) -- printed rather than
   ! silently assumed, so a real regression here doesn't hide invisibly
   ! behind otherwise-passing correctness tests.
   block
      integer :: n_pattern_allocated, kk
      n_pattern_allocated = 0
      do kk = 1, n_cache_entries
         if (allocated(cache_entries(kk)%pattern)) n_pattern_allocated = n_pattern_allocated + 1
      enddo
      write(*,'(A,I0,A,I0,A)') '  (pattern bytes locally stored for ',&
      &n_pattern_allocated, ' of ', n_cache_entries, ' resident cache'//&
      &' slot(s) -- 0 expected under cache_eviction_policy=belady)'
   end block
   write(*,'(A,I0,A)') 'CLEANed ', n_pixels_done, ' pixels.'
   if (n_stopped_niter+n_stopped_abs_flux+n_stopped_auto_nsigma .gt. 0_8) then
      block
         real(dp) :: n_total_stopped_dp
         n_total_stopped_dp = real(n_stopped_niter+n_stopped_abs_flux+&
         &n_stopped_auto_nsigma, dp)
         write(*,'(A)') 'CLEAN stop-reason summary:'
         write(*,'(A,I0,A,F6.2,A)') '  hit niter cap:        ',&
         &n_stopped_niter, ' (',&
         &100.0_dp*real(n_stopped_niter,dp)/n_total_stopped_dp, '%)'
         write(*,'(A,I0,A,F6.2,A)') '  stopped at abs_flux:  ',&
         &n_stopped_abs_flux, ' (',&
         &100.0_dp*real(n_stopped_abs_flux,dp)/n_total_stopped_dp, '%)'
         write(*,'(A,I0,A,F6.2,A)') '  stopped at auto_nsigma:',&
         &n_stopped_auto_nsigma, ' (',&
         &100.0_dp*real(n_stopped_auto_nsigma,dp)/n_total_stopped_dp, '%)'
         write(*,'(A,F0.2,A,I0,A,I0,A,I0,A)') '  n_iter_used: mean=',&
         &real(n_iter_used_sum,dp)/n_total_stopped_dp,&
         &' min=', n_iter_used_min, ' max=', n_iter_used_max,&
         &' (niter cap=', niter, ')'
      end block
   endif

   ! Join any tile writes still in flight -- at most the last two tiles'
   ! writes can reach here undispatched-for-join (each slot's write is
   ! only joined when that slot is reused two tiles later, so the final
   ! tile, and the one before it in the other slot, never gets that
   ! second chance).
   call timer_start(t_stage)
   if (io_overlap) then
      if (write_pending(0)) then
         call tile_write_join(write_thread_id(0))
         write_pending(0) = .false.
      endif
      if (write_pending(1)) then
         call tile_write_join(write_thread_id(1))
         write_pending(1) = .false.
      endif
   endif
   call timer_stop('tile_write_join', t_stage)

   call destroy_fourier_interp_plan(restore_plan_fwd, restore_plan_bwd)

   do k = 1, n_outputs
      call close_output_cube(k)
   enddo

   write(*,'(A)') 'OK: rmclean_cubes complete.'
   call timer_report_summary()
   call log_message('info', 'finalize', 'rmclean_cubes run completed')

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
      have_abs_flux_floor = .false.
      abs_flux_floor = 0.0_sp
      abs_flux_floor_raw_value = 0.0_sp
      abs_flux_floor_unit_scale = -1.0_sp
      have_auto_nsigma = .false.
      auto_nsigma_mult = 0.0_sp
      min_samples_per_fwhm = 2.0_sp
      refine_nsigma = 3.0_sp
      trace_ix = 0
      trace_iy = 0
      log_every = 50
      table_oversample = 20
      have_restore_fwhm = .false.
      restore_fwhm_override = 0.0_sp
      lsq_ref_report_mode_sel = lsq_ref_report_intrinsic
      have_lsq_ref_report_value = .false.
      lsq_ref_report_value = 0.0_sp
      lsq_ref_compute_mode = 'mid'
      have_lsq_ref_compute_value = .false.
      lsq_ref_compute_value = 0.0_sp
      mask_pattern_cache_max = 4096
      cache_eviction_policy = 'belady'
      tile_ra = 0
      tile_dec = 0
      tile_auto = .true.
      mem_frac_ram = 0.25_sp
      min_valid_chan_frac = 0.0_sp
      io_read_threads = 1
      nwriters = 1
      io_overlap = .false.
      log_level = 'info'
      timing_enabled = .false.
      log_output_file = ' '
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
      ! abs_flux_floor=/auto_nsigma= are independently optional and MAY
      ! be combined (whichever fires first wins, per-pixel) -- niter is
      ! always the backstop. Neither given is also valid (niter-only,
      ! same precedent as the old thresh=0.0 -- see planning/
      ! RMCLEAN_INTEGRATION_PLAN.md's "Choosing parameters" section) but
      ! worth flagging, since it's easy to omit both by accident.
      if (.not.have_abs_flux_floor .and. .not.have_auto_nsigma) then
         write(*,'(A)') 'NOTE: neither abs_flux_floor= nor auto_nsigma='//&
         &' given -- CLEAN will run every pixel for the full niter'//&
         &' iterations with no early stop.'
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
      case ('abs_flux_floor')
         call parse_flux_value(val, abs_flux_floor_raw_value,&
         &abs_flux_floor_unit_scale, ios)
         if (ios.ne.0) then
            write(*,*) 'ERROR: abs_flux_floor must be a number, optionally'//&
            &' suffixed with a unit (Jy, mJy, or uJy) with no space,'//&
            &' e.g. abs_flux_floor=0.01 or abs_flux_floor=10mJy'
            status = -1
            return
         endif
         have_abs_flux_floor = .true.
      case ('auto_nsigma')
         read(val, *, iostat=ios) auto_nsigma_mult
         if (ios.ne.0 .or. auto_nsigma_mult.le.0.0_sp) then
            write(*,*) 'ERROR: auto_nsigma must be a positive number'
            status = -1
            return
         endif
         have_auto_nsigma = .true.
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
      case ('trace_ix')
         read(val, *, iostat=ios) trace_ix
         if (ios.ne.0 .or. trace_ix.lt.0) then
            write(*,*) 'ERROR: trace_ix must be a non-negative integer'
            status = -1
            return
         endif
      case ('trace_iy')
         read(val, *, iostat=ios) trace_iy
         if (ios.ne.0 .or. trace_iy.lt.0) then
            write(*,*) 'ERROR: trace_iy must be a non-negative integer'
            status = -1
            return
         endif
      case ('log_every')
         read(val, *, iostat=ios) log_every
         if (ios.ne.0 .or. log_every.lt.1) then
            write(*,*) 'ERROR: log_every must be a positive integer'
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
      case ('cache_eviction_policy')
         select case (trim(val))
         case ('belady', 'hitcount')
            cache_eviction_policy = trim(val)
         case default
            write(*,*) 'ERROR: cache_eviction_policy must be one of'//&
            &' belady|hitcount'
            status = -1
            return
         end select
      case ('mem_frac_ram')
         read(val, *, iostat=ios) mem_frac_ram
         if (ios.ne.0 .or. mem_frac_ram.le.0.0_sp .or. mem_frac_ram.gt.0.95_sp) then
            write(*,*) 'ERROR: mem_frac_ram must be > 0 and <= 0.95'
            status = -1
            return
         endif
      case ('min_valid_chan_frac')
         read(val, *, iostat=ios) min_valid_chan_frac
         if (ios.ne.0 .or. min_valid_chan_frac.lt.0.0_sp .or. min_valid_chan_frac.gt.1.0_sp) then
            write(*,*) 'ERROR: min_valid_chan_frac must be >= 0 and <= 1'
            status = -1
            return
         endif
      case ('tile_ra')
         read(val, *, iostat=ios) tile_ra
         if (ios.ne.0 .or. tile_ra.lt.0) then
            write(*,*) 'ERROR: tile_ra must be a non-negative integer'
            status = -1
            return
         endif
      case ('tile_dec')
         read(val, *, iostat=ios) tile_dec
         if (ios.ne.0 .or. tile_dec.lt.0) then
            write(*,*) 'ERROR: tile_dec must be a non-negative integer'
            status = -1
            return
         endif
      case ('tile_auto')
         tile_auto = flag_from_value_rmclean(val)
      case ('io_read_threads')
         read(val, *, iostat=ios) io_read_threads
         if (ios.ne.0 .or. io_read_threads.lt.1) then
            write(*,*) 'ERROR: io_read_threads must be an integer >= 1'
            status = -1
            return
         endif
      case ('nwriters')
         read(val, *, iostat=ios) nwriters
         if (ios.ne.0 .or. nwriters.lt.1) then
            write(*,*) 'ERROR: nwriters must be an integer >= 1'
            status = -1
            return
         endif
      case ('io_overlap')
         io_overlap = flag_from_value_rmclean(val)
      case ('log_level')
         log_level = trim(val)
      case ('timing_enabled')
         timing_enabled = flag_from_value_logging(val)
      case ('log_output_file')
         log_output_file = trim(val)
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
      &' maskfile=<f> outfile=<base>'
      write(*,'(A)') '    [abs_flux_floor=<v>] [auto_nsigma=<n>]'
      write(*,'(A)') '    [niter=<n>] [gain=<g>] [min_samples_per_fwhm=<f>]'//&
      &' [refine_nsigma=<f>] [table_oversample=<n>] [restore_fwhm=<v>]'
      write(*,'(A)') '    [trace_ix=<n>] [trace_iy=<n>] [log_every=<n>]'
      write(*,'(A)') '    [lsq_ref_report_mode=intrinsic|centroid|min|max'//&
      &'|mid|fixed] [lsq_ref_report_value=<v>]'
      write(*,'(A)') '    [lsq_ref_compute_mode=native|zero|centroid|min'//&
      &'|max|mid|fixed] [lsq_ref_compute_value=<v>]'
      write(*,'(A)') '    [mask_pattern_cache_max=<n>]'//&
      &' [cache_eviction_policy=belady|hitcount] [mem_frac_ram=<f>]'//&
      &' [min_valid_chan_frac=<f>]'//&
      &' [tile_ra=<n>] [tile_dec=<n>] [tile_auto=y|n]'//&
      &' [io_read_threads=<n>] [nwriters=<n>] [io_overlap=y|n]'
      write(*,'(A)') '   or: rmclean_cubes --config <cfgfile>'
      write(*,'(A)') '   or: rmclean_cubes --help | -h'
      write(*,'(A)') ''
      write(*,'(A)') 'ampfile/phafile: the dirty AMP.RMCUBE.FITS/'//&
      &'PHA.RMCUBE.FITS pair rm_synthesis itself wrote.'
      write(*,'(A)') 'maskfile: the matching MASK.CUBE.FITS (with its own'//&
      &' CHANFREQ binary table).'
      write(*,'(A)') 'outfile: base name for the 6 output cubes'//&
      &' (<outfile>.CLEAN/.RESID/.RESTORED.AMP/PHA.RMCUBE.FITS).'
      write(*,'(A)') 'CLEAN stopping criteria -- three independent,'//&
      &' unambiguously-named conditions (planning/RMCLEAN_INTEGRATION_'//&
      &'PLAN.md T8), checked in this order every iteration; niter is'//&
      &' always the hard backstop, the other two are each optional and'//&
      &' MAY be combined (whichever fires first wins for that pixel;'//&
      &' neither given is valid too -- niter alone governs).'
      write(*,'(A)') 'abs_flux_floor: stop the instant a pixel''s own peak'//&
      &' amplitude drops to/below this literal fixed flux value -- native'//&
      &' AMP-cube units for a bare number, or convert from Jy/mJy/uJy'//&
      &' with no space (e.g. abs_flux_floor=10mJy) against the AMP'//&
      &' cube''s own BUNIT. No noise/baseline adjustment -- a pure'//&
      &' "brightest remaining feature below X" comparison.'
      write(*,'(A)') 'auto_nsigma: stop the instant a pixel''s own peak'//&
      &' amplitude drops to/below auto_nsigma x that SAME pixel''s own'//&
      &' noise sigma. The sigma is estimated ONCE per pixel (not'//&
      &' recomputed per iteration) from that pixel''s own FULL dirty'//&
      &' amplitude spectrum''s interquartile range, converted to sigma'//&
      &' via the Rayleigh-distribution analytic relation (dirty'//&
      &' amplitude is Rayleigh-, not Gaussian-, distributed).'
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
      write(*,'(A)') 'trace_ix/trace_iy (default 0, disabled): 1-indexed'//&
      &' GLOBAL pixel to log a per-iteration CLEAN trend for (peak_val/'//&
      &'rms_val/cumulative cleaned flux, plus the exact stop_reason --'//&
      &' niter/abs_flux/auto_nsigma). Debugging/tuning aid, not meant'//&
      &' for whole-cube production runs (one traced pixel only).'
      write(*,'(A)') 'log_every (default 50): throttles trace_ix/trace_iy'//&
      &' iteration lines to every log_every-th iteration; iteration 1'//&
      &' (pre-CLEAN state) and the final iteration are always logged'//&
      &' regardless.'
      write(*,'(A)') 'restore_fwhm (optional): override the restoring'//&
      &' beam FWHM (rad/m^2); default derived from'//&
      &' compute_rmsf_fwhm_multiband.'
      write(*,'(A)') 'lsq_ref_report_mode (default intrinsic, i.e.'//&
      &' lsq_ref_report=0.0): where to report the derotated chi0/'//&
      &' restored phase -- a safe, independent post-processing choice.'
      write(*,'(A)') 'lsq_ref_compute_mode (default mid, the recommended'//&
      &' choice): the reference this program''s own RMSF table/CLEAN'//&
      &' computation uses -- an exact, free choice (derotate_to_lsq_ref,'//&
      &' costs nothing in accuracy), and mid minimizes the per-pixel'//&
      &' escalation search cost (refine_peak_matched_filter''s own'//&
      &' max_offset-driven m_search) -- see this file''s own top comment'//&
      &' for the full story, including why this does NOT change the RM'//&
      &' grid size itself (fixed by whatever rm_synthesis already wrote).'
      write(*,'(A)') 'mask_pattern_cache_max (default 4096): pixels'//&
      &' sharing the same valid-channel mask pattern share one'//&
      &' rmsf_table_t, built once during a serial pre-scan; past this'//&
      &' many DISTINCT patterns, additional patterns fall back to a'//&
      &' one-off table per pixel (safety valve, not a correctness'//&
      &' issue -- just loses the reuse benefit).'
      write(*,'(A)') 'cache_eviction_policy (default belady): which'//&
      &' pattern to evict once mask_pattern_cache_max is reached rather'//&
      &' than simply refusing new entries. belady: offline-optimal,'//&
      &' informed by a one-time whole-mask-cube pre-scan (evicts'//&
      &' whichever cached pattern is not needed again for the longest'//&
      &' time, or ever) -- no tuning parameter, provably the best'//&
      &' possible for a given cache size. hitcount: a plain'//&
      &' least-frequently-used fallback needing no pre-scan at all --'//&
      &' cheaper to start, but not competitive with belady; intended'//&
      &' as an opt-out, not the recommended choice.'
      write(*,'(A)') 'mem_frac_ram (default 0.25): fraction of total'//&
      &' system RAM budgeted for one tile''s own read+compute+write'//&
      &' working set (2 input + 6 output RM-depth arrays per pixel) --'//&
      &' same scheme/key name as rm_synthesis''s own mem_frac_ram'//&
      &' (planning ticket T4).'
      write(*,'(A)') 'min_valid_chan_frac (default 0.0, i.e. off): pixels'//&
      &' whose own valid-channel fraction (summed across ALL bands,'//&
      &' out of the full channel count) falls below this are never'//&
      &' CLEANed -- output NaN across every RM-bin instead. 0.0 means'//&
      &' every pixel with at least 1 valid channel is CLEANed (today''s'//&
      &' original behaviour); e.g. 0.7 skips any pixel with less than'//&
      &' 70% of channels valid. Sparse-coverage pixels have a'//&
      &' coarser/noisier RMSF regardless of compute spent on them --'//&
      &' this also shrinks Pass 0''s own pattern registry and the'//&
      &' runtime cache, since skipped pixels are never registered or'//&
      &' cached at all (docs/dev/RMCLEAN_INTEGRATION_PLAN.md T17).'
      write(*,'(A)') 'tile_ra/tile_dec (default 0/0, i.e. auto): manual'//&
      &' tile size override (pixels); ignored unless tile_auto=n. Auto'//&
      &' policy (tile_auto=y, the default) packs full-RA Dec strips --'//&
      &' tile_ra=nx, as many Dec rows as mem_frac_ram allows -- falling'//&
      &' back to RA-subdivided single rows only if a single full-RA row'//&
      &' does not fit; same policy as rm_synthesis''s own plan_tile.'
      write(*,'(A)') 'tile_auto (default y): y = ignore tile_ra/tile_dec'//&
      &' and derive them from mem_frac_ram; n = use tile_ra/tile_dec as'//&
      &' given (still clamped to the image size).'
      write(*,'(A)') 'io_read_threads (default 1): parallel readonly'//&
      &' chunked reads of each tile''s AMP/PHA RM-depth range -- same'//&
      &' scheme as rm_synthesis''s own io_read_threads; clamped to'//&
      &' min(io_read_threads, OMP thread count, tile''s own nrm).'
      write(*,'(A)') 'nwriters (default 1): parallel chunked'//&
      &' writes of each tile''s 6 output cubes -- same scheme as'//&
      &' rm_synthesis''s own nwriters: >1 switches to raw'//&
      &' stream writes at computed byte offsets, bypassing CFITSIO''s'//&
      &' own ftpsse for the pixel data (concurrent CFITSIO read-write'//&
      &' handles on one file are unsafe -- see write_output_tile''s own'//&
      &' comment); clamped like io_read_threads.'
      write(*,'(A)') 'io_overlap (default n): y = write each tile''s 6'//&
      &' output cubes on a background thread while the NEXT tile''s'//&
      &' read/compute proceeds -- same scheme as rm_synthesis''s own'//&
      &' io_overlap; works with any nwriters setting.'
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

   function flag_from_value_rmclean(val) result(flag)
      !! Same convention as rm_synthesis_mod.f90's own flag_from_value:
      !! first non-blank character '1'/'y'/'Y'/'t'/'T' -> true, anything
      !! else (including blank) -> false.
      character(len=*), intent(in) :: val
      logical :: flag
      character(len=64) :: t
      integer :: i

      t = adjustl(val)
      do i = 1, len_trim(t)
         t(i:i) = achar(iachar(t(i:i)) - merge(32, 0, t(i:i).ge.'A'.and.t(i:i).le.'Z'))
      enddo
      flag = .false.
      if (len_trim(t).eq.0) return
      if (t(1:1).eq.'1' .or. t(1:1).eq.'y' .or. t(1:1).eq.'t') flag = .true.
   end function flag_from_value_rmclean

   subroutine parse_flux_value(val, raw_value, unit_scale, ios)
      !! Splits an abs_flux_floor= value into its leading numeric part and an
      !! optional trailing unit suffix (Jy/mJy/uJy, case-insensitive, no
      !! space -- e.g. "10mJy", "0.01Jy", or a bare "0.01"). Scans
      !! BACKWARD from the end of the string while the character is a
      !! letter -- the boundary between the (purely alphabetic) unit
      !! suffix and the numeric part is the first non-letter found this
      !! way, so an exponent's own 'e'/'E' (embedded, never trailing --
      !! a valid number never ENDS in a bare e/E) is never mistaken for
      !! part of a unit token. unit_scale is left at its -1.0 sentinel
      !! (caller-provided default, meaning "no unit given -- use the
      !! value as-is in the AMP cube's own native units", the historical
      !! bare-number behaviour) when there is no trailing letter run.
      character(len=*), intent(in) :: val
      real(sp), intent(out) :: raw_value
      real(sp), intent(inout) :: unit_scale
      integer, intent(out) :: ios
      character(len=len(val)) :: v
      character(len=len(val)) :: numeric_part, unit_part
      integer :: i, nchar_v, unit_len
      character :: c

      ios = 0
      raw_value = 0.0_sp
      v = adjustl(val)
      nchar_v = len_trim(v)
      if (nchar_v.eq.0) then
         ios = -1
         return
      endif

      unit_len = 0
      do i = nchar_v, 1, -1
         c = v(i:i)
         if ((c.ge.'a'.and.c.le.'z') .or. (c.ge.'A'.and.c.le.'Z')) then
            unit_len = unit_len + 1
         else
            exit
         endif
      enddo
      if (unit_len.ge.nchar_v) then
         ios = -1
         return
      endif

      numeric_part = v(1:nchar_v-unit_len)
      read(numeric_part, *, iostat=ios) raw_value
      if (ios.ne.0) return

      if (unit_len.gt.0) then
         unit_part = v(nchar_v-unit_len+1:nchar_v)
         unit_scale = flux_unit_scale(unit_part(1:unit_len), ios)
      endif
   end subroutine parse_flux_value

   function flux_unit_scale(unit_str, ios) result(scale)
      !! Jy-relative scale factor for a recognised flux-unit token
      !! (case-insensitive). ios is set nonzero for anything else.
      character(len=*), intent(in) :: unit_str
      integer, intent(out) :: ios
      real(sp) :: scale
      character(len=len(unit_str)) :: u
      integer :: i

      u = unit_str
      do i = 1, len_trim(u)
         u(i:i) = achar(iachar(u(i:i)) + merge(32, 0, u(i:i).ge.'A'.and.u(i:i).le.'Z'))
      enddo
      ios = 0
      select case (trim(u))
      case ('jy')
         scale = 1.0_sp
      case ('mjy')
         scale = 1.0e-3_sp
      case ('ujy')
         scale = 1.0e-6_sp
      case default
         scale = 1.0_sp
         ios = -1
      end select
   end function flux_unit_scale

   subroutine read_bunit_keyword(filename, bunit_out, found)
      !! Reads the AMP cube's own BUNIT keyword (rm_synthesis passes this
      !! through from the input Q/U cube, rm_synthesis.f90:2953/2958) --
      !! used to resolve an explicit abs_flux_floor= unit against the
      !! cube's OWN native unit. Typical ASKAP-style
      !! values look like "Jy/beam" -- only the leading flux-unit token
      !! (up to the first '/', if any) is meaningful for the Jy/mJy/uJy
      !! scale comparison here; the per-beam/per-pixel denominator is the
      !! same on both sides of any ratio this program computes, so it is
      !! read but not otherwise interpreted.
      character(len=*), intent(in) :: filename
      character(len=*), intent(out) :: bunit_out
      logical, intent(out) :: found
      integer :: unit, blocksize, fitsstat
      character(len=80) :: comment

      bunit_out = ' '
      found = .false.
      fitsstat = 0
      ! safe_ftopen (fitsio_unit_mod.f90, FTGIOU+FTOPEN) -- FTOPEN
      ! requires a caller-assigned CFITSIO unit number, it does NOT
      ! auto-allocate one (no newunit= support, unlike plain Fortran
      ! open). An uninitialized `unit` here (the original bug, before
      ! this whole file moved to FTGIOU-based allocation) happened to run
      ! without visibly failing on the tiny 32x32 test fixture but
      ! corrupted CFITSIO's internal unit table on a real cube,
      ! SIGSEGVing deep inside libcfitsio -- found via the real ~46GB
      ! Jennifer end-to-end run.
      call safe_ftopen(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         call free_fits_unit(unit)
         return
      endif
      call FTGKYS(unit, 'BUNIT', bunit_out, comment, fitsstat)
      if (fitsstat.eq.0) found = .true.
      fitsstat = 0
      call safe_ftclos(unit, fitsstat)
   end subroutine read_bunit_keyword

   subroutine resolve_abs_flux_floor(status)
      !! Finalises the module-level `abs_flux_floor` (native AMP-cube
      !! units) from abs_flux_floor='s given value+unit -- called once,
      !! after read_amp_pha_geometry (nx/ny/nrm known) and before the main
      !! tile loop. A no-op (status=0, abs_flux_floor left at its init
      !! value) when have_abs_flux_floor is .false. -- auto_nsigma's own
      !! per-pixel percentile-sigma estimate (estimate_percentile_sigma,
      !! src/rmclean.f90) needs no whole-cube pre-scan or resolution
      !! step, unlike the old threshold_snr's noise_floor pre-scan this
      !! subroutine used to also perform (planning/RMCLEAN_INTEGRATION_
      !! PLAN.md T8/T9 -- that whole mechanism is retired, not just
      !! renamed).
      !! Prints exactly what it resolved to, so a run's own stdout log
      !! (scripts/run_pipeline.sh's provenance capture) records the
      !! actual number used, not just the cfg's own request.
      integer, intent(out) :: status
      character(len=80) :: bunit, bunit_token
      logical :: have_bunit
      real(sp) :: native_scale
      integer :: ios, slash_pos

      status = 0
      if (.not.have_abs_flux_floor) return

      call read_bunit_keyword(ampfile, bunit, have_bunit)
      native_scale = 1.0_sp
      if (have_bunit .and. len_trim(bunit).gt.0) then
         bunit_token = adjustl(bunit)
         slash_pos = index(bunit_token, '/')
         if (slash_pos.gt.0) bunit_token = bunit_token(1:slash_pos-1)
         native_scale = flux_unit_scale(trim(bunit_token), ios)
         if (ios.ne.0) then
            write(*,'(A)') 'WARNING: '//trim(ampfile)//' BUNIT="'//&
            &trim(bunit)//'" not recognised as Jy/mJy/uJy -- assuming Jy'//&
            &' for abs_flux_floor unit conversion.'
            native_scale = 1.0_sp
         endif
      else
         write(*,'(A)') 'WARNING: no BUNIT keyword in '//trim(ampfile)//&
         &' -- assuming Jy for abs_flux_floor unit conversion.'
      endif

      if (abs_flux_floor_unit_scale.gt.0.0_sp) then
         abs_flux_floor = abs_flux_floor_raw_value *&
         &(abs_flux_floor_unit_scale/native_scale)
         write(*,'(A,F0.6,A,F0.6,A)') 'abs_flux_floor: ',&
         &abs_flux_floor_raw_value, ' (given unit) -> ', abs_flux_floor,&
         &' (native AMP-cube units, BUNIT="'//trim(bunit)//'").'
      else
         abs_flux_floor = abs_flux_floor_raw_value
      endif
   end subroutine resolve_abs_flux_floor

   subroutine get_mem_total_kb(mem_total_kb)
      !! Total system RAM in kB, from /proc/meminfo's MemTotal -- verbatim
      !! adaptation of reproject_cubes.f90's own get_mem_total_kb (itself
      !! matching rm_synthesis_mod.f90's own tile planner exactly): budget
      !! against TOTAL RAM (deterministic tile size for a given cube/
      !! mem_frac_ram, reproducible across runs) rather than instantaneous
      !! free RAM. 4 GiB fallback if /proc/meminfo can't be read.
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

   subroutine plan_rmclean_tile()
      !! T4a (planning/RMCLEAN_INTEGRATION_PLAN.md): sets tile_ra/tile_dec,
      !! same RA-strips-first auto-tiling + safety-shrink policy as
      !! rm_synthesis_mod.f90's own plan_tile (rm_synthesis_mod.f90:
      !! 3200-3247), applied to rmclean_cubes' own per-tile-pixel byte
      !! budget: 2 input arrays (re,im), single-buffered, + 6 output
      !! arrays (clean/resid/restored re,im), each 4*nrm bytes/pixel --
      !! 14 array-widths total, NOT 8. Found as a real bug (a full-cube
      !! run's own cgroup memory reached 39.9G against a nominal
      !! mem_frac_ram=0.25 budget of ~15.7G, planning/
      !! RMCLEAN_INTEGRATION_PLAN.md T9's own UPDATE): this comment
      !! originally said the 6 output arrays' own double-buffering
      !! (io_read_threads/nwriters/io_overlap, T4d) would be
      !! "added on top of this same budget, not computed here" once
      !! implemented, but T4d actually allocates them permanently
      !! double-buffered (clean_re_buf(...,2) etc., below) REGARDLESS of
      !! whether io_overlap is even on -- the budget arithmetic was never
      !! updated to match, silently under-budgeting every mem_frac_ram=
      !! run by 1.75x (8 assumed vs 14 actual array-widths). RA (NAXIS1)
      !! is the contiguous axis on disk
      !! (FTGSVE/FTPSSE's own natural layout, same reasoning as
      !! reproject_cubes.f90's own tiling comment) -- keeping tile_ra=nx
      !! and packing Dec rows makes every tile read/write one contiguous
      !! run per RM-plane; only an extremely wide image (a single full-RA
      !! Dec row exceeding the budget) falls back to subdividing RA.
      integer(kind=8) :: mem_total_kb, bytes_per_tile_pixel, mem_safe_bytes
      integer(kind=8) :: tile_pixels_max, image_pixels_total, tile_bytes_est

      ! 2 input arrays (single-buffered) + 6 output arrays (each
      ! double-buffered, unconditionally -- see this subroutine's own
      ! top comment) = 2 + 6*2 = 14 array-widths, not 8, each 4 bytes/
      ! voxel x nrm voxels/pixel. T10b: the mask tile (read_mask_tile,
      ! replacing the old whole-cube-resident read_mask_cube) is now
      ! ALSO part of this same per-tile budget -- 1 byte/voxel x nchan
      ! voxels/pixel, single-buffered (read-only working data, never
      ! written out, so no double-buffering needed unlike the outputs).
      bytes_per_tile_pixel = int(4,8) * int(2 + 6*2,8) * int(nrm,8) +&
      &int(nchan,8)
      call get_mem_total_kb(mem_total_kb)
      mem_safe_bytes = int(real(mem_frac_ram,8) * real(mem_total_kb,8) *&
      &1024.0d0, 8)
      if (mem_safe_bytes.le.bytes_per_tile_pixel) mem_safe_bytes = bytes_per_tile_pixel
      tile_pixels_max = mem_safe_bytes / bytes_per_tile_pixel
      if (tile_pixels_max.lt.1_8) tile_pixels_max = 1_8
      image_pixels_total = int(nx,8) * int(ny,8)

      if (tile_auto .or. tile_ra.le.0 .or. tile_dec.le.0) then
         if (tile_pixels_max.ge.image_pixels_total) then
            tile_ra = nx
            tile_dec = ny
         else if (tile_pixels_max.ge.int(nx,8)) then
            tile_ra = nx
            tile_dec = int(tile_pixels_max / int(nx,8))
            if (tile_dec.lt.1) tile_dec = 1
            if (tile_dec.gt.ny) tile_dec = ny
         else
            tile_dec = 1
            tile_ra = int(tile_pixels_max)
            if (tile_ra.lt.1) tile_ra = 1
            if (tile_ra.gt.nx) tile_ra = nx
         endif
      else
         tile_ra = max(1, min(tile_ra, nx))
         tile_dec = max(1, min(tile_dec, ny))
      endif

      ! Safety shrink: reduce Dec rows first (keep full RA contiguous);
      ! only shrink RA once the strip is already a single Dec row -- same
      ! order as rm_synthesis_mod.f90's own plan_tile.
      tile_bytes_est = int(tile_ra,8) * int(tile_dec,8) * bytes_per_tile_pixel
      do while (tile_bytes_est.gt.mem_safe_bytes .and.&
      &(tile_ra.gt.1 .or. tile_dec.gt.1))
         if (tile_dec.gt.1) then
            tile_dec = max(1, tile_dec/2)
         else
            tile_ra = max(1, tile_ra/2)
         endif
         tile_bytes_est = int(tile_ra,8) * int(tile_dec,8) * bytes_per_tile_pixel
      enddo

      n_blocks_total = ((nx + tile_ra - 1) / tile_ra) *&
      &((ny + tile_dec - 1) / tile_dec)

      write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A)') 'Tile plan: tile_ra x tile_dec = ',&
      &tile_ra, ' x ', tile_dec, ' px (image ', nx, ' x ', ny,&
      &' px, mem_frac_ram budget) -- ', n_blocks_total, ' block(s) total.'
   end subroutine plan_rmclean_tile

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
      call safe_ftopen(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(filename)
         status = -1
         call free_fits_unit(unit)
         return
      endif

      call FTMNHD(unit, -1, 'CHANFREQ', 0, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: no CHANFREQ binary table in: ', trim(filename)
         call safe_ftclos(unit, fitsstat)
         status = -1
         return
      endif
      call FTGKYJ(unit, 'NAXIS2', nchan_out, comment, fitsstat)
      if (fitsstat.ne.0 .or. nchan_out.lt.1) then
         write(*,*) 'ERROR: bad or missing NAXIS2 on CHANFREQ table in: ',&
         &trim(filename)
         call safe_ftclos(unit, fitsstat)
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
         call safe_ftclos(unit, fitsstat)
         status = -1
         return
      endif
      call safe_ftclos(unit, fitsstat)
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
      call safe_ftopen(unit, trim(ampfile), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(ampfile)
         status = -1
         call free_fits_unit(unit)
         return
      endif
      call FTGKYJ(unit, 'NAXIS', naxis, comment, fitsstat)
      if (fitsstat.ne.0 .or. naxis.ne.3) then
         write(*,*) 'ERROR: expected NAXIS=3 in: ', trim(ampfile)
         call safe_ftclos(unit, fitsstat)
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
         call safe_ftclos(unit, fitsstat)
         status = -1
         return
      endif
      call safe_ftclos(unit, fitsstat)
      fitsstat = 0
      call safe_ftopen(unit, trim(phafile), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(phafile)
         status = -1
         call free_fits_unit(unit)
         return
      endif
      call FTGKYJ(unit, 'NAXIS1', nx2, comment, fitsstat)
      call FTGKYJ(unit, 'NAXIS2', ny2, comment, fitsstat)
      call FTGKYJ(unit, 'NAXIS3', nrm2, comment, fitsstat)
      call FTGKYD(unit, 'CDELT3', cdelt3_2, comment, fitsstat)
      call FTGKYD(unit, 'CRVAL3', crval3_2, comment, fitsstat)
      call safe_ftclos(unit, fitsstat)
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
      call safe_ftopen(unit, trim(filename), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         call free_fits_unit(unit)
         return
      endif
      fitsstat = 0
      call FTGKYD(unit, 'LSQREF', lsq_ref_dp, comment, fitsstat)
      if (fitsstat.eq.0) then
         found = .true.
         lsq_ref_out = real(lsq_ref_dp, sp)
      endif
      fitsstat = 0
      call safe_ftclos(unit, fitsstat)
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

   subroutine split_range_rmclean(n_total, n_threads, base, rem)
      !! Same convention as rm_synthesis_mod.f90's own
      !! split_channels_across_threads: n_threads threads each get `base`
      !! items, the first `rem` of them get one extra, covering n_total
      !! exactly with contiguous, disjoint ranges.
      integer, intent(in) :: n_total, n_threads
      integer, intent(out) :: base, rem
      base = n_total / n_threads
      rem = mod(n_total, n_threads)
   end subroutine split_range_rmclean

   subroutine read_amp_pha_tile(ampfile, phafile, ix0, iy0, tx, ty, re_out,&
   &im_out, status)
      !! T4a/T4b: reads one tile's AMP/PHA subregion (full RM depth,
      !! [ix0:ix0+tx-1, iy0:iy0+ty-1, 1:nrm]) via FTGSVE, converts to
      !! re/im -- same amp*cos(pha)/amp*sin(pha) forward convention the
      !! old whole-cube read_cube used (matching rm_synthesis's own
      !! p_tile_arr=sqrt(re^2+im^2)/phi_tile_arr=atan2(im,re) convention,
      !! rm_synthesis_mod.f90:1462-1465). Always-natural axis order (this
      !! program's own AMP/PHA cubes are always plain NAXIS=3 with
      !! axis1=RA,axis2=Dec,axis3=RM -- no reproject_cubes-style
      !! axis-order permutation needed).
      !!
      !! T4b: when io_read_threads_eff>1, the tile's own RM range is split
      !! (split_range_rmclean, same base/rem convention as rm_synthesis's
      !! own split_channels_across_threads) across that many threads, each
      !! opening its OWN readonly FTOPEN handle (a genuinely distinct unit
      !! number from fitsio_unit_mod's safe_ftopen) and reading its own
      !! disjoint RM-slice via FTGSVE -- safe because readonly CFITSIO
      !! opens on the same file are exempt from the same-file handle-
      !! aliasing hazard that makes concurrent READ-WRITE handles on one
      !! file unsafe (see write_output_tile's own T4c comment, and
      !! rm_synthesis's own io_read_threads, which this is a direct port
      !! of). Separately, and regardless of same-file-ness: the actual
      !! FTOPEN/FTCLOS calls themselves must go through safe_ftopen/
      !! safe_ftclos, not bare FTOPEN/FTCLOS with only FTGIOU/FTFIOU
      !! critical-guarded -- see read_amp_pha_chunk's own comment for the
      !! real SIGSEGV this caught.
      character(len=*), intent(in) :: ampfile, phafile
      integer, intent(in) :: ix0, iy0, tx, ty
      real(sp), intent(out) :: re_out(tx,ty,nrm), im_out(tx,ty,nrm)
      integer, intent(out) :: status
      integer :: group
      integer :: base_chan, rem_chan, ith, rm_beg, rm_len
      integer :: status_par

      status = 0
      status_par = 0
      group = 1
      call split_range_rmclean(nrm, io_read_threads_eff, base_chan, rem_chan)

      !$omp parallel do schedule(static) default(shared)&
      !$omp& private(ith,rm_beg,rm_len) num_threads(io_read_threads_eff)
      do ith = 0, io_read_threads_eff-1
         if (ith.lt.rem_chan) then
            rm_len = base_chan + 1
            rm_beg = ith*(base_chan+1) + 1
         else
            rm_len = base_chan
            rm_beg = rem_chan*(base_chan+1) + (ith-rem_chan)*base_chan + 1
         endif
         if (rm_len.gt.0) then
            call read_amp_pha_chunk(ampfile, phafile, ix0, iy0, tx, ty,&
            &rm_beg, rm_len, ith, re_out(:,:,rm_beg:rm_beg+rm_len-1),&
            &im_out(:,:,rm_beg:rm_beg+rm_len-1), status_par)
         endif
      enddo
      !$omp end parallel do
      if (status_par.ne.0) status = -1
   end subroutine read_amp_pha_tile

   subroutine read_amp_pha_chunk(ampfile, phafile, ix0, iy0, tx, ty,&
   &rm_beg, rm_len, thread_id, re_out, im_out, status_par)
      !! One io_read_threads worker's own disjoint RM-chunk of one tile --
      !! see read_amp_pha_tile's own comment. Each unit comes from
      !! fitsio_unit_mod's safe_ftopen (FTGIOU+FTOPEN inside one critical
      !! section, FTCLOS+FTFIOU likewise via safe_ftclos) -- NOT bare
      !! FTOPEN/FTCLOS with only the unit-number bookkeeping guarded.
      !! Found the hard way while porting this same fix to
      !! reproject_cubes.f90's own concurrently-called load_wcs: an
      !! isolated 8-thread reproducer doing concurrent FTGIOU+FTOPEN+
      !! (read)+FTCLOS+FTFIOU cycles crashed reliably (SIGSEGV) within a
      !! few hundred iterations unless FTOPEN/FTCLOS themselves share the
      !! SAME critical section as FTGIOU/FTFIOU -- their internal
      !! bookkeeping evidently shares mutable global state. This
      !! subroutine is the highest-concurrency CFITSIO call site in this
      !! program (up to io_read_threads_eff threads, every tile), so it
      !! would have been the first to hit that crash under real
      !! io_read_threads>1 use.
      character(len=*), intent(in) :: ampfile, phafile
      integer, intent(in) :: ix0, iy0, tx, ty, rm_beg, rm_len, thread_id
      real(sp), intent(out) :: re_out(tx,ty,rm_len), im_out(tx,ty,rm_len)
      integer, intent(inout) :: status_par
      integer :: unit, blocksize, fitsstat, group
      integer :: fpixel(3), lpixel(3), incs(3), naxes_full(3)
      logical :: anyflag
      real(sp), allocatable :: amp_chunk(:,:,:), pha_chunk(:,:,:)

      allocate(amp_chunk(tx,ty,rm_len), pha_chunk(tx,ty,rm_len))
      group = 1
      naxes_full = (/ nx, ny, nrm /)
      fpixel = (/ ix0, iy0, rm_beg /)
      lpixel = (/ ix0+tx-1, iy0+ty-1, rm_beg+rm_len-1 /)
      incs = (/ 1, 1, 1 /)

      fitsstat = 0
      call safe_ftopen(unit, trim(ampfile), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(ampfile)
         !$omp atomic write
         status_par = -1
         call free_fits_unit(unit)
         deallocate(amp_chunk, pha_chunk)
         return
      endif
      call FTGSVE(unit, group, 3, naxes_full, fpixel, lpixel, incs, 0.0_sp,&
      &amp_chunk, anyflag, fitsstat)
      call safe_ftclos(unit, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to read tile data from: ', trim(ampfile)
         !$omp atomic write
         status_par = -1
         deallocate(amp_chunk, pha_chunk)
         return
      endif

      fitsstat = 0
      call safe_ftopen(unit, trim(phafile), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(phafile)
         !$omp atomic write
         status_par = -1
         call free_fits_unit(unit)
         deallocate(amp_chunk, pha_chunk)
         return
      endif
      call FTGSVE(unit, group, 3, naxes_full, fpixel, lpixel, incs, 0.0_sp,&
      &pha_chunk, anyflag, fitsstat)
      call safe_ftclos(unit, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: failed to read tile data from: ', trim(phafile)
         !$omp atomic write
         status_par = -1
         deallocate(amp_chunk, pha_chunk)
         return
      endif

      re_out = amp_chunk*cos(pha_chunk)
      im_out = amp_chunk*sin(pha_chunk)
      deallocate(amp_chunk, pha_chunk)
   end subroutine read_amp_pha_chunk

   subroutine open_output_cube(template_file, outname, nx_in, ny_in, nrm_in,&
   &idx, status)
      !! T4a/T4c: creates outname and writes its header ONLY (no pixel
      !! data -- that now comes tile-by-tile via write_output_tile below,
      !! since the whole cube is never resident in memory at once).
      !! Header copied verbatim from template_file (always ampfile: same
      !! WCS on every axis, since this program never resamples anything --
      !! Gate 0 validates the existing RM grid rather than changing it).
      !!
      !! T4c: if nwriters_eff==1, leaves out_unit(idx) OPEN -- the
      !! caller keeps it open across the whole tile loop and closes it
      !! itself via close_output_cube once every tile has been written; a
      !! still-empty (never-written) FITS data segment is a well-defined
      !! intermediate state that later FTPSSE subset writes fill in. If
      !! nwriters_eff>1, fetches this HDU's pixel-data byte offset
      !! (FTGHAD) into out_datastart(idx) and closes the handle
      !! IMMEDIATELY -- this FTCLOS is what makes CFITSIO actually define/
      !! flush the HDU to its full declared NAXIS extent on disk (ported
      !! lesson from rm_synthesis.f90's own T6 postmortem: closing late
      !! risks CFITSIO flushing a stale internal buffer over raw-written
      !! bytes at ffclos time, since CFITSIO has no idea the raw writes
      !! happened; closing right after FTGHAD makes that race impossible
      !! rather than merely unlikely, and only after this close is the
      !! file's on-disk size guaranteed to already span every byte offset
      !! write_output_tile's raw-write path will ever compute).
      character(len=*), intent(in) :: template_file, outname
      integer, intent(in) :: nx_in, ny_in, nrm_in
      integer, intent(in) :: idx
      integer, intent(out) :: status
      integer :: src_unit, fitsstat, blocksize
      integer :: naxes_out(3)
      logical :: simple, extend
      integer(kind=8) :: headstart, dataend, local_datastart
      integer :: local_unit

      status = 0
      out_path(idx) = outname
      naxes_out(1) = nx_in
      naxes_out(2) = ny_in
      naxes_out(3) = nrm_in
      simple = .true.
      extend = .false.

      fitsstat = 0
      call safe_ftopen(src_unit, trim(template_file), 0, blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot open FITS file: ', trim(template_file)
         status = -1
         call free_fits_unit(src_unit)
         return
      endif

      ! safe_ftinit's FTINIT fails (nonzero fitsstat) if outname already
      ! exists on disk -- checked and bailed out on IMMEDIATELY, before
      ! any further CFITSIO call on out_unit(idx): every call after a
      ! failed FTINIT operates on a unit CFITSIO never actually set up,
      ! which is not a clean no-op but undefined behaviour (confirmed
      ! directly: this previously crashed with a SIGSEGV inside CFITSIO
      ! when an output file from an earlier run was left on disk --
      ! tests/run_tests.sh must therefore also clean up its own rmc_*
      ! outputs before each run, which it does).
      fitsstat = 0
      call safe_ftinit(out_unit(idx), trim(outname), blocksize, fitsstat)
      if (fitsstat.ne.0) then
         write(*,*) 'ERROR: cannot create (already exists?): ', trim(outname)
         status = -1
         call safe_ftclos(src_unit, fitsstat)
         call free_fits_unit(out_unit(idx))
         return
      endif
      call FTPHPR(out_unit(idx), simple, -32, 3, naxes_out, 0, 1, extend, fitsstat)
      call copy_generic_header_rmclean(src_unit, out_unit(idx), status)
      if (fitsstat.ne.0 .or. status.ne.0) then
         write(*,*) 'ERROR: failed to write header for: ', trim(outname)
         status = -1
         call safe_ftclos(src_unit, fitsstat)
         return
      endif
      call safe_ftclos(src_unit, fitsstat)
      out_is_open(idx) = .true.
      if (nwriters_eff.gt.1) then
         fitsstat = 0
         local_unit = out_unit(idx)
         ! Zero-init before the call: this system's installed libcfitsio
         ! FTGHAD writes only the LOWER 32 bits of its 3 output arguments
         ! (a genuine Fortran-wrapper/library ABI truncation, confirmed
         ! directly with a minimal standalone reproducer -- the upper 32
         ! bits of an integer(kind=8) receiving variable are left
         ! completely untouched by the call). Leaving them at whatever an
         ! uninitialized automatic variable happens to contain produced a
         ! real, intermittent (stack-content-dependent) corrupted
         ! out_datastart and a genuinely wrong CLEAN.AMP output --
         ! zero-initializing first makes the untouched upper bits always
         ! read back as zero, which is exactly correct for any offset
         ! that (like every FITS header offset in practice) fits in 32
         ! bits.
         headstart = 0_8
         local_datastart = 0_8
         dataend = 0_8
         call FTGHAD(local_unit, headstart, local_datastart, dataend, fitsstat)
         out_datastart(idx) = local_datastart
         if (fitsstat.ne.0) then
            write(*,*) 'ERROR: FTGHAD (data-start offset) failed for: ', trim(outname)
            status = -1
            call safe_ftclos(out_unit(idx), fitsstat)
            out_is_open(idx) = .false.
            return
         endif
         fitsstat = 0
         call safe_ftclos(out_unit(idx), fitsstat)
         out_is_open(idx) = .false.
      endif
   end subroutine open_output_cube

   subroutine write_output_tile(idx_amp, idx_pha, nx_in, ny_in, nrm_in,&
   &ix0, iy0, tx, ty, re_in, im_in, status)
      !! T4a/T4c: re/im -> amp/pha (inverse of read_amp_pha_tile's own
      !! forward convention). nwriters_eff==1 (default): written
      !! via FTPSSE (subset write) at this tile's own [ix0:ix0+tx-1,
      !! iy0:iy0+ty-1, 1:nrm_in] window into the already-open
      !! out_unit(idx_amp/idx_pha) (opened once by open_output_cube,
      !! closed once by close_output_cube after every tile has been
      !! written). nwriters_eff>1: this tile's own nrm_in range is
      !! split (split_range_rmclean) across that many threads, each
      !! writing its own disjoint RM-chunk via write_rm_chunk_raw_rmclean
      !! -- raw stream writes at computed byte offsets, bypassing
      !! CFITSIO's ftpsse for the pixel data entirely, since N CFITSIO
      !! handles opened read-write on the SAME file alias onto one shared
      !! internal buffer (cfitsio's own fits_already_open() contract) and
      !! concurrent writes through them corrupt it -- confirmed by
      !! rm_synthesis's own T4 postmortem (real SIGSEGV in memmove inside
      !! libcfitsio), which is why nwriters>1 never opens extra
      !! CFITSIO handles at all (see open_output_cube's own early-close).
      integer, intent(in) :: idx_amp, idx_pha
      integer, intent(in) :: nx_in, ny_in, nrm_in, ix0, iy0, tx, ty
      real(sp), intent(in) :: re_in(tx,ty,nrm_in), im_in(tx,ty,nrm_in)
      integer, intent(out) :: status
      integer :: fitsstat
      integer :: naxes_out(3), fpixel(3), lpixel(3)
      integer :: base_chan, rem_chan, ith, rm_beg, rm_len
      real(sp), allocatable :: amp_tile(:,:,:), pha_tile(:,:,:)

      status = 0
      allocate(amp_tile(tx,ty,nrm_in), pha_tile(tx,ty,nrm_in))
      amp_tile = sqrt(re_in**2 + im_in**2)
      pha_tile = atan2(im_in, re_in)

      if (nwriters_eff.le.1) then
         naxes_out = (/ nx_in, ny_in, nrm_in /)
         fpixel = (/ ix0, iy0, 1 /)
         lpixel = (/ ix0+tx-1, iy0+ty-1, nrm_in /)

         fitsstat = 0
         call FTPSSE(out_unit(idx_amp), 1, 3, naxes_out, fpixel, lpixel,&
         &amp_tile, fitsstat)
         if (fitsstat.ne.0) then
            write(*,*) 'ERROR: failed to write AMP tile at (', ix0, ',', iy0, ')'
            status = -1
         endif
         fitsstat = 0
         call FTPSSE(out_unit(idx_pha), 1, 3, naxes_out, fpixel, lpixel,&
         &pha_tile, fitsstat)
         if (fitsstat.ne.0) then
            write(*,*) 'ERROR: failed to write PHA tile at (', ix0, ',', iy0, ')'
            status = -1
         endif
      else
         call split_range_rmclean(nrm_in, nwriters_eff, base_chan, rem_chan)
         !$omp parallel do schedule(static) default(shared)&
         !$omp& private(ith,rm_beg,rm_len) num_threads(nwriters_eff)
         do ith = 0, nwriters_eff-1
            if (ith.lt.rem_chan) then
               rm_len = base_chan + 1
               rm_beg = ith*(base_chan+1) + 1
            else
               rm_len = base_chan
               rm_beg = rem_chan*(base_chan+1) + (ith-rem_chan)*base_chan + 1
            endif
            if (rm_len.gt.0) then
               call write_rm_chunk_raw_rmclean(out_path(idx_amp),&
               &out_datastart(idx_amp), nx_in, ny_in, ix0, ix0+tx-1,&
               &iy0, iy0+ty-1, rm_beg, rm_beg+rm_len-1,&
               &amp_tile(:,:,rm_beg:rm_beg+rm_len-1))
               call write_rm_chunk_raw_rmclean(out_path(idx_pha),&
               &out_datastart(idx_pha), nx_in, ny_in, ix0, ix0+tx-1,&
               &iy0, iy0+ty-1, rm_beg, rm_beg+rm_len-1,&
               &pha_tile(:,:,rm_beg:rm_beg+rm_len-1))
            endif
         enddo
         !$omp end parallel do
      endif
      deallocate(amp_tile, pha_tile)
   end subroutine write_output_tile

   subroutine do_tile_write(job)
      !! T4d: the actual "write one tile's 6 output cubes" payload --
      !! callable either inline (io_overlap=n) or as a pthread entry
      !! point's own payload (io_overlap=y); identical logic either way,
      !! so output is bit-for-bit the same regardless of which mode
      !! dispatched it (same invariant as rm_synthesis_mod.f90's own
      !! do_tile_write). Just calls the already-existing, already-tested
      !! write_output_tile (T4a/T4c) three times -- once per (clean,
      !! resid, restored) pair -- using this job's own tile geometry and
      !! the 6 output-array pointers it was populated with. Any write
      !! failure is reported (write_output_tile's own error prints) but
      !! not propagated as a return status: once dispatched onto a
      !! background thread, there is no synchronous caller left to hand a
      !! status to -- same convention rm_synthesis's own do_tile_write
      !! uses (log the error, don't silently drop the write, don't try to
      !! halt the program from a background thread).
      type(tile_write_job_t), intent(inout) :: job
      integer :: status_local

      call write_output_tile(idx_clean_amp, idx_clean_pha, nx, ny, nrm,&
      &job%ix_tile_beg, job%iy_tile_beg, job%tx, job%ty,&
      &job%re_ptr(1)%p, job%im_ptr(1)%p, status_local)
      call write_output_tile(idx_resid_amp, idx_resid_pha, nx, ny, nrm,&
      &job%ix_tile_beg, job%iy_tile_beg, job%tx, job%ty,&
      &job%re_ptr(2)%p, job%im_ptr(2)%p, status_local)
      call write_output_tile(idx_restored_amp, idx_restored_pha, nx, ny, nrm,&
      &job%ix_tile_beg, job%iy_tile_beg, job%tx, job%ty,&
      &job%re_ptr(3)%p, job%im_ptr(3)%p, status_local)
   end subroutine do_tile_write

   function tile_write_thread_entry(arg) bind(C) result(res)
      !! pthread start routine. Unpacks the opaque context pointer back
      !! into the tile_write_job_t it was created from (same process,
      !! same build, so the round-trip through c_loc/c_f_pointer is safe
      !! even though tile_write_job_t is not a bind(C) type) and runs the
      !! write. Verbatim port of rm_synthesis_mod.f90's own
      !! tile_write_thread_entry.
      type(c_ptr), value :: arg
      type(c_ptr) :: res
      type(tile_write_job_t), pointer :: job

      call c_f_pointer(arg, job)
      call do_tile_write(job)
      res = c_null_ptr
   end function tile_write_thread_entry

   subroutine tile_write_dispatch_async(job, thread_id, dispatched)
      !! Launches do_tile_write(job) on a background pthread. `job` must
      !! have the TARGET attribute at the call site (write_job(0:1)'s own
      !! declaration) and must remain valid -- untouched and
      !! undeallocated -- until tile_write_join(thread_id) has returned.
      !! On pthread_create failure this runs the write synchronously right
      !! here instead (safe fallback: a write is never silently dropped),
      !! and reports dispatched=.false. so the caller knows there is
      !! nothing to join later. Verbatim port of rm_synthesis_mod.f90's
      !! own tile_write_dispatch_async.
      type(tile_write_job_t), intent(inout), target :: job
      integer(c_long), intent(out) :: thread_id
      logical, intent(out) :: dispatched
      integer(c_int) :: rc

      rc = c_pthread_create(thread_id, c_null_ptr,&
      &c_funloc(tile_write_thread_entry), c_loc(job))
      if (rc.ne.0_c_int) then
         write(*,*) 'WARNING: pthread_create failed for async tile write;'//&
         &' running inline'
         call do_tile_write(job)
         dispatched = .false.
      else
         dispatched = .true.
      endif
   end subroutine tile_write_dispatch_async

   subroutine tile_write_join(thread_id)
      integer(c_long), intent(in) :: thread_id
      integer(c_int) :: rc
      rc = c_pthread_join(thread_id, c_null_ptr)
   end subroutine tile_write_join

   subroutine close_output_cube(idx)
      integer, intent(in) :: idx
      integer :: fitsstat
      if (.not. out_is_open(idx)) return
      fitsstat = 0
      call safe_ftclos(out_unit(idx), fitsstat)
      out_is_open(idx) = .false.
   end subroutine close_output_cube

   logical function host_is_big_endian_rmclean() result(is_be)
      !! Verbatim port of rm_synthesis_mod.f90's own host_is_big_endian --
      !! see there for the full rationale (FITS mandates big-endian; every
      !! realistic deployment target is little-endian, so this is a
      !! runtime check rather than a hard-coded assumption to stay
      !! correct, not just fast, on a big-endian host).
      integer(kind=4) :: probe
      integer(kind=1) :: bytes(4)
      probe = 1_4
      bytes = transfer(probe, bytes)
      is_be = (bytes(1).eq.0_1)
   end function host_is_big_endian_rmclean

   subroutine swap_bytes_r4_inplace_rmclean(buf, n)
      !! Verbatim port of rm_synthesis_mod.f90's own swap_bytes_r4_inplace.
      integer(kind=8), intent(in) :: n
      real(sp), intent(inout) :: buf(n)
      integer(kind=8) :: i
      integer(kind=1) :: b(4), t
      do i = 1_8, n
         b = transfer(buf(i), b)
         t = b(1); b(1) = b(4); b(4) = t
         t = b(2); b(2) = b(3); b(3) = t
         buf(i) = transfer(b, buf(i))
      enddo
   end subroutine swap_bytes_r4_inplace_rmclean

   subroutine write_rm_chunk_raw_rmclean(file_path, datastart, nx_out,&
   &ny_out, ix_out_beg, ix_out_end, iy_out_beg, iy_out_end, rm_beg,&
   &rm_end, data)
      !! Adaptation of rm_synthesis_mod.f90's own write_rm_chunk_raw --
      !! see there for the full byte-offset-math/two-write-pattern/
      !! endianness rationale. Writes one contiguous RM-bin range
      !! [rm_beg,rm_end] of one tile's output data directly to disk via
      !! plain Fortran stream I/O, bypassing CFITSIO's ftpsse for the
      !! pixel data entirely (see write_output_tile's own comment for
      !! why).
      !!
      !! This routine is reached concurrently (once per nwriters
      !! worker, each opening the SAME file path for its own disjoint
      !! byte range) from write_output_tile's own `!$omp parallel do`.
      !! The Fortran standard does not guarantee any I/O statement --
      !! open() included, newunit= or not -- is safe to call concurrently
      !! without explicit synchronization; a first attempt at fixing a
      !! real, observed corruption here used fixed, pre-assigned
      !! per-thread unit numbers instead, which fixed it but only via
      !! manual, cross-file bookkeeping (every other unit-number range in
      !! this codebase has to stay disjoint from it by inspection, not
      !! anything the compiler/runtime enforces). The named critical
      !! section below makes unit allocation itself atomic instead:
      !! genuinely unique (real newunit= semantics, no manual range to
      !! maintain), while the write itself below still runs fully in
      !! parallel per thread -- only the brief "grab a unit" step is
      !! serialized.
      character(len=*), intent(in) :: file_path
      integer(kind=8), intent(in) :: datastart
      integer, intent(in) :: nx_out, ny_out
      integer, intent(in) :: ix_out_beg, ix_out_end, iy_out_beg, iy_out_end
      integer, intent(in) :: rm_beg, rm_end
      real(sp), intent(in) :: data(*)

      integer :: u, ios
      logical :: full_width, need_swap
      integer(kind=8) :: row_len, n_rows, plane_elems, plane_stride
      integer(kind=8) :: irm, iy_local, mem_off, byte_pos
      real(sp), allocatable :: plane_buf(:)

      full_width = (ix_out_beg.eq.1 .and. ix_out_end.eq.nx_out)
      need_swap = .not. host_is_big_endian_rmclean()
      row_len = int(ix_out_end-ix_out_beg+1, 8)
      n_rows = int(iy_out_end-iy_out_beg+1, 8)
      plane_elems = row_len * n_rows
      plane_stride = int(nx_out, 8) * int(ny_out, 8)

      !$omp critical (raw_write_open_lock)
      open(newunit=u, file=trim(file_path), access='stream',&
      &form='unformatted', status='old', action='write', iostat=ios)
      !$omp end critical (raw_write_open_lock)
      if (ios.ne.0) then
         write(*,*) 'ERROR: write_rm_chunk_raw_rmclean: failed to open ',&
         &trim(file_path)
         return
      endif

      allocate(plane_buf(plane_elems))
      mem_off = 1_8
      do irm = int(rm_beg,8), int(rm_end,8)
         plane_buf = data(mem_off:mem_off+plane_elems-1_8)
         if (need_swap) call swap_bytes_r4_inplace_rmclean(plane_buf, plane_elems)

         if (full_width) then
            byte_pos = datastart + (irm-1_8)*plane_stride*4_8&
            &+ int(iy_out_beg-1,8)*int(nx_out,8)*4_8 + 1_8
            write(u, pos=byte_pos, iostat=ios) plane_buf
            if (ios.ne.0) write(*,*) 'ERROR: write_rm_chunk_raw_rmclean:'//&
            &' write failed (full-width) for ', trim(file_path)
         else
            do iy_local = 0_8, n_rows-1_8
               byte_pos = datastart + (irm-1_8)*plane_stride*4_8&
               &+ (int(iy_out_beg,8)-1_8+iy_local)*int(nx_out,8)*4_8&
               &+ int(ix_out_beg-1,8)*4_8 + 1_8
               write(u, pos=byte_pos, iostat=ios)&
               &plane_buf(iy_local*row_len+1_8:iy_local*row_len+row_len)
               if (ios.ne.0) write(*,*) 'ERROR: write_rm_chunk_raw_rmclean:'//&
               &' write failed (row) for ', trim(file_path)
            enddo
         endif
         mem_off = mem_off + plane_elems
      enddo
      deallocate(plane_buf)
      close(u)
   end subroutine write_rm_chunk_raw_rmclean

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

   function cached_pattern_matches(cidx, pattern) result(matches)
      !! T14 increment 9: whether cache slot cidx's own pattern equals
      !! pattern -- the single place both cache-lookup sites
      !! (cache_lookup_readonly and update_mask_pattern_cache_for_tile's
      !! own inline lookup) go to compare a cached entry's pattern,
      !! since WHERE that pattern's bytes actually live depends on
      !! cache_eviction_policy (see table_cache_entry_t's own pattern(:)
      !! comment): 'hitcount' reads cache_entries(cidx)%pattern
      !! directly (the only place those bytes are stored); 'belady'
      !! reads pattern_registry%entries(cache_entries(cidx)%registry_id)
      !! %pattern instead -- the SAME bytes already live there (Pass 0
      !! put them there), so the cache slot never duplicates them.
      !! PURE READ access (cache_eviction_policy is fixed for the whole
      !! run by the time this is ever called, and pattern_registry is
      !! never mutated during the parallel CLEAN loop), safe to call
      !! concurrently, same as cache_lookup_readonly itself.
      integer, intent(in) :: cidx
      integer(kind=1), intent(in) :: pattern(:)
      logical :: matches
      integer :: rid
      if (cache_eviction_policy.eq.'belady') then
         rid = cache_entries(cidx)%registry_id
         matches = size(pattern_registry%entries(rid)%pattern).eq.size(pattern)
         if (matches) matches = all(pattern_registry%entries(rid)%pattern.eq.pattern)
      else
         matches = size(cache_entries(cidx)%pattern).eq.size(pattern)
         if (matches) matches = all(cache_entries(cidx)%pattern.eq.pattern)
      endif
   end function cached_pattern_matches

   subroutine cache_lookup_readonly(pattern, entry_idx)
      !! Open-addressing (linear probing) lookup into cache_buckets/
      !! cache_entries -- PURE READ access, safe to call concurrently
      !! from every OMP thread in the main CLEAN loop, since T10b's own
      !! update_mask_pattern_cache_for_tile has already fully populated
      !! the cache for THIS tile serially before that tile's own parallel
      !! region ever starts (no insertion ever happens in here).
      !! entry_idx=0 means "not cached" -- either this exact pattern
      !! never occurred during the pre-scan (impossible, since the
      !! pre-scan visits the same mask_tile this is called against) or
      !! the cache was already at mask_pattern_cache_max when this
      !! pattern was first seen --
      !! either way, the caller's own correct response is to build a
      !! one-off throwaway table (see clean_one_pixel).
      !! T14: cache_buckets(slot)==-1 is a TOMBSTONE (a bucket slot
      !! whose own cache_entries row has since been evicted and
      !! overwritten by a DIFFERENT pattern, cache_eviction_policy=
      !! 'hitcount'/'belady') -- standard open-addressing-with-deletion:
      !! a tombstone must NOT stop the probe (the pattern being looked
      !! up may still exist further along the same probe sequence,
      !! exactly as if this were any other non-matching occupied slot),
      !! and must never be dereferenced as a cache_entries index (-1 is
      !! not a valid one). Only a genuinely empty slot (0) proves "not
      !! cached, stop looking."
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
         if (cache_buckets(slot).ne.-1) then
            if (cache_entries(cache_buckets(slot))%hash.eq.h) then
               if (cached_pattern_matches(cache_buckets(slot), pattern)) then
                  entry_idx = cache_buckets(slot)
                  return
               endif
            endif
         endif
         probe = modulo(probe+1, int(n_cache_buckets, 8))
      enddo
   end subroutine cache_lookup_readonly

   subroutine init_mask_pattern_cache()
      !! T10b: allocates the (empty) cache structures once, before the
      !! tile loop starts -- no scanning here any more (that now happens
      !! per tile, see update_mask_pattern_cache_for_tile below), since
      !! the mask is no longer whole-cube-resident at this point in the
      !! program. Running totals for the final summary print (moved to
      !! after the whole tile loop -- there is no longer one single
      !! pre-scan moment to print it at) start at zero here.
      n_cache_buckets = max(16, 4*mask_pattern_cache_max)
      allocate(cache_buckets(0:n_cache_buckets-1))
      cache_buckets = 0
      allocate(cache_entries(max(1, mask_pattern_cache_max)))
      n_cache_entries = 0
      n_distinct_patterns_total = 0
      n_overflow_pixels_total = 0
   end subroutine init_mask_pattern_cache

   subroutine run_pattern_prescan()
      !! T14 (docs/dev/RMCLEAN_INTEGRATION_PLAN.md): Pass 0 -- scans the
      !! WHOLE mask cube once, before real CLEAN compute starts,
      !! building pattern_registry's own complete future-occurrence
      !! timeline for every distinct channel-validity pattern. Called
      !! only when cache_eviction_policy='belady' (see call site) --
      !! 'hitcount' needs no registry and skips this entirely, no point
      !! paying an extra full mask read for a registry that won't be
      !! used.
      !!
      !! Uses next_tile_extent (rmclean_io_mod) for tile geometry --
      !! the SAME function the real Pass-1 tile loop uses, not a
      !! hand-duplicated rule -- so the two passes can never silently
      !! diverge on scan order, the one property Pass 0's own
      !! scan-position bookkeeping absolutely depends on. mask_tile
      !! itself is the SAME array the real Pass-1 loop will use
      !! (already allocated by the caller, before this call).
      !!
      !! Safety valve: if the number of distinct patterns found so far
      !! exceeds 10% of the total image pixel count (nx*ny), a
      !! lookahead-based cache isn't worth its own memory/IO cost for
      !! this dataset (too little pattern reuse to exploit) -- abort
      !! immediately, warn loudly, and fall back to
      !! cache_eviction_policy='hitcount' for the rest of this run,
      !! rather than finish a scan whose own output wouldn't be worth
      !! keeping anyway.
      integer :: ix_tile_beg, iy_tile_beg, ix_cur, iy_cur, tx, ty
      logical :: done
      integer :: ix_l, iy_l, entry_id
      integer(kind=8) :: scan_pos, safety_threshold
      integer :: status_prescan

      write(*,'(A)') 'Pass 0: pre-scanning the whole mask cube to'//&
      &' build the pattern registry (cache_eviction_policy=belady)...'

      call registry_init(pattern_registry)
      safety_threshold = (int(nx, 8) * int(ny, 8)) / 10_8

      scan_pos = 0_8
      ix_tile_beg = 1
      iy_tile_beg = 1
      do
         ix_cur = ix_tile_beg
         iy_cur = iy_tile_beg
         call next_tile_extent(nx, ny, tile_ra, tile_dec, ix_tile_beg,&
         &iy_tile_beg, tx, ty, done)
         if (done) exit

         call read_mask_tile(maskfile, nx, ny, nchan, ix_cur, iy_cur, tx, ty,&
         &io_read_threads_eff, mask_tile(1:tx,1:ty,:), status_prescan)
         if (status_prescan.ne.0) then
            write(*,*) 'ERROR: Pass 0 failed to read mask tile at (',&
            &ix_cur, ',', iy_cur, ')'
            stop 1
         endif

         do iy_l = 1, ty
            do ix_l = 1, tx
               if (all(mask_tile(ix_l,iy_l,:).eq.0_1)) cycle
               ! T17: same threshold, same skip, as
               ! update_mask_pattern_cache_for_tile/clean_one_pixel --
               ! a pixel Pass 1 will never CLEAN must never be
               ! registered here either, or Pass 0's own scan-position
               ! bookkeeping would silently diverge from Pass 1's own
               ! actual visit order (the one property T14's own Belady
               ! eviction decisions absolutely depend on).
               if (min_valid_chan_frac.gt.0.0_sp .and.&
               &real(count(mask_tile(ix_l,iy_l,:).ne.0_1), sp).lt.&
               &min_valid_chan_frac*real(nchan, sp)) cycle
               scan_pos = scan_pos + 1_8
               call registry_lookup_or_insert(pattern_registry,&
               &mask_tile(ix_l,iy_l,:), nchan, scan_pos, entry_id)
            enddo
         enddo

         if (int(pattern_registry%n_entries, 8).ge.safety_threshold) then
            write(*,'(A,I0,A)') 'WARNING: Pass 0 found ',&
            &pattern_registry%n_entries, ' distinct patterns, already'//&
            &' past 10% of the total image pixel count -- too little'//&
            &' pattern reuse for a lookahead-based cache to be worth'//&
            &' its own memory/IO cost. Aborting Pass 0 and falling'//&
            &' back to cache_eviction_policy=hitcount for the rest of'//&
            &' this run.'
            cache_eviction_policy = 'hitcount'
            call registry_init(pattern_registry) ! discard what was built, free its memory
            return
         endif
      end do

      write(*,'(A,I0,A,I0,A)') 'Pass 0: done -- ',&
      &pattern_registry%n_entries, ' distinct pattern(s), ', scan_pos,&
      &' valid pixel(s) scanned.'

      call report_build_time_advisory()
   end subroutine run_pattern_prescan

   subroutine report_build_time_advisory()
      !! T15 (docs/dev/RMCLEAN_INTEGRATION_PLAN.md): predicts, per
      !! block and in total, how much wall-clock time Pass 1 will spend
      !! GENERATING new RMSF tables (build_rmsf_offset_table) -- never
      !! CLEAN-iteration time, which depends on each pixel's own real
      !! SNR/stopping behaviour, not on caching, and would answer a
      !! different question than the one this advisory exists to
      !! answer (how much of the old no-eviction disaster T14 fixes is
      !! actually recoverable). Printed to stdout only -- already
      !! captured in the run's own provenance log by the caller's
      !! existing `> logfile 2>&1` redirection, same convention every
      !! other Pass 0/startup message in this program already relies
      !! on; no separate logging mechanism needed.
      !!
      !! A single, cheap, read-only second scan of the mask cube
      !! (registry_lookup, not insert -- the registry is already
      !! complete) does three things together:
      !!
      !! (1) A full BELADY SIMULATION at the ACTUAL configured
      !! mask_pattern_cache_max: replays the EXACT same admission-
      !! control eviction logic update_mask_pattern_cache_for_tile's
      !! own 'belady' branch uses in real Pass 1, against a throwaway
      !! simulated cache (integer bookkeeping only, no RMSF tables
      !! built). This is NOT optional detail -- an earlier version of
      !! this advisory counted only each pattern's first-ever
      !! occurrence in the whole cube (i.e. assumed an INFINITELY large
      !! cache) and undercounted block 3's own real table-generation
      !! time by ~7x purely because the real 4096-slot cache is far
      !! smaller than that block's own local diversity, forcing genuine
      !! repeated rebuilds that simpler model couldn't see.
      !!
      !! (2) A CACHE-SIZE SWEEP: the SAME simulation, run in parallel
      !! (same single scan, same shared per-pattern occurrence data) at
      !! several LARGER candidate cache sizes (doubling from the actual
      !! configured size up to n_entries -- a cache that size can hold
      !! every distinct pattern simultaneously, the best any cache size
      !! could ever do), reusing each block's own already-measured
      !! ms/channel rate (a property of that block's own DATA, not of
      !! cache policy) rather than re-timing per candidate. Answers a
      !! real, asked-for question directly from data already collected:
      !! how much memory would it take to meaningfully cut this
      !! predicted build time -- without ever building a second,
      !! separately-managed cache (mathematically pointless: Belady is
      !! already optimal for whatever total capacity it's given, so two
      !! independently-managed pools of the same combined size can
      !! never beat one unified cache of that size -- see this ticket's
      !! own progress notes for the full reasoning).
      !!
      !! (3) A TRUE-SINGLETON vs RECURS-BUT-DECLINED split for the
      !! actual cache size's own declined-oneoff pixels: a decline
      !! means "not worth evicting anything for right now," which
      !! covers two different situations -- the pattern's own next
      !! occurrence is huge() (never again, genuinely unique data no
      !! cache size could ever help) versus a real, finite future
      !! position (it WILL recur, just not soon enough to beat whatever
      !! else is cached at that moment -- exactly the case a bigger
      !! cache, not a smarter policy, could capture). If declines are
      !! mostly the latter, the cache-size sweep above is the
      !! actionable lever; if mostly the former, no cache size will
      !! help much.
      !!
      !! All three reuse the SAME per-pattern next-occurrence walk
      !! (sim_next_occ_ptr, SEPARATE from pattern_registry's own
      !! entries(:)%next_occ_ptr, which must stay at its untouched,
      !! post-Pass-0 state for Pass 1's own real use) -- every
      !! occurrence of a pattern advances this pointer exactly once,
      !! regardless of which candidate cache size is being evaluated,
      !! so it is computed once per pixel and shared across all
      !! candidates rather than duplicated.
      !!
      !! Cost: a second full mask-cube read (I/O, comparable to Pass
      !! 0's own first read) plus, per miss, a linear scan bounded by
      !! each candidate's own cache size (exactly the same per-miss
      !! cost real Pass 1 itself pays for its own real eviction
      !! decisions, just repeated for a handful of candidates) -- no
      !! RMSF tables are built except the actual cache size's own
      !! per-block timing samples, so this stays far cheaper than
      !! actual Pass 1 compute.
      integer, parameter :: n_sample_per_block = 5
      integer, parameter :: max_candidates = 8
      integer, allocatable :: pattern_nvalid(:)
      integer(kind=8), allocatable :: sim_next_occ_ptr(:)
      integer :: n_entries_tot

      integer :: n_candidates
      integer :: candidate_cache_max(max_candidates)
      integer, allocatable :: sim_pattern_slot_c(:,:), sim_slot_pattern_c(:,:)
      integer :: sim_n_resident_c(max_candidates)
      integer(kind=8), allocatable :: victim_next_occs_sim(:)
      integer(kind=8) :: cand_hits(max_candidates), cand_admitted(max_candidates),&
      &cand_evicted(max_candidates), cand_declined(max_candidates)
      integer(kind=8) :: cand_declined_singleton, cand_declined_recurs
      real(dp) :: cand_predicted_total_min(max_candidates)
      integer(kind=8) :: cand_block_miss_nvalid_sum(max_candidates)

      integer :: ix_tile_beg, iy_tile_beg, ix_cur, iy_cur, tx, ty
      logical :: done
      integer :: ix_l, iy_l, entry_id, target_slot
      integer :: status_prescan
      integer(kind=8) :: new_next_occ
      integer :: n_blocks, block_cap
      integer, allocatable :: block_hits(:), block_admitted(:),&
      &block_evicted(:), block_declined(:), block_declined_singleton(:),&
      &block_declined_recurs(:)
      integer(kind=8), allocatable :: block_miss_nvalid_sum(:)
      real(dp), allocatable :: block_ms_per_channel(:)
      integer, allocatable :: tmp_i(:)
      integer(kind=8), allocatable :: tmp_i8(:)
      real(dp), allocatable :: tmp_r8(:)
      integer(kind=8) :: tot_hits, tot_admitted, tot_evicted, tot_declined,&
      &tot_declined_singleton, tot_declined_recurs
      real(dp) :: predicted_block_min, predicted_total_min
      integer :: b, k, k2, kc

      integer :: block_sample_ids(n_sample_per_block)
      integer :: block_sample_filled, block_sample_seen
      real(sp) :: rnd
      integer :: rnd_slot
      integer, allocatable :: valid_idx_s(:)
      real(sp), allocatable :: l_sq_valid_s(:)
      type(rmsf_table_t) :: table_s
      integer(kind=8) :: c0, c1, crate
      real(dp) :: sample_total_ms
      integer(kind=8) :: sample_total_channels
      integer :: nvalid_s
      integer :: saved_max_threads
      real(dp) :: bytes_per_table
      logical :: have_bytes_per_table

      ! T15 follow-up: <outfile>.advisory.csv, one row per block, so a
      ! standalone tool (scripts/plot_rmclean_advisory.py) can plot the
      ! predicted table-GENERATION time per block without re-parsing
      ! this program's own stdout text. Kept as a plain CSV rather than
      ! auto-invoking a plotting script from here -- this binary has no
      ! Python/matplotlib dependency anywhere else, and shelling out for
      ! a diagnostic plot is a separate, explicitly deferred decision
      ! (docs/dev/RMCLEAN_INTEGRATION_PLAN.md T15).
      integer :: csv_unit, ios_csv
      character(len=600) :: csv_path

      if (pattern_registry%n_entries.lt.1) return
      n_entries_tot = pattern_registry%n_entries

      ! Candidate cache sizes: the actual configured size, then
      ! doubling, capped at n_entries_tot (a cache that big can hold
      ! every distinct pattern at once -- the best ANY size could do,
      ! a natural ceiling data point). max_candidates=8 is generous
      ! headroom; real datasets need far fewer doublings to reach the
      ! ceiling (e.g. ~4096 -> ~62000 takes 5 steps).
      n_candidates = 0
      k = mask_pattern_cache_max
      do while (n_candidates.lt.max_candidates)
         n_candidates = n_candidates + 1
         candidate_cache_max(n_candidates) = min(k, n_entries_tot)
         if (k.ge.n_entries_tot) exit
         k = k * 2
      enddo

      ! build_rmsf_offset_table has its own internal !$omp parallel do
      ! (src/rmclean.f90:439). In REAL production use it's called from
      ! clean_one_pixel, itself already inside the outer per-pixel
      ! !$omp parallel do -- OMP nesting is off by default, so that
      ! inner parallel region is a no-op there, and the call runs
      ! serially on whichever single thread already owns that pixel.
      ! This whole simulation is called from run_pattern_prescan, which
      ! is NOT nested inside any other parallel region -- without
      ! forcing single-threaded execution for the timing samples below,
      ! the inner parallel do would get REAL multi-thread speedup
      ! production never sees, making every measurement artificially
      ! fast (confirmed directly: an earlier version measured 0.062
      ! ms/channel with 16 threads available, versus 0.222 ms/channel
      ! from a genuinely single-threaded standalone microbenchmark of
      ! the same call -- a ~3.6x discrepancy, not a rounding
      ! difference). Set once for this whole subroutine, restored once
      ! at the end.
      saved_max_threads = omp_get_max_threads()
      call omp_set_num_threads(1)

      allocate(pattern_nvalid(n_entries_tot))
      do k = 1, n_entries_tot
         pattern_nvalid(k) = count(pattern_registry%entries(k)%pattern.ne.0_1)
      enddo
      allocate(sim_next_occ_ptr(n_entries_tot))
      sim_next_occ_ptr = 1_8
      allocate(sim_pattern_slot_c(n_entries_tot, n_candidates))
      sim_pattern_slot_c = 0
      allocate(sim_slot_pattern_c(n_entries_tot, n_candidates)) ! upper bound: largest candidate
      sim_slot_pattern_c = 0
      sim_n_resident_c(1:n_candidates) = 0
      allocate(victim_next_occs_sim(n_entries_tot)) ! upper bound: largest candidate
      allocate(valid_idx_s(nchan), l_sq_valid_s(nchan))
      cand_hits = 0_8; cand_admitted = 0_8; cand_evicted = 0_8; cand_declined = 0_8
      cand_predicted_total_min = 0.0_dp
      cand_declined_singleton = 0_8
      cand_declined_recurs = 0_8
      have_bytes_per_table = .false.
      bytes_per_table = 0.0_dp

      block_cap = 64
      allocate(block_hits(block_cap), block_admitted(block_cap),&
      &block_evicted(block_cap), block_declined(block_cap),&
      &block_declined_singleton(block_cap), block_declined_recurs(block_cap),&
      &block_miss_nvalid_sum(block_cap), block_ms_per_channel(block_cap))
      n_blocks = 0

      ix_tile_beg = 1
      iy_tile_beg = 1
      do
         ix_cur = ix_tile_beg
         iy_cur = iy_tile_beg
         call next_tile_extent(nx, ny, tile_ra, tile_dec, ix_tile_beg,&
         &iy_tile_beg, tx, ty, done)
         if (done) exit

         call read_mask_tile(maskfile, nx, ny, nchan, ix_cur, iy_cur, tx, ty,&
         &io_read_threads_eff, mask_tile(1:tx,1:ty,:), status_prescan)
         if (status_prescan.ne.0) then
            write(*,*) 'WARNING: T15 Belady simulation failed to'//&
            &' re-read a mask tile -- skipping the build-time advisory.'
            call omp_set_num_threads(saved_max_threads)
            return
         endif

         n_blocks = n_blocks + 1
         if (n_blocks.gt.block_cap) then
            block_cap = block_cap * 2
            allocate(tmp_i(block_cap)); tmp_i(1:n_blocks-1)=block_hits
            call move_alloc(tmp_i, block_hits)
            allocate(tmp_i(block_cap)); tmp_i(1:n_blocks-1)=block_admitted
            call move_alloc(tmp_i, block_admitted)
            allocate(tmp_i(block_cap)); tmp_i(1:n_blocks-1)=block_evicted
            call move_alloc(tmp_i, block_evicted)
            allocate(tmp_i(block_cap)); tmp_i(1:n_blocks-1)=block_declined
            call move_alloc(tmp_i, block_declined)
            allocate(tmp_i(block_cap)); tmp_i(1:n_blocks-1)=block_declined_singleton
            call move_alloc(tmp_i, block_declined_singleton)
            allocate(tmp_i(block_cap)); tmp_i(1:n_blocks-1)=block_declined_recurs
            call move_alloc(tmp_i, block_declined_recurs)
            allocate(tmp_i8(block_cap)); tmp_i8(1:n_blocks-1)=block_miss_nvalid_sum
            call move_alloc(tmp_i8, block_miss_nvalid_sum)
            allocate(tmp_r8(block_cap)); tmp_r8(1:n_blocks-1)=block_ms_per_channel
            call move_alloc(tmp_r8, block_ms_per_channel)
         endif
         block_hits(n_blocks) = 0
         block_admitted(n_blocks) = 0
         block_evicted(n_blocks) = 0
         block_declined(n_blocks) = 0
         block_declined_singleton(n_blocks) = 0
         block_declined_recurs(n_blocks) = 0
         block_miss_nvalid_sum(n_blocks) = 0_8
         block_sample_filled = 0
         block_sample_seen = 0
         cand_block_miss_nvalid_sum(1:n_candidates) = 0_8

         do iy_l = 1, ty
            do ix_l = 1, tx
               if (all(mask_tile(ix_l,iy_l,:).eq.0_1)) cycle
               ! T17: same threshold as run_pattern_prescan/Pass 1 --
               ! a below-threshold pixel was never registered, so
               ! looking it up here would hit the "not found" case
               ! below for a reason that has nothing to do with a real
               ! bug.
               if (min_valid_chan_frac.gt.0.0_sp .and.&
               &real(count(mask_tile(ix_l,iy_l,:).ne.0_1), sp).lt.&
               &min_valid_chan_frac*real(nchan, sp)) cycle
               call registry_lookup(pattern_registry, mask_tile(ix_l,iy_l,:),&
               &nchan, entry_id)
               if (entry_id.lt.1) cycle ! defensive: should never happen

               ! Shared across all candidates: this pattern's own
               ! occurrence timeline doesn't depend on cache size, only
               ! on how many times it's been encountered so far.
               if (sim_next_occ_ptr(entry_id).le.&
               &pattern_registry%entries(entry_id)%n_occurrences) then
                  sim_next_occ_ptr(entry_id) = sim_next_occ_ptr(entry_id) + 1_8
               endif
               if (sim_next_occ_ptr(entry_id).gt.&
               &pattern_registry%entries(entry_id)%n_occurrences) then
                  new_next_occ = huge(new_next_occ)
               else
                  new_next_occ = pattern_registry%entries(entry_id)%&
                  &occurrences(sim_next_occ_ptr(entry_id))
               endif

               ! Reservoir-sample this pixel's own pattern for the
               ! per-block timing test below, regardless of hit/miss at
               ! the actual cache size -- standard reservoir sampling
               ! (Algorithm R): the k-th candidate seen is kept with
               ! probability n_sample_per_block/k, replacing a
               ! uniformly random existing slot if so.
               block_sample_seen = block_sample_seen + 1
               if (block_sample_filled.lt.n_sample_per_block) then
                  block_sample_filled = block_sample_filled + 1
                  block_sample_ids(block_sample_filled) = entry_id
               else
                  call random_number(rnd)
                  rnd_slot = 1 + int(rnd*real(block_sample_seen, sp))
                  if (rnd_slot.le.n_sample_per_block) then
                     block_sample_ids(rnd_slot) = entry_id
                  endif
               endif

               do kc = 1, n_candidates
                  if (sim_pattern_slot_c(entry_id,kc).gt.0) then
                     ! HIT for this candidate
                     cand_hits(kc) = cand_hits(kc) + 1_8
                     if (kc.eq.1) block_hits(n_blocks) = block_hits(n_blocks) + 1
                     cycle
                  endif

                  ! MISS for this candidate
                  cand_block_miss_nvalid_sum(kc) = cand_block_miss_nvalid_sum(kc)&
                  &+ int(pattern_nvalid(entry_id), 8)
                  if (kc.eq.1) block_miss_nvalid_sum(n_blocks) =&
                  &block_miss_nvalid_sum(n_blocks) + int(pattern_nvalid(entry_id), 8)

                  if (sim_n_resident_c(kc).lt.candidate_cache_max(kc)) then
                     sim_n_resident_c(kc) = sim_n_resident_c(kc) + 1
                     target_slot = sim_n_resident_c(kc)
                     sim_slot_pattern_c(target_slot,kc) = entry_id
                     sim_pattern_slot_c(entry_id,kc) = target_slot
                     cand_admitted(kc) = cand_admitted(kc) + 1_8
                     if (kc.eq.1) block_admitted(n_blocks) = block_admitted(n_blocks) + 1
                  else
                     do k2 = 1, sim_n_resident_c(kc)
                        if (sim_next_occ_ptr(sim_slot_pattern_c(k2,kc)).gt.&
                        &pattern_registry%entries(sim_slot_pattern_c(k2,kc))%&
                        &n_occurrences) then
                           victim_next_occs_sim(k2) = huge(new_next_occ)
                        else
                           victim_next_occs_sim(k2) = pattern_registry%&
                           &entries(sim_slot_pattern_c(k2,kc))%occurrences(&
                           &sim_next_occ_ptr(sim_slot_pattern_c(k2,kc)))
                        endif
                     enddo
                     target_slot = linear_scan_extreme(&
                     &victim_next_occs_sim(1:sim_n_resident_c(kc)),&
                     &sim_n_resident_c(kc), .true.)
                     if (new_next_occ.lt.victim_next_occs_sim(target_slot)) then
                        sim_pattern_slot_c(sim_slot_pattern_c(target_slot,kc),kc) = 0
                        sim_slot_pattern_c(target_slot,kc) = entry_id
                        sim_pattern_slot_c(entry_id,kc) = target_slot
                        cand_admitted(kc) = cand_admitted(kc) + 1_8
                        cand_evicted(kc) = cand_evicted(kc) + 1_8
                        if (kc.eq.1) then
                           block_admitted(n_blocks) = block_admitted(n_blocks) + 1
                           block_evicted(n_blocks) = block_evicted(n_blocks) + 1
                        endif
                     else
                        cand_declined(kc) = cand_declined(kc) + 1_8
                        if (kc.eq.1) then
                           block_declined(n_blocks) = block_declined(n_blocks) + 1
                           ! T15: true singleton (never needed again --
                           ! no cache size could ever help) vs recurs
                           ! later but declined anyway (a real,
                           ! cache-size-recoverable tradeoff).
                           if (new_next_occ.eq.huge(new_next_occ)) then
                              block_declined_singleton(n_blocks) =&
                              &block_declined_singleton(n_blocks) + 1
                              cand_declined_singleton = cand_declined_singleton + 1_8
                           else
                              block_declined_recurs(n_blocks) =&
                              &block_declined_recurs(n_blocks) + 1
                              cand_declined_recurs = cand_declined_recurs + 1_8
                           endif
                        endif
                     endif
                  endif
               enddo
            enddo
         enddo

         ! Per-block self-timing: time build_rmsf_offset_table on this
         ! block's own reservoir-sampled patterns -- 0 if this block
         ! had no misses at all at the actual cache size (all hits;
         ! e.g. a block whose own pixels are all patterns already
         ! resident from earlier blocks). Reused for EVERY candidate's
         ! own predicted time below (build cost is a property of the
         ! DATA's own channel coverage, not of cache policy).
         if (block_sample_filled.lt.1) then
            block_ms_per_channel(n_blocks) = 0.0_dp
         else
            call system_clock(count_rate=crate)
            sample_total_ms = 0.0_dp
            sample_total_channels = 0_8
            do k = 1, block_sample_filled
               nvalid_s = 0
               do k2 = 1, nchan
                  if (pattern_registry%entries(block_sample_ids(k))%&
                  &pattern(k2).ne.0_1) then
                     nvalid_s = nvalid_s + 1
                     valid_idx_s(nvalid_s) = k2
                  endif
               enddo
               if (nvalid_s.lt.1) cycle
               l_sq_valid_s(1:nvalid_s) = l_sq(valid_idx_s(1:nvalid_s))

               call system_clock(c0)
               call build_rmsf_offset_table(l_sq_valid_s(1:nvalid_s), nvalid_s,&
               &lsq_ref_compute, rm_samp(nrm)-rm_samp(1), real(cdelt3_amp, sp),&
               &table_oversample, table_s)
               call system_clock(c1)

               sample_total_ms = sample_total_ms +&
               &real(c1-c0, dp)/real(crate, dp)*1000.0_dp
               sample_total_channels = sample_total_channels + int(nvalid_s, 8)
               if (.not.have_bytes_per_table) then
                  bytes_per_table = 2.0_dp*real(table_s%n_fine, dp)*8.0_dp
                  have_bytes_per_table = .true.
               endif
            enddo
            if (sample_total_channels.lt.1_8) then
               block_ms_per_channel(n_blocks) = 0.0_dp
            else
               block_ms_per_channel(n_blocks) =&
               &sample_total_ms/real(sample_total_channels, dp)
            endif
         endif

         do kc = 1, n_candidates
            cand_predicted_total_min(kc) = cand_predicted_total_min(kc) +&
            &real(cand_block_miss_nvalid_sum(kc), dp)*block_ms_per_channel(n_blocks)&
            &/1000.0_dp/60.0_dp
         enddo
      end do
      deallocate(pattern_nvalid, sim_next_occ_ptr, sim_pattern_slot_c,&
      &sim_slot_pattern_c, victim_next_occs_sim, valid_idx_s, l_sq_valid_s)
      call omp_set_num_threads(saved_max_threads)

      write(*,'(A)') 'Pass 0 build-time advisory (table-GENERATION'//&
      &' only, not CLEAN -- see docs/dev/RMCLEAN_INTEGRATION_PLAN.md'//&
      &' T15). Hit/miss counts below are from a full Belady SIMULATION'//&
      &' (real admission-control eviction logic replayed against a'//&
      &' throwaway simulated cache, per-block build-cost sampled from'//&
      &' that block''s own real patterns), not the infinite-cache-'//&
      &'assumption first-occurrence count an earlier version of this'//&
      &' advisory used -- keep this log for retrospective analysis.'
      tot_hits = 0_8; tot_admitted = 0_8; tot_evicted = 0_8; tot_declined = 0_8
      tot_declined_singleton = 0_8; tot_declined_recurs = 0_8
      predicted_total_min = 0.0_dp

      csv_path = trim(outfile)//'.advisory.csv'
      open(newunit=csv_unit, file=trim(csv_path), status='replace',&
      &action='write', iostat=ios_csv)
      if (ios_csv.ne.0) then
         write(*,'(A)') 'WARNING: could not open '//trim(csv_path)//&
         &' for the T15 advisory CSV -- continuing without it.'
      else
         write(csv_unit,'(A)') 'block,hits,admitted,evicted_to_admit,'//&
         &'declined_oneoff,declined_singleton,declined_recurs,'//&
         &'ms_per_channel,predicted_build_min'
      endif

      do b = 1, n_blocks
         tot_hits = tot_hits + int(block_hits(b), 8)
         tot_admitted = tot_admitted + int(block_admitted(b), 8)
         tot_evicted = tot_evicted + int(block_evicted(b), 8)
         tot_declined = tot_declined + int(block_declined(b), 8)
         tot_declined_singleton = tot_declined_singleton + int(block_declined_singleton(b), 8)
         tot_declined_recurs = tot_declined_recurs + int(block_declined_recurs(b), 8)
         predicted_block_min = real(block_miss_nvalid_sum(b), dp)*&
         &block_ms_per_channel(b)/1000.0_dp/60.0_dp
         predicted_total_min = predicted_total_min + predicted_block_min
         write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,F0.4,A,F0.2,A)') '  block ',&
         &b, ': hits=', block_hits(b), ' admitted=', block_admitted(b),&
         &' (evicted-to-admit=', block_evicted(b), ') declined-oneoff=',&
         &block_declined(b), ' (true-singleton=', block_declined_singleton(b),&
         &' recurs-but-declined=', block_declined_recurs(b),&
         &') measured=', block_ms_per_channel(b),&
         &' ms/channel -> predicted build time ~', predicted_block_min,&
         &' min.'
         if (ios_csv.eq.0) then
            write(csv_unit,'(I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,F0.4,A,F0.4)')&
            &b, ',', block_hits(b), ',', block_admitted(b), ',',&
            &block_evicted(b), ',', block_declined(b), ',',&
            &block_declined_singleton(b), ',', block_declined_recurs(b),&
            &',', block_ms_per_channel(b), ',', predicted_block_min
         endif
      enddo
      if (ios_csv.eq.0) then
         close(csv_unit)
         write(*,'(A)') 'Wrote '//trim(csv_path)//&
         &' (per-block predicted table-generation time -- see'//&
         &' scripts/plot_rmclean_advisory.py to plot it).'
      endif
      write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') '  TOTAL (actual'//&
      &' mask_pattern_cache_max=', mask_pattern_cache_max, '): hits=', tot_hits,&
      &' admitted=', tot_admitted, ' (evicted-to-admit=', tot_evicted,&
      &') declined-oneoff=', tot_declined, ' (true-singleton=',&
      &tot_declined_singleton
      write(*,'(A,I0,A)') '    recurs-but-declined=', tot_declined_recurs,&
      &' -- these are the ones a BIGGER cache could actually recover,'//&
      &' not the true singletons, which no cache size can help.'
      write(*,'(A,F0.2,A)') '  predicted TOTAL build time (actual cache'//&
      &' size): ~', predicted_total_min, ' min (table generation only'//&
      &' -- actual wall time also includes CLEAN itself, which this'//&
      &' advisory deliberately does not estimate).'

      write(*,'(A)') '  Cache-size sweep (same simulation, larger'//&
      &' candidate cache sizes, reusing each block''s own already-'//&
      &' measured ms/channel rate -- NOT a case for a second, separately'//&
      &' -managed cache: Belady is already optimal for whatever total'//&
      &' capacity it is given, so two independent pools of the same'//&
      &' combined size can never beat one unified cache that size; this'//&
      &' sweep is for sizing the ONE existing cache, mask_pattern_'//&
      &'cache_max):'
      do kc = 1, n_candidates
         write(*,'(A,I0,A,I0,A,I0,A,F0.2,A)') '    cache_max=',&
         &candidate_cache_max(kc), ': admitted=', cand_admitted(kc),&
         &' declined-oneoff=', cand_declined(kc),&
         &' -> predicted build time ~', cand_predicted_total_min(kc),&
         &' min'
         if (have_bytes_per_table) then
            write(*,'(A,F0.2,A)') '      (approx RMSF-table memory at'//&
            &' this size: ~', real(candidate_cache_max(kc),dp)*&
            &bytes_per_table/1024.0_dp/1024.0_dp, ' MiB)'
         endif
      enddo
   end subroutine report_build_time_advisory

   subroutine evict_cache_slot(victim_idx)
      !! T14 (docs/dev/RMCLEAN_INTEGRATION_PLAN.md): marks victim_idx's
      !! OWN bucket entry as a TOMBSTONE (cache_buckets value -1) --
      !! standard open-addressing-with-deletion. Relocates that bucket
      !! slot by re-running victim_idx's own original insertion probe
      !! (starting from its own stored hash) and stopping at the first
      !! slot whose value literally equals victim_idx -- correct and
      !! sufficient because that exact (slot -> victim_idx) association
      !! was established once, at victim_idx's own original insertion,
      !! and open-addressing linear probing is deterministic given a
      !! fixed hash and an unchanged occupancy history; later tombstones
      !! elsewhere don't affect this, since lookups already treat
      !! tombstones as "keep probing," not "stop." Does NOT touch
      !! cache_entries(victim_idx) itself -- the caller overwrites its
      !! hash/pattern/table/hit_count immediately after this returns.
      integer, intent(in) :: victim_idx
      integer(kind=8) :: h, probe
      integer :: tries, slot
      h = cache_entries(victim_idx)%hash
      probe = modulo(h, int(n_cache_buckets, 8))
      do tries = 1, n_cache_buckets
         slot = int(probe)
         if (cache_buckets(slot).eq.victim_idx) then
            cache_buckets(slot) = -1
            return
         endif
         if (cache_buckets(slot).eq.0) exit ! shouldn't happen if invariants hold; defensive only
         probe = modulo(probe+1, int(n_cache_buckets, 8))
      enddo
   end subroutine evict_cache_slot

   subroutine update_mask_pattern_cache_for_tile(tx_in, ty_in)
      !! T10b: serial pre-scan of ONE tile's own mask_tile(1:tx_in,
      !! 1:ty_in,:) subrange, called once per tile, right after that
      !! tile's own read_mask_tile call and BEFORE that tile's own
      !! parallel per-pixel CLEAN loop starts -- same safety argument as
      !! the old one-upfront-global-pass design (cache_entries/
      !! cache_buckets/n_cache_entries are only ever written here,
      !! serially, never during the parallel loop), just repeated once
      !! per tile instead of once for the whole image. A pattern
      !! discovered while scanning an EARLIER tile remains cached and
      !! reusable by any LATER tile (cache_entries/cache_buckets persist
      !! across calls, initialised once by init_mask_pattern_cache, not
      !! reset here).
      !! T14: once the cache is full, a genuinely new pattern either
      !! evicts an existing entry (cache_eviction_policy='hitcount':
      !! least-hit-so-far; 'belady': offline-optimal admission-
      !! controlled eviction, increment 8 -- evicts whichever cached
      !! pattern's own next occurrence is farthest away, but ONLY if the
      !! new pattern's own next occurrence is strictly sooner; otherwise
      !! the cache is left untouched and this pixel's table is built as
      !! a one-off) or, for 'belady' when admission is declined, is
      !! simply not cached. Either way, target_slot is where the new
      !! pattern's own table gets built -- a genuinely NEW cache row
      !! (n_cache_entries incremented) or an EVICTED, now-reused one;
      !! the actual insert logic below is identical either way.
      integer, intent(in) :: tx_in, ty_in
      integer :: ix_l, iy_l, entry_idx, target_slot
      integer(kind=8) :: h, probe
      integer :: tries, slot
      integer :: nvalid_l
      integer, allocatable :: valid_idx_l(:)
      real(sp), allocatable :: l_sq_valid_l(:)
      logical :: do_insert
      ! Debug-only, per-table-build timing (real production cost, not
      ! T15's own artificial self-timing sample) -- log_message('debug',
      ! ...) is a no-op below the configured log_level (single integer
      ! compare before the critical section, see logging_mod.f90), so
      ! this costs nothing when log_level != debug.
      real(dp) :: t_table_build_start_l, t_table_build_ms_l
      character(len=160) :: table_build_msg_l
      ! T14 increment 7: this pixel's own pattern_registry entry id,
      ! looked up fresh only when this pixel's pattern is NOT already a
      ! cache hit (see cache_entries(entry_idx)%registry_id's own
      ! comment for why hits skip this).
      integer :: registry_id_l
      ! T14 increment 8: the Belady eviction decision's own working
      ! state -- new_next_occ is the NEW pattern's own next occurrence
      ! (after the current pixel is consumed); victim_next_occs(k) is
      ! cached slot k's own next occurrence, queried fresh every time a
      ! decision is needed (see this subroutine's own Belady branch for
      ! why a linear scan, not a maintained priority queue).
      integer(kind=8) :: new_next_occ
      integer(kind=8), allocatable :: victim_next_occs(:)
      integer :: k2

      allocate(valid_idx_l(nchan), l_sq_valid_l(nchan))
      allocate(victim_next_occs(max(1, mask_pattern_cache_max)))

      do iy_l = 1, ty_in
         do ix_l = 1, tx_in
            if (all(mask_tile(ix_l,iy_l,:).eq.0_1)) cycle
            ! T17: below the configured coverage threshold -- never
            ! cached, never counted as a distinct pattern; this pixel
            ! will be skipped entirely in clean_one_pixel too (which
            ! applies the exact same threshold independently), so
            ! building/caching a table for it here would be pure waste.
            if (min_valid_chan_frac.gt.0.0_sp .and.&
            &real(count(mask_tile(ix_l,iy_l,:).ne.0_1), sp).lt.&
            &min_valid_chan_frac*real(nchan, sp)) cycle

            h = fnv1a_hash(mask_tile(ix_l,iy_l,:), nchan)
            probe = modulo(h, int(n_cache_buckets, 8))
            entry_idx = -1
            do tries = 1, n_cache_buckets
               slot = int(probe)
               if (cache_buckets(slot).eq.0) exit
               if (cache_buckets(slot).ne.-1) then
                  if (cache_entries(cache_buckets(slot))%hash.eq.h) then
                     if (cached_pattern_matches(cache_buckets(slot), mask_tile(ix_l,iy_l,:))) then
                        entry_idx = cache_buckets(slot)
                        exit
                     endif
                  endif
               endif
               probe = modulo(probe+1, int(n_cache_buckets, 8))
            enddo
            if (entry_idx.ne.-1) then
               cache_entries(entry_idx)%hit_count = cache_entries(entry_idx)%hit_count + 1_8
               ! T14 increment 7: this pattern's own registry timeline
               ! advances regardless of caching -- a cache HIT still
               ! means this pixel's occurrence of the pattern has now
               ! been consumed. registry_id was set once, at this
               ! entry's own original insertion; no fresh lookup needed
               ! here.
               if (cache_eviction_policy.eq.'belady') then
                  call registry_advance(pattern_registry, cache_entries(entry_idx)%registry_id)
               endif
               cycle ! already cached
            endif

            ! T14 increment 7/8: this pattern's own registry consult --
            ! lookup + advance -- happens exactly ONCE per pixel miss,
            ! here, regardless of which branch below ends up handling
            ! it (genuinely-new-slot, evicted-and-reused, or declined
            ! one-off): the registry timeline is a fact about the
            ! PATTERN, not about how this particular miss gets
            ! resolved. new_next_occ (the pattern's own NEXT future
            ! need, after this occurrence) is only meaningful for the
            ! 'belady' eviction decision below, but computing it is
            ! cheap and harmless when unused (e.g. cache not yet full).
            if (cache_eviction_policy.eq.'belady') then
               call registry_lookup(pattern_registry, mask_tile(ix_l,iy_l,:), nchan, registry_id_l)
               call registry_advance(pattern_registry, registry_id_l)
               new_next_occ = registry_next_occurrence(pattern_registry, registry_id_l)
            endif

            do_insert = .true.
            if (n_cache_entries.ge.mask_pattern_cache_max) then
               if (cache_eviction_policy.eq.'hitcount') then
                  target_slot = linear_scan_extreme(&
                  &cache_entries(1:n_cache_entries)%hit_count, n_cache_entries, .false.)
                  call evict_cache_slot(target_slot)
                  n_distinct_patterns_total = n_distinct_patterns_total + 1
               else
                  ! T14 increment 8: 'belady', cache full -- the actual
                  ! offline-optimal eviction decision, with admission
                  ! control. Whichever CURRENTLY CACHED pattern has the
                  ! farthest-away next occurrence is Belady's own
                  ! textbook eviction candidate -- registry_next_
                  ! occurrence's huge() sentinel for "never needed
                  ! again" always wins find_max=.true., exactly
                  ! matching Belady's own "never again always evicts
                  ! first" rule.
                  do k2 = 1, n_cache_entries
                     victim_next_occs(k2) =&
                     &registry_next_occurrence(pattern_registry, cache_entries(k2)%registry_id)
                  enddo
                  target_slot = linear_scan_extreme(victim_next_occs(1:n_cache_entries),&
                  &n_cache_entries, .true.)

                  ! Admission control (the refinement over textbook
                  ! Belady): unlike a CPU cache, a miss here does NOT
                  ! force residency -- this pixel's table gets built
                  ! either way. Only evict-and-admit if the NEW
                  ! pattern's own next need is STRICTLY sooner than the
                  ! farthest victim's; otherwise leave the cache
                  ! untouched and build this pixel's table as a one-off
                  ! throwaway (the pre-T14 overflow behaviour).
                  if (new_next_occ.lt.victim_next_occs(target_slot)) then
                     call evict_cache_slot(target_slot)
                     n_distinct_patterns_total = n_distinct_patterns_total + 1
                  else
                     n_overflow_pixels_total = n_overflow_pixels_total + 1
                     do_insert = .false.
                  endif
               endif
            else
               n_cache_entries = n_cache_entries + 1
               target_slot = n_cache_entries
               n_distinct_patterns_total = n_distinct_patterns_total + 1
            endif
            if (.not. do_insert) cycle

            cache_entries(target_slot)%hash = h
            if (cache_eviction_policy.eq.'belady') then
               ! T14 increment 9: no local pattern-bytes copy under
               ! belady -- this slot's own pattern lives in the
               ! registry instead (Pass 0 already put it there),
               ! reachable via registry_id; cached_pattern_matches
               ! reads from there. See table_cache_entry_t's own
               ! pattern(:) comment for the full story.
               cache_entries(target_slot)%registry_id = registry_id_l
            else
               if (allocated(cache_entries(target_slot)%pattern)) deallocate(cache_entries(target_slot)%pattern)
               allocate(cache_entries(target_slot)%pattern(nchan))
               cache_entries(target_slot)%pattern = mask_tile(ix_l,iy_l,:)
            endif
            cache_entries(target_slot)%hit_count = 1_8

            ! Bucket slot for target_slot: reuse either a genuinely
            ! empty slot OR a tombstone left by an earlier eviction --
            ! either is valid, per evict_cache_slot's own comment.
            probe = modulo(h, int(n_cache_buckets, 8))
            do tries = 1, n_cache_buckets
               slot = int(probe)
               if (cache_buckets(slot).eq.0 .or. cache_buckets(slot).eq.-1) exit
               probe = modulo(probe+1, int(n_cache_buckets, 8))
            enddo
            cache_buckets(slot) = target_slot

            nvalid_l = 0
            do k = 1, nchan
               if (mask_tile(ix_l,iy_l,k).ne.0) then
                  nvalid_l = nvalid_l + 1
                  valid_idx_l(nvalid_l) = k
               endif
            enddo
            l_sq_valid_l(1:nvalid_l) = l_sq(valid_idx_l(1:nvalid_l))
            t_table_build_start_l = omp_get_wtime()
            call build_rmsf_offset_table(l_sq_valid_l(1:nvalid_l), nvalid_l,&
            &lsq_ref_compute, rm_samp(nrm)-rm_samp(1), real(cdelt3_amp, sp),&
            &table_oversample, cache_entries(target_slot)%table)
            t_table_build_ms_l = (omp_get_wtime()-t_table_build_start_l)*1000.0_dp
            write(table_build_msg_l,'(A,I0,A,I0,A,F0.3)')&
            &'table_build stage=cache_populate block=', tile_seq,&
            &' nvalid=', nvalid_l, ' dur_ms=', t_table_build_ms_l
            call log_message('debug','table_build',trim(table_build_msg_l))
         enddo
      enddo

      deallocate(valid_idx_l, l_sq_valid_l, victim_next_occs)
   end subroutine update_mask_pattern_cache_for_tile

   subroutine clean_one_pixel(ix_l, iy_l, ix0, iy0)
      !! One pixel's full CLEAN+restore, called from the main program's
      !! own per-tile `!$omp parallel do` over the tile's own local
      !! (ix_l,iy_l) in [1,tx]x[1,ty]. (ix0,iy0) is the tile's own origin
      !! (ix_tile_beg,iy_tile_beg) -- T10b: mask_tile is now tiled just
      !! like re_tile/im_tile, addressed at the SAME local (ix_l,iy_l),
      !! not the global (ix_g,iy_g) the old whole-cube-resident mask_cube
      !! needed; ix_g/iy_g are still computed and used, but now only for
      !! trace_ix/trace_iy matching (is_traced_p below), not for mask
      !! addressing. Every array declared here is a genuine LOCAL
      !! (automatic per-call) variable -- Fortran gives each concurrent
      !! call its own independent storage for these, with no `save` and
      !! no module-level state written here, so this subroutine is
      !! thread-safe purely by construction, no explicit `private()`
      !! bookkeeping needed for them. Everything this subroutine reads via
      !! host association (re_tile, mask_tile, rm_samp, l_sq, cache_
      !! entries/cache_buckets, restore_plan_fwd/bwd, ...) is READ-ONLY
      !! here; every array it WRITES (the 6 tile-output arrays) is written
      !! only at (ix_l,iy_l,:), a disjoint location per call -- no race
      !! either way.
      integer, intent(in) :: ix_l, iy_l, ix0, iy0
      integer :: ix_g, iy_g
      integer :: nvalid_p, n_iter_used_p, k_p, entry_idx_p
      character(len=16) :: stop_reason_p
      logical :: is_traced_p
      real(sp) :: trace_peak_p(niter), trace_rms_p(niter), trace_flux_p(niter)
      integer :: trace_iter_p
      logical :: trace_log_this_iter_p
      character(len=256) :: trace_msg_p
      real(sp) :: pixel_cleaned_flux_p, pixel_peak_resid_p
      integer, allocatable :: valid_idx_p(:)
      real(sp), allocatable :: l_sq_valid_p(:)
      real(sp) :: dirty_re_p(nrm), dirty_im_p(nrm)
      real(sp) :: comp_re_p(nrm), comp_im_p(nrm)
      real(sp) :: resid_re_p(nrm), resid_im_p(nrm)
      real(sp) :: comp_rm_refined_p(nrm)
      real(sp) :: restored_re_p(nrm), restored_im_p(nrm)
      type(rmsf_table_t) :: throwaway_table
      logical :: used_throwaway
      ! Debug-only per-table-build timing for the one-off/throwaway
      ! path -- see update_mask_pattern_cache_for_tile's own comment on
      ! why this is free when log_level != debug. This call site runs
      ! NESTED inside the outer per-pixel !$omp parallel do, so
      ! build_rmsf_offset_table's own internal parallel do is a no-op
      ! here (nesting disabled) -- this measures that genuinely
      ! single-threaded-in-practice cost, distinct from cache_populate's
      ! own (not nested, really multi-threaded) cost above.
      real(dp) :: t_table_build_start_p, t_table_build_ms_p
      character(len=160) :: table_build_msg_p

      ix_g = ix0 + ix_l - 1
      iy_g = iy0 + iy_l - 1

      nvalid_p = 0
      allocate(valid_idx_p(nchan), l_sq_valid_p(nchan))
      do k_p = 1, nchan
         if (mask_tile(ix_l,iy_l,k_p).ne.0) then
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
         clean_re_tile(ix_l,iy_l,:) = re_tile(ix_l,iy_l,:)
         clean_im_tile(ix_l,iy_l,:) = im_tile(ix_l,iy_l,:)
         resid_re_tile(ix_l,iy_l,:) = re_tile(ix_l,iy_l,:)
         resid_im_tile(ix_l,iy_l,:) = im_tile(ix_l,iy_l,:)
         restored_re_tile(ix_l,iy_l,:) = re_tile(ix_l,iy_l,:)
         restored_im_tile(ix_l,iy_l,:) = im_tile(ix_l,iy_l,:)
         deallocate(valid_idx_p, l_sq_valid_p)
         return
      endif

      if (min_valid_chan_frac.gt.0.0_sp .and.&
      &real(nvalid_p, sp).lt.min_valid_chan_frac*real(nchan, sp)) then
         ! T17 (docs/dev/RMCLEAN_INTEGRATION_PLAN.md): below the
         ! configured coverage threshold -- UNLIKE the nvalid_p<1 case
         ! above, the dirty spectrum here is real, non-NaN data (this
         ! pixel does have SOME valid channels); passing it through
         ! would silently present un-CLEANed dirty data as if it were a
         ! real CLEAN/RESID/RESTORED result. Write explicit NaN
         ! instead, same runtime-0.0/0.0 IEEE NaN convention
         ! rm_synthesis_mod.f90 itself already uses (zero_val_p, a
         ! genuine runtime variable -- a literal 0.0_sp/0.0_sp would be
         ! rejected at compile time).
         block
            real(sp) :: zero_val_p, nan_val_p
            zero_val_p = 0.0_sp
            nan_val_p = zero_val_p/zero_val_p
            clean_re_tile(ix_l,iy_l,:) = nan_val_p
            clean_im_tile(ix_l,iy_l,:) = nan_val_p
            resid_re_tile(ix_l,iy_l,:) = nan_val_p
            resid_im_tile(ix_l,iy_l,:) = nan_val_p
            restored_re_tile(ix_l,iy_l,:) = nan_val_p
            restored_im_tile(ix_l,iy_l,:) = nan_val_p
         end block
         deallocate(valid_idx_p, l_sq_valid_p)
         return
      endif

      ! l_sq_valid_p is filled unconditionally now (not just for the
      ! throwaway-table branch): clean_complex's own T3 refinement step
      ! (refine_peak_matched_filter) needs the valid-channel l_sq list
      ! regardless of which table (cached or throwaway) was used to
      ! build the SUBTRACTION beam.
      l_sq_valid_p(1:nvalid_p) = l_sq(valid_idx_p(1:nvalid_p))

      is_traced_p = (trace_ix.gt.0 .and. ix_g.eq.trace_ix .and. iy_g.eq.trace_iy)

      call cache_lookup_readonly(mask_tile(ix_l,iy_l,:), entry_idx_p)
      used_throwaway = (entry_idx_p.eq.0)
      if (used_throwaway) then
         t_table_build_start_p = omp_get_wtime()
         call build_rmsf_offset_table(l_sq_valid_p(1:nvalid_p), nvalid_p,&
         &lsq_ref_compute, rm_samp(nrm)-rm_samp(1), real(cdelt3_amp, sp),&
         &table_oversample, throwaway_table)
         t_table_build_ms_p = (omp_get_wtime()-t_table_build_start_p)*1000.0_dp
         write(table_build_msg_p,'(A,I0,A,I0,A,I0,A,F0.3)')&
         &'table_build stage=oneoff block=', tile_seq, ' tid=',&
         &omp_get_thread_num(), ' nvalid=', nvalid_p,&
         &' dur_ms=', t_table_build_ms_p
         call log_message('debug','table_build',trim(table_build_msg_p))
      endif

      dirty_re_p = re_tile(ix_l,iy_l,:)
      dirty_im_p = im_tile(ix_l,iy_l,:)

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
         if (is_traced_p) then
            call clean_complex(l_sq_valid_p(1:nvalid_p), nvalid_p,&
            &lsq_ref_compute, rm_samp, nrm, dirty_re_p, dirty_im_p,&
            &throwaway_table, niter, gain,&
            &have_abs_flux_floor, abs_flux_floor,&
            &have_auto_nsigma, auto_nsigma_mult,&
            &comp_re_p, comp_im_p, resid_re_p, resid_im_p, n_iter_used_p,&
            &stop_reason_p, comp_rm_refined_p, nsigma_refine=refine_nsigma,&
            &trace_peak_val=trace_peak_p, trace_rms_val=trace_rms_p,&
            &trace_flux_val=trace_flux_p)
         else
            call clean_complex(l_sq_valid_p(1:nvalid_p), nvalid_p,&
            &lsq_ref_compute, rm_samp, nrm, dirty_re_p, dirty_im_p,&
            &throwaway_table, niter, gain,&
            &have_abs_flux_floor, abs_flux_floor,&
            &have_auto_nsigma, auto_nsigma_mult,&
            &comp_re_p, comp_im_p, resid_re_p, resid_im_p, n_iter_used_p,&
            &stop_reason_p, comp_rm_refined_p, nsigma_refine=refine_nsigma)
         endif
         call destroy_rmsf_offset_table(throwaway_table)
      else
         if (is_traced_p) then
            call clean_complex(l_sq_valid_p(1:nvalid_p), nvalid_p,&
            &lsq_ref_compute, rm_samp, nrm, dirty_re_p, dirty_im_p,&
            &cache_entries(entry_idx_p)%table, niter, gain,&
            &have_abs_flux_floor, abs_flux_floor,&
            &have_auto_nsigma, auto_nsigma_mult,&
            &comp_re_p, comp_im_p, resid_re_p, resid_im_p, n_iter_used_p,&
            &stop_reason_p, comp_rm_refined_p,&
            &nsigma_refine=refine_nsigma, trace_peak_val=trace_peak_p,&
            &trace_rms_val=trace_rms_p, trace_flux_val=trace_flux_p)
         else
            call clean_complex(l_sq_valid_p(1:nvalid_p), nvalid_p,&
            &lsq_ref_compute, rm_samp, nrm, dirty_re_p, dirty_im_p,&
            &cache_entries(entry_idx_p)%table, niter, gain,&
            &have_abs_flux_floor, abs_flux_floor,&
            &have_auto_nsigma, auto_nsigma_mult,&
            &comp_re_p, comp_im_p, resid_re_p, resid_im_p, n_iter_used_p,&
            &stop_reason_p, comp_rm_refined_p, nsigma_refine=refine_nsigma)
         endif
      endif

      ! Per-pixel derived diagnostics: computed straight from clean_
      ! complex's own OUTPUT arrays (no new outputs needed from that
      ! module for these two) -- total flux CLEAN actually extracted
      ! into components, and the peak left behind in the residual once
      ! it stopped.
      pixel_cleaned_flux_p = sum(sqrt(comp_re_p**2 + comp_im_p**2))
      pixel_peak_resid_p = maxval(sqrt(resid_re_p**2 + resid_im_p**2))

      if (is_traced_p) then
         do trace_iter_p = 1, n_iter_used_p
            ! Always log iter=1 (pre-CLEAN state, see clean_complex's own
            ! trace_flux_val comment) and the final iteration; otherwise
            ! throttle to every log_every-th iteration.
            trace_log_this_iter_p = (trace_iter_p.eq.1) .or.&
            &(trace_iter_p.eq.n_iter_used_p) .or.&
            &(mod(trace_iter_p, log_every).eq.0)
            if (.not.trace_log_this_iter_p) cycle
            write(trace_msg_p,'(A,I0,A,I0,A,I0,A,ES14.6,A,ES14.6,A,ES14.6)')&
            &'CLEAN trace ix=', ix_g, ' iy=', iy_g, ' iter=', trace_iter_p,&
            &' peak_val=', trace_peak_p(trace_iter_p), ' rms_val=',&
            &trace_rms_p(trace_iter_p), ' cleaned_flux_so_far=',&
            &trace_flux_p(trace_iter_p)
            call log_message('info','clean_trace',trim(trace_msg_p))
         enddo
         write(trace_msg_p,'(A,I0,A,I0,A,I0,A,I0,A,A)')&
         &'CLEAN trace ix=', ix_g, ' iy=', iy_g, ' STOPPED after ',&
         &n_iter_used_p, ' of ', niter, ' iterations, stop_reason=',&
         &trim(stop_reason_p)
         call log_message('info','clean_trace',trim(trace_msg_p))
         write(trace_msg_p,'(A,I0,A,I0,A,ES14.6,A,ES14.6)')&
         &'CLEAN trace ix=', ix_g, ' iy=', iy_g, ' total_cleaned_flux=',&
         &pixel_cleaned_flux_p, ' peak_residual=', pixel_peak_resid_p
         call log_message('info','clean_trace',trim(trace_msg_p))
      endif

      !$omp atomic
      n_iter_used_sum = n_iter_used_sum + int(n_iter_used_p, 8)
      !$omp atomic
      n_iter_used_min = min(n_iter_used_min, n_iter_used_p)
      !$omp atomic
      n_iter_used_max = max(n_iter_used_max, n_iter_used_p)
      !$omp atomic
      total_cleaned_flux_sum = total_cleaned_flux_sum + real(pixel_cleaned_flux_p, dp)
      !$omp atomic
      peak_residual_sum = peak_residual_sum + real(pixel_peak_resid_p, dp)
      !$omp atomic
      peak_residual_max = max(peak_residual_max, pixel_peak_resid_p)
      !$omp atomic
      peak_residual_min = min(peak_residual_min, pixel_peak_resid_p)
      select case (trim(stop_reason_p))
      case ('abs_flux')
         !$omp atomic
         n_stopped_abs_flux = n_stopped_abs_flux + 1_8
      case ('auto_nsigma')
         !$omp atomic
         n_stopped_auto_nsigma = n_stopped_auto_nsigma + 1_8
      case default
         !$omp atomic
         n_stopped_niter = n_stopped_niter + 1_8
      end select

      call restore_clean(rm_samp, nrm, comp_re_p, comp_im_p, resid_re_p,&
      &resid_im_p, fwhm_rm, restore_plan_fwd, restore_plan_bwd,&
      &restored_re_p, restored_im_p)

      if (lsq_ref_report.ne.lsq_ref_compute) then
         ! comp_re_p/comp_im_p: each bin's flux was accumulated at its own
         ! FLUX-WEIGHTED sub-pixel location (comp_rm_refined_p), not
         ! rm_samp(j) exactly (clean_complex's own comp_rm_refined
         ! comment) -- derotate_to_lsq_ref's rotation angle is location-
         ! dependent, so using rm_samp(j) here silently reintroduces the
         ! same ~dRM/2-scale chi0 error comp_rm_refined was added to
         ! eliminate (tests/test_rmclean_lsqref_flex.f90's own
         ! ipeak/rm_found pattern is the reference this mirrors). Confirmed
         ! this only matters when lsq_ref_compute.ne.lsq_ref_report (this
         ! branch's own guard) -- moot under the old lsq_ref_compute_mode=
         ! native default (ref_diff was usually 0), live now that mid is
         ! the default.
         call derotate_to_lsq_ref(comp_rm_refined_p, nrm, lsq_ref_compute,&
         &lsq_ref_report, comp_re_p, comp_im_p, comp_re_p, comp_im_p)
         ! resid_re_p/resid_im_p and restored_re_p/restored_im_p: both are
         ! genuine regular-grid functions (resid: compute_dirty_rmbeam
         ! evaluates the beam AT rm_samp(j) exactly for every j; restored:
         ! restore_clean's own FFT convolution produces a value that
         ! genuinely lives AT rm_samp(j) after smoothing) -- rm_samp is the
         ! correct location array for these two, unlike comp_re_p/comp_im_p
         ! above.
         call derotate_to_lsq_ref(rm_samp, nrm, lsq_ref_compute,&
         &lsq_ref_report, resid_re_p, resid_im_p, resid_re_p, resid_im_p)
         call derotate_to_lsq_ref(rm_samp, nrm, lsq_ref_compute,&
         &lsq_ref_report, restored_re_p, restored_im_p, restored_re_p,&
         &restored_im_p)
      endif

      clean_re_tile(ix_l,iy_l,:) = comp_re_p
      clean_im_tile(ix_l,iy_l,:) = comp_im_p
      resid_re_tile(ix_l,iy_l,:) = resid_re_p
      resid_im_tile(ix_l,iy_l,:) = resid_im_p
      restored_re_tile(ix_l,iy_l,:) = restored_re_p
      restored_im_tile(ix_l,iy_l,:) = restored_im_p

      deallocate(valid_idx_p, l_sq_valid_p)
      !$omp atomic
      n_pixels_done = n_pixels_done + 1
   end subroutine clean_one_pixel

end program rmclean_cubes
