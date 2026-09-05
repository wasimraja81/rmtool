# Changelog

All notable changes to this project are documented in this file.

## [Unreleased] - RM-CLEAN integration

RM-CLEAN made practically usable end-to-end: a new standalone tool,
`rmclean_cubes`, drives the existing `rmclean_mod` core (Högbom-style
complex CLEAN, Gaussian restore) against REAL dirty AMP/PHA cubes
`rm_synthesis` itself writes, rather than only the synthetic spectra its
own unit test programs built in memory. Full design rationale, decisions
recorded with the user, and ticket-by-ticket verification evidence lives
in `docs/dev/RMCLEAN_INTEGRATION_PLAN.md` (tickets T0-T10b done; T11 not
yet started) — this entry is a summary, not a replacement for that
record.

### Added — `rmclean_cubes` standalone tool (ticket T2)
- New standalone program `src/rmclean_cubes.f90` (own Makefile target
  `make rmclean_cubes`, `bin/rmclean_cubes`), following
  `reproject_cubes`/`convolve_cubes`'s own conventions (own key=value/
  `--config` parser, `--help`, own Makefile block, wired into
  `scripts/make_all.sh` and `docker/dockerfile`).
- **Gate 0**: since this tool cannot resample the RM axis (`CDELT3`/`nrm`
  are fixed by whatever `rm_synthesis` already wrote), it validates the
  existing grid against `get_drm`'s own sampling bound instead — refuses
  to proceed (clear error) rather than silently CLEANing on an
  undersampled grid.
- `lsq_ref_compute` (this tool's own RMSF-table/CLEAN reference) is a
  free, independent choice (`lsq_ref_compute_mode=native|zero|mid|
  centroid|min|max|fixed`), applied via an exact `derotate_to_lsq_ref`
  phase rotation of the already-sampled dirty spectrum — no accuracy
  cost, though also no compute-cost benefit once the RM grid already
  exists (that benefit only applies upstream, at `rm_synthesis`'s own
  new `lsq_ref_mode`, below). Gate 0 always validates against the
  cube's own actual reference (`lsq_ref_native`, read from a new
  `LSQREF` header keyword — see below), never against this choice.
- **Mask-pattern cache**: pixels sharing the same valid-channel mask
  pattern share one RMSF table, built once during a serial pre-scan and
  looked up via a hash-bucketed (open-addressing, collision-safe) table
  for O(1) amortized lookup; capped at `mask_pattern_cache_max` (default
  4096) distinct patterns, with a one-off throwaway table past the cap.
- **OpenMP parallelism**: one thread per pixel (embarrassingly parallel
  along the RM axis); verified bit-identical output across thread
  counts and across cache configurations (generous/disabled/overflowing).
- `restore_fwhm`, `threshold`/`niter`/`gain`/`oversample`/
  `table_oversample`, `lsq_ref_report_mode`/`_value` round out the config
  surface; `--help` documents all of it.

### Added — `rm_synthesis`'s own configurable phase reference (ticket T2b)
- New `lsq_ref_mode` (`zero`|`mid`|`centroid`|`min`|`max`|`fixed`, plus
  `lsq_ref_fixed_value`) config keys, mirroring `rmclean_mod`'s own
  `get_lsq_ref_compute` mode+value strategy (duplicated logic, not a new
  module dependency — this module is the older, heavily
  production-tested core). Default `zero` preserves every existing cfg
  file's behaviour exactly, unaffected unless set explicitly.
- The actual phase reference used is now recorded on both AMP/PHA cubes
  as a new `LSQREF` header keyword, so `rmclean_cubes` (or any other
  downstream reader) never has to assume a value; falls back to 0.0
  (this project's historical convention) with a printed warning for
  cubes written before this keyword existed.

### Added — model-based matched-filter peak refinement (ticket T3)
- `clean_complex`'s own sub-pixel peak refinement no longer fits a local
  parabola directly through the raw, stored `Re`/`Im` samples
  (`peak_interp_parabolic`'s old role) — that fit was only valid if the
  RM axis already satisfied `get_drm`'s demanding `max_offset`-based
  bound (real numbers, thesis P-band: ~22-44x oversampling relative to
  resolution `fwhm` at `lsq_ref=0`), because those samples carry a fast,
  `lsq_ref`-dependent carrier a low-order local fit can't resolve.
- New `rmsf_point_direct` (exact, O(nchan), single-offset RMSF
  evaluation) and `refine_peak_matched_filter` (a local matched-filter
  search against that exact model, using only the nearest stored coarse
  samples as anchors) replace it. The RMSF is analytically known from
  channel λ² alone — no raw per-channel data needed — so this removes
  the outer grid's dependence on `get_drm`'s bound entirely; only
  ordinary resolution-level sampling is needed now.
- Landed in 3 independently-tested phases (isolated validation, additive
  new subroutines, then the `clean_complex` swap), each gated on its own
  passing test, per the project's "validate all claims before
  production" policy. `clean_complex`'s own signature gained
  `l_sq`/`nchan`/`lsq_ref_compute` (needed by the new refinement, not by
  the old parabola) — updated at all 5 call sites.
- A real, measured performance cost was found and addressed during
  validation: the first working tuning (5 anchors, 200 samples/carrier-
  cycle) took `rmclean_cubes` from sub-second to ~1m40s on a 1024-pixel
  test cube. Retuned with evidence (3 vs. 5 anchors were bit-identical
  in every noiseless test; the real failure boundary was found ~7-8
  samples/cycle, with 50/cycle chosen as a checked >=6x margin): ~31.6s
  on the same cube (~3.2x). Flagged to the user as not fully closed —
  a coarse-to-fine two-stage search (ticket T3c, not yet started) should
  close the remaining gap for genuinely large cubes.
### Added — tiered peak refinement + resolution-based Gate 0 (tickets T3b/T3c)
- `refine_peak_matched_filter` is now TIERED (T3c, designed jointly with
  the user): a cheap fast path runs every CLEAN iteration — peak
  location from the log-magnitude parabola, complex amplitude solved
  closed-form against the ≥2 nearest stored anchors at that fixed
  location — and its leftover misfit doubles as a self-consistency
  diagnostic (deliberately NOT a single-point division, which always
  "fits" and can never signal a wrong location). Only when that misfit
  exceeds `nsigma × noise_rms` (noise: the same per-iteration
  `rms_about_mean` the stopping criterion already uses) does the full
  T3 local search run for that iteration. `nsigma` is user-facing
  (`nsigma_refine` optional argument on `clean_complex`,
  `refine_nsigma=` key on `rmclean_cubes`, default 3.0, confirmed
  against a deterministic noisy test scenario). Search constants
  re-anchored physically (statistic oscillates at 2× the carrier rate;
  `m_floor` cut 100→20, tied to `table_oversample`'s own sub-cell
  meaning). Measured: the 1024-pixel test cube went 31.6s → **0.46s**
  (~217x faster than T3's first working version), sources recovered
  unchanged.
- `rmclean_cubes`' Gate 0 recast (T3b) as a pure RM-RESOLUTION
  criterion: `|CDELT3| ≤ fwhm/min_samples_per_fwhm` (fwhm from
  `compute_rmsf_fwhm_multiband`, `lsq_ref`-independent; always the
  data's own fwhm, never the `restore_fwhm` override). New
  `min_samples_per_fwhm=` key (default 2, hard floor 1) replaces the
  retired carrier-based `oversample=` key. Verified both directions
  (real cube passes; CDELT3 forged coarse is refused at default,
  accepted at the floor).
- A real input-contract boundary found by the suite itself: at a
  band-centroid reference the carrier and resolution scales coincide,
  so `test_drm_floor.f90`'s `oversample=1` case is 0.5 samples/fwhm —
  below resolution Nyquist, outside the new contract — and the tiered
  default correctly cannot be relied on there (the sidelobe-dominated
  early-iteration rms makes escalation too lenient to catch it). The
  test now expects failure again with the NEW reason documented
  (sub-resolution input — exactly what Gate 0 refuses in production),
  and its `oversample=2` (1 sample/fwhm) case confirms the default-2
  gate carries genuine margin. `get_drm` itself is unchanged — still
  correct for its remaining jobs (T3c's internal search sizing,
  `rm_synthesis`-side grid planning).

### Fixed
- Two stale doc comments (`src/rmclean.f90`, `tests/thesis_scenario_
  rmclean.f90`) incorrectly claimed `rm_synthesis_mod.f90`'s own
  `extract_general_setup` references the band's own mean λ² — it is, and
  always was, unconditionally at `lambda_sq=0`.
- A real robustness bug in `rmclean_cubes` found via testing, not
  inspection: `FTOPEN`/`FTINIT` return status was not checked before
  subsequent CFITSIO calls in `read_cube`/`read_mask_cube`/
  `write_output_cube` — a failed open/create (e.g. re-running the tool
  without cleaning up a previous run's own output files) left later
  calls operating on a CFITSIO unit that was never actually set up,
  reproduced as a real SIGSEGV. Fixed by checking status immediately
  after every `FTOPEN`/`FTINIT` and bailing out cleanly.

### Verification
- Full regression 85/85 (`tests/run_tests.sh` sections 23-29 cover all
  of T1/T2/T2b/T3): section 28, `RM-CLEAN matched-filter peak
  refinement`, validates the new subroutines against 2 independent
  point-source scenarios at 4 sampling densities each; section 29,
  `rmclean_cubes` end-to-end against a real `lsq_ref_mode=mid` cube,
  covering the `LSQREF` header round-trip, Gate 0, both known injected
  point sources recovered via `check_rm_peak.py`, and mask-pattern cache
  correctness.
- Manually verified: known injected point sources recovered with exact
  RM, ~0.1-0.2% amplitude error, ~0.05° intrinsic-angle error, across
  both the native and an explicitly-derotated `lsq_ref_compute`; 1-thread
  vs. 4-thread and cache-big/zero/small runs all bit-identical.

### Added — memory-budgeted, threaded block I/O for `rmclean_cubes` (ticket T4)
- `rmclean_cubes` moves off "whole-cube-in-memory" onto the SAME
  scheme `rm_synthesis` already uses in production: spatial tiles
  (`tile_ra`/`tile_dec`, RA-strips-first auto policy budgeted by
  `mem_frac_ram`), parallel readonly chunked tile reads
  (`io_read_threads`), raw-stream-write tile output bypassing CFITSIO
  (`io_write_threads`), and pthread-based double-buffered write-behind
  (`io_overlap`) — same key names/semantics as `rm_synthesis`,
  generalized from its 2 named AMP/PHA outputs to `rmclean_cubes`'s 2
  inputs + 6 outputs (CLEAN/RESID/RESTORED x AMP/PHA). The mask cube
  and its mask-pattern -> RMSF-table cache deliberately stay
  whole-cube-resident (tiny relative to the float cubes; the cache
  needs one global pre-scan before any pixel is CLEANed). **Superseded
  by ticket T10b below**: the mask cube is now tiled and the cache
  builds incrementally instead.
- Found and root-caused, not assumed: a forced small-tile run differs
  numerically (not just byte-for-byte) from the default single-tile
  run — confirmed to be `gfortran -O3 -march=native`'s own
  floating-point reassociation of `clean_complex`'s tiered-refinement
  threshold comparison (T3c), sensitive to the runtime memory alignment
  of its own stack arguments, which genuinely differs between tile
  sizes; bit-identical at `-O0`. Not a tiling bug — the regression test
  (`tests/check_tile_consistency.py`) therefore compares with tolerance
  (AMP), not byte-identity, unlike the purely linear reproject_cubes/
  convolve_cubes.

### Fixed
- Concurrent `open(newunit=...)` from different OpenMP threads (each
  opening the same output file path for its own disjoint byte range)
  intermittently produced a corrupted output cube — the Fortran
  standard does not guarantee any I/O statement is safe to call
  concurrently without explicit synchronization. Fixed (in both
  `rmclean_cubes.f90` and `rm_synthesis_mod.f90`'s own
  `write_rm_chunk_raw`, which has the identical pattern) by wrapping
  just the `open(newunit=...)` call in a named `!$omp critical` —
  genuinely unique, no manual unit-range bookkeeping to maintain,
  only the brief allocation step serialized. 20 repeated
  `io_write_threads=4` runs against each tool, 0 mismatches.
- **This system's installed libcfitsio's `FTGHAD` writes only the lower
  32 bits of its 3 output arguments**, leaving an `integer(kind=8)`
  receiving variable's upper 32 bits at whatever was already there —
  confirmed with a minimal standalone reproducer outside this codebase.
  Caused an intermittent (stack-content-dependent), roughly-50%-of-runs
  wrong byte offset in `rmclean_cubes`'s own `io_write_threads>1` path.
  Fixed by zero-initializing the 3 receiving variables immediately
  before every `FTGHAD` call — in `rmclean_cubes.f90` AND in
  `rm_synthesis.f90`'s own pre-existing `io_write_threads` FTGHAD call,
  which has the identical latent exposure (never observed to misbehave
  there only because its receiving variables happen to land in
  zero-initialized static storage on this platform/compiler, not the
  stack — incidental, not guaranteed).

### Verification
- Full regression 93/93 (`tests/run_tests.sh` section 29): forced
  small-tile run (tolerance comparison + known-source RM recovery),
  `io_read_threads=1,2,4` (bit-identical), `io_write_threads=1,2,4`
  (bit-identical, 5 reps each — the FTGHAD bug was probabilistic, so a
  single pass would not have reliably caught it), `io_overlap=y` alone
  (bit-identical, 5 reps) and combined with forced small tiles +
  `io_read_threads=3` + `io_write_threads=3` together (bit-identical,
  5 reps) — the actual combined stress case, not each mechanism only
  in isolation.

### Added — CLEAN convergence/stop-reason logging + gain-tuning workflow (ticket T7)
- `clean_complex` now returns a `stop_reason` string and, for one
  traced pixel (`trace_ix=`/`trace_iy=`, throttled by `log_every=`,
  default 50), a full per-iteration trend (`peak_val`/`rms_val`/
  cumulative cleaned flux) — a debugging/tuning aid, not for
  whole-cube production runs.
- Per-thread block-progress logging (`m of n blocks processed`) and a
  run-end aggregate stop-reason summary (counts and percentages per
  criterion, plus `n_iter_used` mean/min/max), so a real production run
  reports how it actually converged, not just that it finished.
- A reusable subimage-based gain-tuning workflow (`rm_synthesis`'s own
  `subim=y` cutout, not a separate extraction tool) surfaced a real,
  actionable finding about the OLD stopping-criterion mechanism (see
  ticket T8 below) — superseded, not contradicted, by that redesign.

### Added — CLEAN stopping-criteria redesign (ticket T8)
- The old `threshold=`/`threshold_snr=` pair funneled into one
  overloaded, always-multiplicative `thresh` internal variable that
  silently broke `threshold=`'s absolute-flux mode (compared `flux^2`
  against `flux`) and made `threshold_snr=`'s own n-sigma mode a
  self-referential moving target (compared against the CURRENT,
  shrinking residual's own rms every iteration, not a fixed noise
  floor). Found via the T7 gain-sweep workflow, not by inspection.
- Replaced entirely (not aliased) by three independent, unambiguously-
  named, freely-combinable criteria: `niter` (hard cap, unchanged),
  `abs_flux_floor=` (a genuine fixed-flux comparison), `auto_nsigma=`
  (multiplier x a per-pixel noise sigma, estimated ONCE per pixel from
  its own dirty spectrum -- not a whole-cube pre-scan, not recomputed
  per iteration). `threshold=`/`threshold_snr=`/`noise_nlos=`/
  `noise_percentile=`/`noise_seed=` are all retired.
- Fixed T7's own over-cleaning finding at the root: re-running T7's
  subimage with `auto_nsigma=5.0` (the old `threshold_snr=5.0`
  production value) took tail RMS from 6.89 uJy (3.7x below the true
  ~25.5 uJy floor) to 25.24 uJy (matching it), converging in a mean of
  3.89 iterations instead of 233.7.

### Added — `auto_nsigma` correctness: full-spectrum IQR sigma (ticket T9)
- A second, independent bug found the same day: `auto_nsigma` compared
  `peak_val - avg_abs`, not `peak_val` directly -- `avg_abs` (the whole
  residual spectrum's own mean absolute value) is recomputed from the
  CURRENT, shrinking residual every iteration, the same class of
  self-referential bug T8 had just fixed for `thresh`. For a large-scale
  coherent structure (peak roughly equal to average almost everywhere)
  this collapsed to ~0 and stopped CLEAN early regardless of the true
  noise level. Fixed by comparing `peak_val` directly, no subtraction.
- The per-pixel noise-sigma estimator itself went through two more
  iterations, each measured against real data rather than assumed: a
  peak-relative exclusion window (T8's own original design) failed for
  multi-component/extended sightlines; a percentile-subset variant
  (lowest 20% of bins, truncated-Rayleigh-corrected) was tried, then
  measured directly against ten real pixels' own independently-known
  true noise (far-RM-tail ground truth) and found to underperform the
  simplest option -- the pixel's own FULL dirty amplitude spectrum's
  interquartile range, converted to sigma via the analytic untruncated
  Rayleigh-distribution relation (`estimate_iqr_sigma`, fixed factor
  `IQR/0.90656`). Simpler and more accurate; shipped.
- cfg defaults finalised, each justified per-file rather than copied:
  both criteria off unless explicitly opted into; `auto_nsigma`'s own
  default value is `1.0` (was `5.0` under the old design, no longer
  needed once the sigma estimate itself is trustworthy);
  `abs_flux_floor` values always carry an explicit unit suffix now.

### Fixed — 32-bit integer overflow in `read_mask_cube` (ticket T10a)
- A real, serious, silent data-corruption bug: the mask cube was read
  in one `FTGPVB` call whose `nelem` argument was computed as
  `nx*ny*nchan` in ordinary 32-bit integer arithmetic. Any mask cube
  exceeding 2^31 total elements overflows this silently -- a real ASKAP
  dataset's own 4501x4501x288 mask cube is 2.7x over that limit,
  silently truncating the valid-channel read at ~76 of 288 channels for
  EVERY pixel in the image.
- This single bug was the root cause of three things independently
  investigated as separate mysteries on a real full-cube run: a
  specific pixel appearing to never converge in isolation, the
  run's own 99.05% niter-cap-hit rate across all 20,259,001 pixels, and
  spurious CLEAN components with residual exceeding dirty amplitude at
  supposedly signal-free RM planes.
- Fixed by switching to `FTGPVBLL`, CFITSIO's own genuinely-64-bit
  entry point for the identical underlying C call (confirmed directly
  against this project's own bundled cfitsio-4.3.1 source, not
  assumed) -- correct at any dataset size, not just today's.
- `read_mask_cube` was extracted out of `rmclean_cubes`'s own program
  scope into a new module, `src/rmclean_io_mod.f90`, purely so it could
  be independently unit-tested (a fixture exceeding 2^31 elements is
  unavoidably ~2.4GB on disk for a byte-typed cube; kept fast via an
  OS-level sparse file -- only 3 bytes actually written).
- A full real-data confirmation run (replacing the flawed pre-fix
  output): niter-cap hit rate 99.05% -> **0.00%** (not one pixel of
  20,259,001), mean `n_iter_used` 496.72 -> 55.29, total wall time
  ~3h39m -> ~81min (~2.7x faster, entirely as a side effect of CLEAN no
  longer fighting wrong RMSF tables for most of the image), and the
  signal-free-RM-plane finding reversed to the correct direction
  (residual now below dirty at both edge planes, was 1890-2443 uJy vs
  290-361 uJy before the fix).

### Added — RAM-aware, tiled mask handling (ticket T10b)
- The mask cube is now read through the exact same tiled,
  `io_read_threads`-aware mechanism as the AMP/PHA float cubes
  (`read_mask_tile`/`read_mask_chunk`, `FTGSVB`-based -- replacing
  T10a's `read_mask_cube` outright once its only call site was gone),
  rather than held whole in memory regardless of image size. `FTGSVB`
  has no flat pre-multiplied element-count argument at all, so the
  overflow class T10a fixed cannot recur here by construction.
- The mask-pattern -> RMSF-table cache moved from one upfront serial
  pre-scan (needing the whole mask resident first) to an INCREMENTAL
  scheme: scanned once per tile, serially, right before that tile's own
  parallel CLEAN loop starts -- no locks needed at all, since the tile
  loop's own strict sequencing already guarantees a pattern is never
  looked up before it's inserted.
- `plan_rmclean_tile`'s own memory-budget formula gained one more term
  for the mask's own per-pixel footprint (1 byte/voxel x `nchan`), so a
  big machine still gets a tile large enough to recover the old
  whole-cube-resident behaviour for free when it fits, and a small
  machine gets smaller tiles automatically instead of running out of
  memory.
- Verified via the strongest test available: a full real-data
  confirmation run is BYTE-FOR-BYTE IDENTICAL to T10a's own
  confirmation run across all 6 output cubes, with an exactly-matching
  aggregate stop-reason summary -- proof this is a pure architectural
  improvement with zero behavioural change.

### Verification (T7-T10b)
- Full suite green throughout, growing from 121 to 122 (a new sparse-
  fixture overflow regression test for T10a) and back to 121 (that test
  retired once T10b made its own vulnerability class structurally
  impossible, superseded by the existing forced-tile-vs-default test's
  own coverage of the new tiled mask path).
- The `--help`-text-only fix (documenting previously-real-but-
  undocumented `log_level`/`timing_enabled`/`log_output_file`/
  `io_overlap` keys on `reproject_cubes`/`convolve_cubes`/`match_cubes`,
  found during this documentation pass) is unrelated to RM-CLEAN
  specifically but shipped alongside it; full suite unaffected
  (121/121, unchanged).

### Added — Multi-band preprocessing: end-to-end test coverage (ticket T15, `docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md`)
- Closed a real, unaddressed gap found while writing user-facing
  documentation: existing multi-band tests either exercised
  `rm_synthesis`'s own schema/geometry-mismatch REJECTION in isolation,
  or the preprocessing tools' own internal consistency -- none of them
  took genuinely mismatched multi-band data, fixed it with
  `reproject_cubes`/`convolve_cubes`/`match_cubes`, and confirmed
  `rm_synthesis` then recovered the correct answer. Three new
  permanent sections now do exactly that, each isolating one failure
  mode (grid-only, resolution-only, both together via `match_cubes
  stages=both`), each ending in a real `check_rm_peak.py` assertion
  against the known injected sources.
- New dedicated fixture (`tests/make_test_cubes.py`):
  `TEST_BAND2_UNMATCHED.Q/U.FITSCUBE` plus two flat per-channel beam
  logs at genuinely different resolutions (10"/20"), so `convolve_cubes`
  has real smoothing work to do, not a no-op.
- Generalized `strip_fits_ext` (`src/reproject_cubes.f90`,
  `src/convolve_cubes.f90`, `src/match_cubes.f90`) to strip any trailing
  filename extension instead of a hardcoded `.fits`/`.FITS` check. Found
  while building the tests above: this project's own `.FITSCUBE` test
  fixtures weren't recognised as having an extension at all, so output
  filenames doubled up (`name.FITSCUBE_REPROJ.FITS`) instead of swapping
  cleanly (`name_REPROJ.FITS`) -- cosmetic only, never affected real
  users with standard `.fits` inputs, but fixed properly rather than
  left as a quirk. Full suite re-verified at 127/127 after the change.
- Full suite: 127/127 (up from 121).

## [5.0] - `5.0-rc.1` tagged on `develop`; real-scale validation pending before `main`

Multi-band Faraday tomography milestone — by far the largest single body
of work this project has shipped at once. Four parts, each usable on its
own but designed to work as a pipeline: `rm_synthesis` itself can now
merge frequency channels from several input files into one RM synthesis
run; three standalone tools (`reproject_cubes`, `convolve_cubes`, and
`match_cubes`, which consolidates the other two with optional in-memory
chaining) prepare real, mismatched-geometry/mismatched-resolution bands to
actually be combined that way; and beam metadata (`BMAJ`/`BMIN`/`BPA`,
`CASAMBM`/`BEAMS`) is now faithfully carried through the whole chain —
`rm_synthesis`'s own outputs and the entire preprocessing toolchain alike
— rather than silently dropped at any stage. Full design rationale,
decisions recorded with the user, and ticket-by-ticket verification
evidence for all of this lives in
`docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md` (tickets T0-T14) — this entry is
a summary, not a replacement for that record.

### Added — multi-band RM synthesis (`rm_synthesis`, tickets T1-T9)
- Comma-separated-list config schema for every per-band key (`infileQ`,
  `infileU`, `resiQ`/`slopeQ`/`resiU`/`slopeU`, `infileI`/`path_I`,
  `badchan_file`, and — new this release — `chan_blc`/`chan_trc`/
  `chan_inc`): band count is derived from list length, no separate
  `nbands` key. A config with no commas anywhere behaves exactly as
  before — this is additive, not a breaking change to any existing cfg.
- Unified N-band geometry validation (RA/Dec WCS, NAXIS, frequency-axis
  index) against a `reference_band`, loudly refusing before any compute
  on mismatch — the same exact-equality philosophy the existing Q-vs-U
  check already used, generalized rather than replaced.
- Multi-band frequency/λ² merge: every band's channels concatenated into
  one merged spectrum (no deduplication in overlaps, no sort required —
  the DFT kernel is order-independent) and run through the existing
  single-band RM-synthesis compute path unchanged.
- `δRM`/`max RM scale`/per-band `ΔRM` diagnostic, logged (not
  auto-applied) for multi-band runs, since `use_auto_rm_range=1`'s
  existing heuristic is unsafe across bands — verified against a
  thesis-published table (Table 6.1) to within ~1%.
- Multi-tile multi-band runs (previously single-tile only) — verified
  bit-identical to the single-tile result of the same data.
- Per-band channel sub-range selection (`chan_blc`/`chan_trc`/
  `chan_inc`), independent per band — e.g. reject bad edge channels or
  hand-pick a good sub-range per band.
- Per-band bad-channel files — each band flags its own bad channels via
  its own required file (same required-key convention as `infileQ`).
- GPU offload for multi-band (no compute-path changes needed — the same
  kernel already used by CPU and GPU; verified on real GPU hardware,
  including the two-level VRAM staging path).
- `io_read_threads>1`/`io_overlap` enabled for multi-band (previously
  blocked entirely — found by code inspection, on direct challenge, that
  the restriction was unnecessarily conservative, the same pattern as
  several tickets before it).
- A real thesis-scenario regression (`tests/check_thesis_scenario.py`):
  point-source recovery, Faraday-thick component reveal (P+L combined
  ~9x the P-alone peak amplitude in its own RM window), and F2/F3
  resolved-vs-blended behaviour, reproduced from published Table 6.1
  bands (P: 300/30 MHz, L: 1200/120 MHz).

### Added — cross-band preprocessing toolchain (tickets T10-T11)
- `reproject_cubes`: new standalone tool (own binary, `make
  reproject_cubes`) reprojecting two or more FITS cubes onto one common
  sky grid via Starlink AST + `astResampleR`, with three footprint modes
  (`intersection`/`union`/`reference`), full WCS/header propagation
  (including `CROTA`/`PCi_j`/`CDi_j` sky rotation), `mem_frac_ram`-budgeted
  block I/O, and OpenMP parallelism across planes.
- `gaussft_mod` (`src/gaussft.f90`): new pure-computation module for
  elliptical-Gaussian FFT-domain beam-matching convolution (deconvolve
  from a source PSF, reconvolve to a target PSF), thread-safe for OpenMP
  via a plan-once/execute-many split (FFTW's planner is not thread-safe;
  a single plan's "new-array execute" form is, verified directly).
- `commonbeam_mod` (`src/commonbeam.f90`): new module finding the
  smallest common beam every one of N per-channel PSFs can be
  deconvolved from (convex hull + minimum-volume-enclosing-ellipse via
  Khachiyan's algorithm + Sault/MIRIAD deconvolvability validation),
  verified against the `radio_beam` Python package on real ASKAP data.
- `convolve_cubes`: new standalone tool (own binary, `make
  convolve_cubes`) driving `gaussft_mod`/`commonbeam_mod` to convolve
  cubes — across one or several input files together — to one common
  angular resolution. Reads per-channel beams from a CASA-style `BEAMS`
  binary table or a portable ASCII/CSV beam log (`cfg/
  example_beamLog.txt`/`.csv`, ready-to-adapt examples included); a
  channel is bad if missing from the beam source or listed with BMAJ or
  BMIN equal to 0. `max_common_bmaj` guards against silently convolving
  to an unexpectedly coarse auto-derived resolution.

### Added — consolidated in-memory reproject+convolve toolchain (ticket T13)
- `match_cubes`: new standalone tool (own binary, `make match_cubes`)
  consolidating `reproject_cubes` and `convolve_cubes` into one program
  that can run either stage alone, or both CHAINED THROUGH MEMORY with no
  intermediate FITS file (`stages=reproject|convolve|both`,
  `order=convolve_reproject|reproject_convolve`, default
  `convolve_reproject`) — avoids a full extra disk read/write round-trip
  on real 200GB+ cubes. Neither existing standalone tool is modified;
  `match_cubes.f90` adapts their logic instead, a deliberate
  zero-regression-risk tradeoff. Verified via "chaining equivalence"
  (both orders reproduce the two OLD tools' disk-based pipeline
  bit-for-bit), which directly caught two real bugs before release: a
  degenerate-axis header-copy loss, and an FFTW plan-size mismatch in
  `reproject_convolve` order that crashed with heap corruption.

### Added — beam-metadata propagation (tickets T12, T14)
- `rm_synthesis` now propagates `BMAJ`/`BMIN`/`BPA` to all 8 output
  products (previously propagated none at all). If the input has
  `CASAMBM=T` (a genuinely per-channel-varying beam not yet run through
  `convolve_cubes`), the flux-derived outputs (AMP/PHA, and the
  PEAK/RMPEAK/ANGPEAK/SNR maps when `cubestat=y`) additionally get
  `CASAMBM=T` plus the input's own real per-channel `BEAMS` table
  attached as an extension, plus a `HISTORY` note — deliberately not
  MASK/NVALID, which are validity bookkeeping, not flux data. In
  multi-band mode, every band's own beam metadata is now cross-checked
  against the reference band's, with a runtime warning on mismatch.
- `reproject_cubes` (and `match_cubes` with `stages=reproject`) now
  propagate `CASAMBM`/`BEAMS` too: reprojection never touches the beam
  itself, so a genuine per-channel `BEAMS` table on the input is copied
  through to the output unchanged (previously silently dropped — the
  scalar `BMAJ`/`BMIN`/`BPA` alone survived, but not the real per-channel
  table `CASAMBM=T` refers to).
- `convolve_cubes` (and `match_cubes` whenever `convolve` is active) now
  always attach `CASAMBM=T` plus a freshly synthesized `BEAMS` table on
  output, regardless of whether the input had one — one row per channel,
  the common target beam for every channel actually convolved, and the
  same degenerate sentinel CASA itself uses (`tiny(1.0)`, ~1.18e-38) for
  a bad/skipped channel, so a downstream reader can tell exactly which
  channels reached the common resolution rather than a single scalar
  that would misrepresent a bad/NaN channel as sharing it.

### Fixed
- `rm_synthesis` opened its own Q/U/I/mask input cubes `READWRITE`
  despite never writing to any of them (confirmed by grep: no write-type
  CFITSIO call anywhere in the file targets those units), an unnecessary
  risk to irreplaceable input data. Now opened `READONLY`, matching how
  this file's own parallel tile-reader threads for the same files
  already worked.
- `convolve_cubes`' bad-channel detection (both the CASA `BEAMS`-table
  and ASCII/CSV readers) only checked BMAJ for a degenerate (zero) beam
  entry; a channel with BMAJ present but BMIN equal to 0 was silently
  treated as good. Now checks both.
- A latent bug in the per-band channel-count bookkeeping, surfaced (not
  triggered — `subim` was blocked outright for multi-band until the same
  ticket that found it) while adding per-band channel sub-range
  selection: the reference band's own selected-channel count was being
  computed from its raw NAXIS3 rather than its actual selected range.

### Validation
- All 4 build flavours (`scratch/make_all.sh`) clean; full
  `tests/run_tests.sh` 49/49 pass (up from 28 at the start of this
  branch), re-run clean after every change in this release.
- Multi-band `rm_synthesis`: `nbands=1` bit-identical sweep (140/140
  FITS outputs) held after every single ticket in this release, with no
  exceptions — the explicit bar this whole effort was held to throughout.
  Multi-tile-vs-single-tile, per-band-channel-sub-range, and per-band-
  bad-channel-file all verified bit-identical against known-good
  references, not merely "doesn't crash".
- `reproject_cubes`: byte-identical to independently-computed
  (Python/astropy) ground truth at spot-checked pixels; a real
  `FTGSVE` axis-order bug caught by a non-adjacent-sky-axis fixture and
  fixed; 25 repeated stress runs, no failures.
- `gaussft_mod`: identity round-trip and asymmetric-beam cross-check
  against an independent Python implementation; 16-OpenMP-thread
  shared-plan concurrency test bit-identical to serial.
- `commonbeam_mod`: matches `radio_beam` 0.3.9 on a real 286-channel
  ASKAP `BEAMS` table to within 0.003 arcsec (BMAJ/BMIN) and mod-180
  degrees (PA); independently confirmed deconvolvable from all 286 real
  beams.
- `convolve_cubes`: bit-exact identity check (target beam == a
  channel's own native beam reproduces that channel's input exactly,
  validating the SKY-to-PIXEL BPA convention conversion end-to-end);
  smoke-tested against a real cutout of ASKAP data with no NaN/Inf.
- `rm_synthesis` beam propagation: injected real BMAJ/BMIN/BPA and
  `CASAMBM=T`/`BEAMS` cases and confirmed exact propagation to the
  correct output subset only; injected mismatched per-band beams in a
  multi-band run and confirmed the cross-band warning fires correctly
  (and stays silent when bands genuinely match); fed a real
  `convolve_cubes`-produced NaN bad-channel plane into `rm_synthesis`
  with no `badchan_file` and confirmed automatic exclusion via existing
  NaN detection.
- `match_cubes`: single-stage equivalence (`stages=reproject`/`convolve`
  alone byte-identical to the corresponding standalone tool) and
  chaining equivalence (both orders bit-identical to the two old
  standalone tools run back-to-back through a real disk intermediate),
  on a genuine 2-band scenario with offset grids, per-channel beams, a
  real bad channel per band, and an intersection-mode footprint crop.
- `reproject_cubes`/`convolve_cubes`/`match_cubes` `CASAMBM`/`BEAMS`
  propagation: verified against the real 5-column CASA `BEAMS` layout
  (`BMAJ`/`BMIN`/`BPA`/`CHAN`/`POL`) confirmed on a real ASKAP cube;
  `reproject_cubes` output `BEAMS` table byte-identical to the input's
  own; `convolve_cubes` output `BEAMS` table confirmed to carry the
  common target beam on every good channel and the degenerate sentinel
  on the known bad channel, matching the convolved data's own all-NaN
  bad-channel plane; `match_cubes` confirmed identical to both
  standalone tools and to the two-step disk pipeline in both chain
  orders.

### Not yet done
- A full run of the preprocessing toolchain against the complete 23GB
  real ASKAP cube this work targets (only cutouts and synthetic data
  verified so far) — required before an actual `main` release.
- T13 (`match_cubes`) and T14 (`CASAMBM`/`BEAMS` propagation for the
  toolchain) landed on `multi-band-tomography` after `5.0-rc.1` was
  tagged and `develop` merged; not yet merged back to `develop` or
  re-tagged.

## [4.1] - 2026-07-20

Diagnostics milestone: closes a real gap in the run-timing picture found
while looking into why a real large-cube run's swim-lane plot showed an
oddly long first read and first write, and drops a second build path that
had quietly stopped working.

### Added
- Two new timed stages, `io_read_init` and `io_write_init`, cover work
  that always ran but was never counted anywhere: reading the input
  cubes' true dimensions before anything else can happen, and (when
  `io_write_threads>1`) staking out the output files' full size on disk
  the first time they're closed. Both now get their own timer, their own
  line in the run log (using the same read/write categories and colours
  the swim-lane plot already draws, so existing plots need no changes),
  and their own `io_read_init_sec`/`io_write_init_sec` columns in the
  timing CSV — with real byte counts, not placeholders: true FITS header
  bytes read (via `FTGHSP`) for the read side, and the true declared
  AMP+PHA size for the write side.
- The printed "Timing summary" and its "Macro timing breakdown" now fold
  these two stages into the read/write totals, so the reported total
  finally reconciles against the sum of its own named stages instead of
  leaving an unexplained gap.

### Removed
- A second, CMake-based way of building the project (`CMakeLists.txt`,
  `cmake_build.sh`) that had drifted out of sync with how the code is
  actually built (it never learned the preprocessor switch the real
  build now depends on) and no longer worked at all.
- Three duplicate fixed-form source files (`src/myfits_info.f`,
  `src/printerror.f`, `src/rm_synthesis.f`) left over from before the
  current free-form build path existed; only their `.f90` counterparts
  were ever actually compiled.

### Validation
- Clean 4-variant rebuild (0 new warnings), full 28/28 test suite,
  including the two runs that specifically re-check bit-identical output
  under `io_overlap` and `io_write_threads>1` — the exact code path the
  new write-init timing sits next to.
- New log markers hand-checked on a real run: they appear first, in
  correct time order, ahead of tile 1's own read/write, with byte counts
  that check out exactly against the cube's declared dimensions.

## [4.0] - 2026-07-19

Maintainability and documentation milestone rather than a new-capability
one: the two fortran files carrying almost all of rmtool's logic were
restructured around derived types with zero change in observable
behaviour, a real swim-lane plotter bug was found and fixed, and the
README gained a new "Motivation" section explaining rmtool's parallelism
model for a research (not HPC-expert) audience. 

### Added
- README "Motivation" section: explains tiling for memory, the
  disk-layout-driven tile shape, the CPU-side "corner turn" (and why the
  GPU path skips it), read/write/compute overlap, parallel I/O channels,
  and GPU offload — in plain language.

### Changed
- All ~56 config values (paths, tiling/memory-planning knobs, RM sampling,
  masking, GPU, I/O-parallelism, logging/timing keys) are now bundled into
  one `rmsynth_config_t` derived type and read directly as `cfg%field` at
  every use site, replacing a flat scope of loose local variables that
  `read_cfg_keyval` used to populate through a ~56-argument signature (now
  just `(cfgfile, cfg, status)`). One value, one place it lives, instead
  of a config-parsed copy and a separately-mutated local that could in
  principle drift apart.
- The RAM/VRAM tile-size planner (auto-tiling policy, safety-shrink loop,
  VRAM sub-block sizing) is now a named `plan_tile` routine operating on a
  `tile_plan_t` bundle, instead of ~150 lines of inline arithmetic in the
  middle of the main tile loop.
- The per-tile write dispatch's field assembly is now a named
  `populate_write_job` call; the synchronisation code around it (the
  join-before-reuse/join-before-dispatch invariants documented in
  `docs/user/ARCHITECTURE.md`, the two safeguards behind a real historical
  production SIGSEGV) is untouched, in its exact original order.
- Read-side byte-count and thread-split arithmetic in the tile loop moved
  into two small named helpers (`compute_tile_read_bytes`,
  `split_channels_across_threads`).

### Validation
- Every ticket independently gated: clean 4-variant rebuild (0 errors, 0
  new warnings), 28/28 tests, and a bit-identical sweep of all 140
  archived FITS outputs against a frozen baseline.
- Additionally validated on a real production-scale run (`io_overlap=y`,
  `io_read_threads=4`, `io_write_threads=2`, real ASKAP-style cube) —
  confirmed data integrity and passed the structural
  no-overlapping-tile-writes check, the same invariant the historical
  postmortem in `docs/user/ARCHITECTURE.md` was written to guard.

### Fixed
- Swim-lane plotter (`scripts/plot_tile_async_swimlane.py`): the
  legend/info side panel could silently vanish entirely on plots with
  enough content (e.g. a run with both GPU staging slots and a full
  info block) that no candidate layout in a fixed search grid happened
  to fit — every artist got removed and nothing was drawn, with no error.
  A related case rendered both boxes but let their rounded corners
  visually overlap by a few pixels, because the layout check measured
  the info box's raw text extent and missed the padded border box drawn
  around it. Replaced the search entirely with a deterministic layout:
  the side panel is split into a fixed bottom region (legend, narrow
  column count so it can't dwarf the info box) and a top region (info
  box, continuously font-sized to fill whatever the legend left behind).
  The two regions are stacked, not independently placed, so they cannot
  collide by construction, and the info box now uses the space it's
  given instead of stopping at an arbitrary integer point size.

## [3.0] - 2026-07-18

IO-efficiency milestone: parallel reads, genuine parallel writes, async
tile-write overlap, and the crash/correctness work that came with
building them. All planned tickets (T0-T6) are done and validated
end-to-end on real Setonix production hardware, including T6's actual
write-throughput gain (see Validation below). See
`docs/dev/ARCHIVED/RELEASE_NOTES_3.0.md` for the full writeup.

### Added
- `io_read_threads` cfg key: N independent read-only CFITSIO handles per
  input cube, reading disjoint channel ranges concurrently.
- `io_write_threads` cfg key: N independent Fortran STREAM I/O units
  write disjoint RM-bin byte ranges of the AMP/PHA output cubes directly,
  bypassing CFITSIO's `ftpsse`/handle machinery for pixel data entirely
  (T6) -- see Fixed/Known limitations below for the two-stage history of
  why this replaced an earlier, unsafe design rather than being the
  first approach shipped.
- `io_overlap` cfg key: tile N's write runs on a background POSIX thread
  concurrent with tile N+1's read/mask/prep/compute/cubestat. Uses a raw
  pthread rather than an OpenMP task specifically so it cannot silently
  nest/collapse the existing OpenMP parallel regions used by
  `io_read_threads` and the compute kernel.
- Auto-tiler RAM planning is aware of `io_overlap`'s doubled output-side
  buffers, so `tile_dec` is planned smaller automatically under the same
  `mem_frac_ram` -- no new user-facing memory configuration.
- Swim-lane plotter: I/O read and I/O write render as separate lanes
  (previously one shared lane distinguished only by colour), and every
  plot now includes a stage-time-totals bar panel (seconds and % of wall
  time per stage, largest first).
- `tile_read`/`tile_write` log lines now carry a `bytes=<N>` field, and
  the swim-lane plotter renders a new I/O throughput (MB/s) panel from
  it -- stacked directly below the swim-lane/thread panel, sharing its
  time axis so a dip or spike lines up with the Gantt bar above it.
  Absent (no empty panel) for logs predating this field.
- New regression tests: bit-identical `io_overlap=n` vs `y` comparison, a
  structural "no two tile writes ever overlap" invariant check, and a
  bit-identical `io_write_threads=1` vs `=4` comparison across all 8
  output products (`tests/run_tests.sh` §13-14).
- Full, sectioned cfg reference in `README.md`: every key the parser
  accepts, marked required/required-if/optional with its real default,
  cross-checked against `read_cfg_keyval`'s case statements rather than
  written from memory.

### Fixed
- int64-safe flattened tile indices throughout the tile/scatter/mask
  loops (`ipix_tile`, `pix_base`, `src_idx`, `dst_idx`, `ipix_full`,
  `ipix_sub`, `idx_wts`, and the equivalents in `prepare_cpu_data`/
  `prepare_gpu_data`/`tile_extract_gpu_rm_blocked`/
  `cubestat_tail_quantile_maps`) -- same INT32_MAX overflow class as the
  2.0-era `allocate()` fix, previously missed in the runtime index
  arithmetic that actually reads/writes those buffers.
- `io_write_threads>1` root-caused as unsafe (original design): CFITSIO
  aliases repeat read-write opens of an already-open file onto one shared
  internal buffer (`fits_already_open()`), so the "independent" handles
  corrupt each other under concurrent `ftpsse` calls. Produced a real
  SIGSEGV in testing. Hard-clamped to 1 as an interim measure, then fixed
  permanently by replacing the mechanism entirely (see T6 below).
- `io_overlap` write-vs-write race: an initial version only guarded
  buffer-slot reuse (two tiles apart), not whether the *immediately
  preceding* write (a different slot) had finished -- both share the same
  single FITS handle, so two pthreads could call `ftpsse` on it
  concurrently. Produced a real SIGSEGV on a production-scale run (a
  small leftover tile at the bottom of an image whose height wasn't an
  exact multiple of the tile size raced ahead of the previous tile's
  write). Fixed by unconditionally joining any outstanding write before
  dispatching a new one.
- **T6: genuine write-throughput parallelism.** `io_write_threads>1` no
  longer hard-clamped -- each RM-chunk now writes via an independent
  Fortran STREAM I/O unit directly to its byte offset (from `FTGHAD`),
  bypassing CFITSIO for AMP/PHA pixel data entirely, relying only on the
  POSIX guarantee that concurrent writes to disjoint byte ranges of one
  file are safe. Found and fixed a second bug during this work, before it
  ever shipped: leaving CFITSIO's handle for these files open until
  program exit (as the serial path always had) caused CFITSIO's own
  data-fill-check, at final close, to treat the raw-written pixel data as
  "past its own stale end-of-file" bookkeeping (never updated, since the
  raw writer bypassed CFITSIO) and zero-fill over it -- silent data loss,
  not a crash, and only visible in the final on-disk state after close.
  Fixed by closing CFITSIO's handle for AMP/PHA immediately after
  fetching the byte offset it provides, before any raw write happens.
- Swim-lane plotter: the "CPU stage" row (CPU thread-detail view) was
  silently missing its compute segment -- filtered out on the assumption
  that the per-thread lanes above already covered it, which left the row
  summing to less than the tile's actual non-I/O time and a legend entry
  ("CPU compute") that never had a corresponding bar. Restored; the two
  views aren't redundant (the stage row shows *when* the stage ran as a
  whole, the thread lanes show *how* it was parallelised).
- Swim-lane plotter: the GPU pipeline view's synchronous-fallback path
  (a tile that fits in one VRAM sub-block, so there's no async
  double-buffering) expected an old `send N/M` log format the Fortran
  code no longer emits -- current single-shot GPU compute logs plain
  `gpu send`/`gpu recv` notes instead, so non-staged GPU runs were
  silently rendering with no GPU lane at all. Fixed to recognize the
  current format.

### Validation
- Full build matrix remained successful (`OMP/GPU` combinations, zero
  compiler warnings).
- Test suite green (`28/28`), including a full bit-for-bit
  `io_write_threads=1` vs `=4` output comparison and manual verification
  that `io_write_threads>1` combined with `io_overlap=y` also produces
  bit-identical output.
- End-to-end production-scale validation on real Setonix hardware and
  ASKAP/EMU data (13308x11870, 288 channels): the exact case that
  originally crashed now completes without error, `io_read_threads`/
  `io_overlap` confirmed on the real workload, and T6's write-throughput
  gain measured directly -- `io_write_threads=8` dropped write from
  2479.9s (96% of wall time) to 108.3s (6%), a ~23x reduction, taking
  total wall time from 2586.7s to 1945.4s (~25% faster end-to-end). Full
  before/after table in `docs/dev/ARCHIVED/RELEASE_NOTES_3.0.md`.

## [2.0] - 2026-07-17

### Added
- Formalized release-cycle documentation:
  - Added this `docs/dev/ARCHIVED/CHANGELOG.md`.
  - Added `docs/dev/ARCHIVED/RELEASE_NOTES_2.0.md`.

### Changed
- Tile planner memory accounting split in `src/rm_synthesis.f90`:
  - Host RAM tile budget now uses `bytes_per_tile_pixel_ram`.
  - GPU VRAM sub-block budget now uses `bytes_per_vram_pixel`.
- Updated docs to reflect planner behaviour and measured outcomes:
  - CPU full-image benchmark improved after planner split.
  - GPU path remained correct but showed a slight regression on the tested environment.
- Extended swim-lane interpretation notes for CPU thread-detail view:
  - Clarified that single-RM-chunk runs (`nrm_out <= nrm_block_size`) show only odd/non-hatched `cpu_extract` traces.

### Validation
- Full build matrix remained successful (`OMP/GPU` combinations).
- Test suite remained green (`22/22`) during these updates.

## [1.1] - 2026-07-16

### Added
- Expanded observability and timeline diagnostics:
  - Structured stage/tile logging enhancements.
  - Swim-lane plotting workflow and example artifacts.

### Changed
- Host OpenMP performance improvements:
  - Staged gather/scatter loop parallelization in `src/rm_synthesis.f90`.
  - Pack/copy parallelization in `src/rm_synthesis_mod.f90` with nested-region guard (`omp_in_parallel`).
- Documentation refresh across build, parallelism, and design notes.

### Fixed
- Docker tag-build helper compatibility for shell invocation:
  - `docker/build_push_from_tag.sh` now re-execs under bash when invoked with `sh`.
