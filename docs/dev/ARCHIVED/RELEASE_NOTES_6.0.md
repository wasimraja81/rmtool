# Release Notes 6.0

Status: **T0-T10b done on the `rmclean-integration` branch; T11 (CLEAN
divergence/noise-wasting-compute stopping criteria) planned, not yet
started. Not yet merged to `develop`/`main`, and no git tag has been cut
— this document is prepared ahead of that decision, not a record of one
already made.** See `docs/dev/ARCHIVED/CHANGELOG.md`'s own `[Unreleased] - RM-CLEAN
integration` entry and `docs/dev/RMCLEAN_INTEGRATION_PLAN.md` (the full
ticket-by-ticket record, including every dead end and how it was found)
for the living record this document summarizes.

## Summary

6.0 makes RM-CLEAN a real, production-usable part of this package for
the first time. RM synthesis alone leaves a dirty Faraday dispersion
function — sidelobes from the RMSF (Rotation Measure Spread Function)
smeared across every real polarized source. RM-CLEAN is the standard
deconvolution step that removes them, and until this release it existed
only as thesis-era code never wired up against real FITS cubes. This
release ships a standalone tool, `rmclean_cubes`, that drives the
existing pure-computation CLEAN core against the real dirty AMP/PHA
cubes `rm_synthesis` already writes — and, over the course of building
and stress-testing it against a real, moderately large-scale ~46GB
ASKAP test dataset, surfaced and fixed a genuinely serious bug: a
silent 32-bit integer overflow that was corrupting CLEAN's own input
data for any sufficiently large real image, unnoticed until a real
full-cube run at that scale exposed it.

Every numeric claim in this document was independently verified against
either a real dataset or an analytically-derived reference, not
asserted — see Validation below, and in particular the standing project
rule this release reinforced directly: never assert a
statistical/chance claim, or declare a bug fixed, without actually
computing the comparison.

## Highlights

### `rmclean_cubes` — standalone RM-CLEAN against real dirty cubes

A new standalone program (own binary, `make rmclean_cubes`, own
`key=value`/`--config` parser and `--help`, following
`reproject_cubes`/`convolve_cubes`/`match_cubes`'s own conventions)
drives `rmclean_mod` (`src/rmclean.f90`, pure Högbom-style complex CLEAN
+ Gaussian restore, no FITS I/O of its own) against a real dirty
AMP/PHA cube pair plus its matching mask cube. Since this tool cannot
resample the RM axis — it's fixed by whatever `rm_synthesis` already
wrote — it validates the existing grid instead (Gate 0: refuses to
proceed if the grid doesn't resolve the RMSF FWHM at
`min_samples_per_fwhm` samples) rather than silently CLEANing on an
undersampled one. Sub-pixel peak location is found via a tiered
matched-filter refinement against the analytically-known RMSF model
(a cheap fixed-location fit, escalating to a full local search only
when the misfit exceeds a data-driven threshold) rather than
interpolating the stored samples — this means the dirty cube only ever
needs ordinary resolution-level sampling, not the much finer,
carrier-rate sampling an earlier, since-superseded design would have
required. Parallelised one OpenMP thread per pixel, with a
mask-pattern cache sharing one RMSF table across every pixel that
shares the same valid-channel set.

### Three independent, combinable CLEAN stopping criteria

CLEAN stops on whichever of `niter` (a hard iteration cap, always
active), `abs_flux_floor` (stop at a literal fixed flux value — native
units or `Jy`/`mJy`/`uJy`), and `auto_nsigma` (stop at `n`× that same
pixel's own noise sigma, estimated once per pixel from its own full
dirty amplitude spectrum) fires first. This replaced an earlier
`threshold=`/`threshold_snr=` design that turned out, on inspection
prompted by a direct user question ("do we not have a way to stop at
absolute flux? are you sure?"), to funnel both modes into one
overloaded, always-multiplicative internal variable — silently breaking
the absolute-flux mode entirely and making the SNR mode a
self-referential moving target (compared against the CURRENT, shrinking
residual's own noise every iteration, not a fixed floor). Re-running a
real subimage with the old production SNR value under the new design
took the tail residual RMS from 6.89 µJy (3.7x *below* the true ~25.5
µJy noise floor — i.e. actively over-cleaning) to 25.24 µJy (matching
it), converging in a mean of 3.89 iterations instead of 233.7. The
per-pixel noise-sigma estimator itself went through two further,
data-driven revisions before landing on the simplest option that
actually measured best against ten real pixels' own independently-known
true noise: the pixel's own full dirty amplitude spectrum's
interquartile range, converted to sigma via the analytic (untruncated)
Rayleigh-distribution relation.

### A real data-corruption bug, found by a real full-cube run

A complete, real 4501×4501-pixel CLEAN run of a ~46GB ASKAP test dataset
produced a result sharply inconsistent with every small-scale test that
preceded it: 99.05% of all 20,259,001 pixels exhausted the full 500-
iteration budget (mean 496.72 iterations used), versus 0% in every
prior subimage test. Investigated directly rather than dismissed —
bit-for-bit comparison of the real dirty spectrum fed into CLEAN for a
specific pixel, both bit-identical against a small cutout of the same
data — isolated the fault to `read_mask_cube`'s own read of the mask
cube: a single CFITSIO call whose element count was computed as
`nx*ny*nchan` in ordinary 32-bit integer arithmetic, silently
overflowing (wrapping, not erroring) for any mask cube exceeding 2^31
elements. This real dataset's own 4501×4501×288 mask cube is 2.7x over
that limit, silently truncating the valid-channel read at ~76 of 288
channels for **every pixel in the image** — corrupting the RMSF table
and derotation used to CLEAN nearly the entire dataset. This single bug
turned out to be the root cause of three things investigated this
release as apparently separate mysteries: the niter-exhaustion rate
above, a specific pixel that appeared to never converge in isolation,
and spurious CLEAN components (with residual amplitude exceeding the
original dirty amplitude) at RM planes that should have been
signal-free. Fixed by switching to `FTGPVBLL`, CFITSIO's own genuinely
64-bit entry point for the identical underlying C call — confirmed
directly against this project's own bundled CFITSIO source, not
assumed.

### RAM-aware, tiled mask handling

Fixing the immediate bug raised a sharper question: could an
arbitrarily large `mem_frac_ram`/tile-size choice, on a big enough
machine or dataset, still ask a single CFITSIO call to move more
elements than a 32-bit count can hold, even after the fix? Answered
directly rather than assumed safe: `FTGSVB`/`FTGSVE` (the tiled
subsection-read entry points already used for the AMP/PHA float cubes)
have no equivalent exposure at any size — confirmed by reading their
own C implementation, which accumulates the total element count in
genuine 64-bit arithmetic internally. The mask cube's read was
therefore unified with the exact same tiled, `io_read_threads`-aware
mechanism the AMP/PHA cubes already used, rather than patching the old
whole-cube-resident design in place. This also answered a second
question the user posed directly: why did the mask cube need a
separate whole-cube-resident code path at all, when the float cubes
were already tiled regardless of image size? The answer was the
mask-pattern cache's own one-time global pre-scan — reworked to build
INCREMENTALLY instead, once per tile, right before that tile's own
parallel CLEAN loop starts, with no locking needed at all (the tile
loop's own strict sequencing already guarantees a pattern is never
looked up before it's cached). A big machine still gets a tile large
enough to recover the old whole-cube-resident behaviour for free when
it fits the budget; a small machine gets smaller tiles automatically
instead of running out of memory.

## Validation

- Full regression suite green throughout this release: 85/85 after the
  core algorithm/tool/peak-refinement tickets (T1/T2/T2b/T3), 93/93
  after the memory-budgeted tiled I/O port (T4), growing to 122/122
  with a (since-retired) sparse-fixture overflow test added for the
  mask-read bug (T10a), and 121/121 at the point this document was
  written (T10b made that specific test's own vulnerability class
  structurally impossible, so it was retired rather than kept as dead
  weight).
- `read_mask_tile` (the new tiled mask reader) verified directly against
  an independent numpy/astropy reference: a real FITS fixture, reading
  an *offset* 2×2 subregion (not the whole file) split across 3
  `io_read_threads` workers, matched exactly for all 4 sub-pixels.
- The mask-read overflow fix verified two ways: (1) a purpose-built
  sparse-file regression test (an OS-level sparse ~2.4GB fixture — only
  3 bytes actually written — proving the fix reads all the way to the
  last of ~2.4 billion elements, and that it genuinely fails against
  the old code, confirmed by temporarily reverting it); (2) a full
  real-data confirmation run: niter-exhaustion rate 99.05% → **0.00%**
  (not one pixel of 20,259,001), mean iterations used 496.72 → 55.29,
  total wall time ~3h39m → ~81min (~2.7x faster, entirely a side effect
  of CLEAN no longer fighting corrupted RMSF tables for most of the
  image), and the signal-free-RM-plane finding reversed to the correct
  direction (residual now below dirty amplitude at both edge planes,
  was 1890–2443 µJy vs. 290–361 µJy before the fix).
- The RAM-aware retiling (T10b) verified against the strongest test
  available for a change with no *intended* behavioural effect: a
  complete fresh real-data run with the new tiled mask design is
  BYTE-FOR-BYTE IDENTICAL to the T10a confirmation run above across all
  6 output cubes, with an exactly-matching aggregate stop-reason
  summary.
- Every quantitative claim about "expected" or "surprising" statistics
  in this release was independently computed or simulated, never
  asserted from intuition — a standing project discipline reinforced
  directly during this release after an early, corrected instance of
  the opposite.

## Compatibility and behaviour notes

- **Old RM-CLEAN cfg keys are retired, not aliased.** `threshold=`,
  `threshold_snr=`, `noise_percentile=`, `noise_nlos=`, and
  `noise_seed=` have all been removed from `rmclean_cubes`. A cfg file
  using any of them will now fail to parse (a clear "unknown key" error)
  rather than silently misbehaving — replace with `abs_flux_floor=`/
  `auto_nsigma=`. This is a deliberate choice: the old design's own
  correctness bugs mean a silently-accepted old key would be actively
  misleading, not just outdated.
- Every other tool (`rm_synthesis`, `reproject_cubes`, `convolve_cubes`,
  `match_cubes`) is unaffected — no cfg keys changed, removed, or
  renamed for any of them. Three of them (`reproject_cubes`,
  `convolve_cubes`, `match_cubes`) had their own `--help` text extended
  to document `log_level`/`timing_enabled`/`log_output_file`/
  `io_overlap` keys that already existed and worked but weren't
  previously listed — a documentation completeness fix, not a new
  feature or behaviour change.
- Anyone with an existing dirty AMP/PHA/mask cube from a prior
  `rm_synthesis` run can CLEAN it directly with the new
  `rmclean_cubes` — no re-synthesis needed, and no header/format
  changes to the dirty cubes themselves.
- Output filenames are unchanged: `<outfile>.CLEAN/.RESID/.RESTORED.
  AMP/PHA.RMCUBE.FITS`.

## What shipped in this release

- `src/rmclean.f90` (`rmclean_mod`, pure CLEAN computation) and
  `src/rmclean_cubes.f90` (the standalone tool) — CLEAN algorithm,
  Gate 0 grid validation, tiered matched-filter peak refinement, the
  three-criteria stopping design, memory-budgeted/tiled/threaded I/O.
- `src/rmclean_io_mod.f90` — mask-cube I/O extracted into its own
  module (`read_mask_tile`/`read_mask_chunk`), both for independent
  testability and as the vehicle for the T10a/T10b fixes.
- `cfg/rmclean-example.cfg`, `cfg/rmclean-jennifer.e2e.cfg`,
  `cfg/rmclean-e2e-smalltest.cfg` — annotated, current example configs.
- `tests/data/rmclean_pixel65_65/` — a small, permanent, fast-running
  regression fixture built directly from the real dataset that exposed
  the mask-read bug.
- Full documentation pass: `README.md`, `QUICKSTART.md`,
  `docs/user/ARCHITECTURE.md`, `docs/user/PARALLELISM.md`, `docs/dev/ARCHIVED/CHANGELOG.md`, and
  three new docs — `docs/user/APP_REFERENCE.md` (complete parameter
  reference for all 5 tools), `docs/user/TUTORIAL.md` (a step-by-step
  walkthrough), and `docs/user/EXAMPLES.md` (a scenario cookbook: multi-band
  with matched/mismatched grid or resolution, choosing RM-CLEAN
  stopping criteria, memory/IO tuning, GPU vs. CPU, subimage
  extraction) — RM-CLEAN had previously been documented only in the
  dev-facing planning ticket doc, never brought into the user-facing
  docs.

## What's next (beyond this release)

- Ticket T11 (CLEAN divergence / noise-wasting-compute stopping
  criteria): detecting a peak residual that grows without bound
  (divergence) or plateaus/oscillates without real progress
  (stagnation/wasted compute) — neither of today's three stopping
  criteria catches either case, since `abs_flux_floor`/`auto_nsigma`
  only fire when the residual drops low enough, never when it grows or
  stalls. Deliberately deferred, not dropped: the specific symptom that
  originally motivated this ticket (a pixel appearing to never
  converge) turned out to be fully explained by the T10a mask-read bug,
  so the general case remains worth designing for, just no longer
  urgent.
- Merging this branch to `develop`/`main` and cutting an actual tag —
  not yet done as of this document.
- GPU support for `rmclean_cubes` remains explicitly out of scope for
  this release, deferred to a later, separate effort.
