# R1.0 — Confluent Brahmaputra

**Status: first public release of this package.**

## What this package does

This package turns multi-frequency radio polarization (Stokes Q/U)
cubes into Faraday rotation measure (RM) maps, and cleans them. Three
things, usable independently or as one pipeline:

- **RM synthesis** (`rm_synthesis`): the core Fourier-domain transform
  from frequency/λ² space to Faraday depth space, producing a dirty RM
  cube (amplitude + phase) from one or more Stokes Q/U input cubes.
- **Multi-band preprocessing** (`reproject_cubes`, `convolve_cubes`,
  `match_cubes`): prepares bands observed at different pointings and/or
  different angular resolutions so they can genuinely be combined —
  matching sky grid, matching resolution, or both at once — before
  `rm_synthesis` merges their frequency channels into one RM synthesis
  run.
- **RM-CLEAN** (`rmclean_cubes`): deconvolves the RMSF (Rotation
  Measure Spread Function) sidelobes out of a dirty RM cube, the
  standard cleanup step between raw RM synthesis and a publication-ready
  result.

Every tool is built to run unmodified from a laptop to an HPC node: all
five share the same memory-budgeting philosophy (`mem_frac_ram`) —
process data in RAM-budgeted chunks rather than requiring the whole
cube to fit in memory at once, so a run scales down to fit a small
machine's RAM and scales up to make full use of a large one, with no
code change either way.

## Highlights

### Single-band RM synthesis

The core transform: dirty AMP/PHA RM cubes from Stokes Q/U spectra,
with memory-budgeted tiled I/O, OpenMP and optional GPU offload,
configurable RM sampling and channel handling, and per-tile
read/write I/O parallelism tunable independently of compute
parallelism.

### Multi-band Faraday tomography

`rm_synthesis` merges frequency channels from several input bands into
one RM synthesis run — more channels means finer Faraday-depth
resolution and better sensitivity to complex (multi-component) Faraday
structure. Getting there usually needs preparation first, since
different bands are rarely observed on the same sky grid at the same
resolution:

- `reproject_cubes` resamples mismatched bands onto one common sky
  grid (Starlink AST).
- `convolve_cubes` convolves mismatched bands to one common angular
  resolution (auto-derives the smallest common beam every input can
  actually be deconvolved from, or accepts an explicit target).
- `match_cubes` chains both through memory when a band needs both
  fixes — no intermediate FITS file written to disk, which matters at
  real multi-hundred-GB cube sizes.
- Beam metadata (`BMAJ`/`BMIN`/`BPA`) is carried faithfully through the
  whole chain, not silently dropped at any stage.

### RM-CLEAN deconvolution

`rmclean_cubes` drives a Högbom-style complex CLEAN core against the
real dirty AMP/PHA cubes `rm_synthesis` writes. CLEAN stops per pixel
on whichever of three independent, combinable criteria fires first:
`niter` (a hard iteration cap), `abs_flux_floor` (a fixed flux value,
native units or `Jy`/`mJy`/`uJy`), or `auto_nsigma` (n× that pixel's
own estimated noise sigma). Sub-pixel peak location uses a tiered
matched-filter refinement against the analytically-known RMSF model,
so the dirty cube only needs ordinary resolution-level sampling.
Building and stress-testing this against real data surfaced and fixed
a genuine data-corruption bug (see Validation) — not a hypothetical
edge case, something that was silently producing wrong answers at
scale until it was found.

## Validation

The two feature areas have different validation stories, stated
plainly rather than blended into one blanket claim:

- **RM-CLEAN has been validated at real, meaningful scale.** Beyond
  the full regression suite, it was stress-tested against a real,
  moderately large-scale ~46GB ASKAP test dataset — a complete
  4501×4501-pixel run (20,259,001 pixels). That run surfaced a real
  32-bit integer overflow silently corrupting CLEAN's mask-cube input
  for any sufficiently large image. Fixed and reconfirmed on the same
  real data: the fraction of pixels wrongly exhausting their full
  iteration budget dropped from 99.05% to 0.00%, and a follow-up
  memory-management redesign was confirmed byte-for-byte identical to
  that fix on a full repeat run.
- **Multi-band preprocessing is proven correct end-to-end, but only
  against synthetic test data so far.** A dedicated test fixture with
  deliberately mismatched sky grid and resolution is run through
  `reproject_cubes`/`convolve_cubes`/`match_cubes`, and the known
  injected sources are confirmed recovered at the correct RM afterward
  — this exercises the real toolchain against a real mismatch, not a
  no-op. What it has **not** yet done is process genuinely large real
  multi-band data end-to-end. That run is planned as a near-term
  follow-up (real ASKAP multi-band survey cubes, several hundred GB
  total), deferred so far purely because of the disk space a real
  output run needs — tracked as a known next step, not silently
  skipped.
- Full automated regression suite: 127/127, covering both feature
  areas plus every tool's own parameter handling, I/O parallelism
  consistency, and GPU/CPU numerical agreement.

## Getting started

- [../../QUICKSTART.md](../../QUICKSTART.md) — build the tools, quick
  reference for every tool.
- [TUTORIAL.md](TUTORIAL.md) — a full walkthrough from a fresh checkout
  to a cleaned RM cube, every command actually run and verified.
- [EXAMPLES.md](EXAMPLES.md) — recipes for real scenarios: single-band
  vs. multi-band, choosing RM-CLEAN stopping criteria, memory/IO
  tuning, GPU vs. CPU.
- [APP_REFERENCE.md](APP_REFERENCE.md) — complete parameter reference
  for all five tools.

## What's next

- Real-scale validation of the multi-band preprocessing toolchain
  against genuine, large ASKAP multi-band survey data (see Validation
  above) — the immediate next priority.
- Further CLEAN stopping-criteria work: detecting a peak residual that
  diverges or stalls without real progress, a case none of today's
  three stopping criteria catches on its own.
- GPU support for `rmclean_cubes` (currently CPU-only; `rm_synthesis`
  already supports GPU offload).
