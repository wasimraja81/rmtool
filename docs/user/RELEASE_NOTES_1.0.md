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
  `match_cubes`): matches sky grid and/or angular resolution before
  `rm_synthesis` merges frequency channels into one run. Grid-matching
  is a cross-band concern; resolution-matching isn't — a single band's
  own channels commonly need it too, since native beam size varies
  with frequency.
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
- `convolve_cubes` convolves channels to one common angular resolution
  (auto-derives the smallest common beam every input can actually be
  deconvolved from, or accepts an explicit target). This applies
  within a single band too — its own channels' native beam varies with
  frequency — not only when combining bands.
- `match_cubes` chains both through memory when a band needs both
  fixes — no intermediate FITS file written to disk, which matters at
  real multi-hundred-GB cube sizes.
- Beam metadata (`BMAJ`/`BMIN`/`BPA`) is carried faithfully through the
  whole chain, not silently dropped at any stage.
- Both the resolution-matching step (convolution) and the grid-matching
  step (resampling onto a common sky grid) process one image plane at a
  time, using several CPU cores to speed up that single plane rather
  than handing whole planes to separate cores in parallel — so peak
  memory use no longer grows with how many cores a run is given. That's
  what makes a real multi-hundred-GB survey run practical on ordinary
  hardware, and keeps the same code path usable on much smaller
  machines too, down to something like a Raspberry Pi.
- When combining bands that don't share the same sky pointing,
  `match_cubes` tells you if the band you picked to define the output
  grid isn't the most efficient choice, and explains why — matching a
  band onto its own kind of pointing is fast, matching it onto a
  genuinely different one is much slower, so the choice matters.
  Nothing is changed automatically; you decide.
- Choosing the reference band well genuinely matters: measured on a
  genuine two-survey dataset (a 144-channel, 11559×11655-pixel band and
  a 288-channel, 14300×12395-pixel band, different pointings and
  resolutions, ~563GB combined), the grid-matching step alone took
  about 3h40m with a poorly-chosen reference band versus about 1h53m
  with a well-chosen one — nearly half the time, from this one choice
  alone.

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
scale until it was found. A separate real-data finding fixed a large,
avoidable delay at the very start of a run: creating each output cube
was wasting a long stretch of time writing empty placeholder space
that the real cleaned results were always going to overwrite anyway —
for the largest real cubes tested, this alone was costing over half an
hour of otherwise-idle run time before any actual cleaning could even
begin. Output now starts flowing essentially immediately instead.

## Validation

The two feature areas have different validation stories, stated
plainly rather than blended into one blanket claim:

- **RM-CLEAN has been validated at real, meaningful scale, twice
  over.** Beyond the full regression suite, it was first stress-tested
  against a real, moderately large-scale ~46GB ASKAP test dataset — a
  complete 4501×4501-pixel run (20,259,001 pixels). That run surfaced a
  real 32-bit integer overflow silently corrupting CLEAN's mask-cube
  input for any sufficiently large image. Fixed and reconfirmed on the
  same real data: the fraction of pixels wrongly exhausting their full
  iteration budget dropped from 99.05% to 0.00%, and a follow-up
  memory-management redesign was confirmed byte-for-byte identical to
  that fix on a full repeat run. More recently, RM-CLEAN was also run to
  completion against the much larger genuine multi-survey dataset
  described below — an 11715×6664-pixel cube, 59,334,010 pixels
  actually cleaned, no errors.
- **Multi-band preprocessing is proven correct end-to-end against
  synthetic test data, and has now completed a genuine real-scale run
  through the entire pipeline, start to finish.** A dedicated test
  fixture with deliberately mismatched sky grid and resolution is run
  through `reproject_cubes`/`convolve_cubes`/`match_cubes`, and the
  known injected sources are confirmed recovered at the correct RM
  afterward — this exercises the real toolchain against a real
  mismatch, not a no-op. Beyond that, the full chain has now also been
  run on genuine ASKAP survey data: two real bands — a 144-channel,
  11559×11655-pixel WALLABY cube and a 288-channel, 14300×12395-pixel
  EMU cube, ~563GB of Stokes Q/U input between them — with different
  sky footprints and different resolutions, matched onto one common
  grid/beam by `match_cubes` (the matched working cube: 11715×6664
  pixels, the overlap region of both surveys' footprints), merged by
  `rm_synthesis` into a single RM synthesis run across all channels,
  and then cleaned by `rmclean_cubes` — the complete
  preprocessing-through-CLEAN pipeline, run on real multi-survey data,
  completed successfully end to end with no manual intervention.
- Full automated regression suite: 166/166, covering both feature
  areas plus every tool's own parameter handling, I/O parallelism
  consistency, and GPU/CPU numerical agreement.

That same real WALLABY/EMU run is what the README's landing-page demo
GIF is generated from — see [../../README.md](../../README.md#demo).

## Getting started

- [QUICKSTART.md](../../QUICKSTART.md) — build the tools, quick
  reference for every tool.
- [TUTORIAL.md](TUTORIAL.md) — a full walkthrough from a fresh checkout
  to a cleaned RM cube, every command actually run and verified.
- [EXAMPLES.md](EXAMPLES.md) — recipes for real scenarios: single-band
  vs. multi-band, choosing RM-CLEAN stopping criteria, memory/IO
  tuning, GPU vs. CPU.
- [APP_REFERENCE.md](APP_REFERENCE.md) — complete parameter reference
  for all five tools.

## What's next

- Further CLEAN stopping-criteria work: detecting a peak residual that
  diverges or stalls without real progress, a case none of today's
  three stopping criteria catches on its own.
- GPU support for `rmclean_cubes` (currently CPU-only; `rm_synthesis`
  already supports GPU offload).
- Spreading a single large multi-band run across more than one
  machine, so a run that would otherwise take hours on one computer can
  be split up and finish faster — early design work has started, real
  implementation has not.
