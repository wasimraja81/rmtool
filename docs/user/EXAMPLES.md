# Examples: Choosing and Combining Parameters for Real Scenarios

[docs/user/TUTORIAL.md](TUTORIAL.md) walks through one linear path end to
end. This document is different: a set of independent recipes for
specific situations — "my bands don't share a sky grid," "I don't know
whether to stop CLEAN on flux or on sigma," "how much memory should I
give this." Jump straight to whichever one matches what you're facing.

Every command below has been run against this repository's own
synthetic test fixtures, and the output shown is what it actually
produces. A few scenarios need data this repo doesn't ship (e.g. two
bands with different angular resolution) — those are noted
where they come up.

## Contents

1. [Single-band quickstart](#1-single-band-quickstart)
2. [Multi-band: which preprocessing do I need?](#2-multi-band-which-preprocessing-do-i-need)
3. [Choosing RM-CLEAN stopping criteria](#3-choosing-rm-clean-stopping-criteria)
4. [Memory tuning (`mem_frac_ram`)](#4-memory-tuning-mem_frac_ram)
5. [I/O parallelism: quick picks](#5-io-parallelism-quick-picks)
6. [GPU vs. CPU: which build should I use?](#6-gpu-vs-cpu-which-build-should-i-use)
7. [Subimage extraction: iterate fast before committing to a full run](#7-subimage-extraction-iterate-fast-before-committing-to-a-full-run)

---

## 1. Single-band quickstart

If you have exactly one Q cube and one U cube, with all channels
already at the same angular resolution, you don't need any of the
preprocessing tools below — just `rm_synthesis` directly. If your
channels' native beam varies with frequency (common even for a single
band), see [§2c](#2c-resolution-mismatched-sky-grid-already-matched--convolve_cubes-only)
— `convolve_cubes` works the same way with one input file as with
several. See
[docs/user/TUTORIAL.md](TUTORIAL.md) for the full walkthrough (build,
generate a test cube, run, inspect the output) and
[docs/user/APP_REFERENCE.md](APP_REFERENCE.md#1-rm_synthesis) for every
config key. The rest of this document assumes you've already done that
once and are now facing a more specific situation.

---

## 2. Multi-band: which preprocessing do I need?

Real multi-band data comes in four combinations: your bands might
already share a sky grid and resolution, or be missing one, the other,
or both. Each needs a different amount of preprocessing before
`rm_synthesis` will accept them together — doing more than necessary
just wastes time, and doing less means `rm_synthesis` will (correctly)
refuse to run rather than silently combine misaligned data. The
resolution half of this ([§2c](#2c-resolution-mismatched-sky-grid-already-matched--convolve_cubes-only))
applies equally to a single band whose own channels don't share one
resolution — sky-grid mismatch is the only cross-band-only
case here.

The four recipes below are covered by regression tests
(`tests/run_tests.sh` sections 38-40): each takes mismatched synthetic
data, runs it through the recipe shown, and confirms `rm_synthesis`
recovers the known injected sources.

### 2a. Bands already matched — no preprocessing needed

If every band is already on the same sky grid and the same angular
resolution, just list them as comma-separated values directly —
`rm_synthesis` merges the channels itself, no separate tool needed:

```cfg
infileQ = bandA.fits,bandB.fits
infileU = bandA_U.fits,bandB_U.fits
```

Every other per-band key (`resiQ`/`slopeQ`/`resiU`/`slopeU`/
`badchan_file`/`subim_chan_blc`/`subim_chan_trc`/`subim_chan_inc`) also
accepts a comma list, one entry per band, in the same order. Run
against this repo's own two-band test fixture (`tests/data/
TEST.Q.FITSCUBE` + `tests/data/TEST_BAND2.Q.FITSCUBE`, which share the
same sky grid by construction):

```bash
bin/rm_synthesis_release_cpu_omp cfg/rmsynth-e2e-multiband-matched-smalltest.cfg
```

completes cleanly, no preprocessing step involved.

### 2b. Sky grid mismatched, resolution already matched — `reproject_cubes` only

If your bands were observed at different pointings (or just have
slightly different pixel grids for any reason) but the same angular
resolution, feeding them to `rm_synthesis` directly gets refused
loudly rather than silently misaligned — this is deliberate, not a bug
to work around:

```
ERROR: RA WCS mismatch for band            2
Reference CRVAL1/CRPIX1/CDELT1 =    180.000000              17  -1.00000005E-03
Band            2  CRVAL1/CRPIX1/CDELT1 =    185.000000              17  -1.00000005E-03
Quitting now...
```

Fix it with `reproject_cubes` first — reproject the mismatched band
(both its Q and U files) onto whichever band you're treating as the
reference:

```bash
make reproject_cubes
bin/reproject_cubes mode=reference reffile=bandA.Q.fits infiles=bandB.Q.fits,bandB.U.fits mem_frac_ram=0.25
```

This writes `bandB.Q_REPROJ.FITS`/`bandB.U_REPROJ.FITS` on the
reference's own grid. Point `rm_synthesis` at those instead of the
originals:

```cfg
infileQ = bandA.Q.fits,bandB.Q_REPROJ.FITS
infileU = bandA.U.fits,bandB.U_REPROJ.FITS
```

`mode=reference` pins the output grid to the reference file's own
extent (what was done above); use `mode=intersection`/`mode=union` if
you'd rather shrink to the overlap or grow to cover every input
instead — see [docs/user/APP_REFERENCE.md](APP_REFERENCE.md#2-reproject_cubes).

### 2c. Resolution mismatched, sky grid already matched — `convolve_cubes` only

If your bands share a sky grid but were taken at different angular
resolutions (a very common real case — different frequencies naturally
have different native beams), convolve every band down to one common
resolution before combining them. The same tool applies just as much
to a *single* band whose own channels don't share one resolution —
this is the normal case for any cube with a real per-channel `BEAMS`
table (`CASAMBM=T`) or ASCII beam log, since the native beam varies
channel-to-channel even within one band; just pass one `infiles=`
entry instead of several:

```bash
make convolve_cubes
bin/convolve_cubes infiles=bandA.Q.fits,bandB.Q.fits mem_frac_ram=0.25
```

With no `target_bmaj`/`target_bmin`/`target_bpa` given, the common
target beam is auto-derived as the smallest beam every good channel of
every input can actually be deconvolved from — you don't need to work
this out by hand. This reads per-channel beams from each input's own
CASA-style `BEAMS` table automatically; if your data doesn't carry one,
pass `beamfiles=path/to/beamlog.txt,path/to/beamlog2.txt` instead (see
`cfg/example_beamLog.txt`/`.csv` for the plain-text format). Output:
`bandA.Q_CONV.FITS`, `bandB.Q_CONV.FITS` (do the same for the U files),
each now sharing one common beam, ready for `rm_synthesis`. Full
parameter list: [docs/user/APP_REFERENCE.md](APP_REFERENCE.md#3-convolve_cubes).

### 2d. Both mismatched — `match_cubes`, chained through memory

If your bands differ in *both* sky grid and resolution, you need both
steps — and `match_cubes` runs them back-to-back without writing an
intermediate FITS file to disk at all, which matters once your cubes
are tens or hundreds of GB each:

```bash
make match_cubes
bin/match_cubes stages=both order=convolve_reproject \
  footprint_mode=reference reffile=bandA.Q.fits \
  infiles=bandB.Q.fits mem_frac_ram=0.25
```

`order=convolve_reproject` (the default) convolves to common
resolution *before* resampling onto the common grid — low-pass
filtering ahead of interpolation avoids baking in aliasing error, and
is usually cheaper too. Run Q and U through the same call together
where possible so both land on the identical output grid by
construction. Full detail:
[docs/user/APP_REFERENCE.md](APP_REFERENCE.md#4-match_cubes).

**Choosing `reffile=` when combining more than two bands:** if your
bands split into more than one distinct sky pointing, pick the
reference from whichever pointing has the most total channels — every
band sharing that SAME pointing gets a fast, near-identity resample,
while every band from a different pointing pays the full
(much slower) resampling cost regardless of which specific file within
it you chose. `match_cubes` checks this for you automatically and
prints an `ADVISORY:` message if your chosen `reffile=` looks like the
less efficient pick — it never changes your choice for you, since
`reffile=` also fixes the output grid.

**Rule of thumb:** if you're not sure which situation you're in, just
try `rm_synthesis` directly first (2a) — it will refuse loudly and
tell you exactly what mismatched (WCS or otherwise), rather than
silently producing a wrong answer. Let that error message tell you
whether you need 2b, 2c, or 2d.

---

## 3. Choosing RM-CLEAN stopping criteria

`rmclean_cubes` gives you three independent, combinable ways to decide
when a pixel is "clean enough." Which one(s) to use depends on what you
actually know about your data:

| Your situation | Use | Why |
|---|---|---|
| You know your noise floor in physical flux units (e.g. from a previous run, or the survey's own sensitivity documentation) | `abs_flux_floor=<value>` alone | A direct, unambiguous "stop below this flux" — no per-pixel estimation involved, so it can't be fooled by an unusual sightline. |
| You don't know the noise floor ahead of time, or it varies significantly across the image | `auto_nsigma=<n>` alone | Estimates each pixel's own noise from its own dirty spectrum — adapts automatically, no prior knowledge needed. |
| Production runs on real data (recommended default) | **Both together** | Whichever fires first wins per pixel. `auto_nsigma` handles the typical case; `abs_flux_floor` is a safety net for the pixels where the per-pixel sigma estimate itself is untrustworthy (it can be biased for individual pixels even though it's accurate on average). |
| You're debugging, or want to see the FULL iteration history regardless of convergence | Neither — `niter` alone | CLEAN runs every pixel for the full `niter` budget; combine with `trace_ix`/`trace_iy`/`log_every` to inspect one pixel's own per-iteration trend. |

**A worked example**, anchored to one dataset's own measured noise
floor, not an invented number:

```cfg
# abs_flux_floor anchored to this dataset's own independently-measured
# ~25.5 uJy/beam noise floor, with a small margin below it
abs_flux_floor = 20uJy
auto_nsigma    = 1.0
niter          = 500
gain           = 0.1
```

**How `abs_flux_floor` values work:** a bare number is native AMP-cube
units (whatever `BUNIT` your dirty cube carries — often Jy/beam); or
give it with a `Jy`/`mJy`/`uJy` suffix, **no space**
(`abs_flux_floor=10mJy`), and it's converted against the cube's own
`BUNIT` automatically.

**On `gain`:** the default `0.1` is a conventional, conservative Hogbom
CLEAN loop gain — larger values (`0.2`-`0.3`) converge in fewer
iterations but risk overshooting for blended/extended structure;
smaller values (`0.05`) are more careful but need a correspondingly
larger `niter`. If a real run's own stop-reason summary shows an
unexpectedly large fraction of pixels hitting the `niter` cap, that's
usually a sign either `niter` is too low or `gain` is too conservative
for the structure in your data — check the summary rmclean_cubes
prints at the end of every run before assuming your data itself is the
problem.

**If pixels are exhausting `niter`** (the run-end summary shows a large
`hit niter cap` percentage): either raise `niter`, or check whether
`abs_flux_floor`/`auto_nsigma` are set sensibly for your data's own
noise level — a floor set too low, or an `auto_nsigma` multiplier set
too small, asks CLEAN to keep going past the point where there's any
real signal left to remove.

---

## 4. Memory tuning (`mem_frac_ram`)

`mem_frac_ram` is a fraction (0, 0.95] of your machine's **total**
system RAM (not whatever happens to be free at that moment — see
[docs/user/APP_REFERENCE.md](APP_REFERENCE.md#shared-parameters) for why
that distinction matters on a shared machine), budgeted for one chunk's
own working set. **All
five tools in this package share the same key, same (0, 0.95] range,
and the same underlying philosophy** — load data in memory-budgeted
chunks so a run scales down to fit a small machine and scales up to use
a big one, never hardcoding "must fit the whole cube in RAM." This
applies just as much to the three preprocessing tools as it does to
`rm_synthesis`/`rmclean_cubes`; it's easy to miss because they don't
show up together in one place elsewhere in the docs — `mem_frac_ram`
(default `0.25`, same validation range) is defined independently in
every one of `rm_synthesis_mod.f90`,
`rmclean_cubes.f90`, `reproject_cubes.f90`, `convolve_cubes.f90`, and
`match_cubes.f90`).

The two tool families chunk along different axes, though, because
their per-pixel work is different in shape:

- **`rm_synthesis`/`rmclean_cubes`** chunk **spatially** — RA/Dec
  tiles, each holding the *full* channel/RM depth for its pixels
  (`tile_ra`/`tile_dec`, auto-tiled full-width strips shrinking further
  only if a single full-width row doesn't fit). Per-pixel cost: input
  Q/U spectra plus whichever outputs you've requested for
  `rm_synthesis` (AMP/PHA or REAL/IMAG, mask, nvalid, cubestat maps);
  2 input + 6 double-buffered output arrays plus the mask cube's own
  channel depth for `rmclean_cubes` — more per-pixel memory
  than `rm_synthesis` for a comparable image size.
- **`reproject_cubes`/`convolve_cubes`/`match_cubes`** chunk along the
  **frequency/plane axis instead** — each block holds the *full*
  spatial (RA/Dec) extent for a handful of channels
  (`block_planes`, auto-computed from `mem_frac_ram` × total RAM ÷
  bytes-per-plane). This makes sense for these tools because
  reprojection/convolution work independently per channel, so there's
  no need to also split spatially. One consequence worth knowing:
  these three tools process blocks strictly one after another and
  parallelise (OpenMP) only *within* a block — so an unnecessarily
  small `mem_frac_ram` here doesn't just mean more CFITSIO calls, it
  also throws away most of the available speedup (each tool prints a
  WARNING if this happens).

**Practical guidance (applies to all five tools):**

- **Default (`0.25`) is a reasonable starting point** on a shared or
  general-purpose machine — leaves headroom for the OS, other
  processes, and CFITSIO/page-cache overhead that isn't counted in the
  budget itself.
- **On a dedicated machine running one job at a time**, raising it
  (`0.5`-`0.7`) lets a smaller image fit in a single tile/block —
  check with `dry_run=y` first on `rm_synthesis` (below) rather than
  guessing; the other four tools print their computed tile/block size
  to the log on every run, so re-run and read it back if unsure.
- **On a memory-constrained machine, or if you see host out-of-memory
  errors**, lower it (`0.15` or below). This doesn't fail the run — it
  just makes tiles/blocks smaller and increases their count, trading a
  bit more per-chunk overhead for a safer memory ceiling. For
  `reproject_cubes`/`convolve_cubes`/`match_cubes` specifically, don't
  go lower than needed — a `block_planes` that's too small also costs
  OpenMP speedup (see above), not just extra I/O calls.
- **`dry_run=y`** exists only on `rm_synthesis`: it reads headers, runs
  the tile planner, writes `tile_autotune.cfg` (a suggested
  `tile_ra`/`tile_dec` you can copy back into your real cfg) and
  `runtime_estimate.txt`, then exits — no pixel data touched, no output
  written. See [docs/user/APP_REFERENCE.md](APP_REFERENCE.md#1-rm_synthesis).
- **GPU runs** (`rm_synthesis`/`rmclean_cubes` only — the three
  preprocessing tools have no GPU path) budget device memory separately
  via `mem_frac_vram`/`gpu_vram_mib` — unrelated to `mem_frac_ram`,
  which still governs the *host*-side tile size on a GPU run.

For `reproject_cubes`/`convolve_cubes`/`match_cubes`, the equivalent
per-tool detail (their own `mem_frac_ram`/`io_overlap` keys and exact
defaults) is in each tool's own section of
[docs/user/APP_REFERENCE.md](APP_REFERENCE.md) — this EXAMPLES.md section is deliberately
just the decision guide, not the full mechanism.

---

## 5. I/O parallelism: quick picks

Three independent keys — `io_read_threads`, `nwriters`,
`io_overlap` — all default to fully serial. Quick picks, not the full
reasoning:

- **Single local disk, small-to-medium cube**: leave everything at
  default. Parallel I/O mostly helps parallel/shared storage.
- **Lustre, multi-server NFS, or cloud block storage**: try
  `io_read_threads`/`nwriters` set to roughly your storage's
  stripe count (often 4-16) — a single I/O stream rarely saturates
  that kind of storage.
- **RAM-constrained** *and* not on parallel storage: leave
  `io_overlap=n` — it doubles the per-tile output buffer RAM for a
  benefit that mostly doesn't exist on a single physical drive anyway.
- **Dedicated machine, parallel storage, RAM to spare**: this is where
  `io_overlap=y` helps — hides tile N's write behind tile
  N+1's read/compute.
- **Always give `OMP_NUM_THREADS` all your cores** — it's the only
  CPU-bound thread count here; the I/O thread keys are mostly
  waiting on disk/network, not competing for compute.

---

## 6. GPU vs. CPU: which build should I use?

- **CPU (OpenMP)** — `make OMP=1 GPU=0`, then any binary without a
  `gpu_offload` tag. Use this unless you specifically have a GPU
  available and a large enough cube that offload is worth the setup.
  This is what [docs/user/TUTORIAL.md](TUTORIAL.md) uses throughout.
- **GPU offload** — `make GPU=1`, `use_gpu=y` in your cfg. On the
  hardware validated to date, GPU throughput is bounded by host-device
  transfer bandwidth more than raw GPU compute — worth it for large
  cubes, less clear-cut for small ones where transfer overhead
  dominates. `rmclean_cubes` has no GPU support at all yet (CPU-only,
  OpenMP one-thread-per-pixel) — this choice only applies to
  `rm_synthesis`.
- **Same cfg file works on both** — leave `use_gpu=n` on a CPU-only
  binary (default behaviour), or leave it `y` anyway: a CPU-only binary
  prints a warning and falls back to CPU rather than failing.
- Full detail: [README.md](../../README.md#gpu-support-work-in-progress)'s
  "GPU Support" section.

---

## 7. Subimage extraction: iterate fast before committing to a full run

Before running RM synthesis on a large cube, it's often worth
testing your config on a small spatial cutout first — a wrong
`beg_rm`/`end_rm`/`nrm` choice, or a config typo, is much cheaper to
discover on a 9×9-pixel subimage than after an hour-long full run.
`subim=y` plus a pixel range does this without any separate cutout
tool:

```cfg
subim         = y
subim_ra_blc  = 8
subim_ra_trc  = 16
subim_dec_blc = 6
subim_dec_trc = 14
```

(1-indexed, inclusive on both ends — the example above extracts a
9×9-pixel region.) Running this against a 32×32 test cube produces
exactly a 9×9×(RM depth) output cube. Channel-axis subsetting works the
same way via
`subim_chan_blc`/`subim_chan_trc`/`subim_chan_inc` — useful for
excluding a known-bad edge of the band without a separate
`badchan_file` entry per channel. Once you're satisfied with the
result on the cutout, set `subim=n` (or delete the `subim_*` keys
entirely — they're ignored when `subim=n`) and run the full cube.

Full key list: [docs/user/APP_REFERENCE.md](APP_REFERENCE.md#1-rm_synthesis).
