# RM-Synthesis Quick Start Guide

## Contents
1. [Build variants and binary capabilities](#1-build-variants-and-binary-capabilities)
2. [Compile and clean recipes](#2-compile-and-clean-recipes)
3. [Validation test suite](#3-validation-test-suite)
4. [Running on real data](#4-running-on-real-data)
5. [Requirements](#5-requirements)
6. [Swim-lane plots](#6-swim-lane-plots)
7. [Further reading](#7-further-reading)
8. [Preprocessing: reproject_cubes, convolve_cubes, match_cubes](#8-preprocessing-reproject_cubes-convolve_cubes-match_cubes)
9. [RM-CLEAN: rmclean_cubes](#9-rm-clean-rmclean_cubes)

---

## 1. Build variants and binary capabilities

Four independent binaries can be produced; each lives under `bin/` and is also
copied to `bin/rm_synthesis` (last build wins the copy).

| Make flags | Binary produced | What it can do |
|---|---|---|
| `GPU=0 OMP=0` | `bin/rm_synthesis_release_cpu_serial` | **Serial CPU** – one thread, most portable, reference baseline |
| `OMP=1 GPU=0` | `bin/rm_synthesis_release_cpu_omp` | **OpenMP CPU** – multi-threaded parallel DFT loop; use this for large images on any Linux box |
| `GPU=1 OMP=0` | `bin/rm_synthesis_release_gpu_offload` | **GPU offload** – OpenMP target offload build without host-OMP flavor tag |
| `GPU=1 OMP=1` | `bin/rm_synthesis_release_gpu_offload_hostomp` | **GPU offload + host-OMP flavor** |

Key points:
- `MODE=release` (the default) enables `-O3 -march=native`; never compile production runs with `MODE=debug`.
- The GPU binary always defines `-DUSE_GPU`; the `gfortran`-offload path additionally adds `-ffast-math` (the default, auto-selected-first `nvfortran` path does not). Setting `use_gpu=n` in the config at runtime makes it behave like the serial CPU binary.
- GPU and OMP can be enabled together at compile time (`GPU=1 OMP=1`).
- `OMP_NUM_THREADS` controls thread count at runtime for the OMP binary.

---

## 2. Compile and clean recipes

Run all commands from the **repository root**.

### 2a. First-time build (all three variants)

```bash
# Serial CPU  (reference / most portable)
make GPU=0 OMP=0

# OpenMP CPU  (recommended for multi-core production runs)
make OMP=1 GPU=0

# GPU offload  (auto-selects nvfortran → gfortran; can override with GPU_FC=...)
make GPU=1

# --- or with an explicit GPU compiler ---
make GPU=1 GPU_FC=nvfortran     # NVIDIA HPC SDK
make GPU=1 GPU_FC=gfortran      # GNU OpenMP offload (libgomp)
```

All three can coexist — each writes to its own `build/release_omp<N>_gpu<N>/` tree.

### 2b. Partial clean (one variant only)

```bash
make clean GPU=0 OMP=0   # Remove serial artifacts only
make clean OMP=1 GPU=0   # Remove OMP artifacts only
make clean GPU=1         # Remove GPU artifacts only
```

### 2c. Full clean (everything)

```bash
make clean-all
```

### 2d. Debug build (for development only — slow at runtime)

```bash
make MODE=debug GPU=0 OMP=0    # Serial debug
make MODE=debug OMP=1 GPU=0    # OMP debug
```

### 2e. Custom CFITSIO location

```bash
CFITSIO_LIB="-L/opt/cfitsio/lib -lcfitsio" make OMP=1 GPU=0
```

### 2f. Install / uninstall system-wide

```bash
make install           # Copies bin/rm_synthesis → /usr/local/bin/
make uninstall
```

---

## 3. Validation test suite

The test suite lives in `tests/`. It generates synthetic Q/U FITS cubes containing
two known point sources (RM = -5 and +22 rad/m^2), then runs a multi-stage
validation workflow that currently includes:

- CPU serial and CPU OpenMP builds/runs
- GPU and GPU+HostOMP builds/runs when available
- RM peak checks and cube comparisons vs serial reference
- GPU staging-path activation and staged-vs-nonstaged comparison
- auto-tiling shape check (full-RA Dec strips)
- bad-channel masking checks (serial/OMP/GPU)
- cubestat output map checks
- timing summary/CSV emission checks
- `io_overlap` bit-identical comparison plus a structural "no two tile
  writes ever overlap" invariant check
- `nwriters>1` bit-identical comparison across all 8 output products

### 3a. One-shot run (builds + tests everything)

```bash
bash tests/run_tests.sh
```

Representative output (abridged):

```
5. Serial binary – RM peak validation
[OK] src_A: expected RM=-5.0, found RM=-5.00 (err=0.00, tol=1.00)
[OK] src_B: expected RM=+22.0, found RM=+22.00 (err=0.00, tol=1.00)
[PASS] Serial: RM peaks at correct positions

6. OMP binary – bit-identical comparison with serial
[PASS] OMP AMP: matches serial within rtol=1e-4 (FP reassociation)

7. GPU binary – tolerance comparison with serial
[PASS] GPU: RM peaks at correct positions
[PASS] GPU AMP: matches serial within rtol=2e-3

12. Timing report + CSV validation
[PASS] Timing markers present (serial)

Test Summary
Total : ...
Pass  : ...
Fail  : ...
Skip  : ...
RESULT: ALL PASSED
```

`Pass/Fail/Skip` counts vary by platform. In particular, GPU sections are
explicitly skipped when a GPU-capable binary cannot be built.

The small OMP/GPU differences (~1e-4 to ~1e-3 rel.) are normal floating-point
reassociation from parallel reductions and `-ffast-math`; the RM peaks themselves
are exact.

### 3b. Run only the peak-check step on an existing cube

```bash
python3 tests/check_rm_peak.py  <path/to/output.AMP.RMCUBE.FITS> \
                                 tests/data/truth.json
```

### 3c. Compare two output cubes directly

```bash
# Bit-exact comparison (serial vs serial re-run)
python3 tests/compare_cubes.py cube_a.AMP.RMCUBE.FITS cube_b.AMP.RMCUBE.FITS --exact

# Relative-tolerance comparison (e.g. GPU vs serial)
python3 tests/compare_cubes.py cube_a.AMP.RMCUBE.FITS cube_b.AMP.RMCUBE.FITS --rtol 2e-3
```

---

## 4. Running on real data

The launcher script `scratch/run_rmsynthesis_test.sh` handles binary selection,
OMP environment, timing, and output-file checking automatically.

**Usage:**
```
bash scratch/run_rmsynthesis_test.sh  <config>  [num_threads]  [backend]

  <config>        Config file name (relative to cfg/) or absolute path
  [num_threads]   OMP thread count for CPU backend  (default: 6)
  [backend]       auto | cpu | gpu                  (default: auto)
                  auto: reads use_gpu= from cfg to pick the binary
```

---

### 4a. Getting a test object

No Q/U data of your own yet? `scripts/casda_fetch.py` checks the
CSIRO ASKAP Science Data Archive (CASDA) for calibrated POSSUM Stokes
Q/U cubes around a sky position, fetches a small cutout, and generates
a ready-to-run pipeline config from whatever it fetched — enough to
try single/multi-band tomography without downloading a full survey
field or hand-assembling a cfg yourself.

Needs a free CASDA/OPAL account (self-register at `opal.atnf.csiro.au`)
and `astroquery` (`~/venv/rmtool/bin/pip install astroquery`, or your
own venv's pip). The first run prompts for your password once and
remembers it via your OS keyring; every later run is silent.

**Step 1 — check what's there** (the default mode, `--run-mode=dry`;
downloads nothing, safe to run any time):
```bash
python3 scripts/casda_fetch.py --target dancingghosts --username you@example.org
```
`--target` takes either one of this project's own presets (`msh15-56`,
`kes27`, `cena`, `dancingghosts` — objects already investigated for
this project, see `docs/dev/CASDA_FETCH_PLAN.md` §4 for how each was
found) or, failing that, any other name resolvable the way SIMBAD/NED
resolve it (e.g. `--target "Kes 27"`); use `--ra 323.56667 --dec
-53.61528 --radius 8` instead for a position not resolvable by name at
all. This step prints how many observations exist, at how many
frequency bands, with a complete Stokes Q/U pair present, and states
directly whether multi-band synthesis is possible at that position or
only single-band. A band only counts as usable when both Q and U are
present as `.conv` (beam-wise matched across ASKAP's 36 PAF beams) — a
band missing that is reported but excluded, since combining
beam-mismatched data would silently corrupt the result.

**Step 2 — fetch:**
```bash
python3 scripts/casda_fetch.py --target dancingghosts --username you@example.org \
    --run-mode=auto --outdir ./dancingghosts_data
```
`--run-mode=auto` fetches every usable band found, the latest-observed
SBID within each, and the latest reprocessed version of each file —
for when you just want the most complete result CASDA currently
supports, with no manual choices. (`--run-mode=select
--select="1,2:51797"` fetches exactly the band/SBID/version
combination you name instead, using band numbers from step 1's own
report — see `--help` for the full syntax.) Each usable band's I/Q/U/V
cubes (I to visually confirm signal, V as a noise-floor check) land in
their own `<outdir>/SB<n>/` subdirectory; a band that doesn't qualify
is reported `SKIPPED`, not silently dropped. If a fetched cube doesn't
already carry a per-channel BEAMS table of its own, fetching also
pulls and curates the observatory's own beamlog for it automatically —
no manual beam bookkeeping needed either way.

**Step 3 — run it:** step 2 also generates a ready-to-run pipeline cfg
chaining `match_cubes` into `rm_synthesis` (`beamfiles=` already filled
in either way, per the previous paragraph), printed at the end of its
own output along with the exact command to run it:
```bash
./scripts/run_pipeline.sh ./dancingghosts_data/dancingghosts_pipeline.cfg
```
See [TUTORIAL.md §6](docs/user/TUTORIAL.md#6-running-the-full-pipeline-in-one-command)
for what that script does. A starter `rmclean_cubes` cfg is generated
alongside it too, but left out of the pipeline cfg's own `stages=` by
default — CLEAN's stopping criteria are a choice to make after looking
at the dirty cube, not a safe default (pass `--chain-rmclean` at fetch
time in step 2 to opt in anyway; it prints an explicit caveat about
what wasn't checked before doing so). See
[EXAMPLES.md §3](docs/user/EXAMPLES.md#3-choosing-rm-clean-stopping-criteria)
for how to choose those criteria once you're ready.

**Re-running against the same `--outdir`:** `--if-exists` controls
what happens when a previous run's output is still there — `check`
(default) reports what's present and refuses to guess; `reuse` uses it
as-is without re-downloading; `refetch` overwrites it; `clean` wipes
that run's whole output (fetched cubes, generated cfgs, and anything
`run_pipeline.sh` produced from them) before fetching fresh. Full
detail on every flag and a longer worked explanation of the archive
concepts above (SBIDs, reprocessing versions, `.conv` beam-matching):
`python3 scripts/casda_fetch.py --help`.

---

### 4b. Generic full-image run

Canonical config: `cfg/rmsynth.cfg`  
Edit the following keys first for your local data:
- `path`
- `infileQ`
- `infileU`
- `outfile`

**CPU (OpenMP, 8 threads):**
```bash
# Build first if not already done
make OMP=1 GPU=0

bash scratch/run_rmsynthesis_test.sh  cfg/rmsynth.cfg  8  cpu
```

**GPU:**
```bash
# Build GPU binary first
make GPU=1

bash scratch/run_rmsynthesis_test.sh  cfg/rmsynth.cfg  1  gpu
```

**Outputs written to `scratch/`:**
```
RMSYNTH_OUTPUT.AMP.RMCUBE.FITS    # |P(RM)|
RMSYNTH_OUTPUT.PHA.RMCUBE.FITS    # Phase angle (rad)
RMSYNTH_OUTPUT.NVALID.MAP.FITS    # Valid channel count per pixel
```

**Dry-run** (checks tile memory estimates without touching data):
```bash
# Edit cfg to set dry_run=y, then:
bash scratch/run_rmsynthesis_test.sh  cfg/rmsynth.cfg  1  cpu
# Reads tile_autotune.cfg and runtime_estimate.txt in scratch/
```

**Tip - memory tuning for large cubes:**  
The config has `mem_frac_ram=0.30`, which uses 30% of available **host RAM** per
read block. If you see host out-of-memory errors, lower it:
```
mem_frac_ram=0.15
```
For **GPU** runs, the device-memory footprint is controlled separately by
`mem_frac_vram` (fraction of VRAM used per offload block) and `gpu_vram_mib`
(VRAM size in MiB; 0 = auto-detect, else override). If you hit a GPU
out-of-memory (`nvptx_alloc error`), lower `mem_frac_vram` (e.g. 0.4) or set
`gpu_vram_mib` to your card's size. `io_overlap=y` opts into running tile
N's write on a background thread concurrent with tile N+1's read/compute
(CFITSIO is already built reentrant -- nothing extra to install; default
`n`). Doubles the per-tile output buffer RAM -- helps most on
dedicated machines with parallel/fast storage and RAM to spare, less
on a single disk or a RAM-constrained box (see
[docs/user/EXAMPLES.md](docs/user/EXAMPLES.md#5-io-parallelism-quick-picks)
for the decision guide). `io_read_threads`/`nwriters` parallelise the
read/write themselves (independent of `io_overlap`).
Or run a dry-run first to read the auto-tuned tile hint:
```bash
# Temporarily set dry_run=y in cfg, then:
bash scratch/run_rmsynthesis_test.sh  cfg/rmsynth.cfg  1  cpu
cat scratch/tile_autotune.cfg          # copy tile_ra / tile_dec back into cfg
cat scratch/runtime_estimate.txt       # wall-time estimate
```

---

### 4c. Environment variables (advanced)

| Variable | Default (run script) | Effect |
|---|---|---|
| `OMP_NUM_THREADS` | arg 2 (default 6) | CPU thread count |
| `OMP_PROC_BIND` | `close` | Thread affinity (CPU mode) |
| `OMP_PLACES` | `cores` | Thread placement (CPU mode) |
| `OMP_TARGET_OFFLOAD` | `MANDATORY` | GPU mode: `MANDATORY` aborts if no GPU; set `DISABLED` for host-fallback testing |
| `OMP_DEFAULT_DEVICE` | `0` | GPU device index (multi-GPU systems) |

---

## 5. Requirements

| Package | Notes |
|---|---|
| Fortran compiler, CFITSIO, GPU compiler, Starlink AST + FFTW3 | See [BUILD.md](BUILD.md#requirements) for exact install commands per platform and which tool needs which dependency (`rmclean_cubes` needs only FFTW3 of that last group, not Starlink AST). |
| Python 3 + astropy + numpy | `pip install astropy numpy` |
| matplotlib + ffmpeg (diagnostic/plotting scripts only, not needed to build or run any Fortran tool) | `pip install matplotlib`, plus `ffmpeg` via your package manager -- used by `scripts/plot_tile_async_swimlane.py`, `scripts/plot_rmclean_advisory.py`, `scripts/animate_fits_cube.py`, `scripts/benchmark_omp.py` |

```bash
# Minimal Ubuntu install (core build only -- see BUILD.md for the rest)
sudo apt-get install gfortran libcfitsio-dev make
pip install astropy numpy
```

---

## 6. Swim-lane plots

Use `scripts/plot_tile_async_swimlane.py` to visualize overlap across I/O, CPU,
and GPU lanes from the run log.

---

## 7. Further reading

Internal architecture and design documentation (memory/tiling
strategy, RM chunking, CPU/GPU timeline diagnostics, I/O parallelism
postmortems) lives under `docs/dev/`.

---

## 8. Preprocessing: reproject_cubes, convolve_cubes, match_cubes

Real multi-band data rarely arrives on the same sky grid or at the same
angular resolution — and even a single band's own channels often don't
share one resolution, since native beam size varies with frequency.
Three standalone tools (own binaries, independent of the main
`rm_synthesis` build) close that gap before a run — either as two
separate calls:

```bash
# 1. Match angular resolution (all input files together, one call --
#    works the same with one file for a single band or several for multi-band)
make convolve_cubes
bin/convolve_cubes infiles=a.fits,b.fits mem_frac_ram=0.25

# 2. Align sky grids
make reproject_cubes
bin/reproject_cubes mode=intersection reffile=ref.fits infiles=a_CONV.fits,b_CONV.fits

# 3. Run multi-band RM synthesis on the now-matched inputs
bin/rm_synthesis your_multiband.cfg
```

or as one `match_cubes` call that chains both stages through memory, with
no intermediate FITS file written — a real saving for 200GB+ cubes:

```bash
make match_cubes
bin/match_cubes stages=both order=convolve_reproject \
  footprint_mode=intersection reffile=ref.fits infiles=a.fits,b.fits \
  mem_frac_ram=0.25
```

`--help` on any of the three tools gives its full option list.
`convolve_cubes`'/`match_cubes`' per-channel beam input (a CASA-style
`BEAMS` table, or a portable ASCII/CSV beam log — see `cfg/
example_beamLog.txt`/`.csv`) and target-beam derivation are covered in
[docs/user/APP_REFERENCE.md](docs/user/APP_REFERENCE.md). All three
tools propagate `CASAMBM`/`BEAMS` beam metadata through to their
outputs, the same way `rm_synthesis` does for its own.

```bash
python scripts/plot_tile_async_swimlane.py \
  --log scratch/RMSYNTH_OUTPUT.run.log \
  --out scratch/tile_async_swimlane.png \
  --run latest \
  --time-axis absolute
```

Useful flags:

- `--run latest|first|N`: choose which run block from the log.
- `--time-axis absolute|relative`: absolute clock vs relative seconds.

Validation scope:

- Tested in this repo workflow on `NVIDIA GeForce RTX 3050 (6 GiB)`.
- Not yet validated on AMD/ROCm, Intel GPU offload paths, or older NVIDIA
  offload stacks with differing OpenMP target support.

---

## 9. RM-CLEAN: rmclean_cubes

`rmclean_cubes` runs standalone RM-CLEAN (Högbom-style deconvolution)
against a dirty AMP/PHA cube pair `rm_synthesis` already wrote,
producing CLEAN component, residual, and restored cubes.

```bash
make rmclean_cubes

bin/rmclean_cubes ampfile=out.AMP.RMCUBE.FITS phafile=out.PHA.RMCUBE.FITS \
  maskfile=out.MASK.CUBE.FITS outfile=out_cleaned \
  abs_flux_floor=20uJy auto_nsigma=1.0 niter=500 gain=0.1
```

CLEAN stops on whichever of three independent criteria fires first:
`niter` (hard iteration cap, always on), `abs_flux_floor=<v>` (stop at
a literal flux value — accepts a `Jy`/`mJy`/`uJy` suffix, e.g. `10mJy`),
and `auto_nsigma=<n>` (stop at `n`× each pixel's own noise sigma,
estimated once per pixel from its own dirty spectrum). Both are
opt-in — give either, both, or neither.

**Outputs written to `<outfile>`:**
```
<outfile>.CLEAN.AMP.RMCUBE.FITS       # CLEAN component amplitude
<outfile>.CLEAN.PHA.RMCUBE.FITS       # CLEAN component phase
<outfile>.RESID.AMP.RMCUBE.FITS       # Residual amplitude
<outfile>.RESID.PHA.RMCUBE.FITS       # Residual phase
<outfile>.RESTORED.AMP.RMCUBE.FITS    # Restored (components + residual) amplitude
<outfile>.RESTORED.PHA.RMCUBE.FITS    # Restored phase
```

Memory/tiling and I/O parallelism (`mem_frac_ram`/`tile_ra`/`tile_dec`/
`io_read_threads`/`nwriters`/`io_overlap`) work exactly like
`rm_synthesis`'s own keys of the same name (see section 4 above) — a
cube far larger than available RAM CLEANs in bounded memory on any
machine, from a laptop to an HPC node.

**If you're migrating an older cfg**: `threshold=`/`threshold_snr=`/
`noise_percentile=`/`noise_nlos=`/`noise_seed=` have been removed
entirely — replace with `abs_flux_floor=`/`auto_nsigma=` above.

Full parameter reference (every key, default, and meaning) in
[docs/user/APP_REFERENCE.md](docs/user/APP_REFERENCE.md#5-rmclean_cubes); a
fully annotated example config in `cfg/rmclean-example.cfg`; a
step-by-step walkthrough tying `rm_synthesis` and `rmclean_cubes`
together in [docs/user/TUTORIAL.md](docs/user/TUTORIAL.md); a decision guide for
choosing between `abs_flux_floor`/`auto_nsigma` in
[docs/user/EXAMPLES.md](docs/user/EXAMPLES.md#3-choosing-rm-clean-stopping-criteria).
