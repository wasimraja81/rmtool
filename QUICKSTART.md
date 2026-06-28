# RM-Synthesis Quick Start Guide

## Contents
1. [Build variants and binary capabilities](#1-build-variants-and-binary-capabilities)
2. [Compile and clean recipes](#2-compile-and-clean-recipes)
3. [Validation test suite](#3-validation-test-suite)
4. [Running on real data](#4-running-on-real-data)
5. [Requirements](#5-requirements)

---

## 1. Build variants and binary capabilities

Three independent binaries can be produced; each lives under `bin/` and is also
symlinked to `bin/rm_synthesis` (last build wins the symlink).

| Make flags | Binary produced | What it can do |
|---|---|---|
| `GPU=0 OMP=0` | `bin/rm_synthesis_release_omp0_gpu0` | **Serial CPU** – one thread, most portable, reference baseline |
| `OMP=1 GPU=0` | `bin/rm_synthesis_release_omp1_gpu0` | **OpenMP CPU** – multi-threaded parallel DFT loop; use this for large images on any Linux box |
| `GPU=1` | `bin/rm_synthesis_release_omp0_gpu1` | **GPU offload** – OpenMP target offload; falls back to host when `use_gpu=n` in cfg, so the same binary is valid for CPU-only dry-runs too |

Key points:
- `MODE=release` (the default) enables `-O3 -march=native`; never compile production runs with `MODE=debug`.
- The GPU binary is built with `-ffast-math -DUSE_GPU`. Setting `use_gpu=n` in the config at runtime makes it behave like the serial CPU binary.
- OMP and GPU are mutually exclusive at compile time (`GPU=1` forces OMP off).
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
two known point sources (RM = −5 and +22 rad/m²), runs all three binaries, checks
that the RM peaks land at the correct positions, and cross-compares numerical outputs.

### 3a. One-shot run (builds + tests everything)

```bash
bash tests/run_tests.sh
```

Expected output (abridged):

```
5. Serial binary – RM peak validation
[OK] src_A: expected RM=-5.0, found RM=-5.00 (err=0.00, tol=1.00)
[OK] src_B: expected RM=+22.0, found RM=+22.00 (err=0.00, tol=1.00)
[PASS] Serial: RM peaks at correct positions

6. OMP binary – bit-identical comparison with serial
[PASS] OMP AMP: matches serial within rtol=1e-4 (FP reassociation)

7. GPU binary – tolerance comparison with serial
[PASS] GPU: RM peaks at correct positions
[PASS] GPU AMP: matches serial within rtol=2e-4

Test Summary: 8 Pass, 0 Fail, 0 Skip
RESULT: ALL PASSED
```

The small OMP/GPU differences (~1e-4 rel.) are normal floating-point
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
python3 tests/compare_cubes.py cube_a.AMP.RMCUBE.FITS cube_b.AMP.RMCUBE.FITS --rtol 2e-4
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

### 4a. GMRT / CASA full-image run

Config: `cfg/rmsynth-casa.fullim.cfg`  
Data: `/home/wasim/softwares/CURR_DEVEL/fitsio_utils/myfitsio.1.0/DATA/`  
RM range: −200 to +200 rad/m², nrm=201, ofac=4  

**CPU (OpenMP, 8 threads):**
```bash
# Build first if not already done
make OMP=1 GPU=0

bash scratch/run_rmsynthesis_test.sh  cfg/rmsynth-casa.fullim.cfg  8  cpu
```

**GPU:**
```bash
# Build GPU binary first
make GPU=1

bash scratch/run_rmsynthesis_test.sh  cfg/rmsynth-casa.fullim.cfg  1  gpu
```

**Outputs written to `scratch/`:**
```
MY_CASA_RMSYNTH_FULLIM_TEST.AMP.RMCUBE.FITS    # |P(RM)|
MY_CASA_RMSYNTH_FULLIM_TEST.PHA.RMCUBE.FITS    # Phase angle (rad)
MY_CASA_RMSYNTH_FULLIM_TEST.NVALID.MAP.FITS    # Valid channel count per pixel
```

**Dry-run** (checks tile memory estimates without touching data):
```bash
# Edit cfg to set dry_run=y, then:
bash scratch/run_rmsynthesis_test.sh  cfg/rmsynth-casa.fullim.cfg  1  cpu
# Reads tile_autotune.cfg and runtime_estimate.txt in scratch/
```

Expected peak: ~5 rad/m² (RL-corrected GMRT 610 MHz data).

---

### 4b. ASKAP / Jennifer full-image run

Config: `cfg/rmsynth-jennifer.fullim.cfg`  
Data: `/data1/tmp/`  
RM range: −500 to +500 rad/m², nrm=101, ofac=1  
Bad channels: `cfg/askap_nan_channels.burdies`  

**CPU (OpenMP, 12 threads):**
```bash
make OMP=1 GPU=0

bash scratch/run_rmsynthesis_test.sh  cfg/rmsynth-jennifer.fullim.cfg  12  cpu
```

**GPU:**
```bash
make GPU=1

bash scratch/run_rmsynthesis_test.sh  cfg/rmsynth-jennifer.fullim.cfg  1  gpu
```

**Outputs written to `scratch/`:**
```
JENNIFER_TOO_FULLIM_TEST.AMP.RMCUBE.FITS
JENNIFER_TOO_FULLIM_TEST.PHA.RMCUBE.FITS
JENNIFER_TOO_FULLIM_TEST.NVALID.MAP.FITS
```

**Tip — memory tuning for large ASKAP cubes:**  
The config has `mem_frac_ram=0.30`, which uses 30% of available **host RAM** per
read block. On machines with ≥64 GB RAM this is usually fine. If you see host
out-of-memory errors, lower it:
```
mem_frac_ram=0.15
```
For **GPU** runs, the device-memory footprint is controlled separately by
`mem_frac_vram` (fraction of VRAM used per offload block) and `gpu_vram_mib`
(VRAM size in MiB; 0 = auto-detect, else override). If you hit a GPU
out-of-memory (`nvptx_alloc error`), lower `mem_frac_vram` (e.g. 0.4) or set
`gpu_vram_mib` to your card's size. `io_overlap=y` opts into overlapped
read/compute (requires a reentrant libcfitsio build; default `n`).
Or run a dry-run first to read the auto-tuned tile hint:
```bash
# Temporarily set dry_run=y in cfg, then:
bash scratch/run_rmsynthesis_test.sh  cfg/rmsynth-jennifer.fullim.cfg  1  cpu
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

| Package | Ubuntu/Debian | macOS (Homebrew) |
|---|---|---|
| Fortran compiler | `gfortran` | `brew install gcc` |
| CFITSIO | `libcfitsio-dev` | `brew install cfitsio` |
| Python 3 + astropy + numpy | `pip install astropy numpy` | `pip install astropy numpy` |
| GPU compiler (optional) | `nvfortran` (NVIDIA HPC SDK) or `gfortran ≥ 14` with libgomp offload | same |

```bash
# Minimal Ubuntu install
sudo apt-get install gfortran libcfitsio-dev make
pip install astropy numpy
```
