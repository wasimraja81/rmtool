# rmtool


An HPC package for conducting Faraday Tomography (RM-Synthesis) on radio spectro-polarimetric data. The package is built for all machines - scaling from Low-RAM Desktop PCs to HPC Clusters, with GPU acceleration integration underway.

## Features

- **Config-driven RM synthesis** — Use flexible KEY=VALUE configuration files
- **Subimage support** — Extract and process spatial/spectral subsets via config parameters
- **FITS I/O** — Native support for FITS format using CFITSIO
- **GPU offload toggle** — Enable OpenMP target offload with `use_gpu=y` in config
- **Dual build systems** — Makefile for development, CMake for distribution
- **Fortran 77/90** — High-performance numerical code

## Quick Start

### Prerequisites

- **gfortran** (Fortran compiler)
- **CFITSIO** library (`libcfitsio-dev` on Debian/Ubuntu)
- **GNU Make** or CMake

### Build

```bash
# Simple Makefile build (recommended)
make

# Run executable
./bin/rm_synthesis cfg/your_config.cfg
```

Common explicit build variants:

| Build command | Binary produced |
|---|---|
| `make OMP=0 GPU=0` | `bin/rm_synthesis_release_cpu_serial` |
| `make OMP=1 GPU=0` | `bin/rm_synthesis_release_cpu_omp` |
| `make OMP=0 GPU=1` | `bin/rm_synthesis_release_gpu_offload` |
| `make OMP=1 GPU=1` | `bin/rm_synthesis_release_gpu_offload_hostomp` |

The build commands are unchanged; only the binary naming is now clearer.

See [QUICKSTART.md](QUICKSTART.md) for detailed build instructions.

## Documentation

- **[QUICKSTART.md](QUICKSTART.md)** — Quick reference and build overview
- **[BUILD.md](BUILD.md)** — Comprehensive build system documentation
- **[cfg/CONFIG_README.md](cfg/CONFIG_README.md)** — Configuration file reference
- **[docs/DESIGN_CPU_GPU_TIMELINE_AND_RM_BLOCKING.md](docs/DESIGN_CPU_GPU_TIMELINE_AND_RM_BLOCKING.md)** — Architecture rationale: tiling, RM blocking, CPU/GPU parallelization, offload strategy

## Configuration

RM synthesis is controlled via configuration files in the `cfg/` directory. Every configuration must provide required KEY=VALUE pairs for input/output paths, processing parameters, and RM sampling.

**Example:**

```cfg
# Paths and files
path = /path/to/data/
infileQ = Q_cube.fits
infileU = U_cube.fits
outfile = my_rm_synthesis

# RM sampling (linear grid from beg_rm to end_rm with nrm points)
beg_rm = -500
end_rm = 500
nrm = 101
use_auto_rm_range = 0      # 0=manual range, 1=auto from data
ofac = 4                    # Oversampling factor
fac = 3.14159265358979      # Pi for lambda^2 calculations

# Processing options
remove_badchan = n          # Remove channels with RFI (y/n)
badchan_file = unused.txt   # File listing bad channels
rem_mean = 0                # Remove mean Q/U (0 or 1)
remove_qu_bias = n          # Remove I-based bias from Q/U (y/n)
use_gpu = n                 # GPU offload request (y/n). Alias: use_gpus
output_mode = ap            # Output format: ap (amp+phase) or ri (real+imag)
ap_angle_mode = phase       # Phase mode: phase (arg) or pol (0.5*arg)

# Residuals for bias correction
resiQ = 0.0
slopeQ = 0.0
resiU = 0.0
slopeU = 0.0

# Subimage extraction (optional)
subim = y
subim_ra_blc = 1            # RA first pixel
subim_ra_trc = 256          # RA last pixel (0 = max)
subim_ra_inc = 1            # RA step
subim_dec_blc = 1           # Dec first pixel
subim_dec_trc = 256         # Dec last pixel (0 = max)
subim_dec_inc = 1           # Dec step
subim_chan_blc = 1          # Channel first (0 = first)
subim_chan_trc = 0          # Channel last (0 = all)
subim_chan_inc = 1          # Channel step

# Bias correction inputs (required if remove_qu_bias = y)
path_I = /path/to/data/
infileI = I_cube.fits
```

For complete documentation, see [cfg/CONFIG_README.md](cfg/CONFIG_README.md).

## Timing And Benchmark CSV Output

The runtime logger can emit a human-readable timing summary and an optional CSV
row for automation/benchmark tracking.

Add these keys to your config:

```cfg
# Optional timing controls
log_level = info                  # error|warn|info|debug
timing_enabled = y                # master timing switch
timing_tile_enabled = y           # include tile-level stage timers
timing_io_enabled = y             # include I/O stage timers
log_output_file =                 # empty => stdout, else append to file
timing_csv_file = ./timing.csv    # optional: append one CSV row per run
```

Logging behavior:
- `log_level` controls structured log lines emitted via `log_message`.
	- `error`: errors only
	- `warn`: warnings + errors
	- `info`: run lifecycle messages (recommended default)
	- `debug`: reserved for future verbose diagnostics
- `log_output_file` controls destination for both structured log lines and
	timing summary blocks.
	- empty: output goes to stdout
	- non-empty: output is appended to one consolidated run log file
	- every emitted line in this consolidated log is ISO-8601 timestamped

The consolidated log file includes ISO-8601 local timestamps on structured log
entries, for example:

```text
2026-07-12T14:03:09+10:00 [info] [startup] rm_synthesis run started
2026-07-12T14:03:09+10:00 [info] [startup] binary_flavor=gpu_offload
...
2026-07-12T14:03:11+10:00 [info] [finalize] rm_synthesis run completed
```

When enabled, the run prints:
- `Run summary:` (binary flavor and GPU requested/active state)
- `Timing summary (seconds):` (stage totals and percentages)
- `Macro timing breakdown:` (read I/O, compute RM, compute cubestat, write I/O, overhead)

The optional CSV output writes a header (once) and one row per run with:
- run id and mode
- cube and tile dimensions
- stage timings
- process-level I/O counters

Example run:

```bash
bin/rm_synthesis_release_cpu_serial cfg/rmsynth.cfg
```

Example GPU run:

```bash
bin/rm_synthesis_release_gpu_offload cfg/rmsynth.cfg
```

## GPU Acceleration

rmtool supports OpenMP target offload builds for GPU execution.

Build and run:

```bash
# Build GPU/offload-capable binary
make GPU=1

# Run with a GPU-enabled config (use_gpu=y in cfg/rmsynth.cfg)
bin/rm_synthesis_release_gpu_offload_hostomp cfg/rmsynth.cfg
```

Useful runtime environment variables:

```bash
OMP_TARGET_OFFLOAD=MANDATORY   # fail if offload cannot run
OMP_DEFAULT_DEVICE=0           # choose target device index
```

## Swim-Lane Plot Generation

Use the swim-lane script to visualize overlap between I/O, CPU staging, and GPU
compute from a consolidated run log.

Generate a plot from a run log:

```bash
python scripts/plot_tile_async_swimlane.py \
	--log scratch/RMSYNTH_OUTPUT.run.log \
	--out scratch/tile_async_swimlane.png \
	--run latest \
	--time-axis absolute
```

Key options:

- `--run latest|first|N` selects which detected run block from the log to plot.
- `--time-axis absolute|relative` chooses wall-clock vs seconds-from-run-start.
- `--out` controls output PNG path.

The script also prints summary metrics (interval count, window seconds,
GPU-GPU overlap, CPU-GPU overlap) to stdout.

Design rationale and diagnostic interpretation notes are documented in
[docs/DESIGN_CPU_GPU_TIMELINE_AND_RM_BLOCKING.md](docs/DESIGN_CPU_GPU_TIMELINE_AND_RM_BLOCKING.md).

Example swim-lane plots:

Pipeline/stage-overlap view (long-run async example):

![Swim-lane GPU async example](docs/images/swimlane_gpu_example.png)

Pipeline/stage-overlap view (sync-fallback/smaller-run example):

![Swim-lane pipeline example](docs/images/swimlane_pipeline_example.png)

CPU thread-detail view:

![Swim-lane CPU thread example](docs/images/swimlane_cpu_thread_example.png)

### GPU Validation Scope For Swim-Lane Diagnostics

- Tested on: `NVIDIA GeForce RTX 3050 (6 GiB VRAM)` with GNU OpenMP offload
	(`nvptx`) in this repository's current workflow.
- Not yet validated: AMD ROCm offload targets, Intel GPU offload targets, and
	very old NVIDIA GPUs/toolchains where OpenMP target offload support differs.
- If your platform is unvalidated, treat swim-lane diagnostics as experimental
	and confirm with a short controlled run before production execution.

## GPU Runtime Behavior

- `use_gpu=n` runs host execution.
- `use_gpu=y` requests GPU execution when running a GPU-capable binary (`make GPU=1`).
- If `use_gpu=y` is used with a CPU-only binary, the run prints a warning and falls back to CPU.

Useful runtime env vars for GPU runs:

```bash
OMP_TARGET_OFFLOAD=MANDATORY   # fail if offload cannot run
OMP_DEFAULT_DEVICE=0           # choose target device
```

Build examples:

```bash
# CPU-only OpenMP binary
make GPU=0 OMP=1

# GPU/offload-capable binary
make GPU=1
```

**Output Files:**
- `AP + phase mode`: `OUTBASE.AMP.RMCUBE.FITS` and `OUTBASE.PHA.RMCUBE.FITS`
- `AP + pol mode`: `OUTBASE.AMP.RMCUBE.FITS` and `OUTBASE.POLA.RMCUBE.FITS`
- `RI mode`: `OUTBASE.REAL.RMCUBE.FITS` and `OUTBASE.IMAG.RMCUBE.FITS`
- `Common diagnostics`: `OUTBASE.NVALID.MAP.FITS` (valid-channel count map)
- `When cubestat=y`:
	- `OUTBASE.PEAK.MAP.FITS`
	- `OUTBASE.RM_PEAK.MAP.FITS`
	- `OUTBASE.ANG_PEAK.MAP.FITS`
	- `OUTBASE.SNR.MAP.FITS`

## Project Structure

```
rmtool/
├── src/                  Source code (Fortran 77/90)
│   ├── rm_synthesis.f    Main program
│   ├── rm_synthesis_mod.f90  Config parser module
│   └── legacy/           Legacy tools
├── cfg/                  Configuration files and examples
├── bin/                  Compiled executables
├── build/                Build artifacts (Makefile)
├── Makefile              Build configuration
├── CMakeLists.txt        CMake build configuration
└── BUILD.md              Build documentation
```

## Development

### Branch Structure

- **main** — Stable, production-ready releases
- **develop** — Active development branch

### Release Tags

- Formal release tags use `MAJOR.MINOR` format (for example: `1.0`, `1.1`, `2.0`).
- The first formal release tag is `1.0` on `main`.
- Milestone-style tags can still exist for internal checkpoints, but official releases use the numeric format above.

### Building

```bash
# Makefile (development)
make              # Build release
make MODE=debug   # Build with symbols
make clean        # Clean artifacts

# CMake (distribution)
./cmake_build.sh
cd build_cmake && sudo cmake --install .
```

## License

See [LICENSE](LICENSE) file for details.

## Contact

For questions or contributions, please open an issue or contact the maintainers.
