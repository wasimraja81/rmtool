# RM-Synthesis Build System

## Quick Start

Run every command on this page from `$RMTOOL_HOME` -- the directory
you cloned this repository into (see
[QUICKSTART.md](QUICKSTART.md#where-things-go) for the full
convention). Binaries always land in `$RMTOOL_HOME/bin/`.

### Release Tagging Policy

- Public release tags use `R<major>.<minor>` format (for example: `R1.0`,
  `R1.1`, `R2.0`). Full releases are to carry a short codename; `R1.0`,
  not yet released, will be "Confluent Brahmaputra". See
  [RELEASE_NOTES.md](RELEASE_NOTES.md) for the current release.
- `R1.0-beta` is the package's first release: a beta shared for
  feedback ahead of the first full public release, `R1.0` itself.
  Earlier `1.0`–`6.0` version numbers were this project's own internal
  development progression, never tagged as a public release, and are
  not part of the public release line.

### Quick Build Script

```bash
./build.sh
```

Checks dependencies (gfortran, CFITSIO, Starlink AST, FFTW3) with clear
install instructions if anything's missing, then builds everything
[scripts/run_pipeline.sh](scripts/run_pipeline.sh) needs in one go:
`match_cubes`, `reproject_cubes`, `convolve_cubes`, `rmclean_cubes`,
and `rm_synthesis`'s CPU/OpenMP variant. Also attempts `rm_synthesis`'s
GPU variant, best-effort (not required -- reported and skipped if no
working GPU compiler is found). This is the recommended way to build.
For debug builds or any other OMP=/GPU= combination, use `make`
directly (see below).

`./build.sh --cuda` additionally builds the optional CUDA GPU variant
(see "CUDA back end" below). It is skipped with a message if the CUDA
libraries are missing, and it never replaces the default GPU variant.

### Building with Make

For debug builds, a single specific tool, or a specific OMP=/GPU=
combination -- `./build.sh` above already covers the common case.

```bash
# Build (plain `make` with no OMP=/GPU= defaults to the serial CPU variant)
make

# Build with debug symbols
make MODE=debug

# Clean the build artifacts of one flavour (the OMP=/GPU= given)
make clean OMP=1 GPU=0

# Clean every flavour and every tool's build artifacts
make clean-all

# Install rm_synthesis and the cfg/ templates system-wide (optional --
# not needed to run binaries from $RMTOOL_HOME/bin/)
sudo make install

# Show help
make help
```

### Build Variants and Binary Names

Build commands are unchanged (`make OMP=... GPU=...`), but binary names use
clear capability labels:

| Build command | Binary produced | Semantics |
|---|---|---|
| `make OMP=0 GPU=0` | `bin/rm_synthesis_release_cpu_serial` | Pure serial CPU, no parallelism |
| `make OMP=1 GPU=0` | `bin/rm_synthesis_release_cpu_omp` | CPU with OpenMP parallelization |
| `make OMP=0 GPU=1` | `bin/rm_synthesis_release_gpu_offload` | GPU offload, serial host prep |
| `make OMP=1 GPU=1` | `bin/rm_synthesis_release_gpu_offload_hostomp` | GPU offload + host OpenMP |
| `make OMP=0 GPU=1 GPU_BACKEND=cuda` | `bin/rm_synthesis_release_gpu_offload_cuda` | CUDA back end, serial host prep (opt-in) |
| `make OMP=1 GPU=1 GPU_BACKEND=cuda` | `bin/rm_synthesis_release_gpu_offload_hostomp_cuda` | CUDA back end + host OpenMP (opt-in) |

The last two rows are the CUDA flavours described under "CUDA back end".

#### HOST_OMP Policy

Each variant is compiled with a compile-time `HOST_OMP` macro that controls
whether host-side preprocessing loops use OpenMP parallelization:

- **HOST_OMP=0** (`cpu_serial`, `gpu_offload`): Host loops remain serial
- **HOST_OMP=1** (`cpu_omp`, `gpu_offload_hostomp`): Host loops parallelized via `!$omp parallel do`

The `HOST_OMP` macro gates preprocessor directives in the main program 
(`src/rm_synthesis.f90`) and runtime variable checks in free-form Fortran 
(`src/rm_synthesis_mod.f90`), ensuring correct semantic behaviour for each variant.

#### Use Cases

- **cpu_serial**: Baseline reference for benchmarking
- **cpu_omp**: Measure CPU parallelism efficiency vs. serial
- **gpu_offload**: GPU kernel speedup with minimal host overhead
- **gpu_offload_hostomp**: Maximum parallelization (GPU kernel + host prep)
- **gpu_offload_cuda**, **gpu_offload_hostomp_cuda**: the same two, with the
  RM transform done by cuBLAS (opt-in, see below)

#### CUDA back end

`GPU_BACKEND` selects how `GPU=1` builds do the RM transform:

- `GPU_BACKEND=openmp` (default): OpenMP target offload, as above. Needs
  only a GPU-offload-capable compiler.
- `GPU_BACKEND=cuda`: the transform for each block of RM planes is one
  cuBLAS matrix product. Needs `libcublas.so` and `libcudart.so`
  (Debian/Ubuntu: package `nvidia-cuda-dev`, which installs them in
  `/usr/lib/x86_64-linux-gnu`). Set `CUDA_LIBDIR=<directory>` if they are
  elsewhere. `make` stops with a message naming the directory if they are
  missing. It has been built and run with `GPU_FC=gfortran` only.

```bash
make OMP=1 GPU=1 GPU_BACKEND=cuda GPU_FC=gfortran
bin/rm_synthesis_release_gpu_offload_hostomp_cuda cfg/rmsynth.cfg
```

The CUDA flavours are never copied to `bin/rm_synthesis`, and
`make clean ... GPU_BACKEND=cuda` leaves `bin/rm_synthesis` in place. Run the
binary by name. The cfg file is the same as for the OpenMP back end
(`use_gpu=y`). `GPU_BACKEND` has no effect when `GPU=0`.

The startup log of a CUDA binary shows `binary_flavor=gpu_offload_cuda` or
`binary_flavor=gpu_offload_hostomp_cuda`, and prints
"attempting offload with the CUDA back end (cuBLAS)".

## Requirements

### Fortran Compiler
```bash
# Debian/Ubuntu
sudo apt-get install gfortran

# macOS
brew install gcc

# Or any Fortran 2008+ compiler (ifort, pgfortran, etc.)
```

### CFITSIO Library
```bash
# Debian/Ubuntu
sudo apt-get install libcfitsio-dev

# macOS
brew install cfitsio

# From source: https://heasarc.gsfc.nasa.gov/fitsio/
```

### GPU offload (optional, `rm_synthesis` only)
`make GPU=1` with `GPU_FC=gfortran` needs the accelerator compiler for
NVIDIA GPUs, which `gfortran` alone does not include:
```bash
# Debian/Ubuntu
sudo apt-get install gcc-offload-nvptx nvptx-tools
```
A binary built this
way needs, on the machine that runs it, the libgomp NVPTX plugin
(`libgomp-plugin-nvptx1`, installed with `gcc-offload-nvptx`) and the
NVIDIA driver; the CUDA toolkit is not required for the default OpenMP
back end. The opt-in CUDA back end additionally needs the cuBLAS and
CUDA runtime libraries (see "CUDA back end" above). GPU builds are
supported and tested with gfortran and the NVPTX offload compiler only. The
Makefile also contains an untested nvfortran path, which it prefers when
nvfortran is installed, so pass `GPU_FC=gfortran` explicitly.

The container (`docker/dockerfile`) builds the OpenMP GPU offload variants
with gfortran and the NVPTX offload compiler (`gcc-offload-nvptx`). To run
them on a GPU, the host needs an NVIDIA GPU and driver that are made
available to the container. The CUDA back end is not included.

### Starlink AST and FFTW3 (required by `./build.sh`)
The main `rm_synthesis` build links neither: it needs only gfortran + CFITSIO
above (`make OMP=1 GPU=0` builds just that). The standard workflow
(`scripts/casda_fetch.py`, then `scripts/run_pipeline.sh`) also runs the
tools below, so `./build.sh` requires both libraries. Starlink AST is linked
by `reproject_cubes` and `match_cubes`; FFTW3 by `convolve_cubes`,
`match_cubes` and `rmclean_cubes`.
`reproject_cubes`, `convolve_cubes`, `match_cubes` (which consolidates
the other two, with an option to chain them through memory with no
intermediate FITS file), and `rmclean_cubes` (standalone RM-CLEAN,
driving `rmclean_mod` against dirty AMP/PHA cubes `rm_synthesis`
itself wrote) are independent standalone tools (own binaries, own build
targets, not linked into `rm_synthesis`) with their own extra
dependencies:
```bash
# Debian/Ubuntu
sudo apt-get install libstarlink-ast-dev libstarlink-ast-err9 \
    libstarlink-ast-grf3d9 libstarlink-pal-dev libfftw3-dev
```
`rmclean_cubes` needs only FFTW3 + CFITSIO from the above (no Starlink
AST dependency -- it never resamples anything), and `convolve_cubes` needs
FFTW3 but not AST. Both libraries are packaged into
`docker/dockerfile` already, if building via the container is easier
than installing these directly.

### Python 3
Only needed for the test suite, `scripts/casda_fetch.py`, and the
plotting tools -- none of the five compiled binaries above need Python
at all. But `./build.sh` also attempts to set up this project's own
Python venv
(`scripts/setup_venv.sh`, best-effort -- a missing or failed venv never
blocks the binaries above), and that needs `python3` itself already
present, with the standard library's `venv` module available -- on
Debian/Ubuntu, `python3 -m venv` fails without a separate package for
that even when `python3` itself is installed:
```bash
# Debian/Ubuntu
sudo apt-get install python3 python3-venv python3-pip
```
`docker/dockerfile` already includes these.

`scripts/setup_venv.sh` installs `requirements-lock.txt`, not
`requirements.txt` directly -- an exact pinned snapshot of every
Python dependency's version, so the venv, `docker/dockerfile`'s
container, and CI's verify run all install the identical set rather
than each resolving independently against whatever is current on the
day they're built. `requirements.txt` is the unpinned package list to
edit when adding or removing a dependency; regenerate the lock file
afterward from a venv that has it installed:
```bash
~/venv/rmtool/bin/pip install -r requirements.txt
~/venv/rmtool/bin/pip freeze > requirements-lock.txt
```

## Build Targets

| Target | Purpose |
|--------|---------|
| `make` | Build release executable |
| `make MODE=debug` | Build with debugging info |
| `make reproject_cubes` | Build the cross-band sky-grid alignment tool (`bin/reproject_cubes`) |
| `make convolve_cubes` | Build the angular-resolution-matching tool, within a band or across bands (`bin/convolve_cubes`) |
| `make match_cubes` | Build the consolidated reproject+convolve tool, chainable through memory (`bin/match_cubes`) |
| `make rmclean_cubes` | Build the standalone RM-CLEAN tool (`bin/rmclean_cubes`) |
| `make clean` | Remove the build artifacts and `rm_synthesis` binary of the flavour selected by `OMP=`/`GPU=` |
| `make clean-all` | Remove `build/` (every flavour and every tool) and every `rm_synthesis` binary; the other tools rebuild on their next `make` |
| `make install` | Install `rm_synthesis` to `/usr/local/bin` and the `cfg/` templates to `/usr/local/share/rm_synthesis` |
| `make uninstall` | Remove installation |
| `make help` | Show help message |

## Directory Structure

```
rmtool/
├── src/
│   ├── rm_synthesis_mod.f90      # Modern Fortran module
│   ├── rm_synthesis.f90          # Main program; `include`s myfits_info.f90/
│   │                             # printerror.f90 below at compile time
│   ├── myfits_info.f90           # FITS utilities
│   ├── printerror.f90            # Error handling
│   ├── reproject_cubes.f90       # Standalone: cross-band sky-grid alignment
│   ├── gaussft.f90               # gaussft_mod: beam-matching convolution (pure)
│   ├── commonbeam.f90            # commonbeam_mod: smallest common beam
│   ├── convolve_cubes.f90        # Standalone: angular-resolution matching,
│   │                              # within a band or across bands
│   ├── match_cubes.f90           # Standalone: reproject_cubes + convolve_cubes
│   │                              # consolidated, chainable through memory
│   ├── rmclean.f90               # rmclean_mod: RM-CLEAN core (pure computation)
│   └── rmclean_cubes.f90         # Standalone: RM-CLEAN driven against
│                                  # dirty AMP/PHA cubes rm_synthesis wrote
├── build/                        # Build artifacts (Makefile)
│   ├── release_cpu_omp/          # one directory per rm_synthesis flavour
│   ├── release_cpu_serial/       #   (<MODE>_<flavour>), objects and modules
│   ├── reproject_cubes/          # reproject_cubes' own build artifacts
│   ├── convolve_cubes/           # convolve_cubes' own build artifacts
│   ├── match_cubes/              # match_cubes' own build artifacts
│   └── rmclean_cubes/            # rmclean_cubes' own build artifacts
├── bin/                          # Final executables (rm_synthesis, reproject_cubes,
│                                  # convolve_cubes, match_cubes, rmclean_cubes)
├── Makefile                      # Simple build
├── build.sh                      # Quick build script
└── cfg/                          # Configuration files, incl. example_beamLog.txt/.csv
```

## Compiler Options

### GFortran (Default)

```bash
# Release (optimised)
make MODE=release

# Debug (with bounds checking)
make MODE=debug
```

### Other compilers

The compiler is set in the Makefile (`FC := gfortran`), so an
environment variable does not change it; pass it on the `make` command
line instead:

```bash
make FC=ifort OMP=1 GPU=0
```

The Makefile's flags are written for gfortran (for example
`-fallow-argument-mismatch`), so another compiler may reject them. Only
gfortran has been tested.

## Troubleshooting

### CFITSIO Not Found

```bash
# Check if installed
pkg-config --modversion cfitsio

# Or manually specify location
CFITSIO_LIB="-L/usr/local/lib -lcfitsio" make
```

### Compilation Errors

**Error:** `Symbol 'sp' has no type`
- Solution: Ensure `rm_synthesis_mod.f90` compiles first (it defines `sp` kind)

**Error:** `undefined reference to 'ftopen'`
- Solution: Install CFITSIO development package

### AST Link Errors (`reproject_cubes`/`match_cubes` only)

**Error:** `undefined reference to 'astG...'` (e.g. `astGAttr`,
`astGLine`, `astGText`)

`reproject_cubes` and `match_cubes` compile a small hand-written stub,
`src/ast_grf_stub.c`, providing no-op implementations of AST's "GRF"
graphics-callback interface (AST's `Plot` class expects these to exist
at link time, even though this project never creates a `Plot` object,
so they're never called at runtime). `convolve_cubes` and
`rmclean_cubes` don't link AST at all and are unaffected.

The stub's function signatures were transcribed from one system's own
`/usr/include/grf.h`, not from a version-independent contract -- AST's
GRF interface has changed across major versions (v1.0 vs.
v2.0), so a sufficiently different AST installation could need a
different stub.

**What's verified:** build success *and* functional
correctness (`reproject_cubes`/`match_cubes` exercised through
`tests/run_tests.sh`, including resampled values checked against
analytic ground truth, not just a link check) on **Ubuntu 24.04** specifically
-- both this project's own development and its CI/Docker verification
(`docker/dockerfile: FROM ubuntu:24.04`) use that same OS version and
the same package:

```bash
sudo apt-get install libstarlink-ast-dev libstarlink-ast-err9 \
    libstarlink-ast-grf3d9 libstarlink-pal-dev
```

Containers for other base systems are in progress. If you hit this
error on a different OS/AST version in the meantime: compare your
system's own `/usr/include/grf.h` against `src/ast_grf_stub.c` and
adjust the stub's signatures/function set to match, or use the Ubuntu
24.04 Docker image as the known-working path.

### Python Venv Setup Failed

**Error:** `./build.sh` prints `[MISSING] python3 not found -- skipping
venv setup`
- Solution: install `python3`, `python3-venv`, and `python3-pip` (see
  [Requirements](#python-3)
  above), then run `scripts/setup_venv.sh` yourself. The compiled
  binaries are unaffected either way -- this only blocks the test
  suite, `scripts/casda_fetch.py`, and the plotting tools.

**Error:** `scripts/setup_venv.sh` fails partway through with a pip
error, and no network-related message is obviously the cause
- Solution: `python3` was found but the `venv` standard-library module
  wasn't -- on Debian/Ubuntu this needs the separate `python3-venv`
  package even when `python3` itself is already installed. Install it,
  then re-run `scripts/setup_venv.sh`.

## Performance

Release builds (`MODE=release`, the default) already compile with
`-O3 -march=native`, so they are optimised for the machine that builds
them:

```bash
make MODE=release OMP=1 GPU=0
```

Do not set `FFLAGS` to tune this. As an environment variable it is
ignored; on the `make` command line it replaces every flag, including
`-fopenmp` and the preprocessor flags the build needs.

## Timing And Benchmark Logging

All build variants support runtime timing diagnostics through config keys.

Add to your cfg:

```cfg
timing_enabled = y
timing_tile_enabled = y
timing_io_enabled = y
log_output_file =
timing_csv_file = ./timing.csv
```

`log_output_file`:
- empty value writes timing to stdout
- non-empty value appends timing logs to that file

`timing_csv_file`:
- appends one CSV row per run
- useful for scripted benchmark sweeps across `cpu_serial`, `cpu_omp`,
  `gpu_offload`, and `gpu_offload_hostomp`

Example benchmark commands:

```bash
make OMP=0 GPU=0
bin/rm_synthesis_release_cpu_serial cfg/rmsynth.cfg

make OMP=1 GPU=0
bin/rm_synthesis_release_cpu_omp cfg/rmsynth.cfg

make OMP=0 GPU=1
bin/rm_synthesis_release_gpu_offload cfg/rmsynth.cfg
```

## Support

See individual file comments:
- `src/rm_synthesis_mod.f90` - Module documentation
- `src/rm_synthesis.f90` - Main program documentation
- `src/myfits_info.f90` - FITS interface
