# RM-Synthesis Build System

## Quick Start

### Release Tagging Policy

- Official release tags use `MAJOR.MINOR` format (for example: `1.0`, `1.1`, `2.0`, `3.0`).
- Current formal release: `4.1` (on `main`).

### Building with Make

```bash
# Build
make

# Build with debug symbols
make MODE=debug

# Clean build artifacts
make clean

# Install to system
sudo make install

# Show help
make help
```

**Advantages:**
- No dependencies (just Make)
- Fast incremental builds
- Direct control over compilation

### Build Variants and Binary Names

Build commands are unchanged (`make OMP=... GPU=...`), but binary names use
clear capability labels:

| Build command | Binary produced | Semantics |
|---|---|---|
| `make OMP=0 GPU=0` | `bin/rm_synthesis_release_cpu_serial` | Pure serial CPU, no parallelism |
| `make OMP=1 GPU=0` | `bin/rm_synthesis_release_cpu_omp` | CPU with OpenMP parallelization |
| `make OMP=0 GPU=1` | `bin/rm_synthesis_release_gpu_offload` | GPU offload, serial host prep |
| `make OMP=1 GPU=1` | `bin/rm_synthesis_release_gpu_offload_hostomp` | GPU offload + host OpenMP |

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

### Quick Build Script

```bash
# Automated checks + build
./build.sh              # Build release
./build.sh debug        # Build debug
```

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

### Starlink AST and FFTW3 (only needed for `reproject_cubes`/`convolve_cubes`)
The main `rm_synthesis` build needs only gfortran + CFITSIO above.
`reproject_cubes` and `convolve_cubes` are independent standalone tools
(own binaries, own build targets, not linked into `rm_synthesis`) with
their own extra dependencies:
```bash
# Debian/Ubuntu
sudo apt-get install libstarlink-ast-dev libstarlink-ast-err9 \
    libstarlink-ast-grf3d9 libstarlink-pal-dev libfftw3-dev
```
Both packaged into `docker/dockerfile` already, if building via the
container is easier than installing these directly.

## Build Targets

| Target | Purpose |
|--------|---------|
| `make` | Build release executable |
| `make MODE=debug` | Build with debugging info |
| `make reproject_cubes` | Build the cross-band sky-grid alignment tool (`bin/reproject_cubes`) |
| `make convolve_cubes` | Build the cross-band resolution-matching tool (`bin/convolve_cubes`) |
| `make clean` | Remove build artifacts |
| `make install` | Install to /usr/local/bin |
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
│   └── convolve_cubes.f90        # Standalone: cross-band resolution matching
├── build/                        # Build artifacts (Makefile)
│   ├── modules/                  # Compiled .mod files
│   ├── reproject_cubes/          # reproject_cubes' own build artifacts
│   ├── convolve_cubes/           # convolve_cubes' own build artifacts
│   └── *.o                       # Object files
├── bin/                          # Final executables (rm_synthesis, reproject_cubes, convolve_cubes)
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

### Intel Fortran

Modify `Makefile` to use `ifort`:

```bash
FC=ifort make
```

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

## Performance

For production use:

```bash
# Optimize for your CPU
make MODE=release

# Or with aggressive optimisation
FFLAGS="-O3 -march=native -ffast-math" make
```

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

## Distribution

For distributing the package:

1. Create source tarball: `tar czf rm_synthesis-4.1.tar.gz .` (match the current release tag)
2. Users build with: `make && sudo make install`

## Support

See individual file comments:
- `src/rm_synthesis_mod.f90` - Module documentation
- `src/rm_synthesis.f90` - Main program documentation
- `src/myfits_info.f90` - FITS interface
