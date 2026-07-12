# RM-Synthesis Build System

This package supports two modern build approaches for maximum flexibility.

## Quick Start

### Option 1: Simple Makefile (Recommended for Development)

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

The `HOST_OMP` macro gates preprocessor directives in fixed-form Fortran 
(`src/rm_synthesis.f`) and runtime variable checks in free-form Fortran 
(`src/rm_synthesis_mod.f90`), ensuring correct semantic behavior for each variant.

#### Use Cases

- **cpu_serial**: Baseline reference for benchmarking
- **cpu_omp**: Measure CPU parallelism efficiency vs. serial
- **gpu_offload**: GPU kernel speedup with minimal host overhead
- **gpu_offload_hostomp**: Maximum parallelization (GPU kernel + host prep)

### Option 2: CMake (Recommended for Distribution)

```bash
# Build
./cmake_build.sh

# Build debug version
./cmake_build.sh debug

# Install
cd build_cmake
sudo cmake --install .
```

**Advantages:**
- Cross-platform support (Windows, macOS, Linux)
- Automatic dependency detection
- Package configuration files
- Industry standard for C/Fortran projects

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

## Build Targets

| Target | Purpose |
|--------|---------|
| `make` | Build release executable |
| `make MODE=debug` | Build with debugging info |
| `make clean` | Remove build artifacts |
| `make install` | Install to /usr/local/bin |
| `make uninstall` | Remove installation |
| `make help` | Show help message |

## Directory Structure

```
rmtool/
├── src/
│   ├── rm_synthesis_mod.f90      # Modern Fortran module
│   ├── rm_synthesis.f            # Main program
│   ├── myfits_info.f             # FITS utilities
│   └── printerror.f              # Error handling
├── build/                        # Build artifacts (Makefile)
│   ├── modules/                  # Compiled .mod files
│   └── *.o                       # Object files
├── build_cmake/                  # Build artifacts (CMake)
├── bin/                          # Final executable
├── Makefile                      # Simple build
├── CMakeLists.txt                # CMake configuration
├── build.sh                      # Quick build script
├── cmake_build.sh                # CMake build script
└── cfg/                          # Configuration files
```

## Compiler Options

### GFortran (Default)

```bash
# Release (optimized)
make MODE=release

# Debug (with bounds checking)
make MODE=debug
```

### Intel Fortran

Modify `CMakeLists.txt` or `Makefile` to use `ifort`:

```bash
# Makefile
FC=ifort make

# CMake
cmake -DCMAKE_Fortran_COMPILER=ifort ..
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

# Or with aggressive optimization
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
bin/rm_synthesis_release_cpu_serial cfg/rmsynth-casa.fullim.cfg

make OMP=1 GPU=0
bin/rm_synthesis_release_cpu_omp cfg/rmsynth-casa.fullim.cfg

make OMP=0 GPU=1
bin/rm_synthesis_release_gpu_offload cfg/rmsynth-casa.fullim.cfg
```

## Distribution

For distributing the package:

1. Use CMake for configuration
2. Create source tarball: `tar czf rm_synthesis-1.0.tar.gz .`
3. Users build with: `cmake . && make && make install`

## Support

See individual file comments:
- `src/rm_synthesis_mod.f90` - Module documentation
- `src/rm_synthesis.f` - Main program documentation
- `src/myfits_info.f` - FITS interface
