# rmtool

RM synthesis tools for astronomical radio spectro-polarimetry data processing. This package implements Rotation Measure synthesis algorithms for analyzing polarized radio emissions.

## Features

- **Config-driven RM synthesis** — Use flexible KEY=VALUE configuration files
- **Subimage support** — Extract and process spatial/spectral subsets via config parameters
- **FITS I/O** — Native support for FITS format using CFITSIO
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

See [QUICKSTART.md](QUICKSTART.md) for detailed build instructions.

## Documentation

- **[QUICKSTART.md](QUICKSTART.md)** — Quick reference and build overview
- **[BUILD.md](BUILD.md)** — Comprehensive build system documentation
- **[cfg/CONFIG_README.md](cfg/CONFIG_README.md)** — Configuration file reference

## Configuration

RM synthesis is controlled via configuration files in the `cfg/` directory. Key parameters include:

```cfg
# Input FITS cubes
q_fits_file = path/to/Q.fits
u_fits_file = path/to/U.fits

# RM synthesis parameters
rm_min = -500
rm_max = 500
rm_step = 5

# Subimage extraction (optional)
subim = y
subim_ra_blc = 1
subim_ra_trc = 100
subim_dec_blc = 1
subim_dec_trc = 100
```

See [cfg/CONFIG_README.md](cfg/CONFIG_README.md) for complete parameter list.

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
