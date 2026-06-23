# rmtool

RM (Rotation Measure) synthesis tools for analyzing polarized radio observations. This package implements Rotation Measure synthesis algorithms to decompose linear polarization into RM spectra for radio spectro-polarimetry data.

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

**Output Files:**
- `AP + phase mode`: `OUTBASE.AMP.RMCUBE.FITS` and `OUTBASE.PHA.RMCUBE.FITS`
- `AP + pol mode`: `OUTBASE.AMP.RMCUBE.FITS` and `OUTBASE.POLA.RMCUBE.FITS`
- `RI mode`: `OUTBASE.REAL.RMCUBE.FITS` and `OUTBASE.IMAG.RMCUBE.FITS`

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
