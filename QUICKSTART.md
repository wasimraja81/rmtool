#!/bin/bash
# Quick reference card for rm_synthesis build system
# Source this or run: cat QUICKSTART.md

cat << 'EOF'

╔═══════════════════════════════════════════════════════════════════════╗
║                  RM-SYNTHESIS BUILD QUICK START                      ║
╚═══════════════════════════════════════════════════════════════════════╝

📋 DIRECTORY STRUCTURE
─────────────────────────────────────────────────────────────────────
rmtool/
├── src/                      Source code (.f90, .f files)
├── bin/                      Compiled executable
├── build/                    Build artifacts (Makefile mode)
├── build_cmake/              Build artifacts (CMake mode)
├── Makefile                  ← Simple, fast builds
├── CMakeLists.txt            ← Cross-platform, distributable
├── build.sh                  ← Automated script
├── cmake_build.sh            ← CMake automated script
├── BUILD.md                  ← Full documentation
└── cfg/                      Configuration files

🔨 OPTION 1: MAKEFILE (QUICK & SIMPLE)
─────────────────────────────────────────────────────────────────────
Build:
  make                        Build release executable
  make MODE=debug             Build with debugging
  
Clean:
  make clean                  Remove build artifacts
  make install                Install to /usr/local/bin
  make uninstall              Remove installation
  
Help:
  make help                   Show all targets

✅ Status: WORKING ✓

🛠️  OPTION 2: CMAKE (PROFESSIONAL & PORTABLE)
─────────────────────────────────────────────────────────────────────
Automatic build:
  ./cmake_build.sh            Build release
  ./cmake_build.sh debug      Build debug

Manual:
  mkdir build_cmake && cd build_cmake
  cmake -DCMAKE_BUILD_TYPE=Release ..
  cmake --build .
  sudo cmake --install .

✅ Status: READY (create build_cmake/ first)

📝 EXAMPLES
─────────────────────────────────────────────────────────────────────
# Development cycle (fast)
make clean && make MODE=debug
./bin/rm_synthesis <config_file>

# Optimized production build
make MODE=release
make install

# Full clean rebuild
make clean
make MODE=release
./bin/rm_synthesis cfg/your_config.cfg

# Build with custom CFITSIO location
CFITSIO_LIB="-L/opt/cfitsio/lib -lcfitsio" make

✨ SPECIAL FEATURES
─────────────────────────────────────────────────────────────────────
✓ Automatic dependency tracking
✓ Modular compilation (module first, then main)
✓ Two compilation modes (release & debug)
✓ CFITSIO auto-detection
✓ Cross-platform support (via CMake)
✓ Installation to system (/usr/local/bin)
✓ No external dependencies (except CFITSIO)

⚙️  REQUIREMENTS
─────────────────────────────────────────────────────────────────────
1. Fortran compiler (gfortran, ifort, etc.)
2. CFITSIO library (libcfitsio-dev)
3. Make or CMake

Install on Ubuntu/Debian:
  sudo apt-get install gfortran libcfitsio-dev cmake make

Install on macOS:
  brew install gcc cfitsio cmake

🎯 NEXT STEPS
─────────────────────────────────────────────────────────────────────
1. Build: make
2. Test:  ./bin/rm_synthesis cfg/test.cfg
3. For distribution/sharing: use CMake build system

📚 FULL DOCUMENTATION
─────────────────────────────────────────────────────────────────────
See BUILD.md for complete information:
  • Compiler options
  • Troubleshooting
  • Performance tuning
  • Distribution guidelines

EOF
