#!/bin/bash
# Quick build script for rm_synthesis using Makefile
# Usage: ./build.sh [debug|release]

set -e

BUILD_MODE=${1:-release}

echo "================================================"
echo "RM-Synthesis Build Script"
echo "================================================"
echo "Build mode: $BUILD_MODE"
echo "Compiler: gfortran"
echo "Library: CFITSIO"
echo ""

# Check CFITSIO
echo "Checking CFITSIO library..."
if pkg-config --exists cfitsio 2>/dev/null; then
    CFITSIO_VERSION=$(pkg-config --modversion cfitsio)
    echo "✓ CFITSIO found: $CFITSIO_VERSION"
else
    echo "⚠ CFITSIO pkg-config not found, trying system library..."
    if ldconfig -p | grep -q libcfitsio; then
        echo "✓ CFITSIO library found in system"
    else
        echo "✗ CFITSIO not found. Install with:"
        echo "  sudo apt-get install libcfitsio-dev  # Debian/Ubuntu"
        echo "  brew install cfitsio                 # macOS"
        exit 1
    fi
fi

# Check gfortran
echo "Checking gfortran compiler..."
if ! command -v gfortran &> /dev/null; then
    echo "✗ gfortran not found. Install with:"
    echo "  sudo apt-get install gfortran  # Debian/Ubuntu"
    echo "  brew install gcc               # macOS"
    exit 1
fi
GFORTRAN_VERSION=$(gfortran --version | head -1)
echo "✓ $GFORTRAN_VERSION"
echo ""

# Build
echo "Building ($BUILD_MODE mode)..."
make MODE=$BUILD_MODE clean
make MODE=$BUILD_MODE

echo ""
echo "================================================"
echo "✓ Build complete!"
echo "Executable: bin/rm_synthesis"
echo ""
echo "Next steps:"
echo "  ./bin/rm_synthesis --help"
echo "  make install      # Install to /usr/local/bin"
echo "================================================"
