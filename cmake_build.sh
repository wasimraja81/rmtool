#!/bin/bash
# CMake build script for rm_synthesis
# Usage: ./cmake_build.sh [debug|release]

set -e

BUILD_MODE=${1:-Release}

echo "================================================"
echo "RM-Synthesis CMake Build"
echo "================================================"
echo "Build mode: $BUILD_MODE"
echo ""

# Check CMake
if ! command -v cmake &> /dev/null; then
    echo "✗ CMake not found. Install with:"
    echo "  sudo apt-get install cmake  # Debian/Ubuntu"
    echo "  brew install cmake          # macOS"
    exit 1
fi

echo "✓ CMake version: $(cmake --version | head -1)"
echo "✓ Compiler: $(gfortran --version | head -1)"
echo ""

# Create build directory
BUILD_DIR="build_cmake"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Run CMake
echo "Configuring CMake..."
cmake -DCMAKE_BUILD_TYPE=$BUILD_MODE -DCMAKE_INSTALL_PREFIX=/usr/local ..

# Build
echo ""
echo "Building..."
cmake --build . --config $BUILD_MODE

echo ""
echo "================================================"
echo "✓ CMake build complete!"
echo "Executable: build_cmake/rm_synthesis"
echo ""
echo "Next steps:"
echo "  cd build_cmake"
echo "  sudo cmake --install .  # Install to system"
echo "================================================"
