#!/bin/bash
# One-command build for the full rmtool pipeline: every tool
# scripts/run_pipeline.sh needs -- match_cubes, reproject_cubes,
# convolve_cubes, rmclean_cubes, and rm_synthesis's CPU/OpenMP variant
# -- plus a best-effort rm_synthesis GPU build. Checks dependencies
# first, with the same install commands BUILD.md's own Requirements
# section documents, rather than surfacing a raw compiler/linker
# error partway through.
#
# Dependency checks do a real trivial link against each library
# (exactly what the actual build does) rather than searching
# ldconfig's shared-library cache -- that cache misses static
# libraries and any library outside its scanned paths even when the
# compiler/linker can find it fine, which is a false negative, not a
# safe default.
#
# rm_synthesis's GPU variant (OMP=1 GPU=1) is the only piece not
# guaranteed: it needs nvfortran or a GPU-offload-capable gfortran
# (the Makefile auto-selects between them, preferring nvfortran -- see
# `make help`), neither of which is guaranteed present. Its failure is
# reported and skipped, not fatal -- everything else here is expected
# to build on any machine with the dependencies below installed.
#
# For debug builds, every other OMP=/GPU= combination, or anything
# else non-default, use `make` directly -- see BUILD.md.

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${ROOT_DIR}"

log() { echo "[build] $*"; }

echo "================================================"
echo "rmtool build"
echo "================================================"
echo ""

MISSING=0
DEPCHECK_DIR="${ROOT_DIR}/build/.depcheck"
mkdir -p "${DEPCHECK_DIR}"
trap 'rm -rf "${DEPCHECK_DIR}"' EXIT
cat > "${DEPCHECK_DIR}/t.f90" <<'EOF'
program t
end program t
EOF

# check_cmd <human-name> <command> <install-hint-lines...>
check_cmd() {
  local name="$1" cmd="$2"; shift 2
  if command -v "${cmd}" >/dev/null 2>&1; then
    echo "  [OK] ${name} found: $("${cmd}" --version 2>&1 | head -1)"
  else
    echo "  [MISSING] ${name} not found."
    printf '      %s\n' "$@"
    MISSING=1
  fi
}

# check_lib <human-name> <link-flags-as-one-string> <install-hint-lines...>
# Real link test against gfortran -- the same mechanism the actual
# build uses (Makefile's CFITSIO_LIB/AST_LIBS/FFTW_LIBS are plain -l
# flags, no -L path, so "can gfortran link this" IS the real check).
check_lib() {
  local name="$1" linkflags="$2"; shift 2
  if gfortran "${DEPCHECK_DIR}/t.f90" ${linkflags} -o "${DEPCHECK_DIR}/t" >/dev/null 2>&1; then
    echo "  [OK] ${name} found"
  else
    echo "  [MISSING] ${name} -- gfortran could not link against it."
    printf '      %s\n' "$@"
    MISSING=1
  fi
}

log "Checking dependencies (see BUILD.md's Requirements section for full detail)..."
check_cmd "gfortran" gfortran \
  "Debian/Ubuntu: sudo apt-get install gfortran" \
  "macOS:         brew install gcc"
check_lib "CFITSIO" "-lcfitsio" \
  "Debian/Ubuntu: sudo apt-get install libcfitsio-dev" \
  "macOS:         brew install cfitsio"
check_lib "Starlink AST" "-lstarlink_ast -lstarlink_ast_err -lstarlink_ast_grf3d" \
  "Debian/Ubuntu: sudo apt-get install libstarlink-ast-dev libstarlink-ast-err9 libstarlink-ast-grf3d9 libstarlink-pal-dev" \
  "(needed by reproject_cubes/convolve_cubes/match_cubes; rmclean_cubes doesn't need it)"
check_lib "FFTW3" "-lfftw3 -lfftw3f -lfftw3f_omp" \
  "Debian/Ubuntu: sudo apt-get install libfftw3-dev" \
  "(needed by convolve_cubes/match_cubes/rmclean_cubes)"

if [[ "${MISSING}" -eq 1 ]]; then
  echo ""
  echo "Missing dependencies listed above -- install them and re-run, or build via Docker instead (docker/dockerfile bundles everything already)."
  exit 1
fi
echo ""

log "Building rm_synthesis, CPU/OpenMP (OMP=1 GPU=0) -- required, any machine"
make clean OMP=1 GPU=0
make OMP=1 GPU=0

log "Building ancillary tools (match_cubes, reproject_cubes, convolve_cubes, rmclean_cubes)"
rm -rf build/reproject_cubes build/convolve_cubes build/match_cubes build/rmclean_cubes
rm -f bin/reproject_cubes bin/convolve_cubes bin/match_cubes bin/rmclean_cubes
make reproject_cubes
make convolve_cubes
make match_cubes
make rmclean_cubes

log "Attempting rm_synthesis, GPU (OMP=1 GPU=1) -- best effort, not required"
GPU_BUILT=0
if make clean OMP=1 GPU=1 && make OMP=1 GPU=1; then
  GPU_BUILT=1
else
  echo "  GPU variant not built -- no working GPU compiler found (everything above is unaffected). See BUILD.md for GPU setup, or run 'make OMP=1 GPU=1 GPU_FC=<nvfortran|gfortran>' yourself for a detailed error."
fi

echo ""
echo "================================================"
echo "Build complete"
echo "================================================"
echo "  bin/match_cubes"
echo "  bin/reproject_cubes"
echo "  bin/convolve_cubes"
echo "  bin/rmclean_cubes"
echo "  bin/rm_synthesis_release_cpu_omp"
if [[ "${GPU_BUILT}" -eq 1 ]]; then
  echo "  bin/rm_synthesis_release_gpu_offload_hostomp"
fi
echo ""
echo "Next steps:"
echo "  scripts/run_pipeline.sh <pipeline_cfg>   # see cfg/pipeline-example.cfg"
echo "  bin/rm_synthesis_release_cpu_omp --help"
echo "  sudo make install                        # install rm_synthesis to /usr/local/bin"
echo "================================================"
