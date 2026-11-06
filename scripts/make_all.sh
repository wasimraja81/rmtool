#!/bin/bash
# Usage: scripts/make_all.sh [--cuda]
#   --cuda  also build the optional cuBLAS GPU flavour (needs libcublas, libcudart)
RMTOOL_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${RMTOOL_DIR}"

WITH_CUDA=0
case "${1:-}" in
  "") ;;
  --cuda) WITH_CUDA=1 ;;
  *) echo "Usage: scripts/make_all.sh [--cuda]" >&2; exit 1 ;;
esac

# clean build 
make clean OMP=0 GPU=0
make clean OMP=0 GPU=1
make clean OMP=1 GPU=0
make clean OMP=1 GPU=1

# Now build them
make OMP=0 GPU=0
make OMP=0 GPU=1
make OMP=1 GPU=0
make OMP=1 GPU=1

# Opt-in CUDA flavour; never replaces bin/rm_synthesis
if [[ "${WITH_CUDA}" -eq 1 ]]; then
  make clean OMP=1 GPU=1 GPU_BACKEND=cuda GPU_FC=gfortran
  make OMP=1 GPU=1 GPU_BACKEND=cuda GPU_FC=gfortran
fi

# Ancillary standalone tools (reproject_cubes, convolve_cubes,
# match_cubes, rmclean_cubes): each has exactly ONE build flavor --
# always CPU+OpenMP, no GPU offload variant. Their Makefile targets
# hard-code CPU_OPTFLAGS/CPU_OMPFLAGS directly and never read the OMP=/
# GPU= variables at all, so there is no 4-way matrix to rebuild here,
# unlike rm_synthesis above.
rm -rf build/reproject_cubes build/convolve_cubes build/match_cubes build/rmclean_cubes
rm -f bin/reproject_cubes bin/convolve_cubes bin/match_cubes bin/rmclean_cubes
make reproject_cubes
make convolve_cubes
make match_cubes
make rmclean_cubes
