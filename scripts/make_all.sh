#!/bin/bash 
#
RMTOOL_DIR=/home/wasim/softwares/rmtool
cd ${RMTOOL_DIR}

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
