# Release Notes

This is the R1.0-beta release of rmtool -- a beta ahead of the first
full release, shared for feedback.

rmtool provides command-line tools for Faraday rotation measure (RM)
synthesis of radio polarization data cubes, and helper scripts to
fetch and process data from CASDA.

## Tools

- `rm_synthesis` -- RM synthesis on Stokes Q/U cubes (single- or
  multi-band), producing polarized-intensity and angle cubes as a
  function of Faraday depth. Optionally runs the transform on an NVIDIA
  GPU: by OpenMP target offload (`make GPU=1`, `use_gpu=y`, with
  gfortran and the NVPTX offload compiler), or, as an opt-in build, with
  a cuBLAS matrix product (`make GPU=1 GPU_BACKEND=cuda`, needs the
  cuBLAS libraries). On one 5401 x 5401 x 288 channel cube with 168 RM
  planes, the full run took 132 s with the OpenMP back end and 131 s with
  the CUDA back end on a Tesla P100 (36-thread host), and 484 s and 457 s
  on a GeForce RTX 3050 (16-thread host), against 451 s and 930 s for the
  CPU build on the same two machines (one run each). These runs are
  limited by reading and writing the cubes, not by the GPU; the README
  lists the machines and the stage times. GPU builds now accept
  `rem_mean` above 0, and
  staged GPU runs (cubes larger than the GPU memory budget) work under
  `OMP_TARGET_OFFLOAD=MANDATORY`. See [BUILD.md](BUILD.md) for the GPU
  requirements.
- `rmclean_cubes` -- RM-CLEAN deconvolution of the dirty cubes
  `rm_synthesis` produces.
- `reproject_cubes` -- reprojects multiple cubes onto a common sky
  grid, either an existing reference file or a fresh grid built from
  the inputs themselves (`target_system=icrs|fk5|fk4|galactic`,
  `target_projection=tan|sin|car|ait|hpx`, including a true
  HEALPix-aligned grid via `target_projection=hpx`). 12 interpolation
  kernels (`interp=`), plus `interp=exact`: flux-conserving
  reprojection by spherical overlap area, for when getting the total
  right matters more than resolving structure within a pixel. Refuses
  a reprojection onto cells too coarse to sample an input's own beam,
  and checks free disk space before writing.
- `convolve_cubes` -- convolves multiple cubes to a common angular
  resolution. The FFT plans it measures (also `match_cubes`, when it
  convolves) are saved under `wisdom/` in the current directory, or
  wherever `RMTOOL_WISDOM_DIR` points, so later runs at the same image
  size and thread count start faster and give identical results.
- `match_cubes` -- combines reprojection and convolution in one step,
  chained through memory with no intermediate file; the same
  reference-grid/fresh-grid and `interp=` options as `reproject_cubes`.

## Helpers

- `casda_fetch.py` -- checks CASDA for usable POSSUM Stokes Q/U cubes
  at a given sky position, and fetches them. Accepts a position as
  `--ra`/`--dec` (equatorial frames -- `--frame=icrs|fk5|fk4`) or
  `--lon`/`--lat` (any frame, including `--frame=galactic`), in plain
  decimal degrees, unit-lettered sexagesimal (`21h34m16s`), or
  colon-separated sexagesimal (`21:34:16`). Waits until CASDA has
  finished writing each cutout, then checks every downloaded cube
  against CASDA's MD5 checksum.

## Getting started

`./build.sh` builds all five tools and automatically sets up this
project's own Python environment too (needed for the test suite,
`casda_fetch.py`, and the plotting tools -- not for the compiled
binaries themselves). See [QUICKSTART.md](QUICKSTART.md),
[docs/user/EXAMPLES.md](docs/user/EXAMPLES.md), and
[docs/user/TUTORIAL.md](docs/user/TUTORIAL.md) for more.

_Cut: 2026-11-06T10:58:10Z_
