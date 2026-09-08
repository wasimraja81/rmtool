# Release Notes

This is the R1.0-beta release of rmtool -- a beta ahead of the first
full release, shared for feedback.

rmtool provides command-line tools for Faraday rotation measure (RM)
synthesis of radio polarization data cubes, and helper scripts to
fetch and process data from CASDA.

## Tools

- `rm_synthesis` -- RM synthesis on Stokes Q/U cubes (single- or
  multi-band), producing polarized-intensity and angle cubes as a
  function of Faraday depth.
- `rmclean_cubes` -- RM-CLEAN deconvolution of the dirty cubes
  `rm_synthesis` produces.
- `reproject_cubes` -- reprojects multiple cubes onto a common sky
  grid.
- `convolve_cubes` -- convolves multiple cubes to a common angular
  resolution.
- `match_cubes` -- combines reprojection and convolution in one step.

## Helpers

- `casda_fetch.py` -- checks CASDA for usable POSSUM Stokes Q/U cubes
  at a given sky position, and fetches them.

## Getting started

See [QUICKSTART.md](QUICKSTART.md), [docs/user/EXAMPLES.md](docs/user/EXAMPLES.md),
and [docs/user/TUTORIAL.md](docs/user/TUTORIAL.md).

_Cut: 2026-09-08T01:38:53Z_
