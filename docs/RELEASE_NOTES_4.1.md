# Release Notes 4.1

Status: **tagged `4.1` on `main`, 2026-07-20** (merged from `develop`).

## Summary

4.1 is a diagnostics milestone. It closes a real gap found while looking
into a real large-cube production run: the swim-lane plot's very first
read bar and very first write bar looked oddly slow next to every tile
after them, and it turned out that wasn't disk speed at all -- it was
one-time setup work (reading the input cubes' true dimensions, and
staking out the full size of the output files on disk) that always ran,
but was never timed or logged anywhere, so it silently rode along inside
whatever the first tile's numbers looked like. Alongside that, a second,
long-unused way of building the project was found to have quietly
stopped working, and was removed along with a few leftover duplicate
source files.

## Highlights

### Two new timed stages: `io_read_init` and `io_write_init`

- **What they cover.** Before the tile loop can even begin, the code has
  to open the Q and U input cubes once just to learn their true size
  (`myfits_info`), and -- only when `io_write_threads>1`, where output
  pixel data is written via raw stream I/O rather than through CFITSIO
  -- it closes the AMP/PHA output handles for the first time right after
  fetching their data-start offsets, which is the moment CFITSIO commits
  those files to their full declared size on disk. Both of these were
  always happening. Neither was ever counted in any of the run's named
  timing stages.
- **What changed.** Both now get their own `timer_start`/`timer_stop`
  pair (`STAGE_IO_READ_INIT`, `STAGE_IO_WRITE_INIT`), their own line in
  the run log via the existing `log_tile_bounds` call (same `tile_read`/
  `tile_write` categories and colours the swim-lane plotter already
  recognises, so the plot picks these up as one more read/write bar with
  no code changes on the plotting side), and their own
  `io_read_init_sec`/`io_write_init_sec` columns in the timing CSV.
- **Real byte counts, not placeholders.** The read side reports the true
  number of FITS header bytes actually read (via `FTGHSP`'s existing-key
  count, rounded up to whole 2880-byte header blocks) -- this required
  adding a `header_bytes` output argument to `myfits_info` (now on
  `src/myfits_info.f90`, the file that's actually compiled; see
  "Removed" below). The write side reports the true declared size of the
  AMP and PHA output cubes (`nx_out * ny_out * nrm_out * 4 bytes`, each),
  matching the file-size growth an operator would actually see on disk.
- **The printed timing summary now reconciles.** `stage_total_sec` was
  already an independent wall-clock measurement (from the program's very
  first timer call to its last, not a sum of the named stages), so a run
  could show, say, an unexplained ~350-second gap between the reported
  total and the sum of the 8 previously-named stages. That gap was real
  work, just uninstrumented. The "Macro timing breakdown" printed at the
  end of a run now folds both new stages into the read/write totals, so
  the numbers finally add up to something an operator can check rather
  than take on faith.

### Build cleanup: removed a build path that had quietly broken

While tracing through the code for the above, a second, CMake-based way
of building the project (`CMakeLists.txt`, `cmake_build.sh`) was found to
no longer work at all -- it never learned about the `-cpp` preprocessor
flag the real build (the plain `Makefile`, used by every actual release
and by the Docker image) has depended on for a while, so `#if`/`#ifdef`
conditionals in `rm_synthesis_mod.f90` were being fed to the Fortran
compiler as literal source text instead of being evaluated. Confirmed by
actually running `cmake . && make`, which fails immediately. Both files
are removed, along with three duplicate fixed-form source files
(`src/myfits_info.f`, `src/printerror.f`, `src/rm_synthesis.f`) that only
existed because the project's `.f90` free-form versions had superseded
them without anyone deleting the originals -- only the `.f90` files were
ever actually compiled (pulled into `rm_synthesis.f90` via a plain
Fortran `include`). Every doc that referenced the removed files or the
CMake path (`README.md`, `BUILD.md`, `docs/ARCHITECTURE.md`,
`cfg/CONFIG_README.md`) has been corrected.

## Validation

- Clean 4-variant rebuild (0 errors, 0 new warnings) after every source
  change in this release.
- Full 28/28 test suite, including the two tests that specifically
  re-check bit-identical output under `io_overlap` and
  `io_write_threads>1` -- the exact code path the new write-init timing
  sits next to.
- New log markers hand-verified on a real run with debug logging and
  `io_write_threads=4`: both new intervals appear first, in correct time
  order, ahead of tile 1's own read and write, and their reported byte
  counts check out exactly against the test cube's declared dimensions.
- Confirmed the Docker build recipe (`docker/dockerfile`) never
  referenced CMake or the removed files -- it drives the same `Makefile`
  as every other build path, so this cleanup has no effect on it.

## Compatibility and behaviour notes

- No change to any numerical FITS output -- confirmed via the standard
  bit-identical regression checks.
- Two new CSV columns (`io_read_init_sec`, `io_write_init_sec`) appended
  to the end of the timing CSV's existing column list; no existing
  column was renamed, reordered, or removed.
- No cfg keys were added, removed, or renamed.
- If `use_input_mask=y`, the optional mask-cube dimension read still
  happens exactly as before, but is not yet folded into `io_read_init`'s
  timing or byte count -- a small, known gap, not a correctness issue.

## What shipped in this tag

- `STAGE_IO_READ_INIT` / `STAGE_IO_WRITE_INIT` timing stages, their CSV
  columns, and their swim-lane-compatible log markers.
- A `header_bytes` output argument on `myfits_info` (`src/myfits_info.f90`).
- Removal of the CMake build path and three duplicate `.f` source files,
  and the doc corrections that go with it.

## What's next (beyond this tag)

- Folding the optional input-mask-cube dimension read into `io_read_init`
  if that code path turns out to matter for a real run.
- The original diagnostics question that motivated this release --
  understanding the full shape of a real large-cube run's I/O timing --
  continues on the `diagnostics-heuristics` branch.
