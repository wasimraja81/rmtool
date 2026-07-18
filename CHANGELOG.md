# Changelog

All notable changes to this project are documented in this file.

## [Unreleased] - working toward 3.0

IO-efficiency milestone: parallel reads, (safely-capped) parallel writes,
async tile-write overlap, and the crash/correctness work that came with
building them. Merged to `develop` as a checkpoint; not yet tagged as a
formal `3.0` release -- T6 (genuine write-throughput parallelism) is the
next major piece expected before that tag. See
`docs/RELEASE_NOTES_3.0.md` for the full writeup.

### Added
- `io_read_threads` cfg key: N independent read-only CFITSIO handles per
  input cube, reading disjoint channel ranges concurrently.
- `io_write_threads` cfg key: implemented, but hard-clamped to 1 at
  runtime with a startup warning if a higher value is requested -- see
  Fixed/Known limitations below.
- `io_overlap` cfg key: tile N's write runs on a background POSIX thread
  concurrent with tile N+1's read/mask/prep/compute/cubestat. Uses a raw
  pthread rather than an OpenMP task specifically so it cannot silently
  nest/collapse the existing OpenMP parallel regions used by
  `io_read_threads` and the compute kernel.
- Auto-tiler RAM planning is aware of `io_overlap`'s doubled output-side
  buffers, so `tile_dec` is planned smaller automatically under the same
  `mem_frac_ram` -- no new user-facing memory configuration.
- Swim-lane plotter: I/O read and I/O write render as separate lanes
  (previously one shared lane distinguished only by colour), and every
  plot now includes a stage-time-totals bar panel (seconds and % of wall
  time per stage, largest first).
- New regression tests: bit-identical `io_overlap=n` vs `y` comparison,
  a structural "no two tile writes ever overlap" invariant check, and an
  `io_write_threads=4` safety-clamp check (`tests/run_tests.sh` §13-14).

### Fixed
- int64-safe flattened tile indices throughout the tile/scatter/mask
  loops (`ipix_tile`, `pix_base`, `src_idx`, `dst_idx`, `ipix_full`,
  `ipix_sub`, `idx_wts`, and the equivalents in `prepare_cpu_data`/
  `prepare_gpu_data`/`tile_extract_gpu_rm_blocked`/
  `cubestat_tail_quantile_maps`) -- same INT32_MAX overflow class as the
  2.0-era `allocate()` fix, previously missed in the runtime index
  arithmetic that actually reads/writes those buffers.
- `io_write_threads>1` root-caused as unsafe: CFITSIO aliases repeat
  read-write opens of an already-open file onto one shared internal
  buffer (`fits_already_open()`), so the "independent" handles corrupt
  each other under concurrent `ftpsse` calls. Produced a real SIGSEGV in
  testing. Now hard-clamped to 1 at runtime regardless of the cfg value.
- `io_overlap` write-vs-write race: an initial version only guarded
  buffer-slot reuse (two tiles apart), not whether the *immediately
  preceding* write (a different slot) had finished -- both share the same
  single FITS handle, so two pthreads could call `ftpsse` on it
  concurrently. Produced a real SIGSEGV on a production-scale run (a
  small leftover tile at the bottom of an image whose height wasn't an
  exact multiple of the tile size raced ahead of the previous tile's
  write). Fixed by unconditionally joining any outstanding write before
  dispatching a new one.

### Known limitations / next steps
- `io_write_threads` cannot currently be used to speed up writes (capped
  at 1); genuine intra-write parallelism needs either multi-process
  writers or bypassing CFITSIO with raw `pwrite()` at pre-computed byte
  offsets -- tracked as T6 in `planning/IO_PARALLEL_OPTIMISATION_PLAN.md`.
- Confirmed on real Setonix production data (13308x11870, 288 channels):
  once `io_read_threads` makes reads fast, output write becomes the
  dominant cost (96% of wall time in one measured run) precisely because
  it's stuck at a single serial handle -- making T6 the highest-value
  remaining work, not a nice-to-have.

### Validation
- Full build matrix remained successful (`OMP/GPU` combinations, zero
  compiler warnings).
- Test suite green (`28/28`).
- End-to-end production-scale validation on real Setonix hardware and
  ASKAP/EMU data: the exact case that originally crashed now completes
  without error.

## [2.0] - 2026-07-17

### Added
- Formalized release-cycle documentation:
  - Added this `CHANGELOG.md`.
  - Added `docs/RELEASE_NOTES_2.0.md`.

### Changed
- Tile planner memory accounting split in `src/rm_synthesis.f90`:
  - Host RAM tile budget now uses `bytes_per_tile_pixel_ram`.
  - GPU VRAM sub-block budget now uses `bytes_per_vram_pixel`.
- Updated docs to reflect planner behaviour and measured outcomes:
  - CPU full-image Jennifer benchmark improved after planner split.
  - GPU path remained correct but showed a slight regression on the tested environment.
- Extended swim-lane interpretation notes for CPU thread-detail view:
  - Clarified that single-RM-chunk runs (`nrm_out <= nrm_block_size`) show only odd/non-hatched `cpu_extract` traces.

### Validation
- Full build matrix remained successful (`OMP/GPU` combinations).
- Test suite remained green (`22/22`) during these updates.

## [1.1] - 2026-07-16

### Added
- Expanded observability and timeline diagnostics:
  - Structured stage/tile logging enhancements.
  - Swim-lane plotting workflow and example artifacts.

### Changed
- Host OpenMP performance improvements:
  - Staged gather/scatter loop parallelization in `src/rm_synthesis.f90`.
  - Pack/copy parallelization in `src/rm_synthesis_mod.f90` with nested-region guard (`omp_in_parallel`).
- Documentation refresh across build, parallelism, and design notes.

### Fixed
- Docker tag-build helper compatibility for shell invocation:
  - `docker/build_push_from_tag.sh` now re-execs under bash when invoked with `sh`.
