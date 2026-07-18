# Changelog

All notable changes to this project are documented in this file.

## [3.0] - 2026-07-18

IO-efficiency milestone: parallel reads, genuine parallel writes, async
tile-write overlap, and the crash/correctness work that came with
building them. All planned tickets (T0-T6) are done and validated
end-to-end on real Setonix production hardware, including T6's actual
write-throughput gain (see Validation below). See
`docs/RELEASE_NOTES_3.0.md` for the full writeup.

### Added
- `io_read_threads` cfg key: N independent read-only CFITSIO handles per
  input cube, reading disjoint channel ranges concurrently.
- `io_write_threads` cfg key: N independent Fortran STREAM I/O units
  write disjoint RM-bin byte ranges of the AMP/PHA output cubes directly,
  bypassing CFITSIO's `ftpsse`/handle machinery for pixel data entirely
  (T6) -- see Fixed/Known limitations below for the two-stage history of
  why this replaced an earlier, unsafe design rather than being the
  first approach shipped.
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
- `tile_read`/`tile_write` log lines now carry a `bytes=<N>` field, and
  the swim-lane plotter renders a new I/O throughput (MB/s) panel from
  it -- stacked directly below the swim-lane/thread panel, sharing its
  time axis so a dip or spike lines up with the Gantt bar above it.
  Absent (no empty panel) for logs predating this field.
- New regression tests: bit-identical `io_overlap=n` vs `y` comparison, a
  structural "no two tile writes ever overlap" invariant check, and a
  bit-identical `io_write_threads=1` vs `=4` comparison across all 8
  output products (`tests/run_tests.sh` §13-14).
- Full, sectioned cfg reference in `README.md`: every key the parser
  accepts, marked required/required-if/optional with its real default,
  cross-checked against `read_cfg_keyval`'s case statements rather than
  written from memory.

### Fixed
- int64-safe flattened tile indices throughout the tile/scatter/mask
  loops (`ipix_tile`, `pix_base`, `src_idx`, `dst_idx`, `ipix_full`,
  `ipix_sub`, `idx_wts`, and the equivalents in `prepare_cpu_data`/
  `prepare_gpu_data`/`tile_extract_gpu_rm_blocked`/
  `cubestat_tail_quantile_maps`) -- same INT32_MAX overflow class as the
  2.0-era `allocate()` fix, previously missed in the runtime index
  arithmetic that actually reads/writes those buffers.
- `io_write_threads>1` root-caused as unsafe (original design): CFITSIO
  aliases repeat read-write opens of an already-open file onto one shared
  internal buffer (`fits_already_open()`), so the "independent" handles
  corrupt each other under concurrent `ftpsse` calls. Produced a real
  SIGSEGV in testing. Hard-clamped to 1 as an interim measure, then fixed
  permanently by replacing the mechanism entirely (see T6 below).
- `io_overlap` write-vs-write race: an initial version only guarded
  buffer-slot reuse (two tiles apart), not whether the *immediately
  preceding* write (a different slot) had finished -- both share the same
  single FITS handle, so two pthreads could call `ftpsse` on it
  concurrently. Produced a real SIGSEGV on a production-scale run (a
  small leftover tile at the bottom of an image whose height wasn't an
  exact multiple of the tile size raced ahead of the previous tile's
  write). Fixed by unconditionally joining any outstanding write before
  dispatching a new one.
- **T6: genuine write-throughput parallelism.** `io_write_threads>1` no
  longer hard-clamped -- each RM-chunk now writes via an independent
  Fortran STREAM I/O unit directly to its byte offset (from `FTGHAD`),
  bypassing CFITSIO for AMP/PHA pixel data entirely, relying only on the
  POSIX guarantee that concurrent writes to disjoint byte ranges of one
  file are safe. Found and fixed a second bug during this work, before it
  ever shipped: leaving CFITSIO's handle for these files open until
  program exit (as the serial path always had) caused CFITSIO's own
  data-fill-check, at final close, to treat the raw-written pixel data as
  "past its own stale end-of-file" bookkeeping (never updated, since the
  raw writer bypassed CFITSIO) and zero-fill over it -- silent data loss,
  not a crash, and only visible in the final on-disk state after close.
  Fixed by closing CFITSIO's handle for AMP/PHA immediately after
  fetching the byte offset it provides, before any raw write happens.
- Swim-lane plotter: the "CPU stage" row (CPU thread-detail view) was
  silently missing its compute segment -- filtered out on the assumption
  that the per-thread lanes above already covered it, which left the row
  summing to less than the tile's actual non-I/O time and a legend entry
  ("CPU compute") that never had a corresponding bar. Restored; the two
  views aren't redundant (the stage row shows *when* the stage ran as a
  whole, the thread lanes show *how* it was parallelised).
- Swim-lane plotter: the GPU pipeline view's synchronous-fallback path
  (a tile that fits in one VRAM sub-block, so there's no async
  double-buffering) expected an old `send N/M` log format the Fortran
  code no longer emits -- current single-shot GPU compute logs plain
  `gpu send`/`gpu recv` notes instead, so non-staged GPU runs were
  silently rendering with no GPU lane at all. Fixed to recognize the
  current format.

### Validation
- Full build matrix remained successful (`OMP/GPU` combinations, zero
  compiler warnings).
- Test suite green (`28/28`), including a full bit-for-bit
  `io_write_threads=1` vs `=4` output comparison and manual verification
  that `io_write_threads>1` combined with `io_overlap=y` also produces
  bit-identical output.
- End-to-end production-scale validation on real Setonix hardware and
  ASKAP/EMU data (13308x11870, 288 channels): the exact case that
  originally crashed now completes without error, `io_read_threads`/
  `io_overlap` confirmed on the real workload, and T6's write-throughput
  gain measured directly -- `io_write_threads=8` dropped write from
  2479.9s (96% of wall time) to 108.3s (6%), a ~23x reduction, taking
  total wall time from 2586.7s to 1945.4s (~25% faster end-to-end). Full
  before/after table in `docs/RELEASE_NOTES_3.0.md`.

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
