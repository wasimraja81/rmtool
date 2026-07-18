# Release Notes 3.0 (draft — in progress)

Status: **milestone checkpoint on `develop`, not yet tagged**. This document
describes the IO-efficiency work anticipated to ship as `3.0`. T6 (genuine
write-throughput parallelism) is expected to land before the formal tag.

## Summary

rmtool 3.0 targets IO efficiency for large cubes that don't fit in memory:
parallel channel reads, a parallel-write mechanism (subsequently found unsafe
above 1 thread and hard-capped there), and asynchronous tile-write overlap so
a tile's write can proceed on a background thread concurrently with the next
tile's read/compute. Two real, production-scale bugs were found and fixed
during this cycle rather than in the field — both are documented in detail
in `docs/ARCHITECTURE.md` and `planning/IO_PARALLEL_OPTIMISATION_PLAN.md`,
since they carry lessons (about CFITSIO's threading model and about what
bit-identical output tests can and can't catch) that matter for anyone
extending this code later.

## Highlights

- **Parallel reads** (`io_read_threads`): N independent read-only CFITSIO
  handles per input cube, each reading a disjoint channel range
  concurrently. Safe by construction — CFITSIO does not alias read-only
  opens the way it does read-write opens (confirmed against the CFITSIO
  source). Default (`io_read_threads=1`) is bit-identical to the prior
  serial path.
- **Parallel writes** (`io_write_threads`): implemented, but found to be
  unsafe above 1 during this cycle — CFITSIO aliases repeat read-write
  opens of an already-open file onto one shared internal buffer, so the
  "independent" handles corrupt each other under concurrent writes. This
  produced a real SIGSEGV during testing. **Hard-clamped to 1 at runtime**
  regardless of the cfg value, with a startup warning if a higher value is
  requested. Fixing this for real (T6) needs either genuine multi-process
  writers or bypassing CFITSIO with raw `pwrite()` at pre-computed byte
  offsets.
- **Async tile-write overlap** (`io_overlap`): a tile's write runs on a
  background POSIX thread (`pthread_create`/`pthread_join`), concurrent
  with the next tile's read/mask/prep/compute/cubestat. Deliberately not
  an OpenMP task — a task cannot outlive its parallel region, and keeping
  one open across the next tile's own `!$omp parallel do` calls would
  silently nest and collapse them (this codebase never configures
  `OMP_NESTED`/`omp_set_max_active_levels`). Output tile buffers are
  double-buffered (ping-ponged between two physical slots) so the write
  and the next tile's compute never touch the same memory; default
  (`io_overlap=n`) allocates only one slot and is bit-identical to the
  pre-existing behaviour.
- **int64-safe tile indexing**: the flattened `(x, y, channel/RM)` index
  arithmetic throughout the tile, mask-build, and scatter loops is now
  int64, closing the same INT32_MAX overflow class fixed for `allocate()`
  sizes in 2.0 — this time in the runtime arithmetic that actually
  reads/writes those buffers, which had been missed.
- **Swim-lane plotter improvements**: I/O read/write render as separate
  lanes (previously one shared lane, distinguishable only by colour —
  made concurrent reads/writes hard to see as concurrent); every plot now
  includes a stage-time-totals bar panel; the side-panel thread list was
  trimmed to a count.

## Two production bugs, found and fixed during this cycle

1. **`io_write_threads>1` handle aliasing.** CFITSIO's `fits_already_open()`
   explicitly aliases repeat read-write opens of an already-open file onto
   one shared buffer (by CFITSIO's own design — read-only opens are
   exempted for the opposite reason). Concurrent `ftpsse` calls on the
   resulting "independent" handles corrupted that shared buffer, producing
   a SIGSEGV inside `memmove` deep in `libcfitsio.so`. Root-caused against
   the vendored CFITSIO source; fixed by hard-clamping to one handle.
2. **`io_overlap` write-vs-write race.** The first version only guarded
   buffer-slot reuse (joining a slot's previous write when that same slot
   was needed again, two tiles later) — it never guaranteed the
   *immediately preceding* write, in a different slot, had finished, and
   every tile's write shares the same single FITS handle regardless of
   slot. This crashed on a real 4501×4501 production image: the leftover
   partial tile at the bottom (image height not an exact multiple of the
   tile size) raced ahead of the previous tile's still-running write.
   Fixed by unconditionally joining any outstanding write before
   dispatching a new one — a hard `pthread_join()` barrier, not a
   probabilistic race avoidance, so "at most one write in flight" is now
   a structural guarantee independent of tile timing.

Neither bug was caught by bit-identical output comparisons on small/fast
test data, because both are timing-dependent races invisible when every
stage finishes in milliseconds regardless. The test suite now also
includes structural checks (an explicit "no two writes overlap in time"
assertion, and an `io_write_threads=4` safety-clamp check) that verify the
actual concurrency invariants, not just output correctness.

## Real-world validation

First full production-scale run after the fixes: real ASKAP/EMU Q/U cubes
(13308×11870 pixels, 288 channels) on Setonix, `io_read_threads=8`,
`io_overlap=y`. Completed end-to-end with no error. Per-stage timing from
that run:

| Stage | Total (16 tiles) | % of wall time |
|---|---|---|
| I/O write | 2479.9s | 96% |
| CPU prep | 733.5s | 28% |
| CPU compute | 464.7s | 18% |
| I/O read | 364.1s | 14% |
| CPU cubestat | 28.9s | 1% |
| CPU mask | 12.1s | 0% |

(Percentages sum past 100% because stages genuinely overlap in wall time —
see `docs/ARCHITECTURE.md` for the full breakdown.) Despite write now
dominating wall time, `io_overlap` still delivered a measured **~37% wall-
time reduction** versus a fully-serial equivalent, in the same range
predicted from earlier serial-read profiling — just achieved by hiding
write behind compute instead of behind read, since `io_read_threads` made
read fast enough on its own to no longer be the bottleneck doing the
hiding.

This result is also the concrete argument for prioritizing T6 next: write
is now overwhelmingly the dominant remaining cost, specifically because
it's capped at a single serial handle.

## Compatibility and behaviour notes

- All new cfg keys (`io_read_threads`, `io_write_threads`, `io_overlap`)
  default to values that reproduce the pre-existing serial behaviour
  exactly (`1`, `1`, `n` respectively) — existing configs are unaffected
  until a user opts in.
- RM synthesis numerical output is unchanged; all of this cycle's work is
  IO orchestration, not compute-kernel changes (Non-Negotiable Guardrail
  preserved throughout).
- `io_write_threads>1` remains accepted in a cfg file without erroring —
  it just has no effect beyond a startup warning, so existing configs
  that set it are safe but get no speed-up from it.

## Validation snapshot

- Build matrix: CPU serial, CPU OMP, GPU offload, GPU offload + host OMP
  built successfully with zero compiler warnings.
- Tests: full suite passing (`28/28`), including the new structural
  concurrency-invariant checks.
- Production: end-to-end validated on real Setonix hardware against real
  ASKAP/EMU data — the exact case that originally crashed now completes
  without error.

## What's next (before formal `3.0` tag)

- T6: genuine write-throughput parallelism (multi-process writers, or
  bypassing CFITSIO with raw `pwrite()` at pre-computed byte offsets) —
  see `planning/IO_PARALLEL_OPTIMISATION_PLAN.md`.
- Benchmark `io_read_threads=N` against the Lustre stripe count on
  Setonix to tune the recommended default.
