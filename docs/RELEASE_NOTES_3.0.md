# Release Notes 3.0 (draft — in progress)

Status: **milestone checkpoint on `optimise-io`, not yet tagged**. This
document describes the IO-efficiency work anticipated to ship as `3.0`.
All planned tickets (T0-T6) are now done; what remains before a formal tag
is production-scale Setonix validation of T6's actual throughput gain (see
"What's next" below) and merging to `develop`.

## Summary

rmtool 3.0 targets IO efficiency for large cubes that don't fit in memory:
parallel channel reads, genuine parallel writes via raw stream I/O,
asynchronous tile-write overlap so a tile's write can proceed on a
background thread concurrently with the next tile's read/compute, and the
crash/correctness work that came with building all three. Three real,
production-scale bugs were found and fixed during this cycle rather than
in the field — all are documented in detail in `docs/ARCHITECTURE.md` and
`planning/IO_PARALLEL_OPTIMISATION_PLAN.md`, since they carry lessons
(about CFITSIO's threading model, its own internal bookkeeping, and what
bit-identical output tests can and can't catch) that matter for anyone
extending this code later.

## Highlights

- **Parallel reads** (`io_read_threads`): N independent read-only CFITSIO
  handles per input cube, each reading a disjoint channel range
  concurrently. Safe by construction — CFITSIO does not alias read-only
  opens the way it does read-write opens (confirmed against the CFITSIO
  source). Default (`io_read_threads=1`) is bit-identical to the prior
  serial path.
- **Parallel writes** (`io_write_threads`): N independent Fortran STREAM
  I/O units write disjoint RM-bin byte ranges of the AMP/PHA output cubes
  directly, bypassing CFITSIO's `ftpsse`/handle machinery for pixel data
  entirely — relying only on the POSIX guarantee that concurrent writes
  to disjoint byte ranges of one file are safe. The original design (N
  independent read-write CFITSIO handles) was found unsafe during this
  cycle — CFITSIO aliases repeat read-write opens of an already-open file
  onto one shared internal buffer, producing a real SIGSEGV during
  testing — and was hard-clamped to 1 as an interim measure. That
  mechanism is now gone entirely: CFITSIO's handle for these two files is
  closed as soon as the byte offset it's needed for is fetched, before
  any raw write happens, so there's no CFITSIO state left for the raw
  writes to conflict with.
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

## Three production bugs, found and fixed during this cycle

1. **`io_write_threads>1` handle aliasing.** CFITSIO's `fits_already_open()`
   explicitly aliases repeat read-write opens of an already-open file onto
   one shared buffer (by CFITSIO's own design — read-only opens are
   exempted for the opposite reason). Concurrent `ftpsse` calls on the
   resulting "independent" handles corrupted that shared buffer, producing
   a SIGSEGV inside `memmove` deep in `libcfitsio.so`. Root-caused against
   the vendored CFITSIO source; fixed by hard-clamping to one handle as an
   interim measure, then permanently by the T6 raw-write mechanism below
   (which never opens more than one CFITSIO handle onto the file at all).
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
3. **Stale CFITSIO buffer flush zeroed T6's raw-written data.** Found
   during T6 development, before it ever shipped. The first raw-write
   implementation kept the CFITSIO handle for the AMP/PHA files open
   until final program cleanup (as the serial path always had); every raw
   write reported success and read back correctly immediately afterward,
   but every pixel was zero once the program exited. Root cause: CFITSIO
   tracks its own "how big is this file" bookkeeping, updated only by its
   own write calls — since the raw writer never went through CFITSIO,
   that bookkeeping stayed pinned at "header only" all run. At final
   close, CFITSIO's data-fill-check saw the real (raw-written) data as
   "past its own stale end of file," treated that as uninitialized space
   needing zero-fill, and overwrote it. Fixed by closing CFITSIO's handle
   for these two files immediately after fetching the byte offset it
   provides, before any raw write happens — retiring CFITSIO's
   bookkeeping before there's anything left for it to disagree with.

Bugs 1 and 2 were timing-dependent races, invisible on small/fast test
data where every stage finishes in milliseconds regardless of ordering —
bit-identical output comparisons alone could not have caught them. Bug 3
was deterministic (reproduced on every run) but only visible in the
*final* on-disk state, after a call ("close the file") that has nothing
to do with pixel data on its face — every write-time signal available
(iostat, an immediate read-back) reported success. The test suite now
includes structural checks alongside output comparisons: an explicit
"no two writes overlap in time" assertion for bug 2, and a full
bit-for-bit `io_write_threads=1` vs `=4` comparison across all 8 output
products (run to completion, including final close) for bug 3.

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
- `io_write_threads>1` now actually parallelises AMP/PHA writes rather
  than being accepted-but-ignored; existing configs that set it will see
  a behaviour change (real speed-up instead of a silent no-op), though
  output remains bit-identical to `io_write_threads=1`.

## Validation snapshot

- Build matrix: CPU serial, CPU OMP, GPU offload, GPU offload + host OMP
  built successfully with zero compiler warnings.
- Tests: full suite passing (`28/28`), including the structural
  concurrency-invariant check for `io_overlap` and the bit-for-bit
  `io_write_threads=1` vs `=4` comparison for T6.
- Production (`io_read_threads`, `io_overlap`): end-to-end validated on
  real Setonix hardware against real ASKAP/EMU data — the exact case that
  originally crashed now completes without error.
- T6 (`io_write_threads>1` raw-write mechanism): validated bit-identical
  on dev-machine test data only so far, including combined with
  `io_overlap=y`. Production-scale Setonix validation is the main
  remaining item before this ships as a tagged `3.0` — see "What's next".

## What's next (before formal `3.0` tag)

- Benchmark `io_write_threads=N` on real Setonix-scale data to measure
  T6's actual write-throughput improvement — the Setonix run that
  motivated building T6 in the first place (write at 96% of wall time,
  see `docs/ARCHITECTURE.md`) is the case this should be measured
  against. Everything above is correctness validation, not yet a
  production timing result.
- Benchmark `io_read_threads=N` against the Lustre stripe count on
  Setonix to tune the recommended default.
- Merge `optimise-io` to `develop` and tag `3.0` once the above lands.
