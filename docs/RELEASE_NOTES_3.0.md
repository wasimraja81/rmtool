# Release Notes 3.0

Status: **tagged `3.0` on `main`, 2026-07-18** (merged from `develop`). All
planned IO-efficiency tickets (T0-T6) are done and validated end-to-end on
real Setonix production hardware — see "Real-world validation" below for
T6's measured write-throughput gain, the result that was the last open
item before this tag.

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
  trimmed to a count; `tile_read`/`tile_write` log lines now carry a
  `bytes=<N>` field, rendered as a new I/O throughput (MB/s) panel stacked
  directly below the swim-lane/thread panel, sharing its time axis so a
  dip or spike lines up with the Gantt bar above it; the "CPU stage" row
  (CPU thread-detail view) now includes the compute stage segment
  alongside mask/prep/cubestat, so the row sums to the tile's non-I/O time
  instead of silently omitting the largest of the four; and the GPU
  pipeline view's synchronous-fallback path (a tile that fits in one VRAM
  sub-block, so there's no async double-buffering to show) was fixed to
  recognize the log format the Fortran code actually emits — it had gone
  stale relative to an earlier logging format and was silently rendering
  no GPU lane at all for that case.
- **Full cfg reference in README**: every key the parser accepts,
  sectioned by purpose, each marked required/required-if/optional with
  its real default — cross-checked directly against
  `read_cfg_keyval`'s case statements, not written from memory.

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

### T6 validation: the same production case, re-run with `io_write_threads=8`

Same real ASKAP/EMU Q/U cubes, same 13308×11870×288 workload, same 16
tiles, now with `io_write_threads=8` added on top of the settings above:

| Stage | Old (`io_write_threads=1`) | New (`io_write_threads=8`) | Δ |
|---|---|---|---|
| I/O write | 2479.9s (96%) | 108.3s (6%) | **−95.6%** |
| CPU prep | 733.5s (28%) | 704.8s (36%) | −3.9% |
| CPU compute | 464.7s (18%) | 390.7s (20%) | −2.5%* |
| I/O read | 364.1s (14%) | 398.6s (20%) | +9.5% |
| CPU cubestat | 28.9s (1%) | 26.1s (1%) | −10% |
| **Total wall time** | **2586.7s** | **1945.4s** | **−24.9%** |

(*The plotter's own derived "CPU compute" figure runs slightly higher than
this — a long-standing, pre-existing gap between log-timestamp-derived
duration and the Fortran binary's internal stopwatch, unrelated to T6; see
`docs/ARCHITECTURE.md`.)

Write dropped by **~23x** — from 96% of wall time down to 6%, exactly the
result T6 was built for. `tile_prep`'s absolute cost barely moved (733.5s
→ 704.8s, if anything slightly lower), directly ruling out the plausible
worry that 8 background write threads would meaningfully steal cycles
from the concurrent compute pool during `io_overlap`'s overlap window —
read/write threads are I/O-bound (mostly blocked on the actual disk
operation), so the extra threads cost little even stacked on a
fully-subscribed compute pool (see `docs/PARALLELISM.md`, "Thread-pool
interplay"). `tile_prep`'s *share* rose only because the denominator
(total wall time) shrank so much everywhere else.

One consequence worth flagging for whoever picks up the next optimisation
cycle: `tile_prep` is now the single largest cost (36%), a new signal that
was invisible before because write dwarfed everything else. It doesn't
even appear as its own line in the binary's own "Macro timing breakdown"
summary — that summary only tracks four named stages (read/compute/
cubestat/write) and silently folds prep into a catch-all "other overhead"
bucket, which is now over 50% of wall time and mostly *is* prep. Not
something this release needed to fix, but worth knowing before someone
reads that summary and wonders where half the run went.

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
  on dev-machine test data (including combined with `io_overlap=y`) *and*
  on real Setonix production hardware — see "T6 validation" above for the
  measured ~23x write-time reduction and ~25% wall-time reduction on the
  same 13308×11870×288 ASKAP/EMU workload used throughout this cycle.

## What shipped in this tag

- `io_read_threads`, `io_write_threads`, `io_overlap` cfg keys (T0-T6),
  all defaulting to the pre-existing serial behaviour.
- The three production bugs above, found and fixed before or shortly
  after reaching the field.
- Swim-lane plotter: separate I/O read/write lanes, stage-totals bar
  panel, I/O throughput (MB/s) panel, complete CPU-stage row, working
  GPU sync-fallback rendering.
- A full, sectioned cfg reference in `README.md`.

## What's next (beyond this tag)

- `tile_prep` is now the largest single cost (36% of wall time, see "T6
  validation" above) — the natural next optimisation target, now that
  read and write have both been addressed. Not scoped for this release.
- Benchmark `io_read_threads=N` against the Lustre stripe count on
  Setonix to tune the recommended default (io_read_threads=8 was used
  throughout this cycle without a systematic sweep against stripe count).
- Add `tile_prep` as its own line in the Fortran binary's "Macro timing
  breakdown" summary, rather than leaving it folded into "other
  overhead" (see "T6 validation" above).
