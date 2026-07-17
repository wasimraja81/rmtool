# IO Parallel Optimisation Plan

## Context from Setonix swim-lane and runtime logs
- Pipeline timeline views show `I/O read` and `I/O write` taking a substantial
  portion of wall time on large runs.
- CPU thread-detail views indicate extraction threads are mostly productive once
  compute starts, but tile-boundary IO introduces visible serial gaps.
- Current tile loop behaviour is serial for FITS IO in
  `src/rm_synthesis.f90`:
  - tile read uses `FTGSVE` calls,
  - tile write uses `ftpsse`/`ftpssi`/`ftpssb` calls.

Conclusion: easiest wins are in IO overlap/scheduling, not in RM math kernels.

## Non-Negotiable Guardrails
1. Do not change numerical compute kernels.
2. Do not break staging, RM chunking, or tile semantics.
3. Keep new behaviour behind explicit runtime flags with safe serial fallback.

## Definition of Done
1. Correctness
   - `tests/run_tests.sh` remains fully passing.
   - Output sanity checks stay valid (AMP/PHA/NVALID and optional cubestat maps).
2. Isolation
   - IO changes remain in IO orchestration/helpers only.
   - No algorithmic change in extraction or GPU/CPU compute kernels.
3. Observability
   - Existing `tile_read` and `tile_write` timings remain visible.
   - New logs identify serial vs overlapped/parallel IO mode.
4. Performance
   - Setonix benchmark wall time improves versus baseline.
   - `tile_read + tile_write` impact is reduced (absolute or through overlap).
5. Operability
   - Runtime switch can disable new IO path immediately.
   - Fallback path is tested and documented.
6. Setonix Closure
   - Prior failing large-tile case no longer crashes in tile read.
   - Results remain numerically consistent with baseline outputs.
   - IO timings are reported before/after on identical workload.

## Setonix bug discovered during testing
- Symptom: crash in `FTGSVE` tile-read path on first large subset read in
  Setonix CPU run.
- Trigger profile: very large single-call read extent where per-call element
  count exceeds wrapper-level signed 32-bit safety bounds.
- Impact: run abort before compute completes, even when host memory is
  sufficient.
- Scope lock: compute kernels remain unchanged; remediation is IO interface
  and orchestration only.
- Agreed direction on this branch:
  1. implement a 64-bit-capable read interface path first,
  2. use call-level read chunking only as fallback if ABI/symbol constraints
     block direct 64-bit interface usage,
  3. keep serial fallback switchable at runtime.

## Work Tickets

### T0 - Baseline Lock ✓ DONE
- Objective: capture baseline timing and correctness on `optimise-io`.
- Scope: measurement only.
- Evidence captured:
  - Jennifer full-image run on Setonix used as baseline reference.
  - Stage timings in `tile_read`/`tile_write` CSV and log confirmed IO-dominant profile.

### T1 - IO Isolation Layer
- Objective: isolate read/write call sequences into helper routines.
- Status: deferred — T2 and T3/T4 were implemented directly; isolation layer adds
  refactor cost with limited incremental benefit given the current approach.

### T2 - Setonix read-interface hardening ✓ DONE
- Objective: remove large-read overflow risk without changing compute path.
- Commit: `0ac331f` (`large-tile-read-safe` tag)
- Delivered:
  - All `allocate()` size expressions use `int(...,kind=int64)` casts.
  - FTGSVE single-call path preserved; correct with CFITSIO 4.x (Ubuntu 24.04).
  - Crash on 13 308 × 734 × 288 tile eliminated.

### T3 - Parallel Channel Reads ✓ DONE
- Objective: reduce read wall time via parallel CFITSIO handles.
- Commit: `a072c6e`
- Delivered:
  - `io_read_threads` cfg key (default 1, backward-compatible).
  - N independent read-only handles per input file; OMP parallel do over channel
    chunks; each thread writes a contiguous slice of the pre-allocated tile buffer.
  - Clamped to `omp_get_max_threads()` and `nz_out`; serial binary forced to 1.
  - 22/22 tests pass.

### T4 - Parallel RM-bin Writes ✓ DONE
- Objective: reduce write wall time via parallel CFITSIO handles on output cubes.
- Commits: `885f21e` (feature), `904f75d` (OMP ceiling check)
- Delivered:
  - `io_write_threads` cfg key (default 1, backward-compatible).
  - N independent readwrite handles per output cube (AMP/PHA); OMP parallel do over
    RM-bin chunks; non-overlapping pwrite is POSIX-safe.
  - 2D map outputs (NVALID, MASK, cubestat) remain serial.
  - Clamped to `omp_get_max_threads()` and `nrm_out`; serial binary forced to 1.
  - 22/22 tests pass.

### T5 - Async Read/Write Overlap (next)
- Objective: overlap tile N+1 read with tile N write and tile N+1 compute.
- Scope: IO orchestration policy only.
- Plan:
  1. Double-buffer the tile output arrays (p_tile_arr / phi_tile_arr) so tile N
     write can proceed concurrently with tile N+1 read.
  2. Producer-consumer pipeline: read → compute → write in overlapping phases.
  3. Barrier before tile N+1 compute to ensure tile N write is complete.
- Evidence to capture:
  - long-run stability,
  - switchable fallback to serial path,
  - DoD performance/correctness gates met.

## Answers to kickoff questions
1. Can we read FITS cubes in parallel?
   - Yes — implemented via `io_read_threads`.
   - N independent read-only handles per input file read disjoint channel ranges
     simultaneously into the pre-allocated tile buffer.
   - Channel ranges map to different Lustre OSTs; true parallel IO throughput.
   - `io_read_threads=1` (default) is identical to the prior serial single-call path.
2. Can we write FITS in parallel?
   - Yes — implemented via `io_write_threads`.
   - N independent readwrite handles per output cube write disjoint RM-bin ranges
     simultaneously from the scatter output buffer.
   - Non-overlapping `pwrite()` is POSIX-safe; CFITSIO metadata updated only at close.
   - `io_write_threads=1` (default) is identical to the prior serial write path.
   - 2D map outputs (NVALID, MASK, cubestat) remain serial (negligible cost).
3. How do we avoid breaking compute/staging?
   - Compute and staging routines are unchanged; IO change is orchestration-only.
   - Both cfg keys default to 1 (serial); existing behaviour is the default.
   - Both keys are clamped to `omp_get_max_threads()` and their data-depth limit.
   - Serial binary (HOST_OMP=0) forces both to 1 automatically.

## Thread-count design rationale

`io_read_threads` and `io_write_threads` are separate from `OMP_NUM_THREADS`
because compute and IO are constrained by different resources:

| | Compute (`OMP_NUM_THREADS`) | IO threads |
|---|---|---|
| Bottleneck | CPU cores | Lustre OST stripe count |
| Optimal N | All available (e.g. 128) | File stripe count (4–16 typical) |
| Over-count cost | Negligible | Excess handles, metadata pressure |

Recommended values: set `io_read_threads` and `io_write_threads` to the
Lustre stripe count of the data files (`lfs getstripe <file>` on Setonix).

## Current status and next steps

| Ticket | Status | Commit |
|---|---|---|
| T0 Baseline lock | Done | — |
| T1 IO isolation layer | Deferred | — |
| T2 Large-tile read hardening | Done | `0ac331f` (`large-tile-read-safe`) |
| T3 Parallel channel reads | Done | `a072c6e`, `904f75d` |
| T4 Parallel RM-bin writes | Done | `885f21e`, `904f75d` |
| T5 Async read/write overlap | Next | — |

Pending validation: benchmark on Setonix with `io_read_threads=N` and
`io_write_threads=N` (where N = Lustre stripe count) vs baseline
`io_read_threads=1` to quantify wall-time improvement.
