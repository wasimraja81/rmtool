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

## Work Tickets

### T0 - Baseline Lock
- Objective: capture baseline timing and correctness on `optimise-io`.
- Scope: measurement only.
- Evidence to capture:
  - benchmark command + environment,
  - stage timings including `tile_read`/`tile_write`,
  - tests and output sanity.

### T1 - IO Isolation Layer (No Behaviour Change)
- Objective: isolate read/write call sequences into helper routines.
- Scope: refactor only; preserve behaviour.
- Evidence to capture:
  - unchanged outputs for same workload,
  - full test pass,
  - unchanged stage marker semantics.

### T2 - Parallel Read MVP
- Objective: reduce read bottleneck safely.
- Scope: read path only.
- Plan:
  1. first implement read-ahead overlap (producer-consumer double buffer),
  2. then optionally parallelize independent Q/U/(I/mask) reads with separate
     file handles if runtime/library validation is clean.
- Evidence to capture:
  - no corruption/races,
  - reduced read bottleneck in timings.

### T3 - Parallel/Async Write MVP
- Objective: reduce write bottleneck safely.
- Scope: write path only.
- Plan:
  1. first implement async write overlap with compute,
  2. optionally parallelize writes across independent output files while keeping
     each file stream ordered.
- Evidence to capture:
  - valid FITS outputs,
  - deterministic results,
  - improved end-to-end runtime.

### T4 - Integrated IO Overlap Policy
- Objective: combine read-ahead + async write with robust fallback controls.
- Scope: IO orchestration policy only.
- Evidence to capture:
  - long-run stability,
  - switchable fallback to serial path,
  - DoD performance/correctness gates met.

## Answers to kickoff questions
1. Can we read FITS cubes in parallel?
   - Yes, with caveats.
   - Safest immediate gain: read-ahead overlap by tile.
   - Parallel Q/U/I/mask reads are possible with independent handles and
     validated thread-safety on deployed CFITSIO build.
   - Avoid concurrent operations on the same file handle.
2. Can we write FITS in parallel?
   - Yes, with caveats.
   - Safest immediate gain: async write overlap with compute.
   - Parallel writes are safest across independent output products; keep
     in-file ordering deterministic.
3. How do we avoid breaking compute/staging?
   - Keep compute and staging routines unchanged.
   - Implement IO-only orchestration wrappers and feature flags.
   - Require full regressions and benchmark evidence before default enable.

## Easy-gain priority order
1. Read-ahead next tile during current tile compute.
2. Deferred/async write of completed tile during next tile compute.
3. Parallel Q/U (and optional I/mask) reads with separate handles.
4. Parallel writes across independent output files.
5. Extra IO telemetry: bytes, queue depth, and wait reasons.
