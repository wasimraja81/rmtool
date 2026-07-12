# TODO: Real Double-Buffer Overlap (CPU Prep + GPU Compute)

## Goal
Implement true double buffering in the staged GPU path so CPU prep/scatter and GPU compute overlap in time for consecutive sub-blocks.

## Current Problem Statement
From timing traces (swimlane), CPU and GPU are still largely serialized:
- CPU prep/scatter waits for GPU compute completion.
- Overlap is not visible in the first row-band timeline.

## Accountability Rules
- No task is considered complete without:
  - Code change committed in `src/rm_synthesis.f`
  - Successful build with `make OMP=0 GPU=1`
  - Timing evidence showing expected overlap behavior
- Any temporary safety fence (barrier/taskwait) must be documented with reason and removal plan.

## Scope (This TODO)
- Staged two-level GPU path only (`use_staging` path).
- Focus on one real scenario first:
  - `use_gpu_actual = .true.`
  - `n_subblocks_tile > 1`

## Non-Goals (for now)
- Full redesign of non-staging CPU-only path.
- Multi-GPU scheduling policy changes.
- Global algorithmic changes to RM extraction math.

## Execution Steps (Primary Tracking)

### S1. Remove Serialization Points
- [ ] S1.1 Remove in-loop blocking wait that prevents CPU run-ahead.
- [ ] S1.2 Keep only a safe synchronization boundary at end-of-tile (or equivalent correctness boundary).
- [ ] S1.3 Verify no hidden host-side waits remain in gather/compute/scatter control flow.

### S2. Make Buffer Ownership Truly Double-Buffered
- [ ] S2.1 Ensure input staging is slot-local for both current and next slot (`slot_idx`, `next_slot`).
- [ ] S2.2 Make output staging slot-local for in-flight overlap safety (no shared output hazards across sub-blocks).
- [ ] S2.3 Ensure scatter reads only from the slot that is logically complete.

### S3. Make Dependency Wiring Functional (Not Cosmetic)
- [ ] S3.1 `use_async_pipeline` branches must have different behavior for async vs sync mode.
- [ ] S3.2 H2D token chain (`dep_h2d`) must represent readiness of per-slot input.
- [ ] S3.3 Kernel token chain (`dep_kern`) must represent completion of per-slot compute.
- [ ] S3.4 D2H/scatter token chain (`dep_d2h`) must gate slot reuse correctly.
- [ ] S3.5 No if/else branch may contain identical logic under `use_async_pipeline`.

### S4. Validate Execution Semantics
- [ ] S4.1 Confirm OpenMP tasks execute in a valid tasking region (not serialized by context).
- [ ] S4.2 Confirm offload calls used in async path are not forcing immediate host synchronization.
- [ ] S4.3 Confirm slot reuse cannot occur before dependent stages complete.

### S5. Build and Correctness Gates
- [ ] S5.1 Build passes: `make OMP=0 GPU=1`
- [ ] S5.2 No fixed-form formatting regressions (continuation column correctness, no line truncation bugs).
- [ ] S5.3 Existing tests/smoke run for output sanity pass (or clearly documented if unavailable).

### S6. Performance Evidence Gates
- [ ] S6.1 Produce updated swimlane timeline for same workload as baseline.
- [ ] S6.2 Demonstrate visible overlap between CPU prep/scatter and GPU compute.
- [ ] S6.3 Quantify overlap metric (example):
  - overlap_ratio = overlap_time / gpu_compute_time
- [ ] S6.4 Compare before vs after wall time for first row-band and whole run.

## Legacy Checklist Mapping

Legacy labels were kept for historical context only. Use S1-S6 above as the canonical tracker.

### A. Remove Serialization Points (maps to S1)
- [ ] A1. Remove in-loop blocking wait that prevents CPU run-ahead.
- [ ] A2. Keep only a safe synchronization boundary at end-of-tile (or equivalent correctness boundary).
- [ ] A3. Verify no hidden host-side waits remain in gather/compute/scatter control flow.

### B. Make Buffer Ownership Truly Double-Buffered (maps to S2)
- [ ] B1. Ensure input staging is slot-local for both current and next slot (`slot_idx`, `next_slot`).
- [ ] B2. Make output staging slot-local for in-flight overlap safety (no shared output hazards across sub-blocks).
- [ ] B3. Ensure scatter reads only from the slot that is logically complete.

### C. Dependency Wiring Must Be Functional, Not Cosmetic (maps to S3)
- [ ] C1. `use_async_pipeline` branches must have different behavior for async vs sync mode.
- [ ] C2. H2D token chain (`dep_h2d`) must represent readiness of per-slot input.
- [ ] C3. Kernel token chain (`dep_kern`) must represent completion of per-slot compute.
- [ ] C4. D2H/scatter token chain (`dep_d2h`) must gate slot reuse correctly.
- [ ] C5. No if/else branch may contain identical logic under `use_async_pipeline`.

### D. Execution Semantics Validation (maps to S4)
- [ ] D1. Confirm OpenMP tasks execute in a valid tasking region (not serialized by context).
- [ ] D2. Confirm offload calls used in async path are not forcing immediate host synchronization.
- [ ] D3. Confirm slot reuse cannot occur before dependent stages complete.

### E. Build and Correctness Gates (maps to S5)
- [ ] E1. Build passes: `make OMP=0 GPU=1`
- [ ] E2. No fixed-form formatting regressions (continuation column correctness, no line truncation bugs).
- [ ] E3. Existing tests/smoke run for output sanity pass (or clearly documented if unavailable).

### F. Performance Evidence Gates (maps to S6)
- [ ] F1. Produce updated swimlane timeline for same workload as baseline.
- [ ] F2. Demonstrate visible overlap between CPU prep/scatter and GPU compute.
- [ ] F3. Quantify overlap metric (example):
  - overlap_ratio = overlap_time / gpu_compute_time
- [ ] F4. Compare before vs after wall time for first row-band and whole run.

## Definition of Done
Done means all of the following are true:
1. Steps S1-S6 completed.
2. Updated timeline clearly shows non-trivial CPU/GPU overlap.
3. Build is green with `make OMP=0 GPU=1`.
4. No correctness regression observed in sanity outputs.

## Evidence Log Template
Use this section for each implementation step:

### Step <N>: <title>
- Change summary:
- Files touched:
- Build result:
- Timing result:
- Correctness check result:
- Remaining risk:

## Baseline Snapshot (to beat)
- Symptom: CPU and GPU lanes alternate with little/no overlap in first row-band.
- Approximate first row-band duration: ~223 s (from supplied chart).

## Risk Register
- R1: Hidden synchronization in offload runtime may serialize tasks.
- R2: Shared buffers may still cause conservative waits.
- R3: Fixed-form formatting mistakes can silently alter logic.

## Immediate Next Step
- Implement S1.1 + S2.2 in the smallest safe patch, then rebuild and re-profile.
