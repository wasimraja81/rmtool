# overlap-host-device-compute: Execution Tickets

## Non-Negotiable Acceptance Criteria
1. Overlap exists on GPU path (CPU prep/scatter overlaps GPU compute for sub-block pipeline).
2. End-to-end efficiency improves on target workload versus baseline.
3. Numerical correctness is preserved (no regression in tests and output sanity checks).

A stage is failed if any one of the three criteria above is violated.

## Working Rules
- Keep each change set minimal and auditable.
- Every stage must run regression checks before closing.
- No stage is marked done without evidence captured in this file.
- No bluffing: if evidence is missing, status stays OPEN.

## Ticket Format
- ID:
- Objective:
- Scope:
- Change Set:
- Definition of Done (DoD):
- Required Regression:
- Evidence:
- Status: OPEN | IN_PROGRESS | DONE | BLOCKED

## Micro-Step Format (Use under each ticket)
- Step:
- What this is meant to achieve:
- Evidence produced:
- What next is needed to meet final goal:
- Result: PASS | FAIL | PARTIAL

---

## T0 - Baseline Lock
- ID: T0
- Objective: Freeze a reproducible baseline for correctness and timing.
- Scope: Current branch HEAD, no behavior change.
- Change Set: None (measurement only).
- Runbook (Operator View):
  - Step: Build baseline GPU binary (OMP=0, GPU=1).
  - What this is meant to achieve: lock a reproducible baseline executable.
  - Evidence produced: build log.
  - What next is needed to meet final goal: run regressions and collect baseline timing/sanity.
  - Step: Run regression and sanity checks.
  - What this is meant to achieve: confirm baseline correctness before optimization.
  - Evidence produced: tests/run_tests.sh summary and AMP/PHA/NVALID sanity report.
  - What next is needed to meet final goal: capture baseline wall time for later speedup comparison.
  - Step: Run baseline workload with timing capture.
  - What this is meant to achieve: establish reference runtime.
  - Evidence produced: scratch/T0.baseline.omp0.gpu1.log metrics.
  - What next is needed to meet final goal: execute T1 overlap-enabling changes and compare.
- Definition of Done (DoD):
  - Baseline build command and run command are recorded.
  - Baseline timing for target workload is recorded.
  - Baseline correctness checks recorded.
- Required Regression:
  - Build: make OMP=0 GPU=1
  - Tests: tests/run_tests.sh
  - Output sanity: finite/nonzero/min-max checks for AMP/PHA/NVALID
- Evidence:
  - Build: PASS (make OMP=0 GPU=1)
  - Tests: PASS (tests/run_tests.sh => Total 22, Pass 22, Fail 0)
  - Timing: PASS (scratch/T0.baseline.omp0.gpu1.log => wall 3:41.17, user 132.34s, sys 88.88s, max RSS 5336596 KB)
  - Sanity: PASS (AMP/PHA/NVALID finite and nonzero; AMP max 0.3648194, PHA range [-pi, +pi], NVALID min=max=167)
- Status: DONE

## T1 - Async Preconditions and Tasking Validity
- ID: T1
- Objective: Ensure async path executes in a valid host tasking context and is not accidentally serialized.
- Scope: Task region setup and task dependency plumbing only.
- Change Set: src/rm_synthesis.f90
- Runbook (Operator View):
  - Step: Remove per-subblock wait in rm_synthesis.f90:3018.
  - What this is meant to achieve: allows CPU prep/scatter of the next subblock to proceed while prior GPU work is still running, which is required for overlap.
  - Evidence produced: code diff and successful compile.
  - What next is needed to meet final goal: verify correctness did not regress.
  - Step: Keep only end-of-tile synchronization in rm_synthesis.f90:3074.
  - What this is meant to achieve: preserves safety boundary while avoiding unnecessary serialization inside the pipeline.
  - Evidence produced: code diff and no race/corruption observed in outputs.
  - What next is needed to meet final goal: rebuild all binaries from clean state.
  - Step: Cold rebuild all four flavors.
  - What this is meant to achieve: guarantees test and timing results are from current code only.
  - Evidence produced: clean-all plus rebuild logs showing recompilation of main sources.
  - What next is needed to meet final goal: run full regression suite.
  - Step: Run full regression tests.
  - What this is meant to achieve: proves overlap change did not break correctness.
  - Evidence produced: 22/22 pass report.
  - What next is needed to meet final goal: run real workload timing for overlap and speedup proof.
  - Step: Run real workload timing and overlap analysis.
  - What this is meant to achieve: proves whether overlap is actually happening and whether runtime improved vs T0 baseline.
  - Evidence produced: wall time comparison and overlap markers from logs.
  - What next is needed to meet final goal: if speedup is insufficient, iterate on dependency/wait placement; if speedup and correctness both pass, goal is met.
- Definition of Done (DoD):
  - Async path executes with real OpenMP task dependencies.
  - No unconditional in-loop wait that kills run-ahead.
  - End-of-tile synchronization remains for safety.
- Required Regression:
  - Build: make OMP=0 GPU=1 and make OMP=1 GPU=1
  - Tests: tests/run_tests.sh
  - Output sanity: AMP/PHA/NVALID checks
- Evidence:
  - Build: PASS (cold rebuild all flavors after T1.S3 refactor)
  - Tests: PASS (tests/run_tests.sh => Total 22, Pass 22, Fail 0)
  - Timing: PASS (scratch/T1.S3.slotlocal.omp1.gpu1.log => wall 2:58.04, exit status 0)
  - Sanity: PASS (AMP/PHA/NVALID finite/nonzero; ranges consistent with baseline)
- Status: IN_PROGRESS

### Micro-Steps
#### T1.S1
- Step: Establish valid OpenMP tasking lexical region around async sub-block loop.
- What this is meant to achieve: ensure task directives are executed in a valid tasking context so async behavior is structurally possible instead of implicitly serialized by invalid placement.
- Evidence produced:
  - Code update in src/rm_synthesis.f90 introducing omp parallel + omp single region around the sub-block task loop.
  - Cold rebuild of all four flavors succeeded after change:
    - make OMP=0 GPU=0
    - make OMP=1 GPU=0
    - make OMP=0 GPU=1
    - make OMP=1 GPU=1
  - Regression run succeeded: tests/run_tests.sh => Total 22, Pass 22, Fail 0.
- What next is needed to meet final goal: remove per-subblock taskwait to permit run-ahead overlap while preserving end-of-tile safety wait; then rerun cold rebuild + full regression + timing proof.
- Result: PASS

#### T1.S2 (Planned next)
- Step: Remove in-loop per-subblock taskwait and retain only end-of-tile synchronization boundary.
- What this is meant to achieve: allow CPU prep/scatter and GPU compute to overlap across adjacent sub-blocks (core overlap requirement).
- Evidence produced:
  - Code diff applied in src/rm_synthesis.f90.
  - Cold rebuild all flavors passed.
  - Full regression passed: tests/run_tests.sh => Total 22, Pass 22, Fail 0.
  - Real workload failed with runtime GPU errors in scratch/T1.after_wait_removal.omp1.gpu1.log:
    - libgomp/cuCtxSynchronize illegal memory access
    - libgomp/cuMemcpyDtoH_v2 illegal memory access
    - Exit status 1
- What next is needed to meet final goal: restore safety wait, revalidate known-good state, then redesign overlap to avoid reusing/freeing device-side resources while async tasks are still in flight.
- Result: FAIL

#### T1.S2R (Recovery)
- Step: Reintroduce per-subblock taskwait safety gate after failed real-run attempt.
- What this is meant to achieve: restore correctness and eliminate illegal memory access while keeping branch in a stable state.
- Evidence produced:
  - Safety wait restored in src/rm_synthesis.f90.
  - Cold rebuild all four flavors succeeded.
  - Full regression passed: tests/run_tests.sh => Total 22, Pass 22, Fail 0.
- What next is needed to meet final goal: implement overlap with correct resource lifetime boundaries (defer slot-local deallocation/reuse until dependent compute+copy/scatter completion), then repeat rebuild+regression+real timing.
- Result: PASS

#### T1.S3 (Executed)
- Step: Move resource release/reuse boundary so slot-local buffers are not touched until dependency-complete, then retry controlled overlap.
- What this is meant to achieve: prevent illegal memory access while enabling true run-ahead overlap.
- Evidence produced:
  - Refactor in src/rm_synthesis.f90 to use slot-local staging GPU temporaries (slot1/slot2) for prepare, compute, and deallocation paths.
  - Cold rebuild all four flavors succeeded.
  - Full regression passed: tests/run_tests.sh => Total 22, Pass 22, Fail 0.
  - Real workload completed successfully (scratch/T1.S3.slotlocal.omp1.gpu1.log):
    - no illegal memory access messages
    - exit status 0
    - wall time 2:58.04
  - Output sanity checks passed for MY_CASA_RMSYNTH_FULLIM_TEST.{AMP,PHA,NVALID}.
- What next is needed to meet final goal: retry controlled removal of in-loop wait using the slot-local buffers as safety base, then repeat cold rebuild + regression + real-run validation and finally quantify overlap/speedup versus T0.
- Result: PARTIAL

#### T1.S4 (Next)
- Step: Re-attempt in-loop wait relaxation on top of slot-local buffers, with immediate rollback on any real-run failure.
- What this is meant to achieve: enable actual run-ahead overlap while preserving correctness.
- Evidence produced: pending.
- What next is needed to meet final goal: if stable, extract overlap markers and speedup; if unstable, redesign dependency boundaries before further attempts.
- Result: PARTIAL

#### T1.SD1 (Diagnostics)
- Step: Add async lifecycle logging (enqueue/start/done/wait/dealloc) and run real workload with current safe path.
- What this is meant to achieve: replace guesswork with exact ordering evidence for slot lifecycle and synchronization boundaries.
- Evidence produced:
  - Instrumentation added in src/rm_synthesis.f90 under `tile_async` debug category.
  - Build + regression pass after instrumentation: tests/run_tests.sh => Total 22, Pass 22, Fail 0.
  - Real run completed successfully: scratch/T1.diag.async_lifecycle.log => Exit status 0, wall 2:57.46.
  - Ordering evidence from MY_CASA_RMSYNTH_FULLIM_TEST.run.log shows:
    - `async wait begin` occurs immediately after compute enqueue for each slot.
    - `async wait end` occurs only after compute task completion.
    - deallocation/scatter execute after wait.
    - next subblock send for slot 2 starts only after scatter of slot 1 completes.
- What next is needed to meet final goal: redesign task graph so per-slot dealloc/scatter are dependency-driven (not wait-driven) and slot reuse waits on slot completion token, then re-run full gates.
- Result: PASS

#### T1.SD2 (Dependency-Chain + Task-Local)
- Step: Convert async path to dependency-ordered compute->scatter/dealloc tasks and localize task loop variables.
- What this is meant to achieve: remove per-subblock blocking wait while preventing illegal memory access from host/device lifetime races.
- Evidence produced:
  - Refactor in src/rm_synthesis.f90:
    - compute task uses `depend(in:dep_h2d(slot))` -> `depend(out:dep_kern(slot))`
    - scatter/dealloc task uses `depend(in:dep_kern(slot))` -> `depend(out:dep_h2d(slot))`
    - task loop/control vars made task-local via `private(...)`
  - Full cold rebuild + regression pass: tests/run_tests.sh => Total 22, Pass 22, Fail 0.
  - Real workload run pass: scratch/T1.SD2.depchain_tasklocal.log => Exit status 0, wall 2:59.09, no illegal-memory-access messages.
  - Output sanity pass for MY_CASA_RMSYNTH_FULLIM_TEST.{AMP,PHA,NVALID}.
  - Lifecycle logs in MY_CASA_RMSYNTH_FULLIM_TEST.run.log show true overlap ordering:
    - slot2 compute starts before slot1 scatter/dealloc completes.
- What next is needed to meet final goal:
  - Verify numerical equivalence against known-good safe path on real workload (not just finite/nonzero sanity).
  - Fix or reinterpret timing-stage accounting under async tasks (tile_compute dropped to implausibly low value while wall time stayed high).
  - Demonstrate net speedup versus T0 baseline; currently wall time did not improve.
- Result: PARTIAL

## T2 - Correct Double-Buffer Ownership
- ID: T2
- Objective: Guarantee per-slot ownership for staged input/output and safe slot reuse.
- Scope: slot-local arrays and token transitions only.
- Change Set: src/rm_synthesis.f90
- Runbook (Operator View):
  - Step: Enforce slot-local ownership for in-flight input/output buffers.
  - What this is meant to achieve: prevent cross-slot overwrite hazards during overlap.
  - Evidence produced: code diff and slot-indexed dataflow review.
  - What next is needed to meet final goal: verify slot reuse waits for completion token.
  - Step: Gate slot reuse on completion dependency.
  - What this is meant to achieve: avoid reading/writing a slot before compute/scatter is done.
  - Evidence produced: dependency-chain diff and passing regression.
  - What next is needed to meet final goal: measure runtime impact in T4.
- Definition of Done (DoD):
  - Input/output buffers are slot-local for in-flight overlap.
  - Scatter reads only from completed slot.
  - Slot reuse is gated by completion dependency.
- Required Regression:
  - Build: make OMP=0 GPU=1 and make OMP=1 GPU=1
  - Tests: tests/run_tests.sh
  - Output sanity: AMP/PHA/NVALID checks
- Evidence:
  - Build:
  - Tests:
  - Timing:
  - Sanity:
- Status: OPEN

### Micro-Steps
- None logged yet.

## T3 - Enable Overlap Policy on GPU Path
- ID: T3
- Objective: Align gating policy with branch goal so GPU path can overlap when capability exists.
- Scope: use_async_pipeline policy and guard conditions.
- Change Set: src/rm_synthesis.f90, Makefile (if macro policy needs adjustment)
- Runbook (Operator View):
  - Step: Change async gate from flavor-only policy to capability-based policy.
  - What this is meant to achieve: permit overlap path when runtime conditions support it.
  - Evidence produced: gate-condition diff and mode coverage test.
  - What next is needed to meet final goal: confirm fallback path remains correct.
  - Step: Validate fallback behavior when overlap conditions are not met.
  - What this is meant to achieve: preserve correctness across all build/runtime variants.
  - Evidence produced: full regression pass across variants.
  - What next is needed to meet final goal: performance proof in T4.
- Definition of Done (DoD):
  - Gate reflects capability, not conservative flavor naming only.
  - GPU run with sub-blocks can enter overlap path when safe.
  - Fallback path remains correct when overlap conditions are not met.
- Required Regression:
  - Build: make OMP=0 GPU=1 and make OMP=1 GPU=1
  - Tests: tests/run_tests.sh
  - Output sanity: AMP/PHA/NVALID checks
- Evidence:
  - Build:
  - Tests:
  - Timing:
  - Sanity:
- Status: OPEN

### Micro-Steps
- None logged yet.

## T4 - Performance Proof
- ID: T4
- Objective: Prove measurable overlap and speedup.
- Scope: instrumentation/log parsing and before-vs-after comparison.
- Change Set: src/rm_synthesis.f90 and/or scripts used for analysis
- Runbook (Operator View):
  - Step: Run target workload with overlap-enabled build.
  - What this is meant to achieve: collect comparable timing data against T0.
  - Evidence produced: run log with wall/user/sys and stage markers.
  - What next is needed to meet final goal: quantify overlap metric.
  - Step: Compute overlap metric and speedup vs baseline.
  - What this is meant to achieve: prove both overlap existence and efficiency gain.
  - Evidence produced: overlap ratio and runtime delta summary.
  - What next is needed to meet final goal: ensure correctness still passes after performance run.
- Definition of Done (DoD):
  - Timeline/log shows CPU prep/scatter overlapping GPU compute.
  - Overlap metric is reported.
  - End-to-end runtime improvement is reported on same workload.
- Required Regression:
  - Build: make OMP=0 GPU=1 and make OMP=1 GPU=1
  - Tests: tests/run_tests.sh
  - Output sanity: AMP/PHA/NVALID checks
- Evidence:
  - Build:
  - Tests:
  - Timing:
  - Sanity:
  - Overlap metric:
  - Speedup:
- Status: OPEN

### Micro-Steps
- None logged yet.

## T5 - Stabilization and Merge Readiness
- ID: T5
- Objective: Final hardening with clean evidence trail.
- Scope: docs, logs, and final rerun.
- Change Set: TODO.double-buffer-overlap.md + minimal code/doc deltas
- Runbook (Operator View):
  - Step: Final cold rebuild and final full regression rerun.
  - What this is meant to achieve: confirm reproducibility of correctness after all changes.
  - Evidence produced: clean rebuild logs and final test summary.
  - What next is needed to meet final goal: finalize risk register and merge notes.
  - Step: Final timing rerun and closeout report.
  - What this is meant to achieve: confirm speedup claim is reproducible.
  - Evidence produced: final timing comparison and residual risk list.
  - What next is needed to meet final goal: mark all tickets DONE and prepare merge.
- Definition of Done (DoD):
  - All prior tickets are DONE with evidence.
  - Final rerun reproduces correctness and speedup.
  - Remaining risks are explicitly listed.
- Required Regression:
  - Build: make OMP=0 GPU=1 and make OMP=1 GPU=1
  - Tests: tests/run_tests.sh
  - Output sanity: AMP/PHA/NVALID checks
- Evidence:
  - Build:
  - Tests:
  - Timing:
  - Sanity:
  - Residual risks:
- Status: OPEN

### Micro-Steps
- None logged yet.

---

## Stage-by-Stage Execution Order
1. T0 baseline lock
2. T1 tasking validity
3. T2 buffer ownership safety
4. T3 gating policy alignment
5. T4 overlap + speed proof
6. T5 stabilization

## Failure Policy
- If any regression fails, ticket status changes to BLOCKED.
- Root cause and rollback/fix decision must be documented before next ticket starts.
