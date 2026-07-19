# Encapsulation Refactor Plan — `rm_synthesis.f90` / `rm_synthesis_mod.f90`

Branch: `encapsulate-rmsynth` (from `main` @ tag `3.1`)

## Context

The codebase works and is well-tested (28 tests, including bit-identical
output checks and structural concurrency-invariant checks), but two files
carry almost all the logic and have grown unwieldy: `src/rm_synthesis.f90`
(4273 lines, 237 top-level declarations in one flat scope) and
`src/rm_synthesis_mod.f90` (3205 lines). The goal is to hide complexity
behind derived types and module boundaries — modern Fortran supports this
well — **without changing any functional behaviour**. This is an
incremental encapsulation effort, not a rewrite: every ticket is scoped so
its own correctness can be verified before moving to the next, using the
existing test suite as the primary gate, reinforced with bit-identical
output re-checks where a ticket's risk warrants it.

**Hard constraint carried through every ticket below:** zero change in
observable behaviour (bit-identical output, identical log content where
it's user-facing e.g. `tile_autotune.cfg`, identical pass/fail outcome of
all 28 tests). A ticket that can't be verified to meet this is rolled back,
not pushed through.

**Numerical kernels are out of scope, permanently**: `extract_general*`,
`tile_extract_gpu_rm_blocked`, `cubestat_tail_quantile_maps` in
`rm_synthesis_mod.f90` are not touched by any ticket below. They're
already validated and encapsulating them further buys nothing but risk.

**Scope decisions (confirmed with user):**
- GPU staging/async pipeline (dependency-token/slot bookkeeping,
  `rm_synthesis.f90:3284-3826`) is **deferred entirely, out of scope for
  this effort**. This code has a documented history of real production
  SIGSEGVs during its original development (`TODO/TODO.double-buffer-overlap.md`,
  local-only) and three separate past P0 data-corruption bugs in this
  exact region (`TODO/CRITICAL_ASSESSMENT.md`, local-only), and today's
  test coverage for it is weaker than for the I/O write path — no
  structural race-detection test exists for it analogous to `run_tests.sh`
  §13. Revisit only later, and only after building a dedicated structural
  test for it first.
- Deep cfg threading (replacing all ~237 declarations and ~500+ usage
  sites across `rm_synthesis.f90` with `cfg%` field access) **is
  included**, as the final, lowest-priority ticket (T5), done only after
  T1-T3 are stable, split into sub-tickets by variable cluster.

## Tickets

Ticket format follows this repo's existing convention
(`planning/IO_PARALLEL_OPTIMISATION_PLAN.md`): Objective / Scope / Change
Set / Correctness Gate / Rollback Criteria / Effort.

---

### T0 — Baseline Lock — DONE
- **Objective:** Freeze an exact reference before any code changes, so every later ticket has something concrete to diff against.
- **Scope:** Measurement only.
- **Change Set:** None.
- **Correctness Gate:**
  - `bash scratch/make_all.sh` clean build, actual warning count recorded (don't assume zero — record what's really there today).
  - `bash tests/run_tests.sh` full run, PASS/FAIL/SKIP counts recorded verbatim.
  - For every `tests/*.cfg` (and any real-workload cfg available), run once and archive the full output FITS set plus run log under `scratch/baseline_encapsulation/` — the bit-identical reference for every later ticket.
- **Rollback Criteria:** N/A.
- **Effort:** 0.5 session.
- **Evidence (2026-07-19, commit `043c7cf`):** Build clean, 0 compiler
  errors/warnings, 4 pre-existing GPU-offload linker warnings (unchanged
  from every prior build this cycle). `tests/run_tests.sh`: 28/28 pass.
  174 output files (FITS/cfg/log/csv, 42M) archived to
  `scratch/baseline_encapsulation/tests_output/` (gitignored, local-only
  — see `scratch/baseline_encapsulation/T0_MANIFEST.md` for the full
  record and exact commit hash, reproducible from there). No real
  Setonix-scale workload available locally — T3b's large-case re-run
  will need to happen on Setonix directly.

---

### T1 — Config Encapsulation (shallow: call-site wrapper only) — DONE
- **Objective:** `read_cfg_keyval` returns one `type(rmsynth_config_t)` instead of populating ~56 separate `intent(inout)`/`intent(out)` arguments. The struct is unpacked into the existing loose locals immediately after the call, in the same place the call already is — nothing downstream in the 4000+ lines changes.
- **Scope:** `rm_synthesis_mod.f90:1519-2630` (`read_cfg_keyval` signature and its internal assignments only — the `select case` body, duplicate detection, and cross-key validation at `2580-2628` are untouched, just renamed `var` → `cfg%var` inside the subroutine). Call site `rm_synthesis.f90:447-471` gets a mechanical unpack block right after the call.
- **Change Set:** New `rmsynth_config_t` derived type, one field per current argument, grouped by the same clusters the args are already grouped in (paths, subimage bounds, tile/mem planner inputs, RM range, output-mode flags, mask/cubestat flags, IO-thread/logging/timing flags).
- **Correctness Gate:**
  - Behavior-preserving by construction (same values, same locals downstream) — gate is compile + full `run_tests.sh` (28/28).
  - One bit-identical sweep against T0 archives as cheap insurance, not the primary gate.
  - Zero new compiler warnings across all 4 binary variants (unused-argument/typo risk in the unpack block is the one new failure class here).
- **Rollback Criteria:** Any output byte differs from T0; any of the 28 tests newly fail; any new compiler warning. This ticket is mechanical — a failure means a transcription bug, fix and re-verify rather than reconsider the approach.
- **Effort:** 1 session.
- **Evidence (2026-07-19):** New `rmsynth_config_t` (56 fields, `rm_synthesis_mod.f90`), `read_cfg_keyval` signature reduced to `(cfgfile, cfg, status)` — internal `select case`/duplicate-detection/cross-validation logic untouched, only `var` → `cfg%var` renamed via a quote-aware scripted substitution (154 substitutions; string literals like `'Duplicate key ... path'` and `case ('tile_ra')` selectors correctly left untouched, verified by direct inspection). Call site (`rm_synthesis.f90`) reduced to the single call plus a 56-line mechanical unpack block. Clean 4-variant rebuild: 0 errors, 0 new warnings (same 4 pre-existing GPU linker warnings as T0). `tests/run_tests.sh`: 28/28 pass. Full bit-identical sweep against all 140 T0-archived FITS outputs: 134 exact matches, 6 flagged by `compare_cubes.py --exact` (`badchan_{serial,omp,gpu}.{AMP,PHA}.RMCUBE.FITS`) — investigated directly: all 6 have exactly 201 NaN values in both T0 and new output, at identical locations (the intentionally-NaN fully-masked test pixel), and all 205,623 non-NaN elements bit-for-bit identical — a NaN-vs-NaN comparison artifact in the tool's `--exact` mode, not a regression.

---

### T2 — Tile / VRAM Planner Encapsulation — DONE
- **Objective:** Move the RAM/VRAM budget math into a dedicated routine/derived type — it's pure arithmetic, no I/O, no concurrency, the lowest-risk block in the file to relocate. Also resolves the coupling between cfg and the planner: `tile_ra`/`tile_dec`/`mem_frac_ram`/`mem_frac_vram`/`gpu_vram_mib` are simultaneously cfg-parsed *and* auto-tuner-overwritten today; this ticket makes the planner the single place both happen.
- **Scope:** `rm_synthesis.f90:1717-1963` exactly (through the `dry_run` autotune-cfg writer). Does **not** touch the tile loop itself (`2966+`).
- **Change Set:** `tile_plan_t` derived type — inputs (`nz_out`, `nrm_out`, `nx_out`, `ny_out`, `rem_mean`, `use_input_mask`, `need_icube`, `cubestat`, `io_overlap`, `use_gpu_actual`, `mem_frac_ram`, `mem_frac_vram`, `gpu_vram_mib`, `tile_ra_in`, `tile_dec_in`, `tile_auto`) and outputs (`bytes_per_tile_pixel_ram`, `bytes_per_tile_pixel_ram_out`, `bytes_per_vram_pixel`, `tile_pixels_max`, `tile_bytes_est`, `tile_ra`, `tile_dec`, `gpu_vram_mib_eff`, `ny_sub`, `inflight_slots_planned`, `use_staging`, `mem_frac_vram_per_slot`). Body of the arithmetic moves verbatim.
- **Correctness Gate:**
  - `run_tests.sh` §8 ("Auto tiling shape") is the most directly targeted test; §9 (VRAM staging bit-identical) indirectly depends on `ny_sub`/`use_staging`.
  - **Full bit-identical re-check across all T0 cfgs is required** (not optional) — this block's output determines tile geometry for the entire rest of the run, so a subtle transcription error here (dropped `int64` cast, reordered boundary `if`) could silently change downstream FP-reduction order.
  - Explicitly re-verify `tile_autotune.cfg`/`runtime_estimate.txt` output text is byte-identical (user-facing, easy to reformat by accident while moving the block).
- **Rollback Criteria:** §8 or §9 fails; any `tile_autotune.cfg` diff; any AMP/PHA diff beyond the existing accepted FP-reassociation tolerance (test 6/7's OMP-vs-serial class). Cheap to revert — nothing else depends on this ticket having landed.
- **Effort:** 1 session.
- **Evidence (2026-07-19):** New `tile_plan_t` + `plan_tile` (`rm_synthesis_mod.f90`) — the RAM byte-budget accounting, auto-tiling policy, safety-shrink loop, and VRAM sub-block planning (incl. the `gpu_vram_mib_eff` cfg/env/default resolution) moved verbatim, only `var` → `plan%var` renamed. Line numbers had shifted from the plan's original citation post-T1; re-confirmed current scope by grepping `dry_run`/diagnostic-print markers before editing. The `/proc/meminfo` read (genuine file I/O) and the two `write(*,*)`/`open(...)` diagnostic and dry-run-writer blocks stay in `rm_synthesis.f90` as caller concerns, per the objective's "pure arithmetic, no I/O" framing; `mem_avail_kb` is resolved by the caller and passed in as a `plan` input. Five planner-internal outputs (`bytes_per_tile_pixel_ram[_out]`, `bytes_per_vram_pixel`, `mem_safe_bytes`, `tile_pixels_max`) turned out to be write-only in the caller (never read downstream, confirmed by grep) — their now-dead unpack lines and declarations were removed rather than kept for cosmetic parity, along with the VRAM-resolution locals (`env_vram`, `vram_safe_bytes`, etc.) that moved fully into `plan_tile`. Clean 4-variant rebuild: 0 errors, 0 new warnings (same 4 pre-existing GPU linker warnings as T0). `tests/run_tests.sh`: 28/28 pass, including §8 and §9. Full bit-identical sweep of all 140 T0-archived FITS outputs: same result class as T1 — 134 exact, 6 flagged (`badchan_{serial,omp,gpu}.{AMP,PHA}.RMCUBE.FITS`), all confirmed to be the same pre-existing NaN-vs-NaN `--exact` comparison artifact, not a regression. `tile_autotune.cfg`/`runtime_estimate.txt`: no T0-archived copies existed to diff against (not captured at T0 time), so verified directly instead — stashed the T2 diff, rebuilt the pre-T2 serial binary, re-ran test §8's exact dry-run cfg to regenerate "before" copies, then restored T2 and diffed: both files byte-identical before/after.

---

### T3a — I/O Orchestration Cleanup (read-side + stage bookkeeping, low risk) — OPEN
- **Objective:** Clean up the plumbing around the already-modular read/timer/log calls in the tile loop — reduce repeated `timer_start`/`timer_stop`/`log_tile_bounds` pairs and the read-thread-splitting arithmetic into named helper calls, without changing what CFITSIO calls happen or in what order.
- **Scope:** `rm_synthesis.f90:2960-3180` (progress counters, tile-loop entry logging, parallel-read dispatch, mask build). Explicitly excludes the `io_overlap` ping-pong block (`2977-3014`) and everything from `write_job` assembly onward — that's T3b.
- **Change Set:** Extract byte-count computation and read-thread-split index arithmetic into named module-level helpers called from the tile loop. The `!$omp parallel do` loop body itself is not touched (risk of silently changing firstprivate/shared semantics) — only what surrounds it.
- **Correctness Gate:** Standard `run_tests.sh` pass (tests 5-7, 10, 13 all exercise this path every run); one bit-identical sweep against T0 as cheap insurance; manual grep-diff of a sample log's `tile_read`/`tile_write` `bytes=` field values against T0 (would otherwise silently break swim-lane plotting, a real but non-test-covered regression class).
- **Rollback Criteria:** Any test failure; any `bytes=` field diff in logs.
- **Effort:** 1 session.

---

### T3b — I/O Orchestration Cleanup (`io_overlap` ping-pong + write-job dispatch, elevated caution) — OPEN
- **Objective:** Same spirit as T3a, for the concurrency-adjacent code: the buffer-slot ping-pong and the `write_job(cur_slot)%...` assembly block (confirmed: ~65 lines of mechanical `%` assignment at `rm_synthesis.f90:3884-3928`, directly analogous to what T1 does for cfg — this is real, existing precedent in this codebase for the exact "bundle into a derived type" move).
- **Scope:** `rm_synthesis.f90:2977-3014` (ping-pong) and `3884-3928` (write-job field assignment only).
- **Why elevated caution:** two real production SIGSEGVs happened one function-call away from this exact code (CFITSIO handle aliasing, write-vs-write race — both documented in `docs/ARCHITECTURE.md`'s postmortems). Both fixes are timing-dependent invariants (join-before-reuse, join-before-dispatch) that are trivially easy to reorder while "cleaning up," with no compiler error and — per both postmortems explicitly — no small/fast-test failure to catch it.
- **Change Set:** Extract *only* the data assembly (`write_job(cur_slot)%field = local`) into `call populate_write_job(write_job(cur_slot), <locals>)`. The synchronization statements immediately after (`tile_write_join` calls, `write_pending`/`cur_slot` bookkeeping, the dispatch call itself, lines ~3930 onward) stay **completely untouched, in their exact original order**, in the main program. Encapsulate the data, never the control flow.
- **Also fix while here (found during plan verification, zero-risk comment-only fix):** the comment at `rm_synthesis.f90:3935` ("`io_write_threads_eff is hard-clamped to 1`") is stale — it predates T6 and describes the pre-T6 unsafe design as if still current. The join-before-dispatch rule it's justifying is still correct and still unconditional, but the *reasoning* given is out of date (T6 replaced the ftpsse-handle-sharing mechanism for `io_write_threads>1`; the rule is now kept defensively, not because of the old hard clamp). Update the comment to match what's already documented in `docs/ARCHITECTURE.md`/`docs/PARALLELISM.md`. This can be done as a standalone one-line fix independent of the rest of T3b if you want it landed sooner.
- **Correctness Gate (non-negotiable, not just "run the suite"):**
  - `run_tests.sh` §13 (`require_no_overlapping_tile_writes`) and §14 (io_write_threads bit-identical) are the primary signal specifically because they exist due to past breakage in this exact area — everything else is secondary.
  - Both postmortems explicitly state the in-suite tests are "necessary but not sufficient" (invisible on tiny/fast test data) — if a larger/slower test case is available (more tiles, slower storage, or the historical 4501×4501 case from the postmortem), re-run it end-to-end after this change to re-exercise the actual timing window that caused the original races.
  - Full bit-identical sweep across all T0 cfgs is **required**, not optional, for this ticket.
- **Rollback Criteria:** Any §13/§14 failure; any crash/hang on a larger re-run; `git diff` showing synchronization statements moved relative to each other (stop and re-derive even before running tests — per the postmortems, these bugs don't reliably show up as test failures).
- **Effort:** 1.5-2 sessions (mostly careful reading/verification; the actual diff should be small).

---

### T5 — Deep Config Threading (cfg% everywhere) — final, lowest-priority ticket — OPEN
- **Objective:** Replace all ~237 declaration lines and ~500+ usage sites in `rm_synthesis.f90` with `cfg%tile_ra` etc. throughout, removing the "unpack into locals" step T1 leaves in place.
- **Scope:** Whole file, but done only after T1-T3b are landed and stable (so `cfg%` isn't threaded through code about to be restructured again).
- **Change Set:** Split into sub-tickets by variable cluster, each independently gated by the full suite — do not attempt as one large diff:
  - T5a: `cfg%tile_ra`/`cfg%tile_dec` and related planner inputs (highest usage count: `tile_ra` alone appears at 69 sites, `tile_dec` at 64).
  - T5b: `cfg%io_*` (io_overlap: 22 sites, io_write_threads: 16, io_read_threads similar).
  - T5c: `cfg%cubestat`/output-mode/mask flags (cubestat: 16 sites, rem_mean: 16, mem_frac_vram: 14, plus remaining smaller-fanout keys).
  - Continue clustering remaining keys by usage count as they're encountered.
- **Correctness Gate:** Full `run_tests.sh` + bit-identical sweep against T0, per sub-ticket (not deferred to the end of all of T5).
- **Rollback Criteria:** Any test/output diff on a given sub-ticket — revert just that sub-ticket's diff, the others are independent.
- **Effort:** 3-5 sessions total, spread across sub-tickets, each its own session-scale unit of work.

---

## Explicitly Out of Scope

**GPU staging/async pipeline** (`rm_synthesis.f90:3284-3826`: `dep_h2d`/`dep_kern`/`dep_d2h`, `slot_idx`/`slot_subid`, `use_async_pipeline`). Deferred per the scope decision above. If revisited later: build a structural race-detection test for it first (analogous to §13's `require_no_overlapping_tile_writes`), then use the same T1.S1-T1.SD2 micro-ticket granularity already established in `TODO/TODO.double-buffer-overlap.md` (local-only) — that's the granularity that actually worked last time this exact code was touched.

## Summary

| Ticket | Scope | Risk | Bit-identical sweep | Effort |
|---|---|---|---|---|
| T0 | none (baseline) | none | N/A | 0.5 session |
| T1 | cfg call site + signature | low | cheap insurance | 1 session |
| T2 | tile/VRAM planner | low (pure arithmetic) | **required** | 1 session |
| T3a | read/mask orchestration | low | cheap insurance | 1 session |
| T3b | io_overlap + write-job assembly | elevated | **required** + large-case re-run | 1.5-2 sessions |
| T5 (a/b/c...) | cfg% everywhere | low-per-edit, high-volume | required per sub-ticket | 3-5 sessions |
| GPU staging | — | high | — | **out of scope** |

Recommended order: T0 → T1 → T2 → T3a → T3b → T5 (sub-tickets, whenever convenient after the above are stable). Total for T0-T3b: ~5-6 sessions.

## Verification (every ticket)

1. `bash scratch/make_all.sh` — clean rebuild, 4 binary variants, no new compiler warnings vs. T0 baseline.
2. `bash tests/run_tests.sh` — full 28-test suite, must remain 28/28 pass.
3. Where the ticket's gate calls for it: `python3 tests/compare_cubes.py <T0 archive> <new output> --exact` across all archived T0 cfgs.
4. Track progress by updating each ticket's status (OPEN → IN_PROGRESS → DONE, with evidence) directly in this plan file or a promoted copy under `planning/`, matching the convention already used in `planning/IO_PARALLEL_OPTIMISATION_PLAN.md`.
