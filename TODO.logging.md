# TODO: rmtool logging and timing instrumentation

Branch: rmtool-logging
Goal: add low-overhead, configurable timing and run diagnostics for CPU and GPU paths.

## Logging and timing conventions

- [x] Use local wall-clock timestamps in ISO-8601 format with timezone offset.
  - Example format: 2026-07-12T14:03:09+10:00
  - Decision: yes, ISO-8601 is a good default (human-readable and machine-parseable).
- [x] Add config keys:
  - log_level = error|warn|info|debug
  - timing_enabled = y|n
  - timing_tile_enabled = y|n
  - timing_io_enabled = y|n
  - timing_output_file = optional path (stdout if empty)
- [x] Keep instrumentation overhead small when disabled (single boolean guard at callsites).

## Phase 1: Minimal logger and timer primitives

- [x] Add logging/timing helper routines to src/rm_synthesis_mod.f90.
  - Suggested new utilities:
    - log_message(level, stage, message)
    - timer_start(stage_id)
    - timer_stop(stage_id)
    - timer_add(stage_id, dt)
    - timer_report_summary()
  - Use wall time for durations (prefer omp_get_wtime when available, fallback to system_clock).
- [x] Add stage identifiers for major sections (fixed set, stable names).
- [x] Ensure logger supports both stdout and optional file output.

Acceptance criteria:
- [x] Build succeeds for OMP=0/1 and GPU=0/1.
- [x] timing_enabled=n adds no meaningful runtime overhead.

## Phase 2: Macroscopic timing in main pipeline

### Insertion points identified in src/rm_synthesis.f

- [x] Config and startup section (around lines 330-430)
  - Time config file locate + parse (call read_cfg_keyval)
  - Log selected mode summary (use_gpu requested/actual, tile settings, dry_run)
- [x] Input/output FITS open/init section (around lines 939-1200)
  - Time FTOPEN for Q/U/I/mask
  - Time output file existence checks and ftinit creation
- [x] Header propagation section (around lines 1900-2460)
  - Time ftphpr and header keyword propagation block
- [x] Main tile loop section (around lines 2503-2900)
  - Per-tile wall time total
  - Separate timers per tile for:
    - tile read (FTGSVE calls)
    - mask build
    - data layout prep (prepare_gpu_data / prepare_cpu_data)
    - compute kernel(s) (tile_extract_gpu_rm_blocked loops)
    - cubestat_tail_quantile_maps
    - tile writes (ftpsse/ftpssb/ftpssi)
- [x] Final close/deallocation section (around lines 2900-3082)
  - Time FTCLOS blocks and final cleanup
  - Print one summary table at end

Acceptance criteria:
- [x] End-of-run summary includes total time and stage breakdown percentages.
- [x] Tile timing can be turned off separately (timing_tile_enabled).

## Phase 3: Disk I/O macro counters

- [x] Capture process-level I/O counters at run start and run end.
  - Linux source: /proc/self/io
  - Track at least:
    - read_bytes
    - write_bytes
    - syscr (read syscalls)
    - syscw (write syscalls)
- [x] Report deltas in final summary.
- [x] Add fallback message if /proc/self/io unavailable.

Acceptance criteria:
- [x] CPU and GPU runs report comparable macro I/O deltas for same config/cube.

## Phase 4: GPU offload transfer timing

- [ ] Add timing around host-side prep and offload compute calls in tile loop.
- [ ] Split GPU-path timings into:
  - H2D/data prep bucket (prepare_gpu_data + mapping overhead)
  - device compute bucket (tile_extract_gpu_rm_blocked loop)
  - D2H/scatter bucket (scatter back to p_tile_arr/phi_tile_arr and nvalid)
- [ ] Add counters for number of RM blocks and sub-blocks processed.

Notes on current code hotspots:
- Single-level path: around lines 2588-2675
- Two-level staged path: around lines 2676-2798

Acceptance criteria:
- [ ] Summary clearly shows whether GPU path is transfer-bound or compute-bound.

## Phase 5: Compute-time accounting (DFT + stats)

- [ ] Report aggregate compute time for RM extraction kernels.
- [ ] Report cubestat compute time separately.
- [ ] Report percentage split:
  - read I/O
  - compute RM
  - compute cubestat
  - output write I/O
  - other overhead

Acceptance criteria:
- [ ] Percentages sum to ~100% (small residual allowed for instrumentation overhead).

## Phase 6: Output format and usability

- [ ] Emit a concise summary block for humans.
- [ ] Emit an optional CSV line/file for benchmarking automation.
  - Columns: run_id, mode, cube dims, tile dims, stage times, io counters
- [ ] Add one line identifying binary flavor (cpu_serial, cpu_omp, gpu_offload, gpu_offload_hostomp).

## Phase 7: Tests and validation

- [ ] Add a small test config with timing enabled and tiny subimage.
- [ ] Extend tests/run_tests.sh to verify:
  - run completes with timing enabled
  - summary markers appear in output
  - no regressions in existing outputs
- [ ] Compare CPU vs GPU timing reports on same synthetic dataset.

## Nice-to-have items (likely missed otherwise)

- [ ] Add a run UUID and include it in every log line.
- [ ] Add log throttling for per-tile debug lines on very large cubes.
- [ ] Add compile-time switch to fully compile out debug logging if needed.
- [ ] Record page-fault counters from /proc/self/stat at start/end (minor/major deltas).

## Definition of done

- [ ] Logging/timing is config-driven, default-safe (off or info-level low-noise).
- [ ] CPU and GPU runs both produce reliable stage breakdowns.
- [ ] Data transfer vs compute bottleneck is directly visible in final report.
- [ ] Documentation updated in README or BUILD with usage examples.
