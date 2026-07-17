# rmtool Architecture (Master Document)

## Purpose
This is the master architecture document for rmtool.

It records implemented architecture and accepted design decisions in released
versions.

## System Overview
rmtool performs RM synthesis on large FITS spectral cubes using a tiled,
memory-bounded workflow designed for both CPU and GPU execution.

High-level flow:

1. Read a spatial tile from FITS inputs.
2. Build/merge validity masks for that tile.
3. Execute RM extraction on CPU or GPU.
4. Compute optional cubestat maps.
5. Write tile outputs back to FITS products.

## Codebase Map

### Core Runtime
- `src/rm_synthesis.f90`: top-level orchestration, tile loop, IO stages,
  backend selection, and execution control flow.
- `src/rm_synthesis_mod.f90`: shared module utilities, timers/logging,
  configuration parsing, and helper routines.

### Build and Delivery
- `Makefile`: primary development build matrix (OMP/GPU variants).
- `CMakeLists.txt`: distribution-oriented build entry.
- `docker/`: container build and release helpers.

### Configuration and Operations
- `cfg/`: runtime configuration files and schema guidance.
- `scratch/slurm/`: HPC job script examples and operational wrappers.

### Validation and Diagnostics
- `tests/`: regression and behaviour checks.
- `scripts/`: tooling including swim-lane plotting and benchmark helpers.

## Detailed Architecture: Memory, Parallelisation, Offload, and Diagnostics

For focused,
maintained deep-dives, see:
- [PARALLELISM.md](PARALLELISM.md) for memory and execution decomposition diagrams and
  examples.
- [DESIGN_CPU_GPU_TIMELINE_AND_RM_BLOCKING.md](DESIGN_CPU_GPU_TIMELINE_AND_RM_BLOCKING.md) for timeline, RM chunking,
  and swim-lane interpretation details.

### Purpose
This section records architecture-level design choices in rmtool:
- how data is tiled and staged in memory,
- how CPU and GPU execution are parallelised,
- why RM chunking is used,
- and how runtime diagnostics map to those design choices.

The primary driver is scalable processing for very large cubes under
constrained host RAM and device VRAM.

### Design Goals
- Process cubes larger than available host RAM using tiled reads and writes.
- Process cubes larger than available device VRAM using sub-block staging and
  bounded offload windows.
- Keep a single scientific kernel structure that works across CPU and GPU
  builds.
- Preserve numerical correctness while enabling practical throughput on
  desktops and HPC nodes.

### Data Movement and Memory Strategy

#### Host RAM tiling
- Input FITS data is read in tile windows (RA x Dec x full channel span for the
  tile).
- Tile dimensions are selected by user config or auto-tuner and bounded by
  `mem_frac_ram`.
- This limits peak host memory while keeping enough work per tile to amortize
  IO overhead.

#### Device VRAM staging
- For GPU-enabled runs, tiles can be further partitioned into Dec-strip
  sub-blocks sized by `mem_frac_vram` and effective VRAM.
- This staging allows datasets larger than VRAM to run by streaming bounded
  sub-blocks.
- In host-OMP-enabled GPU runs, a two-slot pipeline can overlap
  prep/compute/scatter phases with dependency ordering.

#### Why this is necessary
- Large data products can exceed either host RAM or VRAM if processed
  monolithically.
- Tiling and staging convert memory capacity limits into scheduling and
  streaming decisions.

### Compute Kernel Strategy

#### RM chunking
- The RM synthesis axis is processed in blocks (`nrm_block_size`) rather than
  all RM bins at once.
- This bounds active template footprint and provides controllable work
  granularity.

#### Why RM chunking exists
- GPU: required for bounded offload and kernel launch sizing under VRAM
  constraints.
- CPU: retained intentionally for code-path alignment and for potential cache
  and translation lookaside buffer (TLB) benefits when `nrm` is large.
- Cross-platform consistency: same high-level decomposition helps validation
  and maintenance.

#### CPU-specific note
- CPU path uses RM chunking too; this is algorithmically valid and intentional.
- Benefit is workload-dependent:
  - large `nrm`: can improve locality and scheduling behaviour,
  - small `nrm` (single block): mostly neutral overhead.

### Parallelisation Strategy

#### CPU execution
- OpenMP host parallelism is applied over collapsed loop dimensions in compute
  kernels.
- Thread-level work is data-parallel for (pixel, RM-in-block) combinations.
- Data packing loops in `prepare_cpu_data` and `prepare_gpu_data` are also
  host-parallelised in HOST_OMP builds, with a guard to avoid nested OpenMP
  oversubscription when already inside an active parallel region.

#### GPU execution
- OpenMP target offload is used for device kernels.
- For staging mode, host-side orchestration can pipeline sub-block phases when
  host OMP threads are available.
- In staging mode with HOST_OMP enabled, host gather and scatter loops run with
  OpenMP loop parallelism (`parallel do`/`taskloop`) inside the existing
  dependency-ordered slot pipeline.

#### Synchronisation model
- Dependency-ordered task sequencing ensures slot reuse safety for staged GPU
  operation.
- Non-staging paths remain simpler and synchronous.

### Offload Activation Semantics
- Config key `use_gpu` requests GPU execution.
- GPU-capable binaries can offload when runtime/toolchain/device support is
  present.
- CPU-only binaries warn and fall back to CPU behaviour when `use_gpu` is
  requested.

### IO Architecture Section
- Tile IO is stage-bounded around compute and currently orchestrated in the
  main tile loop.
- IO optimisation is treated as an orchestration concern and must not change RM
  numerical kernels.

### Observability and Diagnostics

#### Structured logs
- Stage logs include `tile_read`, `tile_compute`, `tile_scatter`,
  `tile_write` lifecycle markers.
- GPU async paths emit `tile_async` markers for enqueue/start/done transitions.
- CPU extraction emits `tile_thread` timing markers (thread id, RM block,
  duration).

#### Swim-lane plotting behaviour
- Plotting is an observability layer on top of the runtime logs.
- Mode selection is data-driven:
  - CPU thread detail view for `tile_thread` timing streams,
  - pipeline timeline view for macro-stage overlap.
- For GPU-enabled logs without `tile_async` events, `tile_compute` send/done
  boundaries are used as a synchronous fallback GPU interval.

#### How to read the swim-lane plots

##### Two view modes
- `Process timeline` + `View: Pipeline timeline`
  - Use when GPU/CPU/IO stage overlap is the main question.
  - Lanes are coarse (`GPU`, `CPU`, `I/O`) and bars show stage intervals.
- `Process timeline` + `View: CPU thread detail`
  - Use when thread balance, per-thread gaps, and extraction distribution are
    the main question.
  - Lanes include `CPU stage`, optional `I/O`, and one lane per OpenMP thread
    (`T<tid>`).

##### Legend semantics (pipeline view)
- `I/O read`, `I/O write`
  - FITS or intermediate read/write intervals.
- `CPU prep (odd blocks)`, `CPU prep (even blocks)`
  - Host-side staging/prep intervals for alternating RM blocks.
  - Even blocks use hatch pattern to visually separate alternating block groups.
- `CPU scatter`
  - Host-side scatter/merge of computed data back into tile outputs.
- `GPU compute async slot 1`, `GPU compute async slot 2`
  - True async device compute intervals from `tile_async start/done` markers.
  - Presence of both slots usually indicates double-buffered pipelining.
- `GPU compute (synchronous fallback)`
  - Used when async markers are absent in a GPU-enabled run.
  - Constructed from `tile_compute send -> done` boundaries.
  - Indicates inferred synchronous/proxy GPU timing, not explicit async slot
    timing.

##### Legend semantics (CPU thread detail view)
- `cpu_extract rm_block odd`, `cpu_extract rm_block even`
  - Per-thread extraction intervals from `tile_thread` markers.
  - Even blocks are hatched for parity separation.
  - If a run executes only one RM block (`nrm_out <= nrm_block_size`), only
    the odd/non-hatched series appears; hatched-even traces are not expected.
- Optional stage overlays: `CPU mask`, `CPU prep`, `CPU compute`,
  `CPU cubestat`.
- Optional IO overlays: `I/O read`, `I/O write`.

##### Side-panel metadata semantics
- `Run log file`
  - Basename of the log used for parsing.
- `Run selector`
  - Selected run index mode (`latest`, `first`, or numeric selector).
- `Plot date`
  - Local timestamp when the figure was generated.
- `Total wall time (s)`
  - Duration from first to last parsed event inside the selected run window.
- `Execution context`
  - `GPU run inferred` when startup marker `GPU requested and enabled` is
    present.
  - `CPU only run inferred` otherwise.
- `GPU startup marker`
  - Explicitly reports `found` or `not found` for the GPU-enabled startup
    marker.
- `GPU async s1/s2/sync-fb`
  - Count summary of GPU compute intervals by source:
    - `s1`: async slot 1 intervals,
    - `s2`: async slot 2 intervals,
    - `sync-fb`: synchronous fallback intervals (`send -> done` proxy).
- `Overlap metrics`
  - `GPU-GPU overlap (s)`: overlap among GPU compute intervals.
  - `CPU-GPU overlap (s)`: overlap between CPU stage intervals and GPU compute
    intervals.

##### Interpreting common patterns
- `s1/s2 > 0` and `sync-fb = 0`
  - Async slot markers are present; timeline reflects explicit pipelined GPU
    activity.
- `s1 = s2 = 0` and `sync-fb > 0`
  - No async slot markers were logged; GPU intervals are inferred from
    synchronous boundaries.
- Large gaps in CPU thread view with non-empty `CPU stage` lane
  - Usually indicates non-extract phases dominating that time window.
- High `CPU-GPU overlap (s)`
  - CPU host work and GPU compute are overlapping; expected in effective
    pipelined runs.

### Design Trade-offs
- A unified decomposition (tile + sub-block + RM-chunk) simplifies verification
  across backends but can introduce small overhead in trivial workloads.
- Dynamic scheduling and conservative memory bounds improve robustness across
  machines, sometimes at the cost of peak idealized throughput.

### Practical Implications
- For very large datasets, memory-bound decomposition is essential for
  successful completion.
- For CPU-only small-RM jobs, RM chunking is mostly structural; performance
  gains are not guaranteed.
- For GPU jobs, RM and spatial blocking are core to fitting and streaming work
  through limited VRAM.
- Staging-loop host OpenMP improvements apply only when staging is active,
  i.e. GPU-active runs where `use_staging` is true. CPU-only runs with
  `use_staging=false` do not execute those staging gather/scatter loops.

### Released Design Baseline (2.0)
- Staged gather/scatter host loops are OpenMP-parallel in HOST_OMP builds while
  preserving dependency ordering.
- Spectral pack/copy loops in `prepare_cpu_data` and `prepare_gpu_data` use
  host OpenMP parallelisation with `omp_in_parallel` protection against nested
  oversubscription.
- Tile-memory planning uses separate host RAM and device VRAM accounting:
  - `bytes_per_tile_pixel_ram` for host tile sizing.
  - `bytes_per_vram_pixel` for VRAM sub-block sizing.

