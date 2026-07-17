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

### IO Architecture

#### Tile IO stage structure

Each tile passes through five sequential stages before the next tile begins:

```
tile_read → tile_mask → tile_prep → tile_compute → tile_scatter → tile_write
```

IO is bounded to the `tile_read` and `tile_write` stages. Compute kernels are
not affected by any IO change.

#### FITS disk layout and access pattern

The output cube is stored RA-fastest on disk (FITS NAXIS ordering):

```
chan=1 / rm=1:  [ dec=1, ra=1..nx ] [ dec=2, ra=1..nx ] … [ dec=ny, ra=1..nx ]
chan=2 / rm=2:  …
…
```

For a tile read of `[ra=1..nx, dec=d1..d2, chan=1..nz]`, each channel plane
requires one contiguous read followed by a seek over the remaining Dec rows:

| per tile, per channel | size |
|---|---|
| Contiguous read (Dec strip) | `tile_dec × tile_ra × 4 B` |
| Seek gap to next channel | `(ny − tile_dec) × tile_ra × 4 B` |

With a typical Setonix cube (13 308 × 734 Dec strip × 288 channels), the
read efficiency is `734/11870 ≈ 6%` — the remaining 94% is seek gaps.  This
is intrinsic to the FITS layout vs the RM synthesis access pattern (all
channels needed simultaneously per pixel); it cannot be eliminated without
re-ordering the cube on disk.

#### Parallel read — `io_read_threads`

Channel-chunked parallel reads using `N` independent CFITSIO handles per
input file.  Each thread reads a disjoint channel range of the same spatial
tile into a contiguous slice of the pre-allocated buffer.

```
Serial (io_read_threads=1):
  thread 0:  seek→read chan 1   seek→read chan 2   …   seek→read chan 288
  wall time ∝  288 × (t_seek + t_read)

Parallel (io_read_threads=N):
  thread 0:  seek→read chan   1..288/N  ─┐
  thread 1:  seek→read chan 288/N+1..   ├─ simultaneous
  …                                      ─┘
  wall time ∝  (288/N) × (t_seek + t_read)
```

Channel ranges map to widely-separated byte regions in the file, landing on
different Lustre OSTs and enabling true parallel filesystem throughput.

| Property | Value |
|---|---|
| Decomposition axis | Frequency (channel) |
| Seeks per thread | `nz_out / N` (same total as serial) |
| Read per seek | `tile_dec × tile_ra × 4 B` (unchanged — no spatial loss) |
| Spatial tile footprint | Unchanged — full `tile_ra × tile_dec` per call |
| Buffer layout | Contiguous per thread; no interleaving |
| Lustre benefit | Channel ranges → different byte offsets → different OSTs |

For `io_read_threads=1` (default): `par_unit_Q(1)=21`, `par_unit_U(1)=22`
etc. alias the existing open handles.  No file re-opens, no overhead, loop
runs once — bit-identical to the prior serial path.

#### Parallel write — `io_write_threads`

RM-bin-chunked parallel writes using `N` independent readwrite CFITSIO
handles per output cube (AMP and PHA).  Each thread writes a disjoint RM-bin
range of the tile output buffer to non-overlapping byte regions of the file.

```
Serial (io_write_threads=1):
  thread 0:  write rm 1..nrm (AMP)   then write rm 1..nrm (PHA)

Parallel (io_write_threads=N):
  thread 0:  write rm   1..nrm/N  (AMP + PHA)  ─┐
  thread 1:  write rm nrm/N+1..   (AMP + PHA)  ├─ simultaneous
  …                                              ─┘
```

| Property | Value |
|---|---|
| Decomposition axis | RM (slowest axis in output cube) |
| Output files parallelised | AMP cube (unit 41) and PHA cube (unit 42) |
| 2D map outputs (NVALID, MASK, cubestat) | Remain serial — negligible cost |
| POSIX safety | Concurrent `pwrite()` to non-overlapping byte regions is guaranteed safe |
| CFITSIO metadata safety | File-level metadata updated only at `FTCLOS`, not during subset puts |

For `io_write_threads=1` (default): `par_wunit_amp(1)=41`,
`par_wunit_pha(1)=42` alias existing handles — zero overhead.

#### Why `io_read_threads` and `io_write_threads` are separate from `OMP_NUM_THREADS`

Compute threads and IO threads are constrained by different resources:

| | Compute (`OMP_NUM_THREADS`) | IO threads (`io_read/write_threads`) |
|---|---|---|
| Bottleneck | CPU cores | Lustre OST stripe count, CFITSIO handle overhead |
| Optimal N on Setonix | All available cores (e.g. 128) | File stripe count (typically 4–16) |
| Cost of over-counting | Idle threads, negligible | Excess handles, metadata pressure, diminishing throughput |
| Typical value | `OMP_NUM_THREADS` environment variable | `lfs getstripe <file>` stripe count |

Setting `io_read_threads=128` on a 128-core node would open 128 CFITSIO
handles per input file.  Beyond the Lustre stripe count there is no additional
throughput and metadata server load increases.  The recommended value is the
Lustre stripe count of the input/output files.

#### Thread ceiling enforcement

| Build | Behaviour when `io_*_threads > OMP thread pool` |
|---|---|
| `HOST_OMP=1` | Clamped to `omp_get_max_threads()`; warning printed at startup |
| `HOST_OMP=0` (serial) | Forced to 1 with a note; N>1 in serial binary runs handles sequentially (worse than a single call) |

Additional clamp: `io_read_threads_eff ≤ nz_out` and
`io_write_threads_eff ≤ nrm_out` (can't have more threads than data partitions).

#### Runtime behavior table

| `io_read_threads` | `OMP_NUM_THREADS` | Effective IO threads | Extra file handles |
|---|---|---|---|
| 1 (default) | any | 1 | 0 — aliases units 21/22/40/45 |
| 4 | ≥ 4 | 4 | 3 per input file |
| 8 | 4 | 4 (clamped, warning) | 3 per input file |
| 4 | any | 1 (serial binary) | 0 (forced, note printed) |

Same table applies to `io_write_threads` / units 41/42 / `nrm_out`.

#### Design constraints preserved
- Compute kernels are unchanged.
- Tile spatial footprint is set by the RAM planner; parallel IO does not
  alter it.
- All tests pass at `io_read_threads=1`, `io_write_threads=1` (default).
- Switching to serial is always possible by omitting both keys from the cfg.

### FITS Read Correctness for Large Tiles

#### Background
FITS subset reads are performed via the CFITSIO Fortran wrapper `FTGSVE`.
Internally `FTGSVE` computes a total element count `nelem` as the product of
the per-axis extents being read.  In CFITSIO versions prior to 3.47 this
product was a 32-bit signed integer; in CFITSIO 4.x it is `LONGLONG` (64-bit)
throughout the call chain (`ffgsve` → `ffgr4b` → `ffgbyt`).

#### Bug
When the RAM planner selects a large tile (e.g. 13 308 × 734 × 288 channels =
2.8 × 10⁹ elements), two independent 32-bit integer overflows occur:

1. `allocate(specQ(tile_ra*tile_dec*nz_out))` — size computed in default-integer
   (32-bit) arithmetic wraps silently, producing an undersized buffer.
2. Older CFITSIO (< 3.47) computes `nelem` as a 32-bit int inside `ffgsve`;
   the wrapped negative value reaches `ffgbyt` as an invalid byte count,
   causing a `SIGSEGV` (backtrace: `ffgbyt ← ffgr4b ← ffgsve ← ftgsve_`).

#### Possible fixes
- **Reduce tile size** in the planner — defeats the purpose of RAM-driven tiling.
- **Channel-batch the FTGSVE calls** — keeps the tile intact but splits the read
  along the frequency axis into chunks ≤ INT32_MAX elements.  Correct for any
  CFITSIO version, but adds IO calls and complexity.
- **Fix the integer arithmetic** — correct `allocate()` sizes to use int64 and
  rely on CFITSIO 4.x (which uses `LONGLONG` for `nelem` throughout its internal
  call chain) for single large reads.

#### Adopted fix
The integer arithmetic fix was adopted (`src/rm_synthesis.f90`):

- All `allocate()` size expressions use `int(..., kind=int64)` casts; buffers
  are correctly sized for any tile the RAM planner selects.
- FITS reads remain single direct `FTGSVE` calls; no tile size reduction.

Integer-safety after fix:

| Layer | Safe up to | Mechanism |
|---|---|---|
| Array allocation | 2^63 − 1 elements | `int(...,int64)` in all `allocate()` size expressions |
| CFITSIO `nelem` | 2^63 − 1 elements | CFITSIO 4.x uses `LONGLONG` throughout `ffgsve → ffgr4b → ffgbyt` |
| Coordinate arrays `fpixels` / `lpixels` | 2^31 − 1 per axis | Pixel index values; no current survey axis approaches this limit |

#### Not covered — backward compatibility with CFITSIO < 3.47
The single-call path requires CFITSIO ≥ 3.47.  On an older library, reads of
tiles where `tile_ra × tile_dec × nz_out > INT32_MAX` will crash as above.
The correct mitigation is to update CFITSIO to a current release, not to
reduce tile size or add call-splitting logic.

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

