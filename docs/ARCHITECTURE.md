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

RM-bin-chunked parallel writes for the AMP and PHA output cubes. Each
thread writes a disjoint RM-bin range of the tile output buffer directly
to its byte range on disk via independent Fortran STREAM I/O
(`write_rm_chunk_raw` in `rm_synthesis_mod.f90`), bypassing CFITSIO's
`ftpsse`/handle machinery for the pixel data entirely.

```
Serial (io_write_threads=1):
  thread 0:  ftpsse rm 1..nrm (AMP)   then ftpsse rm 1..nrm (PHA)

Parallel (io_write_threads=N):
  thread 0:  raw-write rm   1..nrm/N  (AMP + PHA)  ─┐
  thread 1:  raw-write rm nrm/N+1..   (AMP + PHA)  ├─ simultaneous
  …                                                  ─┘
```

| Property | Value |
|---|---|
| Decomposition axis | RM (slowest axis in output cube) |
| Output files parallelised | AMP cube and PHA cube |
| Mechanism (N>1) | Independent `newunit=` STREAM units, one `open`/`write`/`close` cycle per thread per RM-chunk, byte position computed from `datastart + pixel_offset*4` |
| 2D map outputs (NVALID, MASK, cubestat) | Remain serial via `ftpsse`/`ftpssb`/`ftpssi` — negligible cost |

For `io_write_threads=1` (default): a single `ftpsse` call per tile through
the one CFITSIO handle already open for the file, unchanged from the
original serial path — bit-identical, zero overhead.

For `io_write_threads>1`: **CFITSIO's handle for the AMP/PHA files
(units 41/42) is closed immediately** after `FTGHAD` fetches the
data-start byte offset (`datastart_amp`/`datastart_pha` in
`rm_synthesis.f90`), *before* any tile write happens — see the postmortem
below for why this early close is load-bearing, not just tidiness.
`write_rm_chunk_raw` computes every subsequent byte position by pure
arithmetic from that one fetched offset:

```
byte_pos = datastart + (irm-1)*ny_out*nx_out*4 + (iy_out_beg-1)*nx_out*4 + 1
```

Each RM-chunk's plane data is copied into a scratch buffer, byte-swapped
to big-endian in place (`host_is_big_endian`/`swap_bytes_r4_inplace` —
checked at runtime, not assumed, so this is also correct on a big-endian
host), and written with a single `write(u, pos=byte_pos) plane_buf` per
RM-plane when the tile spans the cube's full RA width (the recommended
tiling mode — one contiguous run per plane), or one `write` per row when
it doesn't (rare: only reached via an explicit `tile_ra` override or the
single-Dec-row auto-tiler fallback).

> **History: `io_write_threads>1` was unsafe from 2.0 through most of the
> 3.0 cycle, before this raw-write mechanism replaced the CFITSIO-handle
> approach.** The original design opened N independent
> `FTOPEN(...,1,...)` (read-write) handles onto the same file, on the
> (false) assumption that this gave N independent CFITSIO buffers.
> CFITSIO's `fits_already_open()` (`cfitsio-4.3.1/cfileio.c:1512-1520,1653`)
> explicitly aliases repeat read-write opens of an already-open file onto
> one shared `FITSfile` buffer — by CFITSIO's own design ("the file MUST
> only be physically opened once"). Read-only opens are exempted from
> this aliasing for the opposite reason CFITSIO gives: *"2 different
> threads cannot share the same FITSfile pointer"* — which is exactly why
> `io_read_threads` was safe from the start and the original
> `io_write_threads` design was not. Concurrent `ftpsse` calls on the "N"
> write handles from `!$omp parallel do` corrupted CFITSIO's shared
> buffer bookkeeping, confirmed against a real crash: a SIGSEGV inside
> `memmove`, deep in unsymbolised `libcfitsio.so` frames, immediately
> after the first tile's compute finished with `io_write_threads=4` in
> the log. Hard-clamped to 1 as an interim measure until the raw-write
> mechanism above (T6) replaced the whole approach — no CFITSIO handle is
> ever opened more than once for these files now, so there is nothing
> left for `fits_already_open()` to alias.

#### Postmortem: stale CFITSIO buffer flush zeroed raw-written pixel data (found while implementing T6)

The first version of the raw-write mechanism kept the CFITSIO handle for
units 41/42 open until the very end of the program (as the serial path
always had), only using it for header keywords and calling `FTCLOS` at
final cleanup alongside every other output file. This produced silent
data loss, not a crash: every raw write reported `iostat=0` and read back
correctly immediately afterward, `datastart`/`dataend` were confirmed
correct, and the file was the right size — but after the program exited,
every AMP/PHA pixel was `0.0`.

Bisected by reading the file's own bytes from inside the running program,
directly before and after the final `FTCLOS(41,...)` call: the pixel data
was intact right up to that call and zero immediately after it, with
`FTGHAD` re-queried at that same point confirming `datastart`/`dataend`
hadn't shifted (ruling out a header-growth/data-relocation explanation).
Tracing CFITSIO 4.3.1's own source confirmed the mechanism: `ffclos` →
`ffchdu` → `ffpdfl` ("insure correct data fill values") checks/pads the
tail of the data unit against CFITSIO's own `Fptr->filesize`/
`Fptr->logfilesize` bookkeeping — fields updated *only* by CFITSIO's own
write calls. Since the raw stream writer never went through CFITSIO, that
bookkeeping never advanced past "header only" for the whole run. When
`ffpdfl` (via `ffmbyt`/`ffldrc` in `buffers.c`) tried to check the real
end of the pixel data, it saw a position past its own stale `filesize`,
concluded "past EOF", and fabricated an all-zero in-memory buffer marked
dirty — which then got flushed to disk by `ffflsh`, overwriting the real
data with zeros.

**Fix:** close the CFITSIO handle for AMP/PHA right after `FTGHAD`, before
any raw write happens, so CFITSIO's bookkeeping for that file is retired
(and its buffers freed) before there is anything to disagree with. The
final cleanup section skips re-closing units 41/42 when this early close
already happened (`ampha_handles_closed_early` flag in
`rm_synthesis.f90`), avoiding a double-close error.

**Why this wasn't caught immediately:** the bug is invisible to anything
that checks "did the write report an error" or "is the data on disk right
after I wrote it" — both were true. It only shows up in the *final*
on-disk state, after a CFITSIO call that has nothing to do with pixel
data on its face ("close the file"). `tests/run_tests.sh` §14 now
compares full AMP/PHA/MASK/NVALID/PEAK/RM_PEAK/ANG_PEAK/SNR output
bit-for-bit between `io_write_threads=1` and `=4` specifically because a
mid-write `iostat` check would not have caught this class of bug.

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

#### Async read/write overlap — `io_overlap`

Even with `io_read_threads`/`io_write_threads` tuned, the tile loop is
still fully serial *between* tiles: `write(N)` finishes before `read(N+1)`
even starts, though the two touch entirely different files (output
AMP/PHA vs. input Q/U) with no data dependency between them. Real Setonix
timing (`JENNIFER_TOO_FULLIM.run.cpu.setonix.log`, 27 tiles) showed read
(~124s/tile) and write (~95-110s/tile) together at 84% of wall time —
overlapping them was the largest single lever available.

```
Serial (io_overlap=n, default):
  tile N:   read → mask → prep → compute → cubestat → write(N)
  tile N+1:                                                    read → mask → ...
                                                          ^ waits for write(N)

Overlapping (io_overlap=y):
  tile N:   read → mask → prep → compute → cubestat → write(N) ─────────┐
  tile N+1:                                            read → mask → prep → compute → cubestat → write(N+1)
                                                                              ^
                                                             write(N) hidden behind tile N+1's own pipeline
```

**Why a pthread, not an `!$omp task`.** An OpenMP task can never outlive
its enclosing parallel region — the region's exit is always an implicit
barrier. Keeping the write's region open across tile N+1's own
`!$omp parallel do` calls (read chunking, mask build, prep, compute) would
make every one of them a *nested* parallel region, and this codebase never
configures `OMP_NESTED`/`omp_set_max_active_levels` — libgomp silently
collapses nested regions to one thread by default. `io_read_threads` and
the compute kernel's thread count would quietly stop being parallel the
moment `io_overlap` was turned on, with no error. A raw POSIX thread
(`pthread_create`/`pthread_join` via `iso_c_binding`; see
`tile_write_job_t` and `do_tile_write` in `rm_synthesis_mod.f90`) has no
lifetime coupling to OpenMP's team model, so read/mask/prep/compute keep
using their existing parallel regions completely undisturbed. read and
write always touch different files (input Q/U vs. output AMP/PHA), so
they were never at risk of the `io_write_threads>1` handle-aliasing
hazard above — but a *separate* write-vs-write hazard was, and did,
happen; see the postmortem below.

**Thread-pool cost of this design.** The independence from OpenMP's team
model that makes the pthread safe also means it doesn't *share* anything
with the main thread's OpenMP teams — it competes with them for cores
instead. libgomp's worker-thread pool is owned by whichever host thread
first encounters a parallel region, not shared globally; verified
empirically (a background pthread and the main thread each entering
their own 4-thread `omp parallel` region at the same instant produced 8
distinct OS threads running fully concurrently, not 4 shared between
them). Concretely: the write pthread — spawned fresh per tile by
`tile_write_dispatch_async` — is, itself, always +1 OS thread on top of
whatever the main thread is doing. If `io_write_threads>1` as well, that
pthread's own `!$omp parallel do` over RM-chunks spins up a *second*,
separate team of that many OS threads (fresh every tile, since the
pthread itself is fresh every tile), fully additive to the main thread's
concurrent read or compute team. This is usually cheap in practice —
read/write threads are I/O-bound (mostly blocked on the actual
disk/network operation) rather than CPU-bound like compute, so the OS
scheduler gives compute the core time whenever a read/write thread is
blocked — but it is genuine core oversubscription, not something the
runtime silently avoids. See `docs/PARALLELISM.md` ("Thread-pool
interplay") for the full breakdown and a thread-count rule of thumb.

**Double buffering.** `p_tile_arr`, `phi_tile_arr`, `mask_tile_arr`,
`nvalid_tile_arr`, and the four cubestat maps (when `cubestat=y`) are
`pointer, contiguous` in `rm_synthesis.f90`, each backed by two physical
targets (`..._s0`/`..._s1`). The tile loop ping-pongs which target the
pointers reference each tile; a tile's write reads whichever slot it was
given while the next tile's mask/prep/compute writes the other slot, so
compute and write are never touching the same memory at once. With
`io_overlap=n` only slot 0 is ever allocated and the pointers never move —
bit-identical to a plain allocatable array, zero behaviour change. The RAM
auto-tiler doubles just the output-side bytes-per-pixel term when
`io_overlap=y`, so `tile_dec` is planned smaller automatically under the
same `mem_frac_ram` — nothing new for the user to configure.

**Synchronisation.** Two independent rules apply, guarding two different
hazards:
1. *Buffer safety* — a slot's pending write is joined when that same slot
   is about to be reused (two tiles later), so the next tile's
   mask/prep/compute never overwrites memory its predecessor's write is
   still reading.
2. *Handle safety* — before **any** new write is dispatched, whichever
   write is currently outstanding (either slot) is joined first,
   unconditionally, regardless of `io_write_threads`. This predates the
   T6 raw-write mechanism, from when every tile's write shared the same
   two FITS handles regardless of slot and two pthreads calling `ftpsse`
   on the same handle at once was unsafe no matter how well the buffers
   were separated (see postmortem below) — that hazard is what motivated
   the rule. It remains in force unconditionally today, including under
   `io_write_threads>1`'s raw-write path (each tile's own RM-chunks
   already write to genuinely disjoint byte ranges via independent
   stream units, so two *different* tiles' writes could in principle also
   overlap safely by the same POSIX guarantee — but the rule doesn't
   special-case that, since it costs nothing: a tile's write almost
   always finishes before the next tile's own read/mask/prep/compute
   pipeline does anyway, so the join is rarely the thing actually being
   waited on).

Both slots are also joined unconditionally after the tile loop ends,
since the last two tiles' writes never get a "reused two tiles later"
trigger to join them within the loop.

#### Postmortem: tile-write race between consecutive tiles (found on real Setonix-scale data)

The first shipped version of `io_overlap` only implemented synchronisation
rule 1 above (buffer safety) and asserted, incorrectly, that this was
sufficient — see the corrected claim earlier in this section. It is not:
rule 1 only guarantees the write that *previously owned this same slot*
has finished; it says nothing about whatever write is running in the
*other* slot at that moment.

This crashed a real run: a 4501×4501 image tiled into 9 Dec-strips of 561
rows each, except the last, a 13-row leftover (4501 is not a multiple of
561). Tile 8's entire read→mask→prep→compute→cubestat pipeline finished
in ~4s — far faster than a full tile's write (~4-27s) — so tile 8's write
was dispatched while tile 7's write (a *different* buffer slot, so not
caught by rule 1) was still running. Both pthreads ended up calling
`ftpsse` on the same AMP/PHA handles concurrently, producing a SIGSEGV
inside `memmove`, deep in unsymbolised `libcfitsio.so` frames, inside a
thread stack (`start_thread`/`clone3` in the backtrace) — i.e. inside the
write pthread, not the main thread.

**Fix:** added synchronisation rule 2 above (join any outstanding write,
from either slot, before dispatching a new one) in `rm_synthesis.f90`.
Because `pthread_join()` is a hard block rather than a probabilistic race
avoidance, this makes "at most one write in flight" a structural
guarantee rather than a timing-dependent hope — true regardless of how
fast or slow any individual tile happens to be.

**Why the original test suite didn't catch this:** bit-identical output
comparisons validate data correctness, not concurrency safety. A race
that only manifests when one tile's write outlasts the next tile's entire
compute pipeline is invisible on tiny/fast test data, where every stage
finishes in milliseconds regardless. `tests/run_tests.sh` §13 now also
parses the `tile_write` start/done log markers and asserts no two writes'
time windows ever overlap (`require_no_overlapping_tile_writes`) — a
structural check of the actual invariant, not a proxy for it via output
comparison. Re-validated end-to-end on the original crashing case (same
4501×4501 image, `io_read_threads=4`, `io_overlap=y`) — completed without
error.

| Property | Value |
|---|---|
| Mechanism | `pthread_create`/`pthread_join`, independent of the OpenMP runtime |
| Buffers doubled | p/phi/mask/nvalid + cubestat maps (if enabled); not specQ/specU (never touched by write) |
| Default (`io_overlap=n`) | Single buffer, inline write — identical to pre-T5 behaviour |
| Scope | Non-Negotiable Guardrail: read/mask/prep/compute kernels unchanged |
| GPU staging path | Not specially handled, needs none — its own async sub-block pipeline fully closes its OpenMP region before tile-write is reached |
| Platform assumption | `pthread_t` representable as C `long` (glibc/x86_64 Linux) |

Verified bit-identical (`io_overlap=n` vs `y`) across a 7-tile run
(`tests/run_tests.sh` §13) and manually at 32 tiles and on the serial
(non-OMP) binary, plus the no-overlapping-writes structural check
(§13) and an end-to-end production-scale rerun (see postmortem above).

#### Real Setonix production result: write is now the dominant, and largely unhidden, cost

First full-scale production run after the fix: real ASKAP/EMU Q/U cubes
(13308×11870, 288 channels), `io_read_threads=8`, `io_overlap=y`, 16
Dec-strip tiles. Completed end-to-end with no error — the crash is fixed.
The per-stage numbers are worth recording because they reveal a limit of
`io_overlap` that the original (serial-read) profiling data didn't
surface:

| Stage | Total (16 tiles) | % of 2586.7s wall time |
|---|---|---|
| I/O write | 2479.9s | 96% |
| CPU prep | 733.5s | 28% |
| CPU compute | 464.7s | 18% |
| I/O read | 364.1s | 14% |
| CPU cubestat | 28.9s | 1% |
| CPU mask | 12.1s | 0% |

(Percentages sum past 100% because stages genuinely overlap in wall
time — that's the point of `io_overlap`. See the stage-totals bar panel
now included in the swim-lane plots, which renders exactly this table.)

Steady-state write is ~140-150s/tile (one 394s outlier on the very first
tile — a one-off cost from creating/extending the output file on Lustre;
subsequent writes only update already-allocated sub-regions of the
existing cube, hence the much shorter steady state). The rest of the
pipeline (read+mask+prep+compute+cubestat) totals only ~95-115s/tile.
Because `io_read_threads=8` made read fast (~22-30s vs. ~124s in the
original serial-read profile), the "budget" available to hide write
behind is now *smaller* than write itself — so most of write's cost is
exposed rather than hidden. Concretely: the join-before-dispatch fix
(see postmortem above) runs on the main thread, and `read(N+2)` is the
very next statement after that join in program order — so whenever
write(N) outlasts tile N+1's own pipeline, the main thread blocks there,
and read(N+2) waits too. The concurrency window is real but bounded to
one tile's worth of work; it does not compound across multiple tiles.

Despite that, `io_overlap` still delivered a real reduction: a fully
serial equivalent (sum of every stage: 364+12+733+465+29+2480 ≈ 4083s)
would have taken ~68 minutes; the actual overlapped run took 2586.7s
(~43 minutes) — **≈37% wall-time reduction**, in the same range as the
~34% predicted from the original serial-read profile, just achieved
differently (mostly by hiding write behind *compute* now, since read
got fast enough on its own to no longer be the thing doing the hiding).

**Implication that motivated T6:** `io_read_threads` and `io_overlap` had
both done their job — write was overwhelmingly the dominant cost (96% of
wall time) specifically *because* it was stuck at a single serial handle
(`io_write_threads` hard-clamped to 1 at the time this run was measured).
That made T6 (genuine write-throughput parallelism, described above) the
highest-value remaining work rather than a nice-to-have: unlike the
read/overlap work, which was fighting a bottleneck that was already a
*minority* of wall time, T6 was attacking the piece that was, by a wide
margin, most of it. T6 has since landed (see "Parallel write —
`io_write_threads`" above); a production-scale Setonix rerun to measure
the actual wall-time improvement from raw-write parallelism is still
pending — the validation so far is bit-identical correctness on
dev-machine test data, not a production throughput measurement.

#### When `io_overlap=y` can be detrimental

Two independent levers determine whether overlap helps, does nothing, or
actively hurts — they don't move together, so reasoning about "big vs.
small dataset" alone is not enough.

- **RAM size → how much smaller `io_overlap` forces tiles to be.** The
  auto-tiler doubles the output-side byte estimate under `io_overlap=y`
  (see the RAM tile planner above), so `tile_dec` shrinks to compensate
  under a fixed `mem_frac_ram`. Generous RAM: doubling a big allowance
  still leaves a big one, tile count barely moves. Tight RAM: doubling can
  eat a large fraction of an already-small budget, pushing tile count up a
  lot — and more tiles means more fixed per-tile overhead (OMP region
  spin-up, mask-build loop launch, log writes, cubestat call), which is
  pure CPU/software cost independent of disk speed.
- **Disk speed *and architecture* → how much write-time there is to hide,
  and whether "concurrent" means genuinely overlapped or contended.**
  This needs a further split:
  - *Single-spindle disk* (one actuator, classic HDD): reading file A and
    writing file B concurrently forces the head to seek back and forth
    between two locations repeatedly, instead of each operation running
    as one long, mostly-sequential pass serially. Overlap can make this
    **actively worse**, not merely less beneficial.
  - *High-latency but parallel* storage (Lustre, multi-server NFS, cloud
    block storage): concurrent read+write can genuinely use independent
    hardware paths, so the latency really is hidden — this is what the
    Setonix profiling data (and the benefit model in the read/write
    overlap section above) assumed.

|  | Fast disk | Slow disk |
|---|---|---|
| **Small RAM** | Worst case: real cost (more, smaller tiles → more fixed overhead), almost no benefit (fast disk ⇒ little write-time to hide). Likely net *negative*. | Depends on architecture. Parallel FS: probably still a net win, blunted by tile fragmentation and by the OS having less spare page-cache to buffer both streams at once. Single-spindle: likely net *negative* — seek-thrashing stacks on top of the fragmentation cost. |
| **Large RAM** | Neutral: doubling the buffer costs nothing, but there's little write-time to hide either (likely compute- or read-bound already). Harmless to leave on, no real win. | Best case: lots of write-time to hide, RAM absorbs the buffer doubling without meaningfully shrinking tiles. This is the regime the feature targets, closest to the Setonix profile. |

The single sharpest predictor of "detrimental" is **small RAM combined
with either a fast disk or a single-spindle slow disk** — both make the
cost (tile fragmentation) outweigh the benefit (not much I/O time to hide,
or the "hiding" itself backfiring via seek contention). Before enabling
`io_overlap=y` on a RAM-constrained machine, check whether the target
storage is parallel (Lustre/NFS/cloud) or a single physical drive — that
answer flips the recommendation. Given how architecture-dependent this is,
treat the table as a starting hypothesis, not a substitute for a short
head-to-head timing run on the actual target machine (the swim-lane
plotter's separate I/O read/write lanes make this cheap to check).

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
- `tile_read`/`tile_write` start/done lines also carry a `bytes=<N>` field
  -- the actual payload size for that call (read: Q+U plus optional input
  mask/I cube; write: AMP+PHA plus whichever of MASK/NVALID/PEAK/RM_PEAK/
  ANG_PEAK/SNR are enabled), letting the swim-lane plotter compute real
  MB/s per interval instead of only duration.
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
  - Lanes are coarse (`GPU`, `CPU`, `I/O read`, `I/O write` — read and
    write are separate lanes specifically so overlapping reads/writes
    under `io_overlap` render as genuinely concurrent bars instead of
    colliding in one shared lane).
- `Process timeline` + `View: CPU thread detail`
  - Use when thread balance, per-thread gaps, and extraction distribution are
    the main question.
  - Lanes include `CPU stage`, `I/O read`/`I/O write`, and one lane per
    OpenMP thread (`T<tid>`).
- `Stage time totals` (bottom panel, both views)
  - Horizontal bar chart of total wall-clock seconds per stage, largest on
    top, labelled with seconds and % of the run's total wall time.
  - A bar rather than a pie: real runs are often extremely skewed (one
    stage at >90% of wall time with everything else in single digits),
    which a pie renders as one slice and an unreadable sliver soup.
  - Percentages can sum past 100% — that's expected, not a bug: stages
    genuinely overlap in wall time when `io_overlap=y` (write(N)
    concurrent with tile N+1's read/compute), so their individual shares
    of the wall clock aren't mutually exclusive slices of a whole.
  - Side-panel `Thread IDs` were dropped in favour of just `Threads
    active` (a count) — the full ID list added noise without adding
    diagnostic value once thread counts got into the teens.
- `I/O throughput (MB/s)` (middle panel, both views, only when the log's
  `tile_read`/`tile_write` lines carry a `bytes=` field)
  - One flat horizontal segment per I/O interval, at that interval's
    average MB/s (`bytes / duration`) — a per-operation average, not a
    continuously sampled signal, so segments are flat for their duration
    with gaps in between where nothing was in flight.
  - Stacked directly below the swim-lane/thread panel, sharing its time
    axis (unlike the stage-totals bar, which is categorical) — a dip or
    spike here lines up visually with the Gantt bar directly above it.
  - Absent entirely (no empty panel shown) for older logs predating the
    `bytes=` field, or for any run where every I/O interval lacks a byte
    count.

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

### Released Design Baseline (3.0)
Full detail in the "IO Architecture" section above; summary here for the
same at-a-glance baseline record the 2.0 entry provides:
- `io_read_threads`: N independent read-only CFITSIO handles per input
  cube, reading disjoint channel ranges concurrently.
- `io_write_threads`: N independent Fortran STREAM I/O units write
  disjoint RM-bin byte ranges of the AMP/PHA output cubes directly,
  bypassing CFITSIO for pixel data entirely (the original CFITSIO-handle
  design was found unsafe -- see the T4 postmortem -- and is gone, not
  merely clamped).
- `io_overlap`: a tile's write runs on a background pthread concurrent
  with the next tile's read/mask/prep/compute/cubestat, double-buffered,
  serialized against itself via a hard `pthread_join()` barrier.
- All three default to the pre-existing serial behaviour; existing
  configs are unaffected until a user opts in.
- Real Setonix production validation: write dropped from 96% to 6% of
  wall time (~23x reduction) with `io_write_threads=8`; total wall time
  fell ~25% end-to-end on a 13308×11870×288 ASKAP/EMU workload. Full
  before/after breakdown in `docs/RELEASE_NOTES_3.0.md`.
- Swim-lane plotter: separate I/O read/write lanes, stage-totals bar
  panel, I/O throughput (MB/s) panel, complete CPU-stage row (mask/prep/
  compute/cubestat), working GPU synchronous-fallback rendering.
- Test suite grown from 22 to 28, adding structural concurrency-invariant
  checks (no-overlapping-writes for `io_overlap`, bit-identical
  `io_write_threads=1` vs `=4`) alongside the existing output-correctness
  checks.

