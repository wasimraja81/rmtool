# rmtool

An HPC package for conducting Faraday Tomography (RM-Synthesis) on radio spectro-polarimetric data. The package is built for all machines - scaling from Low-RAM Desktop PCs to HPC Clusters, with GPU acceleration via OpenMP target offload.

## Motivation

Modern spectro-polarimetric surveys routinely produce data cubes bigger
than the memory available on any one computer — a laptop, a workstation,
or a single machine in a larger computing cluster. That is simply the
normal size of the data now, and rmtool was designed with this in mind
from the outset.

Many existing tools force a choice that made sense when datasets were
smaller: try to load the whole cube into memory and hope it fits, or
process it one piece at a time on a single processor while the rest of
the computer's capacity goes unused. rmtool does neither. It adapts to
whatever computer it is given — using many processors and plenty of
memory when they are available, and working just as reliably on a
single processor with limited memory when they are not:

- **Processes the cube in memory-sized tiles, shaped to match how the
  data actually sits on disk — and never chops it up more than it has
  to.** rmtool is never *forced* to hold the full cube in RAM, but it
  isn't forced to fragment it either: if the cube is small enough, or the
  machine has enough memory, to fit the whole thing within budget, rmtool
  processes it as a single tile and only subdivides further when the
  image genuinely doesn't fit. That budget is a user-set fraction of the
  machine's *total* memory — deliberately not whatever happens to be
  free at that moment — so the tile size for a given cube and config is
  reproducible on the same machine regardless of what else is running on
  it at the time. (On a busy shared node, that also means the budget
  isn't automatically reduced for other jobs' usage, so a large fraction
  is worth setting conservatively there.) The same configuration scales
  from a modest workstation to a large HPC node unchanged. The tile
  shape itself is chosen for read speed too, when tiling is needed: a FITS
  cube stores each frequency channel as one contiguous RA/Dec plane, with
  RA varying fastest, so rmtool tiles as full-width, multi-row strips.
  Every tile read is then one contiguous block on disk, not a scatter of
  small fragments.
- **Reorganizes each tile once in memory for CPU compute, not on disk —
  and skips that step entirely on the GPU path.** The channel-by-channel
  layout that makes reads fast is not the layout the CPU's per-pixel RM
  synthesis wants, which needs each pixel's full frequency spectrum
  contiguous for an efficient inner loop. rmtool performs that
  reorganization — sometimes called a "corner turn" — once per tile, in
  memory, in parallel, rather than paying for it as slow, scattered disk
  access. The GPU path doesn't pay this cost at all: its preferred
  layout already matches the order data is read in, so tile preparation
  there is a masked copy, not a transpose. If an input cube instead
  arrived already stored spectrum-first per pixel, the CPU path's
  in-memory step could be skipped too; that is not how imaging pipelines
  produce FITS cubes today, but it is a concrete, known place with
  further speed to gain if that convention changes.
- **Overlaps a tile's write with the next tile's read and compute,
  instead of running every stage in strict sequence.** While one tile's
  results are being written to disk, the next tile is already being read
  and processed — concurrently, on a background thread, rather than
  waiting for the write to finish first. A tile's write and the next
  tile's read/compute have nothing to do with each other, so there's no
  reason to make one wait on the other — keeping storage and compute
  both busy at once instead of idling in turn.
- **Reads and writes over multiple parallel channels.** A single I/O
  stream rarely saturates the bandwidth available on shared HPC storage
  (Lustre and similar filesystems); rmtool can open several read and
  write channels concurrently to make fuller use of it.
- **Supports GPU acceleration, with more headroom still to exploit.**
  Enabling GPU offload moves the core computation onto the graphics
  card. The same config file works unmodified on the CPU-only build
  too — leave `use_gpu=n`, or leave it `y` anyway and the CPU-only build
  prints a warning and proceeds on CPU rather than failing. (Running a
  GPU-capable build with `use_gpu=y` on a machine with no physical GPU at
  all is untested territory — use the CPU-only build there.) On the
  hardware validated to date, GPU throughput is bounded by host-device
  transfer bandwidth (PCIe) more than by the GPU's own compute capacity —
  a faster PCIe link and a higher-end GPU than our test hardware should
  see a larger gain than we've measured so far. This is an active area
  for further tuning, not a finished ceiling.
- **Most of the above is configuration, not code.** Parallel I/O and
  I/O/compute overlap are opt-in settings in a plain-text config file,
  no recompilation either way. GPU offload is also a config toggle — but
  which of the four build variants to run is still a one-time, per-machine
  build choice; the config format itself doesn't change between them.

The result is one configuration format that scales itself to whatever
machine it's run on, from a single workstation to an HPC
facility, without requiring the user to reason about memory budgets,
concurrency, or hardware-specific tuning.

## Features

- **Config-driven RM synthesis** — Use flexible KEY=VALUE configuration files
- **Multi-band support** — merge frequency channels from several input files
  (different pointings, epochs, or telescope bands) into one RM synthesis run
- **Subimage support** — Extract and process spatial/spectral subsets via config parameters
- **FITS I/O** — Native support for FITS format using CFITSIO
- **GPU offload toggle** — Enable OpenMP target offload with `use_gpu=y` in config
- **Fortran 77/90** — High-performance numerical code
- **Cross-band preprocessing toolchain** — `reproject_cubes` (grid alignment)
  and `convolve_cubes` (common angular resolution) prepare mismatched bands
  for a multi-band run; see [Multi-Band Preprocessing Toolchain](#multi-band-preprocessing-toolchain) below

## Quick Start

### Prerequisites

- **gfortran** (Fortran compiler)
- **CFITSIO** library (`libcfitsio-dev` on Debian/Ubuntu)
- **GNU Make**

### Build

```bash
# Simple Makefile build (recommended)
make

# Run executable
./bin/rm_synthesis cfg/your_config.cfg
```

Common explicit build variants:

| Build command | Binary produced |
|---|---|
| `make OMP=0 GPU=0` | `bin/rm_synthesis_release_cpu_serial` |
| `make OMP=1 GPU=0` | `bin/rm_synthesis_release_cpu_omp` |
| `make OMP=0 GPU=1` | `bin/rm_synthesis_release_gpu_offload` |
| `make OMP=1 GPU=1` | `bin/rm_synthesis_release_gpu_offload_hostomp` |

The build commands are unchanged; only the binary naming is now clearer.

See [QUICKSTART.md](QUICKSTART.md) for detailed build instructions.

## Documentation

- **[QUICKSTART.md](QUICKSTART.md)** — Quick reference and build overview
- **[BUILD.md](BUILD.md)** — Comprehensive build system documentation
- **[cfg/CONFIG_README.md](cfg/CONFIG_README.md)** — Configuration file reference
- **[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md)** — Master architecture document for implemented codebase design
- **[docs/PARALLELISM.md](docs/PARALLELISM.md)** — Parallelism and memory decomposition deep-dive
- **[docs/DESIGN_CPU_GPU_TIMELINE_AND_RM_BLOCKING.md](docs/DESIGN_CPU_GPU_TIMELINE_AND_RM_BLOCKING.md)** — Architecture rationale: tiling, RM chunking, CPU/GPU parallelization, offload strategy
- **[planning/IO_PARALLEL_OPTIMISATION_PLAN.md](planning/IO_PARALLEL_OPTIMISATION_PLAN.md)** — IO optimisation plan: parallel read/write, async overlap, and genuine write-throughput parallelism (T0-T6 all adopted)
- **[planning/ENCAPSULATION_REFACTOR_PLAN.md](planning/ENCAPSULATION_REFACTOR_PLAN.md)** — Encapsulation refactor plan: config/tile-planner/IO-orchestration derived types, ticket-by-ticket (T0-T5, all adopted)
- **[planning/MULTI_BAND_TOMOGRAPHY_PLAN.md](planning/MULTI_BAND_TOMOGRAPHY_PLAN.md)** — Multi-band Faraday tomography plan: config schema, frequency merge, cross-band geometry/resolution matching (`reproject_cubes`/`convolve_cubes`), beam-metadata propagation, ticket-by-ticket (T0-T12, all adopted)
- **[CHANGELOG.md](CHANGELOG.md)** — Release history and key changes by version
- **[docs/RELEASE_NOTES_2.0.md](docs/RELEASE_NOTES_2.0.md)** — Detailed release notes for tag 2.0
- **[docs/RELEASE_NOTES_3.0.md](docs/RELEASE_NOTES_3.0.md)** — Detailed release notes for tag 3.0 (IO-efficiency milestone)
- **[docs/RELEASE_NOTES_4.0.md](docs/RELEASE_NOTES_4.0.md)** — Detailed release notes for tag 4.0 (maintainability/documentation milestone)
- **[docs/RELEASE_NOTES_4.1.md](docs/RELEASE_NOTES_4.1.md)** — Detailed release notes for tag 4.1 (diagnostics milestone)
- **[docs/RELEASE_NOTES_5.0.md](docs/RELEASE_NOTES_5.0.md)** — Detailed release notes for 5.0 (multi-band Faraday tomography milestone; in preparation, not yet tagged)

## Configuration

RM synthesis is controlled via configuration files in the `cfg/` directory: plain
`KEY=VALUE` pairs, one per line, `#` for comments. The parser is strict —
unknown keys, duplicate keys, and unparsable values are all rejected outright
(see [cfg/CONFIG_README.md](cfg/CONFIG_README.md) for the exact rules), so a
config that loads at all is already validated in that sense.

Keys marked **(required)** below must always be present. Keys marked
**(required if ...)** are conditionally required. Everything else is optional
and defaults to the value shown if omitted — an omitted optional key is not an
error, unlike an omitted required one.

**Full annotated example, sectioned by purpose:**

```cfg
# --- Input / output paths (required) ---
path    = /path/to/data/        # directory containing both input cubes
infileQ = Q_cube.fits           # Stokes Q input cube, relative to `path`
infileU = U_cube.fits           # Stokes U input cube, relative to `path`
outfile = my_rm_synthesis       # output basename; cube-type suffixes are
                                 # appended automatically (see "Output format"
                                 # below for exactly which files this produces)

# --- Bad-channel handling (required) ---
remove_badchan = n              # y: drop channels listed in the file below
badchan_file   = unused.txt     # one channel index per line; alias: global_badchan_file
                                 # (the key must be present even when
                                 # remove_badchan=n, but the file is only
                                 # ever opened when remove_badchan=y, so any
                                 # placeholder path is fine -- it doesn't
                                 # need to exist)

# --- Subimage extraction (required: subim; rest optional, default = full cube) ---
subim = n                       # y: only process the pixel/channel ranges below;
                                 # n: process the entire cube (ranges below are ignored)
subim_ra_blc   = 1              # RA first pixel (>= 1)
subim_ra_trc   = 0              # RA last pixel (0 = max)
subim_ra_inc   = 1              # RA step (>= 1)
subim_dec_blc  = 1              # Dec first pixel (>= 1)
subim_dec_trc  = 0              # Dec last pixel (0 = max)
subim_dec_inc  = 1              # Dec step (>= 1)
subim_chan_blc = 0              # channel first pixel (0 = first)
subim_chan_trc = 0              # channel last pixel (0 = max)
subim_chan_inc = 1              # channel step (>= 1)

# --- Q/U processing & bias correction (required) ---
rem_mean       = 0              # 1: subtract the per-pixel mean from Q/U before
                                 # synthesis; 0: use Q/U as-is
remove_qu_bias = n              # y: apply an I-cube-based bias correction to Q/U
                                 # (requires path_I/infileI below); n: skip it
resiQ  = 0.0                    # bias-correction residual/slope terms -- only
slopeQ = 0.0                    # meaningful when remove_qu_bias=y; leave at 0.0
resiU  = 0.0                    # otherwise
slopeU = 0.0

# path_I/infileI (required if remove_qu_bias = y)
path_I  = /path/to/data/
infileI = I_cube.fits

# --- RM synthesis sampling (required: ofac, fac, use_auto_rm_range) ---
use_auto_rm_range = 0           # 0: use beg_rm/end_rm/nrm below (manual grid);
                                 # 1: derive the RM range from the data instead
                                 # (beg_rm/end_rm/nrm become optional overrides)
ofac = 4                        # oversampling factor (>= 1); nrm_out = nrm * ofac
fac  = 3.14159265358979         # pi, for lambda^2 conversion -- leave as-is

# beg_rm/end_rm/nrm (required if use_auto_rm_range = 0)
beg_rm = -500                   # RM grid start (rad/m^2)
end_rm =  500                   # RM grid end (rad/m^2). Alias: max_rm
nrm    =  101                   # number of RM samples before oversampling. Alias: nrm_out

# --- Output format (optional) ---
output_mode   = ap              # ap: write amplitude+phase cubes (default);
                                 # ri: write real+imaginary cubes instead
ap_angle_mode = phase           # (output_mode=ap only) phase: PHA cube is arg(F);
                                 # pol: PHA cube is 0.5*arg(F) (polarization angle)
                                 # ap+phase -> OUTBASE.AMP.RMCUBE.FITS + OUTBASE.PHA.RMCUBE.FITS
                                 # ap+pol   -> OUTBASE.AMP.RMCUBE.FITS + OUTBASE.POLA.RMCUBE.FITS
                                 # ri       -> OUTBASE.REAL.RMCUBE.FITS + OUTBASE.IMAG.RMCUBE.FITS

# --- Masking & optional outputs (optional) ---
mask_cube_file       =          # empty: mask cube (if write_mask_output=y) is
                                 # written to OUTBASE.MASK.CUBE.FITS; non-empty:
                                 # overrides that output path instead
mask_input_cube_file =          # empty: no external input mask; non-empty: path
                                 # to a FITS mask applied on top of NaN/Inf and
                                 # global-bad-channel masking
mask_trust_mode      = safe     # safe: tolerate minor mask/data shape mismatches;
                                 # strict: reject on any mismatch
write_mask_output    = y        # y: write the per-channel MASK cube (default on)
write_nvalid_output  = y        # y: write the per-pixel NVALID map (default on)

# --- Cubestat / peak maps (optional) ---
cubestat = n                    # y: also write PEAK/RM_PEAK/ANG_PEAK/SNR 2D maps
                                 # (per-pixel peak amplitude, its RM, angle at
                                 # peak, and S/N) alongside the RM cubes

# --- GPU (optional) ---
use_gpu      = n                # y: request GPU offload (alias: use_gpus).
                                 # On a gpu_offload* binary this enables offload;
                                 # on a CPU-only binary it prints a warning and
                                 # falls back to CPU (not an error, not silent)
gpu_vram_mib = 0                # 0: auto-detect available VRAM; non-zero: override
                                 # the detected size (MiB) used for sub-block planning
mem_frac_vram = 0.70            # fraction of VRAM to budget per compute sub-block

# --- Tile memory planning & I/O parallelism (optional) ---
# Full explanation, defaults, and tuning guidance in "Tile Memory Planning
# and I/O Parallelism" below -- these are usually left at their defaults
# unless you're on a large cube or a specific HPC filesystem.
tile_auto        = y            # y: auto-size tiles from mem_frac_ram (recommended)
tile_ra          = 0            # manual override (0 = auto); ignored when tile_auto=y
tile_dec         = 0            # manual override (0 = auto); ignored when tile_auto=y
mem_frac_ram     = 0.25         # fraction of total system RAM to budget per tile
io_read_threads  = 1            # N independent read-only FITS handles per input cube
io_write_threads = 1            # N-way parallel AMP/PHA writes (raw stream I/O)
io_overlap       = n            # y: overlap tile N's write with tile N+1's compute

# --- Logging & timing (optional) ---
# Full explanation in "Timing And Benchmark CSV Output" below.
log_level            = info     # error|warn|info|debug
log_output_file      =          # empty: stdout; non-empty: append to this file
timing_enabled       = n        # y: print the Timing summary / Macro breakdown
timing_tile_enabled  = n        # y: include tile-level stage timers
timing_io_enabled    = n        # y: include I/O stage timers
timing_csv_file      =          # empty: no CSV; non-empty: append one row per run

# --- Misc (optional) ---
dry_run = n                     # y: read input cube headers, run the tile
                                 # planner, write tile_autotune.cfg (suggested
                                 # tile_ra/tile_dec/mem_frac_* KEY=VALUE lines
                                 # to copy into a real cfg) and
                                 # runtime_estimate.txt, then exit -- no pixel
                                 # data is read and no output cubes are created
```

For the parser's exact validation rules (which combinations are required,
range checks, duplicate/unknown-key rejection), see
[cfg/CONFIG_README.md](cfg/CONFIG_README.md).

## Timing And Benchmark CSV Output

The runtime logger can emit a human-readable timing summary and an optional CSV
row for automation/benchmark tracking.

Add these keys to your config:

```cfg
# Optional timing controls
log_level = info                  # error|warn|info|debug
timing_enabled = y                # master timing switch
timing_tile_enabled = y           # include tile-level stage timers
timing_io_enabled = y             # include I/O stage timers
log_output_file =                 # empty => stdout, else append to file
timing_csv_file = ./timing.csv    # optional: append one CSV row per run
```

Logging behaviour:
- `log_level` controls structured log lines emitted via `log_message`.
	- `error`: errors only
	- `warn`: warnings + errors
	- `info`: run lifecycle messages (recommended default)
	- `debug`: reserved for future verbose diagnostics
- `log_output_file` controls destination for both structured log lines and
	timing summary blocks.
	- empty: output goes to stdout
	- non-empty: output is appended to one consolidated run log file
	- every emitted line in this consolidated log is ISO-8601 timestamped

The consolidated log file includes ISO-8601 local timestamps on structured log
entries, for example:

```text
2026-07-12T14:03:09+10:00 [info] [startup] rm_synthesis run started
2026-07-12T14:03:09+10:00 [info] [startup] binary_flavor=gpu_offload
...
2026-07-12T14:03:11+10:00 [info] [finalize] rm_synthesis run completed
```

When enabled, the run prints:
- `Run summary:` (binary flavor and GPU requested/active state)
- `Timing summary (seconds):` (stage totals and percentages)
- `Macro timing breakdown:` (read I/O, compute RM, compute cubestat, write I/O, overhead)

The optional CSV output writes a header (once) and one row per run with:
- run id and mode
- cube and tile dimensions
- stage timings
- process-level I/O counters

Example run:

```bash
bin/rm_synthesis_release_cpu_serial cfg/rmsynth.cfg
```

Example GPU run:

```bash
bin/rm_synthesis_release_gpu_offload cfg/rmsynth.cfg
```

## Tile Memory Planning and I/O Parallelism

For cubes too large to fit in RAM, the image is auto-tiled into full-RA
Dec strips sized to a fraction of system RAM, with optional parallel
reads/writes and background-thread write overlap. All of these are opt-in
and default to the pre-existing serial behaviour.

```cfg
# Tile memory planning
tile_auto = y                # y: auto-size tiles from mem_frac_ram (recommended)
tile_ra = 0                  # manual override (0 = auto); ignored when tile_auto=y
tile_dec = 0                 # manual override (0 = auto); ignored when tile_auto=y
mem_frac_ram = 0.30          # fraction of total system RAM to budget per tile
mem_frac_vram = 0.70         # fraction of GPU VRAM to budget per sub-block (GPU only)

# I/O parallelism
io_read_threads = 4          # N independent read-only FITS handles per input
                              # cube, each reading a disjoint channel range.
                              # Safe to increase (try your Lustre stripe count).
io_write_threads = 1         # N-way parallel AMP/PHA writes via raw stream
                              # I/O (bypasses CFITSIO for pixel data). Safe
                              # to increase (try your Lustre stripe count).
io_overlap = n               # y: overlap tile N's write with tile N+1's
                              # read/mask/prep/compute on a background thread.
```

**`io_write_threads>1` bypasses CFITSIO for AMP/PHA pixel writes.**
Setting it higher used to open multiple read-write handles onto the same
output file for parallel RM-chunked writes, but CFITSIO aliases repeat
read-write opens of an already-open file onto one shared internal
buffer — the "independent" handles weren't independent, and concurrent
writes through them corrupted that shared state badly enough to crash
with a segfault on a real run. Rather than clamp this to 1 forever, each
RM-chunk is now written by an independent Fortran STREAM unit directly to
its byte offset (computed once via `FTGHAD`), relying only on the POSIX
guarantee that concurrent writes to disjoint byte ranges of one file are
safe — CFITSIO itself is closed for these two files as soon as that byte
offset is fetched, before any tile write happens, so there's no shared
handle left to corrupt. See `docs/ARCHITECTURE.md` ("Parallel write —
`io_write_threads`") for the full root cause and the design.

**`io_overlap`'s writes are still serialized against each other, by
design — independent of `io_write_threads`.** `io_overlap` guarantees at
most one tile's write is ever in flight — enforced with a blocking
`pthread_join()` before dispatching the next one, not a timing-dependent
assumption. (An earlier version of this feature only checked this per
double-buffer slot, which left a gap: a small/fast tile immediately after
a large/slow one — e.g. the leftover partial tile at the bottom of an
image whose height isn't an exact multiple of the tile size — could
dispatch its write before the previous tile's write had finished,
corrupting CFITSIO's shared handle state at the time. Fixed and
re-validated end-to-end on the exact case that crashed; see the T5
postmortem in `docs/ARCHITECTURE.md` / `planning/IO_PARALLEL_OPTIMISATION_PLAN.md`.)
This doesn't reduce the actual overlap benefit — tile N+1's
read/mask/prep/compute/cubestat already run fully concurrently with
write(N) regardless; only *dispatch* of write(N+1) waits for write(N) to
finish, which write(N) has almost always already done by then. Both
`io_write_threads>1` and `io_overlap=y` together are validated
bit-identical to the fully serial path.

**`io_write_threads` and `io_overlap` are independent keys — use either
alone, or both.** `io_write_threads=N>1` with `io_overlap=n` still gives
each tile's write N-way parallelism; it just means the main thread waits
for that write to finish before starting the next tile's read (no
hide-behind-the-next-tile benefit, but also no doubled RAM buffers). Set
`io_overlap=y` on top only once you also have the RAM and the parallel
storage to benefit from it (see below).

**Thread budget: give compute all your cores; read/write don't need
their own.** `OMP_NUM_THREADS` should be all available cores — it's the
only genuinely CPU-bound thread type here. `io_read_threads` never
competes with compute (they run sequentially on the main thread, same
tile), so it's essentially free to set to your storage's stripe count.
`io_write_threads` *does* run concurrently with the next tile's compute
when `io_overlap=y`, but read/write threads spend most of their time
blocked on the actual disk/network operation rather than burning CPU, so
a modest value (stripe count, typically 4–16) costs little even stacked
on top of a fully-subscribed compute pool — just don't set it close to
your full core count. Full mechanics (which threads share libgomp's pool
and which don't) in `docs/PARALLELISM.md` ("Thread-pool interplay").

**`io_overlap` is not a free win — check your RAM and disk before turning
it on.** It doubles the RAM used by the per-tile output buffers (to let a
background thread write tile N while tile N+1 is computed), and the
benefit depends entirely on how much of your wall time is actually spent
waiting on I/O:

| | Fast disk | Slow disk |
|---|---|---|
| **Small RAM** | Likely **worse**: doubled buffers shrink tiles a lot (more per-tile overhead), fast disk means little write time to hide anyway. | Depends: helps on parallel storage (Lustre/NFS); likely **worse** on a single physical drive (concurrent read+write causes seek-thrashing instead of hiding latency). |
| **Large RAM** | Harmless but pointless: tile count barely changes, but there's little I/O time to hide either. | Best case — this is what the feature targets. |

Rule of thumb: if you're RAM-constrained *and* not on a parallel
filesystem (Lustre, multi-server NFS, cloud block storage), leave
`io_overlap=n`. Full reasoning in `docs/ARCHITECTURE.md` under "When
`io_overlap=y` can be detrimental" — when in doubt, time a short run both
ways on your actual target machine; the swim-lane plotter (below) renders
I/O read and I/O write as separate lanes specifically to make this cheap
to check.

## Recent Performance Enhancements

`3.0`'s IO-efficiency milestone (T0-T6, see `docs/RELEASE_NOTES_3.0.md` and
`planning/IO_PARALLEL_OPTIMISATION_PLAN.md` for the full history):

- **Parallel reads** (`io_read_threads`): N independent read-only CFITSIO
  handles per input cube, each reading a disjoint channel range concurrently.
- **Genuine parallel writes** (`io_write_threads`, T6): N independent
  Fortran STREAM I/O units write disjoint RM-bin byte ranges of the
  AMP/PHA output cubes directly, bypassing CFITSIO for pixel data entirely.
- **Async tile-write overlap** (`io_overlap`): a tile's write runs on a
  background thread concurrent with the next tile's read/compute.
- int64-safe tile indexing throughout, closing an INT32_MAX overflow class
  in the runtime arithmetic that reads/writes tile buffers.

All three cfg keys default to the pre-existing serial behaviour; existing
configs are unaffected until a user opts in.

Real Setonix production validation (13308×11870 pixels, 288 channels,
16 tiles), before vs. after this milestone:

| Stage | Before | After | Δ |
|---|---|---|---|
| I/O write | 2479.9s (96%) | 108.3s (6%) | **−95.6%** |
| Total wall time | 2586.7s | 1945.4s | **−24.9%** |

Write dropped ~23x — from 96% of wall time to 6% — by far the largest
single lever in the tool to date. Full before/after breakdown by stage in
`docs/RELEASE_NOTES_3.0.md`.

Validation summary:
- Build matrix: all four `OMP/GPU` variants compile with zero warnings.
- Tests: `28/28` passing.
- Production: end-to-end validated on real Setonix hardware against real
  ASKAP/EMU data — the exact case that originally crashed (see the T4/T5
  postmortems in `docs/ARCHITECTURE.md`) now completes without error.

## GPU Acceleration

rmtool supports OpenMP target offload builds for GPU execution.

Build and run:

```bash
# Build GPU/offload-capable binary
make GPU=1

# Run with a GPU-enabled config (use_gpu=y in cfg/rmsynth.cfg)
bin/rm_synthesis_release_gpu_offload_hostomp cfg/rmsynth.cfg
```

Useful runtime environment variables:

```bash
OMP_TARGET_OFFLOAD=MANDATORY   # fail if offload cannot run
OMP_DEFAULT_DEVICE=0           # choose target device index
```

## Swim-Lane Plot Generation

Use the swim-lane script to visualize overlap between I/O, CPU staging, and GPU
compute from a consolidated run log.

Generate a plot from a run log:

```bash
python scripts/plot_tile_async_swimlane.py \
	--log scratch/RMSYNTH_OUTPUT.run.log \
	--out scratch/tile_async_swimlane.png \
	--run latest \
	--time-axis absolute
```

Key options:

- `--run latest|first|N` selects which detected run block from the log to plot.
- `--time-axis absolute|relative` chooses wall-clock vs seconds-from-run-start.
- `--out` controls output PNG path.

The script also prints summary metrics (interval count, window seconds,
GPU-GPU overlap, CPU-GPU overlap) to stdout.

Every plot now includes a **stage time totals** bar panel underneath the
timeline: total wall-clock seconds per stage, sorted largest-first, with
seconds and % of total run wall time labelled on each bar. A bar chart
rather than a pie, since real runs are often extremely skewed (one stage
at >90% of wall time) -- a pie would render that as one slice and an
unreadable sliver soup. Percentages can add up to more than 100%; that's
expected when stages overlap in wall time (e.g. `io_overlap=y`), not a
bug. The side panel's `Thread IDs` line (CPU thread detail view) was
dropped in favour of just `Threads active` (a count) -- the full ID list
stopped being useful information once thread counts got into the teens.

Design rationale and diagnostic interpretation notes are documented in
[docs/DESIGN_CPU_GPU_TIMELINE_AND_RM_BLOCKING.md](docs/DESIGN_CPU_GPU_TIMELINE_AND_RM_BLOCKING.md).

Example swim-lane plots, current as of the I/O throughput (MB/s) panel and
the "CPU stage" row's compute segment (both described above):

Pipeline/stage-overlap view (async, double-buffered GPU dispatch --
`gpu_offload_hostomp` binary, staged VRAM sub-blocks):

![Swim-lane GPU async example](docs/images/swimlane_gpu_example.png)

Pipeline/stage-overlap view (synchronous-fallback GPU dispatch -- a tile
that fits in one VRAM sub-block, so there's no double-buffering to overlap):

![Swim-lane pipeline example](docs/images/swimlane_pipeline_example.png)

CPU thread-detail view (`io_read_threads`/`io_write_threads`/`io_overlap`
all active):

![Swim-lane CPU thread example](docs/images/swimlane_cpu_thread_example.png)

### GPU Validation Scope For Swim-Lane Diagnostics

- Tested on: `NVIDIA GeForce RTX 3050 (6 GiB VRAM)` with GNU OpenMP offload
	(`nvptx`) in this repository's current workflow.
- Not yet validated: AMD ROCm offload targets, Intel GPU offload targets, and
	very old NVIDIA GPUs/toolchains where OpenMP target offload support differs.
- If your platform is unvalidated, treat swim-lane diagnostics as experimental
	and confirm with a short controlled run before production execution.

## GPU Runtime Behaviour

- `use_gpu=n` runs host execution.
- `use_gpu=y` requests GPU execution when running a GPU-capable binary (`make GPU=1`).
- If `use_gpu=y` is used with a CPU-only binary, the run prints a warning and falls back to CPU.

Useful runtime env vars for GPU runs:

```bash
OMP_TARGET_OFFLOAD=MANDATORY   # fail if offload cannot run
OMP_DEFAULT_DEVICE=0           # choose target device
```

Build examples:

```bash
# CPU-only OpenMP binary
make GPU=0 OMP=1

# GPU/offload-capable binary
make GPU=1
```

**Output Files:**
- `AP + phase mode`: `OUTBASE.AMP.RMCUBE.FITS` and `OUTBASE.PHA.RMCUBE.FITS`
- `AP + pol mode`: `OUTBASE.AMP.RMCUBE.FITS` and `OUTBASE.POLA.RMCUBE.FITS`
- `RI mode`: `OUTBASE.REAL.RMCUBE.FITS` and `OUTBASE.IMAG.RMCUBE.FITS`
- `Common diagnostics` (on by default; disable with `write_mask_output=n`/
  `write_nvalid_output=n`):
	- `OUTBASE.MASK.CUBE.FITS` (per-channel validity mask cube)
	- `OUTBASE.NVALID.MAP.FITS` (valid-channel count map)
- `When cubestat=y`:
	- `OUTBASE.PEAK.MAP.FITS`
	- `OUTBASE.RM_PEAK.MAP.FITS`
	- `OUTBASE.ANG_PEAK.MAP.FITS`
	- `OUTBASE.SNR.MAP.FITS`

**Beam metadata:** if the input Q cube carries `BMAJ`/`BMIN`/`BPA`, every
output above carries the same values, unchanged. If the input instead has
`CASAMBM=T` (a genuinely per-channel-varying restoring beam, e.g. an
un-convolved CASA multi-beam cube), that propagated scalar is only the
input's own nominal/reference value and means nothing on its own — so the
AMP/PHA cubes and the PEAK/RMPEAK/ANGPEAK/SNR maps (not MASK or NVALID,
which are validity bookkeeping, not flux data) also get `CASAMBM=T` plus
the input's own real per-channel `BEAMS` table attached as an extension,
plus a `HISTORY` note explaining why. Run `convolve_cubes` first (below) if
a single, well-defined resolution is required. In multi-band mode, a
mismatch between bands' own beam metadata is reported as a runtime warning
(not a hard error — RM synthesis itself doesn't depend on beam metadata for
correctness).

## Multi-Band Preprocessing Toolchain

Real bands rarely arrive on the same sky grid or at the same angular
resolution. Two standalone tools (own binaries, own build targets,
independent of the main `rm_synthesis` build graph) close that gap before
a multi-band run:

1. **`reproject_cubes`** — reprojects two or more FITS cubes onto one
   common sky grid, using Starlink AST for WCS handling and `astResampleR`
   for resampling. Three footprint modes:

   ```bash
   make reproject_cubes
   bin/reproject_cubes mode=intersection reffile=ref.fits infiles=a.fits,b.fits
   # mode=union | reference also available; --help for the full option list
   ```

2. **`convolve_cubes`** — convolves all channels, across all input files,
   to one common angular resolution (or an explicit target). Reads
   per-channel beams from a CASA-style `BEAMS` binary table (auto-detected)
   or a portable ASCII/CSV beam log (see `cfg/example_beamLog.txt` and
   `cfg/example_beamLog.csv` for the format — nobody should have to
   reinvent it from scratch). A channel is treated as bad — its output
   plane written as all-NaN, not convolved, and automatically excluded by
   `rm_synthesis`'s own NaN detection later, no extra config needed — if
   it's missing from the beam log entirely, or present with `BMAJ` or
   `BMIN` equal to 0:

   ```bash
   make convolve_cubes
   bin/convolve_cubes infiles=bandA.fits,bandB.fits mem_frac_ram=0.25
   # target beam is auto-derived (smallest beam every good channel of
   # every input can be deconvolved from) unless target_bmaj/target_bmin/
   # target_bpa are given explicitly; --help for the full option list
   ```

Typical order for a genuinely mismatched multi-band dataset:
`reproject_cubes` (align grids) → `convolve_cubes` (match resolution,
across all bands together in one call) → `rm_synthesis` (multi-band RM
synthesis on the now grid- and resolution-matched inputs).

Full design detail, verification evidence, and the underlying computation
modules (`src/gaussft.f90`, `src/commonbeam.f90`) are documented in
[planning/MULTI_BAND_TOMOGRAPHY_PLAN.md](planning/MULTI_BAND_TOMOGRAPHY_PLAN.md)
(tickets T10-T12) and in each source file's own header comment.

## Project Structure

```
rmtool/
├── src/                       Source code (Fortran 77/90)
│   ├── rm_synthesis.f90       Main program (free-form F90); `include`s
│   │                          myfits_info.f90/printerror.f90 below at compile time
│   ├── rm_synthesis_mod.f90   Shared module: config parser, timers/logging, helpers
│   ├── myfits_info.f90, printerror.f90   Free-form F90 helpers, pulled into
│   │                                      rm_synthesis.f90 via `include`
│   ├── reproject_cubes.f90    Standalone: cross-band sky-grid alignment (own binary)
│   ├── gaussft.f90            gaussft_mod: pure elliptical-Gaussian FFT-domain
│   │                          deconvolve/reconvolve computation, no I/O
│   ├── commonbeam.f90         commonbeam_mod: smallest common beam across N PSFs
│   ├── convolve_cubes.f90     Standalone: cross-band resolution matching (own binary),
│   │                          drives gaussft_mod + commonbeam_mod
│   └── legacy/                Older standalone FITS utilities, not part of the build
├── cfg/                        Configuration files, examples, and ARCHIVED/ (63 historical configs);
│                                example_beamLog.txt/.csv for convolve_cubes' ASCII beam format
├── docs/                       Architecture, parallelism, and design deep-dives; release notes
├── planning/                   IO optimisation plan and ticket history
├── scripts/                    Swim-lane plotting and benchmark tooling
├── tests/                      Regression suite (tests/run_tests.sh)
├── TODO/                       Historical development logs and assessments
├── docker/                     Container build/release helpers
├── scratch/                    Ad-hoc run outputs, example logs/plots (gitignored)
├── bin/                        Compiled executables
├── build/                      Build artifacts (Makefile)
├── Makefile                    Primary development build (OMP/GPU variants)
├── build.sh                    Quick build wrapper script
└── BUILD.md, QUICKSTART.md     Build and quick-start documentation
```

## Development

### Branch Structure

- **main** — Stable, production-ready releases
- **develop** — Active development branch

### Release Tags

- Formal release tags use `MAJOR.MINOR` format (for example: `1.0`, `1.1`, `2.0`).
- The first formal release tag is `1.0` on `main`.
- Milestone-style tags can still exist for internal checkpoints, but official releases use the numeric format above.

### Building

```bash
make              # Build release
make MODE=debug   # Build with symbols
make clean        # Clean artifacts
```

## License

See [LICENSE](LICENSE) file for details.

## Contact

For questions or contributions, please open an issue or contact the maintainers.
