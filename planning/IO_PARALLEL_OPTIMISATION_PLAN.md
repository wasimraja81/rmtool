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

#### T4 postmortem: `io_write_threads>1` crashes on a real Setonix run — ⚠ unresolved
The "non-overlapping pwrite is POSIX-safe" assumption above does not hold
for CFITSIO. `io_write_threads>1` opens N independent `FTOPEN(...,1,...)`
(read-write) handles onto the *same* output file, but CFITSIO's own
`fits_already_open()` (`cfitsio-4.3.1/cfileio.c:1512-1520,1653`) explicitly
aliases repeat read-write opens of an already-open file onto one shared
underlying `FITSfile` buffer/state — the comment there says as much:
*"If the file is opened/reopened with write access, then the file MUST
only be physically opened once."* Read-only opens are exempted from this
aliasing for exactly this reason (*"2 different threads cannot share the
same FITSfile pointer"*) — which is why `io_read_threads` is safe and
`io_write_threads` is not: the N "independent" write handles are not
independent at all, and concurrent unsynchronised `ftpsse` calls on them
from an `!$omp parallel do` corrupt CFITSIO's shared buffer bookkeeping.
Root-caused from a real crash: SIGSEGV inside `memmove`, deep in
unsymbolised `libcfitsio.so` frames, occurring immediately after the first
tile's compute finished and the parallel write section fired — log showed
`io_write_threads=4` at the time. `io_write_threads=1` (default) sidesteps
this entirely (aliases the single pre-existing handle, no parallel region
is entered). **Action for now: do not set `io_write_threads>1`** until a
follow-up ticket implements one of: (a) genuine multi-process writers
(CFITSIO's alias table is per-process, so separate processes each get real
independent buffers), or (b) bypass CFITSIO for the pixel write entirely —
pre-compute the byte offset via `fits_get_hduaddrll`, write big-endian
bytes directly via `pwrite()` on an independently-opened OS file
descriptor, and let CFITSIO only own the header and final checksum.

### T5 - Async Read/Write Overlap ✓ DONE
- Objective: overlap tile N write with tile N+1 read/mask/prep/compute/cubestat.
- Scope: IO orchestration policy only; no compute-kernel changes.
- Cfg key: `io_overlap` (default `n`, pre-existing placeholder from T2-T4 wired
  up by this ticket).

#### Why a pthread instead of an OMP task
An `!$omp task` cannot outlive its enclosing parallel region (region exit is
always an implicit barrier). Making the write's region span the next tile's
read/mask/prep/compute would turn every one of *their* existing
`!$omp parallel do` calls into **nested** parallel regions — and this
codebase never sets `OMP_NESTED`/`omp_set_max_active_levels`, so libgomp
silently collapses nested regions to one thread by default. That would
silently disable `io_read_threads` and the compute kernel's thread count
the moment `io_overlap` was turned on, with no error or warning.
A raw POSIX thread (`pthread_create`/`pthread_join` via `iso_c_binding`,
`tile_write_job_t` in `rm_synthesis_mod.f90`) has no lifetime coupling to
OpenMP's team model, so read/mask/prep/compute keep using their existing
parallel regions completely undisturbed. CFITSIO is already built
`_REENTRANT` (confirmed via the vendored `cfitsio-4.3.1` source), and the
write pthread only ever touches the output AMP/PHA/MASK/NVALID/stat-map
handles while the main thread's read touches the *different* input Q/U
handles — disjoint files, so no risk of the same-file handle-aliasing
hazard that made `io_write_threads>1` unsafe (see T4 note / crash
postmortem below).

#### Double buffering
`p_tile_arr`, `phi_tile_arr`, `mask_tile_arr`, `nvalid_tile_arr`, and (when
`cubestat=y`) the four cubestat maps are declared as `pointer, contiguous`
in `rm_synthesis.f90`, each backed by two physical target arrays
(`..._s0`/`..._s1`). The tile loop ping-pongs the pointers each tile
(`cur_slot = mod(tile_seq,2)`); a tile's write reads whichever slot it was
given, while the *next* tile's mask/prep/compute writes into the other
slot, so there is never a moment where the same memory is being written by
one tile's compute and read by the previous tile's write. `io_overlap=n`
(default) allocates only slot 0 and never repoints the pointers — bit-
identical to the pre-T5 plain-allocatable-array behaviour, zero regression
risk. The auto-tiler's RAM byte-per-pixel estimate doubles just the
output-side term when `io_overlap=y`, so `tile_dec` is planned smaller
automatically under a fixed `mem_frac_ram` — the user still only sets a
memory fraction.

#### Synchronisation
Each buffer slot's pending write is joined (`pthread_join`) only when that
*same* slot is about to be reused, two tiles later — not before every
tile, which would serialise everything and defeat the overlap. Both slots
are joined unconditionally after the tile loop ends, since the last two
tiles' writes never get a "two tiles later" trigger within the loop.

#### Correctness evidence
- All 22 pre-existing tests still pass with `io_overlap=n` (default),
  confirming no behavioural change to the existing path.
- New test (`tests/run_tests.sh` section 13): identical 7-tile run with
  `io_overlap=n` vs `io_overlap=y` (`tile_auto=n`, `tile_dec=5`,
  `cubestat=y`) produces bit-identical AMP/PHA/MASK/NVALID/PEAK/RM_PEAK/
  ANG_PEAK/SNR outputs. Manually also verified at 32 tiles (single-Dec-row
  tiles, stresses many join cycles) and on the serial (non-OMP) binary
  (the pthread mechanism does not depend on `HOST_OMP`).

#### Scope not covered by this ticket
- The GPU two-level (VRAM sub-block staging) tile path is not specially
  handled, but needs none: its own async sub-block pipeline fully closes
  its OpenMP region before the tile-write step is reached, so there is no
  structural conflict, and it reuses the same `p_tile_arr` et al. pointers
  transparently.
- `io_write_threads>1` remains unsafe for the reasons in the T4 postmortem
  below; `io_overlap` does not change that risk (still only one thread
  ever touches a write handle at a time in the shipped default).
  Genuine intra-write parallelism (bypassing CFITSIO with raw `pwrite()`
  at pre-computed byte offsets) is a separate, not-yet-started ticket.

## Answers to kickoff questions
1. Can we read FITS cubes in parallel?
   - Yes — implemented via `io_read_threads`.
   - N independent read-only handles per input file read disjoint channel ranges
     simultaneously into the pre-allocated tile buffer.
   - Channel ranges map to different Lustre OSTs; true parallel IO throughput.
   - `io_read_threads=1` (default) is identical to the prior serial single-call path.
2. Can we write FITS in parallel?
   - `io_write_threads` is implemented but **unsafe above 1** — see the T4
     postmortem above. The "non-overlapping pwrite() is POSIX-safe"
     assumption this ticket shipped with does not hold: CFITSIO aliases
     repeat read-write opens of the same file onto one shared buffer, so
     the N "independent" handles corrupt each other under concurrent
     `ftpsse` calls. Keep `io_write_threads=1` (default) until a follow-up
     ticket lands genuine multi-process or raw-`pwrite()` write parallelism.
   - `io_write_threads=1` (default) is identical to the prior serial write path.
   - 2D map outputs (NVALID, MASK, cubestat) remain serial (negligible cost).
   - Separately, `io_overlap` (T5) overlaps the *whole* write step for tile
     N with tile N+1's read/mask/prep/compute — this is safe today since it
     keeps `io_write_threads` at 1 handle, just moves which thread runs it.
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
| T4 Parallel RM-bin writes | Done, `io_write_threads>1` unsafe (see postmortem) | `885f21e`, `904f75d` |
| T5 Async read/write overlap (`io_overlap`) | Done | — |
| T6 Safe intra-write parallelism (multi-process or raw `pwrite()`) | Not started | — |

Pending validation: benchmark `io_overlap=y` on Setonix against the
serial-read baseline log (`JENNIFER_TOO_FULLIM.run.cpu.setonix.log`:
read ~124s/tile, write ~95-110s/tile) to confirm the wall-time reduction
predicted from that data (~34% at serial read, scaling down as
`io_read_threads` increases — see swim-lane analysis). Also benchmark
`io_read_threads=N` (N = Lustre stripe count) vs baseline `io_read_threads=1`
to quantify read-side wall-time improvement; do not set
`io_write_threads>1` until T6 lands.
