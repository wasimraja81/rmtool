# Critical Assessment of rm_synthesis

*Based on code review and debugging sessions, 2026-07-04/05. Updated through commit `a86138a`.*

---

## What the code does well

### Memory architecture is sound
The two-level tiling (RAM tile → VRAM Dec-strip) is the right design for
large-image RM synthesis. Tile size is RAM-budget-driven (`mem_frac_ram`), the
auto-tuner writes `tile_autotune.cfg` with a suggested subimage region, and
the staging path correctly bounds VRAM usage. For production runs on 25k×25k
images this architecture is what keeps the code viable.

### Pre-computed templates
`cos_arr`/`sin_arr` are computed once before the tile loop and held resident
for all tiles. No trigonometry is repeated per pixel or per RM block. This is
the single biggest algorithmic win in the hot path.

### Unified mask in a single pass
Three separate masking sources — global bad channels, NaN/Inf in Q or U, and
an optional input mask FITS — are reduced to one `int8` array in one parallel
loop per tile. Clean and cache-friendly compared to the previous three nested
loops.

### Test suite
Ten tests covering serial, OMP, GPU, staging, and auto-tiling paths with
synthetic ground truth, bit-identical CPU checks, and GPU tolerance checks.
This is what caught every regression during the refactoring sessions.

---

## What is bad

### 1 — ~~Latent data-corruption bug: `stMaskOut` reads uninitialised memory~~ ✅ Fixed in `022e7e8`

`stMaskOut` is allocated, **never written**, and then read in the staging
scatter:

```fortran
! src/rm_synthesis.f, staging scatter loop
mask_tile_arr(dst_idx) = stMaskOut(src_idx)   ! stMaskOut is uninitialised
```

The correct source is `stMask_tile_arr`, which was gathered from
`mask_tile_arr` before the kernel call. As written, every run with staging
active and `write_mask_output = y` writes a garbage mask FITS cube. The test
suite never exercises that combination, so this has been invisible.

**Severity:** data-corruption bug in a code path that exists specifically to
support large images — the very runs most likely to use staging.

---

### 2 — ~~Duplicate DFT implementations: the root cause of the sign-reversal bug~~ ✅ Fully fixed in `36eb833`
*(Sign reversal fixed in `5e0f9ea`; `tile_extract_cpu` deleted in `36eb833`; CPU path now calls `prepare_gpu_data` + `tile_extract_gpu_rm_blocked`)*

There are two complete, independent implementations of the same DFT:

| Kernel | Path | Where |
|---|---|---|
| `tile_extract_cpu` (was `tile_extract_gpu`) | CPU | `rm_synthesis_mod.f90` |
| `tile_extract_gpu_rm_blocked` | GPU (and CPU fallback) | `rm_synthesis_mod.f90` |

They must track each other's invariants — channel ordering, L_sq convention,
normalisation, masking logic — but nothing enforces that. The sign-reversal
bug fixed in this session was caused by exactly this duplication: the L_sq
loop direction changed in `rm_synthesis.f` as part of the mask-consolidation
refactoring, only the GPU kernel's direct-index convention happened to survive
correctly, and the CPU kernel's old reversal `iz = nz_out - cnt2 + 1` silently
produced wrong RM signs for every serial and OMP run.

The fix was three lines. The structural problem remains. The CPU path should
call `prepare_gpu_data` + `tile_extract_gpu_rm_blocked` with
`use_gpu_actual=false` — the `collapse(2)` directive works on CPU threads —
and `tile_extract_cpu` should be deleted.

*Fixed in `36eb833`: `tile_extract_cpu` deleted; CPU path now identical to GPU
path. `use_staging` guard fixed in `a73d639` to prevent staging being entered
for CPU runs. All 10 tests pass.*

**Severity:** structural defect that has already caused one silent correctness
regression and will cause another.

---

### 3 — ~~`tile_extract_gpu` is the CPU kernel~~ ✅ Fixed in `022e7e8`, renamed in `HEAD`
*(Module header comment corrected in `022e7e8`; function renamed `tile_extract_cpu` in follow-up commit)*

The kernel (then named `tile_extract_gpu`, now renamed `tile_extract_cpu` in `505f829`) was called exclusively on the **CPU path**.
The module header compounds the confusion:

```fortran
public :: tile_extract_gpu  ! Legacy: kept for compatibility (now wraps tile_extract_gpu_rm_blocked)
```

It does not wrap anything. It reimplements the DFT independently. The comment
is false. During this debugging session the misleading name sent analysis in
the wrong direction multiple times.

---

### 4 — ~~Comment contradicts code on L_sq ordering~~ ✅ Fixed in `022e7e8`

```fortran
! 1) Build L_sq for ALL channels (good and bad) in ascending lambda_sq order
cnt2 = 0
do i = zpix_beg, zpix_end, incs(freq_axis)   ! ascending frequency = descending L_sq
```

The comment says ascending L_sq. The loop produces descending L_sq (ascending
frequency → L_sq = (c/f)² decreases). A comment that contradicts the code is
worse than no comment: it caused a full debugging session to chase the wrong
hypothesis.

---

### 5 — ~~`sampled_freq.txt` displays mismatched columns~~ ✅ Fixed in `022e7e8`

The diagnostic file is written as:

```fortran
do i = 1, nz_out
   write(78,*) zval(zpix_end - (i-1)*incs(freq_axis)), L_sq(i), flag_arr_out(i)
```

`zval(zpix_end - (i-1)*incs)` and `L_sq(i)` do **not** correspond to the same
physical channel. The frequency column counts down from the highest channel
while L_sq counts up from index 1. The file looks like paired data but is not.
Any user or developer reading it to check channel ordering will draw the wrong
conclusion.

---

### 6 — ~~`wsum` recomputed O(nrm_out) times per pixel in the GPU kernel~~ ✅ Fixed in `d1509ea`

Inside `tile_extract_gpu_rm_blocked`, `wsum` is accumulated inside the
`collapse(2)` loop — once per `(pixel, RM)` pair:

```fortran
do ipix = 1, npix
  do i_rm_local = 1, nrm_block_now          ! <-- wsum recomputed here
    wsum = 0.0
    do iz = 1, nz_out
      wsum = wsum + wts_gpu(ipix, iz)        ! identical for every i_rm_local
    ...
    rc_cor = rc_cor / wsum
```

`wsum` depends only on `ipix`, not on `i_rm_local`. Recomputing it
`nrm_block_now` times per pixel wastes `npix × nrm_out × nz_out` additions
where `npix × nz_out` would suffice. For a 512×512-pixel tile with 201 RM
bins and 236 channels this is ~12 billion redundant FMA operations per tile
just to compute a normalisation constant. `wsum` should be precomputed per
pixel in `prepare_gpu_data` (which already loops over the same data) and
passed as a `wsum_gpu(npix)` array.

*Fixed in `d1509ea`: `wsum_gpu(npix)` precomputed in `prepare_gpu_data` via one
parallel pass over `wts_gpu`; kernel uses `wsum_gpu(ipix)` directly.*

---

### 6b — ~~GPU path never populated `nvalid_tile_arr`~~ ✅ Fixed in `475a74c`

The per-pixel valid-channel count (`nvalid_tile_arr`) is written to a FITS map
when `write_nvalid_output = y`. `tile_extract_cpu` sets it during extraction.
`tile_extract_gpu_rm_blocked` never touched it — so any run with `use_gpu = y`
and `write_nvalid_output = y` wrote uninitialised data to the NVALID FITS map.

With `wsum_gpu(npix)` now available after `prepare_gpu_data` returns, the fix
required no extra data passes:

```fortran
do ipix_tile = 1, nx_tile*ny_tile
  nvalid_tile_arr(ipix_tile) = int(wsum_gpu(ipix_tile), kind=2)
end do
```

Applied to both the single-level and staging GPU paths.

---

### 6c — ~~`use_staging` triggered on CPU runs~~ ✅ Fixed in `a73d639`

`use_staging = (ny_sub < tile_dec)` was independent of `use_gpu_actual`.
Staging is a VRAM management mechanism and makes no sense on CPU; the staging
block allocates and calls GPU-path functions that require `use_gpu_actual=true`.
On a system where `gpu_vram_mib` was set to a small value in the config, a CPU
run would enter staging and crash or produce wrong output.

Fixed by making staging conditional on GPU being active:
```fortran
use_staging = (ny_sub.lt.tile_dec) .and. use_gpu_actual
```

---

### 7 — ~~Mask build is serial while N−1 cores sit idle~~ ✅ Fixed in `af85709`

```fortran
!$omp parallel do default(shared) private(idx_wts, iz)
do idx_wts = 1, nx_tile*ny_tile*nz_out
  mask_tile_arr(idx_wts) = 1
  iz = (idx_wts - 1) / (nx_tile*ny_tile) + 1
  if (flag_arr_out(iz) == 0) mask_tile_arr(idx_wts) = 0
  if (nan_check_on) then
    if (specQ(idx_wts) /= specQ(idx_wts) ...) mask_tile_arr(idx_wts) = 0
  end if
  if (use_input_mask) then
    if (specMask(idx_wts) <= 0.5) mask_tile_arr(idx_wts) = 0
  end if
end do
```

For a 500×500-pixel tile with 236 channels this was 59 million iterations of
non-trivial logic including integer division and three conditional branches
running on one core. Fixed with a single `!$omp parallel do` directive.

---

### 8 — No I/O / compute overlap

The tile loop is strictly serial: `read → compute → write`. For any
reasonably fast NVMe the FITS read of tile T+1 could be overlapped with the
extraction of tile T using a double-buffer pattern (one I/O thread, one
compute thread). This is a standard HPC technique and would hide the FITS read
latency at no algorithmic cost.

---

### 9 — The outer tile loops cannot be parallelised

`FTGSVE` and `FTPSSE` use hardcoded Fortran unit numbers (21, 22, 41, 42),
which are global mutable state. You cannot open the same file on two units
simultaneously to read different tiles concurrently. No matter how many nodes,
sockets, or cores are available, the tile loop is permanently serial.

Moving to an I/O abstraction that opens and closes per tile (or using
thread-local unit numbers) would unlock outer-level tile parallelism, which
for a 25k×25k image processed in many small tiles would give near-linear
scaling with core count.

---

### 10 — ~~Vestigial dead variables in `tile_extract_cpu`~~ ✅ Fixed in `022e7e8`

After the mask-consolidation refactoring, `tile_extract_cpu` (formerly `tile_extract_gpu`) declared and
never used: `pix_base`, `iz`, `kk`, `per_pix_valid`, `mask_val`. These are
leftovers from the old dense-packing approach. They generate compiler warnings,
inflate the `private(...)` clause of the OpenMP directive with phantom names,
and signal to any reader that the refactoring was not completed.

---

### 11 — Fixed-form Fortran in the main program

`rm_synthesis.f` is fixed-form Fortran 77: 6-character statement labels,
72-column line limit, `-` continuation markers. The module `rm_synthesis_mod.f90`
is free-form F90. Mixing them forces the reader to context-switch between two
syntaxes. The fixed-form constraints prevent using modern Fortran features
(`associate`, `block`, named `do`-loop constructs, `error stop`) in the main
program without rewriting it. The fixed-form code is not wrong, but it is a
maintenance liability that grows over time.

---

### 12 — Data copies in `prepare_cpu_data` / `prepare_gpu_data` double tile memory

Both prepare functions allocate new 2D arrays and copy the flat FITS data:
- CPU binary: `specQ_cpu(nz_out, npix)` — stride-1 layout needed for cache efficiency
- GPU binary: `specQ_gpu(npix, nz_out)` — coalesced layout for warp access

The copy is the necessary cost of the layout fix (`bf76380`). Without it,
the inner DFT loop suffers stride-5.9 MB cache misses (the 4.3× slowdown
we fixed). The 2× tile RAM overhead is the deliberate trade-off.

The alternative — reading FITS data directly in the right layout — would
require either reading one channel at a time (expensive I/O) or accepting
the old cache-unfriendly layout. The copy remains open as a P3 item but
it is a design trade-off, not a bug.

---

### 13 — `use_gpu` / `use_gpu_actual` is two flags where one would do

```fortran
use_gpu_actual = .false.
if (use_gpu) then
#ifdef USE_GPU
  use_gpu_actual = .true.
#endif
end if
```

Whether the binary was compiled with `USE_GPU` is a compile-time constant.
`use_gpu_actual` could simply be a preprocessor macro or a single Fortran
`parameter`. Instead, both flags are declared, set at runtime, logged, and
threaded through every call site. The split exists because `use_gpu` is
user-configurable; the `#ifdef` guard is not. Exposing both as runtime
variables conflates these two orthogonal concerns.

---

## Priority order for fixes

| Priority | Issue | Impact | Status |
|---|---|---|---|
| **P0** | `stMaskOut` reads uninitialised memory | Corrupt mask FITS output | ✅ `022e7e8` |
| **P0** | Duplicate DFT kernels — sign reversal | Caused sign-reversal bug; structural | ✅ Sign `5e0f9ea`; deleted `36eb833` |
| **P0** | GPU path never set `nvalid_tile_arr` | Corrupt NVALID FITS on GPU runs | ✅ `475a74c` |
| **P0** | `use_staging` triggered on CPU runs | Crash/wrong output on CPU + small VRAM cfg | ✅ `a73d639` |
| **P1** | `wsum` recomputed per RM in GPU kernel | ~200× wasted work in hot path | ✅ `d1509ea` |
| **P1** | Serial mask build | Wastes N−1 cores before every tile | ✅ `af85709` |
| **P2** | Misleading names and false comments | Caused multi-hour debugging sessions | ✅ `022e7e8`, `505f829` |
| **P2** | Dead variables in `tile_extract_cpu` | Compiler warnings, incomplete refactor | ✅ `022e7e8` |
| **P2** | Hardcoded speed of light | 0.069% systematic error on all L_sq | ✅ `2ead708` |
| **P3** | Fixed-form main program | Maintenance liability | 🔲 Open |
| **P3** | Data copies in prepare functions | 2× tile RAM; necessary layout trade-off | 🔲 Open (by design) |
| **P4** | No I/O / compute overlap | Performance opportunity | 🔲 Open |
| **P4** | Serial tile loop (global unit numbers) | Scalability ceiling | 🔲 Open |

---

## Efficiency profile (current architecture)

*Assessed at `a86138a`. Figures use representative parameters:
tile 512×512 pixels, nz\_out=236 channels, nrm\_out=201 RM bins, N=6 cores.*

### Hot path per tile

| Phase | Parallelism | Operations | Notes |
|---|---|---|---|
| FITS read | serial | O(npix × nz) float reads | I/O bound |
| Mask build | **parallel** (`af85709`) | O(npix × nz) with int-div + 3 branches | ✅ fixed |
| `prepare_cpu_data`/`prepare_gpu_data` copy | OMP (wsum loop only; copy serial) | O(npix × nz) reshape into target layout | necessary for cache/coalesce |
| DFT (`tile_extract_gpu_rm_blocked`) | OMP collapse(2) over (npix × nrm\_block) | O(npix × nrm × nz) FMA | peak compute |
| FITS write | serial | O(npix × nrm) float writes | I/O bound |

### DFT arithmetic intensity

For the test cube (512×512 × 236ch × 201RM):
- Work: npix × nrm × nz × 8 FP ops (4 FMA per channel × rc/rs/ic/is)  
  = 262,144 × 201 × 236 × 8 ≈ **100 GFLOPs per tile**
- Read traffic: specQ\_gpu + specU\_gpu = 2 × 262,144 × 236 × 4B ≈ **494 MB**
- Template traffic: cos\_arr + sin\_arr = 2 × 236 × 201 × 4B ≈ **380 KB** (fits in L2/L3)
- Arithmetic intensity: ~100 GFLOPs / 494 MB ≈ **~200 FLOP/byte** → compute-bound on CPU

### Estimated serial bottleneck breakdown (before profiling)

```
Mask build:         O(npix × nz)   parallel  ~  59 M iterations    ✅ `af85709`
prepare_gpu_data:   O(npix × nz)   parallel  ~ 247 MB read+write
DFT kernel:         O(npix×nrm×nz) collapse  ~  49 B FMA          ← dominant
FITS read:          O(npix × nz)   serial    ~ 247 MB disk read
FITS write:         O(npix × nrm)  serial    ~ 211 MB disk write
```

At 8 cores and ~50 GFLOPs/core peak: DFT budget ≈ **2.5 s/tile** (compute).  
Mask build on 1 core: ~59 M ÷ ~500 M iter/s ≈ **0.1 s** — small but avoidable.  
FITS I/O at 500 MB/s: read+write ≈ **~1 s/tile** — overlappable with compute.

### Validated result: compile-time layout selection (`bf76380`)

The `(npix, nz_out)` GPU layout gave stride = 5.9 MB per channel step in the
inner DFT iz-loop on a 4501×4501×288 run (L3 = 20 MB). Every channel access
was a DRAM miss: 288 misses × 101 RM bins = 29,088 DRAM misses per pixel.

Fixed by `prepare_cpu_data` with `(nz_out, npix)` layout (stride-1 channel
access) for CPU binaries; `prepare_gpu_data` keeps `(npix, nz_out)` for GPU.

**CASA fullim benchmark (6 cores, data in page cache):**

| | Regressed CPU path | Fixed CPU path |
|---|---|---|
| Wall time | ~5 min | **25 s** |
| CPU-seconds | ~1200 s | **108 s** |
| Minor page faults | millions/tile | 730 K (first-touch only) |
| Peak RSS | 9+ GB | 2.93 GB |
| CPU utilisation | 544% / 600% | 440% / 600% |

12× wall-time speedup on CASA fullim (data in page cache). The remaining
unused CPU capacity is CFITSIO serial processing overhead (byte-swapping
46 GB through a serial API), confirmed on Jennifer fullim — see below.

**Jennifer fullim benchmark (6 cores, 4501×4501×288, data from disk):**

| | Regressed CPU path | Fixed CPU path |
|---|---|---|
| Wall time | 1:13:52 | **13:59** |
| CPU% | 544% | 204% |
| Bottleneck | DRAM cache misses | Disk I/O (142 MB/s) |

**Jennifer fullim (back-to-back, data in page cache):**
- CPU%: 211% — Amdahl serial fraction = **37%**
- 37% × 839s = **309s serial** = CFITSIO byte-swap overhead processing 47 GB
- Even from page cache, CFITSIO processes each float serially (FITS big-endian → x86 little-endian)
- True compute ceiling with current CFITSIO API: **~270% CPU** on 6 cores

### Remaining open performance work (in priority order)

1. ~~**Parallelise mask build**~~ ✅ Done (`af85709`)
2. **Profile DFT inner loop** — use `perf stat` / `gprof` to confirm cache miss rate and FP utilisation. The `collapse(2)` schedule with default `static` may leave cores idle on the last partial block.
3. **Overlap FITS I/O** — double-buffer: start reading tile T+1 while processing tile T.
4. **Eliminate `prepare_cpu_data` copy** — thread the flat index formula directly into the kernel; saves 2× tile RAM and one memory pass.
