# Critical Assessment of rm_synthesis

*Based on code review and debugging session, 2026-07-04.*

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
an optional input mask FITS — are reduced to one `int8` array in one serial
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

### 2 — Duplicate DFT implementations: the root cause of the sign-reversal bug
*(Sign reversal fixed in `5e0f9ea`; structural duplication remains — see priority table)*

There are two complete, independent implementations of the same DFT:

| Kernel | Path | Where |
|---|---|---|
| `tile_extract_gpu` | CPU | `rm_synthesis_mod.f90` |
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
and `tile_extract_gpu` should be deleted.

**Severity:** structural defect that has already caused one silent correctness
regression and will cause another.

---

### 3 — ~~`tile_extract_gpu` is the CPU kernel~~ ✅ Fixed in `022e7e8`, renamed in `HEAD`
*(Module header comment corrected in `022e7e8`; function renamed `tile_extract_cpu` in follow-up commit)*

The kernel named `tile_extract_gpu` is called exclusively on the **CPU path**.
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

### 6 — `wsum` recomputed O(nrm_out) times per pixel in the GPU kernel

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

---

### 7 — Mask build is serial while N−1 cores sit idle

```fortran
do idx_wts = 1, nx_tile*ny_tile*nz_out    ! serial
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

For a 500×500-pixel tile with 236 channels this is 59 million iterations of
non-trivial logic including integer division and three conditional branches.
It runs on one core while all others wait. Adding `!$omp parallel do` here
costs two lines and gives full thread utilisation before the extraction begins.

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

### 10 — ~~Vestigial dead variables in `tile_extract_gpu`~~ ✅ Fixed in `022e7e8`

After the mask-consolidation refactoring, `tile_extract_gpu` declares and
never uses: `pix_base`, `iz`, `kk`, `per_pix_valid`, `mask_val`. These are
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

### 12 — `prepare_gpu_data` doubles the tile memory footprint

It copies `specQ`/`specU` into `specQ_gpu`/`specU_gpu` — a full duplicate of
the tile data. The GPU kernel could read from `specQ` directly using a
computed flat index with no semantic change. The copy is justifiable only if
the memory layout genuinely matters for VRAM cache-line coalescing on real GPU
hardware. On CPU with `OMP_TARGET_OFFLOAD=DISABLED` it is pure overhead and
doubles the effective tile RAM requirement.

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

| Priority | Issue | Impact |
|---|---|---|
| Priority | Issue | Impact | Status |
|---|---|---|---|
| **P0** | `stMaskOut` reads uninitialised memory | Corrupt mask FITS output | ✅ `022e7e8` |
| **P0** | Duplicate DFT kernels — sign reversal | Caused sign-reversal bug; will recur | ⚠️ Sign fixed `5e0f9ea`; duplication open |
| **P1** | `wsum` recomputed per RM in GPU kernel | ~200× wasted work in hot path | 🔲 Open |
| **P1** | Serial mask build | Wastes N−1 cores before every tile | 🔲 Open |
| **P2** | Misleading names and false comments | Caused multi-hour debugging sessions | ✅ `022e7e8` |
| **P2** | Dead variables in `tile_extract_gpu` | Compiler warnings, incomplete refactor | ✅ `022e7e8` |
| **P3** | Fixed-form main program | Maintenance liability | 🔲 Open |
| **P3** | `prepare_gpu_data` memory copy | 2× tile RAM, pure overhead on CPU | 🔲 Open |
| **P4** | No I/O / compute overlap | Performance opportunity | 🔲 Open |
| **P4** | Serial tile loop (global unit numbers) | Scalability ceiling | 🔲 Open |
