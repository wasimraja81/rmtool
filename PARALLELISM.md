# RM-Synthesis Parallelism and Memory Architecture

This document describes how `rm_synthesis` tiles the sky image, loads it into
RAM/VRAM, and distributes work across CPU cores or GPU threads.

---

## 1 — Serial FITS Load into RAM (Tiled I/O)

The full sky image is too large to hold in RAM at once (e.g. 25k×25k×236
channels). The code reads it in **spatial tiles** chosen to fit within a
user-controlled fraction of available RAM (`mem_frac_ram`).

```
DISK  ─────────────────────────────────────────────────────────────────────────
  Q.FITS  [nx_total × ny_total × nz_out]   (e.g. 25600 × 25600 × 236 ch)
  U.FITS  [nx_total × ny_total × nz_out]
            │
            │  FTGSVE  (CFITSIO strided sub-image read)  ← serial, one tile
            ▼
RAM  ──────────────────────────────────────────────────────────────────────────
  specQ  [tile_ra × tile_dec × nz_out]     flat float32 array
  specU  [tile_ra × tile_dec × nz_out]     flat float32 array
  mask_tile_arr [tile_ra × tile_dec × nz_out]  int8, built in single pass:
        ├─ global bad channels  (flag_arr_out)
        ├─ NaN/Inf in Q or U
        └─ input mask FITS       (if provided)
```

The outer loops over tiles are **purely serial** — tiles are processed one at a
time, each written back to the output FITS before the next tile is read:

```
for ix_tile = xpix_beg … xpix_end  step tile_ra          ← serial
  for iy_tile = ypix_beg … ypix_end  step tile_dec        ← serial
    FTGSVE → specQ, specU                                  ← serial disk read
    build mask_tile_arr                                    ← serial single pass
    extract P(RM, pixel)  [CPU or GPU — see §2/§3]
    FTPSSE → output FITS                                   ← serial disk write
```

---

## 2 — CPU Parallelism

The CPU kernel (`tile_extract_gpu`) parallelises over **pixels** using OpenMP.
Each core independently computes the full RM spectrum for its assigned pixels.

```
RAM tile:  specQ(npix, nz_out)   npix = tile_ra × tile_dec
           specU(npix, nz_out)
           mask_tile_arr(npix, nz_out)
           cos_arr(nz_out, nrm_out)   ← read-only, shared by all cores
           sin_arr(nz_out, nrm_out)   ← read-only, shared by all cores

  !$omp parallel do  (over ipix = 1 … npix)
  OpenMP divides npix into N contiguous chunks, one per core:

  ipix:  1 ──── chunk ────► npix/N │ npix/N+1 ──── chunk ────► 2*npix/N │ …
         └──── Core 0 ────┘         └──────── Core 1 ──────┘

  ┌──────────────────────────────────────────────────────────────────┐
  │  Core 0          │ Core 1          │ Core 2    │ … │ Core N-1   │
  │  pix 1…npix/N    │ pix npix/N+1…  │ …         │   │ …npix      │
  │                  │   2*npix/N      │           │   │            │
  │                                                                  │
  │  Each core, for its pixel ipix:                                  │
  │    for i_rm = 1 … nrm_out          (sequential)                 │
  │      for cnt2 = 1 … nz_out         (sequential, SIMD eligible)  │
  │        rc += wt * Q[cnt2] * cos_arr[cnt2, i_rm]                 │
  │        rs += wt * Q[cnt2] * sin_arr[cnt2, i_rm]                 │
  │        ic += wt * U[cnt2] * cos_arr[cnt2, i_rm]                 │
  │        is += wt * U[cnt2] * sin_arr[cnt2, i_rm]                 │
  │      P[ipix, i_rm] = sqrt((rc-is)²+(rs+ic)²) / wsum            │
  │  !$omp end parallel do                                           │
  └──────────────────────────────────────────────────────────────────┘

Output:  p_tile_arr   [npix × nrm_out]  flat float32
         phi_tile_arr [npix × nrm_out]  flat float32
```

**Work partition:**

```mermaid
flowchart LR
    T["RAM Tile\nnpix pixels\nnrm_out RM bins\nnz_out channels"]

    subgraph OMP["OpenMP thread pool"]
        direction TB
        C0["Core 0\npix 0…k"]
        C1["Core 1\npix k+1…2k"]
        C2["Core 2\n…"]
        CN["Core N-1\npix …npix"]
    end

    T -->|"!$omp parallel do\nover pixels"| OMP

    subgraph Work["Each core computes (for its pixels)"]
        direction TB
        W1["outer loop: nrm_out RM bins  (sequential)"]
        W2["inner loop: nz_out channels  (sequential / SIMD)"]
        W3["DFT accumulation + normalise"]
    end

    C0 & C1 & C2 & CN --> Work
    Work --> OUT["p_tile_arr\nφ_tile_arr\n(written to FITS)"]
```

---

## 3 — GPU Parallelism

The GPU kernel (`tile_extract_gpu_rm_blocked`) parallelises over **pixel × RM**
pairs using `collapse(2)`. The channel loop remains sequential per work-item.

### 3a — Single-level (tile fits in VRAM)

```
RAM tile  →  prepare_gpu_data  →  specQ_gpu(npix, nz_out)   float32
                                   specU_gpu(npix, nz_out)   float32
                                   wts_gpu  (npix, nz_out)   float32  (0/1)

Templates (read-only, stays on device across RM blocks):
  cos_arr(nz_out, nrm_out)
  sin_arr(nz_out, nrm_out)

RM-block loop  (CPU, serial):
  for i_rm_block = 1 … nrm_out  step nrm_block_size
    ┌────────────────────────────────────────────────────────────────────┐
    │  !$omp target teams distribute parallel do  collapse(2)           │
    │  [offloaded to GPU; falls back to CPU threads if no GPU]          │
    │                                                                    │
    │  for ipix      = 1 … npix           ┐                             │
    │  for i_rm_loc  = 1 … nrm_block_now  ┘  collapsed → GPU threads   │
    │                                                                    │
    │    Each GPU thread (one per pixel×RM pair):                        │
    │      i_rm_global = i_rm_block + i_rm_loc - 1                      │
    │      for iz = 1 … nz_out          (sequential)                    │
    │        wt = wts_gpu[ipix, iz]                                      │
    │        rc += wt*(Q[ipix,iz]-μQ) * cos_arr[iz, i_rm_global]        │
    │        rs += wt*(Q[ipix,iz]-μQ) * sin_arr[iz, i_rm_global]        │
    │        ic += wt*(U[ipix,iz]-μU) * cos_arr[iz, i_rm_global]        │
    │        is += wt*(U[ipix,iz]-μU) * sin_arr[iz, i_rm_global]        │
    │      P[ipix, i_rm_global] = sqrt(…) / wsum                        │
    │  !$omp end target …                                                │
    └────────────────────────────────────────────────────────────────────┘
```

```mermaid
flowchart TB
    FITS["FITS on Disk\n(nx_total × ny_total × nz_out)"]
    RAM["RAM Tile\nspecQ / specU\n(tile_ra × tile_dec × nz_out)"]
    PREP["prepare_gpu_data\n(CPU, serial)\nreshape + apply mask"]
    VRAM["VRAM\nspecQ_gpu(npix, nz_out)\nspecU_gpu(npix, nz_out)\nwts_gpu (npix, nz_out)\ncos/sin_arr (nz_out, nrm_out)"]

    subgraph RMBLK["RM-block loop  (CPU, serial)"]
        direction TB
        BLK["current block: i_rm_block … i_rm_block+nrm_block_size"]
    end

    subgraph GPU["GPU  —  collapse(2) over (pixel × RM_in_block)"]
        direction TB
        T00["thread (pix=0, rm=0)"]
        T01["thread (pix=0, rm=1)"]
        T10["thread (pix=1, rm=0)"]
        TD["… O(npix × nrm_block)\nthreads total"]
        subgraph INNER["Each thread: sequential over nz_out channels"]
            DFT["DFT accumulation\nrc,rs,ic,is += wt·Q/U·cos/sin"]
        end
    end

    OUT["p_tile_arr / φ_tile_arr\n→ written to FITS"]

    FITS -->|"FTGSVE (serial)"| RAM
    RAM -->|"CPU serial"| PREP
    PREP --> VRAM
    VRAM --> RMBLK
    RMBLK --> GPU
    T00 & T01 & T10 & TD --> INNER
    GPU --> OUT
```

### 3b — Two-level staging (tile too large for VRAM)

When the RAM tile does not fit in VRAM, it is further subdivided into **Dec
strips** (`ny_sub` rows). Each strip is gathered into compact staging buffers,
offloaded, and results scattered back.

```
RAM tile  [tile_ra × tile_dec × nz_out]
  │
  │  for iy_sub_beg = 1 … tile_dec  step ny_sub    ← serial, CPU
  │
  ├──► gather:  stQ/stU/stMask_tile_arr  [tile_ra × ny_sub_now × nz_out]
  │
  ├──► prepare_gpu_data  →  st_Q_gpu / st_U_gpu / st_wts_gpu
  │
  ├──► RM-block loop  →  tile_extract_gpu_rm_blocked  (same as §3a)
  │       GPU parallelism: (tile_ra × ny_sub_now) × nrm_block_size threads
  │
  └──► scatter:  stP / stPhi  back into  p_tile_arr / phi_tile_arr
```

```mermaid
flowchart LR
    subgraph RAM_TILE["RAM Tile\ntile_ra × tile_dec × nz_out"]
        S1["Dec strip 1\ntile_ra × ny_sub"]
        S2["Dec strip 2\ntile_ra × ny_sub"]
        SN["Dec strip N\ntile_ra × ny_sub_last"]
    end

    subgraph VRAM_SUB["VRAM (bounded per strip)"]
        G1["offload strip 1 → GPU"]
        G2["offload strip 2 → GPU"]
        GN["offload strip N → GPU"]
    end

    subgraph OUT_TILE["p_tile_arr / φ_tile_arr"]
        R1["scatter results 1"]
        R2["scatter results 2"]
        RN["scatter results N"]
    end

    S1 -->|"gather+prepare\n(CPU serial)"| G1 -->|"collapse(2) kernel"| R1
    S2 -->|"gather+prepare\n(CPU serial)"| G2 -->|"collapse(2) kernel"| R2
    SN -->|"gather+prepare\n(CPU serial)"| GN -->|"collapse(2) kernel"| RN
```

---

## 4 — Summary: Parallelism Dimensions

| Dimension | CPU path | GPU path |
|---|---|---|
| **Tiles (RA × Dec)** | serial | serial |
| **Pixels within tile** | `!$omp parallel do` — N_cores threads | `collapse(2)` — O(npix × nrm_block) GPU threads |
| **RM bins** | sequential per core | batched per block; collapsed into pixel dimension |
| **Channels (nz_out)** | sequential (SIMD by compiler) | sequential per GPU thread |
| **VRAM staging** | N/A | serial Dec-strip loop when tile > VRAM |

**Key invariant:** `cos_arr` and `sin_arr` are pre-computed once, held
resident in RAM (CPU) or VRAM (GPU), and never recomputed per-pixel or
per-RM-block.
