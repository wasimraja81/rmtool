# Multi-Band Faraday Tomography — Feasibility & Design Plan

Branch: `multi-band-tomography` (from `develop`)

**Status: planning only. No code changes on this branch yet.**

## 1. Motivation

Today rmtool ingests exactly one Q cube and one U cube
(`cfg%infileQ`/`cfg%infileU`), covering a single, contiguous, uniformly-sampled
frequency band, and RM-synthesises them directly. Multi-band tomography means
ingesting **several Q/U cube pairs from different frequency bands** (e.g. a
low-band and a high-band observation of the same field, possibly from
different instruments/pipelines) and RM-synthesising across the combined
frequency coverage in one run. Wider, multi-band λ² coverage is what actually
buys better RM resolution and reduces sidelobes in the RM spread function —
that is the scientific point of this effort.

## 2. Requirements (as given)

1. Ingest multiple Q cubes and multiple U cubes (one Q/U pair per band).
2. Verify that all cubes share a matching geometry before doing anything else
   — if they don't match, refuse to run rather than produce a silently wrong
   result.
2a. **(added, confirmed with user, 2026-07-20)** The new config schema must
   be **backward compatible**: an existing single-band `cfg` file, unedited,
   must keep working with the new code, with no observable behaviour
   change. See §5 for the resolved design.
3. **Decision needed:** attempt on-the-fly geometry matching (reprojection),
   or require the caller to hand rmtool already-matched cubes?
   - 3a. Geometry matching in the on-the-fly case requires RA/Dec overlap
     detection and WCS reprojection.
   - 3b. Cubes must also match in angular resolution (beam) — today rmtool
     assumes one resolution for the whole run and never checks it.

**Decision recorded (confirmed with user, 2026-07-20):** rmtool will
**require pre-matched geometry** for v1. Reprojection/regridding and
resolution-matching (convolution to a common beam) are the caller's
responsibility, done upstream with existing astronomy tools (e.g.
`reproject`, `montage`, a common-resolution convolution step) before handing
cubes to rmtool. rmtool's job is to **validate** that what it was given
actually matches — pixel grid, WCS, and (new) resolution — and refuse to run
otherwise, exactly as it already does for the single Q/U pair today. This is
a much smaller lift than building a reprojection engine, and it keeps the
tiled I/O and parallelism model completely untouched (see §4). It does not
solve the problem for a user whose bands are natively on different pixel
grids — they must regrid first — but that tradeoff was made explicitly,
not by default.

**Follow-up decision recorded (confirmed with user, 2026-07-20), resolving
3b:** the "matching resolution" check will **not** be implemented as a
`BMAJ`/`BMIN`/`BPA` header comparison — those keywords are frequently absent
or unreliable in multi-frequency FITS products, so a header-based beam check
would be brittle in exactly the case it's meant to protect. Instead, the
only automated geometry gate is **exact RA/Dec pixel-grid equality across
all N cubes** (extending today's Q-vs-U check, §3), with a loud warning and
a hard refusal to run on any mismatch — no silent tolerance, matching
requirement 2. Actual angular-resolution matching between bands remains the
caller's responsibility and is not verified by rmtool at all; this is a
narrower automated guarantee than originally scoped in 3b, traded off
explicitly rather than attempted unreliably.

## 3. Current architecture: what "ingest" and "geometry" mean today

Evidence, all from `src/rm_synthesis.f90` and `src/myfits_info.f90` on
`develop`:

- **Ingestion is hardcoded to one Q/U pair.** `cfg%infileQ`/`cfg%infileU`
  are single `character` fields in `rmsynth_config_t`
  (`src/rm_synthesis_mod.f90:109`), read once each via `read_cfg_keyval`
  and enforced as required keys
  (`src/rm_synthesis_mod.f90:2585-2589`). There is no array/list-of-files
  concept anywhere in the config parser.
- **Geometry, as validated today, means exact equality**, not overlap or
  reprojection. `myfits_info` (`src/myfits_info.f90:21-166`) reads per-cube
  `NAXIS`, `NAXES`, and `CRVAL1/2`, `CRPIX1/2`, `CDELT1/2` (RA/Dec) plus
  `CRVAL/CRPIX/CDELT` on whichever axis carries `CTYPE*=FREQ`
  (`myfits_info.f90:105-118`, `146-156`). `rm_synthesis.f90:594-666` then
  requires Q and U to agree on: frequency-axis index, `NAXIS`, and every
  `NAXES(i)` (pixel-for-pixel dimension equality); a second block,
  `rm_synthesis.f90:733-806`, separately enforces bit-identical
  `CRVAL/CRPIX/CDELT` on the RA, Dec, and frequency axes, hard `stop`-ing on
  any mismatch with no tolerance. None of this is a WCS-aware overlap
  computation — it's exact-value equality end to end. There is no rotation
  matrix (`PCi_j`/`CROTA2`) handling anywhere in the codebase — the existing
  check implicitly assumes no differential rotation between cubes, which
  holds trivially today because Q and U are always the same physical cube
  pair.
- **No resolution (beam) information is read or checked at all.** A search
  for `BMAJ`/`BMIN`/`BPA` across `src/` returns nothing. Matching resolution
  across the whole run is an unstated, unenforced assumption today, true only
  because there has only ever been one input band. Per the decision in §2,
  this stays true by design for multi-band too — no beam check will be
  added, since `BMAJ`/`BMIN`/`BPA` are not reliably present in multi-frequency
  FITS headers. The sole automated gate for multi-band geometry is exact
  RA/Dec pixel-grid equality (extending the existing NAXIS/CRVAL/CRPIX/CDELT
  checks below from Q-vs-U to N-cubes-vs-reference), loud-warn-and-refuse on
  mismatch.
- **Frequency sampling is assumed to be one linear grid per cube.** The
  spectral axis is described by a single `(CRVAL, CRPIX, CDELT)` triple per
  cube (`myfits_info.f90:142-156`); `rm_synthesis.f90:1511-1562` derives the
  absolute frequency of every channel, and hence every channel's λ² in
  `L_sq` (`rm_synthesis.f90:1482`, populated inside the same block), from
  that single triple plus `nz_totpix`. Bad/flagged channels are masked via
  `flag_arr`/`flag_arr_out` (weight-zeroed in the DFT), but they are still
  members of the one uniform grid — nothing in the codebase currently
  supports two cubes whose channels are not literally contiguous, evenly
  spaced samples of the same linear axis. Two bands from different receivers
  will essentially never satisfy that; combining them means building one
  merged, non-uniform frequency/λ² list, which today's `L_sq` construction
  does not do.
- **RM-axis sampling parameters (`ofac`, `fac`, `beg_rm`, `end_rm`, `nrm`) are
  single, run-wide scalars** (`rmsynth_config_t` fields,
  `src/rm_synthesis_mod.f90`; consumed at `rm_synthesis.f90:1615-1616` via
  `extract_general_setup`). Multi-band synthesis still wants exactly one RM
  axis for the whole run — this part needs no new concept, just a frequency
  list long enough to cover every band's channels.

## 4. Parallelism framework: what carries over, what doesn't

Reference: `docs/PARALLELISM.md`, `docs/ARCHITECTURE.md`.

**There is no MPI, and no multi-process/multi-node decomposition of any
kind** (`grep -rniE "mpi_init|use mpi|mpif" src/ Makefile build.sh` — zero
hits; the Slurm example at `scratch/slurm/run_rm-synthesis.sbatch:9-11`
requests `--nodes=1 --ntasks=1`). Every parallelism axis described below is
intra-node: OpenMP host threads and/or GPU target-offload. Multi-band
tomography therefore stays within the same single-node execution model —
there is no existing "distribute work across processes" mechanism to either
reuse or work around.

The existing decomposition is:
- **Tiles**: serial loop over 2D RA/Dec strips, sized by `mem_frac_ram`
  (`plan_tile`/`tile_plan_t`, `src/rm_synthesis_mod.f90:168-195`). Each tile
  read pulls the *full* channel span for that spatial footprint.
- **I/O parallelism**: `io_read_threads` splits **channels** of a single
  input file across independent CFITSIO handles; `io_write_threads` splits
  **RM bins** of the output cube across independent STREAM writers;
  `io_overlap` runs one tile's write concurrently with the next tile's
  read/compute on a background pthread.
- **Compute parallelism**: OpenMP/GPU parallelise over **pixels** (and
  pixel×RM-bin on GPU) within a tile. The channel loop is always sequential
  per pixel, and `cos_arr`/`sin_arr` are precomputed once for the whole run
  from the one `nz_out`-length frequency grid.

### What fits without disruption

| Aspect | Fits current framework? | Why |
|---|---|---|
| Spatial tiling (RA/Dec) | **Yes, unchanged** | Tiling is per-sky-position, independent of how many bands feed each pixel's spectrum. A tile is still "this RA/Dec footprint, all channels" — just now "all channels" spans multiple files. |
| Per-pixel OpenMP/GPU parallelism | **Yes, unchanged** | The compute kernel already treats the channel axis as one flat array of length `nz_out` per pixel; it does not care whether the channels came from one file or several, provided the merged spectrum + weights + λ² arrays are assembled correctly before the kernel runs. |
| `io_read_threads` | **Extends naturally** | Read parallelism is already "N independent CFITSIO handles per input file, splitting channels." With K bands, the natural extension is one such handle-set per band (or continue splitting by global channel range, now spanning multiple files) — same mechanism, larger fan-out. |
| `io_write_threads` / `io_overlap` | **Yes, unchanged** | Both operate purely on the *output* AMP/PHA cube and the RM axis, which stays single and run-wide regardless of how many input bands feed it. |
| RM chunking (`nrm_block_size`) | **Yes, unchanged** | Operates on the single merged RM axis; number of input bands is invisible to this stage. |

### What does not fit and needs new design

| Aspect | Fits current framework? | Why not, and what it needs |
|---|---|---|
| Config schema (single `infileQ`/`infileU`) | **No** | Needs a list-of-bands concept (§5) — every place `cfg%infileQ` is read once needs to become "for each band." This is schema and orchestration work, not parallelism work. **Must stay backward compatible** (requirement 2a, §5): existing single-band `cfg` files are required to keep working unedited, with bit-identical output — the new schema is additive (an optional `nbands` key gating a second, per-band key set), never a reinterpretation of today's keys. |
| Geometry validation | **No, must be extended, not reused as-is** | Today's check is Q-vs-U only; needs to become "all N cubes' RA/Dec pixel grid vs. a chosen reference" (band 1, per the §7 decision), still exact-equality, loud-warn-and-refuse on mismatch — **no** resolution/beam check is added (decided against in §2/§3: `BMAJ`/`BMIN`/`BPA` are unreliable in multi-frequency headers). |
| Frequency/λ² grid construction | **Mixed — better news than it first looks** | The DFT template kernel itself (`extract_general_setup`, `src/rm_synthesis_mod.f90:642-728`) computes `cos_arr`/`sin_arr` from each channel's λ² **individually** (a direct sum, not an FFT), so it does not algebraically require uniform λ² spacing — a concatenated, gapped, multi-band λ² array would compute correctly through this kernel essentially as-is. The actual gap is entirely upstream, in how that array is *produced*: `myfits_info` reads one `(CRVAL,CRPIX,CDELT)` triple per cube and `linspace`s it into a uniform ramp (`src/rm_synthesis.f90:1561-1563`) — there is no per-channel frequency table anywhere. Multi-band needs that replaced by **concatenating each band's own linspace-derived channel list** (each band can keep its own internal linear grid) into one merged, sorted array, with per-channel weights/flags carried through per band. This changes `nz_out`'s meaning from "one cube's channel count" to "sum of all bands' channel counts," touching every allocation sized by it (`data_arrQ/U`, `flag_arr_out`, `L_sq`, `cos_arr`/`sin_arr`, tile-local `specQ`/`specU`) — sizing/plumbing work, not a kernel rewrite. **Overlapping-band frequency ranges (§7 decision): no deduplication.** Per-channel weighting today is a uniform 0/1 flag (`flag_arr`/`flag_arr_out`, counted into `wsum` — "count of valid channels", `src/rm_synthesis.f90:2988-2989`), not a noise/sensitivity-based weight, so flat concatenation of every band's good channels already implements "weight by both" for free: an overlap region simply ends up with more equally-weighted channel terms in the same DFT sum, exactly like non-overlapping channels elsewhere in the run. If a noise-based per-channel weight is ever added later, the same merged-list design extends to it unchanged, since the kernel already takes an arbitrary per-channel weight. Separately, **`use_auto_rm_range=1`'s default RM-range heuristic does assume uniform spacing** (`dfreq = (freq_MHz(npts)-freq_MHz(1))/(npts-1)`, `src/rm_synthesis_mod.f90:656-687`) and would silently compute a wrong range/resolution across a multi-band gap; **§7 decision: forbidden outright for `nbands>1`** — see the RM-range diagnostic row below. |
| RM-range/resolution diagnostic for multi-band runs | **New, additive — not a parallelism concern** | With `use_auto_rm_range` forbidden for `nbands>1` (§7), the user must supply `beg_rm`/`end_rm`/`nrm` explicitly, but has no easy way to know what's actually achievable from their specific band combination — non-uniform λ² spacing across bands means the classical uniform-spacing formulas (already used for the existing single-band auto-range heuristic, `src/rm_synthesis_mod.f90:656-687`) don't directly apply. **Decision (confirmed with user, 2026-07-20):** compute and log (stdout + run log) a generalized theoretical RM resolution and un-aliased RM range/span from the merged, non-uniform λ² array — `RM_res` from the **maximum separation in λ²** across the full merged coverage, `RM_span` from the **width of the smallest (finest) individual channel in λ²** — as an informational diagnostic only, guiding the user's choice of `beg_rm`/`end_rm`/`nrm` rather than auto-selecting a range. The exact generalized (non-uniform-spacing) formulas are to be sourced from Wasim Raja's PhD thesis (RRI digital repository) at implementation time; the classical uniform-spacing case already implemented (`src/rm_synthesis_mod.f90:656-687`, Brentjens & de Bruyn-style relations) is the special case this generalizes from. |
| Tile read stage | **Partially — needs a per-band read loop, not a redesign** | `tile_read` currently issues one `FTGSVE` call (or `io_read_threads`-many) per input file. With K bands it becomes a loop of K such call-groups into disjoint slices of one enlarged `specQ`/`specU` buffer — additive complexity, not a new decomposition axis. Bad-channel/mask handling (`flag_arr`) needs to become per-band-aware so a channel flagged bad in band 2 doesn't collide with band 1's indexing. |
| Bias correction / Q-U bias fields (`resiQ`, `slopeQ`, `resiU`, `slopeU`, `infileI`) | **Decided: per-band** | These are currently single scalars/one I-cube for the whole run. **Decision (confirmed with user, 2026-07-20):** Q-U bias correction is physically an instrumental effect and must be done **per band** — each band needs its own `resiQ`/`slopeQ`/`resiU`/`slopeU` and, when `remove_qu_bias=y`, its own Stokes-I cube (`infileI`/`path_I`). This is config-schema and per-band-loop plumbing, not a parallelism concern — the bias correction itself is applied per-channel before the DFT sum, so it composes with the merged-frequency-list design the same way bad-channel flags do. |

### Summary

The tiling, I/O-thread, and RM-chunk decomposition axes are **orthogonal to
"how many bands feed the frequency axis"** and need no rework — and,
better than initially assumed, neither does the DFT template kernel itself,
which is already per-channel rather than assuming a uniform grid. The real
cost of this feature is concentrated in a narrower band than a first pass
would suggest: the **ingestion/validation/frequency-assembly** layer (turning
"one cube's grid" into "N cubes' grids, merged into one list, each channel
still tagged with which band/file it came from for I/O purposes") plus the
**auto-RM-range heuristic**, which does assume uniform spacing and needs an
explicit multi-band guard. Once the merged per-channel frequency/λ² list and
weight array exist, they can be handed to the existing tile/compute/write
pipeline essentially unchanged. All of this stays within the existing
single-node OpenMP/GPU execution model — there is no distributed-computing
dimension to design around.

## 5. Sketch of a config schema (not final, no code written)

Current single-band schema (`cfg/rmsynth.cfg`):
```
infileQ = Q_cube.fits
infileU = U_cube.fits
resiQ = 0.0  slopeQ = 0.0  resiU = 0.0  slopeU = 0.0
infileI = I_cube.fits  path_I = /path/to/data/
```

### Backward compatibility (requirement 2a): two explicit, non-overlapping modes

**Design: `nbands` is a new, optional key, defaulting to 1.** Whether it is
present (and `>1`) selects which of two mutually exclusive key sets is read
— there is no reinterpretation of old keys and no silent merging of the two
forms, so there is no ambiguity to get wrong:

- **`nbands` absent, or `nbands = 1` (legacy/default mode).** The unsuffixed
  keys — `infileQ`, `infileU`, `resiQ`/`slopeQ`/`resiU`/`slopeU`,
  `infileI`/`path_I` — are read exactly as today, through the same
  single-band code path, with **zero observable behaviour change**. Every
  existing `cfg` file in `cfg/`, `tests/*.cfg`, and any user's own configs
  never sets `nbands`, so they fall into this branch automatically and keep
  working unedited. Setting the new suffixed, per-band keys (below) in this
  mode is a config error (ambiguous intent), not silently ignored.
- **`nbands > 1` (multi-band mode).** The unsuffixed keys become invalid in
  this mode (config error if set — same reasoning: avoid silent ambiguity
  about which value wins), replaced by per-band suffixed keys:
  ```
  nbands  = 2
  # reference_band: whose RA/Dec pixel grid every other band is validated
  # against (band 1 by default, per the §7 decision).
  reference_band = 1

  infileQ_1 = low_Q.fits   infileU_1 = low_U.fits
  resiQ_1 = 0.0  slopeQ_1 = 0.0  resiU_1 = 0.0  slopeU_1 = 0.0
  infileI_1 = low_I.fits  path_I_1 = /path/to/data/

  infileQ_2 = high_Q.fits  infileU_2 = high_U.fits
  resiQ_2 = 0.0  slopeQ_2 = 0.0  resiU_2 = 0.0  slopeU_2 = 0.0
  infileI_2 = high_I.fits  path_I_2 = /path/to/data/
  ```
  (or a comma-separated list form — either is a config-parser syntax choice,
  independent of the feasibility questions above and of the backward-compat
  design; picking the exact syntax is an implementation-time decision).
  Bias-correction fields are per-band in this mode (§4, §7 decision).
  `use_auto_rm_range` is rejected at config-validation time whenever
  `nbands>1` (§4, §7 decision) — `beg_rm`/`end_rm`/`nrm` become required
  keys in that case, informed by the logged RM-range/resolution diagnostic.

**Correctness gate for this requirement**, to be carried into the phase-1
ticket (§6): every existing `tests/*.cfg` and `cfg/*.cfg` file, run unedited
against the new code, must produce **bit-identical output** to the current
`develop` baseline — the same gate this repo already applies to structural
refactors (`planning/ENCAPSULATION_REFACTOR_PLAN.md`'s "zero change in
observable behaviour" constraint). Concretely, this likely means the legacy
(`nbands`-absent) path should route to the *same* single-band code that
exists today, rather than a "multi-band code with N always forced to 1" —
the latter risks subtly different floating-point summation order or
allocation sizing even when logically equivalent; the former guarantees
bit-identical output by construction.

## 6. Recommended phasing

1. **N-file ingestion + geometry validation, gated on backward compatibility
   first.** Smallest useful slice: prove the "N files in (config schema,
   §5), RA/Dec pixel grid validated against band 1, loud refusal on
   mismatch, existing pipeline unchanged downstream" shape works, with the
   fewest moving parts. No frequency merge yet — could initially even reject
   overlapping/multi-band frequency coverage to isolate ingestion-plumbing
   risk from frequency-merge risk. **Correctness gate includes requirement
   2a's bit-identical check** (§5): every existing `tests/*.cfg`/`cfg/*.cfg`
   file, unedited, must produce identical output on the new code before this
   phase is considered done — the same bar as this repo's existing
   refactor-correctness gates.
2. **Frequency/λ² merge + per-band bad-channel handling.** The structural
   work identified as the biggest remaining gap in §4 — concatenate each
   band's channel list (no deduplication needed in overlaps, §4/§7) into one
   merged, sorted λ²/weight array sized by the new `nz_out` meaning.
3. **RM-range/resolution diagnostic + `use_auto_rm_range` guard.** Reject
   `use_auto_rm_range=1` whenever `nbands>1`; compute and log the
   generalized (non-uniform-λ²-spacing) `RM_res`/`RM_span` diagnostic from
   Wasim Raja's PhD thesis formulas, to guide the user's explicit
   `beg_rm`/`end_rm`/`nrm` choice.
4. **Per-band bias-correction implementation** (`resiQ`/`slopeQ`/`resiU`/
   `slopeU`/`infileI` per band, §4/§7 decision).
5. **Multi-band-aware diagnostics**: swim-lane/log output currently reports
   one `bytes=` figure per tile read; extending it to show per-band
   breakdown is a nice-to-have, not a blocker.

Each phase should get its own ticket(s) in the style of
`planning/ENCAPSULATION_REFACTOR_PLAN.md` / `IO_PARALLEL_OPTIMISATION_PLAN.md`
(Objective/Scope/Change Set/Correctness Gate/Rollback Criteria/Effort) once
this feasibility plan is agreed, rather than being written speculatively here.

## 7. Decisions recorded (confirmed with user, 2026-07-20)

All five open questions from the original draft of this document are now
resolved:

1. **Overlapping frequency ranges between bands**: no error, no explicit
   merge policy needed — "weight by both" falls out for free from flat
   concatenation, since today's per-channel weighting is a uniform 0/1 flag
   (`wsum`, `src/rm_synthesis.f90:2988-2989`), not a noise/sensitivity
   weight. See §4's frequency/λ² grid construction row.
2. **Bias correction fields** (`resiQ`/`slopeQ`/`resiU`/`slopeU`/`infileI`):
   **per-band**, since Q-U bias is an instrumental effect specific to each
   receiver/band. See §4, §5.
3. **Resolution-mismatch tolerance**: **no `BMAJ`/`BMIN`/`BPA` check at
   all** — those keywords are unreliable in multi-frequency FITS headers, so
   a header-based beam-equality check would be brittle exactly where it's
   needed most. The only automated gate is **exact RA/Dec pixel-grid
   equality**, loud-warn-and-refuse on mismatch. Actual resolution matching
   between bands is the caller's responsibility, unverified by rmtool. See
   §2, §3, §4.
4. **Reference geometry for validating N cubes**: **the first band
   listed** (`reference_band` defaults to band 1, config-overridable — §5)
   — the rationale being that all bands are expected to already share one
   geometry by construction (§2's pre-matched-geometry decision), so which
   one is nominally "the reference" is a validation-order convenience, not a
   meaningful scientific choice.
5. **`use_auto_rm_range` for multi-band runs**: **forbidden outright**
   whenever `nbands>1` — the existing heuristic assumes uniform channel
   spacing (§3, §4) and would silently miscompute across a multi-band gap.
   In its place, rmtool computes and logs (stdout + run log) a **diagnostic
   RM resolution (`RM_res`) and un-aliased RM range/span (`RM_span`)**
   generalized to non-uniform λ² spacing: `RM_res` from the **maximum
   separation in λ²** across the full merged band coverage, `RM_span` from
   the **width of the smallest (finest) individual channel in λ²** — per
   Wasim Raja's PhD thesis (RRI digital repository), generalizing the
   uniform-spacing relations already implemented for the single-band case
   (`src/rm_synthesis_mod.f90:656-687`). This is informational only: rmtool
   does **not** auto-run the full theoretical RM range for multi-band —
   `beg_rm`/`end_rm`/`nrm` remain required, explicit user choices, guided by
   (not overridden by) the logged diagnostic.

## 8. Non-goals for this effort

- WCS reprojection/regridding of mismatched cubes (rejected for v1, §2).
- Convolution to a common resolution (rejected for v1, §2 — resolution
  matching is validated, not performed).
- Any `BMAJ`/`BMIN`/`BPA`-based resolution-equality check (rejected
  outright, §2/§7 — not merely deferred: multi-frequency FITS headers don't
  carry these reliably enough to check).
- Mosaicking (combining cubes covering *different* sky regions into one
  footprint) — out of scope; this effort is about combining different
  *frequency* coverage of the *same* sky region.
- Auto-selecting or auto-running the full theoretical RM range for
  multi-band data (rejected, §7 point 5) — rmtool logs the achievable
  `RM_res`/`RM_span` as a diagnostic; the user still chooses
  `beg_rm`/`end_rm`/`nrm` explicitly.
- Making `use_auto_rm_range=1` itself band-gap-aware (rejected in favor of
  forbidding it outright for `nbands>1`, §7 point 5) — a band-aware version
  of the auto-range heuristic is not being pursued.
- Any change to the numerical RM-synthesis kernel itself
  (`extract_general_setup`, `tile_extract_gpu_rm_blocked`) beyond how many
  channels it's handed — same guardrail this repo already applies in
  `planning/ENCAPSULATION_REFACTOR_PLAN.md`.
