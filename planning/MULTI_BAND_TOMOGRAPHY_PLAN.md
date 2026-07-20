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
2b. **(added, confirmed with user, 2026-07-20)** The ingestion/validation/
   compute pipeline must be **unified, not duplicated**: there is exactly
   one code path, parameterized by `nbands`, and the legacy single-band case
   is equivalent to running that one path with `nbands=1` — not a separate,
   independently-maintained legacy subroutine. See §5 for what this trades
   away (a weaker, empirically-verified rather than by-construction
   bit-identical guarantee) and why that trade was accepted.
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
| Config schema (single `infileQ`/`infileU`) | **No** | Needs a list-of-bands concept (§5) — every place `cfg%infileQ` is read once needs to become "for each band." This is schema and orchestration work, not parallelism work. **Must stay backward compatible, via one unified pipeline, not two** (requirements 2a/2b, §5): `nbands` defaults to 1, legacy unsuffixed keys are parse-time aliases for band 1's keys, and there is exactly one ingestion/validation/compute pipeline for all `nbands` values — bit-identical output for existing configs is guaranteed only by the test-sweep correctness gate (§5), not by construction. |
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

### Backward compatibility (requirements 2a/2b): one unified pipeline, `nbands` defaults to 1

**Revised design (confirmed with user, 2026-07-20, superseding an earlier
two-path draft — see the note at the end of this subsection): `nbands` is a
new, optional key, defaulting to 1, and there is exactly one ingestion →
geometry-validation → frequency-merge → tile/compute/write pipeline for
every value of `nbands`, including 1.** The legacy single-band case is not a
separately maintained branch; it *is* the `nbands=1` case of the general
pipeline.

- **Config-parser level (the only place backward compatibility needs
  explicit handling):** the unsuffixed legacy keys — `infileQ`, `infileU`,
  `resiQ`/`slopeQ`/`resiU`/`slopeU`, `infileI`/`path_I` — remain accepted
  and are resolved, once, at parse time, as **aliases for band 1's suffixed
  keys** (`infileQ` ≡ `infileQ_1`, etc.) into one internal per-band array
  representation (conceptually `band(1:nbands)`, each holding the fields
  `rmsynth_config_t` has today, just multiplied out). Downstream of parsing
  there is only ever "an array of `nbands` bands" — no branch on whether the
  legacy or new key spelling was used.
  ```
  # legacy spelling (nbands absent/1) -- parses to band(1)%infileQ = ...
  infileQ = Q_cube.fits
  infileU = U_cube.fits

  # equivalent explicit spelling, same result
  nbands = 1
  infileQ_1 = Q_cube.fits
  infileU_1 = U_cube.fits
  ```
  Mixing the two spellings for the same band (e.g. both `infileQ` and
  `infileQ_1` set) remains a config-parse error — that guard is about
  avoiding ambiguous *input*, not about maintaining two downstream code
  paths, so it doesn't conflict with unification.
- **Multi-band configs** (`nbands > 1`) use only the suffixed form, since
  there's no single unsuffixed name that could mean "band 2":
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
  (or a comma-separated list form — a config-parser syntax choice,
  independent of everything else here). `use_auto_rm_range` is rejected at
  config-validation time whenever `nbands>1` (§4, §7 decision) —
  `beg_rm`/`end_rm`/`nrm` become required keys in that case, informed by the
  logged RM-range/resolution diagnostic.

**What "unified" costs, and how the correctness gate compensates.** Under
the earlier two-path draft, the legacy case got its bit-identical guarantee
*for free*, by construction — it literally called the untouched existing
subroutine. Unifying gives up that free guarantee: the general pipeline is
new code, and floating-point non-associativity means a loop that is
"logically equivalent" to today's for `nbands=1` (e.g. an outer per-band
loop that happens to run once) is not automatically bit-identical to it —
compiler vectorization/instruction-scheduling/FMA-contraction decisions can
depend on loop shape even when the arithmetic sequence is mathematically the
same. **This was a conscious trade the user made explicitly**, accepting the
added verification burden in exchange for a single maintained pipeline. The
correctness gate is therefore promoted from a backstop to the **primary and
only** mechanism ensuring requirement 2a: every existing `tests/*.cfg` and
`cfg/*.cfg` file, run unedited against the new code, must produce
**bit-identical output** to the current `develop` baseline (the same "zero
change in observable behaviour" bar `planning/ENCAPSULATION_REFACTOR_PLAN.md`
already applies to structural refactors in this codebase) — this is no
longer a sanity check on top of a by-construction guarantee, it is the
guarantee. Two implementation practices worth carrying into the ticket that
does this work, to keep the `nbands=1` case as close to "no extra floating
point operations" as the unification allows:
- The per-band frequency/λ² construction (`myfits_info` → `linspace`, one
  call per band) should stay a straight per-band call with no change to its
  internal arithmetic; only the *assembly* of per-band results into the
  merged list is new.
- For `nbands=1`, the merge/concatenation step (whatever form it takes for
  `nbands>1` — sort, interleave, tag-by-band) should reduce to using band
  1's array directly rather than passing it through generic merge logic
  that happens to be a no-op for one input — a real short-circuit, not
  merely an algorithm that's expected to behave like one.

(**Earlier draft, superseded:** a prior version of this document proposed
two explicit, mutually exclusive code paths gated on `nbands==1` vs. `>1`,
matching the pattern this repo already uses for `io_read_threads`/
`io_write_threads`/`io_overlap` defaulting to their serial/off behaviour via
an unchanged code branch. That remains the lower-risk option if the
verification cost of the unified design proves expensive in practice; it is
recorded here in case a future revision needs to fall back to it.)

## 6. Recommended phasing

1. **Unified N-band ingestion + geometry validation, `nbands` defaulting to
   1.** Smallest useful slice: build the single ingestion/validation
   pipeline (config-parser aliasing of legacy keys to band 1, §5;
   per-band-array internal representation; RA/Dec pixel grid validated
   against `reference_band`) and prove it works for `nbands=1` before
   exercising `nbands>1` at all. No frequency merge yet — could initially
   even reject `nbands>1` outright to isolate ingestion-plumbing risk from
   frequency-merge risk. **This phase carries the full weight of
   requirements 2a/2b's correctness gate** (§5), since the unified design
   means there is no by-construction fallback: every existing
   `tests/*.cfg`/`cfg/*.cfg` file, unedited, must produce bit-identical
   output on the new code before this phase is considered done — the same
   bar as this repo's existing refactor-correctness gates, but now the sole
   mechanism rather than a backstop.
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
(Objective/Scope/Change Set/Correctness Gate/Rollback Criteria/Effort).
**T0 and T1 (covering the start of phase 1) are now written — see §9.** Later
tickets (phase 1's remainder, phases 2-5) are deliberately not written yet;
each should be scoped once the ticket(s) before it have landed and their
Evidence sections are filled in, the same incremental discipline
`ENCAPSULATION_REFACTOR_PLAN.md` used.

## 7. Decisions recorded (confirmed with user, 2026-07-20)

All five open questions from the original draft of this document are now
resolved, plus one follow-up decision (0) reached after the rest:

0. **Backward-compatibility mechanism (requirement 2b, revising an earlier
   draft of §5): unify into one pipeline rather than maintain two code
   paths.** `nbands` defaults to 1; legacy single-band configs are the
   `nbands=1` case of the same general pipeline, not a separately kept
   legacy branch. This was a deliberate trade: it gives up the "bit-identical
   by construction" guarantee two explicit paths would have given the
   legacy case for free, in exchange for not duplicating the
   ingestion/validation/compute logic. The correctness gate (bit-identical
   output on every existing `tests/*.cfg`/`cfg/*.cfg`, §5) is promoted from
   backstop to sole safety mechanism as a result. See §5 for the full
   design and the superseded two-path alternative, kept there as a fallback
   if the unified design's verification cost proves too high in practice.

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
6. **Multi-band test data: synthetic only, not real observational cubes**
   (confirmed with user, 2026-07-20). All multi-band correctness testing —
   geometry-validation accept/reject paths, frequency/λ² merge, the
   `RM_res`/`RM_span` diagnostic — is validated against synthetic Q/U cubes
   with known injected point sources at known RM, extending the existing
   `tests/make_test_cubes.py` generator (already used by `tests/run_tests.sh`
   for single-band tests, validated via `tests/check_rm_peak.py` recovering
   the known-truth RM and `tests/compare_cubes.py` for bit-identical/diff
   comparisons). **This is an explicit, narrower scope than the
   encapsulation-refactor effort accepted for itself**: that project's T3b
   needed a real production-scale Setonix run to validate beyond what
   synthetic in-suite tests could show (see project memory
   `project_jennifer_t3b_validation`); no equivalent real-data validation
   step is planned here. A real multi-band production run remains something
   the user can do separately once this lands, but it is not a gate any
   ticket in this plan depends on passing.

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

## 9. Tickets

Ticket format follows this repo's existing convention
(`planning/ENCAPSULATION_REFACTOR_PLAN.md`, `planning/IO_PARALLEL_OPTIMISATION_PLAN.md`):
Objective / Scope / Change Set / Correctness Gate / Rollback Criteria /
Effort, each getting an **Evidence (...)** section appended once done.

---

### T0 — Baseline Lock

- **Objective:** Freeze an exact, reproducible reference of current
  `develop` behaviour before any multi-band code changes land, so every
  later ticket's bit-identical correctness gate (§5, §7 decision 0) has a
  concrete baseline to diff against — the same purpose T0 served in
  `planning/ENCAPSULATION_REFACTOR_PLAN.md`, adapted for this branch. This
  step matters more here than it did there: that effort's baseline was a
  safety net for a refactor with no algorithmic change intended anywhere;
  this effort's own §7 decision 0 explicitly gave up the "bit-identical by
  construction" guarantee, making this baseline the *primary* correctness
  instrument, not a backstop.
- **Scope:** Measurement only, on `multi-band-tomography` at its current
  tip (commit `a2417f6` at the time this ticket is written — re-confirm the
  exact commit when the ticket is actually executed, since more planning
  commits may land first). No source changes.
- **Change Set:** None.
- **Correctness Gate:**
  - Clean build via `bash scratch/make_all.sh` (or this repo's current
    equivalent full-matrix build script) across all four variants
    (CPU-serial, CPU-OMP, GPU-offload, GPU-offload-hostomp); record the
    actual warning count verbatim, don't assume zero.
  - Full `bash tests/run_tests.sh` run; record PASS/FAIL/SKIP counts
    verbatim.
  - Archive the complete output FITS set (AMP/PHA/MASK/NVALID/etc.) plus run
    logs for every `tests/*.cfg` under a new
    `scratch/baseline_multiband/` directory (gitignored, local-only,
    mirroring `scratch/baseline_encapsulation/`'s existing pattern) — the
    literal byte-for-byte reference every later ticket in this plan diffs
    against via `tests/compare_cubes.py --exact`.
  - Record the exact git commit hash this baseline was built from in a
    `scratch/baseline_multiband/T0_MANIFEST.md`, mirroring
    `scratch/baseline_encapsulation/T0_MANIFEST.md`'s existing format.
- **Rollback Criteria:** N/A (measurement only).
- **Effort:** 0.5 session.

---

### T1 — Unified `nbands`-Parameterized Ingestion + Geometry Validation

- **Objective:** Introduce the `nbands` config key (default 1) and the
  single, unified per-band ingestion/validation pipeline (§5, §7 decision
  0) — proving the `nbands=1` case is bit-identical to the T0 baseline
  before any frequency-merge logic (phase 2) exists at all. This is
  deliberately the narrowest slice that exercises the new schema and the
  new geometry-validation loop end-to-end without touching the RM-synthesis
  numerics.
- **Scope:**
  - `src/rm_synthesis_mod.f90:106-153` (`rmsynth_config_t`): add `nbands`
    (default 1) and `reference_band` (default 1) fields. Add a new small
    per-band derived type (e.g. `band_cfg_t`) holding the per-band fields
    identified in §4/§5 (`infileQ`, `infileU`, `resiQ`, `slopeQ`, `resiU`,
    `slopeU`, `infileI`, `path_I`) and an allocatable array of it, sized to
    `nbands`, inside `rmsynth_config_t`.
  - `read_cfg_keyval` (`src/rm_synthesis_mod.f90:1632` onward): add
    `nbands`/`reference_band` parsing; add the legacy-key-aliasing logic
    (§5) that populates `band(1)` from unsuffixed `infileQ`/`infileU`/etc.
    when `nbands` is absent or 1, and from suffixed `infileQ_1`, `infileQ_2`
    ... when `nbands>1`; reject mixed spellings (both a legacy and a
    suffixed key set for the same run) as a config-parse error, per §5.
    Required-key enforcement (today at `:2585-2589` for `infileQ`/`infileU`)
    extends to every configured band.
  - `src/rm_synthesis.f90:572-667` (today's `myfits_info` calls + Q-vs-U
    NAXIS/dimension check) and `:733-806` (today's Q-vs-U WCS-value exact
    match): generalize from "Q vs U, exactly two cubes" to "every band's Q
    and U vs `reference_band`'s geometry, N cubes" — same exact-equality
    philosophy (§3, §7 decision 3), same loud-refuse-on-mismatch behaviour,
    just looped over `nbands` instead of hardcoded to one pair.
  - For `nbands>1` specifically: after geometry validation passes, stop
    with a clear, explicit "multi-band frequency merge not yet implemented"
    message rather than attempting synthesis — phase 2's job, out of scope
    here. This keeps T1's blast radius to ingestion/validation plumbing
    only.
  - Test fixtures: extend `tests/make_test_cubes.py` to optionally emit a
    second synthetic Q/U band (distinct frequency range from the existing
    550-750 MHz GMRT-like band, e.g. an 800-950 MHz band, same RA/Dec
    geometry and same injected point sources) plus a deliberately
    geometry-mismatched variant (different `CRVAL1`/`CDELT1` or pixel
    count) — enough to exercise both the accept and loud-refuse paths of
    the new N-band geometry validation. Per §7 decision 6, this stays
    synthetic-only; no real multi-band cubes are sourced for this ticket.
- **Correctness Gate:**
  - **Legacy/default path (`nbands` absent, and explicit `nbands=1`):**
    every `tests/*.cfg`/`cfg/*.cfg` file, run unedited, produces
    bit-identical output (`tests/compare_cubes.py --exact`) against the T0
    baseline archive. This is the ticket's central gate, carrying the full
    weight of requirements 2a/2b (§5, §7 decision 0) — there is no
    by-construction fallback if this fails.
  - `tests/run_tests.sh`'s existing pass/fail counts are unchanged from T0.
  - New multi-band fixture tests (added to `tests/run_tests.sh` or a
    parallel script): a matched-geometry two-band config passes validation
    and reaches (and stops cleanly at) the "not yet implemented" message; a
    mismatched-geometry two-band config is loudly refused before any
    compute begins.
  - A config mixing legacy and suffixed keys is rejected with a clear error
    at parse time.
- **Rollback Criteria:** If bit-identical output cannot be achieved for the
  `nbands=1` case within this ticket's effort budget, roll back to before
  this ticket rather than merging a change that silently breaks existing
  users — re-evaluate whether the unified-pipeline decision (§7 decision 0)
  needs revisiting in favour of the superseded two-path design kept in §5.
- **Effort:** 1.5-2 sessions (config-parser rework + new derived type +
  N-band validation loop + synthetic multi-band fixture generation, gated
  by a bit-identical sweep that itself takes real wall-clock time to run
  and diff).
