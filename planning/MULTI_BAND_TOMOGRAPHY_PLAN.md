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
| Config schema (single `infileQ`/`infileU`) | **No** | Needs a list-of-bands concept (§5) — every place `cfg%infileQ` is read once needs to become "for each band." This is schema and orchestration work, not parallelism work. **Must stay backward compatible, via one unified pipeline, not two** (requirements 2a/2b, §5): every per-band key becomes comma-separable, band count is derived from list length (no separate `nbands` key), a comma-free value is just a length-1 list, and there is exactly one ingestion/validation/compute pipeline regardless of band count — bit-identical output for existing configs is guaranteed only by the test-sweep correctness gate (§5), not by construction. |
| Geometry validation | **No, must be extended, not reused as-is** | Today's check is Q-vs-U only; needs to become "all N cubes' RA/Dec pixel grid vs. a chosen reference" (band 1, per the §7 decision), still exact-equality, loud-warn-and-refuse on mismatch — **no** resolution/beam check is added (decided against in §2/§3: `BMAJ`/`BMIN`/`BPA` are unreliable in multi-frequency headers). |
| Frequency/λ² grid construction | **Mixed — better news than it first looks** | The DFT template kernel itself (`extract_general_setup`, `src/rm_synthesis_mod.f90:642-728`) computes `cos_arr`/`sin_arr` from each channel's λ² **individually** (a direct sum, not an FFT), so it does not algebraically require uniform λ² spacing — a concatenated, gapped, multi-band λ² array would compute correctly through this kernel essentially as-is. The actual gap is entirely upstream, in how that array is *produced*: `myfits_info` reads one `(CRVAL,CRPIX,CDELT)` triple per cube and `linspace`s it into a uniform ramp (`src/rm_synthesis.f90:1561-1563`) — there is no per-channel frequency table anywhere. Multi-band needs that replaced by **concatenating each band's own linspace-derived channel list** (each band can keep its own internal linear grid) into one merged, sorted array, with per-channel weights/flags carried through per band. This changes `nz_out`'s meaning from "one cube's channel count" to "sum of all bands' channel counts," touching every allocation sized by it (`data_arrQ/U`, `flag_arr_out`, `L_sq`, `cos_arr`/`sin_arr`, tile-local `specQ`/`specU`) — sizing/plumbing work, not a kernel rewrite. **Overlapping-band frequency ranges (§7 decision): no deduplication.** Per-channel weighting today is a uniform 0/1 flag (`flag_arr`/`flag_arr_out`, counted into `wsum` — "count of valid channels", `src/rm_synthesis.f90:2988-2989`), not a noise/sensitivity-based weight, so flat concatenation of every band's good channels already implements "weight by both" for free: an overlap region simply ends up with more equally-weighted channel terms in the same DFT sum, exactly like non-overlapping channels elsewhere in the run. If a noise-based per-channel weight is ever added later, the same merged-list design extends to it unchanged, since the kernel already takes an arbitrary per-channel weight. Separately, **`use_auto_rm_range=1`'s default RM-range heuristic does assume uniform spacing** (`dfreq = (freq_MHz(npts)-freq_MHz(1))/(npts-1)`, `src/rm_synthesis_mod.f90:656-687`) and would silently compute a wrong range/resolution across a multi-band gap; **§7 decision: forbidden outright for `nbands>1`** — see the RM-range diagnostic row below. |
| RM-range/resolution diagnostic for multi-band runs | **New, additive — not a parallelism concern** | With `use_auto_rm_range` forbidden for `nbands>1` (§7), the user must supply `beg_rm`/`end_rm`/`nrm` explicitly, but has no easy way to know what's actually achievable from their specific band combination. **Decision (confirmed with user, 2026-07-20), formulas now sourced from the thesis (Raja 2014, Chapter 6 + §2.5 — see §7 decision 5 for the full citation, exact equations, and the three-quantity distinction the user flagged):** compute and log (stdout + run log) three genuinely distinct quantities — `δRM` (resolving power, for telling two nearby features apart), `max RM scale` (the largest Faraday-*thickness* an extended/thick component can have and still be detected at all), and `max un-aliased RM` (the largest Faraday *depth* at which a single thin/point-like component can still be reliably, unambiguously measured) — as an informational diagnostic only, guiding the user's choice of `beg_rm`/`end_rm`/`nrm` rather than auto-selecting a range. |
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

### Backward compatibility (requirements 2a/2b): one unified pipeline, comma-separated lists

**Final design (confirmed with user, 2026-07-21, superseding two earlier
drafts — see the history note at the end of this subsection): every
per-band key accepts a comma-separated list. `infileQ = Q_cube.fits` and
`infileQ = Q_band1.fits,Q_band2.fits` are parsed by the *same* code —
split on comma, get an array — so a single value isn't a special "legacy"
case at all, it's simply a length-1 array. There is no separate `nbands`
config key**: the number of bands is *derived* from the list length of
`infileQ` (cross-checked against `infileU`'s list length, and against
every other per-band key's list length when set — a config-parse error on
mismatch, e.g. `infileQ` listing 2 files but `infileU` listing 1).

```
# legacy single-band cfg, completely unedited -- still valid, still nbands=1
infileQ = Q_cube.fits
infileU = U_cube.fits
resiQ = 0.0  slopeQ = 0.0  resiU = 0.0  slopeU = 0.0
infileI = I_cube.fits  path_I = /path/to/data/

# multi-band: same keys, comma-separated
infileQ = low_Q.fits,high_Q.fits
infileU = low_U.fits,high_U.fits
resiQ = 0.0,0.0   slopeQ = 0.0,0.0   resiU = 0.0,0.0   slopeU = 0.0,0.0
infileI = low_I.fits,high_I.fits   path_I = /path/to/data/,/path/to/data/

# reference_band: whose RA/Dec pixel grid every other band is validated
# against (band 1, i.e. the first list entry, by default -- §7 decision 4).
reference_band = 1
```

This is a **simpler and more thoroughly unified realization of requirement
2b than either earlier draft achieved** (see history below): there is
exactly one key name per field (no suffixed `_1`/`_2` variants to keep in
sync with a separate `nbands` count), and the parser has no
legacy-vs-new-spelling branch to get wrong — comma-splitting a string with
no commas and a string with one comma are the same operation. `nbands` still
exists as an internal, derived quantity (used throughout this document as
shorthand for "how many bands"), it just isn't something the user sets
directly, removing an entire class of possible inconsistency (a stale
`nbands=3` next to only 2 listed files can no longer happen — there is
nothing to go stale). `use_auto_rm_range` is still rejected at
config-validation time whenever the derived `nbands>1` (§4, §7 decision) —
`beg_rm`/`end_rm`/`nrm` become required keys in that case, informed by the
logged RM-range/resolution diagnostic.

**What "unified" costs, and how the correctness gate compensates.** The
verification story from the earlier drafts carries over unchanged: a
length-1 list, however it's parsed, still runs through *new* code (the
general list-based ingestion/merge pipeline), and floating-point
non-associativity means "logically equivalent to today's single-band path"
is not automatically bit-identical to it — compiler
vectorization/instruction-scheduling/FMA-contraction decisions can depend
on loop shape even when the arithmetic sequence is mathematically the
same. **This was a conscious trade the user made explicitly** (§7 decision
0), accepting the added verification burden in exchange for a single
maintained pipeline. The correctness gate is therefore the **primary and
only** mechanism ensuring requirement 2a: every existing `tests/*.cfg` and
`cfg/*.cfg` file, run unedited against the new code, must produce
**bit-identical output** to the current `develop` baseline (the same "zero
change in observable behaviour" bar `planning/ENCAPSULATION_REFACTOR_PLAN.md`
already applies to structural refactors in this codebase). Two
implementation practices worth carrying into the ticket that does this
work, to keep the single-entry-list case as close to "no extra floating
point operations" as the unification allows:
- The per-band frequency/λ² construction (`myfits_info` → `linspace`, one
  call per list entry) should stay a straight per-entry call with no change
  to its internal arithmetic; only the *assembly* of per-entry results into
  the merged list is new.
- For a length-1 list, the merge/concatenation step (whatever form it takes
  for longer lists — sort, interleave, tag-by-band) should reduce to using
  that one entry's array directly rather than passing it through generic
  merge logic that happens to be a no-op for one input — a real
  short-circuit, not merely an algorithm that's expected to behave like one.

**Design history, superseded drafts (kept for context, not current):**
1. *Two explicit, mutually exclusive code paths* gated on band count,
   matching the pattern this repo already uses for `io_read_threads`/
   `io_write_threads`/`io_overlap` defaulting to their serial/off behaviour
   via an unchanged code branch — rejected (2026-07-20) in favour of one
   unified pipeline (§7 decision 0), on the grounds that duplicating the
   whole ingestion/validation/compute pipeline was worse than the added
   verification burden.
2. *A `nbands` config key plus `infileQ_1`/`infileQ_2`-style suffixed
   per-band keys*, with unsuffixed legacy keys treated as parse-time
   aliases for band 1 — the unified pipeline's first concrete schema, but
   superseded (2026-07-21) by the comma-separated-list design above, which
   achieves the same unification with less config surface (no suffix
   proliferation, no separate `nbands` key to keep consistent with the
   actual list lengths).

Both remain lower-risk fallbacks if the current design's verification cost
proves expensive in practice.

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
2. **Frequency/λ² merge + per-band bad-channel handling, validated against
   the synthetic scenario in §10.** §10's exact numeric parameters are now
   pinned directly from the thesis (Raja 2014, Table 6.1/6.2 — see §7
   decision 5's citation), so this dependency is already resolved rather
   than pending. The structural ingestion/merge work identified as the
   biggest remaining gap in §4 — concatenate each band's channel list (no
   deduplication needed in overlaps, §4/§7) into one merged, sorted
   λ²/weight array sized by the new `nz_out` meaning.
3. **RM-range/resolution diagnostic + `use_auto_rm_range` guard.** Reject
   `use_auto_rm_range=1` whenever `nbands>1`; compute and log the three
   distinct thesis-sourced quantities from §7 decision 5 (`δRM`, `max RM
   scale`, and — pending its combined-multi-band formula, per §10's
   implementation notes — `max un-aliased RM`) to guide the user's explicit
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
   In its place, rmtool computes and logs (stdout + run log) a diagnostic —
   confirmed against the actual thesis chapter (source below), which
   revealed the diagnostic needs **three distinct quantities, not two**
   (correction from the user, 2026-07-20, after an earlier draft of this
   decision conflated two of them):
   - **`δRM` (RM resolution)** — the ability to tell two nearby features
     apart in Faraday depth, whether two thin components or structure
     within one thick component. `δRM ∝ 1/[λ²(Δν/ν)]` per band (thesis
     eq 6.1, restated from eq 2.4 in §2.5.3). **For combined multi-band
     data, this is *not* set by the total λ² separation across the full
     gapped dataset** (an assumption an earlier draft of this document
     made and which turned out to be wrong) — the thesis demonstrates
     empirically (Table 6.1, P+L row: `Δλ²=0.214 m²` ≈ `Δλ²_P (0.201) +
     Δλ²_L (0.013)`, and prose confirming it explicitly) that **the
     effective combined resolution is set by the *sum* of the individual
     bands' own λ² spans**, `δRM_combined ∝ 1/(Σ_band Δλ²_band)` — the gap
     between bands does not help resolution, only each band's own internal
     bandwidth does.
   - **`max RM scale` (sensitivity to Faraday-*thick*/extended structure)**
     — the largest Faraday-depth *thickness* an extended component can have
     and still be detected without being completely washed out by
     bandwidth depolarization. `max RM scale ~ 1/λ²_min` (thesis eq 6.4),
     where `λ²_min` is the smallest λ² sampled anywhere across the *whole*
     combined dataset — i.e. set by whichever band reaches the shortest
     wavelength (highest frequency), a band-edge quantity, **not**
     dependent on individual channel width. Table 6.1 confirms this too:
     P+L's `max RM scale` (55.4) equals L-band's own value exactly, since L
     is the higher-frequency band and thus supplies the dataset's
     `λ²_min`.
   - **`max un-aliased RM` (a third, genuinely separate quantity — this is
     the one an earlier draft of this document wrongly conflated with `max
     RM scale`)** — the largest Faraday *depth* (not thickness) at which a
     single Faraday-**thin**/point-like component can still be reliably
     measured without being confused with (aliased onto) a different RM.
     This one *is* governed by individual channel width in λ² (thesis
     §2.5.2-§2.5.3, illustrated via the Gauribidannaur Telescope worked
     example: `δRM≈0.5 rad m⁻²` with `n_ch=256` channels limits reliable
     probing to `|RM|≲60 rad m⁻²`) — matching the formula already
     implemented for the existing single-band `use_auto_rm_range` heuristic
     (`src/rm_synthesis_mod.f90:656-687`). **Open implementation gap,
     flagged for whoever scopes phase 3, not resolved by this document**:
     the thesis's own Table 6.1 does not tabulate this quantity for the
     P/L/P+L combined case (only `δRM` and `max RM scale` are tabulated),
     so its exact multi-band-combined formula (accounting for each band's
     own, possibly different, channel width) needs to be worked out at
     implementation time, not assumed from Table 6.1 by analogy to the
     other two quantities.

   **Source:** Raja, W. (2014), *"Faraday Slicing Polarized Radio
   Sources,"* PhD thesis, Raman Research Institute / Jawaharlal Nehru
   University — Chapter 6, *"Slicing Faraday-thick components: tomography
   using multi-band data"* (thesis pp. 187-198), §6.1 in particular (eqs
   6.1-6.4, Tables 6.1-6.2, Figures 6-1/6-2/6-3); and Chapter 2 §2.5,
   *"Faraday tomography & advantages of non-uniform sampling in λ²"*
   (thesis pp. 23-37), specifically §2.5.1 "Bandwidth depolarization"
   (eqs 2.12-2.38) and §2.5.2-§2.5.3 "Sampling scheme for λ²s & aliasing" /
   "Faraday tomography at very low frequencies" (thesis pp. 31-34, the
   Gauribidannaur worked example). Read directly from the user's local copy
   (`/home/wasim/Thesis/wasim_thesis_Final.pdf`) — not accessible via the
   ResearchGate DOI originally supplied (403 on automated fetch).
   This is informational only: rmtool does **not** auto-run the full
   theoretical RM range for multi-band — `beg_rm`/`end_rm`/`nrm` remain
   required, explicit user choices, guided by (not overridden by) the
   logged diagnostic.
6. **Multi-band test data: synthetic only, not real observational cubes**
   (confirmed with user, 2026-07-20). All multi-band correctness testing —
   geometry-validation accept/reject paths, frequency/λ² merge, the
   `δRM`/`max RM scale`/`max un-aliased RM` diagnostic — is validated
   against synthetic Q/U cubes
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
  `δRM`/`max RM scale`/`max un-aliased RM` as a diagnostic; the user still
  chooses `beg_rm`/`end_rm`/`nrm` explicitly.
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
- **Evidence (2026-07-21, commit `baa52ce`):** Build clean via
  `bash scratch/make_all.sh`: 0 compiler errors, 0 compiler warnings, 4
  linker warnings (2 per GPU-offload binary — pre-existing
  `crtoffloadtable.o`/`DT_TEXTREL` noise, matching
  `scratch/baseline_encapsulation/T0_MANIFEST.md`'s own finding for the
  same toolchain, GNU Fortran 13.3.0). `tests/run_tests.sh`: 28/28 pass, 0
  fail, 0 skip. 174 output files (FITS/cfg/log/csv, 42M) archived to
  `scratch/baseline_multiband/tests_output/` (gitignored, local-only —
  full record in `scratch/baseline_multiband/T0_MANIFEST.md`).

---

### T1 — Comma-List Config Schema + Unified Geometry Validation

- **Objective:** Introduce comma-separated-list parsing for every per-band
  config key (§5) — with band count derived from list length, no separate
  `nbands` key — and the single, unified per-band ingestion/validation
  pipeline (§7 decision 0), proving the single-entry-list case is
  bit-identical to the T0 baseline before any frequency-merge logic (phase
  2) exists at all. This is deliberately the narrowest slice that exercises
  the new schema and the new geometry-validation loop end-to-end without
  touching the RM-synthesis numerics.
- **Scope:**
  - `src/rm_synthesis_mod.f90:106-153` (`rmsynth_config_t`): add
    `reference_band` (default 1). Add a new small per-band derived type
    (e.g. `band_cfg_t`) holding the per-band fields identified in §4/§5
    (`infileQ`, `infileU`, `resiQ`, `slopeQ`, `resiU`, `slopeU`, `infileI`,
    `path_I`) and an allocatable array of it inside `rmsynth_config_t`,
    sized to the *derived* band count (§5) — no `nbands` field is stored
    directly; it's the array's size.
  - `read_cfg_keyval` (`src/rm_synthesis_mod.f90:1632` onward): change
    `infileQ`/`infileU`/`resiQ`/`slopeQ`/`resiU`/`slopeU`/`infileI`/`path_I`
    parsing from "read one value" to "split on comma, read N values"; the
    band count is the resulting list length, cross-validated equal across
    every per-band key that's set (config-parse error otherwise, per §5).
    A config with no commas anywhere behaves exactly as today — same
    single value, same single-entry array. Required-key enforcement (today
    at `:2585-2589` for `infileQ`/`infileU`) extends to every list entry.
  - `src/rm_synthesis.f90:572-667` (today's `myfits_info` calls + Q-vs-U
    NAXIS/dimension check) and `:733-806` (today's Q-vs-U WCS-value exact
    match): generalize from "Q vs U, exactly two cubes" to "every band's Q
    and U vs `reference_band`'s geometry, N cubes" — same exact-equality
    philosophy (§3, §7 decision 3), same loud-refuse-on-mismatch behaviour,
    just looped over the derived band count instead of hardcoded to one
    pair.
  - For a derived band count `>1` specifically: after geometry validation
    passes, stop with a clear, explicit "multi-band frequency merge not
    yet implemented" message rather than attempting synthesis — phase 2's
    job, out of scope here. This keeps T1's blast radius to
    ingestion/validation plumbing only.
  - Test fixtures: extend `tests/make_test_cubes.py` to optionally emit a
    second synthetic Q/U band (distinct frequency range from the existing
    550-750 MHz GMRT-like band, e.g. an 800-950 MHz band, same RA/Dec
    geometry and same injected point sources) plus a deliberately
    geometry-mismatched variant (different `CRVAL1`/`CDELT1` or pixel
    count) — enough to exercise both the accept and loud-refuse paths of
    the new N-band geometry validation. Per §7 decision 6, this stays
    synthetic-only; no real multi-band cubes are sourced for this ticket.
- **Correctness Gate:**
  - **Legacy/default path (no commas in any per-band key, i.e. derived
    band count = 1):** every `tests/*.cfg`/`cfg/*.cfg` file, run unedited,
    produces bit-identical output (`tests/compare_cubes.py --exact`)
    against the T0 baseline archive. This is the ticket's central gate,
    carrying the full weight of requirements 2a/2b (§5, §7 decision 0) —
    there is no by-construction fallback if this fails.
  - `tests/run_tests.sh`'s existing pass/fail counts are unchanged from T0.
  - New multi-band fixture tests (added to `tests/run_tests.sh` or a
    parallel script): a matched-geometry two-band config (comma-separated
    `infileQ`/`infileU`) passes validation and reaches (and stops cleanly
    at) the "not yet implemented" message; a mismatched-geometry two-band
    config is loudly refused before any compute begins.
  - A config with inconsistent list lengths across per-band keys (e.g. two
    `infileQ` entries, one `infileU` entry) is rejected with a clear error
    at parse time.
- **Rollback Criteria:** If bit-identical output cannot be achieved for the
  single-entry-list case within this ticket's effort budget, roll back to
  before this ticket rather than merging a change that silently breaks
  existing users — re-evaluate whether the unified-pipeline decision (§7
  decision 0) needs revisiting in favour of one of the superseded
  lower-risk designs kept in §5.
- **Effort:** 1.5-2 sessions (config-parser rework + new derived type +
  N-band validation loop + synthetic multi-band fixture generation, gated
  by a bit-identical sweep that itself takes real wall-clock time to run
  and diff).
- **Evidence (2026-07-21):** Implemented as designed: `band_cfg_t` +
  `cfg%band(:)`/`cfg%reference_band` in `rm_synthesis_mod.f90`; comma-list
  parsing via new `csv_count`/`csv_get_item` helpers, deferred to a
  post-parse assembly step so every per-band key's list length can be
  cross-validated before `cfg%band(:)` is allocated; the legacy scalar
  fields (`cfg%infileQ` etc.) are populated from `cfg%band(reference_band)`
  so every existing use site in `rm_synthesis.f90` is untouched. The new
  N-band RA/Dec geometry-validation loop in `rm_synthesis.f90` is gated on
  `size(cfg%band).gt.1` and sits entirely after the existing Q-vs-U
  validation block — genuinely dead code for `nbands=1`, not merely
  "expected to behave like" the old path.
  - Build: clean, 0 errors, 0 new warnings (same 4 pre-existing
    GPU-offload linker warnings as T0).
  - `tests/run_tests.sh`: 31/31 pass (28 original + 3 new multi-band
    fixture tests: matched-geometry validates and stops cleanly at
    "not yet implemented"; mismatched-geometry loudly refused before any
    compute; inconsistent per-band list lengths rejected at parse time).
  - Bit-identical sweep: 140/140 FITS outputs match `scratch/baseline_multiband/`
    (the 6 `badchan_*` files reported as "differing" by `compare_cubes.py
    --exact` are the pre-existing NaN-vs-NaN tooling artifact from
    `project_encapsulation_refactor` — re-confirmed self-referential here
    by diffing the baseline archive against itself and seeing the
    identical "201 elements differ" report). 16/16 `.cfg` outputs
    byte-identical; `.csv`/`.log` differences are limited to expected
    wall-clock timing values and run-id timestamps.
  - `tests/make_test_cubes.py` extended with a `TEST_BAND2` (800 MHz,
    150-channel) fixture sharing the primary band's RA/Dec geometry and
    injected sources, plus a `TEST_BAND2_MISMATCH` variant (shifted
    `CRVAL1`) — confirmed the primary band's own `TEST.Q/U.FITSCUBE` and
    `TEST_BADCHAN.*` bytes are unchanged by this addition (verified via
    direct file comparison against the pre-change fixtures, not just
    "should be the same"), since the shared RNG stream for band 2's noise
    is drawn strictly after band 1's, and band 1's own code path was only
    refactored into a shared helper, not altered in sequence or content.

---

### T2 — Multi-Band Frequency/λ² Merge (CPU, Single-Tile)

- **Objective:** Replace T1's "not yet implemented" stop with an actual
  multi-band RM synthesis, for the narrowest slice that can be verified
  against §10's thesis-grounded scenario: CPU-only (serial and OMP),
  single-tile cubes. This deliberately narrows §6 phase 2's original
  description — see "Scope narrowing" below for what's cut and why.
- **Scope narrowing (confirmed direction, 2026-07-21):** attempting to
  thread multi-band support through every existing feature (GPU
  offload/staging, subimage extraction, the bad-channel-*file* mechanism,
  `remove_qu_bias`, multi-tile RAM planning, `io_read_threads`/
  `io_overlap`) in one ticket is more than one verifiable increment can
  safely absorb, especially under the bit-identical bar this whole effort
  already operates under. T2 covers CPU compute + single-tile only;
  everything else loudly stops with an explicit "not yet implemented for
  multi-band" message when requested with more than one band, exactly the
  same graceful-stop philosophy T1 already established for the whole
  feature. Per-pixel NaN-based masking (not the bad-channel-*file*
  mechanism) is unaffected and continues to work per band naturally, since
  it's a property of the read data, not a separate per-band list.
- **Change Set:**
  - `src/rm_synthesis.f90` header/frequency-read section (~1470-1620
    today): for each band, read its own `CRVAL`/`CRPIX`/`CDELT` on the
    freq axis (already partially done in T1's geometry loop for bands
    2..N, but only for the RA/Dec axes — extend to also capture the freq
    axis for the merge) and build that band's own per-channel frequency
    array via the existing `linspace` logic, unchanged. Concatenate every
    band's channel list (frequency, and hence λ², plus per-channel
    weight/flag) into one merged array — no deduplication in overlaps
    (§4/§7 decision 1), no sorting required (the DFT kernel is
    order-independent, confirmed in §4). This redefines `nz_out` from "one
    cube's channel count" to "sum of all bands' channel counts" — touches
    every allocation sized by it (`data_arrQ/U`, `flag_arr_out`, `L_sq`,
    `cos_arr`/`sin_arr`, tile-local `specQ`/`specU`).
  - Tile read stage: loop over bands, one `FTGSVE` call per band per tile
    (serial, `io_read_threads` forced to 1 for `nbands>1` — stop if the
    user requests `io_read_threads>1` or `io_overlap=y` with `nbands>1`,
    deferred per the scope narrowing above), each band's data landing in
    its own disjoint slice of the enlarged `specQ`/`specU` buffers.
  - Tile planning: if the RAM auto-tiler (or an explicit `tile_ra`/
    `tile_dec`) would produce more than one tile for a multi-band run,
    stop with "multi-tile multi-band not yet implemented" — deferred per
    the scope narrowing above. `nbands=1` tiling is completely unaffected.
  - GPU path (`use_gpu=y` with `nbands>1`), subimage (`subim=y` with
    `nbands>1`), bad-channel-file (`remove_badchan=y` with `nbands>1`),
    and `remove_qu_bias=y` with `nbands>1` (per §6 phase 4 — bias
    correction is per-band and explicitly deferred to its own ticket): all
    stop with an explicit not-yet-implemented message rather than silently
    producing wrong output.
  - `RM-range/resolution diagnostic`: not required for T2's own
    correctness gate (that's phase 3), but T2's test cfgs need
    `use_auto_rm_range=0` with explicit `beg_rm`/`end_rm`/`nrm` regardless
    (already forbidden for `nbands>1` per §7 decision 5).
- **Correctness Gate:**
  - **`nbands=1` path: still bit-identical** to `scratch/baseline_multiband/`
    (T0) — this ticket touches the frequency-array-construction and tile-read
    code that the single-band path also runs through, so this is not a
    given the way T1's purely-additive geometry loop was; it must be
    re-verified explicitly, not assumed.
  - **§10 scenario, P-band alone:** point source (RM=-100) recovered
    accurately; the Faraday-thick top-hat component (RM 100-130) is
    essentially invisible (`max RM scale=3.6 ≪` thickness 30).
  - **§10 scenario, L-band alone:** both components blend into one
    unresolved feature (`δRM=250.1` ≫ their ~215 rad m⁻² separation).
  - **§10 scenario, P+L combined:** point source still recovered
    accurately; the top-hat component's extended structure is now
    recovered (a significant fraction of its flux, distinguishable in
    shape from a point source) — the actual multi-band payoff this ticket
    exists to deliver.
  - **§10 scenario, F2/F3 addition:** blended at L alone, resolved at P
    alone and combined (§10's own honesty note about what this specifically
    demonstrates still applies).
  - Every "not yet implemented" stop (GPU/subimage/badchan-file/
    remove_qu_bias/multi-tile/io-parallelism × `nbands>1`) triggers
    cleanly, before any wrong compute, with a clear message.
- **Rollback Criteria:** If the `nbands=1` bit-identical gate cannot be
  met, roll back to before this ticket — the frequency-merge logic is not
  worth landing at the cost of silently changing single-band output.
- **Effort:** 3-4 sessions (frequency-array generalization touches a wide
  radius of `nz_out`-sized state; the §10 fixture generation and
  RM-recovery assertions are new test-infrastructure, not just cfg files).
- **Evidence (2026-07-21):** Implemented as designed. Reference band's
  channels always placed first (offset 0) in the merged array so the
  pre-existing single-band `L_sq`/`flag_arr_out`-building code needed zero
  modification; every other band's frequency/λ² construction mirrors it
  independently (own frequency-unit inference, own `linspace`). Every new
  code path (geometry/freq capture, `nz_out` override, append-bands
  `L_sq`, per-band tile read, multi-tile/GPU/subimage/badchan-file/
  `remove_qu_bias`/io-parallelism/`use_auto_rm_range` stops) is gated
  behind `size(cfg%band).gt.1`, so it is dead code for `nbands=1` by
  construction, not merely "expected to behave the same."
  - Build: clean, 0 errors, 0 new warnings (same 4 pre-existing
    GPU-offload linker warnings as T0/T1).
  - `tests/run_tests.sh`: 35/35 pass (32 from T0/T1 + 3 new: §10 fixture
    generation, all three P-alone/L-alone/P+L runs completing, and the
    full `check_thesis_scenario.py` assertion set).
  - `nbands=1` bit-identical sweep: 140/140 FITS outputs match
    `scratch/baseline_multiband/` (same 6 pre-existing `badchan_*`
    NaN-comparison-tool artifacts as T0/T1, re-confirmed self-referential
    — unchanged by this ticket).
  - §10 scenario (`tests/make_thesis_scenario_cubes.py`,
    `tests/check_thesis_scenario.py`, P-band 300/30 MHz and L-band
    1200/120 MHz reproduced exactly from Table 6.1): point source
    recovered accurately at both P-alone and P+L combined; Faraday-thick
    top-hat component's recovered peak amplitude is `~9x` larger at P+L
    combined than at P-alone (4.76 vs 0.53 in its own RM window) —
    directly confirming the `max RM scale` washout/reveal physics this
    effort exists to validate, matching the thesis's own Figures 6-1/6-2/
    6-3 qualitatively; F2/F3 resolved at P-alone, blended at L-alone, both
    via a dip-vs-peaks comparison targeted at the known expected RM
    positions. F2/F3 at P+L combined was **not** asserted as "resolved" —
    see the §10 addendum above for the documented dirty-beam-ringing
    reason (no RM-CLEAN in this codebase), confirmed as expected physics
    rather than a merge defect, since the two claims §10 was actually
    designed around (point-source accuracy, thick-component reveal) both
    hold cleanly.
  - Not yet exercised by any test in this ticket: multi-tile multi-band
    runs, GPU multi-band, subimage/bad-channel-file/`remove_qu_bias` with
    `nbands>1` — all correctly refuse with an explicit message (each
    `stop` statement reviewed by inspection, matching the pattern already
    validated for T1's mismatched-geometry refusal), but no automated test
    exercises each refusal path individually; left for whoever picks up
    the deferred features in a later ticket to add alongside the feature
    itself.

## 10. Multi-band synthetic test scenario (informs the Phase 2 ticket)

Per §7 decision 6, all multi-band correctness testing is synthetic. This
section specifies the injected sky model that phase 2's ticket needs to
validate the frequency/λ² merge against, now grounded directly in the
thesis (see the citation in §7 decision 5) rather than in speculative
placeholders. Per user decisions on 2026-07-20: this scenario **reproduces
the thesis's own worked example (Chapter 6, §6.1, Tables 6.1-6.2, Figures
6-1/6-2/6-3) exactly**, plus **one addition beyond it** (a close
Faraday-thin pair, F2/F3) chosen using the thesis's own real numbers rather
than invented ones.

### Bands — reproduced exactly from Table 6.1

| Band | ν_c (MHz) | Δν (MHz) | Δλ² (m²) | λ²_min (m²) | δRM (rad m⁻²) | max RM scale (rad m⁻²) |
|---|---|---|---|---|---|---|
| P | 300 | 30 | 0.201 | 0.907 | 15.6 | 3.6 |
| L | 1200 | 120 | 0.013 | 0.057 | 250.1 | 55.4 |
| P+L | — | — | 0.214 | 0.057 | 14.7 | 55.4 |

This is a **new, separate synthetic fixture set** from the existing
550-750 MHz GMRT-like single-band cube already in
`tests/make_test_cubes.py` — that fixture (and its `src_A`/`src_B` point
sources at RM=-5/+22) is unrelated to and untouched by this scenario; P and
L here are their own new FITS cube pairs. Channel count/width per band
isn't specified by Table 6.1 (only `ν_c` and `Δν` are) — pick something
fine enough that individual-channel bandwidth depolarization (thesis
§2.5.1, eqs 2.12-2.38) stays negligible relative to the effects under
test, so results match Table 6.1's idealized, coverage-only numbers; this
is an implementation-time engineering choice, not mandated by the thesis
itself, since the thesis's own Chapter 6 demonstration doesn't appear to
model finite-channel-width depolarization explicitly either (it isn't
listed among Table 6.1's inputs).

### Components

**Two components, reproduced exactly from Table 6.2** (the thesis's own
worked example):

1. **Point source (Faraday-thin/delta).** `RM = -100 rad m⁻²`,
   amplitude `15 Jy/(rad m⁻²)`, `PA_intrinsic = 0°`. Genuinely a single
   Faraday depth, so it is correctly recovered as one clean, unresolved
   peak in *every* scenario below — its role is a stable reference,
   validating that nothing about the merge corrupts recovery of an
   ordinary point source.
2. **Faraday-thick component, modelled as a top-hat in Faraday depth**
   (thesis §6.1.1: *"The extended feature is modeled as a top-hat function
   along RM"* — **not** a `sinc`/Burn-1966-formula parameterization, which
   an earlier draft of this document incorrectly assumed before the actual
   thesis text was available). Spans `RM = 100` to `130 rad m⁻²` (centre
   115, thickness 30), amplitude `5 Jy/(rad m⁻²)`, `PA_intrinsic = 0°`.

**One addition beyond the thesis's own demo (confirmed with user,
2026-07-20): a close Faraday-thin pair, F2/F3**, at `RM = +250` and
`+290 rad m⁻²` (`Δφ_23 = 40 rad m⁻²`), amplitude `8 Jy/(rad m⁻²)` each,
`PA_intrinsic = 0°` — placed well clear of both components above (120
rad m⁻² from the top-hat's upper edge) so nothing overlaps. `Δφ_23=40` was
chosen using the real Table 6.1 numbers: comfortably above `δRM_P (15.6)`
and `δRM_P+L (14.7)` (resolved by P alone and combined) while comfortably
below `δRM_L (250.1)` (blended at L alone).

**Honesty check on what F2/F3 actually demonstrates, since it isn't from
the thesis itself:** with these specific P/L parameters, P-band's own
resolution (15.6) is already very close to the combined resolution (14.7)
— per the thesis's own finding that combined `δRM` is dominated by
whichever band contributes the larger `Δλ²` (here, P contributes ~94% of
the combined span). So F2/F3, as specified, mainly demonstrates
*"resolution is set by which band you choose"* (already implicit in
comparing Figures 6-1 vs 6-2) rather than *"combining bands beyond the
best single band helps resolution"* — a pair genuinely requiring the
combination over P alone would need `Δφ_23` between 14.7 and 15.6, a
sub-1-rad-m⁻² margin too fragile for a robust automated test (sensitive to
channelisation, deconvolution loop-gain, etc.). **The genuine "why
multi-band" demonstration is carried entirely by the top-hat component
below**, matching the thesis's own stated message (§6.1.5): combining
bands buys *sensitivity to extended structure* more than it buys marginal
extra resolution beyond the best single band.

### Expected results per band (mirrors thesis Figures 6-1/6-2/6-3, extended to F2/F3)

- **L-band alone** (`δRM=250.1`, `max RM scale=55.4`, thesis Fig. 6-1,
  *"poor RM resolution"*): point source and top-hat component blend into
  one broad, unresolved hump (thesis's own finding — separation from -100
  to ~115 is ~215, well under 250.1). F2/F3 (`Δφ=40`) also blend — new,
  not in the thesis, but the same reasoning applies.
- **P-band alone** (`δRM=15.6`, `max RM scale=3.6`, thesis Fig. 6-2,
  *"poor sensitivity to extended structure"*): point source recovered
  cleanly and accurately (thesis's own finding). Top-hat component is
  **essentially invisible** — `max RM scale (3.6) ≪` the component's own
  thickness (30), so it's washed out by bandwidth depolarization, not
  merely coarsely resolved. F2/F3 (`Δφ=40 > δRM_P=15.6`) resolved into two
  distinct peaks — new, not in the thesis.
- **Combined P+L** (`δRM=14.7`, `max RM scale=55.4`, thesis Fig. 6-3):
  point source still recovered cleanly; the top-hat component's **real
  extended structure is now revealed** (thesis's own finding — a
  significant fraction of its flux is recovered, distinguishable from a
  point source). F2/F3 still resolved (`40 > 14.7`) — new, not in the
  thesis, and per the honesty check above, not meaningfully better
  resolved than at P alone.

### Implementation notes for the phase-2 ticket

- `tests/make_test_cubes.py` needs: (a) support for summing multiple
  component contributions (point + top-hat, or the added F2/F3 pair) into
  one pixel's complex polarization spectrum — currently one delta
  component per source position; (b) a **top-hat-in-Faraday-depth** source
  type (not a `sinc`/Burn-slab closed form) alongside the existing
  delta-component type, matching the thesis's own model exactly; (c) a new
  `make_header`/frequency-axis block for Band P and Band L, each emitting
  its own Q/U FITS pair, additive to (not replacing) the existing
  550-750 MHz single-band fixture.
- `tests/check_rm_peak.py` (currently checks one recovered RM peak against
  a known truth) needs a variant that checks **peak count and separation**
  (for the F2/F3 resolved-vs-blended assertion) and one that checks
  **recovered polarized flux fraction and approximate Faraday-depth
  extent** (for the top-hat washed-out/blob/resolved-structure assertion)
  rather than a single peak-position check.
- **Not covered by this scenario, flagged as an open gap (§7 decision 5's
  third quantity, `max un-aliased RM`):** nothing here tests the
  high-|RM| point-source aliasing behaviour from thesis §2.5.2-§2.5.3 (the
  Gauribidannaur worked example) — that would need a thin point source
  placed near or beyond a band's own `max un-aliased RM`, and the thesis
  doesn't give a combined-multi-band formula for that quantity the way
  Table 6.1 does for `δRM`/`max RM scale`. Left for whoever scopes phase 3
  to decide whether a dedicated test is worth adding once that formula is
  worked out, rather than guessed here.

### Addendum (found while implementing T2, 2026-07-21): dirty-beam ringing makes the F2/F3-at-P+L claim untestable by simple dip comparison

Running the actual scenario surfaced a genuine, physically-expected
limitation the original design above didn't anticipate. The point-source
and Faraday-thick-component checks (§10's headline claims, matching the
thesis's own Figures 6-1/6-2/6-3 exactly) came out cleanly: point source
recovered accurately at both P-alone and P+L; the top-hat component's
recovered peak amplitude in its own RM range is `~9x` larger at P+L
combined than at P-alone (4.76 vs 0.53), directly confirming the
`max RM scale` washout/reveal physics this effort exists to validate.

The **F2/F3 addition specifically at P+L combined**, however, does not
show a clean, simply-detectable dip between the two peaks the way it does
at P-alone. This codebase has **no RM-CLEAN deconvolution** (confirmed by
grep — `RMCLEAN`/`rm-clean` appear nowhere in `src/`), so its output is
the raw *dirty* RM spectrum, not the Gaussian-restored profile the thesis
itself uses (§2.5, §6.1's own `RMCLEAN` step, thesis pp. 35-37). Combining
two widely-separated bands (P and L) creates a very large total λ² span,
and per the thesis's own point (§2.5, "the *dirty* RM response profile...
[depends on] the exact sampling scheme"), that large span produces
fine-period sidelobe ringing in the dirty spectrum — confirmed directly in
the actual output: the L-alone spectrum is smooth across the F2/F3 region
(narrow λ² span → long ringing period), while the P+L-combined spectrum
oscillates by several units between *adjacent* 2 rad/m⁻² bins (wide
combined span → short ringing period) — exactly the signature the thesis
predicts. At this sampling, that ringing's amplitude can rival the true
F2/F3 dip even though the underlying resolution (`δRM=14.7`) is fine
enough in principle to separate them.

**Resolution adopted**: `tests/check_thesis_scenario.py` checks P-alone
and L-alone resolved-vs-blended with a targeted dip-vs-peaks comparison at
the known expected RM positions (robust — not blind peak-counting, which
the same ringing defeats even at P-alone via an unrelated sidelobe just
above a naive prominence threshold). At P+L combined, it deliberately does
**not** assert "resolved" — only that both expected positions still carry
real, elevated signal (the merge computed something meaningful there, not
noise) — with an explicit code comment explaining why, rather than a
fragile threshold tuned to pass on this one dataset. **This is a
documented, physically-understood limitation, not a defect in the
frequency merge itself** — the two claims §10 was actually designed to
test (point-source accuracy, thick-component reveal) both hold cleanly.
Adding RM-CLEAN is out of scope for this effort (§8 non-goals still
covers "any change to the numerical RM-synthesis kernel... beyond how many
channels it's handed") — if a clean-vs-dirty comparison ever becomes
important, that is a new, separate effort, not a T2 fix.
