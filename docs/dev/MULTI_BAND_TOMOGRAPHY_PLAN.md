# Multi-Band Faraday Tomography — Feasibility & Design Plan

Branch: `multi-band-tomography` (from `develop`)

**Status: T0-T14 all merged to `develop` and tagged `5.0-rc.1` (real-scale
validation pending before an actual `main` release -- see
`docs/dev/ARCHIVED/CHANGELOG.md`'s `[5.0]` entry and `docs/dev/ARCHIVED/RELEASE_NOTES_5.0.md`). T13
(`match_cubes`) and T14 (`CASAMBM`/`BEAMS` propagation for the cross-band
toolchain) were implemented and verified on `multi-band-tomography` and
merged into `develop` via `5e2291f` -- confirmed directly against git
history (`git merge-base --is-ancestor 5e2291f develop`), not assumed;
this line previously said "not yet merged," which was stale.**

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

Reference: `docs/user/PARALLELISM.md`, `docs/user/ARCHITECTURE.md`.

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
change in observable behaviour" bar `docs/dev/ENCAPSULATION_REFACTOR_PLAN.md`
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
`docs/dev/ENCAPSULATION_REFACTOR_PLAN.md` / `IO_PARALLEL_OPTIMISATION_PLAN.md`
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
     **Resolved (confirmed with user, 2026-07-22): eq 2.4 (thesis p. 15,
     found on a second, more targeted read — not in the pages originally
     consulted for decision 5), the exact source formula behind the
     Gauribidannaur illustration**:
     ```
     δRM  ~  π/(2λc²) · (νc/Δν)            -- per-band resolution
     ΔRM  ~  n_ch × δRM                     -- per-band un-aliased span
     ```
     (`δ` = resolution, `Δ` = un-aliased span, `λc = c/νc`). Verified
     directly against real numbers: `δRM_P=15.73`, `δRM_L=251.68` (Table
     6.1: 15.6/250.1, within ~1%); Gauribidannaur (`νc=34.5 MHz,
     Δν=1.5 MHz, n_ch=256`) gives `ΔRM/2=61.2`, matching the thesis's own
     quoted *"maximum absolute RM of only about 60 rad m⁻²"* almost
     exactly — confirming `ΔRM` is the *full* span (`−ΔRM/2` to `+ΔRM/2`),
     not the max `|RM|` itself. Algebraically equivalent (linear-bandwidth
     approximation) to the existing single-band `use_auto_rm_range`
     heuristic's `d_nu = fac/t_span` (`src/rm_synthesis_mod.f90:656-687`,
     now documented in-code, `f8a89eb`) when `fac=π` (this codebase's
     default) — the same physics, reparameterized.

     **Decision on the multi-band *combined* value (confirmed with user,
     2026-07-22): NOT computed by this ticket.** Per-band `ΔRM` is logged
     for each band individually. The user's own words: *"I need to think
     what the overall max RM would be when we combine the bands. In FFT
     where all L_sq channels are equal, this is easy. But in DFT the
     un-aliased span should be more than the FFT case."* — i.e. the naive
     multi-band generalization (by analogy with `δRM`/`max RM scale`'s
     already-established combined forms) is *not* simply summing or
     min/max-ing the per-band values, because this codebase's DFT-based
     extraction (§4: order-independent, no FFT/uniform-sampling
     requirement) should, in principle, tolerate a *larger* un-aliased
     span than an equal-spacing FFT treatment would — a genuine open
     research question, not an implementation gap to fill in by analogy.
     Deferred to a future ticket once resolved, rather than guessed here.

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
  `docs/dev/ENCAPSULATION_REFACTOR_PLAN.md`.

## 9. Tickets

Ticket format follows this repo's existing convention
(`docs/dev/ENCAPSULATION_REFACTOR_PLAN.md`, `docs/dev/IO_PARALLEL_OPTIMISATION_PLAN.md`):
Objective / Scope / Change Set / Correctness Gate / Rollback Criteria /
Effort, each getting an **Evidence (...)** section appended once done.

---

### T0 — Baseline Lock

- **Objective:** Freeze an exact, reproducible reference of current
  `develop` behaviour before any multi-band code changes land, so every
  later ticket's bit-identical correctness gate (§5, §7 decision 0) has a
  concrete baseline to diff against — the same purpose T0 served in
  `docs/dev/ENCAPSULATION_REFACTOR_PLAN.md`, adapted for this branch. This
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

---

### T3 — RM-Range/Resolution Diagnostic (`δRM`, `max RM scale`, per-band `ΔRM`)

- **Objective:** Give multi-band users the guidance `use_auto_rm_range=1`
  would otherwise have provided (forbidden for `nbands>1` since T1/T2) by
  computing and logging (stdout + run log) the three §7-decision-5
  quantities: `δRM` per band and combined, `max RM scale` per band and
  combined, and `ΔRM` (un-aliased span) per band only — informational,
  does not gate or auto-select `beg_rm`/`end_rm`/`nrm`, which stay
  required, explicit user choices.
- **Scope:**
  - Pure diagnostic — read-only with respect to every array T2 built
    (`band_czval`/`band_czpix`/`band_zinc`/`band_nz`, already fully
    populated for every band, including the reference band, by the end of
    T1/T2's geometry-validation block in `rm_synthesis.f90`). Touches
    `nz_out`, `L_sq`, `specQ`/`specU`, and every other array T2's
    correctness gate covers **not at all** — this ticket only adds
    `write(*,*)`/`log_message` calls, nothing numerically consumed
    downstream.
  - Per band `i`: compute edge frequencies from `band_czval(i)`/
    `band_czpix(i)`/`band_zinc(i)`/`band_nz(i)` (same `z1`/`zn`
    construction already used for the append-other-bands `L_sq` loop, half
    a channel extended on each edge, mirroring the existing single-band
    `Lsq1`/`Lsq2` edge treatment in `extract_general_setup`), then:
    - `Δλ²_i` = λ² at the low-frequency edge − λ² at the high-frequency
      edge.
    - `δRM_i = cfg%fac / Δλ²_i` (§7 decision 5, thesis eq 2.4 in its
      edge-based form; `cfg%fac` defaults to `π`, matching eq 2.4's
      literal `π/(2λc²)·(νc/Δν)` to the linear-bandwidth approximation).
    - `max_RM_scale_i = cfg%fac / λ²_min_i` (thesis eq 6.4), `λ²_min_i` =
      λ² at band `i`'s own high-frequency edge.
    - `ΔRM_i = band_nz(i) × δRM_i` (thesis eq 2.4's un-aliased-span form);
      log both the full span and `ΔRM_i/2` (max `|RM|`) — the thesis's own
      Gauribidannaur example quotes the latter ("maximum absolute RM"), so
      log both to avoid the ambiguity a future reader could otherwise hit.
  - Combined (all bands together): `δRM_combined = cfg%fac /
    (Σ_i Δλ²_i)`; `max_RM_scale_combined = cfg%fac / (min_i λ²_min_i)`
    (both confirmed against Table 6.1's P+L row in §7 decision 5). **No
    combined `ΔRM`** — explicitly deferred per the 2026-07-22 decision
    recorded in §7 decision 5; log only the per-band values for this
    quantity, with a one-line note in the output explaining why no
    combined figure is shown.
  - Gated behind `n_bands_t2.gt.1` — a no-op for `nbands=1`, where the
    existing single-band `use_auto_rm_range` heuristic already covers this
    ground and remains untouched.
- **Correctness Gate:**
  - `nbands=1` bit-identical sweep unaffected (trivially — this ticket's
    code doesn't execute for `nbands=1` at all).
  - For the §10 scenario (P: 300/30 MHz, L: 1200/120 MHz — real,
    already-known-correct numbers to check against): logged `δRM_P`,
    `δRM_L`, `δRM_combined`, `max_RM_scale_P`, `max_RM_scale_L`,
    `max_RM_scale_combined` match Table 6.1 (15.6/250.1/14.7 and
    3.6/55.4/55.4 respectively) to within the same ~1% tolerance the
    hand-verification above already showed; per-band `ΔRM_P`/`ΔRM_L`
    computed and logged (no combined-`ΔRM` line present).
  - Output FITS files bit-identical to what T2 alone would have produced
    for the same inputs (this ticket adds no new numerical computation
    that reaches any output array).
- **Rollback Criteria:** N/A in the usual sense (no correctness gate can
  meaningfully fail here beyond a formula transcription error, which the
  Table 6.1 cross-check above exists to catch) — if the logged numbers
  don't match Table 6.1 within tolerance, fix the formula before merging
  rather than shipping a diagnostic that misleads users.
- **Effort:** 0.5-1 session (pure diagnostic, no architecture change, no
  new test infrastructure beyond checking logged numbers against Table 6.1
  by hand or a small script).
- **Evidence (2026-07-22):** Implemented as designed, gated behind
  `n_bands_t2.gt.1`, inserted right after the existing `ngood_chan`
  book-keeping loop in `rm_synthesis.f90` — reads only
  `band_czval`/`band_czpix`/`band_zinc`/`band_nz` (already fully
  populated by T1/T2) and `cfg%fac`; writes nothing to any array. Build
  clean, 0 errors, 0 new warnings (same 4 pre-existing GPU-offload linker
  warnings). `tests/run_tests.sh`: 35/35 pass, unchanged from T2. `nbands=1`
  bit-identical sweep: 140/140 FITS match `scratch/baseline_multiband/`
  (same 6 pre-existing NaN-artifact diffs) — expected, since this ticket's
  code doesn't execute at all for `nbands=1`.

  Ran against the §10 scenario (P: 300/30 MHz, L: 1200/120 MHz) and
  checked the logged numbers directly against Table 6.1:

  | | `δRM` (rad m⁻²) | `max RM scale` (rad m⁻²) |
  |---|---|---|
  | P (logged / Table 6.1) | 15.651 / 15.6 | 3.468 / 3.6 |
  | L (logged / Table 6.1) | 250.419 / 250.1 | 55.494 / 55.4 |
  | Combined (logged / Table 6.1) | 14.731 / **14.7** | 55.49 / 55.4 |

  All within ~1% or better (the combined `δRM` essentially exact) — the
  code's exact per-channel-edge computation tracks Table 6.1 more tightly
  than the linear-λc hand-verification done when the formula was resolved
  (§7 decision 5). Per-band `ΔRM` (un-aliased span) logged for both bands
  with no combined figure, per the 2026-07-22 decision — confirmed the
  "no combined value" note prints correctly alongside the per-band ones.

---

### T4 — Multi-Tile Multi-Band Runs

- **Objective:** Remove T2's single-tile restriction, reopening the
  thesis's own actual motivating use case (§6.2: combining Arecibo P+L
  data over the GALFACTS survey field to study diffuse Galactic polarized
  emission — inherently wide-field, not reachable under a single-tile
  ceiling) rather than leaving multi-band tomography usable only on
  images small enough to fit one RAM-budgeted tile.
- **Why this is a small, contained change (confirmed by direct code
  inspection before scoping, 2026-07-22):** T2's per-band tile-read code
  already sits inside the existing RA/Dec tile loop and already reuses
  the tile loop's own live `fpixels`/`lpixels` for its RA/Dec bounds —
  the same mechanism the reference band's own read already relies on
  across many tiles in existing single-band production runs. `specQ`/
  `specU` (sized once, `tile_ra×tile_dec×nz_out`) and the merged
  `L_sq`/`cos_arr`/`sin_arr`/`flag_arr_out` (built once, tile-invariant)
  are already correctly tile-agnostic. Per-band FITS units are opened
  once and reused across tiles exactly like the reference band's own
  21/22. `io_overlap` stays separately blocked for multi-band regardless
  (T2's own restriction, untouched here) — no double-buffering
  interaction to introduce. The only genuinely single-tile-specific code
  is the one blocking check itself (`rm_synthesis.f90:2105-2106`, right
  after `plan_tile`).
- **Scope:**
  - Remove (or relax to a no-op) the `cfg%tile_ra.lt.nx_out .or.
    cfg%tile_dec.lt.ny_out` stop for `nbands>1`.
  - No other source change anticipated — this ticket's real weight is in
    the **Correctness Gate** below, not new code, since the read-stage
    generalization was already written tile-agnostically in T2.
  - GPU multi-band, subimage-for-multi-band, bad-channel-file-for-
    multi-band, `remove_qu_bias`-for-multi-band, and
    `io_read_threads>1`/`io_overlap`-for-multi-band remain out of scope,
    unaffected by this ticket — their own `stop` checks are untouched.
- **Correctness Gate:**
  - `nbands=1` bit-identical sweep unaffected (this check only ever fired
    for `n_bands_t2.gt.1`, so removing it cannot change single-band
    behaviour).
  - **New invariant this ticket exists to prove**: a multi-tile multi-band
    run produces **bit-identical output** to the single-tile run of the
    *same* multi-band data — tiling must not change the scientific
    answer, the same bar the original single-band tiling implementation
    was already held to. Verify by forcing a small `tile_ra`/`tile_dec`
    (via explicit cfg values, not `tile_auto`) on the §10 P+L scenario so
    the run spans multiple tiles, and comparing output FITS against the
    existing single-tile §10 run via `compare_cubes.py --exact`.
  - §10 scenario's own scientific assertions (point-source accuracy,
    thick-component reveal, F2/F3 behaviour) still hold under the
    multi-tile run — re-run `check_thesis_scenario.py` against the
    multi-tile output, not just the tile-count-invariance diff above.
- **Rollback Criteria:** If multi-tile output is not bit-identical to
  single-tile for the same data, do not ship this ticket — investigate
  whether some tile-loop state actually was implicitly single-tile-only
  despite the inspection above (e.g. a variable only initialised before
  the tile loop that should instead be per-tile), rather than assume the
  inspection was complete.
- **Effort:** 0.5-1 session (the change itself is small; most of the
  effort is the multi-tile-vs-single-tile bit-identical verification).
- **Evidence (2026-07-22):** Confirmed by direct code inspection before
  writing any change (user's own architectural reasoning, verified rather
  than assumed): T2's per-band tile-read code already sits inside the
  RA/Dec tile loop (`rm_synthesis.f90:3164` loop, per-band read at
  `:3315`) and already reuses the tile loop's live `fpixels`/`lpixels` —
  confirmed these are set to the *current* tile's bounds earlier in the
  same iteration, exactly like the reference band's own read. The single
  blocking check (`:2105-2106`) was removed outright (not relaxed),
  replaced with an informational tile-count log line.

  Verified the new invariant directly: forced `tile_ra=tile_dec=16` (4
  tiles on the §10 scenario's 32×32 image, `tile_auto=n`) and confirmed
  all 4 output products (AMP/PHA/MASK/NVALID) are **bit-identical** to the
  single-tile P+L run — `tests/compare_cubes.py --exact` on all four,
  wired into `run_tests.sh` as a permanent regression check (T4's own
  correctness gate, not just a one-off manual check). One real bug caught
  by the automated version that the manual check missed: the generated
  multi-tile cfg initially reused the single-tile run's own `outfile`
  path unchanged, colliding with its still-present output files (the
  binary correctly refused to overwrite) — fixed by explicitly
  substituting the `outfile` line rather than copying it verbatim.

  Build clean (0 errors, 0 new warnings, same 4 pre-existing GPU-offload
  linker warnings). `tests/run_tests.sh`: 37/37 pass (35 from T0-T3 + 2
  new: 4-tile-run confirmation, bit-identical-to-single-tile). `nbands=1`
  bit-identical sweep: 140/140 FITS match `scratch/baseline_multiband/`
  (same 6 pre-existing NaN-artifact diffs) — expected, since the removed
  check only ever fired for `nbands>1`.

---

### T5 — Split-Band Identity Test (contiguous split == undivided cube)

- **Objective:** Add the single most direct, unambiguous regression check
  for the multi-band merge mechanism (confirmed with user, 2026-07-22, as
  a distinct and *more fundamental* test than the §10 scientific
  scenario): split one existing single-band test cube's channel axis into
  two **contiguous** halves (no gap — channels 1..N/2 as "band A", N/2+1..N
  as "band B"), feed them back through rmtool as a 2-band multi-band run,
  and check the result against the *same* cube run undivided, single-band.
  Unlike §10 (which validates real multi-band physics against thesis
  numbers, but needs qualitative interpretation for some assertions —
  e.g. the F2/F3-at-P+L ringing caveat), this test has a crisp, purely
  mechanical pass/fail: for a contiguous split with the low-frequency half
  as `reference_band`, the merged channel sequence reconstructs *exactly*
  the original cube's own channel order, so the expected result is not
  merely "close" but **bit-identical** — any future change to the
  frequency-merge/tile-read architecture that breaks this silently would
  be caught immediately, without needing to reason about whether an
  observed difference is expected science or a regression.
- **Scope:**
  - `tests/make_test_cubes.py`: after building the existing primary-band
    `q_cube`/`u_cube` arrays (unchanged), slice them (not regenerate) into
    `TEST_SPLIT_LO.Q/U.FITSCUBE` (channels 1..100) and
    `TEST_SPLIT_HI.Q/U.FITSCUBE` (channels 101..200), each with `CRVAL3`
    adjusted to that half's own first-channel frequency and `NAXIS3=100`,
    `CDELT3` unchanged. Slicing the already-built array (rather than
    re-synthesizing two smaller cubes independently) guarantees
    pixel-for-pixel correspondence to the undivided cube by construction,
    not by re-deriving the same noise/signal separately.
  - `tests/run_tests.sh`: new section running `infileQ =
    TEST_SPLIT_LO.Q.FITSCUBE,TEST_SPLIT_HI.Q.FITSCUBE` (same for
    `infileU`) with the *same* numeric parameters (`beg_rm`/`end_rm`/`nrm`/
    `fac`/`ofac`/etc.) already used by section 5's existing single-band
    `serial` run, so the comparison target is that run's own
    already-produced `tests/output/serial.*.FITS` — no separate "run the
    undivided cube again" step needed.
- **Correctness Gate:**
  - AMP/PHA/MASK/NVALID from the split-band run **bit-identical**
    (`compare_cubes.py --exact`) to `tests/output/serial.*.FITS`. If this
    turns out to be merely close rather than exact, that itself is a
    finding worth recording precisely (which floating-point step
    introduced the difference), not silently loosening the check to
    `--rtol`.
  - `nbands=1` bit-identical sweep and all prior tests unaffected (this is
    a new, additive test; no source change anticipated unless the gate
    above surfaces a real bug).
- **Rollback Criteria:** N/A for the test itself. If the gate fails, that
  is the signal to investigate — do not weaken the test to make it pass.
- **Effort:** 0.5 session (reuses existing fixtures and an existing
  comparison target; no new source code expected, only new test
  fixtures/cfg/assertions).
- **Evidence (2026-07-22):** Confirmed no source change was needed — the
  gate passed on the first run. `tests/make_test_cubes.py` slices the
  already-built primary-band arrays into `TEST_SPLIT_LO`/`TEST_SPLIT_HI`
  (100 channels each, `CRVAL3` adjusted per half, `CDELT3` unchanged);
  confirmed by direct file comparison that every existing fixture
  (`TEST.Q/U`, `TEST_BADCHAN.*`, `TEST_BAND2.*`) is byte-identical to
  before this addition. `run_tests.sh` section 17 runs the 2-band split
  config and diffs AMP/PHA/MASK/NVALID against section 5's existing
  `serial.*.FITS` — **bit-identical on the first attempt**, exactly as
  the architectural reasoning predicted (contiguous split with the
  low-frequency half as `reference_band` reconstructs the exact original
  channel sequence, so the DFT sum performs the identical floating-point
  operations in the identical order). Build unchanged (no source edited
  this ticket; still 0 errors, 0 new warnings). `tests/run_tests.sh`:
  38/38 pass (37 from T0-T4 + this one). `nbands=1` bit-identical sweep
  unaffected (140/140 FITS, 6 expected pre-existing artifacts).

### T6 — Per-Band Channel Sub-Range Selection (`subim` channel axis only)

- **Objective:** Investigation of the four remaining "not yet implemented"
  T2 restrictions (`rem_mean`, `subim`, `remove_badchan`,
  `remove_qu_bias`/`resiQ`/`slopeQ`/`resiU`/`slopeU`, confirmed with user
  2026-07-22) found: `rem_mean` needs no change (pure per-pixel,
  all-channels-loaded computation, confirmed correct as-is);
  `remove_qu_bias`/`resiQ`/`slopeQ`/`resiU`/`slopeU` are dead code even in
  the pre-existing single-band tool -- `cfg%resiQ` etc. are parsed, stored,
  and printed, but never applied to `specQ`/`specU` anywhere in the
  compute pipeline, and the I-cube read into `stI` is never passed to
  `tile_extract_gpu_rm_blocked`. Confirmed with user this is intentional:
  the I-cube read path is being kept alive deliberately as a placeholder
  for future Q/U-vs-I calibration (normalising by the I-spectrum where
  Q/U show non-smooth behaviour across band boundaries), not a bug to fix
  under this project. `remove_badchan` is a genuine gap needing new
  per-band plumbing (own ticket, not in scope here). `subim` splits into
  two cases: RA/Dec-only sub-windowing already works with zero changes
  (shared `fpixels(1)/(2)` apply identically to every band's read, since
  bands are geometry-matched by policy); channel-range sub-windowing
  (`subim_chan_blc/trc/inc`) does not -- confirmed by user as a wanted
  feature in its own right (reject bad edge channels or hand-pick a good
  sub-range independently per band), not merely a restriction to lift.
- **Scope:**
  - `rm_synthesis_mod.f90`: `band_cfg_t` gains `chan_blc`, `chan_trc`,
    `chan_inc` (integers, defaults 0/0/1, same semantics as today's
    `subim_chan_blc/trc/inc`), parsed via the same raw-comma-buffer +
    deferred-assembly pattern as `resiQ`/`slopeQ`/`resiU`/`slopeU`
    (optional keys: only band-count-validated/parsed if the key was seen
    at all, else every band defaults to 0/0/1 -- matching today's
    "key absent" behaviour exactly). Legacy scalars
    `cfg%subim_chan_blc/trc/inc` populate from `cfg%band(reference_band)`,
    so every existing single-band cfg file is unaffected.
  - `rm_synthesis.f90`: new per-band scratch arrays `band_chan_blc(:)`,
    `band_chan_inc(:)` (parallel to `band_czval`/`band_nz`/etc). At the
    point each band's raw NAXIS3 becomes known (reference band ~line 947,
    each other band ~line 1046), resolve that band's own
    `chan_blc/trc/inc` (0 defaults to "full band", mirroring the existing
    single-band `subim` resolution at lines 1645-1679) into
    `band_chan_blc(iband)`/`band_chan_inc(iband)` and a **corrected**
    `band_nz(iband)` = that band's own *selected* channel count (not its
    raw NAXIS3) -- with an explicit bounds check that
    `chan_blc`/`chan_trc` don't exceed that band's actual NAXIS3, stopping
    with a per-band error if they do (mirroring the existing single-band
    bounds check at lines 1705-1718). This also fixes a latent bug the
    investigation surfaced: `band_nz(cfg%reference_band)` was being set
    to the *raw* `naxes(freq_axis)` before `fpixels`/`lpixels` even get
    built from `subim_chan_blc/trc`, so `nz_out_band`/`band_offset`
    (computed from `band_nz` before the existing `subim` block runs) would
    have silently used the wrong, oversized count for the reference band
    the moment channel-range subimaging was allowed for multi-band --
    unreachable today only because `subim` is currently blocked outright
    for `nbands>1`.
  - The T2 "append other bands' L_sq/flag_arr_out" loop (~1875-1905) and
    the T3 diagnostic loop (~1920+): both currently assume a band's first
    contributing channel is raw pixel 1 (`z1_band` derived straight from
    that band's `CRVAL`/`CRPIX`/`CDELT`). Both need `z1_band` shifted by
    `(band_chan_blc(iband)-1)*band_zinc(iband)` and the `zn_band` step
    multiplied by `band_chan_inc(iband)`, so frequency/lambda^2 values are
    computed for the *selected* channels, not always starting at channel 1.
  - The per-band tile-read `FTGSVE` calls (~3326-3354): currently hardcode
    `io_par_fpixels(freq_axis)=1, io_par_lpixels(freq_axis)=band_nz(iband)`
    (always the band's full range) with the shared `incs` array. Needs a
    per-band `fpixels(freq_axis)=band_chan_blc(iband)`,
    `lpixels(freq_axis)=band_chan_blc(iband)+(band_nz(iband)-1)*band_chan_inc(iband)`,
    and a per-band copy of `incs` with `incs(freq_axis)=band_chan_inc(iband)`
    (RA/Dec entries unchanged, shared across all bands as today).
  - Remove the blanket `cfg%subim` stop-check for `nbands>1` (lines
    903-907) entirely -- no replacement guard needed, since channel-range
    subimaging is now band-aware and RA/Dec-only subimaging already worked.
- **Correctness Gate:**
  - `nbands=1` bit-identical sweep, unaffected.
  - New test reusing T5's fixtures with zero new cube generation: 2-band
    run with band 1 = `TEST_SPLIT_LO` (unrestricted) and band 2 = the
    *full* undivided `TEST.Q/U` cube (200 ch) with `chan_blc,chan_inc =
    <blc>,1` restricting it, per-band, down to exactly channels 101-200 --
    i.e. `TEST_SPLIT_HI`'s own range. Must be bit-identical to T5's
    already-passing split-band result (`tests/output/serial.*.FITS`).
    Only passes if the per-band offset/count/z1-shift arithmetic above is
    exactly right, not merely "doesn't crash" -- this is a real exercise
    of the new code, not a restatement of T5.
- **Rollback Criteria:** revert to the blanket `subim` block for `nbands>1`
  if the identity test doesn't come out bit-identical (investigate and
  fix, don't loosen to a tolerance).
- **Effort:** ~1 session (schema + two arithmetic-correction sites + a
  per-band read-loop parameterization + one new test).
- **Evidence (2026-07-22):** Implemented as scoped. `band_cfg_t` gained
  `chan_blc`/`chan_trc`/`chan_inc`; legacy scalars
  `cfg%subim_chan_blc/trc/inc` now populate from
  `cfg%band(reference_band)`, so single-band cfg files are unaffected.
  Fixed the latent bug the investigation surfaced: `band_nz` for both the
  reference band and every other band is now resolved from that band's own
  `chan_blc/trc/inc` (with an explicit per-band NAXIS bounds check) at the
  point each band's header is read, instead of the raw NAXIS3 -- this also
  corrects `nz_out_band`/`band_offset`, which are accumulated from
  `band_nz` before the reference band's own (pre-existing) `subim` block
  runs. `z1_band`/`zn_band` in both the T2 L_sq-append loop and the T3
  diagnostic loop now shift by `(band_chan_blc(iband)-1)*band_zinc(iband)`
  and step by `band_chan_inc(iband)`, so per-band frequency/lambda^2 values
  reflect the *selected* channels, not always channel 1. The per-band
  tile-read `FTGSVE` calls now use `band_chan_blc(iband)` as
  `fpixels(freq_axis)`, `band_chan_blc(iband)+(band_nz(iband)-1)*band_chan_inc(iband)`
  as `lpixels(freq_axis)`, and a per-band incs copy with
  `band_chan_inc(iband)` in the freq slot. The blanket `cfg%subim`
  stop-check for `nbands>1` is removed with no replacement guard.
  Build: 0 errors, 0 new warnings across all 4 variants (still exactly 4
  pre-existing GPU-offload linker warnings). `tests/run_tests.sh`: 39/39
  pass (38 from T0-T5 + new section 18) -- the new test (band 1 =
  `TEST_SPLIT_LO` unrestricted, band 2 = the full undivided 200-channel
  `TEST` cube with `subim_chan_blc,subim_chan_trc = 0,101` / `0,200`
  restricting it down to exactly `TEST_SPLIT_HI`'s own channel range)
  came out **bit-identical on the first attempt** to both T5's split-band
  result and the original undivided single-band run -- confirming the
  per-band offset/count/z1-shift arithmetic is correct, not merely
  non-crashing. `nbands=1` bit-identical sweep unaffected (140/140 FITS
  compared, 134 exact + the same 6 pre-existing NaN-artifact diffs seen in
  every prior sweep).

### T7 — Per-Band Bad-Channel Files

- **Objective:** The one remaining genuine gap identified in the T2
  restriction investigation (see T6): let each band flag its own bad
  channels via its own file, rather than blocking `remove_badchan`
  outright for `nbands>1`. Confirmed with user (2026-07-22): a per-band
  file is wanted, required (exact per-band list, same rule as `infileQ`)
  rather than optional/shared, even though this means a `remove_badchan=n`
  multi-band cfg must still supply one placeholder path per band (a minor
  typing cost inherited from the pre-existing single-band key already
  being unconditionally required regardless of the flag's value).
- **Scope:**
  - `band_cfg_t` gains `badchan_file` (comma-list, always required like
    `infileQ` -- not optional like the T6 keys). Legacy scalar
    `cfg%badchan_file` populates from `cfg%band(reference_band)`, so
    single-band cfg files are unaffected.
  - `rm_synthesis.f90`: each non-reference band's own bad-channel file is
    read fresh inside the existing T2 append loop (not kept as a
    persistent per-band array) -- same open/read-list/close pattern as the
    reference band's own `bad_chan`/`flag_arr`, but sized to *that band's
    own* raw NAXIS3 (`band_naxes(iband,freq_axis)`, already stored, T6).
    List entries are raw pixel indices into that band's own file, same
    semantics as the reference band today -- this is what lets bad-channel
    flagging and T6 channel-subimaging compose correctly (a flagged raw
    channel outside the T6-selected range is simply never sampled).
  - The T2 `flag_arr_out` write for each band changes from unconditionally
    marking every channel good to looking up that band's own resolved
    flag array at the raw index `band_chan_blc(iband)+(i-1)*band_chan_inc(iband)`.
  - `NBADGLOB` output header: extended to the sum of bad-channel counts
    across every band (equals the reference band's own count, unchanged,
    whenever `nbands=1`) -- metadata-only, not part of the correctness
    gate (headers aren't compared by `compare_cubes.py`), fixed anyway
    since leaving it silently reference-band-only would misdocument
    multi-band output. `MASKSRC` needed no change -- inspection found its
    final value already derives from `cfg%remove_badchan` directly (a
    single global flag), not from any per-band count, so it was already
    band-agnostic.
  - Remove the blanket `remove_badchan` stop-check for `nbands>1`.
- **Correctness Gate:**
  - `nbands=1` bit-identical sweep, unaffected.
  - New test: flagging raw channel 150 (of the undivided 200-channel
    `TEST` cube) via a plain single-band `badchan_file` must reproduce,
    bit-identically, a 2-band split run (`TEST_SPLIT_LO`+`TEST_SPLIT_HI`)
    that flags the *same* raw channel via band 2's own `badchan_file`
    (channel 50 in band 2's own numbering, since band 2 starts at original
    channel 101). Exercises the new per-band read/apply code specifically,
    composed with the existing frequency-merge architecture -- not merely
    "doesn't crash".
- **Rollback Criteria:** revert to the blanket `remove_badchan` block for
  `nbands>1` if the identity test doesn't come out bit-identical
  (investigate and fix, don't loosen to a tolerance).
- **Effort:** ~0.5-1 session (schema + one read/apply site inside an
  existing loop + one new test; no new arithmetic-correction sites, unlike
  T6, since bad-channel indices are raw pixel indices independent of any
  subim selection).
- **Evidence (2026-07-22):** Implemented as scoped. `band_cfg_t` gained
  `badchan_file` (required per-band list); legacy scalar
  `cfg%badchan_file` populates from `cfg%band(reference_band)`. Existing
  multi-band test cfgs in `run_tests.sh` (sections 15-18, plus the shared
  `make_thesis_cfg` helper in section 16) updated to supply a per-band
  `global_badchan_file` list (`/dev/null` repeated to match band count)
  since the key is now band-count-validated like `infileQ`. Each
  non-reference band's bad-channel file is read fresh inside the T2 append
  loop, sized to that band's own raw NAXIS3, and applied via a per-band
  flag lookup at the T6-aware raw index. `NBADGLOB` now sums bad-channel
  counts across all bands (unchanged for `nbands=1`); confirmed by
  inspection that `MASKSRC` needed no change since its final value already
  derives from `cfg%remove_badchan` directly, not any per-band count.
  Build: 0 errors, 0 new warnings across all 4 variants (still exactly 4
  pre-existing GPU-offload linker warnings). `tests/run_tests.sh`: 41/41
  pass (39 from T0-T6 + new section 19) -- the new test (raw channel 150
  flagged via a single-band `badchan_file` on the undivided `TEST` cube
  vs. the same raw channel flagged via band 2's own `badchan_file`,
  channel 50 in its own numbering, in a 2-band split run) came out
  **bit-identical on the first attempt**, confirming the per-band
  read/apply code is correct, not merely non-crashing. `nbands=1`
  bit-identical sweep unaffected (140/140 FITS compared, 134 exact + the
  same 6 pre-existing NaN-artifact diffs seen in every prior sweep --
  re-verified directly this ticket, by inspecting the actual data at the
  differing positions, that these are the single deliberately
  fully-masked test pixel (25,25) outputting the identical `NaN` in both
  the baseline and current output at every RM bin; `NaN != NaN` under
  IEEE 754 is why `compare_cubes.py --exact`'s `numpy.array_equal` flags
  them, not a behavioural difference).

With T6 and T7 complete, all four T2 restrictions from the investigation
that opened T6 are now resolved: `rem_mean` needed no change (confirmed
correct as-is), RA/Dec-only `subim` needed no change, channel-range
`subim` is now band-aware (T6), and `remove_badchan` is now band-aware
(T7). `remove_qu_bias`/`resiQ`/`slopeQ`/`resiU`/`slopeU` remain
deliberately unimplemented (dead code, pre-existing even for single-band,
kept as a placeholder for future Q/U-vs-I calibration -- not a multi-band
gap). Next planned topic (per user, 2026-07-22): GPU offload for
multi-band.

### T8 — GPU Offload for Multi-Band

- **Objective:** User's assessment (2026-07-22, confirmed by code
  inspection before agreeing, mirroring the T4 discipline): GPU offload
  should already work for multi-band with no compute-path changes, since
  `tile_extract_gpu_rm_blocked` is the same kernel used by both CPU and
  GPU paths, and `prepare_gpu_data`/`prepare_cpu_data` are parameterized
  purely by `nz_out` (already the correct merged total by the time either
  is called -- the multi-band override runs well before). Also checked the
  user's suggestion that auto-tile logic might need updating for
  multi-band cube sizes/mem-fractions: found this is *already* correct,
  not something to change -- `plan_tile`'s RAM/VRAM byte-budget arithmetic
  is entirely `plan%nz_out`-driven, and `plan%nz_out` is assigned after
  the multi-band override, so it already auto-shrinks tiles for a wider
  merged spectrum with zero changes. The one thing genuinely never
  exercised: the two-level VRAM sub-block *staging* path
  (`plan%use_staging = (ny_sub<tile_dec) .and. use_gpu_actual`), unreached
  while `use_gpu_actual` was forced false for `nbands>1` -- needed a real
  GPU-hardware check, not just inspection.
- **Scope:** Remove the blanket `use_gpu_actual` stop-check for `nbands>1`
  (the one line at the top of the T2 scope-narrowing block). No other
  source change.
- **Correctness Gate:**
  - `nbands=1` bit-identical sweep, unaffected.
  - New test (real GPU hardware, RTX 3050 present in this environment):
    the T5 split-band 2-band config run through the GPU binary.
    Non-staged: `AMP.RMCUBE` within `rtol=2e-3` of the CPU reference
    (`tests/output/split_identity.AMP.RMCUBE.FITS`, itself already proven
    bit-identical to the undivided single-band run by T5) plus RM-peak
    validation via `check_rm_peak.py` -- mirrors test 7's own tolerance
    exactly, not a multi-band-specific relaxation. Staged
    (`gpu_vram_mib=1` forcing VRAM sub-block subdivision, confirmed via
    the `Staging sub-blocks:  T` log line): `AMP.RMCUBE` bit-identical to
    the non-staged multi-band GPU run, mirroring test 9's own pattern
    (`OMP_TARGET_OFFLOAD=DISABLED` for determinism, relying on the same
    already-established fact that this kernel has no cross-thread
    reduction, so host-fallback and real-device dispatch of the same
    `-ffast-math`-compiled code are bit-identical).
  - `PHA.RMCUBE` deliberately excluded from the GPU rtol check, matching
    every existing GPU test in this suite -- confirmed by direct
    comparison during this ticket that single-band `serial.PHA.RMCUBE`
    vs `gpu.PHA.RMCUBE` already exceeds `rtol=2e-3` (1.971e-02), so this
    is pre-existing phase-near-low-amplitude sensitivity to `ffast-math`
    reassociation, not a multi-band regression.
- **Rollback Criteria:** revert to the blanket GPU block for `nbands>1` if
  either GPU test fails to meet its tolerance/bit-identical bar
  (investigate and fix, don't loosen the bar).
- **Effort:** ~0.5 session (one line removed; the correctness case rests
  on code already proven generic by T2/T4/T6/T7, verified here on real
  hardware rather than re-derived).
- **Evidence (2026-07-22):** Implemented as scoped -- the blanket
  `use_gpu_actual` stop-check for `nbands>1` removed, no other source
  change. Build: 0 errors, 0 new warnings across all 4 variants (still
  exactly 4 pre-existing GPU-offload linker warnings). `tests/run_tests.sh`:
  44/44 pass (41 from T0-T7 + new section 20, run against a real NVIDIA
  RTX 3050 present in this environment) -- multi-band GPU non-staged run:
  RM peaks correct, `AMP.RMCUBE` within `rtol=2e-3` of the CPU reference
  (max relative diff 9.376e-05, well inside tolerance); staged run
  (confirmed via the `Staging sub-blocks:  T` log line, 32 VRAM
  sub-blocks processed): bit-identical to the non-staged multi-band GPU
  run. Confirmed by direct comparison that excluding `PHA.RMCUBE` from
  the GPU rtol gate is correct practice, not a multi-band-specific
  loosening: single-band `serial.PHA.RMCUBE` vs `gpu.PHA.RMCUBE` already
  exceeds `rtol=2e-3` (1.971e-02) with zero multi-band involvement, worse
  than the multi-band case's own PHA gap (7.652e-03) -- this is why the
  pre-existing single-band GPU tests (7, 7b, 9) never check `PHA.RMCUBE`
  either. `nbands=1` bit-identical sweep unaffected (140/140 FITS
  compared, 134 exact + the same 6 pre-existing NaN-artifact diffs seen
  in every prior sweep). Both the user's core claim (compute path needs
  no changes) and the one correction to it (auto-tile logic needs no
  update either, since it is already `nz_out`-driven) held up under
  direct verification.

All four tickets opened by the "what else would be zero-effort" line of
investigation (T6 channel subimaging, T7 per-band bad channels, T8 GPU
offload) are now complete, alongside the earlier T0-T5 foundation.
`remove_qu_bias` remains a deliberate, known non-gap (dead code,
pre-existing even for single-band, kept as a placeholder for future
Q/U-vs-I calibration). `io_read_threads>1`/`io_overlap` is a genuine,
still-open gap -- see T9 below.

### T9 — Enable `io_read_threads>1`/`io_overlap` for Multi-Band (MUST DO)

- **Status:** not yet scoped or started. Tracked explicitly here per user
  instruction (2026-07-22) after a miscommunication: the user never
  intended this restriction and was questioning why it existed at all,
  not asking for it to be added -- this is a required fix, not an
  optional enhancement, and matters for large multi-band cubes.
- **Objective:** remove the blanket
  `cfg%io_read_threads.gt.1 .or. cfg%io_overlap` stop-check for
  `nbands>1` ([rm_synthesis.f90:932-936](../../src/rm_synthesis.f90#L932-L936)),
  following the same discipline as T6/T7/T8: investigate by code
  inspection first (is the restriction genuinely necessary, or was it
  over-cautious like channel subimaging/GPU turned out to be?), scope
  concretely, implement, verify via a bit-identical/tolerance test plus
  the standing `nbands=1` sweep, record Evidence.
- **Not yet investigated**: whether the per-band tile-read loop (T2)
  already reuses the parallel-channel-read/async-write infrastructure
  safely for multi-band, or whether genuine new work is needed (e.g. the
  per-band FTGSVE calls added in T2/T6 currently run as plain sequential
  calls, forced by `io_read_threads` being blocked at 1 for `nbands>1` --
  worth checking whether that per-band loop can be parallelized the same
  way the reference band's own channel-chunked read already is).
- **Evidence (2026-07-22):** Investigated by code inspection before
  touching anything, per the user's direct challenge to justify the
  restriction. Found: the `io_read_threads` parallel channel-split read
  ([rm_synthesis.f90:3391-3404](../../src/rm_synthesis.f90#L3391-L3404))
  only ever touches the *reference* band's own channel range and buffer
  offset -- the other bands' reads happen afterward, sequentially,
  entirely outside that parallel region, so there is no data race to
  introduce. The `io_overlap` ping-pong double-buffering
  ([rm_synthesis.f90:2350-2364](../../src/rm_synthesis.f90#L2350-L2364)) is
  sized from `nz_out`, already the correct post-merge multi-band total at
  that point in the program, and the write-dispatch logic
  (`populate_write_job`) operates purely on already-merged output tile
  arrays with no notion of band count at all. Same over-cautious-blanket-
  restriction pattern as T6/T7/T8. Implemented: removed the stop-check,
  no other source change. Build: 0 errors, 0 new warnings across all 4
  variants (still exactly 4 pre-existing GPU-offload linker warnings).
  `tests/run_tests.sh`: 49/49 pass (44 from T0-T8 + two new sections) --
  the T5 split-band config, forced into 7 tiles (uneven remainder, same
  shape as the existing single-band io_overlap/io_read_threads tests):
  with `io_overlap=y`, confirmed 7 tiles ran, no overlapping tile writes,
  and all 4 output products bit-identical to the existing single-tile
  split-band reference; with `io_read_threads=4`, confirmed the 4-handle
  parallel-read path was taken and all 4 output products bit-identical to
  the same reference. `nbands=1` bit-identical sweep unaffected (140/140
  FITS compared, 134 exact + the same 6 pre-existing NaN-artifact diffs
  seen in every prior sweep). Both code-inspection predictions confirmed
  empirically on the first attempt.

### T10 — Cross-Band Geometry Alignment (`reproject_cubes`)

- **Status:** done.
- **Objective:** T0-T9 above assume every band already shares one sky
  grid (same RA/Dec pixel-for-pixel) and, from T3 onward, one angular
  resolution -- real multi-band ASKAP data does not arrive that way.
  Before any of the multi-band merge machinery above can run on genuine
  survey data, the bands need a preprocessing step that puts them on a
  common WCS grid. `reproject_cubes` is that step: a standalone tool
  (own binary, own build target, not linked into `rm_synthesis` at all --
  see the Makefile's own comment on why it's kept off the main build
  graph), built on Starlink AST for WCS handling and `astResampleR` for
  the actual resampling.
- **Scope:** `src/reproject_cubes.f90` (single file). Three footprint
  modes (`intersection`/`union`/`reference`) computed from N input
  files' own sky footprints against a reference file, via AST
  `SkyFrame`-to-`SkyFrame` conversion (not whole-`FrameSet` conversion --
  `astConvert`'s domain search does not recurse into a `CmpFrame`'s
  internal components, so aligning two compound "STOKES-SKY-SPECTRUM"
  frames needs the sky axes picked out first). Output axis layout always
  puts the 2 sky axes first (RA-fastest-on-disk, matching
  `rm_synthesis`'s own tile-read assumption), with full header
  propagation: per-axis WCS keywords with correct CRPIX shift for
  crop/grow, `PCi_j`/`CDi_j` sky rotation, `BUNIT`/`BMAJ`/`BMIN`/`BPA`/
  `OBJECT`/etc. via a generic verbatim header copy. Reads/resamples/
  writes in `mem_frac_ram`-budgeted blocks of planes (same tile-planning
  concept as `rm_synthesis`'s own, see `get_mem_total_kb`), OpenMP-
  parallelised across planes within a block -- each thread builds its
  own private AST pixel-to-pixel Mapping from scratch rather than
  sharing one across threads, since this Fortran AST binding exports no
  `astLock_`/`astUnlock_` for handing an AST Object between threads
  (checked against the actual linked `.so`'s symbol table, not assumed).
- **Correctness gate:** resampled output matches independently-computed
  (Python/astropy) ground truth exactly at spot-checked pixels in both
  reference and intersection/union modes; intersection-mode `CRPIX`
  shift verified against `rm_synthesis`'s own existing subimage-CRPIX
  formula; union-mode uncovered-area NaN count matches the expected
  geometric area exactly; a genuine axis-order bug (`FTGSVE` filling its
  output array in ascending-axis-number order among non-degenerate axes,
  not always sky-first) caught by a non-adjacent-sky-axis fixture,
  fixed, and reverified byte-identical against the pre-fix output on
  every other fixture; 25 repeated stress runs at default thread count,
  no failures.
- **Effort:** built incrementally across several commits (footprint
  computation, resampling, output writing, CLI, OpenMP block I/O,
  header metadata, sky-rotation propagation) -- see `git log
  src/reproject_cubes.f90` for the full commit-by-commit record; not
  reconstructed here since each commit message already documents its
  own change and verification.
- **Note on real-data readiness:** ASKAP data specifically needed the
  `CROTA`/`PCi_j`/`CDi_j` sky-rotation propagation (not just per-axis
  CRVAL/CRPIX/CDELT) -- confirmed present on real ASKAP headers, added
  after being initially flagged as a documented future gap rather than
  assumed absent.

### T11 — Cross-Band Resolution Matching (`gaussft_mod`, `commonbeam_mod`, `convolve_cubes`)

- **Status:** done.
- **Objective:** T10 solves grid alignment; this ticket solves the
  matching problem T3's own `δRM`/resolution diagnostic exists to warn
  about in the first place -- real bands (and even single bands: ASKAP
  per-channel restoring beams vary continuously across a band, not just
  between bands) do not share one angular resolution, and merging
  channels of different resolution into one RM synthesis without first
  convolving them to a common beam is exactly the kind of silent
  correctness gap this project's whole design philosophy (loudly refuse
  or visibly warn, never silently produce a misleading answer) exists to
  close. Three new pieces, each independently reusable, deliberately
  split the way `reproject_cubes` already splits computation from I/O:
  - **`src/gaussft.f90`** (`gaussft_mod`): pure computation, no FITS I/O.
    Given one image plane, its own source elliptical-Gaussian PSF, and a
    target PSF (same pixel grid), returns the plane convolved from one
    to the other via FFT-domain deconvolve-then-reconvolve
    (multiply by target-PSF-FT / source-PSF-FT, inverse-transform).
    Split into `plan_convolution`/`convolve_to_beam`/
    `destroy_convolution_plan` specifically so a caller can parallelise
    across planes with OpenMP: FFTW's planner is not thread-safe and
    must run once, serially, but a single plan's `fftw_execute_dft` (the
    "new-array execute" form) is documented safe for concurrent use by
    multiple threads supplying their own arrays -- verified directly,
    not just trusted: 16 threads sharing one plan across 64 planes with
    distinct per-plane beams reproduces a serial run bit-for-bit.
    Corrected a real amplitude bug versus the original `src/gaussft.f`
    prototype (and its direct Python port, and upstream `racs_tools`,
    which all carry the same bug): the closed-form 2D-Gaussian FT
    amplitude is `2*pi*sigma_x*sigma_y`, not
    `sqrt(2*pi*sigma_x*sigma_y)` -- confirmed via closed-form derivation,
    direct numerical integration to 15 significant figures, and
    cross-checking against MIRIAD's own independent `gaufac` formula
    (`au2.gauss_factor`, also in this repo).
  - **`src/commonbeam.f90`** (`commonbeam_mod`): given N per-channel
    beams, finds the smallest common beam every one of them can be
    deconvolved from. A simpler "just take the largest beam" shortcut is
    not generally correct -- real ASKAP per-channel BPA varies by more
    than 90 degrees across a band, confirmed on a real cube's own BEAMS
    table, so the single largest-major-axis beam does not always
    deconvolve every other channel's beam. Follows the standard approach
    (CASA `ia.commonbeam()`, the `radio_beam` Python package): sample
    each beam's boundary, reduce to the 2D convex hull, fit the minimum-
    volume enclosing ellipse (Khachiyan's algorithm), validate against
    every real input beam via the Sault/MIRIAD "gaupar" deconvolution
    formula (same formula family as `au2.gauss_factor` above), retry
    with a larger safety margin if needed. One deliberate departure from
    `radio_beam`'s own algorithm: since a beam here is pure shape with
    no position, the ellipse fit uses the simpler origin-centred variant
    of Khachiyan's algorithm rather than `radio_beam`'s general
    free-centre one -- verified against `radio_beam` 0.3.9 itself on a
    real 286-channel ASKAP BEAMS table (within 0.003 arcsec on
    BMAJ/BMIN, PA matching mod 180 degrees, and independently confirmed
    deconvolvable from all 286 real beams via `radio_beam`'s own
    `deconvolve_optimized`).
  - **`src/convolve_cubes.f90`**: the main program driving both modules,
    mirroring `reproject_cubes`' own I/O-vs-computation split. Reads
    per-channel PSFs via a CASA-style `BEAMS` binary table extension
    (auto-detected) or a portable ASCII/CSV fallback (`cfg/
    example_beamLog.txt`/`.csv`); a channel is bad -- written as an
    all-NaN plane, not convolved -- if it's missing from the file
    entirely, or present with BMAJ or BMIN equal to 0 (either alone is
    enough; not an AND of both). Pools every good channel across ALL
    input files before calling `commonbeam_mod` once, so multi-band
    support needs no extra machinery: every band gets convolved to the
    exact same shared target. `max_common_bmaj` lets a user cap the
    auto-derived target and refuse to proceed if it comes out coarser
    than expected, rather than silently convolving to an unchecked
    resolution. `mem_frac_ram`-budgeted block I/O + OpenMP, same concept
    as `reproject_cubes`. One correctness-critical piece specific to
    this tool: FITS `BMAJ`/`BMIN`/`BPA` is a SKY-frame convention
    (position angle from North through East), while `gaussft_mod`'s own
    convention is PIXEL-frame -- converted via
    `bpa_pixel = atan2(sign(CDELT2)*cos(theta), sign(CDELT1)*sin(theta))`,
    derived from the local tangent-plane geometry and checked against
    both real ASKAP CDELT signs' special cases (North, East) by hand
    before being verified empirically via a bit-exact identity
    round-trip (target beam set equal to one channel's own native beam
    reproduces that channel's input data to the last bit).
- **Correctness gate:** `gaussft_mod`'s own identity/asymmetric-beam/
  thread-safety tests (above); `commonbeam_mod` against `radio_beam`
  (above); `convolve_cubes`' BEAMS-table and ASCII/CSV readers verified
  against real ASKAP header conventions and cross-checked against each
  other; bad-channel union (degenerate beam, badchan_file, BMAJ-or-BMIN
  zero) verified with exact expected counts; `max_common_bmaj` cutoff
  verified to refuse correctly; full pipeline smoke-tested against a
  genuine cutout of real ASKAP data (`/data1/tmp/cutout-stokesQ.fits`)
  with no NaN/Inf and output stats matching input to high precision for
  already-near-common-resolution channels.
- **Not yet done:** a full run against the complete 23GB real cube
  (only cutouts and synthetic data verified so far).
- **Build/packaging:** own Makefile targets (`make reproject_cubes`,
  `make convolve_cubes`, both independent of the main `rm_synthesis`
  build graph), both packaged into `docker/dockerfile` (`libfftw3-dev`
  added for `convolve_cubes`' FFTW3 dependency, `libstarlink-ast-dev`
  and related packages already present for `reproject_cubes`) --
  verified via a real `docker build` + `docker run`, not just Makefile
  inspection.

### T12 — `rm_synthesis` Beam-Metadata Propagation and Input Safety

- **Status:** done.
- **Objective:** two gaps found while making T10/T11's outputs
  actually usable as `rm_synthesis` inputs, neither specific to
  multi-band but both surfaced by working through this pipeline
  end-to-end.
  1. `rm_synthesis` never propagated `BMAJ`/`BMIN`/`BPA` (or any beam
     metadata at all) to any of its 8 output products (AMP/PHA cubes,
     mask, nvalid, peak/rmpeak/angpeak/snr maps) -- confirmed by grep,
     zero references anywhere before this ticket. A user running
     `rm_synthesis` on a `convolve_cubes`-processed (single, well-
     defined resolution) input had no way to recover that resolution
     from the output files at all.
  2. All of `rm_synthesis`'s own input cubes (Q/U/I/mask, units
     21/22/40/45) were opened `READWRITE` (`rwmode=1`) despite the code
     never writing to any of them -- confirmed by grep (no `FTPKYx`/
     `FTP2Dx`/`FTPPRx`/`FTPCLx`/`FTPHIS`/`FTPSSE` call anywhere in the
     file targets those units) and by the fact that this file's own
     parallel tile-reader threads for the same files already
     independently open `READONLY`. An unnecessary, real (if latent)
     risk to irreplaceable input science data, with no upside.
- **Scope/decisions (confirmed with user):**
  - `BMAJ`/`BMIN`/`BPA` propagated from the input Q cube's primary
    header to all 8 outputs, unchanged from whatever the input carries
    (including when meaningless -- see next point).
  - If the input has `CASAMBM=T` (a genuine per-channel-varying beam,
    e.g. an un-convolved CASA multi-beam cube), the propagated scalar
    BMAJ/BMIN/BPA is only the input's own nominal/reference value and
    means nothing on its own -- but instead of hiding that, the output
    ALSO gets `CASAMBM=T` plus the input's own real per-channel `BEAMS`
    binary table, attached as its own extension HDU, plus `HISTORY`
    cards explaining why. Applies to AMP/PHA and, when `cubestat=y`,
    the PEAK/RMPEAK/ANGPEAK/SNR maps -- every one of these is derived
    from the actual flux-bearing data, so "which beams went into this"
    is a real provenance question for all of them.
  - Deliberately NOT applied to MASK.CUBE.FITS or NVALID.MAP.FITS: both
    are pure per-pixel/per-channel validity bookkeeping (a flag, a
    count), not flux data -- nobody convolves a flag table, and a
    BEAMS extension there would only invite a confusing "why does this
    have a beam?" instead of the intended "have we processed this
    correctly?". Still carry the plain BMAJ/BMIN/BPA scalar, unchanged.
  - Multi-band mode: the above only ever looked at the reference band's
    own Q file -- extended to cross-check every non-reference band's
    own primary header against the reference (BMAJ/BMIN/BPA numeric
    mismatch beyond a small tolerance, presence mismatch, or its own
    `CASAMBM=T`), warning with the actual differing values per band
    rather than silently reflecting only one band's metadata. Not a
    hard error -- this project already hard-stops on genuine geometry
    mismatches (WCS/NAXIS/frequency-axis-index) earlier in the same
    multi-band file-opening loop; a beam-metadata mismatch gets the
    same warn-and-continue treatment as the single-band CASAMBM case,
    since RM synthesis itself does not depend on beam metadata for
    correctness.
  - `rwmode` changed from 1 to 0 for units 21/22/40/45. Also brings
    this file in line with `docs/user/ARCHITECTURE.md`'s own documented
    CFITSIO lesson (a real historical SIGSEGV, see its "History:
    `io_write_threads>1` was unsafe" postmortem): CFITSIO aliases
    repeat `READWRITE` opens of an already-open file onto one shared
    buffer, but exempts `READONLY` opens from that aliasing by design.
- **Correctness gate:** identity check (target beam == one channel's
  own native beam reproduces that channel's input exactly) already
  covers `convolve_cubes`' own math; for T12 specifically -- injected
  real BMAJ/BMIN/BPA into a test cube and confirmed exact propagation
  to all 8 outputs, confirmed clean degradation when absent; injected a
  genuine `CASAMBM=T` + `BEAMS` table and confirmed `CASAMBM`/BEAMS/
  HISTORY land exactly on AMP/PHA/PEAK/RMPEAK/ANGPEAK/SNR and nowhere
  else (MASK/NVALID untouched), BEAMS table content byte-identical to
  the input's own, primary pixel data uncorrupted despite the HDU-
  append-then-return-to-primary sequence; injected mismatched per-band
  beams in a 2-band multi-band run and confirmed the cross-band warning
  fires with the correct differing values, and confirmed it stays
  silent (no false positive) when bands genuinely match; fed a real
  `convolve_cubes`-produced NaN bad-channel plane into `rm_synthesis`
  with no `badchan_file` at all and confirmed automatic exclusion via
  the existing (default-on) NaN-check mechanism -- `NVALID` correctly
  read 5 (not 6) everywhere, AMP cube had no NaN.
- **Build/test:** all 4 build flavours (`scratch/make_all.sh`) clean;
  full `tests/run_tests.sh` 49/49 pass, re-run clean after every
  sub-change in this ticket.

### T13 — `match_cubes`: In-Memory Chaining of Reproject + Convolve

- **Status:** done.
- **Objective:** `reproject_cubes` and `convolve_cubes` (T10/T11) each
  read a full input cube and write a full output cube to disk. Run
  back-to-back -- the recommended pipeline before multi-band
  `rm_synthesis` -- the reprojected intermediate is written in full, then
  immediately re-read in full by `convolve_cubes`, for an artifact nobody
  wants on its own. At the user's real scale (200GB+ cubes) this doubles
  disk I/O and disk space for no reason. `match_cubes` is a new
  standalone tool that can run either stage alone, or both chained
  THROUGH MEMORY with no intermediate FITS file, order configurable
  (default convolve-then-reproject, confirmed with the user as the
  right default -- see below).
- **Order matters, not just for tidiness:** convolving to the common
  target beam before reprojection low-pass-filters the image before
  `astResampleR`'s linear interpolation touches it, so resampling
  operates on smooth, well-sampled data rather than a band's own native
  (possibly only marginally Nyquist-sampled) sharp PSF -- avoiding
  interpolation/aliasing error that convolving afterward cannot undo,
  since the error is already baked into the resampled pixel values by
  then (the same reasoning as an anti-alias filter before downsampling
  in ordinary signal processing). It also usually costs less: in
  `union` footprint mode the reprojected output grid is larger than any
  input's own native grid, so convolving first does the expensive FFT
  work on the smaller native footprint. `reproject_convolve` order is
  not wrong, just not the default, and remains fully selectable.
- **Scope/decisions (confirmed with user):**
  - Neither `reproject_cubes.f90` nor `convolve_cubes.f90` is touched --
    both remain fully independent, already-tested, already-shipped tools
    for anyone who wants just one stage without `match_cubes`.
    `src/match_cubes.f90` therefore duplicates (adapts, not `use`s) the
    subroutines it needs from both, rather than extracting a shared
    module -- a real, accepted maintenance cost in exchange for zero
    regression risk to two tools already in production use.
  - `rm_synthesis` itself is out of scope here -- feeding it directly,
    with no intermediate file at all, is a separate, harder design
    question the user flagged for later (RM synthesis wants full
    per-pixel spectra across the whole frequency axis for a spatial
    tile, a fundamentally different access pattern than reproject/
    convolve's own plane-at-a-time one -- the "corner turn" `rm_synthesis`
    already performs internally exists for exactly this reason).
  - Axis-scope handling is deliberately asymmetric by stage:
    `stages=reproject` alone keeps `reproject_cubes`' own fully general
    N-dimensional "other axes" handling (any number of non-sky axes) --
    no new restriction. `stages=convolve` or `stages=both` adopt
    `convolve_cubes`' own existing restriction instead (exactly 2 sky
    axes + 1 FREQ axis, every other axis degenerate) -- not a new
    limitation, the scope `gaussft_mod`'s own per-channel convolution
    already has.
  - Tool named `match_cubes` (matches the sky grid AND the resolution),
    confirmed with the user over an alternative (`chain_cubes`).
- **Implementation:** one new source file, `src/match_cubes.f90`
  (~2400 lines). Pre-scan phase (once per run, order-independent of each
  other): reproject footprint bounds (adapted from `reproject_cubes`'
  own main-program logic) and convolve beam-pooling/common-beam
  derivation (adapted from `convolve_cubes`' own), each only run when
  its stage is active. Per-file processing dispatches to
  `process_one_file_general` (stages=reproject alone, a near-verbatim
  port of `reproject_cubes`' `write_reprojected_file`/`read_one_block`/
  `write_one_block`) or `process_one_file_restricted` (stages=convolve
  or both -- the new logic: one `!$omp parallel` region per block,
  `!$omp single` read / `!$omp do` per-plane compute / `!$omp single`
  write, coexisting two different thread-safety models in the same
  region -- each thread's own *private* AST pixel-to-pixel Mapping,
  required since AST objects can't be shared across threads without
  lock/unlock this binding doesn't have, alongside ONE *shared*
  `gaussft_mod` FFTW plan built once outside the parallel region, safe
  to share via `fftw_execute_dft`'s "new-array execute" form). A bad
  channel (from the convolve pre-scan's own bad-channel policy) skips
  BOTH stages entirely and writes NaN straight to the final buffer,
  cheaper than running either stage on known-bad data.
- **Two real bugs found by this ticket's own verification, both fixed
  before landing:**
  1. `process_one_file_restricted`'s header writer explicitly copies
     sky-axis WCS (from the reference) and the FREQ axis' own WCS (from
     the input) when reproject is active, but originally never copied
     any OTHER degenerate axis' WCS (e.g. a size-1 STOKES axis) at all --
     `exclude_axis_wcs=true` in the generic header copy correctly
     excludes ALL axis-indexed keywords from the verbatim pass-through,
     but nothing then explicitly restored a degenerate axis outside
     sky1/sky2/freq_axis, silently losing its CTYPE/CRVAL/CRPIX/CDELT.
     Caught by this ticket's own chaining-equivalence diff (below)
     against the two-step disk-based reference, which does carry this
     through via `reproject_cubes`' own general "other axes" handling.
     Fixed: an explicit loop over every axis that isn't sky1/sky2/
     freq_axis, copying its own keywords through unchanged.
  2. A serious one: `gaussft_mod`'s FFTW plan is sized for one specific
     `(nx,ny)` and must never be executed against arrays of a different
     size (FFTW's own "new-array execute" contract requires matching
     dimensions -- using a mismatched-size array is undefined behaviour,
     not a clean error). `order=reproject_convolve` convolves on the
     OUTPUT grid, not the native one -- a genuinely different size
     whenever `footprint_mode` crops or grows it -- but the plan was
     being built once, unconditionally, for the native grid regardless
     of order. Manifested as a real crash (`munmap_chunk(): invalid
     pointer`, heap corruption) the first time `order=reproject_convolve`
     was actually exercised end-to-end, not a subtle numerical
     discrepancy -- exactly the kind of bug the chaining-equivalence
     verification below exists to catch before it reaches real data.
     Fixed: plan sized for whichever grid convolution will actually run
     on for the selected order (native for convolve_reproject or
     convolve-alone, output for reproject_convolve).
- **Correctness gate (all four verified bit-identical, header AND
  data, on a genuine 2-band scenario with offset grids, varying
  per-channel beams including a real bad channel per band, and an
  intersection-mode footprint crop):**
  1. `match_cubes stages=reproject` vs standalone `reproject_cubes` on
     the same inputs: bit-identical.
  2. `match_cubes stages=convolve` vs standalone `convolve_cubes`: 
     bit-identical.
  3. `match_cubes stages=both order=convolve_reproject` vs the two *old*
     standalone tools run back-to-back through a real disk intermediate
     (`convolve_cubes` -> `_CONV.FITS` -> `reproject_cubes` on that
     file): bit-identical -- confirms the in-memory chain changes
     nothing numerically, only removes the disk round-trip.
  4. `match_cubes stages=both order=reproject_convolve` vs the same two
     tools in the other order (`reproject_cubes` -> `_REPROJ.FITS` ->
     `convolve_cubes` on that file): bit-identical, after the two fixes
     above.
  - Incidental finding while constructing gate 3/4's own two-step disk
    reference: `reproject_cubes`' own output never carries a `BEAMS`
    binary table forward (it only ever copies primary-header cards, and
    never touches binary table extensions at all) -- a real,
    pre-existing characteristic of that already-shipped standalone tool,
    not a bug introduced or fixed here. Worked around for the
    verification itself by using an explicit ASCII beam log for both the
    two-step reference and `match_cubes` in the `reproject_convolve`
    gate, rather than relying on BEAMS-table auto-detection on an
    already-reprojected intermediate.
- **Build/packaging:** `make match_cubes` (own Makefile target, needs
  both `AST_LIBS` and `FFTW_LIBS` plus `gaussft_mod`/`commonbeam_mod`/
  `ast_grf_stub.o`, mirroring `reproject_cubes`/`convolve_cubes`' own
  targets exactly); `docker/dockerfile` builds it alongside the other
  two, needing no new apt packages (union of two already-packaged
  dependency sets).
- **Build/test:** all 4 `rm_synthesis` build flavours
  (`scratch/make_all.sh`) clean; `reproject_cubes`/`convolve_cubes`/
  `match_cubes` all build clean; full `tests/run_tests.sh` 49/49 pass,
  unaffected (this ticket touches neither `rm_synthesis` nor the two
  existing standalone tools' own source).

### T14 — `CASAMBM`/`BEAMS` Propagation for `reproject_cubes`/`convolve_cubes`/`match_cubes`

- **Status:** done.
- **Objective:** T12 gave `rm_synthesis` a genuine `CASAMBM`/`BEAMS`
  passthrough for its own outputs, downstream of the cross-band
  preprocessing toolchain. The toolchain itself (`reproject_cubes`,
  `convolve_cubes`, and T13's `match_cubes`) had no equivalent: a
  multi-beam input's real per-channel restoring beam was silently lost
  the moment it passed through either tool, well before `rm_synthesis`
  ever saw it. This ticket closes that gap upstream, using the same
  `CASAMBM=T` + `BEAMS` binary table extension convention T12 already
  established, "much the same way" (the user's own framing) but adapted
  to what each tool actually does to the beam.
- **Decision (confirmed with user):** the two tools' own beam semantics
  differ, so their propagation rule differs too, deliberately:
  - `reproject_cubes` (and `match_cubes` with `stages=reproject`) never
    touch the beam -- reprojection is a pure spatial resample, and
    `BMAJ`/`BMIN`/`BPA` are stored in sky degrees, unaffected by it. So
    the rule mirrors `rm_synthesis`' own T12 rule exactly: if the input
    has `CASAMBM=T`, copy its real `BEAMS` table through to the output
    unchanged.
  - `convolve_cubes` (and `match_cubes` whenever `convolve` is active),
    by contrast, genuinely changes the beam -- every good channel is
    convolved to one common target beam. Copying the input's own BEAMS
    table through unchanged would therefore be actively wrong (it would
    describe the INPUT's pre-convolution beams, not the output's). The
    user was asked directly whether the output should always get a
    synthesized `CASAMBM=T`/`BEAMS` table, or only when the input itself
    had one (mirroring `rm_synthesis`' trigger condition) -- confirmed:
    **always attach**, since `convolve_cubes` always tracks a real
    per-channel good/bad split internally regardless of the input's own
    format (a plain scalar `BMAJ`/`BMIN`/`BPA` header, or an external
    ASCII beam log, carries the same good/bad information a `BEAMS`
    table would), so a scalar-only output would misrepresent every
    bad/NaN (skipped) channel as sharing the common target beam.
- **Real BEAMS table layout confirmed against production data:** read
  directly off `/data1/tmp/cutout-stokesQ.fits` (the same real ASKAP cube
  T10-T12 already used) via `astropy`: `BinTableHDU`, `EXTNAME='BEAMS'`,
  5 columns `[BMAJ,BMIN,BPA,CHAN,POL]`, formats `[1E,1E,1E,1J,1J]`,
  `BMAJ`/`BMIN` in arcsec, `BPA` in deg, `CHAN` 0-indexed, `POL` always 0
  (single Stokes product per file). `convolve_cubes`'/`match_cubes`' own
  new synthesized-table writer matches this layout exactly, including the
  `POL=0` convention.
- **Implementation:**
  - `reproject_cubes.f90`'s `write_reprojected_file` and
    `match_cubes.f90`'s `process_one_file_general`: after the existing
    generic header copy (which already carries the scalar `CASAMBM`
    keyword through verbatim as a raw header card, but cannot reach a
    separate extension HDU), check `CASAMBM` on the freshly-written
    output header; if true, open the INPUT file on its own dedicated unit
    (`reproject_cubes.f90`: 45; `match_cubes.f90`: 45, disjoint from that
    file's other unit numbers -- see each file's own comment), move to
    its `BEAMS` HDU, and `ftcopy` it onto the output, then `ftmahd` back
    to the output's own primary HDU. Missing `BEAMS` despite `CASAMBM=T`
    (a malformed input) is a runtime warning, not a hard error -- matches
    `rm_synthesis`' own T12 handling of the same edge case.
  - `convolve_cubes.f90`'s new `write_beams_table` and
    `match_cubes.f90`'s new `write_beams_table_match` (a verbatim port,
    per this project's standing decision not to share code between
    `match_cubes.f90` and the two standalone tools): builds the 5-column
    table described above via `FTIBIN`/`FTGCNO`/`FTPCLE`/`FTPCLJ`, one
    row per channel -- `(tgt_bmaj, tgt_bmin, tgt_bpa)` for a good channel,
    `(tiny(1.0), tiny(1.0), 0.0)` for a bad one (the same degenerate
    sentinel this project's own `BEAMS`-table readers already treat as
    "no valid beam" -- see T11's own read-side handling), `CHAN=ich-1`,
    `POL=0`. Called once, right before the output file closes (after all
    per-channel data has already been written), alongside an explicit
    `FTPKYL(...,'CASAMBM',.true.,...)` on the primary header --
    `copy_generic_header_convolve`'s/`copy_generic_header_match`'s own
    existing exclusion of `BMAJ`/`BMIN`/`BPA`/`CASAMBM` from the generic
    verbatim copy (needed since those are recomputed, not copied) already
    covered this; only its own comment needed updating, since it
    previously (accurately, at the time) said this program never creates
    a `BEAMS` extension at all.
- **Correctness gate (real per-channel-beam fixtures, one genuine bad
  channel, both good and bad rows exercised):**
  1. `reproject_cubes` alone on a `CASAMBM=T` input: output `BEAMS` table
     byte-identical to the input's own (all 3 good channel values, the
     bad-channel row, in original order).
  2. `convolve_cubes` alone: output `BMAJ`/`BMIN`/`BPA` scalar matches
     the auto-derived common beam; `BEAMS` table has that same common
     beam on every good channel and the `tiny(1.0)` sentinel on the known
     bad channel (channel 4, confirmed absent from the ASCII beam log
     fixture); convolved data plane for the bad channel confirmed
     all-NaN, good-channel planes confirmed non-NaN.
  3. `match_cubes stages=reproject` vs standalone `reproject_cubes`:
     `BEAMS` table (and every other header field, and the data)
     identical -- only the `HISTORY` provenance text differs, expectedly,
     by tool name (`reproject_cubes:` vs `match_cubes:`, a length
     difference that incidentally also changes whether CFITSIO wraps
     that one `HISTORY` card into two -- confirmed benign, unrelated to
     this ticket's own changes, by checking every OTHER header key,
     the data, and the `BEAMS` table are all exactly equal).
  4. `match_cubes stages=convolve` vs standalone `convolve_cubes`: same
     result, `BEAMS` table and data identical.
  5. `match_cubes stages=both order=convolve_reproject` vs the two OLD
     standalone tools run back-to-back through a real disk intermediate
     (`convolve_cubes` -> `reproject_cubes` on that intermediate,
     confirming `reproject_cubes`' own new copy-through logic correctly
     forwards the BEAMS table `convolve_cubes` just synthesized on the
     intermediate file): `BEAMS` table and data identical.
  6. `match_cubes stages=both order=reproject_convolve` vs the same two
     tools in the other order: `BEAMS` table and data identical.
- **Build/test:** `convolve_cubes`/`reproject_cubes`/`match_cubes` all
  build clean; all 4 `rm_synthesis` build flavours (`scratch/make_all.sh`)
  clean; full `tests/run_tests.sh` 49/49 pass, unaffected (this ticket
  touches neither `rm_synthesis` nor its own test fixtures).

### T15 — End-to-End Preprocessing Test: Genuinely Mismatched Bands, Fixed, Correct Answer Recovered

**Gap found, not assumed:** while writing user-facing documentation
(`docs/user/EXAMPLES.md`, a scenario cookbook covering "my bands don't share
a sky grid/resolution"), the user asked directly whether a multi-band
*test* backs the recipes being documented -- prompting a check of
`tests/run_tests.sh`'s actual coverage rather than assuming it existed.
Confirmed: sections 15-22 test `rm_synthesis`'s own multi-band SCHEMA
handling in isolation (comma-list parsing, geometry-mismatch REJECTION,
frequency merge, split-band identity, per-band channel/bad-channel
handling, GPU/IO-parallelism for multi-band) -- all real, all passing,
but every one of them uses bands that are either already matched by
construction, or deliberately mismatched specifically to prove
`rm_synthesis` correctly REFUSES them. Sections 30/33-36 test
`reproject_cubes`/`convolve_cubes`/`match_cubes`'s own internal
correctness (skip-if-already-matched, `io_overlap` consistency) --
section 33's own `convolve_cubes` test in particular applies a fake
uniform beam log to a single copy of the same test cube, not two
genuinely different-resolution inputs. **No test anywhere took
genuinely mismatched multi-band data, ran it through the preprocessing
toolchain to fix the mismatch, and confirmed `rm_synthesis` then
combined the result correctly** -- the actual point of the toolchain,
never verified end to end as one thing.

**New fixture (`tests/make_test_cubes.py`):** `TEST_BAND2_UNMATCHED.Q/
U.FITSCUBE` -- the same underlying injected-source signal as the
existing `TEST_BAND2.Q/U.FITSCUBE` (so the correct post-fix answer is
exactly `truth.json`'s own `SOURCES`), with the same `CRVAL1` shift
`TEST_BAND2_MISMATCH` already used for the grid half of the mismatch,
plus two new flat per-channel ASCII beam logs (`band1_beamlog.txt` at a
fixed 10", `band2_beamlog.txt` at a genuinely different 20",
`convolve_cubes`'s own portable format, same convention section 33
already uses) for the resolution half -- `convolve_cubes` has real
smoothing work to do on band 1, not a no-op. `band2_beamlog.txt`
doubles as the resolution mismatch for the ALREADY-grid-matched
`TEST_BAND2` too, letting the grid and resolution mismatches be tested
independently as well as together.

**Three new permanent sections (`tests/run_tests.sh` 38-40), each
isolating one failure mode, each ending in the same real assertion --
`check_rm_peak.py` against `truth.json` on the final combined
`rm_synthesis` output, not just "the preprocessing step didn't crash":**
1. **Grid-only** (§38): `TEST_BAND2_MISMATCH` (existing fixture, no
   beam mismatch) fixed by `reproject_cubes` alone, then `rm_synthesis`
   recovers `src_A`/`src_B`.
2. **Resolution-only** (§39): `TEST_BAND2` (grid-matched) + the new beam
   logs, fixed by `convolve_cubes` alone (band 1 genuinely smoothed
   10"->~20", confirmed by the derived common beam in the tool's own
   log output), then `rm_synthesis` recovers both sources.
3. **Both together** (§40): `TEST_BAND2_UNMATCHED` + beam logs, fixed by
   one `match_cubes stages=both order=convolve_reproject` call (all 4
   files -- both bands' own Q/U -- passed together, so the derived
   common beam and grid account for every input, not just the
   mismatched one), then `rm_synthesis` recovers both sources.

**A real mechanical finding while building this, not assumed -- found,
then fixed, not left as a documented quirk:**
`reproject_cubes`/`convolve_cubes`/`match_cubes`'s own `strip_fits_ext`
output-naming helper originally only stripped a trailing `.fits`/
`.FITS` (hardcoded 5-character check) -- this project's own test
fixtures use `.FITSCUBE` (9 characters), which didn't match, so the
real output filename was the FULL original filename with
`_REPROJ.FITS`/`_CONV.FITS`/`_MATCHED.FITS` simply appended (e.g.
`TEST_BAND2_MISMATCH.Q.FITSCUBE_REPROJ.FITS`, not an extension swap) --
confirmed directly by reading the helper's own source rather than
assumed from the shorter, `.fits`-suffixed examples used elsewhere in
the docs. Rather than leave this as a documented limitation, all three
copies of `strip_fits_ext` (`src/reproject_cubes.f90`,
`src/convolve_cubes.f90`, `src/match_cubes.f90`) were generalized to
strip whatever the trailing extension actually is (anything after the
last `.` in the filename's own basename, ignoring any path component),
not a hardcoded `.fits` check -- so a `.FITSCUBE` input now produces a
clean `TEST_BAND2_MISMATCH.Q_REPROJ.FITS`, same as a real user's
`.fits` input always did. All bash-side filename references in
`tests/run_tests.sh` that depended on the old double-extension
behaviour (sections 30, 33-36, 38-40) were updated to match, and the
full suite re-verified at 127/127 after the change -- not just the new
T15 sections, since sections 30/33-36 (pre-existing, not part of this
ticket) also construct these tools' own output filenames and would
otherwise have silently broken. Also caught directly during manual
verification (before either was ever run inside the permanent suite): a
`rm -f *MATCHED*` cleanup command deleted `TEST_BAND2_UNMATCHED.*`
itself (the glob matches "UNMATCHED" as well as the intended
`*_MATCHED.FITS` outputs) -- recovered immediately since
`make_test_cubes.py` is fully deterministic, but the permanent test
sections use exact filenames (`${x}_MATCHED.FITS`), not a bare
substring glob, specifically to avoid this class of mistake recurring.

**Verification:** all three scenarios manually verified against real
invocations before being written into the permanent suite (not written
from assumption), each including a genuine mismatch → tool fixes it →
`check_rm_peak.py` confirms both known sources recovered at the correct
RM. Full suite: 127/127 (up from 121) after landing, and 127/127 again
after the later `strip_fits_ext` generalization and its knock-on
filename-reference updates across sections 30/33-36/38-40.

### T16 — `badchan_file` Made Genuinely Per-Band in `convolve_cubes`/`match_cubes`

**Gap found the hard way, while preparing the real WALLABY+EMU
multi-band validation run (the first genuinely large-scale exercise of
this whole toolchain -- see R1.0's own "what's next"):** WALLABY has
one real bad channel (index 1), EMU has two (indices 161, 178),
confirmed directly from each cube's own real CASA BEAMS binary table,
not assumed. The first attempt to pass this to `match_cubes` used a
single combined `badchan_file` listing all three indices -- which,
read directly from `apply_badchan_list`'s own source, applies every
listed index identically to *every* infile, bounded only by that
infile's own channel count. For WALLABY (144 channels) and EMU (288
channels) this meant EMU's own perfectly good channel 1 would also get
excluded, purely as a side effect of sharing one list across bands with
different channel numbering. This was caught and rejected before any
real run happened -- the workaround was judged unacceptable, not a
tolerable rounding error, and the fix was done properly instead of
shipped as a documented quirk.

**A second finding, this one reassuring rather than a bug:** neither
`convolve_cubes` nor `match_cubes` actually needed a manual
`badchan_file` for this specific case at all. `read_beams_table`
already auto-detects a bad channel per file, directly from that file's
own BEAMS table (BMAJ or BMIN < 1e-6 arcsec -- CASA's own ~1.18e-38
sentinel for a flagged channel), independently for every infile, with
no cross-band contamination possible. The bad channel's pixel data is
written as NaN (not garbage), the output's own BEAMS table re-flags it
the same way for any downstream tool, `reproject_cubes` passes NaN as
`astResampleR`'s own `badval` (excluded from the resampling kernel, not
blended into neighbours), and `rm_synthesis` does its own independent
per-pixel/per-channel NaN/Inf detection on the actual Q/U data. The
whole chain already handles real BEAMS-table-flagged bad channels
correctly and safely end to end with zero manual configuration.

**Why fix `badchan_file` anyway, if auto-detection already covers it:**
a CASA BEAMS table is *a* source of bad-channel information, not the
only one -- RFI or other known-bad channels may never show up as a
degenerate beam entry at all. `rm_synthesis` already had a genuine,
independent, per-band manual override for exactly this (`badchan_file`/
`global_badchan_file`, comma-separated, one file per band, count
validated against band count). `convolve_cubes`/`match_cubes` only had
the single-shared-list version -- verified directly by reading both
files' own `parse_args`/`apply_kv` (`character(len=512) :: badchan_file`,
a scalar, not `beamfiles`' own `character(len=512) :: beamfiles
(max_inputs)` array), not assumed from the `--help` text alone, since
the `--help` text itself was ambiguous on this point.

**Fix:** `badchan_file` in `src/convolve_cubes.f90` and
`src/match_cubes.f90` (independent copies, adapted not shared, same
convention as `strip_fits_ext`'s own three copies) changed from a
scalar to a comma-separated, one-entry-per-infile list -- exactly
mirroring `beamfiles`' own existing CSV-per-infile convention (same
`raw_*` staging variable, same post-parse `cfg_csv_count`/
`cfg_csv_get_item` split, same "must list exactly as many entries as
infiles" validation, empty entries allowed for an infile with no
manual list of its own). Each infile's own bad-channel reading moved
from once before the per-infile loop to inside it, keyed on that
infile's own list entry.

**Verification:** manually verified first -- two synthetic bands, each
given a distinct, deliberately different bad-channel index, confirmed
each output cube's own designated channel (and only that one) came out
as an all-NaN plane, for both `convolve_cubes` and `match_cubes`
independently. Then locked in as a permanent regression (`tests/
run_tests.sh` section 41, both tools). Full suite: 129/129 (up from
127).

### T17 — Long `infiles=` CLI Argument Silently Truncated at 512 Chars

**Gap found the hard way, running `scripts/run_pipeline.sh` on the real
WALLABY+EMU data:** the real run failed with a bizarre "failed to open
FITS file: .../match_input_symlinks/im" -- a filename cut off after 2
characters. Root-caused directly, not guessed: `run_pipeline.sh`'s own
symlink-redirection scheme (needed to land `match_cubes`' output on a
different disk than its raw input -- see T18 below) lengthens every
path by the symlink directory's own prefix; the combined `infiles=`
argument for 4 real files came to 565 characters. `parse_args`'s own
`character(len=512) :: this_arg` CLI-token buffer in `convolve_cubes.
f90`/`match_cubes.f90`/`reproject_cubes.f90` silently truncates
anything longer (Fortran raises no error on this), and the cut
happened to land exactly at "...match_input_symlinks/im" -- byte for
byte matching the error. This bug was real but dormant before this
session: every earlier manual invocation used the shorter raw
`/data1/tmp/multi-band/...` paths directly, never crossing 512 chars.

**Fix:** widened `this_arg`/`cli_val` and the `raw_infiles`/
`raw_beamfiles`/`raw_badchan_file` CSV-staging buffers (and the
equivalent cfg-file-reading buffers) from 512 to 16384 characters in
all three tools -- comfortable headroom over the theoretical worst
case (`max_inputs`=50 entries x up to 512 chars each). Per-entry
buffers (`infiles(:)`/`beamfiles(:)`/`badchan_file(:)`, each already
512 chars) were not the problem and are unchanged.

**Verification:** confirmed the exact failure mode first (a 565-char
argument reproduces the identical truncation), then confirmed the fix
with a 515-char and later a 533-char real `infiles=` argument across
all three tools -- both bands processed correctly, full paths intact
throughout. Locked in as a permanent regression (`tests/run_tests.sh`
section 42, all three tools, using a deliberately deep nested fixture
path to exceed 512 chars). Full suite: 132/132 (up from 129).

### T18 — `scripts/run_pipeline.sh`: Any Non-Empty Subset of Stages

**Gap found while setting up the real WALLABY+EMU run's hybrid disk
plan** (matched cubes on a slow, high-capacity disk alongside the raw
input; `rm_synthesis`/`rmclean_cubes` output on a fast NVMe disk):
`run_pipeline.sh` has exactly one `outdir=`, used both as the
destination for `match_cubes`' own symlink-redirected output AND for
`rm_synthesis`/`rmclean_cubes`' own output -- it cannot express two
different disks for one invocation. Compounding this, the script
previously hard-required `rmsynth` in every `stages=` list ("match
alone is not a full pipeline run -- use bin/match_cubes directly for
that"), so there was no way to run `match` in one invocation (`outdir`
on the slow disk) and `rmsynth`(`,rmclean`) in a separate, later
invocation (`outdir` on the fast disk) through the script itself.

**Fix:** removed the hard `rmsynth`-required check. Any non-empty
subset of `{match, rmsynth, rmclean}` is now valid. Whenever a stage's
own inputs would normally come from an earlier stage that is NOT
present in the SAME invocation, that stage's own cfg template's own
path-shaped keys are used exactly as written, instead of being
overridden with a chained value -- `rmsynth`'s `path=`/`infileQ=`/
`infileU=` (already behaved this way whenever `match` was absent) and,
newly, `rmclean`'s `ampfile=`/`phafile=`/`maskfile=`/`outfile=`
whenever `rmsynth` is absent (with existence checks against whatever
the template says, so a stale/missing path fails loudly and by name
rather than a bare "not found").

**A second, previously-masked bug found in the process:**
`scratch/run_rmsynthesis_test.sh`'s own cfg-value extraction
(`awk -F= '/^path=/{...}'`) was anchored to a bare `key=value` with NO
whitespace tolerance, unlike this same script's own `use_gpu`/
`dry_run` extraction a few lines away (`/^use_gpu[[:space:]]*=/`).
This never mattered before because `run_pipeline.sh` always rewrote
`path=`/`infileQ=`/`infileU=`/`outfile=` into tight `key=value` form
via its own `cfg_set_inplace` whenever `match` ran -- but this
project's own common aligned cfg style (`path                = ...`)
is exactly what every `rmsynth_cfg_template` in this repo already uses,
and a `stages=rmsynth`-alone invocation (T18's whole point) feeds such
a template straight through unmodified. Fixed to match the same
`[[:space:]]*=` tolerance already used two lines above it in the same
file.

**Verification:** the existing full-chain smoketests (`cfg/
pipeline-e2e-smalltest.cfg`, `cfg/pipeline-e2e-multiband-smalltest.
cfg`) re-run and confirmed unchanged. The new capability itself
verified with a genuine 3-way split -- `stages=match` alone (real
reprojection of a genuinely-mismatched band, correct manifest,
`outdir` on one disk) -> `stages=rmsynth` alone (template's own
`path=`/`infileQ=`/`infileU=` pointed at the match-only run's real
output, `outdir` on a different disk) -> `stages=rmclean` alone
(template's own `ampfile=`/`phafile=`/`maskfile=` pointed at the
rmsynth-only run's real output) -- all three ran correctly end to end,
each stage's real output consumed correctly by the next, run as three
fully separate `run_pipeline.sh` invocations. Full suite: 132/132
(unaffected -- this script is not part of `tests/run_tests.sh`'s own
automated coverage, verified instead via its existing manual smoketest
cfg convention).

### T19 — `io_overlap=y` Hang: CFITSIO Concurrency Hazard, Fixed with Raw Stream Writes

**Gap found running the real WALLABY+EMU `match` stage (`stages=both`)
for real:** the run stalled mid-write with `io_overlap=y` -- output file
size unchanged, every OS thread asleep (`S` state), zero CPU, zero disk
I/O, indefinitely, well after real progress had already happened (so
not a startup/config problem). Never seen on the synthetic test
fixtures because it is timing-dependent: `io_overlap=y`'s whole point
is to overlap the CURRENT block's write (on a background pthread) with
the NEXT block's read (on the main thread) -- on the small, fast
fixtures used by `tests/run_tests.sh`, one side always finishes before
the other genuinely overlaps in wall-clock time; on the real, slow
`/data1` disk and large real cubes, they do. Root cause: the background
write pthread's own CFITSIO call (`FTPSSE` on `out_unit`) and the main
thread's own next-block read (fresh CFITSIO calls, on a different
input unit, inside `read_freq_block`/`read_one_block`) genuinely
executing concurrently -- CFITSIO is not guaranteed thread-safe across
*any* two of its own calls from different OS threads at the same time,
even on different file units, without a specific reentrant build
(confirmed: this build links the ordinary, non-reentrant `libcfitsio`).
The two threads block on what is almost certainly an internal CFITSIO
lock. A real, structural hazard, not a performance question --
`io_overlap=n` (dropping the background write thread entirely) would
have made it go away, but only by removing the overlap, not the actual
unsafe pattern; rejected as a workaround for exactly that reason.

**Fix:** applied `rm_synthesis_mod.f90`'s own already-working
`io_write_threads>1` design (T6, `docs/dev/
IO_PARALLEL_OPTIMISATION_PLAN.md`) to every affected write path.
CFITSIO is used only once per output file, single-threaded, to create
the file, write every header keyword (including the CASAMBM/BEAMS
extension copy where present), and fetch the data section's own byte
offset via `FTGHAD` -- then the handle is closed immediately.
Every pixel write from that point on goes through plain Fortran stream
I/O (`open(access='stream', form='unformatted')` + `write(u,
pos=byte_pos)` at a computed byte offset, with explicit host-endianness
handling since CFITSIO's own big-endian conversion no longer runs).
CFITSIO is therefore never touched concurrently, because after header
setup it is not touched AT ALL during the write loop -- nothing to
lock, rather than a lock relied on to make concurrent access safe
(the same reasoning `write_reprojected_file`'s own block-batching
comment already used for a related, earlier concurrency concern).
Fixed in all four affected write paths:
- `convolve_cubes.f90`'s `write_convolved_file` (single fixed
  `freq_axis`, contiguous plane-per-plane byte offset --
  `write_freq_block_raw`).
- `match_cubes.f90`'s `process_one_file_restricted` (`stages=convolve`/
  `stages=both`, same single-`freq_axis` design, its own
  `write_freq_block_raw`).
- `reproject_cubes.f90`'s `write_reprojected_file` (the canonical
  general N-other-axes design -- `write_one_block_raw`).
- `match_cubes.f90`'s `process_one_file_general` (`stages=reproject`
  alone, a verbatim port of the above -- `write_one_block_raw_general`,
  reusing `process_one_file_restricted`'s own
  `host_is_big_endian_mc`/`swap_bytes_r4_inplace_mc` rather than a
  third copy).

For the two general-axis paths, the byte offset is not simply
`(chan_start-1)*nx*ny*4` (that only holds for a single fixed freq
axis): output axis 3 (`other_axes(1)`) is the one that varies within a
block, every slower axis (`other_axes(2:n_other)`) is held fixed for
the whole call, and axes 1,2 (sky) are always written at full extent.
Since FITS/Fortran arrays are both axis-1-fastest, this means the
`fixed_offset` contribution from the slower, fixed axes is a constant
per call (`sum over k=2..n_other of (other_idx(k)-1) * stride(axis
2+k)`, with `stride` accumulated as the running product of every
faster axis's own extent), and only the `chan_start` term varies block
to block -- computed once per `write_one_block_raw[_general]` call, not
re-derived per plane.

A Makefile gap was found and fixed in the process:
`$(CONVOLVE_EXECUTABLE)`'s link recipe was missing `$(SRCDIR)/
printerror.f90` (unlike `reproject_cubes`/`match_cubes`'s own recipes,
which already had it) -- pre-existing, masked because the only
`printerror` call site that existed before this fix was apparently
optimized away as unreachable; the new, genuinely-reachable `FTGHAD`
error path exposed it as a link failure.

**Not a hang risk in the first place:** `rmclean_cubes.f90`'s own
`io_overlap` already defaults its background-write path through the
same raw-stream design whenever `nwriters>1` (renamed from
`io_write_threads`, T12, `docs/dev/RMCLEAN_INTEGRATION_PLAN.md`); only
its `nwriters=1` (default) + `io_overlap=y` combination still calls
`FTPSSE` from the background thread and could in principle hit an
equivalent hazard. Out of scope for this ticket (not exercised by the
real WALLABY+EMU run, which reaches `rmclean_cubes` as a separate,
later invocation) -- noted here for whoever picks up `rmclean_cubes`
next, not fixed.

**Verification:** full suite 132/132 throughout (unchanged), including
sections 33/34/35/36's own `io_overlap=y` bit-identical-output checks
for all four write paths -- these only ever verified correctness of
the write itself (byte-for-byte identical to `io_overlap=n`), which the
old `FTPSSE`-based code already got right; they do not, and structurally
cannot, exercise the real hang (small fixtures never genuinely overlap
in wall-clock time). The actual concurrency fix's own confidence comes
from the design itself, not a timing-dependent regression test: CFITSIO
is provably never called concurrently because it is never called at
all once the write loop starts. The real, final verification is the
WALLABY+EMU run itself, relaunched with `io_overlap=y` restored (see
`cfg/pipeline-multiband-wallaby-emu-match.cfg`).

### T20 -- Second Real Deadlock: Concurrent AST Object Creation, Fixed with a Critical Section

**Gap found relaunching the real WALLABY+EMU `match` stage after T19:**
band 1 (WALLABY Q, `stages=both`) completed correctly (T19's own fix
holding up), but the run then deadlocked at the very start of band 2
(EMU Q) -- zero CPU, zero disk I/O, indefinitely, a second real hang
distinct from T19's. A live `gdb` thread dump (`thread apply all bt`
on all 6 threads, non-destructive attach/detach) showed every thread,
including the main one, parked in `pthread_mutex_lock` deep inside
`libstarlink_ast.so.9`, called from `ast_read_` <- `load_wcs` <-
`process_one_file_restricted` -- the per-thread WCS-Mapping setup step
that runs once at the start of processing each file (`stages=both`/
`stages=reproject`, wherever reprojection is part of the run;
`convolve_cubes.f90` never touches AST at all and was never exposed to
this). Two threads waited on each of two distinct internal AST
mutexes, the classic shape of a real deadlock, not just slow
contention. Never triggered on the synthetic test fixtures (a timing-
dependent race needing genuine concurrent `ast_read_` calls, which
only reliably happens with real, larger FITS headers) -- and never
triggered on band 1 either, purely by luck of thread scheduling.

**Why per-thread Mapping creation exists in the first place:** each
thread builds its own private pixel-to-pixel Mapping from scratch
(own `ast_begin` context, own `load_wcs`/`extract_sky_mapping`/
`compose_pix2pix` calls) rather than sharing one Mapping built by the
caller, because AST's own documented multi-threading model (SUN/211
Sec 4.12, "AST Objects within Multi-threaded Applications") requires
`astLock`/`astUnlock` to hand an Object from the thread that created
it to another thread -- and this Fortran binding exports neither
(confirmed directly against `libstarlink_ast.so.9`'s own symbol table:
no `ast_lock_`/`ast_unlock_`/`ast_thread_` at all, only the internal
`astLockId_`/`astUnlockId_`/`astCheckLock_`/`astThreadId_` C-level
helpers, none of them exposed as Fortran-callable wrappers).

**An initial fix attempt was tried and reverted:** build the Mapping
ONCE, single-threaded, then hand each worker thread its own
independent `ast_copy` (confirmed via a standalone test program that
`ast_copy(obj, status)` -- the inherited-status convention every other
AST call in this codebase already uses -- returns a genuinely distinct,
valid object). This is exactly what AST's own docs prescribe for
objects needed by multiple threads (SUN/211 Sec 4.13, "Copying
Objects": "a deep copy of the Object should be taken for each thread,
using astCopy"). It does not work in this Fortran binding: AST enforces
per-thread object OWNERSHIP at runtime regardless of how an object was
obtained, and the docs' own prescribed procedure requires the
`astUnlock`/`astLock` handoff around that copy -- which, per the above,
isn't available. Confirmed the hard way: every worker thread other
than the one that ran `ast_copy` itself raised `astResampler`'s own
"Invalid Object pointer given ... currently owned by another thread"
error the moment it tried to use its "own" copy. Two real regression
tests caught this immediately (T15's own grid-only-mismatch and long-
`infiles=` sections), which is exactly why the test suite is run after
every change, not just at the end.

**Actual fix:** each thread still builds its own private Mapping from
scratch, unchanged -- but the creation step itself
(`ast_begin`/`load_wcs`/`extract_sky_mapping`/`compose_pix2pix`) is now
wrapped in a named `!$omp critical (ast_setup)` section, serializing
just that one step. This directly targets the actual, confirmed root
cause (6 threads calling `ast_read_` genuinely simultaneously,
contending on `libstarlink_ast.so.9`'s own internal object-creation
bookkeeping, shared across all threads regardless of how independent
the resulting objects are) without needing any cross-thread object
sharing at all -- each thread's Mapping remains entirely its own,
squarely inside what this AST binding actually supports. The real
per-plane resampling (`ast_resampler`, already proven safe under
genuine concurrency by band 1's own successful completion) is
untouched and stays fully parallel; only the brief, one-time-per-file
setup window is now one-thread-at-a-time instead of six-at-once.
Applied to all three affected call sites: `reproject_cubes.f90`'s
`write_reprojected_file` (the canonical source), `match_cubes.f90`'s
`process_one_file_general` (`stages=reproject` alone), and
`match_cubes.f90`'s `process_one_file_restricted` (`stages=both`,
inside its own `do_reproject_l`-gated per-block setup).

**Verification:** full suite 132/132 after each of the three call
sites (reproject_cubes fixed and tested first, then match_cubes' two
call sites). The real, load-bearing verification: the actual real
WALLABY+EMU `match` run, relaunched clean from scratch with both T19's
`io_overlap=y` and this fix in place, genuinely completed band 1
(WALLABY Q, `finished:`/`OK: wrote` logged, confirmed via `gdb`/
`/proc/<pid>/io`/`iostat` that the earlier "stall" reports during this
run were real, active work, not a hang) and moved cleanly into band 2
(EMU Q) -- the exact transition that deadlocked every time before this
fix, now passed cleanly. The full real run went on to complete all 4
bands (Q/U x WALLABY/EMU), ~9h15m total, manifest written, pipeline
exit code 0 -- the definitive end-to-end confirmation of both T19 and
T20 together.

### T21 -- `nwriters`: Configurable N-Way Concurrent Block Writers for `convolve_cubes`/`match_cubes`/`reproject_cubes`

**Motivation:** while watching the real WALLABY+EMU `match` run, a live
`perf`/`/proc/<pid>/io` profile of one write-heavy window showed the
majority of sampled CPU cycles inside the block write itself (raw
stream I/O copying a multi-GB buffer into the page cache), with zero
actual bytes reaching the physical block device in the same window --
i.e. real, measured CPU cost, not disk latency, for a single giant
per-block write. `rm_synthesis_mod.f90`/`rmclean_cubes.f90` already
solved exactly this with `nwriters` (T6, renamed from
`io_write_threads`, T7/T12 above): split one block/tile's own planes
into N disjoint concurrent writers instead of one. `convolve_cubes.f90`/
`match_cubes.f90`/`reproject_cubes.f90` had no equivalent option --
only `io_overlap` (whether a block's write overlaps the NEXT block's
read+compute), never how many writers do that write once dispatched.
Explicitly scoped as **configurable, not forced**: on a genuinely
rotational disk (confirmed this machine's own `/data1` earlier this
session), concurrent writes are not obviously a win and may hurt --
the point is giving the user the choice, defaulting to the existing,
already-proven serial behaviour (`nwriters=1`).

**Design, mirrored exactly from `rm_synthesis_mod.f90`'s own
`do_tile_write`:** `nwriters_eff = max(1, min(nwriters,
omp_get_max_threads()))`, further bounded by the current block's own
`chan_len` (never more writers than planes in a single block). When
`nwriters_eff>1`, the block's planes are split into
`nwriters_eff` disjoint, contiguous chunks via the same even-split
(base+remainder) scheme `do_tile_write` already uses, each chunk
written by its own `write_freq_block_raw`/`write_one_block_raw[_general]`
call (T19's own raw-stream writers, already safe for concurrent,
disjoint-byte-range use) inside a plain nested `!$omp parallel do
num_threads(nwriters_eff)` -- reusing OpenMP's own thread-spawning
rather than managing raw pthreads directly. This composes cleanly with
`io_overlap`: when `io_overlap=y`, the whole block-write call already
runs on one background pthread (T19's own `block_write_dispatch_async`),
which itself becomes the "master" of this nested parallel region --
`io_overlap` decides WHEN a block's write runs, `nwriters` decides how
many concurrent writers do it once dispatched, and the two knobs are
otherwise independent, same as in `rm_synthesis`.

**The clamp formula itself was explicitly discussed and left
unchanged** (not switched to a spare-core-headroom formula capping at
`omp_get_num_procs() - omp_get_max_threads()`): in practice the writer
thread(s) piggyback on cores that are genuinely idle at write-dispatch
time (either the physical cores outside this machine's own
`OMP_NUM_THREADS=6` convention, or the OMP worker threads' own idle
window between one block's `!$omp end parallel`/`!$omp end do` and the
next block's first compute region, while the next block's
single-threaded read runs) -- see T7/T12's own entry for the fuller
account of why real contention is lower than a naive thread-count sum
suggests.

**Applied to:** `convolve_cubes.f90`'s `write_convolved_file`/
`do_block_write`; `match_cubes.f90`'s `process_one_file_restricted`
(`stages=convolve`/`both`) and `process_one_file_general`
(`stages=reproject` alone), each with their own `do_block_write`/
`do_block_write_general`; `reproject_cubes.f90`'s
`write_reprojected_file`/`do_block_write` (the general-axis case,
splitting along `other_axes(1)`, the axis varying within a block, with
every slower axis already fixed for the whole write job). Same
`nwriters` cfg key (CLI and `--config`), same default (1), same
`--help` text pattern, across all three -- and now consistent with
`rm_synthesis`/`rmclean_cubes`' own `nwriters` (T7/T12), closing the
naming gap that motivated the rename in the first place.

**Verification:** full suite, 4 new bit-identical regression tests
(`nwriters=4` vs `nwriters=1`, `io_overlap=y`, on the existing small
fixtures) added alongside the existing `io_overlap` ones -- one per
tool/path: `convolve_cubes` (§43), `match_cubes` `stages=convolve`
(§44), `reproject_cubes` (§45), `match_cubes` `stages=reproject`
(§46). Full suite: 136/136.

### T22 -- Per-File Stage-Timing Visibility at INFO Level, and a Genuine Convolve-vs-Reproject Split for `match_cubes`

**Gap:** watching the real WALLABY+EMU `match` run, there was no way
to tell where its own wall-clock time was actually going -- `log_level`/
`timing_enabled` defaulted to `info`/`n`, and at that level `logging_
mod.f90`'s existing `log_message('debug', 'tile_thread', ...)` per-
block/per-thread lines (the only timing detail that existed at all)
are gated behind `debug`, several times noisier than wanted for a
simple "which stage dominates" question, and were never emitted at all
for this run. Separately, `match_cubes.f90`'s `process_one_file_
restricted` (`stages=both`) wraps its ENTIRE per-plane loop -- both
`convolve_to_beam` AND (when `do_reproject_l`) `ast_resampler` -- in
one lumped `'block_convolve'` timer, so even turning on the existing
`timing_enabled` end-of-run summary couldn't have separated "time in
convolution" from "time in reprojection" for exactly the run in
question.

**Fix, two parts:**
1. **`logging_mod.f90`:** added a second, always-on accumulator
   (`file_stage_totals`, entirely separate from the existing
   `stage_totals`/`timing_enabled`-gated one, which keeps its exact
   original behaviour unchanged) plus two new subroutines,
   `timer_reset_file_stages()` (call once before starting a file) and
   `timer_report_file_summary(file_label)` (prints unconditionally at
   `[info] [timing]` level, no-op if nothing was timed). `timer_stop`
   now always updates `file_stage_totals` (one wall-clock read + one
   atomic add -- negligible next to the FFT/AST compute it's timing)
   regardless of `timing_enabled`. Wired into `convolve_cubes.f90`,
   both of `match_cubes.f90`'s write paths, and `reproject_cubes.f90`
   -- reset right before each file's `starting:` log line, reported
   right after its `finished:` line. Deliberately NOT wired into
   `rmclean_cubes.f90`: that tool processes one amp/pha image pair per
   invocation, not a loop over multiple input files, so the "per-file"
   framing doesn't apply -- it already has its own end-of-run summary
   via the existing, unchanged `timer_report_summary()`/
   `timing_enabled`.
2. **`match_cubes.f90`'s `process_one_file_restricted`:** added
   `timer_start`/`timer_stop('convolve_compute', ...)` and
   `('reproject_compute', ...)` directly around each of the file's
   three `convolve_to_beam`/`ast_resampler` call sites (`stages=
   convolve` alone, and both orderings of `stages=both` --
   `convolve_reproject`/`reproject_convolve`), using new per-thread-
   private `t_conv_local`/`t_resamp_local` start times. The existing
   outer `'block_convolve'` timer is untouched (still a real, useful
   bound on the whole per-block parallel loop's own wall-clock time) --
   this adds a finer split alongside it, not a replacement.

**Worth being honest about (not fixed, by design):**
`convolve_compute`/`reproject_compute` are summed ACROSS every
worker thread (each thread's own accumulated time, atomic-added into
one shared total) -- a genuine "which stage costs the most total
CPU-seconds" signal, directly useful for "where did the time go." The
other stages (`block_read`/`block_write`/`block_write_join`/the outer
`block_convolve`) are single-thread wall-clock measurements, taken
outside the parallel region. Both are correct on their own terms, but
they are not on the same basis -- the report's own "pct" column sums
all registered stages together as if they were, so (for example) a
well-parallelised `convolve_compute`+`reproject_compute` total can
legitimately exceed the outer `block_convolve` wall-clock figure when
several threads ran concurrently. Confirmed directly on a real test
run (6+ threads available, no `OMP_NUM_THREADS` cap set): `convolve_
compute`=0.009s + `reproject_compute`=0.132s (thread-summed) against
`block_convolve`=0.069s (wall-clock) for the exact same per-plane
loop -- not a bug, just two different units sharing one column. A
reader wanting genuine wall-clock breakdown should treat `block_*`
entries as the wall-clock bound and `*_compute` entries as "which
compute stage dominates," not literally add every row to 100%.

**Verification:** full suite 136/136, unaffected (pure additive
logging, no behaviour/output change -- confirmed via the existing
`nwriters`/`io_overlap` bit-identical tests still passing bit-for-bit).
Manually confirmed the new per-file summary appears by default (no
flags needed) on a real `stages=both` run of the small mismatched-band
fixture, correctly showing `reproject_compute` (62%) dominating
`convolve_compute` (4%) for that case.

### T23 -- Real Output Is 100% NaN: `convolve_to_beam` Has No NaN Handling (found, root-caused, fixed, verified -- committed)

**Gap found:** after the real WALLABY+EMU `match` run finally completed
end-to-end (all 4 bands, exit code 0 -- see T19/T20's own verification
notes above), a direct check of the real output (`~/venv/rmtool/bin/
python3` + astropy, per [[reference_python_venv]]) found every single
channel of every one of the 4 real output files 100% NaN -- not just
the 1-2 channels genuinely auto-detected as bad. All 4 files (both
Stokes, both bands) deleted afterward to reclaim ~268GB (2026-08-04);
none of this real output is usable.

**Root-caused, not guessed** (per [[feedback_verify_code_behavior_before_asserting]]):
1. Raw-byte inspection of the disk file (`struct.unpack`, bypassing
   astropy) showed the identical canonical big-endian quiet-NaN pattern
   (`7f c0 00 00`) repeated for every pixel, in both a genuinely-bad
   channel and a channel never flagged bad -- ruling out a byte-swap/
   endianness bug (T19's own write path): a corrupted real value would
   produce *varied* bit patterns per pixel, not one repeated canonical
   encoding. The value handed to the writer genuinely is NaN.
2. The real BEAMS table's own values (checked directly) are sane
   (~9-10 arcsec, correct `TUNIT=arcsec` convention, matching what the
   code assumes) -- not a units-mismatch bug.
3. The raw *input* data (before any processing) is genuinely clean for
   non-bad channels -- real flux values, ~7-8% NaN (real ASKAP edge-
   of-primary-beam blanking, a completely normal feature of real survey
   data), not corrupted.
4. A real-scale repro (14 real channels at the full real spatial extent,
   so `nx_pad=11664` matches exactly, built from the actual WALLABY
   file + its real BEAMS table) reproduced 100% NaN output even with
   `OMP_NUM_THREADS=1` -- ruling out any FFTW concurrent-plan-execution
   theory; the bug is deterministic, not timing-dependent.

**Actual root cause:** `convolve_to_beam` (`src/gaussft.f90`) performs
the whole convolution as a single global 2D FFT (`cimg = FFT(image)`,
multiply by the analytically-defined Gaussian-ratio kernel, inverse
FFT) with no NaN handling at all. A 2D FFT is a genuinely global
transform -- every output pixel depends on *every* input pixel -- so a
single NaN anywhere in the plane propagates through the transform and
poisons the entire output, not just nearby pixels. Every synthetic test
fixture used throughout this project (127+ tests) is fully finite, so
this was never exercised; only genuinely real ASKAP data (which is
essentially always partially blanked) triggers it. Confirmed this
predates all of today's other tickets (T19-T22) -- a latent bug in
`gaussft_mod.f90` since it was first written, unrelated to any of
today's changes.

**Design investigation (real-world precedent checked, not
re-derived from scratch):** `rm_synthesis.f90`'s own NaN handling
(`specQ(idx)/=specQ(idx)`, `rm_synthesis.f90:3955-3997`) doesn't
directly transfer -- RM synthesis is a per-pixel transform along the
frequency axis (zero spatial coupling between pixels), so masking a
bad `(pixel,channel)` only drops one term from that one pixel's own
sum. Spatial FFT convolution has no such luxury -- it mixes pixels
together by definition. Checked the actual, real-world precedent
instead: RACS-tools (`AlecThomson/RACS-tools`, `racs_tools/
convolve_uv.py`'s `convolve()`), a peer-used tool solving this exact
problem for the same class of ASKAP data with the same analytic-FFT-
Gaussian technique. Its actual approach (read directly from source,
not guessed): zero-fill NaN pixels, convolve once, then separately
convolve the binary NaN mask through the *same* kernel and re-NaN only
the output pixels where the mask-convolution shows the kernel footprint
was **entirely** filler (`mask_conv >= 1`, within float tolerance) --
partially-contaminated boundary pixels are left with the plain,
uncorrected zero-fill result, with no renormalization attempted.

**Design settled (more rigorous than the RACS-tools precedent, not
just a copy of it):** unlike RACS-tools, which only ever fully rejects
a pixel (`mask_conv >= 1`) and otherwise reports the plain, unnormalized
zero-fill result, this implementation renormalizes every partially-
contaminated pixel exactly, not just the fully-clean or fully-dirty
ones. Zero-fill NaN pixels, convolve once through the existing analytic
kernel to get `C_D`; separately build a 0/1 validity mask, convolve
it through the *same* kernel, and divide by `g_ratio` (the kernel's own
DC gain -- see this module's header comment on why it isn't unit gain)
to get `C_M`, the exact 0..1 fraction of each output pixel's kernel
WEIGHT (not a raw pixel count -- the kernel is a Gaussian, so a NaN
near its centre counts far more than one in its tail) that came from
real data. `image_out = C_D/C_M` where `C_M >= conv_nan_reject_frac`
(a fixed module-level constant, `src/gaussft.f90`, currently `0.5` --
reject if more of the kernel's weight fell on NaN than on valid
pixels), NaN otherwise. Reduces to byte-for-byte the original
NaN-agnostic computation when a plane has no NaN at all (verified: full
139-test suite, all synthetic fixtures are NaN-free, zero regressions).

**Implemented** in `src/gaussft.f90`'s `convolve_to_beam` only --
self-contained inside `gaussft_mod.f90`, so `convolve_cubes.f90` and
both of `match_cubes.f90`'s call sites (convolve-first and
reproject-first paths) get the fix automatically, no caller changes
needed. NaN detection via `ieee_is_nan` (`ieee_arithmetic`, the same
intrinsic already used elsewhere in this project, e.g.
`convolve_cubes.f90`'s own bad-channel writer).

**Verification:** new regression test (`tests/run_tests.sh` §48) --
injects real NaN into a 6x6 spatial corner block (all 200 channels) of
band 1's Q/U test cubes, far from both known sources (`src_A` at
(12,10), `src_B` at (22,20)), then runs the full `convolve_cubes` ->
`rm_synthesis` pipeline (mirroring §39's own structure). Confirms three
things a plane-poisoning regression would each independently break:
(1) output is NOT mostly NaN (measured 2.3%, vs. 100% before this fix
-- the direct regression check for T23's own bug), (2) pixels well
inside the injected block correctly remain NaN (the threshold rejection
is actually engaging, not just accidentally producing *some* non-NaN
output), (3) pixels far from the block are completely clean, and both
`src_A`/`src_B` RM peaks are still recovered downstream through
`rm_synthesis` -- proving the fix is genuinely LOCAL, not merely
"less broken". Full suite: 142/142.

### T24 -- `dry_run`: Disk-Type-Aware `io_overlap`/`nwriters` Advisory for `convolve_cubes`/`match_cubes`/`reproject_cubes`

**Motivation:** raised alongside T21's own `nwriters` work -- given
`nwriters`/`io_overlap`'s correct setting genuinely depends on the
target disk's own physical characteristics (a real, confirmed-rotational
disk vs. SSD/NVMe -- this machine's own `/data1` vs. `/home` being the
concrete example throughout this session), and given `rm_synthesis`
already has an established `dry_run` precedent (planning-only pass,
writes a suggested `tile_autotune.cfg`, touches no real data) for a
different kind of advisory (tile/memory/VRAM sizing) -- extend the same
convention to these three tools for I/O-parallelism settings, rather
than expecting every user to independently rediscover this session's
own disk-type reasoning (T7/T12/T21's own "why nwriters is safe on
spare cores, not free of contention on a spinning disk" discussion).

**Design:** `dry_run` cfg key (same name as `rm_synthesis`' own),
default `n`. When `y`: checks the first `infiles=` entry's own target
disk via `/sys/block/<dev>/queue/rotational`, resolved from the path's
own mount point -- the exact same, already-verified mechanism (`df`
--output=source` + `basename` + NVMe `pN`-suffix-aware parent-device
resolution) this session used interactively to confirm `/data1` is
genuinely rotational and `/home` is NVMe, reused here rather than a
parallel, untested reimplementation. Shells out via
`execute_command_line` (Fortran has no direct block-device-resolution
API) with its stdout captured to a plain file in the current working
directory -- never `/tmp` or any other system directory, applying this
project's own standing scratch-file convention to the tool's own
runtime behaviour too, since HPC compute nodes often restrict or omit
`/tmp` entirely. Prints the detected disk type and a suggestion
(`io_overlap=n, nwriters=1` for spinning; `io_overlap=y, nwriters=2`
for SSD/NVMe), and writes a `<tool_name>_dryrun.cfg` the user can copy
values from or pass directly on the command line -- mirroring
`rm_synthesis`' own `tile_autotune.cfg` shape. Explicitly advisory
only: does not change `nwriters`' own clamp formula (T7/T12's own
already-settled decision), just suggests a starting value within it.
Processes no real file when `dry_run=y` -- returns immediately after
writing the suggestion.

Implemented identically (adapt, don't share, per this project's own
module convention) in `convolve_cubes.f90`, `match_cubes.f90`, and
`reproject_cubes.f90` -- `reproject_cubes.f90`'s own CLI-then-cfg-then-
override dual-parse structure needed the extra `cli_dry_run`/
`cli_seen_dry_run` plumbing the other two tools' shared `apply_kv`
design doesn't.

**Verification:** manually confirmed against real ground truth on this
machine -- `/data1/tmp/multi-band/...` correctly reports "spinning
(rotational=1)"/suggests `io_overlap=n nwriters=1`; a path under
`scratch/` (NVMe) correctly reports "non-rotational (rotational=0)"/
suggests `io_overlap=y nwriters=2` -- for all three tools. New
regression test (`tests/run_tests.sh` §47) locks in that `dry_run=y`
writes a valid suggested cfg (`io_overlap=`/`nwriters=` present) and
genuinely processes no data (checked: the real output file that a live
run would have produced does not exist afterward). Full suite:
139/139.

### T25 -- Explicit `target_bmaj`/`target_bmin`/`target_bpa` Never Validated as a Real Deconvolution Target (found during T23's own investigation, not started)

**Gap found** while confirming a premise for T23's own NaN-handling
design (that `convolve_to_beam`'s kernel is always positive/well-
behaved in real space, which the design relies on): checked both
`convolve_cubes.f90` and `match_cubes.f90` directly for their own
`have_target` (explicit `target_bmaj=`/`target_bmin=`/`target_bpa=`)
branch, and neither calls `deconvolve_is_valid`
(`commonbeam.f90:578`, the ported MIRIAD "gaupar" validity check) on
it at all -- unlike the auto-derived common-beam path
(`find_common_beam`), which validates every candidate before ever
returning it. So a user-supplied explicit target beam that isn't
actually deconvolvable from some channel's real native beam (target
smaller than native in some direction, or otherwise invalid) currently
passes through with no check, and `convolve_to_beam`'s own kernel
exponent (`g_arg - dg_arg`) can go *positive* instead of negative at
high frequencies for that channel -- turning the intended low-pass
smoothing kernel into an unbounded high-frequency amplifier, which
would silently produce garbage (most likely Inf/NaN) output with no
warning at all.

**Not the cause of T23's own real-run NaN bug** -- confirmed directly:
the real WALLABY+EMU run's own log ("Derived common beam from 858 good
channels...") shows it went through the auto-derived (validated) path,
not an explicit target -- these are two genuinely separate gaps, found
back to back but with independent causes and independent fixes.

**Fix (not started):** call `deconvolve_is_valid` on the explicit
target beam against every input channel's own native beam, in both
`convolve_cubes.f90` and `match_cubes.f90`'s `have_target` branch,
matching exactly what `find_common_beam` already does for the
auto-derived path -- fail loudly (a real `ERROR:`, not silent garbage)
naming the specific offending channel/file, before any FFT ever runs.

### T26 -- `mem_frac_ram` Block-Sizing Ignored `convolve_to_beam`'s Own Transient Memory, Causing Real Swap on a 62.7GB Host (found, root-caused, fixed, verified -- committed)

**Gap found:** while investigating whether this workload could be
ported to a RAM-constrained host (a Raspberry Pi 5, 8GB), computed
`convolve_to_beam`'s own peak per-plane working memory directly from
its actual allocate/deallocate sequence (`src/gaussft.f90`) at
WALLABY's real resolution (`nx=11559, ny=11655`, FFT-padded to
`11664x11664` -- both taken directly from the live real run's own
log) -- came to **~7.05GB for a single plane, single thread**. Checked
this against the real, currently-running WALLABY+EMU `match_cubes`
process's actual memory (`/proc/<pid>/status`, `free -h`) rather than
trust the arithmetic alone: `VmPeak=71.65GB` against this machine's
**62.7GB** of physical RAM, with `VmSwap=3.47GB` for the process alone
and **12GB of system-wide swap in active use** -- the run has been
surviving by swapping, not by actually fitting, the entire time.

**Root cause:** `bytes_per_plane` (the formula `block_planes` --
`mem_frac_ram`'s own block-sizing knob -- is derived from, in both
`match_cubes.f90:process_one_file_restricted` and
`convolve_cubes.f90:write_convolved_file`) only ever budgeted the raw
single-precision I/O block buffer (`4 * (nx*ny + 2*nx_out*ny_out)`
bytes). It had **no knowledge at all** of `convolve_to_beam`'s own
NaN-aware-path (T23) transient working set -- a native-size NaN mask,
two padded-size `complex(dp)` FFT buffers alive at once, two
native-size `real(dp)` buffers, plus the caller's own native-size
`plane_native` held for the whole call. That transient cost is
**per OpenMP thread** (all thread-private locals, not shared) and is
entirely separate from, and additive to, the I/O block buffer -- so no
`mem_frac_ram` setting protected against it; the budget calculation
was simply blind to roughly half (or more, at high thread counts) of
the real peak memory a run would actually touch.

**Fix:** derived the exact peak transient cost by tracing every
allocate/deallocate in `convolve_to_beam` (`gaussft.f90`): `20*(nx*ny)
+ 32*(nx_pad*ny_pad)` bytes per thread, where `nx*ny` is the native
grid convolution actually runs on (`nx_in`/`ny_in`, or `nx_out`/`ny_out`
for `order=reproject_convolve`) and `nx_pad*ny_pad` is its FFT-padded
size (`next_fast_fft_size`). In both `match_cubes.f90` and
`convolve_cubes.f90`, this is now computed *before* `bytes_per_plane`
is sized, `omp_get_max_threads()` copies of it are reserved out of the
`mem_frac_ram` budget first, and `block_planes` is sized from whatever
budget remains -- clamping to 1 plane (with a clear `WARNING:`) if the
per-thread transient cost alone already exceeds the budget, rather
than silently proceeding to overcommit RAM regardless of
`mem_frac_ram`. `match_cubes.f90`'s own `conv_nx`/`conv_ny` (and the
`nx_pad`/`ny_pad` they imply) determination was moved earlier in the
subroutine, ahead of the memory budget that now needs it (previously
computed just before `plan_convolution`, well after the old budget
calculation) -- `next_fast_fft_size` is cheap and pure, so computing it
once here and again (identically) inside `plan_convolution` moments
later is harmless redundancy, not a behaviour change.

**Verification:** full suite unchanged, 142/142 (no test fixture is
large enough to exercise this path differently -- the fix only changes
behaviour once the transient cost is large relative to the
`mem_frac_ram` budget, which none of the small synthetic fixtures
trigger). Directly computed what the fix now produces for the real
run's own actual settings (`mem_frac_ram=0.25`, `OMP_NUM_THREADS=6`,
WALLABY's real dimensions): reserved transient cost for 6 threads
(42.29GB) already exceeds the entire 16.64GB budget -- exactly
reproducing the real overshoot that was directly measured via
`VmPeak`/`VmSwap` above, confirming the fix catches the actual
observed failure mode, not just a theoretical one. Same computation
shows the actual safe range on this 62.7GB host at `mem_frac_ram=0.25`
is only 1-2 threads before the transient cost alone stops fitting --
useful, concrete guidance for both this machine's own future real runs
and for the original Raspberry Pi 5 (8GB) sizing question that
prompted this investigation.

**Refined (same ticket, follow-up):** the original fix above was
all-or-nothing -- reserve `omp_get_max_threads()` copies of the
transient cost, and if that doesn't fit, floor straight to 1
plane/serial, even when some smaller thread count in between (e.g. 3
of the 6 requested) would genuinely have fit. Replaced with a plain
decrementing search in both files: start at `omp_get_max_threads()`,
step down while that many threads' own reserved transient cost is
`>=` the `mem_frac_ram` budget, floor at 1. The result
(`omp_threads_cap`) is what actually caps real OpenMP concurrency now
-- `nthreads = min(omp_threads_cap, block_planes)` (previously
`min(omp_get_max_threads(), block_planes)`), which feeds directly into
`!$omp parallel num_threads(nthreads)`, so this is a genuine cap on
threads actually spawned, not just an I/O-buffer sizing knob. Re-ran
the real numbers with T27's own (now smaller, single-precision)
transient cost: on this 62.7GB host at `mem_frac_ram=0.25`,
`omp_threads_cap` comes out to **3** for WALLABY's own dimensions and
**2** for EMU's (EMU is natively larger, `14300x12395` vs WALLABY's
`11559x11655`, so its own transient cost per thread is higher and its
safe thread count lower) -- confirmed these differ per file, not a
single fixed number, by construction (see the "known gap" note below
for why that per-file correctness doesn't yet extend to an *upfront*
multi-file estimate). Full suite re-verified after this refinement,
142/142.

**On the formula's own two `16`s (raised directly, worth being precise
about rather than treating as one coefficient):** `16*(nx*ny) +
16*(nx_pad*ny_pad)` is two DIFFERENT sums that happen to both equal
16, not one halved-and-reused number. Native-size term: `nan_mask`
(4B, `logical`) + `c_d`-or-`c_m` (4B, now `real(sp)` per T27) + the
CALLER's own `plane_native`/`plane_in` (8B, `real(dp)` -- deliberately
NOT converted, so `convolve_to_beam`'s external interface didn't need
to change) = 16. Padded-size term: `g_final` (8B, now `complex(sp)`) +
`cimg`-or-`cmsk` (8B, now `complex(sp)`) = 16. The padded term
genuinely halved (32->16, both contributors converted to single
precision); the native term only partly halved (20->16), held back by
`plane_native` staying double precision -- converting that too (a
caller-side change, not just `gaussft_mod`'s own internals) would drop
the native term further, to 8. Not implemented.

**Known gap, not yet addressed (raised directly by the user, worth
recording precisely rather than glossing over):** `block_planes`/
`omp_threads_cap` are computed **live, once per file**, at the exact
moment `process_one_file_restricted`/`write_convolved_file` begins
that file -- confirmed via source that this uses each file's own true
`naxes_f(i,:)` (populated during this tool's own per-file header
pre-scan), not the reference file's dimensions and not just the first
file in `infiles=` -- so the *live* calculation is genuinely per-file
correct, not a "look at one file and assume" shortcut. What's
genuinely missing: there is no **upfront, whole-batch preview** --
nothing surveys every file in `infiles=` before real processing starts
and reports (or plans around) the per-file thread/memory picture for
the WHOLE run in advance. A user only discovers a given file's own
`omp_threads_cap` (and whether it had to be reduced) via the warning
printed right as that file begins -- potentially hours into a run, for
a file near the end of the processing order. `dry_run=y` (T24) does
not help here either: it only checks the *disk type* of the first
`infiles=` entry for its own `io_overlap`/`nwriters` suggestion, and
never touches the convolve memory/thread question at all. A genuinely
thorough prep step would read every file's header (NAXIS-only, cheap --
same pattern already used for the reference file's own WCS-only read)
during `dry_run=y`, compute each file's own `transient_bytes_per_thread`/
`omp_threads_cap`, and report the full per-file breakdown up front --
deliberately NOT collapsed to one worst-case number for the whole
batch, since that would needlessly throttle the smaller files (WALLABY
at 3 threads) down to the biggest file's own more conservative number
(EMU at 2). Not started.

### T27 -- `gaussft_mod`'s FFT Buffers/Intermediates Now Single Precision Internally (implemented, verified -- committed)

**Motivation:** raised alongside T26's own investigation, while
assessing how to make this workload fit a RAM-constrained host (the
Raspberry Pi 5, 8GB, that prompted T26). `convolve_to_beam`'s NaN-aware
path (T23) carried every array -- the FFT buffers, the validity mask,
the intermediate data/mask-convolution results -- at `real(dp)`/
`complex(dp)`, despite FITS output already being written at
`BITPIX=-32` (`FTPHPR(..., -32, ...)`, both `write_convolved_file` and
`write_matched_file`) -- there was never a real dynamic-range/precision
argument for double precision flux data in the first place, and it
was costing roughly double the memory of the alternative for no
accuracy the pipeline's own output actually keeps.

**Design:** `gaussft_mod`'s FFT buffers (`cimg`, `cmsk`, `g_final`) and
intermediate real arrays (`c_d`, `c_m`) are now `complex(sp)`/`real(sp)`
(`sp => real32`) internally, executed via FFTW's single-precision
legacy Fortran entry points (`sfftw_plan_dft_2d`/`sfftw_execute_dft`/
`sfftw_destroy_plan`, from `libfftw3f` -- a genuinely separate library
from `libfftw3`, confirmed present on this system via `nm -D` before
committing to the plan, and added to `Makefile`'s shared `FFTW_LIBS`).
One deliberate exception, not swept into the blanket precision
reduction: the kernel exponent itself (`g_arg - dg_arg`) is a
DIFFERENCE of two potentially large terms before `exp()` -- exactly
the shape of computation vulnerable to catastrophic cancellation in
single precision -- so `g_arg`/`dg_arg` and their combination stay
`real(dp)` throughout, with only the finished kernel value cast down
to `sp` at the point `g_final` is actually stored. `convolve_to_beam`'s
own external interface (`image`/`image_out`, both dummy arguments)
deliberately stays `real(dp)` -- the precision reduction is entirely
internal to this module, so neither `match_cubes.f90` nor
`convolve_cubes.f90` needed a single line changed to pick this up.

**Verified peak-memory effect** (WALLABY's real dimensions, same
method as T26): per-thread transient cost dropped from **7.05GB to
4.33GB** (-38.5%). On this 62.7GB host at `mem_frac_ram=0.25`, that
moves the actual safe thread count from ~1-2 (T26's own finding,
pre-T27) to ~3 -- a real, usable improvement, not just a memory-safety
margin.

**T26's own budget formula had to be updated too, not left stale:**
its hardcoded byte multipliers (`20*(nx*ny) + 32*(nx_pad*ny_pad)`)
were derived from `gaussft_mod`'s PRE-T27 double-precision array
sizes. Left unchanged, T26's own budget would have kept reserving
memory for a cost that no longer existed -- safe (over-conservative),
but silently wasteful, working against the very optimisation T27 was
for. Re-derived and updated in both `match_cubes.f90` and
`convolve_cubes.f90` to `16*(nx*ny) + 16*(nx_pad*ny_pad)` (see either
file's own comment for the coefficient derivation) -- flagged there
explicitly as NOT computed from `gaussft_mod`'s own type declarations,
so it will drift stale again if those kinds ever change a third time,
unless re-derived by hand alongside any future change.

**Verification:** full suite, 142/142, including one real, EXPECTED
tolerance update rather than a silent pass: `tests/test_gaussft_padding
.f90`'s "target beam == native beam is a no-op" check compares an FFT
round-trip to `image` bit-for-bit up to a tolerance, previously
`1.0e-9` (calibrated for the old double-precision FFT) -- now measures
`~3.1e-8`, well within single precision's own `~1.19e-7` epsilon and
NOT achievable at the old tolerance by design, so the tolerance was
loosened to `1.0e-6` (comfortably above the observed error, still tight
enough to catch a real fault) with a comment explaining why. The
OTHER check in the same test (`auto-pad-inside-convolve_to_beam
matches manual pre-padding`, comparing two invocations of the same new
code against each other rather than against a fixed absolute
tolerance) still matches to the bit (`max|diff|=0.000E+00`) --
confirming the padding logic itself is untouched, only the absolute
precision ceiling changed, exactly as intended.

**Not yet done (assessed, not implemented):** computing `g_final`
(and/or the NaN mask) analytically per-pixel inline, instead of
materialising them as full padded-size arrays, as a further, OPTIONAL
low-RAM code path -- discussed as worthwhile specifically for
severely memory-constrained hosts (the original Pi 5 motivation) but
deliberately not implemented alongside T27, since it would duplicate
`convolve_to_beam`'s own per-pixel loop under a runtime/compile-time
switch rather than being a drop-in replacement. Left as a distinct,
separately-scoped future ticket.

### T28 -- Convolve Parallelization Redesign: Serial-Over-Planes, Threaded-Within-Plane (implemented, verified -- committed)

**Motivation:** raised directly by the user while reviewing T26/T27's
own numbers. The pre-T28 model parallelized `convolve_to_beam` by
running N OpenMP threads *concurrently*, each independently processing
its own whole plane -- so peak transient memory scales linearly with
thread count (T26's own `omp_threads_cap` exists purely to cap that
multiplication). The user's question: why parallelize *across* planes
at all, when the actual bottlenecks for this workload are I/O and RAM,
not FFT compute? A plane's own 2D FFT can instead be parallelized
*internally* (FFTW natively supports multi-threaded execution of a
single transform), while planes themselves are processed one at a
time -- trading "N threads x 1 plane's memory each" for "1 thread's
worth of memory, N-way-threaded". This does not need linear FFT
speedup to be worthwhile, since FFT is one step of the whole pipeline,
not the dominant cost. It also directly serves the original
Raspberry-Pi-5 portability motivation behind T26/T27: a design whose
peak memory no longer multiplies by core count is a genuine, real
claim about running this tool on an 8GB single-board machine, not just
an aspiration.

**Design:** `plan_convolution` (`gaussft.f90`) takes a new `nthreads`
argument and calls FFTW's `sfftw_init_threads`/
`sfftw_plan_with_nthreads(nthreads)` before planning, so the resulting
plan itself executes multi-threaded on every subsequent
`sfftw_execute_dft` call -- linked against `libfftw3f_omp`
(OpenMP-backed), chosen over the pthreads-backed `libfftw3f_threads`
alternative purely for consistency with the rest of the codebase's
`-fopenmp` model; both expose identical `sfftw_*` legacy Fortran
symbols (confirmed via `nm -D` before committing to the choice).
`convolve_to_beam`'s own internals (kernel build, NaN masking, array
-syntax elementwise ops) are now parallelized in-place with
`!$omp parallel do collapse(2)`/`!$omp workshare` rather than relying
on the caller to run multiple concurrent instances. In
`match_cubes.f90`/`convolve_cubes.f90`, the per-plane loop that used to
be `!$omp parallel ... !$omp do schedule(dynamic)` (N planes at once,
each single-threaded) is now a plain serial `do` loop, one
`convolve_to_beam` call at a time, each internally using
`omp_threads_cap` threads. `transient_bytes_per_thread` (T26/T27) is
renamed `transient_bytes_per_plane` and reserved ONCE regardless of
thread count, since only one plane's transient memory is ever alive at
a time now; `omp_threads_cap` reverts from T26's memory-safety search
to a plain performance choice (`omp_get_max_threads()`), since more
threads no longer costs more memory.

**Cross-check performed before implementing (per explicit instruction
to document, then cross-check, then implement):** two issues found,
both resolved before writing code, not discovered mid-implementation:

1. FFTW's `sfftw_plan_with_nthreads` and the old
   concurrent-multiple-threads-each-calling-`sfftw_execute_dft` pattern
   are mutually exclusive -- combining them (a plan already configured
   for N-way internal threading, itself invoked from N concurrent
   OpenMP threads) would oversubscribe by a factor of N. Resolution:
   full replacement of the old parallel-over-planes loop, not a
   layering of the two strategies.
2. A naive single-serial-pass redesign covering BOTH convolve and
   reproject together would have been a real regression for
   reproject-heavy files: `ast_resampler` (via `reproject_cubes`'s own
   `reproject_compute` stage) has no internal threading of its own, so
   reproject's existing parallel-over-planes strategy is the ONLY
   source of its parallelism -- serializing it alongside convolve would
   have silently made reproject-heavy runs slower, trading a real
   existing speedup for one that doesn't apply to that stage.
   Resolution: kept as two independent passes over each block's
   planes, not one combined loop -- convolve always
   serial-over-planes/threaded-within-plane (T28's new strategy,
   regardless of whether it runs first or second), reproject always
   parallel-over-planes exactly as before (T20's critical-section
   AST setup, unchanged). Which pass runs first still follows the
   existing `convolve_first_l`/`order=` config, matching the
   pre-T28 ordering semantics exactly -- only the internal
   parallelization strategy of each stage changed, not the pipeline's
   observable behaviour.

**Implementation notes:**
- `Makefile`'s `FFTW_LIBS` gained `-lfftw3f_omp`; `tests/run_tests.sh`'s
  standalone `gaussft_padding`/`gaussft_threading` build commands
  needed the same library added separately (missed once, caught by a
  link failure, fixed).
- T22's per-file stage timers (`convolve_compute`/`reproject_compute`)
  had to be re-added at all 5 call sites after being dropped during the
  rewrite -- a self-caught regression, not something the two-pass
  split intentionally changed.
- A genuine bug was found (not by the user -- via a failing
  integration test, `match_cubes stages=both order=convolve_reproject`)
  and fixed during implementation: the original single combined
  per-plane loop initialized `t_status = 0` once before AST's
  "inherited status" convention was first used
  (`call ast_begin(t_status)`; AST no-ops every call once its status
  argument is non-zero, so a garbage initial value silently no-ops the
  entire WCS/mapping setup without raising any visible error). Splitting
  the combined loop into three separate `!$omp parallel` regions (T28's
  two-pass structure) dropped this initialization from the two branches
  that do reprojection, leaving `t_status` as uninitialized
  thread-private stack garbage on entry -- confirmed directly via a
  temporary debug print, which showed real AST error code `32` (and,
  on one thread, uninitialized garbage `-311422384`) at the critical
  section, despite both passes' own debug logs showing every plane
  "done". Fixed by restoring `t_status = 0` immediately before each
  branch's own `!$omp critical (ast_setup)` block.

**Verification:** full suite, 143/143 (142 pre-existing + one new:
`tests/test_gaussft_threading.f90`, wired in as section 49). The new
test checks, against a 401x401 synthetic plane exercising the
NaN-handling (`has_nan=.true.`) path: (1) `nthreads=4` output matches
`nthreads=1` (max|diff|=0, i.e. bit-identical, not merely within
tolerance); (2) `nthreads=4` run twice against the same FFTW plan is
bit-identical (confirms the threaded decomposition is genuinely
deterministic run-to-run, not just claimed so by FFTW's own
documentation); (3) output is not degenerate (NaN block stays NaN,
finite pixels present with a real positive peak). The
`match_cubes stages=both order=convolve_reproject` scenario that
surfaced the `t_status` bug above was independently re-verified by
direct reproduction outside the test harness after the fix (exit code
0, all four synthetic multi-band files processed and written).

**TODO -- real-scale re-run of `match_cubes` on the WALLABY+EMU data,
post-T28 (not started, deferred for disk space):** everything verified
above is against synthetic fixtures only. The existing real matched
output under `/data1/tmp/pipeline_multiband_e2e/match_input_symlinks/`
(`*_MATCHED.FITS`, ~267.8GB across the 4 Q/U files) was produced on
2026-08-04, BEFORE T28 was implemented -- it reflects the old
concurrent-multi-thread-per-plane convolve strategy, not the new
serial-per-plane/threaded-within-plane one. Need to re-run
`match_cubes` on the same real WALLABY+EMU input under the current
(T28) code and compare against that existing output: correctness
(should match the pre-T28 result numerically -- T28 changed
parallelization strategy, not the convolution math itself) and real
wall-clock/peak-memory behaviour (the actual point of T28 -- does peak
RSS on this run now stay flat regardless of `OMP_NUM_THREADS`, as the
synthetic-fixture testing implies but has never been checked at this
real 563GB-input scale). Blocked purely on disk space: `/data1` is at
67% (1.2TB free) but the existing pre-T28 `*_MATCHED.FITS` output would
need to coexist with a fresh post-T28 run's own output for the
comparison to mean anything, and `/home`'s own NVMe scratch is already
at 83% (131GB free, currently held by the in-progress real
`rmclean_cubes` run's own pre-allocated output cubes). Deferred until
some of that space is cleared -- not forgotten, not silently dropped.

### T29 -- Multi-Band Merge Implicitly Assumes Bands Are Already Flux-Cross-Calibrated; No Stokes-I Use in the Multi-Band Path (documented, not started)

**Raised by the user, verified against the actual code rather than
taken on trust** ("check and verify, do not believe me since I am a
human with finite memory"): does anything in the multi-band toolchain
use Stokes-I, and does anything cross-calibrate Q/U flux scale between
bands before `rm_synthesis` merges their channels into one run?

**Verified: no, on both counts.**
- `reproject_cubes`/`convolve_cubes`/`match_cubes` -- the entire
  multi-band preprocessing toolchain -- never reference Stokes-I at
  all (grepped `src/{reproject,convolve,match}_cubes.f90`: zero hits
  for `infileI`/any Stokes-I concept). They operate on Q/U only.
- `rm_synthesis` DOES have an existing Stokes-I option
  (`infileI`/`path_I` cfg keys, `rm_synthesis_mod.f90`), but it is used
  for exactly ONE purpose: `remove_qu_bias`, a linear leakage-bias
  correction (`resiQ`/`slopeQ`/`resiU`/`slopeU`, i.e. `Q_corrected = Q
  - (resiQ + slopeQ*I)` and the U equivalent -- values supplied
  directly in the cfg, not derived from I automatically) for
  instrumental Q/U-vs-I leakage, NOT a flux-scale renormalization.
- That one existing use is explicitly BLOCKED for multi-band runs:
  `rm_synthesis.f90`'s own T2-scope-narrowing block (`if
  (n_bands_t2.gt.1)`) hard-stops with `'ERROR: multi-band Q/U bias
  correction is not yet implemented.'` the instant `remove_qu_bias=y`
  is combined with `nbands>1`.
- Net effect: for any actual multi-band run today, Stokes-I is never
  read, never touched, anywhere in the pipeline. `rm_synthesis`
  combines each band's Q/U channels directly, with an entirely
  IMPLICIT assumption that every band's own flux scale is already
  mutually consistent (correctly absolute-calibrated relative to each
  other) before they're combined -- nothing in the code checks,
  verifies, warns about, or corrects for a cross-band flux-scale
  mismatch.

**Future direction (user's own idea, explicitly NOT planned or scoped
yet -- documented here so it isn't lost, not committed to):**
per-pixel Stokes-I (a smoothed version of each band's own I cube, per
pixel) could in principle be used to cross-calibrate Q/U flux scale
BETWEEN bands before merging -- i.e., use each band's own I as a
per-pixel flux-scale reference, rather than trusting the bands were
already calibrated onto a common scale. This itself carries its own
further assumption, worth stating precisely rather than glossing over:
it only works if each band's own Stokes-I is ITSELF absolute-flux-
calibrated (correct in an absolute sense, not just internally
consistent within that band) -- if a band's own I calibration is off,
using it as a Q/U flux-scale reference would just propagate that same
error into the "corrected" Q/U. Not started; no design work done
beyond capturing the idea and its own precondition here.
