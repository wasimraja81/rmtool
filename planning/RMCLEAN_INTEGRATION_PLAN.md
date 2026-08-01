# RM-CLEAN Integration Plan

Branch: `rmclean-integration` (from `develop`)

**Status: T0 (MASK.CUBE.FITS per-channel frequency/band table), T1
(RM-CLEAN core algorithm module, `src/rmclean.f90`), T2 (standalone
tool `rmclean_cubes` + per-pixel mask-pattern caching + OpenMP
parallelism), T2b (`rm_synthesis`'s own `lsq_ref_mode`), T3
(model-based matched-filter peak refinement, replacing
`peak_interp_parabolic`'s own complex-value step), T3b (Gate 0 recast
as an RM-resolution criterion, `min_samples_per_fwhm`), T3c (tiered
fast-path refinement with data-driven `refine_nsigma` escalation --
superseding the originally-sketched coarse-to-fine idea, and closing
T3's own performance gap), and T4 (memory-budgeted, threaded block I/O
-- `plan_rmclean_tile`/`io_read_threads`/`io_write_threads`/
`io_overlap`, the same scheme `rm_synthesis` already uses in
production) all done and verified -- see each ticket's
own Evidence section below. GPU support for RM-CLEAN is explicitly
deferred to a later, separate effort (decision 6 below); everything
here is OpenMP-CPU-only by design. T7 (CLEAN convergence/stop-reason
logging + a reusable subimage-based gain-tuning workflow) done and
verified; its own gain-sweep findings are a real, actionable result
(see T7's own Evidence section) about the OLD stopping-criterion
mechanism -- superseded, not contradicted, by T8. T8 (CLEAN
stopping-criteria redesign: `abs_flux_floor`/`auto_nsigma`/`niter`,
per-pixel RM-tail sigma computed once, exact `stop_reason` logging) done
and verified -- the old `threshold=`/`threshold_snr=` pair funneled
into one overloaded, always-multiplicative `thresh` that silently
broke the absolute-flux mode and made the auto/n-sigma mode a
self-referential moving target; fixing it took this same T7 subimage's
tail RMS from 6.89 uJy (3.7x below the true ~25.5 uJy noise floor) to
25.24 uJy (matching it), converging in a mean of 3.89 iterations instead
of 233.7. `threshold=`/`threshold_snr=`/`noise_nlos=`/`noise_seed=`/
`tail_exclude_nfwhm=`/(the T9-only, now ALSO retired) `noise_percentile=`
are ALL RETIRED (not aliased) -- use `abs_flux_floor=`/`auto_nsigma=`
in any cfg going forward, nothing else. T9 (same day) fixed a second,
independent bug in `auto_nsigma`: `- avg_abs` was subtracted from
`peak_val` before comparison, inherited unexamined from the
pre-existing algorithm and itself a self-referential moving target
(recomputed from the shrinking residual every iteration) that stopped
CLEAN early on large-scale coherent structures -- `auto_nsigma` now
compares `peak_val` directly, no subtraction. T9 ALSO tried a per-pixel
lowest-percentile-fraction sigma estimate (`noise_percentile=`,
Rayleigh-IQR-corrected) as a safer alternative to T8's own peak-
relative exclusion window, but then measured it directly against real
data and found it performed WORSE than simply using the pixel's own
FULL dirty spectrum's IQR (a fixed Rayleigh factor, `IQR/0.90656`, no
subsetting) -- so `noise_percentile=` was retired again, same day it
was introduced, in favour of the simpler full-spectrum estimator
(`estimate_iqr_sigma`). cfg defaults finalised: both criteria off by
default (must opt in explicitly); `auto_nsigma`'s own default value is
`1.0` (not `5.0` -- justified now that the sigma estimate itself is
trustworthy); `abs_flux_floor` values are always given with an explicit
unit suffix; `cfg/rmclean-jennifer.e2e.cfg` sets a REAL, active
`abs_flux_floor=20uJy` (anchored to this dataset's own independently-
measured ~25.5 uJy/beam floor) alongside `auto_nsigma=1.0` as a genuine
safety net, not a placeholder. See T9's own Evidence section (including
its own UPDATE) for the full mechanism, the accuracy comparison that
drove the full-spectrum-over-percentile decision, and per-cfg
justification.**

## Context

The user has a working Fortran RM-CLEAN implementation from their own
thesis (source at `~/softwares/CURR_DEVEL/RM_CLEAN_TESTS/`, fixed-form
`.f`), used for their thesis but never integrated into `rmtool`. The goal
is to port the relevant algorithmic core into `rmtool` as a new option
(the user's framing: "integrate it to rmsynthesis as an option"),
modernized (free-form Fortran, `rmtool`'s own module/kind conventions)
and made robust, not a raw port. Thesis reference for the deferred
bandwidth-depolarization (BW-depol) work: `~/Thesis/wasim_thesis_Final.pdf`,
chapter 2.5.1 "Bandwidth depolarization: resultant modification of the
complex polarized spectra" (printed pp. 25-30, PDF pages 33-38 — a
consistent +8 offset between printed and PDF page numbers in that file).

## Source code inventory (`~/softwares/CURR_DEVEL/RM_CLEAN_TESTS/`)

Read in full before any design discussion. Files to actually port:
`rm_clean.f` (core Högbom-style complex clean), `compute_dirty_rmbeam.f`
(RMSF/dirty-beam construction via matched-filter correlation against the
same `cos_arr`/`sin_arr` templates `rm_synthesis` already builds),
`rm_restore.f` (Gaussian-restore step), `bw_depol_correct.f` +
`bw_depol_correct_setup.f` (analytic bandwidth-depolarization correction,
deferred — see T-future below), `peak_interp.f` + `quad_interp.f`/
`fourier_interp.f`/`fourier_interp_re.f`/`sinc_interp.f` (sub-pixel peak
interpolation), `index_absmax.f`.

Explicitly NOT ported: `extract_general_v3.f`/`extract_general_setup.f`/
`extract_single.f` (the old code's own from-scratch RM synthesis —
redundant, `rm_synthesis_mod.f90:675` already has a modernized
`extract_general_setup` with `cos_arr`/`sin_arr` built once and reused
across the whole run, `rm_synthesis.f90:2167-2173` — RM-CLEAN's own beam
construction should consume those same templates, not duplicate them);
`rm_clean_single.f`/`compute_dirty_rmbeam_single.f` (the newer "clean Q
and U separately" variant exploiting real-DFT Hermitian symmetry —
skipped per the user's own explicit call, since their README already
flags it as "needs more checks", not validated); the interactive
PGPLOT driver programs, `myplot1.f`, `modified_QU.f` + the random-number
generators (simulation/test harness, not runtime code); `fort_lib.f`
(grab-bag utility library — checked against `rm_synthesis_mod.f90` and
`compute_mean`/`dot_product_custom` already exist there, modernized,
matching `fort_lib.f`'s `mean`/`dotproduct` exactly — reuse those
directly rather than duplicating; only a new `compute_rms`, "rms about
the mean", needs adding in the same style, since nothing equivalent
exists yet and the core `rm_clean.f` path needs it for its stopping
criterion. `robust_meanrms.f` is unused scope — only ever referenced by
the shelved "single" variant driver and one dated backup file).

## Design decisions confirmed with the user

1. **Restore order (critical fix vs. the original code):** the original
   `rm_clean.f` sums `ClnFluxQ(i) = ClnFluxQ(i) + ResiQ(i)` (components
   + residual) *before* `rm_restore.f` ever runs, so the Gaussian
   restoring-beam convolution operates on the sum, not just the
   components. The ported version must convolve ONLY the pure
   delta-function clean-component map with the restoring beam, then add
   the (unconvolved) residual afterward — standard Högbom order.
2. **Complex clean only.** `rm_clean_single.f` (separate real-DFT-symmetry
   Q/U clean) is shelved, not ported, until validated independently of
   this work.
3. **RMSF is genuinely RM-dependent, not shift-invariant — confirmed by
   direct derivation, not by assumption.** Without BW-depol, per-channel
   weights are RM-independent, so `R(φ;φ₀) = (1/K)Σ wₖ·exp(−2i(φ−φ₀)λ²ₖ)`
   IS, as an algebraic identity, a pure function of the offset
   `Δ=φ−φ₀` alone — this holds for any λ² sampling, uniform or not (this
   was worked through explicitly in discussion and agreed). What the
   user actually needs recompute-at-exact-peak for is that an
   interpolated (generally off-grid) peak location can't be honored by
   index-shifting a discretely tabulated template without reintroducing
   interpolation error. With BW-depol switched on later, weights become
   genuinely `φ₀`-dependent per-channel (not just a scalar factor of
   `φ₀`), and shift-invariance breaks for real at that point — see
   T-future.
4. **Sub-pixel peak interpolation: interpolate Re and Im (or Q and U)
   channels separately, at the same offset, for every interpolation
   mode** (parabolic/sinc/Fourier-amplitude), not just for the
   Fourier-complex mode (which is the only mode that currently computes
   an interpolated phase at all — the other three modes leave phase at a
   stale/zero value and print `"Incorporate the scheme NOW!!!"` in the
   original code, a self-flagged incomplete state). Determine the offset
   from the magnitude/amplitude array as each mode already does, then
   apply that same interpolation kernel independently to Re and Im and
   evaluate both at that offset.
5. **`extract_general_setup` reuse — checked, compliant as-is.** Already
   in `rm_synthesis_mod.f90:675`, modernized, and its `cos_arr`/`sin_arr`
   output is built once per run and reused across every pixel
   (`rm_synthesis.f90:2167-2173`, `:4261` etc.) — exactly the structure
   RM-CLEAN's own beam-builder needs as an input. No changes needed to
   this subroutine.
6. **Parallelization: CPU/OMP first, one thread per pixel running its own
   full serial CLEAN loop** — pixels are independent, the loop within a
   pixel isn't (iteration i+1 needs iteration i's residual). Directly
   analogous to `gaussft_mod`'s own established "one thread, one full
   serial FFT-convolve per plane" pattern in `convolve_cubes.f90`,
   applied at pixel granularity. GPU is explicitly deferred: highly
   variable per-pixel iteration counts cause real warp/thread-block
   divergence, and this project's existing GPU kernels are straight-through,
   not work-stealing — a harder redesign to justify without real
   profiling showing CLEAN is actually the bottleneck.
7. **`compute_mean`/`dot_product_custom` reuse — checked, compliant.**
   `rm_synthesis_mod.f90:1614` and `:1627` already match `fort_lib.f`'s
   `mean`/`dotproduct` behavior, modernized. Reuse directly. One new
   subroutine needed: `compute_rms` (rms about the mean), same style,
   since nothing equivalent exists yet.
8. **Restoring-beam FWHM constant: keep `π` (not `2√3`)** — the user's
   deliberate, confirmed choice, not something to "fix" toward the
   Brentjens & de Bruyn (2005) eq. 61 constant.
   **UPDATE (found via the user cross-checking a real ASKAP low-band
   restoring-beam value against the RMSF's own expected ~50 rad/m²):**
   the implementation had silently drifted from this decision —
   `compute_rmsf_fwhm`/`compute_rmsf_fwhm_multiband` computed
   `0.5*pi/lsq_span`, not `pi/lsq_span`, an erroneous extra factor of
   0.5 with no basis in this decision or in the original `rm_restore.f`'s
   own definition of its `FWHM_RM` argument. Fixed to plain
   `pi/lsq_span`. Confirmed via full-codebase search that nothing
   compensates for the old factor elsewhere (`rm_synthesis.f90`/
   `rm_synthesis_mod.f90` have zero references to "fwhm" at all; CLEAN's
   own component-finding runs against the exact dirty RMSF, never this
   approximate value) — the fix's effect is confined to `restore_clean`'s
   restoring-beam width (now 2x wider, matching the true theoretical
   resolution), Gate 0's `drm_required` criterion (now correctly
   permissive), and the logged "Restoring beam FWHM" value. Follow-up:
   `tests/test_matched_filter_refine.f90`'s own coarsest resolution tier
   moved from `fwhm/2` to `fwhm/3` (2 real samples/fwhm turns out to
   exceed the fast path's reliable margin even for a noiseless single
   component — a real, previously-unexercised property of the tiered
   design, not a new bug; see that test's own `run_tier_mechanics`
   comment for the measured residual/threshold numbers).
9. **RM-dependent-RMSF optimization (pre-BW-depol): an offset-table keyed
   on `Δ=φ−φ_peak`, built once per dataset, interpolated into at each
   CLEAN iteration instead of recomputing the full `O(N_chan)` sum every
   time** — valid per point 3's derivation. Confined strictly inside
   `compute_dirty_rmbeam`'s own implementation, behind the same general
   `(φ_grid, RM_in, phase_in) → beam array` interface the original code
   already has, so nothing upstream (the CLEAN loop's own control flow)
   needs to know or care how the beam was actually computed — this
   specific acceleration strategy can be swapped out later (when
   BW-depol lands) without touching `rm_clean`'s own code. Must carry
   clear comments documenting this is a pre-BW-depol-only optimization
   and exactly why it stops being valid once `w_k` becomes `RM_in`-dependent
   (see T-future).
10. **Per-pixel valid-channel patterns are NOT guaranteed identical across
    a cube** — a global `badchan_file` is band-wide, but individual pixels
    can carry additional NaN-driven flagged channels (this is exactly
    what `MASK.CUBE.FITS` already records, per-pixel, per-channel). A
    single dataset-wide offset-table would be silently wrong for any
    pixel whose mask deviates from the dominant pattern. Resolved via a
    **pre-processing pass** (not on-the-fly caching during the CLEAN
    loop, which the user correctly flagged as relying on too many
    unstated assumptions about scheduling/locality): scan
    `MASK.CUBE.FITS` block-wise once, enumerate the distinct per-pixel
    valid-channel patterns present (keyed by a canonical hash), build one
    offset-table per distinct pattern; the actual CLEAN pass is then a
    pure lookup per pixel. A safety valve caps the number of distinct
    patterns tracked, falling back to always-recompute past that
    threshold rather than holding an unbounded number of tables — same
    "degrade gracefully, don't silently blow up" instinct as
    `max_common_bmaj` in `convolve_cubes`.
11. **Architecture: standalone tool first, consuming an existing dirty
    cube** (AMP/PHA), mirroring the `reproject_cubes`/`convolve_cubes`
    standing-alone-before-`match_cubes`-chained-them trajectory already
    proven in this project. Decouples CLEAN's very different
    parallelization shape from `rm_synthesis`'s own tile/OMP/GPU
    machinery while the port is still settling, and allows re-running
    CLEAN with different `niter`/`gain`/`thresh` without re-paying for
    the DFT step. An inline `rm_synthesis` option (the user's original
    framing) remains the eventual goal, deferred until the standalone
    version is proven — same trajectory `match_cubes` walked relative to
    the two tools it later chained.
12. **MASK.CUBE.FITS per-channel metadata is a genuine prerequisite, not
    something CLEAN can route around** — checked `rm_synthesis.f90:2793-2805`:
    `MASK.CUBE.FITS`'s own `ctype3`/`crval3`/`cdelt3` header is a single
    linear `FREQ` WCS, sized to the full (possibly multi-band-concatenated)
    channel axis (`naxes_mask(3)=nz_out`, `:2531`) with no `nbands`-specific
    handling found near that block — the same class of gap already
    documented elsewhere in this codebase for `use_auto_rm_range=1` being
    rejected outright for `nbands>1` (a single CRVAL/CDELT can't represent
    a concatenated multi-band channel list). Rather than rely on that
    linear axis at all (for either single- or multi-band), attach an
    explicit per-channel binary table to `MASK.CUBE.FITS` — same
    binary-table-extension pattern already proven for `CASAMBM`/`BEAMS`
    this project. **Promoted to T0 below, at the user's explicit
    request, so it doesn't keep delaying the CLEAN ticket itself.**

## Choosing parameters

`rmclean_mod` (T1) has no config/CLI layer yet — that's T2's job, not
started — so every parameter below is set explicitly by the caller today
(as in `tests/*.f90`). This section is the reference for what T2's own
defaults should be, and for anyone calling the module directly meanwhile.
Added at the user's own request for organized, scope-clean naming and an
explicit decision guide, after `lsq_ref` alone started doing double duty
for two genuinely unrelated concerns.

**Two deliberately separate reference parameters — do not conflate them:**

| Name | Role | Governs | Default recommendation |
|---|---|---|---|
| `lsq_ref_compute` | Internal computation reference (`build_rmsf_offset_table`, `compute_dirty_rmbeam_direct`, `get_drm`) | ONLY grid cost (`get_drm`'s bound) — never chi0 precision, once `comp_rm_refined`/`derotate_to_lsq_ref` are used correctly | `get_lsq_ref_compute(..., mode=lsq_ref_compute_mid, ...)` — cheapest safe grid, RECOMMENDED default |
| `lsq_ref_report` | Reporting convention (`derotate_to_lsq_ref`'s own target reference) | ONLY where chi0 is quoted — a bookkeeping choice, not a computation | `mode=lsq_ref_report_intrinsic` (this project's thesis convention) as the default; `mode=lsq_ref_report_centroid` as an explicit opt-in for statistically decorrelated (RM, chi0) uncertainty pairs (Brentjens & de Bruyn 2005-style) |

These answer different questions and are optimized by different,
non-interchangeable criteria — confirmed, not just asserted:
`mode=lsq_ref_compute_mid` minimizes `max_k|l_sq(k)-lsq_ref_compute|` (a
minimax/Chebyshev-centre problem over the two extreme channels only,
driven by numerical grid-stability, no noise involved); `mode=
lsq_ref_report_centroid` targets the (weighted) mean specifically to
decorrelate a jointly-fitted RM and chi0's own statistical uncertainty (a
linear-regression centering result, driven by noise/covariance, no
grid-cost involved). They coincide only by coincidence for special cases
(a single, symmetric, evenly-sampled band) and diverge meaningfully for
realistic, channel-imbalanced multi-band data (`tests/
test_optimal_lsq_ref.f90`'s own P+L case: `0.5806` vs
`0.3771`/`1.0011`/`0.0626` depending which of the four candidates you
compute). Neither is "more correct" — pick `lsq_ref_compute` for cost,
`lsq_ref_report` for what you want to quote, independently of each other.

**Interface, after three iterations driven directly by the user's own
pushback — each caught a real problem, not bikeshedding:**

1. First pass named the two "give me a good value" helpers
   `optimal_lsq_ref_midpoint` and `lsq_ref_report_centroid` — inconsistent
   with each other (one prefixed, one suffixed) and not self-explanatory
   about which VALUE each one filled in.
2. Second pass tried `suggest_lsq_ref_compute`/`suggest_lsq_ref_report`/
   `suggest_drm` — one consistent prefix, but "suggest" doesn't convey
   that the result is actually COMMITTED for the real computation, only
   that it's advisory — the user's own objection, with real established
   precedent behind the fix: HDF5's `H5Pget_chunk`/`H5Pset_chunk`, and
   cuFFT's own `cufftGetSize1d` (vs. `cufftEstimate1d`, its own
   heuristic-only sibling) both use "Get" specifically for the value
   that gets committed, not "Suggest"/"Estimate".
3. Final design, this section's own: `get_lsq_ref_compute`/
   `get_lsq_ref_report`, each taking a REQUIRED `mode` argument (an
   `integer, parameter` enum, matching this module's own pre-existing
   `fftw_forward`/`fftw_estimate` convention) plus an optional
   `fixed_value` (read only when `mode=..._fixed`) — resolving a THIRD
   objection along the way: passing a raw number through an optional
   argument on a "getter" is really a SETTER wearing a getter's name. The
   fix: a getter only ever computes/derives; a caller who wants a
   specific number they chose themselves just writes ordinary Fortran
   (`lsq_ref_compute = 0.734_sp`) — no module call needed at all, exactly
   how `niter`/`gain` already work. The getters exist only for values
   worth deriving: named, meaningful conventions (`..._intrinsic`,
   `..._centroid`) or simple derived statistics (`..._min`/`..._max`/
   `..._mid`) — plus `..._fixed` for symmetry/self-documentation
   (`call get_lsq_ref_report(l_sq, nchan, mode=lsq_ref_report_fixed,
   lsq_ref_report=x, fixed_value=0.734_sp)` reads more clearly at the
   call site than a bare assignment, even though it does the same thing).

Call sites use explicit `mode=`/keyword-argument syntax throughout (the
user's own explicit style request) — every `tests/*.f90` call reads
self-documenting, e.g. `call get_lsq_ref_compute(l_sq, nchan,
mode=lsq_ref_compute_mid, lsq_ref_compute=lsq_ref_compute)`.

Recognized modes (each set is its own independent `integer, parameter`
namespace — `lsq_ref_compute_mid=1` and `lsq_ref_report_intrinsic=1` do
NOT collide, they're never compared cross-set):

| `get_lsq_ref_compute`'s modes | `get_lsq_ref_report`'s modes |
|---|---|
| `lsq_ref_compute_mid` (=1, RECOMMENDED DEFAULT) | `lsq_ref_report_intrinsic` (=1, default) |
| `lsq_ref_compute_intrinsic` | `lsq_ref_report_centroid` |
| `lsq_ref_compute_centroid` | `lsq_ref_report_min` |
| `lsq_ref_compute_min` | `lsq_ref_report_max` |
| `lsq_ref_compute_max` | `lsq_ref_report_mid` |
| `lsq_ref_compute_fixed` (+ `fixed_value`) | `lsq_ref_report_fixed` (+ `fixed_value`) |

Extensible by design (the user's own framing): a new mode is just one
more named constant plus one more `case` branch, nothing structural.

**Reading a component's own precise RM/chi0 (not the restored map):**
always use `clean_complex`'s own `comp_rm_refined(j)` for the RM value,
never `rm_samp(j)` and never a fresh `peak_interp_parabolic` pass over
the restored map — root-caused empirically (see `derotate_to_lsq_ref`'s
own comment): the restored map's own peak is grid-quantised (up to
`dRM/2` away from the true continuous RM) even when it looks perfectly
clean and symmetric, and no amount of re-interpolating that output,
however exact, can recover information already discarded at the moment
the component was filed into an integer bin. The restored map itself
remains correct for its own separate uses — visual inspection, and
integrating an EXTENDED feature's total flux (`window_flux`) — just not
for a single component's own precise RM/phase.

**There is no raw `dRM` input anywhere in this module's API, by design**
— the user's own explicit simplification, after asking directly "why is
Nyquist even an option?": `lsq_ref_compute`/`lsq_ref_report` are genuine
free choices with no wrong answer, but sampling below the stability bound
is a CORRECTNESS failure, not a preference, so the only knob `get_drm`
exposes is `oversample`, a multiplier ALWAYS applied on top of the
mandatory floor — asking `get_drm` for an unsafe `dRM` is architecturally
impossible, not just discouraged.

`oversample`'s own floor is ENFORCED, not advisory, and is `2`, not the
bare 2-point Nyquist value of `1` — root-caused empirically (`tests/
test_drm_floor.f90`, `tests/test_drm_floor_enforcement.f90`), not assumed:
`oversample=1` does NOT diverge (residual power stays small in every
scenario checked) but recovers a WRONG chi0, not just an imprecise one —
a single point source (RM=50, chi0=0.3 rad, band-mean `lsq_ref_compute`)
came back at chi0=-0.025 rad, while `oversample=2` already recovered
0.3002 rad (matching `oversample=15`'s own 0.3000 to float32 rounding).
`get_drm` stops with a clear error if asked for `oversample<2` — checked
on one scenario only (single point source, no noise), so treat `2` as a
floor with real safety margin built in, not a knife-edge value tuned to
just this case. Default when `oversample` is omitted: `4` (the user's own
suggestion), giving further margin at a fraction of `oversample=15`'s own
grid cost (`nrm=80` vs `nrm=297` in the tested scenario).

**The two same-named "oversample" parameters are unrelated — do not
confuse them:**

| Where | Controls | Cost if increased |
|---|---|---|
| `get_drm`'s `oversample` | The OUTER RM-grid's own fineness relative to the bare stability floor (`nrm`, hence the size of every `nrm`-length array in `clean_complex`) | Expensive — scales the whole grid |
| `build_rmsf_offset_table`'s `oversample` | The INTERNAL 1D RMSF lookup table's own fineness relative to the outer grid's `dRM` | Cheap — a 1D table built once, independent of `nrm` |

Defaults used throughout this project's own tests: `4` for the former
(the enforced floor is `2`; see above), `20` for the latter (comfortably
cheap regardless of value).

**CLEAN loop parameters (`niter`, `gain`, `thresh`):** `gain=0.1` is the
conventional Hogbom choice, no reason found to deviate. `thresh` needs a
documented warning, not a blind default: a moderate `thresh` caused
`clean_complex` to stop after a single iteration on a broad
(Faraday-thick) feature in the P+L thesis scenario, because that
feature's own mean residual sits close to its peak from the very first
iteration — `thresh=0.0` (letting `niter` alone govern convergence) was
needed there. A single fixed default is genuinely risky across source
morphologies; T2 should either default to `thresh=0.0` with a generous
`niter`, or implement a shape-aware stopping heuristic, not just pick a
number.

## Ticket sequence

### T0 — `MASK.CUBE.FITS` per-channel frequency/band binary table

- **Status:** done. Pre-requisite for all of T1+; nothing else in
  this plan can rely on `MASK.CUBE.FITS` metadata until this lands.
- **Objective:** attach an explicit per-channel binary table extension to
  `MASK.CUBE.FITS` (own EXTNAME, e.g. `CHANFREQ`) recording each merged
  channel's true frequency (and, for multi-band runs, which band it came
  from), replacing reliance on the existing linear `ctype3`/`crval3`/
  `cdelt3` `FREQ` axis header — which is likely already wrong for
  multi-band runs (decision 12 above) and which a standalone RM-CLEAN
  tool would otherwise have no reliable way to trust.
- **Scope:** single mechanism for both single-band and multi-band runs
  (no special-casing) — always attach the table, regardless of `nbands`.
  Columns: channel index (0-indexed, matching the `BEAMS`-table `CHAN`
  convention already established), frequency (Hz, matching the existing
  `MASK.CUBE.FITS` `cunit3` convention), band index (only meaningful for
  multi-band; a single-band run just gets band index 0 for every row).
  Does not change `MASK.CUBE.FITS`'s existing linear `FREQ` axis header
  (left as-is for now, for any existing downstream reader relying on it)
  — this ticket is purely additive.
- **Correctness gate:** verify against real single-band and multi-band
  runs — table row count matches `nz_out`/`naxes_mask(3)` exactly; single-band
  table's frequencies match what the existing linear CRVAL/CDELT axis
  would predict (a direct check that the new table isn't just
  *differently* wrong); multi-band table's frequencies match each band's
  own true per-channel frequency (verifiable against the input Q/U cubes'
  own headers directly, not derived); full regression
  (`scratch/make_all.sh` + `tests/run_tests.sh`) stays clean throughout.
- **Implementation:** `rm_synthesis.f90`, new declarations at the
  CASAMBM/BEAMS block's own location (`chanfreq_chan`/`chanfreq_band`/
  `chanfreq_freq`, `chanfreq_colnum`/`chanfreq_fitsstat`/`chanfreq_i`/
  `chanfreq_iband`), and a new `if(out_mask_open)` block placed at the
  very end of the run, after the tile-processing loop and before the
  final `FTCLOS`/deallocate sequence (`L_sq`/`flag_arr_out`/
  `band_offset`/`band_nz` are all still populated and unchanged there).
  `FTIBIN`/`FTGCNO`/`FTPCLJ`/`FTPCLD` build a 3-column (`CHAN` 1J,
  `FREQ` 1D/Hz, `BAND` 1J) `CHANFREQ` binary table extension on unit 43
  (MASK), one row per merged channel, `FREQ` derived from `L_sq`
  (`freq_Hz = c_velocity(Mm/s)*1e6/sqrt(L_sq)`, matching MASK's own
  existing `cunit3='Hz'` convention) rather than duplicating any
  frequency computation. `BAND` filled via `band_offset`/`band_nz`
  range-walking (0 for every channel when `n_bands_t2<=1`). Purely
  additive — MASK's existing linear `ctype3`/`crval3`/`cdelt3` header is
  left untouched for any existing downstream reader.
- **Correctness gate (all passed):**
  - Multi-band run (`TEST.Q/U.FITSCUBE` as the reference band + `TEST_BAND2.Q/U.FITSCUBE`,
    the same fixture pair `tests/run_tests.sh`'s own T1/T2 multi-band test uses):
    `CHANFREQ` row count (350) matches `MASK.CUBE.FITS`'s own `naxis3`
    exactly; `BAND` correctly splits 200 rows to band 0 / 150 rows to
    band 1 (matching each input cube's own `naxis3`); `FREQ` for each
    band's own rows matches that band's own INPUT cube header
    (`crval3`/`cdelt3`) directly — not re-derived from anything already
    written by this program — to `rtol=1e-6` (residual is single-precision
    `L_sq` round-trip, expected and consistent with the rest of this
    codebase's own `real(sp)` convention).
  - Single-band run (`TEST.Q/U.FITSCUBE` alone): `CHANFREQ` row count
    (200) matches `naxis3`; every row's `BAND=0`; `FREQ` matches what
    `MASK.CUBE.FITS`'s own existing linear `crval3`/`cdelt3` axis would
    predict, to `rtol=1e-6` — confirming the new table isn't just
    *differently* wrong, it agrees with the one piece of existing
    metadata that was already trustworthy for the single-band case.
  - `write_mask_output=n`: run completes cleanly, no `MASK.CUBE.FITS`
    (and therefore no `CHANFREQ` table) written, no crash — confirms the
    new block's `out_mask_open` gate holds.
  - All 4 build flavours (`scratch/make_all.sh`) clean; full
    `tests/run_tests.sh` 49/49 pass, unaffected (purely additive change,
    confirmed by the existing bit-identical `compare_cubes.py` checks
    still passing throughout).

### T1 — RM-CLEAN core algorithm module (`rmclean_mod`)

- **Status:** done. Pure computation only, no FITS I/O, no standalone
  binary yet — mirrors `gaussft_mod`'s own split between narrowly-scoped
  computation and a caller that owns I/O/config. New file: `src/rmclean.f90`.
- **Objective:** port the actual CLEAN algorithm
  (`rm_clean`/`compute_dirty_rmbeam`/`rm_restore` from
  `~/softwares/CURR_DEVEL/RM_CLEAN_TESTS/`), modernized and with every
  fix from the design decisions above, into a tested, self-contained
  module — the prerequisite for T2's standalone tool.
- **Scope:** `index_absmax`, `peak_interp_parabolic` (decision 4),
  `plan_fourier_interp`/`fourier_interp_complex` (FFTW-based, replacing
  the original's in-house `fft1d`/`ifft1d`, per explicit instruction),
  `rmsf_table_t`/`build_rmsf_offset_table`/`compute_dirty_rmbeam_direct`/
  `compute_dirty_rmbeam` (decisions 3/9), `clean_complex` (decision 1),
  `compute_rmsf_fwhm`/`restore_clean` (decisions 1/8). `rms_about_mean`:
  a private local copy inside `rmclean_mod` itself (the public
  `compute_rms` added to `rm_synthesis_mod.f90` per decision 7 is for
  T-future's eventual inline-into-`rm_synthesis` path — `rmclean_mod`
  deliberately doesn't `use rm_synthesis_mod` at all, same standalone-module
  precedent as `gaussft_mod`/`commonbeam_mod`, so it needs its own tiny
  copy today). Sub-pixel interpolation scope deliberately reduced from
  the original's 4 modes (parabolic/Fourier-amplitude/Fourier-complex/
  sinc) to parabolic only, since decision 4's fix (interpolate Re/Im
  independently at the magnitude-derived offset) is what actually
  matters and preserving all 4 modes added little beyond what the
  original comment already called out as "(in use)"; `fourier_interp_complex`
  is still built and tested as a general-purpose utility, just not wired
  into `clean_complex`'s own peak-search as an alternative mode.
- **Three real bugs found during the port, all via direct re-derivation
  against the algorithm's own math, not just carried over or guessed at:**
  1. `quad_interp.f`'s own closed-form vertex-value formula
     (`peak_val = beta - 0.25*(alpha-gama)*peak`) used the *absolute*
     fractional bin index where the correct formula requires the *local*
     vertex offset alone — only numerically close in the original because
     it was always called on a small (5-point) window. My own first
     draft of the Re/Im-separate version repeated the same class of
     mistake before being caught by re-deriving the general 3-point
     quadratic (`y(x)=A x^2+Bx+C` via standard central-difference
     coefficients) and evaluating it directly at the externally-supplied
     offset, rather than reusing the vertex-only shortcut.
  2. The original's `ClnFluxQ`/`ClnFluxU` component bookkeeping
     accumulated `frac*ResiQ(imax)`/`frac*ResiU(imax)` — dimensionally
     inconsistent with `frac*re_beam(i)`/`frac*im_beam(i)`, what is
     actually subtracted from the residual every iteration, and unrelated
     to it. Fixed to the standard Hogbom convention: the delta-function
     component's own complex value is `frac*exp(i*phase_val)` (magnitude
     `frac`, at the peak's own phase) — `frac*cos(phase_val)`/
     `frac*sin(phase_val)` here.
  3. The original passed `phase_val` (the residual's own `atan2(Im,Re)`
     angle, already in the doubled `2*chi0` convention) straight into
     `compute_dirty_rmbeam`'s `phase_in`, which itself doubles internally
     (`2*(RM_in*L_sq+phase_in)`, matching the standard
     `P(L_sq)=p*exp(2i*chi0)*exp(2i*RM*L_sq)` Faraday relation) — doubling
     it again. Verified directly against this module's own
     `compute_dirty_rmbeam` (not just derived on paper) that `phase_in=
     phase_val/2` is the self-consistent choice: the beam's own value at
     its exact peak (`delta=0`) must equal `exp(i*phase_val)`, unit
     magnitude at exactly the residual's own phase.
- **FFTW normalization bug (own new code, caught before it shipped):**
  `fourier_interp_complex`'s first draft carried over the original's own
  `norm=nout/npts` scale factor unchanged — correct for the original's
  own `ifft1d`, which (like most hand-rolled FFT/IFFT pairs) normalizes
  internally by `1/N`; wrong for FFTW's `dfftw_execute_dft`, which is
  unnormalized in both directions. Fixed to `norm=1/npts`, the only
  factor needed once FFTW itself does no internal normalization at all.
- **`rm_restore.f`'s own FFT-shift sequence not ported as-is:** the
  original does `fft` -> `fftshift` on both signal and kernel -> multiply
  -> `ifft` -> `fftshift` on the result, an ordering not straightforward
  to verify correct by inspection. Reimplemented as the standard,
  textbook circular-convolution-via-FFT recipe instead: build the
  Gaussian kernel directly in *wrapped* form (index 1 = RM=0, indices
  past the midpoint = the "negative lag" tail wrapped to the array's own
  end — the standard layout circular convolution via FFT requires), FFT
  both arrays with no shifting at all, multiply, inverse FFT, divide by
  `nrm`. Same underlying convolution, easier to verify by construction.
- **Correctness gate (all verified numerically, not just derived):**
  1. `peak_interp_parabolic`: synthetic Gaussian-envelope complex signal
     with a known sub-pixel peak (offset 0.35 bins) — recovered offset
     to `7e-7`; symmetric on-grid peak recovers offset `~0` and the exact
     injected amplitude; degenerate (flat) input returns offset `0`, no
     crash.
  2. `fourier_interp_complex`: pure sinusoid at an exact DFT-basis
     frequency reproduced to `~6e-7` at every interpolated point (the
     defining property of bandlimited interpolation), and exactly (`0`
     error) at the original sample locations.
  3. `compute_dirty_rmbeam` (table-based) vs. `compute_dirty_rmbeam_direct`
     (exact `O(nchan)` reference): cross-validated at 5 off-grid
     `(RM_in,phase_in)` points, agreement to `~1.6e-6`. On-grid
     self-beam sanity check (`RM_in` exactly on the `rm_samp` grid,
     `phase_in=0`) gives `re=1.0,im=0.0` exactly, as physically required
     (a point source's own beam at its own peak).
  4. `restore_clean`: delta-function-input identity check (convolving a
     unit spike with the restoring beam must reproduce the beam itself,
     centred at the spike) — matched the analytic Gaussian to `~1e-6`
     directly, with the one larger residual (`2.8e-4` at the array's far
     edge) confirmed via independent calculation to be the exact,
     expected periodic-wraparound contribution of the circular
     convolution (`exp(-37^2/(2*sigma^2))` matched the measured error to
     6 significant figures) — not a bug. Zero-component/nonzero-residual
     passthrough: output exactly equals the residual, `0` error.
  5. **End-to-end**: synthetic single point source (`RM=23.7` rad/m^2,
     off-grid; intrinsic angle `chi0=0.35` rad; amplitude `5.0`; 200-channel
     700-900 MHz synthetic band; no noise) run through
     `clean_complex`+`restore_clean` in full: recovered dominant-bin
     amplitude `4.9995` (true `5.0`); recovered phase `0.7000001` (true
     `2*chi0=0.7`, `~6e-6` error — direct empirical confirmation that bug
     3's phase-convention fix is correct, not just algebraically argued);
     recovered RM snapped to the nearest grid bin (`24.0` vs true `23.7`,
     within `dRM/2` as expected for delta-function component placement);
     residual power dropped to `2.7e-13` of the dirty map's own power
     (full convergence, as expected for a noise-free single point source
     matching the basis exactly).
- **Build/test:** `src/rmclean.f90` compiles clean standalone (`gfortran`
  + `-fopenmp -lfftw3`); all 4 `rm_synthesis` build flavours
  (`scratch/make_all.sh`) clean; full `tests/run_tests.sh` (49 pre-existing
  + 6 new, see addendum below) pass, unaffected (this ticket adds one new
  self-contained module file plus the additive `compute_rms` in
  `rm_synthesis_mod.f90`, neither consumed by any existing tool yet).

**Addendum — thesis-scenario test, two more bugs found:** at the user's
request, added `tests/thesis_scenario_rmclean.f90`, a single line-of-sight
reproduction of Raja (2014) Chapter 6 Figures 6.1/6.2/6.3 (Table 6.1/6.2
exact: a point source at RM=-100 amplitude 15, a Faraday-thick top-hat
RM 100-130 total amplitude 5, both PA=0; CLEAN parameters niter=1000,
loop-gain=0.1 taken directly from the figures' own panel labels), cleaned
from P-band alone, L-band alone, and P+L combined — no FITS/cube I/O
needed, `rmclean_mod` used directly. Wired into `tests/run_tests.sh` as
section 23 (compiles fresh each run, skips gracefully if FFTW3 isn't
available).

This test caught two further real bugs beyond the three already listed
above, both via comparing against Table 6.1's own published numbers, not
by inspection:

1. `compute_rmsf_fwhm` (T1's own new subroutine, not carried from the
   original) came out **negative** for every case except one. Root
   cause: the underlying `lsq_span` formula's sign depends on which end
   of the array `l_sq(1)`/`l_sq(nchan)` represents, and the original
   `rm_restore.f` this was ported from implicitly assumed the *opposite*
   ordering convention from the one `L_sq` actually has everywhere else
   in this project (`rm_synthesis.f90`'s own comment: descending, i.e.
   `l_sq(1)` is the largest/lowest-frequency channel). A span is a
   magnitude, not an ordering-dependent signed quantity — fixed with
   `abs()`, refactored into a shared `padded_lsq_span` function reused by
   both `compute_rmsf_fwhm` and the new multi-band variant below.
2. Feeding a naively-concatenated P+L channel array into
   `compute_rmsf_fwhm` gave a wildly wrong (order-of-magnitude-too-small)
   FWHM. That subroutine's own span is effectively
   `max(l_sq)-min(l_sq)` across the *whole* array — for two
   widely-separated, non-contiguous bands this is dominated by the large
   *gap* between them, not by either band's own actual channel coverage.
   The thesis states the correct rule directly (Sec 6.1.4): "the
   effective RM-resolution is provided by the **sum** of the
   lambda-squared-spans at the two individual bands" — added a new
   `compute_rmsf_fwhm_multiband(l_sq, nchan, band_offset, band_nz,
   n_bands, fwhm_rm)`, same `band_offset`/`band_nz` per-band segmentation
   convention `rm_synthesis.f90` already uses, summing each band's own
   `padded_lsq_span`. Verified directly against Table 6.1 itself: P alone
   0.2007 vs tabulated 0.201; L alone 0.01255 vs tabulated 0.013; summed
   0.2133 vs the table's own P+L value 0.214 — all three matching to
   within the table's rounding.

**Addendum revision — flux metric, phase-reference convention, and a
genuine CLEAN-divergence bug, all found from further user review:**

1. **Flux metric was wrong.** The first version reported a bare restored
   *peak* inside the `[100,130]` window for the Faraday-thick component.
   Per the user's own correction: a point source's peak already equals
   its total flux (unresolved, standard convention), but an *extended*
   feature's total flux is the *integral* (area) of its restored profile
   over its own extent — exactly how a real interferometric image's
   resolved-source flux is measured, not read off a single peak pixel.
   Fixed with a new `window_flux` (Riemann sum, `Σ spec(j)*dRM` inside
   `[lo,hi]`), then further corrected to divide by the restoring beam's
   own continuous area (`sigma*sqrt(2*pi)`, `sigma=0.42466*FWHM_RM`) —
   the restoring kernel is normalized to unit *peak* height, not unit
   area (`restore_clean`'s own `kernel_c(i)=exp(...)`, height 1 at
   center), so integrating the restored map directly double-counts the
   beam's own area on top of the true flux (caught empirically: an
   early attempt reported L-alone recovering `1526%` of the true flux,
   obviously unphysical, before the beam-area correction).
2. **Phase-reference convention: `lsq_ref` must be an explicit,
   caller-supplied parameter, not silently computed as `mean(l_sq)`
   inside `compute_dirty_rmbeam_direct`/`build_rmsf_offset_table`.** Per
   the user's explicit instruction ("My thesis refers to lambda_sq=0
   for phase reference, and I want to use the same here"), the test now
   injects the sky model and builds the dirty spectrum both at
   `lsq_ref=0` (no subtraction), matching the thesis's own convention,
   rather than this project's usual mean-referencing (used purely for
   numerically smaller phase arguments — a cosmetic convenience, not a
   physical requirement). This is a real correctness constraint, not a
   style choice: `R(Δ)` is only a pure function of the offset
   `Δ=RM-RM_true` for a FIXED reference shared between the injected data
   and the matched-filter template — mismatching references was
   confirmed empirically to distort *relative* amplitudes between
   components, not just add a rotation. `compute_rmsf_fwhm_multiband`
   also needed a genuine fix here, not just this test: the naive
   concatenated-array approach used `max(l_sq)-min(l_sq)` across bands,
   which is dominated by the inter-band gap; fixed by summing each
   band's own `padded_lsq_span` directly (already covered above).
3. **`lsq_ref=0` triggered a genuine CLEAN divergence, root-caused as a
   Nyquist-type sampling violation, not a numerical-precision issue —
   and the fix is now a DERIVED bound, not a tuned constant.** With
   `lsq_ref=0`, the dirty spectrum's Re/Im oscillate at a rate set by
   each channel's *absolute* λ² (~1.1 for P-band), not the small spread
   around a mean (~0.2) that mean-referencing produces. The original
   RM-grid spacing (`dRM=2.0`, `nrm=226`) was too coarse for
   `peak_interp_parabolic`'s 3-point quadratic fit to stay valid across
   that oscillation, and CLEAN's residual grew exponentially
   (`111 -> 8e16` over ~700 iterations) — confirmed via a dedicated A/B
   test (`lsq_ref=0` vs `lsq_ref=mean`, otherwise identical) that isolated
   `lsq_ref=0` as the true trigger, and via an independent finer-grid
   test (`dRM=0.1`) that converged cleanly under otherwise-identical
   conditions. The first fix hardcoded that empirically-found `dRM=0.1`
   directly — flagged by the user as looking ad hoc, correctly. Replaced
   with a genuine sampling-theorem derivation, `rmclean_mod`'s new public
   `suggest_drm(l_sq, nchan, lsq_ref, oversample, drm_max)`:
   the dirty spectrum is `P(φ)=(1/K)Σ_k p_k·exp(-2iφ(l_sq(k)-lsq_ref))`,
   a sum of complex sinusoids in φ. Its *magnitude* envelope `|P(φ)|` is
   exactly independent of `lsq_ref` (a reference change multiplies the
   whole sum by `exp(-2iφ·Δref)`, unit modulus for every φ) and is set
   only by how the channels' own λ² are spread relative to each other —
   the SPAN, `compute_rmsf_fwhm`'s own quantity, the physically
   meaningful resolution (the user's own analogy: a wave packet's GROUP
   velocity/envelope, which carries the actual information content and
   cannot be changed by a mere choice of coordinate origin). But CLEAN
   operates on Re(φ)/Im(φ) *separately* (needs the actual complex value,
   not just the magnitude, to know what to subtract) and each term k is
   its own sinusoid with period `π/|l_sq(k)-lsq_ref|` — the analogue of a
   carrier's PHASE velocity, conveying no extra resolution/information
   by itself, but which the grid must still track faithfully since raw
   Re/Im values are what get sampled, fit, and subtracted. The binding
   constraint is therefore the *farthest* channel from `lsq_ref`:
   `max_k|l_sq(k)-lsq_ref|`. Bare two-point Nyquist for that term gives
   `dRM ≤ π/(2·max_offset)`; an `oversample` factor (used at `15` in the
   test) tightens this further, since `peak_interp_parabolic`'s *local*
   3-point quadratic fit needs several samples per cycle to be valid,
   not merely the two points needed to avoid aliasing a global
   reconstruction (a distinct, weaker requirement). Re-deriving `dRM`
   this way (rather than the hardcoded literal) lands at `dRM≈0.0947`
   (`nrm=4748`) for this scenario — matching the empirically-found
   `0.1` closely, confirming the original number was in fact the right
   order of magnitude, just not derived. All checks still pass
   identically with the derived value. Separately, `peak_interp_parabolic`
   gained a scientifically-grounded zero-magnitude guard (`log10(0)` is a
   genuine domain violation; falls back to the sampled center value,
   mirroring the pre-existing flat-curvature fallback) — a real
   robustness fix, verified not to be what resolved the divergence on
   its own (the dRM fix is the actual root cause).
4. **Pass/fail criteria rebalanced per the user's explicit direction:**
   strict, quantitative checks belong only on the *isolated* Faraday-simple
   point source ("retrieve at its RM location with nearly equal power as
   simulated"); the Faraday-thick component (never fully resolved by
   design, per the thesis's own scenario) uses generous, qualitative
   markers instead ("significant fraction... restored"), not a number
   tuned tight to whatever the code happens to currently produce. Point
   source: within `1.0 rad/m^2` and `15%` of true amplitude (actual
   deviations ~0.1 rad/m^2 / ~2.5%, so the tolerance is not
   loosened-until-passing). Thick component: P+L flux `>5x` P-alone flux
   (actual ~27x) and `>10%` of the true simulated total flux (actual
   ~40%) — both generous relative to the observed margins.
5. **Plotting (`tests/plot_thesis_scenario_rmclean.py`), matching the
   thesis's own Figure 6.1-6.3 panel style directly** (dirty Re/Im/Amp +
   cleaned amplitude, dirty + restored phase, true input model, λ²
   coverage), two more bugs found by inspecting the rendered plots
   against the thesis's own figures:
   - The λ² coverage panel initially tried to detect per-band segments
     from λ²-spacing statistics alone (gap > 10x the global median
     diff) — this silently shattered the P-band into 61 spurious
     single-point "segments", because P-band's own channel-to-channel λ²
     spacing (coarser, near the low-frequency end) is itself ~30x larger
     than L-band's, so no single global threshold can separate a real
     inter-band gap from a coarser band's own ordinary spacing. Fixed by
     having the Fortran side emit explicit band membership (a `band`
     column in `<slug>_lsq.csv`, from the same `band_offset`/`band_nz`
     already used for `compute_rmsf_fwhm_multiband`) rather than
     inferring it — the plot now groups by true band membership,
     confirmed visually to render both the L-band and P-band segments
     correctly in the combined case.
   - Once `lsq_ref=0` was correctly implemented with the tightened,
     `suggest_drm`-derived grid, the dirty Re/Im panels show
     many oscillation
     cycles across the RM range, matching the thesis's own Figure 6.1/6.2
     (previously showed only one cycle, under the old mean-referenced,
     coarser-grid version) — confirmed by visual inspection of the
     regenerated PNGs.
6. **Channel counts (`nchan=61` P-band, `nchan=121` L-band):** not
   independently specified anywhere in Table 6.1 itself; inherited from
   `tests/make_thesis_scenario_cubes.py`'s own pre-existing precedent
   (`BAND_P`/`BAND_L` dicts) for this same thesis scenario, used
   as-is rather than invented fresh for this test.

Final result (all checks passing, current numbers): P-alone point
source `97.6%` of true amplitude / thick-component flux `1.5%` of true
total; L-alone point `86.9%` / thick `7.8%`; P+L combined point `97.8%`
/ thick `39.9%` — matching the thesis's own qualitative conclusion (Sec
6.1.2-6.1.4) that combining bands recovers a significant fraction of the
extended feature's flux that neither band alone can, while the isolated
point source is recovered reliably and to essentially full power in
every case.

**Addendum — build coverage, a real precision gap, and lsq_ref flexibility
(post-commit follow-up, user review):**

1. **`scratch/make_all.sh` now also builds the ancillary standalone tools**
   (`reproject_cubes`, `convolve_cubes`, `match_cubes`) — previously it
   only cycled `rm_synthesis`'s own 4 OMP/GPU flavours, silently never
   rebuilding the other three tools at all. Checked the Makefile
   directly rather than assuming: those three targets hard-code
   `CPU_OPTFLAGS`/`CPU_OMPFLAGS` and never read the `OMP=`/`GPU=`
   variables, so unlike `rm_synthesis` there is no 4-way matrix to
   rebuild — one clean+build per tool suffices. (RM-CLEAN itself still
   has no Makefile target — T2, not yet started.)
2. **Precision audit (`float32` vs `float64`), requested directly —
   found one real gap.** Read the whole module against this question:
   every numerically-delicate part (the RMSF offset table, `compute_
   dirty_rmbeam`'s interpolation, `restore_clean`'s convolution, FFTW's
   own `complex(dp)` arrays matching the double-precision `-lfftw3` it
   links against) already does its internal arithmetic in double
   precision, casting to `real(sp)` only at the public boundary — except
   `clean_complex`'s own accumulator loop, which was `real(sp)`
   throughout. Unlike the PEAK-FINDING step (recomputed fresh from the
   current residual every iteration, so any single-precision rounding
   there is bounded and non-compounding), the residual/component
   SUBTRACTION (`resid = resid - frac*beam`) is applied repeatedly to
   the SAME accumulator across up to `niter` (often 1000+) iterations —
   exactly the case where rounding error compounds. Confirmed this was
   real, not theoretical: a noise-free single-point-source case had
   previously converged to a residual power of only ~`1e-13` of the
   dirty map's own power — suspiciously close to `float32`'s own
   `~1.2e-7` relative-amplitude (squared: `~1.4e-14` power) floor,
   meaning the accumulator's own precision, not the algorithm or the
   data, was capping how deep CLEAN could converge. Fixed by promoting
   `resid_re`/`resid_im`/`comp_re`/`comp_im` to `real(dp)` internally
   for the duration of the loop (a fresh `real(sp)` snapshot is taken
   each iteration for the unchanged peak-finding calls), casting back to
   `real(sp)` only in the final output — the same "double internally,
   single at the boundary" convention every other routine in this module
   already follows, extended to the one place that didn't. Negligible
   cost (a few extra KB per pixel — CLEAN runs one line of sight at a
   time, not per-cube). Real data's own thermal noise floor is almost
   always far above `float32`'s epsilon regardless, so this rarely
   mattered in practice — but it was a real, identifiable, and now
   removed limitation, not a hypothetical one.
3. **`lsq_ref` flexibility + `derotate_to_lsq_zero`, at the user's
   explicit request** ("users should be given a chance to select
   reference lambda-sq... we would still like to report extracted
   intrinsic polarisation angle de-rotated to lambda_sq=0"). `lsq_ref`
   was already a free, explicit parameter to `compute_dirty_rmbeam_
   direct`/`build_rmsf_offset_table` (this session's earlier fix) — the
   missing piece was a way to report results at the STANDARD
   `lambda_sq=0` convention regardless of which `lsq_ref` was used
   internally for the actual computation. New public subroutine,
   `derotate_to_lsq_zero(rm_samp, nrm, lsq_ref, re_in, im_in, re_out,
   im_out)`: the dirty/restored spectrum built at reference `lsq_ref` is
   `P_ref(phi) = P_0(phi)*exp(2i*phi*lsq_ref)` (the same general
   reference-change relation `suggest_drm`'s own comment
   derives), so `P_0(phi) = P_ref(phi)*exp(-2i*phi*lsq_ref)` recovers
   the `lambda_sq=0` convention exactly, pointwise, with no re-synthesis
   needed — a genuine closed-form identity, not an approximation.
   `|P_0|=|P_ref|` exactly (a unit-modulus factor changes only phase),
   directly verified in the new test below.
   - **This connects directly to `suggest_drm`'s own finding**
     (a band's own centroid as `lsq_ref` needs a far coarser, cheaper RM
     grid than `lsq_ref=0` does) — this is very likely the actual
     motivation for wanting `lsq_ref` flexibility at all, so the two
     features were designed to compose: pick whatever `lsq_ref` is
     computationally convenient, then de-rotate to report chi0 at the
     standard convention.
   - **A second, genuine precision coupling found while validating this
     (not by inspection — the first version of the new test failed,
     revealing it):** chi0 = `0.5*phase_val - RM_found*lsq_ref`. At
     `lsq_ref=0` this is exact regardless of any imprecision in
     `RM_found` (the multiplication by zero cancels it) — the reason
     `lsq_ref=0` is uniquely "free" for phase reporting, not a
     coincidence. At a NONZERO `lsq_ref`, the SAME `RM_found`
     imprecision (bounded by grid resolution, roughly half a cell even
     after sub-pixel refinement) gets multiplied by `lsq_ref` before
     reaching chi0 — a real, physical trade-off for choosing a coarser,
     cheaper grid, not a bug. Documented directly in `derotate_to_lsq_
     zero`'s own comment, with the practical mitigation: sub-pixel
     refine `RM_found` first (`peak_interp_parabolic`, the same routine
     `clean_complex` already uses internally every iteration) before
     reading off chi0, rather than settling for the nearest grid point.
4. **New test, `tests/test_rmclean_lsqref_flex.f90`** (wired into
   `tests/run_tests.sh` as section 24): a single point source with a
   DELIBERATELY NONZERO intrinsic angle (chi0=0.3 rad — the existing
   thesis-scenario test's own sky model is all PA=0, which cannot catch
   a sign error in a derotation formula), cleaned once at `lsq_ref=0`
   and once at the band's own mean, confirming: (a) the mean reference
   allows a `>5x` coarser grid (confirmed `10.7x` here); (b) both cases
   recover chi0 within a tolerance DERIVED from each case's own grid
   resolution (`0.02` rad at `lsq_ref=0`, exact; `0.5*dRM*|lsq_ref|` at
   the coarser reference, not an arbitrary loosened number — matches
   finding 3 above directly); (c) `|P|` is unchanged by derotation to
   better than `1e-4` relative, confirming the transform is phase-only.
   All checks pass; full regression (`scratch/make_all.sh` + `tests/
   run_tests.sh`) is 65/65, up from 57/57 (this test's own 7 new checks
   plus the ancillary-tools build coverage above).

**Addendum — the actual root cause of the chi0-precision gap, and a real
fix (`comp_rm_refined`), from a follow-up discussion pressing on whether
a coarse grid can be trusted at all:** the user pushed back hard, and
correctly, on the framing above — derotation is an exact, deterministic
phase rotation and cannot itself change any error, so the `~0.22` rad
discrepancy had to come from somewhere else being conflated with
`lsq_ref` choice. Chased in stages, each checked empirically rather than
argued from theory alone:

1. **Confound #1 (found, fixed): the two cases used different grid
   spacings.** `suggest_drm`'s bound (stability only) differs
   ~10x between the two references; forcing the SAME (fine) grid onto
   both cases dropped the discrepancy from `0.22` to `0.02` rad —
   confirming most of it was about grid resolution, not `lsq_ref` per se.
2. **Confound #2 (found, understood, not the actual fix): Fourier
   interpolation, tried as "the" fix, initially made things much
   worse** (`0.71` rad). Root-caused directly: `fourier_interp_complex`
   does a *global*, whole-array FFT upsample, which implicitly assumes
   the array is periodic; the worst reconstruction error landed exactly
   on the array's own last sample, the classic Gibbs/wraparound
   signature. Confirmed this is real by comparing against an
   independently-computed exact reference: reconstruction near the true
   peak was accurate to `~3e-4` (excellent), but the global array's own
   distant edges corrupted it. Even confining the Fourier upsample to a
   local, genuinely quiet window (`+/-8*FWHM_RM`, edge amplitudes
   `~1e-6`) still gave a bad answer (`0.71` rad) — ruling out "just a
   local window" as the fix, and pointing at something else entirely.
3. **The actual root cause: not an interpolation problem, a discarded-
   information problem.** The coarse grid's own restored profile, right
   at its recovered peak, was perfectly clean and symmetric
   (`9.547, 10.00, 9.547`) — but centred exactly ON the grid point
   `rm(imax)`, itself `~0.22` rad/m^2 from the true continuous RM. No
   interpolation of an already-quantised output, however exact, can
   recover information that isn't there. `clean_complex` already computes
   a precise sub-pixel `peak_loc` every single iteration (via
   `peak_interp_parabolic`, correctly used to build the right beam for
   subtraction) — but discarded it the moment the component was filed
   into `comp_re(imax)`/`comp_im(imax)`, an integer grid bin. Verified by
   patching a scratch copy of `clean_complex` to accumulate a
   flux-weighted average of `peak_loc` for the dominant bin: this
   recovered `RM=50.0000`, `chi0=0.3000` — EXACT to the precision
   checked, at the identical cheap, coarse grid that gave `0.22` rad
   error before.
4. **Fix, made permanent:** `clean_complex` gained a new required output,
   `comp_rm_refined(nrm)` — per-bin (not just the single dominant
   component, so this stays correct for multi-component scenarios like
   the thesis scenario's own point-source-plus-thick-feature case),
   flux-weighted `peak_loc`, falling back to `rm_samp(j)` for bins that
   never received any component flux (moot there, since `comp_re`/
   `comp_im` are exactly zero too, but keeps the array always
   well-defined). `derotate_to_lsq_zero`'s own doc now recommends this
   explicitly: read chi0 off `comp_re`/`comp_im` (the pure component map)
   with `comp_rm_refined(j)` as the RM value, NOT the restored map's own
   peak nor a fresh `peak_interp_parabolic` pass over it. The restored
   map remains correct for its own existing uses (visual inspection,
   integrating an EXTENDED feature's total flux) — this only changes how
   a single component's own precise RM/phase should be read.
   `tests/test_rmclean_lsqref_flex.f90` rewritten to use this directly;
   both `lsq_ref=0` and `lsq_ref`=band-mean cases now recover chi0 to the
   SAME tight tolerance (`0.02` rad) at their own, very different
   (`10.7x`) grid costs. `tests/thesis_scenario_rmclean.f90` also updated
   to report its own point-source RM via `comp_rm_refined` (a real,
   if here modest, precision improvement, unrelated to that test's own
   already-passing checks). Full regression stays 65/65.
5. **Direct answer to "can the RM-grid be coarse (nonzero, band-centroid
   `lsq_ref`) and still trust chi0(lambda_sq=0)?": yes — the coarse grid
   was never the actual obstacle.** The real requirement is using
   `comp_rm_refined`, not a bare grid coordinate, when converting a
   recovered component's phase to the standard convention. `lsq_ref`
   choice is now a genuinely free, purely computational decision
   (grid cost via `suggest_drm`) decoupled from achievable chi0
   precision.

**UPDATE (2026-07-31) — the "once `comp_rm_refined`/`derotate_to_lsq_ref`
are used correctly" precondition above was NOT actually satisfied in
`rmclean_cubes.f90`, found live in the Jennifer v2 verification run.**
`clean_complex` computes `comp_rm_refined_p` (a required output) and
`rmclean_cubes.f90` received it into a local variable — but never read
it again. The program's own three `derotate_to_lsq_ref` calls (compute
frame -> report frame, guarded by `lsq_ref_report.ne.lsq_ref_compute`)
all used `rm_samp` uniformly, including for `comp_re_p`/`comp_im_p` --
silently reintroducing the exact ~dRM/2-scale location error
`comp_rm_refined` was built to eliminate, for every run where
`lsq_ref_compute` genuinely differs from `lsq_ref_report`. This was
dormant under the OLD default (`lsq_ref_compute_mode=native`, which for
most cubes made `lsq_ref_compute` equal `lsq_ref_report`'s own default
of `0`, so `ref_diff=0` and the choice of location array never
mattered) -- it went live the moment this session switched the default
to `mid` (decision 8's own UPDATE above), since `mid` deliberately picks
a nonzero `lsq_ref_compute`.

Found while investigating a user question about the Jennifer v2 run's
own flux normalization (comparing CLEAN/RESTORED amplitude against the
dirty cube): attempting to numerically verify the exact `dirty =
[comp(*)RMSF_dirty] + resid` identity (provable by induction on
`clean_complex`'s own subtraction loop) from the OUTPUT cubes alone kept
failing by a large margin. Root cause traced directly to this gap --
reconvolving `comp_re`/`comp_im` using `rm_samp(imax)` as each
component's location (the same imprecise convention the derotation call
itself used) breaks the delicate near-total cancellation this pixel's
12 overlapping components require to reproduce the (heavily
depolarized) dirty amplitude, since the RMSF is oscillatory enough that
even sub-grid-cell location error matters when many large terms need to
cancel to a small net value.

**Fix**: `rmclean_cubes.f90`'s `derotate_to_lsq_ref` call for
`comp_re_p`/`comp_im_p` now passes `comp_rm_refined_p` as the location
array instead of `rm_samp` -- exactly the pattern
`tests/test_rmclean_lsqref_flex.f90` already validated at the
`rmclean_mod` library level (`ipeak`/`rm_found = comp_rm_refined(ipeak)`),
just never connected into the production pipeline. `resid_re_p`/
`resid_im_p` and `restored_re_p`/`restored_im_p` correctly keep
`rm_samp` -- both are genuine regular-grid functions (the beam is
evaluated AT `rm_samp(j)` for the residual; `restore_clean`'s own FFT
convolution produces a value that genuinely lives AT `rm_samp(j)` after
smoothing), unlike the raw component map, whose per-bin flux is
FLUX-WEIGHTED-AVERAGE-located, not grid-centred.

**Verified**: full regression suite still 109/109 (this branch was
previously untested at nonzero `ref_diff` with a tolerance tight enough
to catch it). Direct before/after comparison on section 29's own
`lsq_ref_mode=mid` test cube (rebuilding pre-fix and post-fix binaries,
running both against the identical input cubes): `CLEAN.PHA` changes by
a real amount (max `0.111` rad, mean `0.0097` rad over 701 nonzero
bins), `CLEAN.AMP` is bit-identical to float32 epsilon (`5.96e-8` --
confirms a pure phase correction, no amplitude side effect), and
`RESTORED.PHA` is EXACTLY unchanged (`0.0` -- confirms the fix is
correctly scoped to `comp_re_p`/`comp_im_p` only, as intended). Not yet
re-run against the full Jennifer dataset (multi-hour real-data job) --
the fix's correctness is established at the mechanism level (matching
the already-validated library-level test pattern) and confirmed to have
a real, correctly-scoped, nonzero effect; a full Jennifer re-run would
only add production-scale confirmation, not new evidence of
correctness.

**Addendum — `suggest_lsq_ref_compute`: the cheapest safe `lsq_ref`
choice, and why a plausible-sounding alternative is actually poor.** The
user asked directly, correctly declining to just accept "band centroid":
is the midpoint of the channel set's own l_sq extent really optimal, and
shouldn't centring on the LOWEST-frequency band's own centroid be
better (that band drives `suggest_drm`'s tightest per-band
constraint)? Checked numerically on the same imbalanced P-band(61ch)/
L-band(121ch) combination used throughout this plan (`suggest_drm`,
`oversample=15`, `rm_span=450`):

| `lsq_ref` choice | value | `dRM≤` | `nrm` |
|---|---|---|---|
| `0` | `0` | `0.0948` | `4748` |
| channel-count-weighted mean | `0.3771` | `0.1440` | `3127` |
| P-band (lowest-freq) centroid | `1.0011` | `0.1109` | `4060` |
| L-band centroid | `0.0626` | `0.1005` | `4479` |
| **midpoint of overall extent** | `0.5806` | `0.1999` | **`2253`** |

Centring on the lowest-frequency band is barely better than `lsq_ref=0`
(`4060` vs `4748`, ~15%) — it only shrinks that ONE band's own offset,
while leaving the OTHER band almost as exposed as `lsq_ref=0` did,
merely relocating the worst case rather than reducing it.
`suggest_drm`'s bound is `max_k|l_sq(k)-lsq_ref|` — a max over
the WHOLE channel set — and for any fixed set of values this maximum is
determined entirely by the two EXTREME values (`min(l_sq)`,
`max(l_sq)`): every other channel, however many there are or wherever
they sit between the extremes, can never exceed whichever extreme is
farther from `lsq_ref`, so it never enters the max. Minimizing
`max(|max(l_sq)-lsq_ref|, |lsq_ref-min(l_sq)|)` is solved exactly by
equalizing the two terms — the midpoint, `(min(l_sq)+max(l_sq))/2` — a
classic 1D minimax/Chebyshev-centre result, not an approximation.
Channel count (per-band or overall) provably does not enter this
criterion at all — confirmed by construction, not just observation:
the weighted mean (`0.3771`) sits far from the midpoint (`0.5806`)
specifically because it's pulled toward L-band's own values by its 2x
channel-count advantage, yet the midpoint still beats it.

Channel count DOES matter for a different, unrelated question — the
actual achievable statistical precision (SNR-driven, the usual
RM-synthesis noise-propagation relations) — but per finding 5 above,
that is now fully decoupled from `lsq_ref` choice once `comp_rm_refined`
and `derotate_to_lsq_zero` are used correctly: `lsq_ref` affects ONLY
grid cost, never precision, so there is no need to trade one against
the other.

New public `suggest_lsq_ref_compute(l_sq, nchan, lsq_ref_opt)`:
`lsq_ref_opt = 0.5*(minval(l_sq)+maxval(l_sq))`. New test, `tests/
test_optimal_lsq_ref.f90` (`run_tests.sh` section 25): confirms the
formula directly, and confirms the midpoint's own `suggest_drm`
bound is at least as cheap as every alternative above, on this
deliberately-imbalanced combination where a wrong intuition looks most
plausible. Full regression: 72/72 (up from 65/65: this test's own 6
checks).

**Addendum — naming cleanup and a generalized reporting reference, at
the user's own explicit request** ("be correct", "be organised — clean
in variable naming and scope", "separate compute efficiency choices from
'reporting to community' choices"). Three changes, all covered by the
new "Choosing parameters" section above:

1. **Renamed `lsq_ref` → `lsq_ref_compute` everywhere it means the
   internal computation reference** (`compute_dirty_rmbeam_direct`,
   `build_rmsf_offset_table`, `suggest_drm`,
   `suggest_lsq_ref_compute`'s own output) — a single bare name was
   starting to do double duty for two unrelated concerns once a second,
   reporting-side reference entered the picture; renamed throughout via
   a word-boundary rename (verified it left `lsq_ref_dp` and similar
   longer identifiers untouched) plus a full rebuild/regression, not a
   partial or comment-only rename.
2. **`derotate_to_lsq_zero` generalized to `derotate_to_lsq_ref(rm_samp,
   nrm, lsq_ref_compute, lsq_ref_report, re_in, im_in, re_out, im_out)`**
   — the same exact pointwise identity as before, just parameterized by
   an arbitrary target `lsq_ref_report` instead of a hardcoded `0.0`
   (`lsq_ref_report=0.0` reproduces the old behaviour exactly). New
   public `suggest_lsq_ref_report(l_sq, nchan, lsq_ref_report)`: the
   Brentjens & de Bruyn (2005)-style reporting reference — the
   (weighted) mean, which this module's own current equal-per-channel
   `(1/nchan)` weighting convention makes a plain arithmetic mean (would
   need to become a genuine weighted mean if this module ever gains
   real per-channel noise weights — flagged directly in its own doc
   comment so that future addition doesn't silently leave it stale).
3. **Why B&dB's centroid and `suggest_lsq_ref_compute`'s own midpoint
   are NOT competing answers to the same question**, worked through with
   the user directly: B&dB's centroid minimizes the STATISTICAL
   COVARIANCE between a jointly-fitted RM and chi0 (a classic linear-
   regression centering result — centering a regressor at its own
   weighted mean decorrelates the fitted intercept and slope) — a
   concern about NOISE, most directly relevant to a QU-fitting-style
   joint estimation, or to what reference to QUOTE final results at for
   statistically clean reported uncertainties. `suggest_lsq_ref_compute`
   minimizes `suggest_drm`'s worst-case grid-stability bound — a
   concern about NUMERICAL COST for the iterative, grid-based CLEAN
   algorithm specifically, with no noise or covariance anywhere in it.
   They coincide only for a single, symmetric, evenly-sampled band;
   `tests/test_optimal_lsq_ref.f90`'s own imbalanced P+L combination
   shows them diverging (`0.5806` vs `0.3771`).

`tests/test_rmclean_lsqref_flex.f90` extended (not just renamed) to
directly exercise the new reporting-reference capability: for each
`lsq_ref_compute` case, chi0 is now checked at BOTH `lambda_sq=0` and at
`suggest_lsq_ref_report`'s own value, against an INDEPENDENTLY derived
closed-form expectation (`chi0_true + RM_true*lsq_ref_report`, not
derived from `derotate_to_lsq_ref` itself, so this is a real test of the
formula, not a tautological confirmation) — confirmed both
`lsq_ref_compute` choices agree with each other and with the closed form
exactly. Full regression: 74/74 (up from 72/72: two new checks per
compute case).

**Addendum — one more naming pass, and a real factual correction, both
from direct user pushback:**

1. **The `suggest_*` naming convention, adopted after the user pointed
   out the previous names weren't self-explanatory as a family:**
   `optimal_lsq_ref_midpoint` → `suggest_lsq_ref_compute`,
   `lsq_ref_report_centroid` → `suggest_lsq_ref_report`,
   `required_drm_nyquist` → `suggest_drm` (the user's own specific
   complaint: "required_drm_nyquist" has no verb and doesn't parse
   cleanly as a name — true of the other two as well, just less
   obviously). All three now share one pattern: "don't know what to set
   `<X>` to? Call `suggest_<X>`." The underlying math (midpoint / mean /
   Nyquist-type bound) moved entirely into each routine's own doc
   comment — a reader only needs the prefix to know a routine is
   optional and fills in a value, not the specific method, to use it
   correctly. Renamed via a word-boundary regex across `src/rmclean.f90`
   and all four `tests/*.f90` files plus this plan, verified with a
   second pass for markdown line-wrap artifacts the regex couldn't
   catch (`required_drm_\nnyquist` split across lines in two places,
   caught and fixed by hand) — a full rebuild and 74/74 regression after
   every rename pass, not just a text substitution.
2. **Factual correction to `suggest_lsq_ref_report`'s own doc comment:**
   it previously claimed "no separate per-channel noise-weight concept
   anywhere in rmclean_mod today" as if masking were entirely absent
   from this project — wrong, and corrected after the user pointed out
   `MASK.CUBE.FITS`'s own per-pixel channel masking exists precisely for
   this. Verified directly against `rm_synthesis_mod.f90`'s own source
   (`wts_gpu(ipix, iz) = real(mask_tile(src_idx), sp)`) rather than
   re-asserting the correction from memory: this project's masking is a
   genuine 0/1 weight, and nowhere in this codebase, checked directly,
   is any GRADUATED (e.g. inverse-noise-variance) weight ever used —
   masking is the entire weighting scheme that exists. A 0/1 mask
   reduces exactly to "average over the surviving channels," which is
   precisely what `suggest_lsq_ref_report` already computes.
   `rmclean_mod` itself still takes no mask argument directly (same as
   every other routine here — it only ever sees `l_sq(nchan)`), so
   correct masking is the CALLER's responsibility (build `l_sq` per
   pixel from only that pixel's own unmasked channels, e.g. from
   `MASK.CUBE.FITS`, before calling in) — but given that,
   `suggest_lsq_ref_report` is ALREADY correct for per-pixel masking
   today, not pending a future addition as the old comment implied.
   Only a genuinely graduated (non-0/1) weight, if this project ever
   adopts one, would require the formula itself to change.

**Addendum — a third naming pass (`suggest_*` → `get_*`, mode-based
interface) and the oversample-floor finding, from continued direct user
pushback:**

1. **`suggest_lsq_ref_compute`/`suggest_lsq_ref_report`/`suggest_drm`
   renamed to `get_lsq_ref_compute`/`get_lsq_ref_report`/`get_drm`.** The
   user's objection: "suggest" doesn't make clear the value is actually
   committed for the real work, only offered as advice. Checked
   established precedent rather than picking a word by feel: HDF5's own
   `H5Pget_chunk`/`H5Pset_chunk` get/set convention for tunable
   parameters, and cuFFT's own `cufftGetSize1d` (the committed answer)
   vs. `cufftEstimate1d` (a quick heuristic) — "Get" is the established
   choice for a value that gets used, not merely suggested.
2. **Dropped the earlier "optional `requested=` override" design
   entirely** — a second objection, equally valid: passing a raw number
   through an optional argument on something named `get_` is really a
   SETTER wearing a getter's name. Resolution: a getter now only ever
   computes/derives a value; a caller who wants a specific number they
   already know just writes it directly in their own code
   (`lsq_ref_compute = 0.734_sp`) — no module call needed, exactly how
   `niter`/`gain` already work, and no ambiguity about what "get" means.
3. **Introduced a `mode` argument (a required `integer, parameter` enum,
   matching this module's own pre-existing `fftw_forward`/`fftw_estimate`
   convention) plus an optional `fixed_value`**, after the user asked
   directly how to let a caller choose between a qualitative/statistical
   convention (`min_lsq`, `max_lsq`, `centroid`, `mid`, `intrinsic`) and
   an arbitrary fixed number, through one uniform interface. `fixed_value`
   is read only when `mode=..._fixed`, enforced with a clear runtime stop
   if missing — the same "refuse loudly" convention already used
   elsewhere in this project. Naming clarified along the way: Fortran's
   `parameter` attribute (a fixed, compile-time-constant integer, never
   assigned to at runtime) is a completely different concept from
   "parameter" meaning "function argument" in general programming usage —
   worth spelling out explicitly, since this was a genuine, reasonable
   point of confusion, not a design flaw.
4. **`get_drm` redesigned around the user's own simpler idea: no raw
   `dRM` input at all, only `oversample`, always applied on top of a
   mandatory, ENFORCED floor** ("why is Nyquist even an option? Should we
   not ever undersample?"). This directly motivated finding the REAL
   floor empirically rather than assuming bare Nyquist (`oversample=1`)
   was safe: it is NOT — `tests/test_drm_floor.f90` shows `oversample=1`
   recovers a wrong chi0 (not merely imprecise), while `oversample=2`
   already matches `oversample=15`'s own precision. Floor enforced at `2`;
   default `4` if `oversample` is omitted (the user's own suggested
   number, confirmed reasonable by this check). `tests/test_drm_floor_
   enforcement.f90` confirms `get_drm` actually stops (nonzero exit) for
   `oversample<2`, wired into `run_tests.sh` section 27 with intentionally
   INVERTED pass logic (a nonzero exit is the correct, expected outcome
   for that one test — a genuine, deliberate exception to this suite's
   own PASS/FAIL convention, documented as such at the call site so it
   isn't mistaken for an oversight later).
5. **`lsq_ref_compute_mid` established as the explicit, documented
   RECOMMENDED DEFAULT** (the user's own instruction) — `get_lsq_ref_
   compute`'s own doc comment now states plainly, with both halves of the
   claim verified rather than asserted: harmless for accuracy (chi0 is
   provably decoupled from `lsq_ref_compute` choice once `comp_rm_refined`
   is used) and actively smart for cost (the same minimax proof as
   before). The other five modes remain available for experimentation or
   matching an external convention, at no safety cost, since `get_drm`
   enforces its own floor against whatever `lsq_ref_compute` value
   actually results, regardless of which mode produced it.

Full regression: 78/78 (up from 74/74: two new test programs, `tests/
test_drm_floor.f90` and `tests/test_drm_floor_enforcement.f90`, sections
26/27).

### T2 — `rmclean_cubes`: standalone RM-CLEAN tool + mask-pattern cache + OpenMP

New standalone program, `src/rmclean_cubes.f90` (own Makefile target
`make rmclean_cubes`, `bin/rmclean_cubes`, wired into `scripts/
make_all.sh` and `docker/dockerfile`), driving `rmclean_mod`'s already-
hardened API (T1) against a REAL dirty AMP/PHA cube pair `rm_synthesis`
itself wrote, plus its `.MASK.CUBE.FITS`/`CHANFREQ` table (T0). Follows
`reproject_cubes.f90`/`convolve_cubes.f90`'s own standalone-tool
conventions (own key=value/`--config` parser, `--help`, own Makefile
block).

1. **Gate 0 (no re-gridding)**: this tool cannot resample the RM axis --
   `CDELT3`/`nrm` are fixed by whatever `rm_synthesis` already wrote.
   Before any CLEAN work, `get_drm(l_sq, nchan, lsq_ref_native, drm_
   required, oversample)` is compared against the cube's own `CDELT3`;
   refuses to proceed (clear error, nonzero exit) if the existing grid
   is too coarse. Uses the FULL per-run channel list (not any one
   pixel's own valid subset) — conservative and safe, since dropping
   channels only ever reduces `get_drm`'s own bound.

2. **`lsq_ref_native` vs `lsq_ref_compute` — a real design correction
   made during implementation** (the user pushed back directly: "even if
   the dirty cubes are at lsq=0, shouldn't we be able to RM-clean at any
   computationally favourable lsq_ref?"): `lsq_ref_native` is the
   reference the cube was ACTUALLY built at (read from its own new
   `LSQREF` header keyword -- see T2b below -- falling back to 0.0, with
   a printed warning, for cubes written before that keyword existed).
   Gate 0 always validates against `lsq_ref_native` specifically (sampling
   adequacy is a property of how the data was generated). `lsq_ref_
   compute` (this tool's OWN RMSF-table/CLEAN reference) is a fully free,
   independent choice (`lsq_ref_compute_mode=native|zero|mid|centroid|
   min|max|fixed`, default `native`), applied via `derotate_to_lsq_ref`
   -- an exact, lossless phase rotation of the already-sampled dirty
   spectrum (verified `|P|`-invariant, T1's own `test_rmclean_lsqref_
   flex.f90`), so choosing a non-native value costs nothing in accuracy.
   It also saves nothing in compute once the RM grid already exists
   (`build_rmsf_offset_table`/`clean_complex`/`restore_clean`'s own cost
   depends only on `nrm`/`niter`/`table_oversample`) -- the grid-size win
   a favourable reference buys only applies UPSTREAM, at `rm_synthesis`'s
   own `lsq_ref_mode` (T2b), not here.

3. **Mask-pattern pre-scan + hash-bucketed `rmsf_table_t` cache**
   (decision 10): `build_mask_pattern_cache` scans the full mask cube
   ONCE, serially, before any parallel region; canonicalizes each
   pixel's own valid-channel bit pattern (the mask row itself, already
   0/1 bytes -- no further packing needed), hashes it (64-bit FNV-1a),
   and builds ONE `rmsf_table_t` per DISTINCT pattern via an open-
   addressing (linear-probing) hash-bucket table for O(1) amortized
   lookup -- a plain linear scan across cache entries was rejected as
   too slow for large pixel counts. Collision-safe: every lookup does a
   full byte-for-byte pattern compare on a hash match, never trusts the
   hash alone. Capped at `mask_pattern_cache_max` (default 4096)
   distinct patterns; past the cap, a pixel falls back to a one-off
   throwaway table (safety valve, not a correctness issue). Verified:
   `mask_pattern_cache_max=4096` (744 distinct patterns cached) vs `=0`
   (every pixel a throwaway) vs `=50` (mixed hit+overflow) all produce
   BIT-IDENTICAL output on a synthetic cube with 744 distinct per-pixel
   masking patterns (`tests/run_tests.sh` section 28).

4. **OpenMP parallelism**: one `!$omp parallel do collapse(2) schedule(
   dynamic)` over `(ix,iy)`, calling `clean_one_pixel` per pixel. Every
   per-pixel scratch array (`dirty_re_p`, `comp_re_p`, etc.) is a genuine
   LOCAL automatic variable inside that subroutine — thread-safe purely
   by Fortran's own call-stack semantics, no explicit `private()`
   bookkeeping needed for them. The shared FFTW restore plan (built ONCE
   outside the parallel region) and the mask-pattern cache (fully built
   serially beforehand) are both read-only inside the parallel region.
   Verified: `OMP_NUM_THREADS=1` vs `=4` produce BIT-IDENTICAL output on
   the same test cube.

5. **A real robustness bug found and fixed via testing, not by
   inspection**: `write_output_cube`/`read_cube`/`read_mask_cube`
   originally called `FTGPVE`/`FTGPVB`/`FTPHPR` immediately after
   `FTOPEN`/`FTINIT` WITHOUT checking the returned status first -- if
   the underlying call failed (e.g. `FTINIT` refusing to create a file
   that already exists on disk, exactly what happens re-running this
   tool without cleaning up a previous run's own outputs), every
   subsequent CFITSIO call operated on a unit CFITSIO never actually set
   up: not a clean no-op but undefined behaviour, reproduced as a real
   SIGSEGV. Fixed by checking status immediately after every `FTOPEN`/
   `FTINIT` and bailing out cleanly. Caught by `tests/run_tests.sh`
   itself segfaulting on a second consecutive run (stale output files
   from the first run) -- section 28 now also cleans up its own `rmc_*`
   outputs up front, and the fix itself is a genuine correctness
   improvement independent of the test.

### T2b — `rm_synthesis`'s own `lsq_ref_mode` option

Added alongside T2 (the user: "instead of hard coding lsq=0, we should
add flexible option for rmsynthesis - following the same strategy as in
rmclean"). `rm_synthesis_mod.f90`'s `extract_general_setup` was
unconditionally at `lambda_sq=0` (Finding A, confirmed directly: its own
`phi_tmp=omega*t(kk)` used raw `L_sq`, no subtraction). Now: `cfg%lsq_ref_
mode` (`zero`|`mid`|`centroid`|`min`|`max`|`fixed`) + `cfg%lsq_ref_fixed_
value`, resolved once via a new `compute_lsq_ref` function (duplicated
from, not `use`-coupled to, `rmclean_mod`'s own `get_lsq_ref_compute` --
this module is the older, heavily production-tested core; the newer,
algorithm-specific `rmclean_mod` should not become a hard dependency of
it). Baked into `extract_general_setup`'s own `cos_arr`/`sin_arr`
template construction (`phi_tmp = omega*(t(kk)-lsq_ref)`) -- the ONLY
call site (`rm_synthesis.f90:2178`), so every downstream consumer
(`extract_general`'s own dot products) automatically inherits whichever
reference was chosen, no other code changes needed. The actual value
used is recorded in a new `LSQREF` header keyword on both AMP/PHA cubes,
so `rmclean_cubes` (T2) never has to assume one. Default `zero` preserves
every existing cfg file's behaviour exactly (verified: full 78-test
suite unchanged before this addition, still green after).

Full regression after T2+T2b: 84/84 (up from 78/78 — new section 28,
`rmclean_cubes` end-to-end against a real `lsq_ref_mode=mid` cube: LSQREF
header round-trip, Gate 0, both known point sources recovered via
`check_rm_peak.py`, mask-pattern cache correctness).

### T3 — model-based matched-filter peak refinement

Grew out of an extended dialogue with the user pressure-testing T2's own
Gate 0: why does `get_drm` demand ~22-44x oversampling relative to
`fwhm` at `lsq_ref=0` (real numbers, thesis P-band, `nu_c=300 MHz`,
`dnu=30 MHz`, 61 channels: `dRM(oversample=4)=0.356 rad/m^2` against
`fwhm=7.826 rad/m^2`)? Root cause, derived and then verified against the
actual code: `clean_complex`'s own `peak_interp_parabolic` fits the CLEAN
peak's complex value (`re_at_peak`/`im_at_peak`) via a local parabola
fit **directly through the raw, stored `Re`/`Im` samples** -- and those
samples carry a fast, `lsq_ref`-dependent carrier (`R(Δ)` decomposes into
a slow, `lsq_ref`-independent envelope times a unit-modulus carrier at
`max_offset = maxₖ|l_sq(k)-lsq_ref|`), so a local low-order fit to them is
only valid if the OUTER grid already satisfies that same demanding bound.
The peak's *location* step (a separate parabola on log-magnitude) was
never the problem -- magnitude carries no carrier.

**The fix**: since `R(Δ)` is analytically known exactly, at any offset,
from `l_sq` alone (no raw per-channel Q/U data needed), replace the
raw-sample parabola with a local matched-filter search: for a fine grid
of trial offsets spanning `[rm_samp(imax)-drm, rm_samp(imax)+drm]`, fit
the one-complex-unknown amplitude `A = Σⱼconj(R(Δⱼ))·dataⱼ / Σⱼ|R(Δⱼ)|²`
against the 3 nearest *stored* coarse samples, keeping the trial that
maximizes `|Σⱼconj(R)·data|²/Σⱼ|R|²`. Two new subroutines, `rmsf_point_
direct` (O(nchan), single-offset exact `R(Δ)` -- NOT `compute_dirty_
rmbeam_direct`, which needs pre-built `cos_arr`/`sin_arr` templates and
always loops over the full output array even for one point) and
`refine_peak_matched_filter` (the search itself), both in `src/
rmclean.f90`.

**Landed in 3 independently-tested phases, per the user's own explicit
"validate all claims before production" policy**:
- Phase 0 (`tests/test_matched_filter_refine.f90`, zero production code
  touched): free-standing prototype, validated against 2 independent
  point-source scenarios (`RM=50/chi0=0.3` and `RM=-27.3/chi0=-1.1`) at
  4 sampling densities each, including sampling far below `get_drm`'s
  own floor. Found and fixed a real sub-issue along the way: the LOCAL
  search grid's own density must ALSO track `max_offset` (not a fixed
  constant) -- a coarse discrete search over `R(Δ)`, itself
  fast-oscillating, can lock onto the wrong nearby local maximum of an
  oscillatory objective (confirmed directly: chi0 error vs. search
  density was non-monotonic in [40,400], only converging reliably past
  ~500 with the original nanchor=5 config).
- Phase 1: ported into `rmclean_mod` as new, additive public
  subroutines; `clean_complex` itself still 100% unchanged; promoted
  test wired into `tests/run_tests.sh` (calls the production
  subroutines directly, confirmed identical to the validated prototype).
- Phase 2: swapped `clean_complex`'s own internals (`peak_interp_
  parabolic` call + `re_at_peak`/`im_at_peak` derivation replaced by one
  `refine_peak_matched_filter` call; nothing else in the Hogbom loop
  changed). Required a signature change (`clean_complex` now takes
  `l_sq`/`nchan`/`lsq_ref_compute`, needed by the refinement step but not
  by the old parabola) -- updated at all 5 call sites (`tests/thesis_
  scenario_rmclean.f90`, `tests/test_drm_floor.f90`, `tests/test_rmclean_
  lsqref_flex.f90`, `src/rmclean_cubes.f90`'s 2 call sites).

**A genuine, expected consequence, not a regression**: `tests/test_drm_
floor.f90` was originally written to demonstrate `peak_interp_parabolic`'s
own Nyquist-floor requirement (`oversample=1` gives a confidently WRONG
chi0). Since T3 removes that specific requirement, `oversample=1` now
ALSO passes -- through the FULL Hogbom CLEAN loop, not just Phase 0's
isolated check. The test's own expectation was inverted (documented in
its own updated header) rather than treated as a failure -- direct,
additional confirmation of the fix, at real production code, not just
the prototype.

**A real, measured performance cost, found and partially addressed**:
first working version (5 anchors, 200 samples/carrier-cycle,
empirically-generous but never load-tested) took `rmclean_cubes` from
sub-second to ~1m40s on a 32x32/1024-pixel test cube (niter=200) -- niter
x npixels calls to a per-call-cheap routine still add up in aggregate.
Retuned with evidence: 3 vs 5 anchors gave BIT-IDENTICAL results in
every noiseless scenario checked (the 2 extra anchors add no
discriminating power); a direct sweep against both Phase-0 scenarios
found the actual failure boundary at ~7-8 samples/cycle -- 50/cycle
(chosen) is a checked >=6x margin above that boundary, not a guess.
Result: 1m40s -> 31.6s on the same test cube (~3.2x), full 8-check
Phase-0 harness and full regression suite still green at the new
settings. **Not fully closed**: 31.6s for 1024 pixels does not scale to
a real million-pixel cube (~8+ hours extrapolated) -- flagged explicitly
to the user as a follow-on ticket (T3c: a coarse-to-fine two-stage
search -- coarse pass at enough points/cycle to avoid aliasing between
carrier cycles, fine pass only within the winning cycle -- should cut
cost from O(M) to roughly O(sqrt(M)), but needs its own validation cycle
before landing, not rushed into this same ticket).

Full regression after T3: run via `tests/run_tests.sh` (new section 28,
`RM-CLEAN matched-filter peak refinement`, ahead of the renumbered
`rmclean_cubes` section 29); all RM-CLEAN sections (23-29) and the
pre-existing full suite green at the retuned settings.

### T3c — tiered fast-path refinement (supersedes the coarse-to-fine sketch)

The design came out of an extended clarification dialogue with the user
that also caught two real defects in T3's first implementation along
the way: (a) the search window was `+/-drm` around the discrete peak
when `+/-drm/2` is all the true peak can occupy (the discrete maximum
bin is at most half a cell from the true peak, else a neighbour would
have won) -- fixed, halving the search for free; (b) the search-density
constants (`samples_per_cycle=50`/`m_floor=100`) were empirically swept
but not physically anchored -- the user's requirement: "choose them
based on real physical considerations, not arbitrarily."

**The tiered design (the user's own, refined jointly):**

- **Fast path, runs every call**: peak LOCATION from the log-magnitude
  parabola (`peak_interp_parabolic`'s own location step -- sound at
  resolution-level sampling, magnitude carries no carrier). At that
  FIXED location, the one-complex-unknown amplitude is solved
  closed-form against the (up to 3) nearest stored anchors -- the same
  matched-filter formula as the search's inner step, shared via one
  internal `eval_trial` helper so the tiers cannot drift apart
  numerically. Crucially this is NOT a single-point division (an
  explicitly rejected earlier design): a single-anchor division always
  "fits" perfectly (2 unknowns, 2 equations) and can never signal a
  wrong location; the >=2-anchor fixed-location fit is over-determined,
  so its LEFTOVER misfit is a genuine self-consistency diagnostic.
- **Escalation criterion (data-driven, user's own choice)**: accept the
  fast path when `sqrt(leftover/nanchor) <= nsigma * noise_rms`, with
  `noise_rms` the same per-iteration `rms_about_mean` estimate
  `clean_complex`'s own stopping criterion already uses (no separate
  estimator), and `nsigma` a user-facing knob (`nsigma_refine` optional
  argument on `clean_complex`, `refine_nsigma=` config key on
  `rmclean_cubes`), default 3.0 -- a conventional n-sigma cut,
  confirmed against the test program's own noisy scenario rather than
  only asserted. Beyond threshold: the full T3 local search runs for
  that iteration, unchanged.
- **Physically-anchored search constants** (replacing the arbitrary
  values): the search samples the matched-filter STATISTIC, whose
  fastest oscillation is TWICE the carrier rate (squared magnitude
  doubles frequency), carrier rate set by `max_offset` exactly as in
  `get_drm` -- so `cycles_in_window = 2*window*max_offset/pi`, at 50
  samples per statistic cycle (validated margin ~12x above the
  empirically-located wrong-local-maximum failure boundary of ~4 per
  statistic cycle). `m_floor` cut from 100 to 20, now tied to
  `table_oversample`'s own default meaning elsewhere in this project
  ("how finely to resolve within one native grid cell") rather than
  being an independent invented constant.

**Evidence** (`tests/test_matched_filter_refine.f90`, rewritten for the
tiers -- per scenario x density it now checks the forced-search path
(`nsigma=-1`), the forced-fast path (`nsigma=huge`), and the tiered
default (`nsigma=3`), all against truth, plus two dedicated mechanics
checks and a noisy-default check):

- The fast path ALONE is already accurate on clean single-component
  data even at `fwhm/2` sampling (chi0 error ~0.01 rad) -- the
  parabola's LOCATION was never the weak point of the old method, its
  raw-sample complex-value fit was; location + model-anchored
  closed-form amplitude is enough.
- Tier mechanics, same threshold and claimed noise floor, only the data
  differing: clean single component -> fast path ACCEPTED (no search);
  two blended components ~0.8 fwhm apart -> escalation FIRED. Both
  confirmed.
- Noisy single component (deterministic LCG noise, sigma=2% of
  amplitude) at default `nsigma=3`: recovery within tolerance, fast
  path taken.
- **Performance, measured**: `rmclean_cubes` on the same 1024-pixel
  test cube (niter=200) -- ~1m40s (T3 first working version) -> 31.6s
  (T3 retuned) -> **0.46s (T3c tiered)**, with recovered source
  RM/amplitude/chi0 unchanged against truth. The originally-sketched
  coarse-to-fine two-stage search is superseded: the fast path skips
  the search entirely on the (overwhelmingly common) clean iterations
  rather than making it cheaper, and the full search remains, exact
  and unchanged, where the data genuinely needs it.
- **A real input-contract boundary, found by the regression suite, not
  by design**: `tests/test_drm_floor.f90`'s own `oversample=1` case
  (band-CENTROID reference) failed again under the tiered default
  after passing under T3's pure search. Root-caused, not papered over:
  at a centroid reference `max_offset ~= span/2`, so the carrier and
  RESOLUTION scales coincide -- `samples/fwhm = oversample/2`, and
  `oversample=1` is only 0.5 samples/fwhm, BELOW resolution Nyquist,
  i.e. outside the very input contract T3b's Gate 0 enforces. The fast
  path's log-magnitude location step has nothing valid to work with on
  such a grid, and early iterations' sidelobe-dominated rms makes the
  escalation criterion too lenient to catch the miss (the pure search,
  being phase-sensitive across the whole +/-dRM/2 window, happened to
  survive it). The test was updated to expect failure again -- with its
  doc comment now stating the NEW reason (sub-resolution input, the
  exact failure mode Gate 0 exists to refuse) rather than the old
  carrier/parabola reason -- and its `oversample=2` case (1 sample/
  fwhm: above Gate 0's hard floor, below its default 2) confirmed
  still recovering correctly, evidence the default carries genuine
  margin. Module-level contract made explicit by this: `clean_complex`
  requires resolution-adequate input; `rmclean_cubes`' Gate 0 is what
  guarantees it in production.

### T3b — Gate 0 recast as an RM-resolution criterion

`rmclean_cubes`'s Gate 0 no longer validates against `get_drm`'s
`max_offset`-based (carrier-driven, `lsq_ref`-dependent) bound -- with
T3/T3c in `clean_complex`, the stored grid never needs to resolve the
fast Re/Im carrier at all. The gate now enforces the one requirement
that genuinely remains: `|CDELT3| <= fwhm/min_samples_per_fwhm`, with
`fwhm` from `compute_rmsf_fwhm_multiband` -- an `lsq_ref`-INDEPENDENT
quantity set by the lambda^2 span alone (always the DATA's own fwhm,
never the `restore_fwhm` override, which only shapes the restoring
beam). New config key `min_samples_per_fwhm` (default 2.0, hard floor
1.0) replaces the retired `oversample=` key entirely. Verified both
ways manually: the real test cube (`CDELT3=0.5`, `fwhm=11.41`, ~22.8
samples/fwhm) passes; the same cube with `CDELT3` forged to 8.0 is
refused at the default (needs `<=5.71`) and accepted at
`min_samples_per_fwhm=1` (needs `<=11.41`) -- and the AMP/PHA
geometry-consistency check was confirmed to fire before Gate 0 when
only one of the pair is altered. `get_drm` itself (and its
`oversample>=2` floor) is UNCHANGED at the module level -- still the
correct tool for its remaining jobs (T3c's own internal search-density
sizing, and `rm_synthesis`-side grid planning for anyone who wants
carrier-resolved dirty cubes); only `rmclean_cubes`'s gate stopped
using it.

### T4 — memory-budgeted, threaded block I/O for `rmclean_cubes`

`rmclean_cubes` was "Stage A" (T2's own scope): whole-cube-in-memory,
correct but unbounded in RAM, with no throughput mechanisms for
production-scale cubes. T4 ports the SAME scheme `rm_synthesis` already
uses in production (`plan_tile`/`io_read_threads`/`io_write_threads`/
`io_overlap` in `rm_synthesis_mod.f90`/`rm_synthesis.f90`), adapted from
rm_synthesis's 2 named pixel-cube outputs (AMP/PHA) to `rmclean_cubes`'s
2 inputs + 6 outputs (CLEAN/RESID/RESTORED x AMP/PHA), done as four
independently-tested sub-tickets, T4a-T4d, in order. All work is in
`src/rmclean_cubes.f90` only (`rmclean_mod`/`src/rmclean.f90` untouched).

**Deliberate divergence (not an oversight):** the mask cube + FNV-hash
mask-pattern -> `rmsf_table_t` cache stay whole-cube-resident, not
tiled -- the cache needs one global pre-scan across the whole image
before any pixel is CLEANed (pixels anywhere can share a table), unlike
rm_synthesis's own per-tile mask. Tiny relative to the float cubes (1
byte/voxel x `nchan` vs. 4 bytes/voxel x `nrm` x 8 arrays), so keeping
it resident costs little.

**T4a (core tile geometry):** `plan_rmclean_tile` ports `plan_tile`'s
RA-strips-first auto-tiling + safety-shrink policy verbatim (same
`mem_frac_ram`/`tile_ra`/`tile_dec`/`tile_auto` key names/defaults),
budgeted against rmclean_cubes' own per-tile-pixel cost (2 input + 6
output RM-depth arrays, `4*8*nrm` bytes/pixel). Read/compute/write per
tile mirrors reproject_cubes.f90's own 3-phase cadence (this ticket's
own first phase, before any of T4b-d's threading lands on top).
`read_amp_pha_tile`/`open_output_cube`/`write_output_tile`/
`close_output_cube` replace the old whole-cube `read_cube`/
`write_output_cube`; `clean_one_pixel` gained a tile-local
`(ix_l,iy_l)` + tile-origin `(ix0,iy0)` signature to address the
now-tiled float cubes locally while still addressing the
whole-cube-resident mask/cache globally.

**Real finding, not anticipated:** comparing a forced small-tile run
against the default (single-tile) run showed small but genuine
numerical differences (~1e-4 in AMP, larger in PHA at near-zero AMP).
Root-caused, not assumed: bit-identical at `-O0` (no
auto-vectorization) across the same two tile configurations, confirmed
directly by rebuilding both configurations at `-O0` and diffing --
`gfortran -O3 -march=native`'s auto-vectorization reassociates
`clean_complex`'s own tiered-refinement threshold comparison
(T3c) differently depending on the runtime memory alignment of its
stack-allocated arguments, which genuinely differs between tile sizes
(different automatic-array footprints in the enclosing call frames) even
though every value going INTO `clean_complex` is bit-identical (verified
directly by instrumenting `clean_one_pixel`). Not a tiling logic bug --
a pre-existing floating-point-reassociation characteristic of the
optimized build that T4a's own tile-shape change was simply the first
thing to expose (the previous test suite never compared two runs that
processed pixels in a different order). Consequence: T4a's own
regression test (`tests/check_tile_consistency.py`) compares with a
tolerance (AMP atol 1e-2, PHA atol 0.05 rad only where AMP clears a
floor), not byte-identity -- the right standard for an iterative,
threshold-branching algorithm, unlike the purely linear resampling
`reproject_cubes`/`convolve_cubes` compare byte-identically.

**T4b (`io_read_threads`):** `split_range_rmclean` (same base/rem
convention as `split_channels_across_threads`) splits a tile's own
`nrm` range across `io_read_threads_eff` threads, each opening its OWN
readonly `FTOPEN` handle (300+t/400+t unit numbers) and reading its own
disjoint RM-slice via `FTGSVE` -- safe because readonly opens are
exempt from the same-file handle-aliasing hazard read-write handles hit
(see T4c). Verified bit-identical to the serial baseline at
`io_read_threads=1,2,4`.

**T4c (`io_write_threads`):** ports `write_rm_chunk_raw`
(`write_rm_chunk_raw_rmclean`) + `host_is_big_endian`/
`swap_bytes_r4_inplace` verbatim -- raw stream writes at computed byte
offsets, bypassing CFITSIO's `ftpsse` for the pixel data entirely,
since N CFITSIO handles opened read-write on the same file alias onto
one shared internal buffer (confirmed via rm_synthesis's own T4
postmortem). `open_output_cube` fetches each output's data-start byte
offset via `FTGHAD` once and closes the CFITSIO handle immediately
(the T6 lesson: closing late risks a stale internal buffer flush over
raw-written bytes).

**Two real, non-obvious bugs found and fixed while validating T4c**
(both caught by repeated-run stress testing, not a single pass --
the first only manifested probabilistically):

1. Concurrent `open(newunit=...)` calls from different OpenMP threads
   (each opening the SAME output file path for its own disjoint byte
   range) intermittently produced a corrupted/all-zero output cube --
   the Fortran standard does not guarantee any I/O statement is safe to
   call concurrently without explicit synchronization, so this was an
   unsynchronized use of a construct never promised thread-safe, not a
   gfortran-specific defect. First fix: fixed, pre-assigned per-thread
   unit numbers instead of `newunit=` -- worked, but only via manual
   cross-file bookkeeping (every other unit-number range in the file
   has to stay disjoint from it by inspection). Superseded (per the
   user's own explicit follow-up request) by wrapping just the
   `open(newunit=u,...)` call itself in a named
   `!$omp critical (raw_write_open_lock)` -- genuinely unique real
   `newunit=` semantics, no manual range to maintain, only the brief
   allocation step serialized (the write/close per thread still runs
   fully in parallel). Same fix applied to `rm_synthesis_mod.f90`'s own
   `write_rm_chunk_raw` (io_write_threads' original implementation,
   which has the identical pattern). 20 repeated `io_write_threads=4`
   runs against each tool post-fix, 0 mismatches.
2. **This system's installed libcfitsio's `FTGHAD` writes only the
   LOWER 32 bits of its 3 output arguments**, leaving the upper 32 bits
   of an `integer(kind=8)` receiving variable completely untouched --
   confirmed with a minimal standalone reproducer entirely outside this
   codebase (a bare `FTINIT`/`FTPHPR`/`FTGHAD` program: initializing the
   receiving variables to `-1_8` before the call reproducibly yields
   `datastart = -2^32`; initializing to `0_8` first yields the correct
   value). `open_output_cube`'s local `headstart`/`local_datastart`/
   `dataend` are automatic (stack) variables with no guaranteed initial
   value, so this showed up as a genuinely intermittent (stack-content-
   dependent) wrong byte offset roughly half the time. Fixed by
   zero-initializing all 3 immediately before every `FTGHAD` call.
   **`rm_synthesis.f90`'s own `FTGHAD` call (io_write_threads' own
   original implementation) has the exact same latent exposure** --
   never observed to misbehave there only because its receiving
   variables are main-program-scope locals that happen to land in
   zero-initialized static storage on this platform/compiler (not the
   stack), which is incidental storage-class behaviour, not a
   guarantee. Fixed there too, same zero-init, confirmed with the full
   regression suite still green afterward. Verified via 20 repeated
   `io_write_threads=4` runs post-fix (0 failures; the same loop
   pre-fix failed roughly half the time) -- `tests/run_tests.sh`'s own
   T4c section repeats this check 5x per thread count rather than
   trusting one pass, given the bug's own probabilistic nature.

**T4d (`io_overlap`):** ports the pthread interface bindings
(`c_pthread_create`/`c_pthread_join`) and double-buffer/dispatch
machinery (`tile_write_job_t`, `tile_write_dispatch_async`,
`tile_write_thread_entry`, `do_tile_write`) from
`rm_synthesis_mod.f90`, generalizing `tile_write_job_t`'s named
`unit_amp`/`unit_pha` fields to a 3-pair array (`re_ptr(3)`/`im_ptr(3)`,
addressing clean/resid/restored in a fixed order, each written via the
already-existing, already-tested `write_output_tile`). The two join
rules ported exactly as rm_synthesis has them (both load-bearing, not
stylistic): buffer-reuse join (before repopulating slot `cur_slot`,
join that slot's own prior pending write) and handle-safety join
(before dispatching ANY new write, join BOTH slots unconditionally
first, since all 6 output files' handles are shared across slots
regardless of which slot is being written). Output-tile storage is
double-buffered via a trailing slot dimension (`clean_re_buf(:,:,:,2)`
etc.) with `clean_re_tile` etc. as POINTERs re-targeted at
`cur_slot+1` every tile, rather than rm_synthesis's own named `_s0`/
`_s1` array pairs. Verified bit-identical to `io_overlap=n` (5 reps),
and — the real gate — the FULL combined stress case (forced small
non-full-width tiles + `io_read_threads=3` + `io_write_threads=3` +
`io_overlap=y` together, not each mechanism only in isolation) bit-
identical to the small-tile-only reference (5 reps).

**Evidence:** `tests/run_tests.sh` section 29 (rmclean_cubes end-to-end)
grew four new checks (T4a small-tile tolerance comparison, T4b
`io_read_threads` sweep, T4c `io_write_threads` sweep x5 reps, T4d
`io_overlap` alone x5 reps + the combined stress case x5 reps) — full
suite 93/93 green.

### T-future — Bandwidth-depolarization in RM-CLEAN

- Explicitly deferred (decision 3/9). Revisit the offset-table
  optimization once `w_k` becomes genuinely `RM_in`-dependent per-channel
  — likely needs either a full `(φ_peak, Δ)` 2D table or, more promisingly,
  a separable moment-expansion exploiting `bw_depol_correct.f`'s own
  Taylor-expansion structure (order `tn` in powers of `λ²ₖ`) — unworked,
  flagged as a direction not a result. Ground this against the user's
  own thesis derivation, chapter 2.5.1 (see Context above), before
  designing.

### T5 — GPU/CPU hybrid acceleration for RM-CLEAN (exploratory, not started)

**Status: design capture only, nothing scheduled.** Preserved here per
the user's own explicit request ("write this up in a plan doc... we
will implement these in a very near-future project") — not a ticket to
start, and merged into this file rather than kept as a separate
document once the user pointed out it's a continuation of this same
integration effort, not a distinct initiative.

**Motivation, measured not assumed.** A moderately-large-scale
verification run (`cfg/rmclean-jennifer.e2e.cfg`, 4501x4501x101, ~20.26M
pixels) measured `tile_compute` at 97.7-98.2% of `rmclean_cubes`'s own
total wall time across two runs (13023s of 13328s total stage time; the
earlier pre-`lsq_ref_compute=mid` run showed 16739s of 17046s).
`tile_read`/`tile_write` combined are under 2.3% in both — unambiguously
compute-bound, not I/O-bound; SSD/NVMe speed has no bearing on it.
Traced the actual per-iteration cost directly: `refine_peak_matched_
filter`'s fast path calls `rmsf_point_direct` (`src/rmclean.f90:491-512`)
up to 3 times per CLEAN iteration for its own closed-form amplitude fit,
and up to `m_search` (tens to hundreds) more times on escalation. Each
call is a direct `O(nchan)` sum of `cos`/`sin` over every good channel
(~286 for the Jennifer band) — transcendental-heavy, not simple
arithmetic. With `niter` up to 500 and ~20M independent pixels, this is
both embarrassingly parallel across pixels and individually
transcendental-heavy within one — a shape GPUs are architecturally
suited to (raw parallelism plus dedicated special-function-unit
throughput for `sin`/`cos`), *if* the already-known risk (warp/wavefront
divergence from variable per-pixel iteration/escalation counts, flagged
and deferred at decision 6 above) can be managed. The dev machine has a
real, testable GPU (`nvidia-smi`: NVIDIA GeForce RTX 3050), and this
project already has a working OMP `target teams distribute parallel do`
GPU-offload build path (`rm_synthesis_mod.f90`, `GPU_NVFLAGS`/
`GPU_GNUFLAGS` in the `Makefile`) — reusing that same offload convention
here, rather than a new GPU programming model, is the natural low-risk
path if/when this is implemented.

**Explicit design goal: generic, self-calibrating, minimal supervision.**
A hybrid CPU/GPU dispatch mechanism tuned to one dataset's own S/N
distribution is not the goal — this tool runs on data from different
telescopes, bands, and noise regimes, with no expectation a human
re-tunes it per run. Direct analogy to FFTW's own planner
(`FFTW_MEASURE`/`FFTW_PATIENT`: benchmark real candidate strategies on
the actual machine and problem at hand, then commit to the measured
winner, rather than hardcode a strategy): calibrate from a small, cheap,
real sample of the actual data at hand, every run, not from a fixed
rule. This project already has exactly this pattern for a different
parameter — `threshold_snr` derives the CLEAN stopping flux from "500
random sightlines," not a hardcoded value — the mechanism below reuses
that already-validated pattern rather than inventing a new one.

**Idea 1 — per-sightline Faraday-complexity metrics from rm_synthesis.**
Raw peak S/N (already output as `PEAK.MAP.FITS`/`SNR.MAP.FITS`) is not a
tight proxy for CLEAN cost: a single bright, isolated, unresolved point
source can converge as fast as a faint one (iteration count is driven
more by model complexity than raw amplitude), while a moderate-S/N
sightline with two close, blended components can escalate on nearly
every iteration and cost far more. Proposed richer metrics, motivated by
the same qualitative questions a human would ask when visually judging a
dirty Faraday spectrum's own complexity: (1) **peak count** — how many
local maxima exceed a noise-relative threshold (one vs. more than one is
the coarsest, cheapest signal); (2) **peak width ratio** — is the
dominant peak's own width close to the theoretical RMSF FWHM
(`compute_rmsf_fwhm_multiband`, already computed once per run), i.e.
unresolved/simple, or genuinely wider (resolved/extended, intrinsically
harder to CLEAN); (3) **peak clustering span** — are significant peaks
within one contiguous Faraday-depth range, or multiple well-separated
clusters (the same diagnostic a person uses to judge "one blended
feature or genuinely separate RM components far apart"). All three are
computable directly from the dirty spectrum `rm_synthesis` already
builds per pixel, at negligible incremental cost alongside the existing
`PEAK.MAP`/`SNR.MAP`/`RM_PEAK.MAP`/`ANG_PEAK.MAP` outputs — candidate new
maps `NPEAKS.MAP.FITS`/`PEAKWIDTH.MAP.FITS`/`PEAKSPAN.MAP.FITS`, read by
`rmclean_cubes` the same way it already reads `MASK.CUBE.FITS`.

**Idea 2 — self-calibrating CPU/GPU dispatch.** A calibration pass, run
once per `rmclean_cubes` invocation before the main tile loop: (1) sample
N random pixels (matching `threshold_snr`'s own convention) from real,
unmasked sightlines; (2) run each through the existing, unchanged CPU
`clean_complex` path, recording its complexity metric(s) (Idea 1, or S/N
as a fallback) alongside its actual `n_iter_used` and escalation count
(both cheap — `n_iter_used` is already a `clean_complex` output); (3)
bucket the sample by complexity metric and measure the *within-bucket
variance* of iteration/escalation count — this, not correlation with the
metric itself, is the real criterion, since divergence only costs GPU
throughput when a dispatched batch has high internal variance,
regardless of whether the metric predicts absolute cost well; (4) if a
stratification exists where some buckets are meaningfully more
homogeneous than others, assign homogeneous buckets to GPU and
heterogeneous ones to CPU, sized against *measured* relative throughput,
not pixel count alone; (5) if no useful stratification is found (flat
complexity, no GPU available, or too small a sample to distinguish
buckets confidently), fall back to CPU-only automatically — no human
decides whether hybrid dispatch is worthwhile for a given dataset, the
calibration pass decides, every run, from that run's own real data.
Architecturally: preserve the existing contiguous spatial-tile I/O
structure (`tile_ra`/`tile_dec`) rather than gathering/scattering pixels
by classification across the whole image — within each tile, split into
a GPU sub-batch and CPU sub-batch by the calibrated criterion, and
dispatch both concurrently (CPU threads work the tile's own hard pixels
while the GPU works its easy ones at the same time), no rewrite of the
existing tile-based buffer/I/O machinery.

**Idea 3 — RM-faceting (genuinely unexplored territory).** The user's
own excitement here is deliberate: if this works, it's a real
algorithmic contribution to RM-CLEAN, not merely an engineering speedup.
Analogy from wide-field radio synthesis imaging: CLEAN is made tractable
(and parallel) by faceting the sky plane into smaller regions and
deconvolving each mostly independently, but for high-dynamic-range
imaging this isn't fully independent — sources *outside* a facet's own
nominal footprint still contaminate it through the dirty beam's own
sidelobes extending past the facet boundary, so accurate deconvolution
needs *joint* correction across facets. Proposed analogy: facet the
**Faraday depth axis** (not the sky axis) for a single sightline — a
facet is a region of Faraday-depth space where real signal exists,
identified via Idea 1's own peak-clustering diagnostic; a sightline with
multiple well-separated peaks could have its own RM range split into
facets, each CLEANed more independently/in parallel, a natural
additional axis of parallelism and plausibly a very GPU-friendly shape
for the case Idea 1 flags as "simple, isolated."

**Joint correction is full-range, not nearest-neighbor — corrected after
review.** An earlier draft argued the 1-D topology of Faraday-depth
facets (each with only two immediate neighbors, unlike an 8-connected
2-D sky grid) made the joint-correction bookkeeping simpler than the
imaging case. That conflated topological adjacency with the actual
physical reach of the effect, and is wrong: `rmsf_point_direct`'s own
R(delta) is a sum of `cos`/`sin` terms across every channel — structurally
a sinc-like kernel (Brentjens & de Bruyn's own RMSF shape), and sinc-like
kernels decay as roughly `1/delta`, not exponentially. A power-law decay
means a bright enough component's own sidelobe stays non-negligible far
from its own facet, exactly the same character imaging PSF sidelobes
have — precisely why wide-field imaging cannot get away with correcting
only adjacent facets either. The correct picture: a facet's own found
components must be predicted and subtracted against the *whole*
Faraday-depth range, not just its immediate neighbors — the same
major-cycle structure wide-field imaging already uses (predict the full
current model, subtract from the full data, re-facet, repeat), not a
cheaper local patch-up. What genuinely does still favour RM-CLEAN over
the imaging case: the RMSF used here (today, pre-bandwidth-depolarization
— see caveat below) is exact and closed-form, a direct function of the
channels' own `l_sq` values, not an empirically-characterized,
direction-dependent PSF the way imaging often has to deal with — so
whatever the joint correction has to do, full-range or not, can be
computed precisely rather than approximated. That is the real structural
advantage over the imaging analogy; the 1-D neighbor count is not.

**Caveat: this section assumes a shift-invariant RMSF, which stops being
true once bandwidth-depolarization lands (T-future above, decision 3/9).**
Without BW-depol, per-channel weights `w_k` are RM-independent, so
`R(φ;φ₀) = (1/K)Σ w_k·exp(−2i(φ−φ₀)λ²ₖ)` is, as an algebraic identity, a
pure function of the offset `Δ=φ−φ₀` alone (shift-invariant, for any
`λ²` sampling — exactly why `rmsf_point_direct` takes only `delta`, never
an absolute Faraday depth, as its shape argument). With BW-depol on,
`w_k` becomes genuinely `φ₀`-dependent per channel, not just a scalar
factor, and shift-invariance breaks for real. T-future's own point 9
already anticipated this for a *different* optimization (the pre-BW-depol
offset-table inside `compute_dirty_rmbeam`) and deliberately confined it
behind `compute_dirty_rmbeam`'s own general `(φ_grid, RM_in, phase_in) →
beam array` interface specifically so it could be swapped out later
without touching `rm_clean`'s own code, rather than hard-coding
shift-invariance into the CLEAN loop itself. RM-faceting's own
joint-correction step must follow the same discipline: today it can use
a single, fixed R(Δ) kernel valid everywhere in Faraday-depth space; once
BW-depol lands, the correction kernel would need to become genuinely
RM-dependent (per T-future's own note, likely a full `(φ_peak, Δ)` 2D
table or a separable moment-expansion of `bw_depol_correct.f`'s own
Taylor structure — unworked, a direction not a result). BW-depol is also
expected to make the per-evaluation RMSF cost itself more expensive than
today's simple `O(nchan)` sum, which would shift the compute profile
Idea 2's own calibration pass is measuring — today's profiling is a
snapshot of the pre-BW-depol state, not necessarily representative once
BW-depol lands. Design implication: keep RM-faceting's own
joint-correction logic behind a similarly swappable beam-evaluation
interface, not hard-coded against shift-invariance, so it does not need
a full redesign when BW-depol arrives.

This connects back to Idea 1 directly: sightlines classified as "one
contiguous cluster" have nothing to facet (no benefit, added complexity
for free); sightlines with multiple well-separated clusters are the
natural candidates where RM-faceting — if it works — could pay off.
Explicitly flagged as unexplored: not attempted in this project, and not
known to the user to have been explored in the RM-CLEAN literature
generally. If it works, this is a genuine algorithmic contribution worth
its own validation and write-up, independent of whatever it happens to
also do for GPU dispatch.

**How the three ideas relate.** Idea 1 (complexity metrics) is the
measurement layer — it doesn't by itself change how CLEAN runs, only
what's known about each sightline going in. Idea 2 (self-calibrating
dispatch) is the near-term payoff — usable as soon as Idea 1's metrics
(or even just existing S/N) exist, no algorithmic change to CLEAN
itself, only where each pixel's own unmodified CLEAN loop executes. Idea
3 (RM-faceting) is the longer-term, higher-risk research direction — a
genuine algorithmic change to CLEAN, motivated by and targeted using
Idea 1's own clustering diagnostic, needing real numerical validation
(does joint cross-facet correction actually recover the same answer as
unfaceted CLEAN?) before it could be trusted on real data.

**Explicitly out of scope for now:** no implementation of any of the
above in this session or the immediately following one; no GPU kernel
code, no new `rm_synthesis` output maps, no dispatch logic; no real
`n_iter_used`-vs-complexity-metric data pulled from an actual run
(deferred at the user's own explicit request, to avoid anchoring the
generic design on one dataset's own characteristics before the mechanism
itself is designed to be dataset-agnostic).

### T6 -- validation.f90: a flexible cube-slice statistics tool (not started)

**Status:** planning capture only, per the user's own explicit request --
"I will start that work in a new session." No code written.

**Motivation.** Verifying the Jennifer v2 run's correctness (the
`comp_rm_refined`/derotation fix above) relied on ad hoc, one-off Python
(astropy + numpy, reading FITS slices directly) to measure things like
"RMS of AMP in the RM tail (-500 to -360 rad/m^2) across a handful of
random pixels" -- real, useful checks, but not reusable, not in the
codebase, and not something a future run's own correctness can be
quickly re-verified against without re-deriving the same ad hoc script.
The user wants this promoted to a real, permanent tool.

**Concrete motivating check (the one that triggered this ticket):**
tail-region (RM far from any real emission) dirty/restored/resid AMP
RMS for 8 random Jennifer v2 sightlines came out ~13-30 uJy/beam
(mean ~22-26 uJy/beam across the three cubes), independently confirmed
by the user to be close to the band-averaged Q/U noise level -- a good
end-to-end sanity check. The user's own follow-up, not yet checked:
**the residual should ALSO drop to this same ~26 uJy/beam noise floor
in the RM channels AT/NEAR the peak emission**, not just far away in
the tail -- if CLEAN successfully removed the real signal, the
residual near the source should look like pure noise too, indistin-
guishable from the tail. This is a genuine "goodness of CLEAN"
diagnostic (did CLEAN converge, or is the threshold/niter/gain
combination leaving real, structured signal behind near the peak) that
today has no automated check anywhere in this project.

**What the tool should do (scope, as discussed, not yet designed in
detail):**
- Extract a slice of any RM cube (AMP/PHA, or the underlying RE/IM
  pair) in flexible ways: an RM-index or RM-value range (e.g. "the
  15 bins nearest a given RM", or "the tail beyond +/-X rad/m^2 from
  the sightline's own peak", or an absolute RM range), a spatial
  region (a box, or a list of specific (x,y) pixels, or N random
  pixels satisfying some validity criterion like `NVALID.MAP > threshold`),
  or some combination (e.g. "the residual cube, in a window centred on
  each sightline's own peak RM, for N random valid pixels").
- Report standard statistics over that slice: RMS, mean, std, min/max --
  matching what the ad hoc Python already computed this session, not a
  novel statistics design.
- Compare-across-cubes mode: given the same slice definition, report
  the same statistic for two or more cubes side by side (dirty vs
  restored vs resid, the exact pattern used this session) -- this is
  what makes "residual near the peak should match the tail's own noise
  floor" a single, repeatable command instead of a bespoke script.
- Likely CLI-driven like every other tool in this project (`bin/validate_cube`
  or similar; the user's own suggested source file name is
  `validation.f90`, naming/scope to be finalized when this work actually
  starts) -- matching `rmclean_cubes`/`rm_synthesis`'s own cfg-or-CLI
  convention rather than inventing a new interface style.

**Explicitly not designed yet:** exact CLI/cfg schema, exact slice
syntax, exact output format (human-readable table vs CSV vs a value
usable in `tests/run_tests.sh` assertions), whether this becomes a true
regression-test building block (e.g. "assert residual-near-peak RMS is
within Nx of tail RMS" as an automated pass/fail, not just a printed
number) -- all open questions for the session that actually picks this
up.

### T7 -- CLEAN convergence/stop-reason logging + subimage gain-tuning workflow (done)

**Objective.** Two things the user asked for in the same session, done
together since the second depends on the first: (1) make it possible to
answer "is CLEAN converging, and why did it stop" from the log alone,
for both a single traced sightline and the whole run; (2) get real,
data-driven evidence for choosing `niter`/`gain`/`threshold_snr` on the
actual Jennifer dataset, without paying for a full 4501x4501 run per
trial.

**Scope/what was built.**
- `clean_complex` (`src/rmclean.f90`) gained a new required output,
  `stopped_by_threshold` (alongside the pre-existing `n_iter_used`):
  disambiguates "threshold criterion fired" from "niter cap exhausted
  without ever satisfying it" -- `n_iter_used==niter` alone cannot tell
  these apart (the threshold criterion can also fire, coincidentally,
  on the literal last allowed iteration). Updated all 5 call sites
  (`rmclean_cubes.f90` x2, plus 3 direct test callers).
- Three new OPTIONAL per-iteration trace outputs on the same
  subroutine -- `trace_peak_val`/`trace_rms_val`/`trace_flux_val`
  (cumulative `sum(|comp|)` so far) -- filled only when the caller
  supplies the arrays, so the hot (non-traced) path pays nothing extra.
  `rmclean_mod` stays "pure computation, no logging" (its own top
  comment's stated design boundary): it only returns more data
  optionally; `rmclean_cubes.f90`'s `clean_one_pixel` decides which
  pixel (`trace_ix`/`trace_iy` cfg keys, 1-indexed global, default 0 =
  off) to request it for and does the actual `log_message` calls.
  `log_every` (default 50) throttles those lines to every Nth
  iteration; iteration 1 (the pre-CLEAN state -- `trace_flux_val(1)`
  is exactly 0 there) and the final iteration are always logged
  regardless, so the trend's start/end are never missed.
- Run-wide aggregate summary (always on, no cfg key -- cheap
  population-level counters, not per-iteration): printed once at run
  end alongside the pre-existing `CLEANed N pixels.` line -- count
  stopped-by-threshold vs hit-niter-cap (with %), `n_iter_used`
  mean/min/max, total cleaned flux summed over every pixel, and
  peak-residual mean/min/max across pixels. Answers "why did CLEAN
  stop" and "how much flux/residual" at the whole-run level without
  per-pixel log spam.
- Block-progress logging (a smaller, earlier ask in the same session):
  `plan_rmclean_tile` now prints the total block count once
  (`-- N block(s) total.`), and the main tile loop prints `Block m of n
  processed.` after each block's CLEAN compute finishes -- previously
  neither was logged at all, so progress was unreadable without
  manually deriving block count from the tile-plan geometry.
- Subimage-based fast-iteration workflow: no new cutout tool needed --
  reused `rm_synthesis`'s own pre-existing `subim=y` CFITSIO-subsection
  read directly against the already-convolved Jennifer Q/U cubes
  (`cfg/rmsynth-jennifer.subim128.cfg`, 128x128 px centred on the
  4501x4501 image, RA/Dec 2187-2314, full 286-channel band, otherwise
  identical RM settings to the production e2e cfg). Runs in ~7s
  (vs. hours for the full cube); confirmed the centre-of-image crop
  genuinely contains real signal (SNR up to ~20, not just noise), so
  it exercises both the tail-noise-floor check and the "residual near a
  real peak" check in one small cube. Working files under
  `scratch/rmclean_tuning_subim128/`.

**Correctness gate.** New standalone test `tests/test_rmclean_stop_reason.f90`
(registered as `run_tests.sh` section 37): a noiseless single-point-
source scenario, Case A (generous niter/gain/threshold) asserts
`stopped_by_threshold=.true.`, `n_iter_used<niter`, a real (>5) initial
trace peak, a weakly-monotonically-DECREASING peak trace, and a weakly-
monotonically-INCREASING flux trace starting at exactly 0; Case B
(niter=3, deliberately too small) asserts `stopped_by_threshold=.false.`
and `n_iter_used==niter` exactly (cap exhausted, not a coincidental
threshold hit on iteration 3). Full regression suite re-run after the
`clean_complex` signature change (affects every caller): 118/118 pass
(109 pre-existing + 7 new), including the existing small-tile/
io_read_threads/io_write_threads/io_overlap bit-identical checks --
confirms the new outputs/logging don't change any existing numerical
result. Also smoke-tested end-to-end against real subimage data
(`trace_ix=112 trace_iy=89`, the subimage's own SNR-peak pixel):
iteration-1 trace peak matched that pixel's own dirty PEAK.MAP value,
`log_every=50` correctly emitted iter=1/50/final(82) and skipped the
rest, and the run-end summary's stop-reason counts matched expectations.

**Evidence -- gain sweep on the 128x128 subimage** (niter=500,
threshold_snr=5.0, `lsq_ref_compute_mode=mid` -- identical to
production except `gain`; OMP_NUM_THREADS=6, OMP_PROC_BIND=close,
OMP_PLACES=cores per this host's own standing thread-count convention):

| gain | stopped@threshold | hit niter cap | n_iter_used mean | tail RMS (uJy, median) | peak-window RMS (uJy, median) | peak/tail ratio |
|------|-------------------:|--------------:|-----------------:|-----------------------:|-------------------------------:|-----------------:|
| 0.1  | 96.81% | 3.19% | 233.7 | 6.89 | 7.31 | 1.061 |
| 0.2  | 99.88% | 0.12% | 132.9 | 5.65 | 5.93 | 1.049 |
| 0.3  | 99.88% | 0.12% |  98.2 | 4.59 | 4.73 | 1.031 |
| 0.5  | 99.87% | 0.13% |  74.4 | 2.86 | 2.85 | 0.995 |
| 0.7  | 98.07% | 1.93% |  89.7 | 1.67 | 1.44 | 0.862 |

Tail RMS on the matching DIRTY cube (no CLEAN at all): median 25.46
uJy/beam, mean 25.62 uJy/beam -- consistent with the ~22-26 uJy/beam
band-averaged Q/U noise floor already established on the full
`jennifer_e2e` cube (T6's own motivating check), confirming this
subimage's noise properties are representative, not an artefact of the
crop.

**Finding (real, not just a tuning preference).** Every gain tested
leaves a CLEANed tail RMS well BELOW the true instrumental noise floor
(6.89 uJy at gain=0.1, i.e. already ~3.7x quieter than the 25.46 uJy
dirty floor) -- because `threshold_snr`'s stopping criterion is
self-referential (`thresh x rms_val`, `rms_val` recomputed from the
CURRENT residual every iteration over the WHOLE spectrum including the
tail), so a more aggressive gain shrinks its own stopping target as it
goes, converging to a genuinely lower absolute residual rather than
just reaching the same target faster. This gets monotonically WORSE as
gain increases, and at `gain=0.7` the peak/tail ratio drops BELOW 1.0
(0.862) -- the region right at a real source ends up QUIETER than blank
sky, a clear over-subtraction signature (removing genuine noise/
structure as if it were flux), not a healthy convergence result.
`gain=0.1` (current production default, used in the v1/v2/v3 Jennifer
runs) has the peak/tail ratio closest to the theoretically-expected
value (slightly ABOVE 1 -- real structure at the peak leaving marginally
more residual than pure tail noise, the correct qualitative sign) of
all 5 gains tested, but pays for it with a 3.19% niter=500 non-
convergence rate (vs. ~0.12% at gain=0.2/0.3/0.5). **Conclusion: `gain=
0.1` is not just a conservative default, it is the best of the 5 tested
for residual noise fidelity -- the ~3% non-convergent pixels should be
fixed by RAISING `niter` (e.g. 800-1000), not by raising `gain`, since
raising gain demonstrably trades residual noise fidelity away.** Whether
niter=500 itself, or the `threshold_snr=5.0` self-referential stopping
criterion's own design, should change is an open question for a future
session -- not resolved here.

**UPDATE (T8, same day):** the self-referential stopping criterion
itself -- not just niter/gain tuning around it -- turned out to be a
real bug, not just a design choice worth revisiting. See T8 below: fixing
it (per-pixel, ONCE-computed RM-tail sigma, `auto_nsigma=`) took this same
subimage's tail RMS from 6.89 uJy (gain=0.1, old `threshold_snr` design --
3.7x quieter than the true ~25.5 uJy noise floor) to 25.24 uJy (matching
the dirty cube's own floor almost exactly), AND converged in a mean of
3.89 iterations instead of 233.7. The gain-fidelity ranking above (0.1
best of 5) was a real, correctly-measured finding about the OLD
mechanism, superseded by T8's fix rather than contradicted by it.

### T8 -- CLEAN stopping-criteria redesign: abs_flux_floor/auto_nsigma/niter, per-pixel tail sigma, exact stop_reason logging (done)

**Objective.** The user asked directly: "why then do we have n-sigma?
... clean up this mess ... we want unambiguous variable names ... stop
based on whichever condition is met first: niter, absolute value
(Jy/mJy/uJy), auto threshold (user multiplier x sigma derived from the
data itself, RM tail etc) ... logs should report exact cause."

**Root cause found.** `clean_complex`'s old `thresh` parameter was fed
by BOTH `rmclean_cubes.f90`'s `threshold=<absolute Jy value>` and
`threshold_snr=<n-sigma multiplier>`, and was UNCONDITIONALLY multiplied
by `rms_val` every iteration (`peak_val-avg_abs <= thresh*rms_val`). So
`threshold=` never actually compared an absolute flux value at all --
it silently became another sigma-multiplier, comparing flux^2 against
flux. And `threshold_snr` computed an absolute number ONCE up front
(`threshold_snr x noise_floor`, from a whole-cube pre-scan of
`noise_nlos` random sightlines), then THAT got multiplied by `rms_val`
again -- `rms_val` recomputed fresh from the CURRENT (shrinking)
residual every single iteration, over the WHOLE spectrum. A
self-referential, ever-tightening target, not a fixed noise floor --
directly explaining T7's own gain-sweep finding (residual RMS well
below the true noise floor, getting worse with higher gain).

**Design (three independent, unambiguously-named criteria, whichever
fires first per pixel, niter always the hard backstop):**
1. `niter` -- unchanged, the loop's own iteration cap.
2. `abs_flux_floor=<v>` -- stop the INSTANT a pixel's own peak amplitude
   drops to/below this literal fixed value (native units, or Jy/mJy/uJy
   suffix). No noise/baseline adjustment -- a pure "brightest remaining
   feature below X" comparison. This is what `threshold=` was always
   supposed to mean and never actually did.
3. `auto_nsigma=<n>` -- stop when `(peak_val-avg_abs) <=
   auto_nsigma x tail_sigma`, where `tail_sigma` is estimated ONCE (new
   `estimate_tail_sigma`, `src/rmclean.f90`) from THIS pixel's own dirty
   spectrum: find its own peak bin, exclude a window of
   `tail_exclude_nfwhm x fwhm_rm` (default 3.0 FWHM) around it, and take
   `rms_about_mean` of whatever RM bins remain (the pixel's own "RM
   tail" -- far from its own peak, presumed noise-dominated regardless
   of whether a real source is present, same reasoning as T6's own
   motivating check). NOT recomputed per iteration -- this is the fix
   for the moving-target bug above. NOT a whole-cube pre-scan either
   (unlike the old `threshold_snr`/`estimate_noise_floor`/`noise_nlos`
   machinery, deleted entirely) -- genuinely per-pixel, matching the
   user's own explicit choice ("literal per-pixel RM-tail sigma", not
   reusing the old whole-cube pre-scan) when asked to disambiguate this
   design fork.

Unlike the old `threshold=`/`threshold_snr=` pair (mutually exclusive,
exactly one required), `abs_flux_floor=`/`auto_nsigma=` are each
independently optional and MAY be combined -- whichever fires first
wins for that pixel. Neither given is also valid (niter-only, same
precedent as the old `thresh=0.0`), with a printed NOTE so it's not
silently missed.

**Logging (the user's 4th, non-negotiable requirement).** `clean_complex`
now returns `stop_reason` (`'niter'`/`'abs_flux'`/`'auto_nsigma'`, not
just a `stopped_by_threshold` boolean, which could not distinguish 3
outcomes) -- disambiguates `n_iter_used==niter` (cap exhausted) from a
criterion coincidentally firing on the literal last allowed iteration.
`rmclean_cubes.f90` tallies all 3 reasons run-wide (extends the T7
aggregate summary: `hit niter cap` / `stopped at abs_flux` / `stopped at
auto_nsigma`, each with a %) and logs the exact `stop_reason` per
iteration/pixel for `trace_ix=`/`trace_iy=`-traced pixels (T7's own
per-iteration trace mechanism, unchanged otherwise).

**Naming discipline (the user's other explicit request).** Every new
parameter name states exactly what it is -- `abs_flux_floor` (a flux
value), `auto_nsigma_mult` (a multiplier), `tail_sigma`/`tail_exclude_
nfwhm` (a sigma estimate and its own exclusion window) -- no quantity is
shared between two different roles. `auto_nsigma_mult` is deliberately
NOT called anything containing bare "nsigma" alone, to avoid collision
with the pre-existing, unrelated `nsigma_refine` (the tiered
peak-refinement escalation threshold in `refine_peak_matched_filter`) --
two genuinely different "n-sigma" concepts in this codebase, now
impossible to confuse by name alone. See
[[feedback_unambiguous_naming]] for the standing memory this created.

**Correctness gate.** `clean_complex`'s public signature changed (new
required `fwhm_rm`/`have_abs_flux_floor`/`abs_flux_floor`/
`have_auto_nsigma`/`auto_nsigma_mult`/`tail_exclude_nfwhm` args, `thresh`
removed, `stopped_by_threshold`->`stop_reason`) -- updated all 5 call
sites (`rmclean_cubes.f90` x2 tables x2 trace-variants = 4, plus 3 direct
test callers: `test_rmclean_lsqref_flex.f90`/`test_drm_floor.f90`
translated to niter-only mode, exact behavioral match for their old
near-zero/n-sigma thresh values in their noiseless synthetic scenarios;
`thesis_scenario_rmclean.f90` translated to niter-only, an EXACT match
for its old `thresh=0.0`). `tests/test_rmclean_stop_reason.f90` rewritten
for 3 cases: Case A (niter-only) asserts `stop_reason=='niter'` and
`n_iter_used==niter` exactly with no early-exit possible by
construction; Case B (`abs_flux_floor`, noiseless) asserts
`stop_reason=='abs_flux'`, early convergence, and the final residual
peak at/near the floor (plus the pre-existing trace-array checks); Case
C (`auto_nsigma`, WITH fixed-seed injected Box-Muller noise -- a
noiseless spectrum's own tail sigma is ~0, which would only test
machine-precision convergence, not the criterion itself) asserts
`stop_reason=='auto_nsigma'` and early convergence. Full regression
suite: 121/121 pass (all pre-existing plus the rewritten section 31 --
`abs_flux_floor`/`auto_nsigma` unit conversion, determinism across two
runs with no seed at all now needed since per-pixel tail-sigma has zero
randomness, valid-combined and valid-neither-given cases replacing the
old mutually-exclusive/required-one checks -- and section 37's own 3
new cases). `cfg/rmclean-e2e-smalltest.cfg`, `cfg/rmclean-jennifer.e2e.cfg`,
`cfg/rmclean-example.cfg` updated for the renamed keys.

**Evidence -- the fix actually fixes T7's finding, not just the logging.**
Re-ran the 128x128 subimage (same data as T7) with `auto_nsigma=5.0`
(matching the OLD `threshold_snr=5.0` production value) under the NEW
design: tail RMS **25.24 uJy (median) / 25.42 uJy (mean)** -- matching
the dirty cube's own true noise floor (25.46/25.62 uJy) almost exactly,
versus the OLD design's 6.89/7.22 uJy (3.7x too quiet) at the same
nominal "5-sigma" setting. Mean iterations to convergence: **3.89**
(min 1, max 19) versus the old design's 233.7 -- the self-referential
moving target was not just numerically wrong, it was also needlessly
slow. Stop-reason summary on this run: 100.00% `auto_nsigma`, 0%
niter/abs_flux, confirming every pixel now stops for the intended,
physically-meaningful reason.

**UPDATE (T9, same day):** the `tail_exclude_nfwhm`-based per-pixel
sigma estimate above was itself replaced -- see T9. It silently fails
for multi-component or extended/Faraday-thick sightlines (a second peak,
or a broad plateau wider than the excluded window, still contaminates
the "tail"); T9's percentile+IQR method (further corrected for the
Rayleigh, not Gaussian, distribution of amplitude data) replaces it.
Separately, T9 also found and fixed a second bug in the `auto_nsigma`
criterion ITSELF (the `- avg_abs` subtraction, inherited unexamined from
this ticket's own carryover of the pre-existing algorithm) that this
ticket's own Evidence section did not catch, because the T7 subimage's
own gain-sweep tuning happened not to exercise a large-scale coherent
structure the way T9's own high-signal-RM-plane test did.

### T9 -- auto_nsigma correctness: percentile+IQR sigma (Rayleigh-corrected), avg_abs removed (done)

**Objective.** T8 shipped `auto_nsigma`'s per-pixel sigma as a peak-
relative exclusion window (`tail_exclude_nfwhm`). The user immediately
flagged the real weakness: "What if there are extended rm features and
multiple peaks?" -- a second peak, or a broad plateau wider than the
excluded window, silently contaminates the "tail" that window design
relies on. Separately, the user asked for a concrete image-plane test
("rms in an image plane... at a high signal RM plane") to see whether
the residual sigma there had actually decreased to noise level -- that
test surfaced a SECOND, independent bug in the criterion itself.

**Fix 1 -- percentile+IQR replaces the exclusion window.** New
`estimate_percentile_sigma` (`src/rmclean.f90`, replacing
`estimate_tail_sigma`): sorts THIS pixel's own dirty AMPLITUDE spectrum,
takes the lowest `noise_percentile` fraction (default 0.20) of bins --
an order-statistic selection that doesn't care where or how many real
features exist, only that a real source is a narrow-in-amplitude-rank
excursion (same reasoning the old, retired, whole-cube `noise_percentile`
pre-scan used, applied per-pixel here instead). `fwhm_rm`/
`tail_exclude_nfwhm` are gone from `clean_complex`'s own signature
entirely (no longer needed -- one fewer required argument, one fewer
concept). cfg key renamed `tail_exclude_nfwhm=` -> `noise_percentile=`
throughout (`rmclean_cubes.f90`, `cfg/rmclean-example.cfg`,
`cfg/rmclean-jennifer.e2e.cfg`, 3 direct test callers, `tests/
test_rmclean_stop_reason.f90`'s own Case C).

**The image-plane test and what it found.** Ran the T7 subimage's
CLEAN-peak pixel through the RMSF FWHM-corrected design, then looked at
the FULL 128x128 image plane at that pixel's own RM_PEAK channel (not
just that one pixel's own RM-axis spectrum). Two real findings:

1. **This subimage's "high signal" isn't a point source.** 3801 of 4317
   high-SNR (>7) pixels (88%) share nearly the same RM_PEAK (~150
   rad/m^2) -- a large-scale coherent/extended polarised structure
   spanning most of the field, not an isolated source. (The
   spatially-widespread sharing of one RM_PEAK does NOT imply each
   individual pixel's own RM-axis spectrum is broad -- that would be a
   real-signal/beam-sidelobe distinction specific to each pixel, not a
   field-wide property; conflating the two was an error corrected
   mid-session, see point 3 below.)
2. Image-plane RMS at that channel: dirty 108.6 uJy, resid only dropped
   to 41.9-87.6 uJy across iterations of this fix (still well above the
   ~24-26 uJy true tail floor at that point) -- prompting a direct
   per-iteration trace of the brightest pixel to find out why.

**Fix 2 -- the REAL bug: `- avg_abs` in the auto_nsigma comparison.**
Tracing the brightest pixel showed CLEAN stopping at `peak_val=37 uJy`
after only 23 (of 500) iterations -- not niter exhaustion (0% of all
16384 pixels hit the cap), and not because `percentile_sigma` itself
was reached. The actual comparison, inherited unexamined from the
pre-existing algorithm, was `(peak_val - avg_abs) <= auto_nsigma_mult *
percentile_sigma`, where `avg_abs` (mean |residual| over the WHOLE
spectrum) is recomputed from the CURRENT residual every iteration --
the exact same self-referential moving-target flaw already fixed once
for the multiplicative side (`rms_val` -> `percentile_sigma`, fixed at
T8), still present on the subtractive side. For a large-scale coherent
structure, peak~average almost everywhere, so `peak-avg_abs` collapses
toward zero long before the peak has actually reached anything close to
the noise floor -- stopping CLEAN early even though the residual peak
was still far above it. Fixed by dropping the subtraction entirely:
`auto_nsigma` now compares `peak_val <= auto_nsigma_mult *
percentile_sigma` directly, matching `abs_flux_floor`'s own (already
subtraction-free) form. Re-traced the same pixel: iterations to
convergence went from 23 -> 63, peak_val at stop from 37 -> 15.5 uJy --
converging genuinely deeper, not just differently.

**Fix 3 -- Rayleigh, not Gaussian, IQR-to-sigma conversion.** Even after
Fix 2, the traced pixel's own `percentile_sigma` (2.84 uJy) was ~8x
below its own independently-measured true tail sigma (22.4 uJy). Root
cause: dirty AMPLITUDE (`sqrt(re^2+im^2)` of iid-Gaussian re/im) is
Rayleigh-distributed, not Gaussian, AND the quartiles are taken WITHIN
the bottom `noise_percentile` slice, not the whole population -- a
doubly-different shape from a plain Gaussian IQR, for which the
standard `IQR/1.349` conversion is analytically wrong. Derived the
correct, closed-form relation: for a Rayleigh(sigma) variate, the p-th
quantile is `x_p = sigma*sqrt(-2*ln(1-p))`; quartiles taken within the
bottom fraction `q` correspond to OVERALL percentiles `0.25q`/`0.75q`,
giving

```
IQR_q(sigma) = sigma * [ sqrt(-2*ln(1-0.75q)) - sqrt(-2*ln(1-0.25q)) ]
```

so `sigma = IQR_measured / IQR_q_factor(q)`, computed once from
`noise_percentile` alone (no new data dependency). For q=0.20,
`IQR_q_factor ~= 0.2497` (vs the Gaussian 1.349 and the untruncated-
Rayleigh-population 0.906 -- both wrong for this specific quantity).
Implemented in `estimate_percentile_sigma`. Caveat, stated explicitly in
the code comment: this still assumes the bottom-q fraction is itself
noise-dominated; an exceptionally bright source's own RMSF-sidelobe
spread CAN still occupy more than `1-q` of a spectrum (confirmed
directly for the one brightest pixel in this subimage: 91 of 101 bins
above 5% of its own peak) -- a real, inherent limit of any percentile
method, not fixed by the conversion factor alone, and not evidence that
"widespread spatial structures have broad per-pixel spectra" in
general (that conflation was corrected mid-session; ordinary sightlines
occupy only a few RM bins, exactly as expected).

**Evidence.** Full regression suite: 121/121 pass after all three fixes
(same test count as T8 -- no new tests needed beyond updating existing
call sites/values for the signature change). Re-ran the T7/T8 subimage
end to end with the fully-fixed design (`auto_nsigma=5.0,
noise_percentile=0.20`):
- Aggregate tail RMS (true noise-only region): **25.26/25.44 uJy**
  (median/mean) -- matching the dirty cube's own floor (25.46/25.62 uJy)
  more closely than either prior attempt, and for the right reason this
  time (verified via the avg_abs/Rayleigh mechanism, not coincidence).
- At the brightest pixel's own peak (RM=170): dirty 165.7 uJy -> resid
  87.2 uJy. **This is the correct, expected outcome, not under-
  convergence** -- `auto_nsigma=5.0` means "stop once the peak is no
  longer >5 sigma above noise," so a properly-converged residual AT A
  REAL SOURCE's peak should settle near `N*sigma` (~5*22=110 uJy here,
  87 uJy is close, modulo discrete gain=0.1 iteration steps), NOT
  collapse to the bare 1-sigma noise floor -- only genuinely empty sky
  (the tail) should reach that. Recognizing this distinction (tail ->
  1-sigma; real-source peak -> N-sigma) corrects an implicit
  expectation this whole investigation had been carrying since T7 that
  "residual near a real peak should equal the tail" -- true only in the
  weak sense of "not dramatically below the noise floor" (T7's own
  actual finding), not "should reach exactly the tail's own value."

**UPDATE (same day):** T9's own percentile-then-IQR design (lowest
`noise_percentile` fraction, default 0.20, truncated-Rayleigh-corrected)
was itself replaced after the user asked directly why restrict to a
subset at all rather than using the full available spectrum. Measured
both methods directly against many real pixels' own independently-known
true noise (same far-RM-tail ground truth as above), split into the two
scenarios that matter:

| scenario | method | ratio to true (median) | ratio to true (mean +/- std) |
|---|---|---:|---:|
| highest-SNR pixel (1 case) | full-spectrum IQR | 1.00 | -- |
| highest-SNR pixel (1 case) | truncated-20% IQR | 0.68 | -- |
| 10 typical high-SNR pixels | full-spectrum IQR | 1.00-1.07 | 1.15-1.16 +/- 0.31-0.32 |
| 10 typical high-SNR pixels | truncated-20% IQR | 0.78-0.79 | 0.77-0.78 +/- 0.20-0.21 |

Full-spectrum wins decisively in BOTH scenarios, not a tradeoff -- the
truncated design's own smaller sample (~20 of 101 bins) combined with
its much larger IQR-to-sigma correction factor (~0.25 for q=0.20 vs
~0.91 for the full population, a ~4x amplification) cost more accuracy
than real per-pixel signal contamination ever saved. That contamination
was ALSO measured directly (not assumed): a first check in fraction-of-
peak units was flawed (it mostly counted ordinary Rayleigh noise
fluctuation for modest-SNR sources, not real contamination); redone in
absolute-sigma units against each pixel's own true noise, real
sidelobe-driven contamination is real but moderate (~5-12% of bins
exceed 3x true sigma for a typical bright pixel, vs ~1 bin expected by
pure chance) -- comfortably under the 25% that would bias a full-
spectrum Q1, EXCEPT for at least one individually-checked pixel where
contamination reached ~30% of bins and full-spectrum IQR overestimated
true sigma by ~1.9x (a real, demonstrated limit of ANY percentile-free
estimator, not fixed by this change).

`estimate_percentile_sigma`/`noise_percentile` retired entirely (not
aliased) -- replaced by `estimate_iqr_sigma` (`src/rmclean.f90`): full
101-bin dirty amplitude spectrum's own IQR, FIXED Rayleigh full-
population factor (`IQR/0.90656`, a constant, not data-dependent).
Simpler code (one fewer cfg key, one fewer function parameter, no
subsetting logic) AND more accurate. Full regression suite: 121/121
pass (updated call sites in `rmclean_cubes.f90`'s 4 `clean_complex`
call sites, 3 direct test callers, `tests/test_rmclean_stop_reason.f90`'s
own 3 cases).

**cfg defaults finalised (same day), each justified per-file, not
copied from the old design:**
- `abs_flux_floor` and `auto_nsigma` both default OFF
  (`have_abs_flux_floor`/`have_auto_nsigma` start `.false.`) -- neither
  criterion silently activates; niter alone governs unless the user
  explicitly opts in to one or both (the user's own explicit
  requirement: "if user does not specify the multiplier... we should
  not trigger this as stopping criterion").
- `auto_nsigma`'s own default VALUE, where used, is `1.0` (not the old
  `5.0`) -- justified directly by the full-spectrum IQR sigma's own
  measured accuracy above (median ratio ~1.00-1.07 vs true noise): a
  1-sigma target is no longer an arbitrary aggressive choice once the
  sigma estimate itself is trustworthy.
- `abs_flux_floor` values are now spelled out with an explicit unit
  suffix (uJy) in every cfg, not a bare native-units number -- makes
  the actual physical value legible in the cfg file itself without
  needing to know that dataset's own BUNIT by heart.
- `cfg/rmclean-e2e-smalltest.cfg`: `abs_flux_floor=10000uJy` (=0.01
  Jy/beam, this fixture's own BUNIT) -- UNCHANGED value, just re-
  expressed with an explicit unit; this is a regression-test
  convenience choice (~1% of the fixture's own ~1 Jy injected-source
  peak, deep enough to exercise this section's own tiling/threading/
  io_overlap/mask-cache bit-identical checks) predating this session,
  not a noise-derived choice -- left as-is since changing it would risk
  invalidating already-verified baselines for no new information.
- `cfg/rmclean-example.cfg` (generic template, unknown target data):
  `abs_flux_floor=0.0uJy` -- a documented but EFFECTIVELY INERT
  reminder (real dirty amplitude is always >0 from noise alone, so
  `peak_val<=0.0` essentially never fires) that the key exists, since
  any nonzero number here would be arbitrary for genuinely unknown
  data; `auto_nsigma=1.0` (uncommented, active) carries the real
  default, since it needs no prior knowledge of the target data.
- `cfg/rmclean-jennifer.e2e.cfg` (real, specific dataset):
  `abs_flux_floor=20uJy` -- a GENUINE, ACTIVE second criterion, not a
  placeholder, anchored to this dataset's own independently-and-
  repeatedly-measured ~25.5 uJy/beam noise floor (T7/T9's own tail-RMS
  checks) with a small margin below it. Deliberately NOT left inert
  here: T9 directly demonstrated `auto_nsigma`'s own per-pixel sigma
  estimate can be biased for individual pixels even with the corrected
  full-spectrum method (the ~1.9x-overestimate pixel above), so a fixed
  floor anchored to data we actually trust is a genuine safety net --
  whichever of the two fires first wins per pixel. Verified on the real
  subimage: a meaningful, non-degenerate split (10.03% of pixels
  stopped via `abs_flux_floor`, 89.97% via `auto_nsigma`), confirming
  both criteria are doing real work, not one dominating.

**T9 UPDATE 2 (same day) -- unrelated pair of bugs found running the
real full Jennifer cube with the curated cfg:**
1. Per-thread `dur_ms` log field (`F10.3`) overflows to `**********`
   once a single tile block's own compute stage runs longer than
   ~16.7 min (only 6 integer digits fit) -- found live, a block
   genuinely took 55 min. Display-only (the timing itself was always
   correct); widened to `F18.3`.
2. `plan_rmclean_tile`'s own `bytes_per_tile_pixel` budget formula
   assumed 8 array-widths (2 input + 6 output, all single-buffered),
   but the 6 output arrays (`clean_re_buf(...,2)` etc.) are allocated
   PERMANENTLY double-buffered for the tile-write-join mechanism,
   regardless of whether `io_overlap` is even on -- 14 array-widths,
   not 8, a 1.75x silent under-budget that predates this session
   (T4d). Found because the real run's own cgroup memory reached
   39.9G against a nominal `mem_frac_ram=0.25` budget of ~15.7G;
   confirmed by recomputing `tile_pixels_max` with the old formula --
   it matched the run's own logged tile size almost exactly. Fixed
   (`2 + 6*2` instead of `8`); the corrected formula accounts for
   ~27.4G of the observed 39.9G, the remainder being CFITSIO/page-
   cache overhead `mem_frac_ram` was never designed to budget for at
   all. Full regression suite: 121/121 pass after both fixes.
