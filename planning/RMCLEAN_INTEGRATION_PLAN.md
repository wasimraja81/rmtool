# RM-CLEAN Integration Plan

Branch: `rmclean-integration` (from `develop`)

**Status: T0 (MASK.CUBE.FITS per-channel frequency/band table) and T1
(RM-CLEAN core algorithm module, `src/rmclean.f90`) done and verified.
T2 (standalone tool + per-pixel mask-pattern caching) scoped from design
discussion but not yet implemented; expect numbering/scope to firm up
once started, per this project's own convention of writing a ticket's
Evidence section only once its own work is done.**

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
   `required_drm_nyquist(l_sq, nchan, lsq_ref, oversample, drm_max)`:
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
     `required_drm_nyquist`-derived grid, the dirty Re/Im panels show
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

### T2 (not yet detailed — standalone tool consuming T1, scoped from
### design discussion above, numbering/scope to firm up once started)

- Per-pixel mask-pattern pre-scan + pattern-keyed offset-table cache
  (decision 10), consuming T0's new per-channel table.
- Standalone tool wiring: CLI/config schema, own binary/build target,
  consuming existing dirty AMP/PHA cubes (decision 11).

### T-future — Bandwidth-depolarization in RM-CLEAN

- Explicitly deferred (decision 3/9). Revisit the offset-table
  optimization once `w_k` becomes genuinely `RM_in`-dependent per-channel
  — likely needs either a full `(φ_peak, Δ)` 2D table or, more promisingly,
  a separable moment-expansion exploiting `bw_depol_correct.f`'s own
  Taylor-expansion structure (order `tn` in powers of `λ²ₖ`) — unworked,
  flagged as a direction not a result. Ground this against the user's
  own thesis derivation, chapter 2.5.1 (see Context above), before
  designing.
