#!/usr/bin/env python3
"""
UR-12/UR-13 (docs/dev/BRANCH_STATUS_update-reproject.md): interpolation
accuracy-vs-speed benchmark for reproject_cubes' interp= schemes, plus
an apples-to-apples comparison against astropy's own reproject package
(reproject_interp/reproject_exact) wherever a genuine equivalent exists.

Ground truth is analytic, never resampled: each test case defines a
source as a closed-form function of CONTINUOUS input-pixel coordinates
(a sum of 2D Gaussians, or a sharp step edge). The actual FITS input is
built by sampling that function onto the input's own integer pixel
grid -- the same discretization any real image has. The ground truth
at a given OUTPUT pixel is the exact value of the SAME continuous
function, evaluated at the exact fractional input-pixel coordinate that
output pixel's own sky position maps to (via astropy.wcs, independent
of any resampling kernel) -- not a value sampled from any resampled
grid. This is what makes the comparison non-circular.

Two different, established operations are being measured, not one
scheme failing where another succeeds -- the distinction geophysical
regridding tools (ESMF/xESMF/CDO) draw between "interpolation" and
"conservative remapping", not a framing specific to this file:

  - INTERPOLATION (every scheme here except interp=exact, which this
    file does not cover -- see check_interp_cost.py): treats the input
    as point samples of a continuous, band-limited signal and
    reconstructs that signal's value AT A POINT, via the
    Shannon-Whittaker sampling theorem a sinc-family kernel
    implements. rms_error below measures exactly this: how close the
    reconstructed value at each output pixel's own coordinate is to
    the true analytic value at that same point.
  - CONSERVATIVE REMAPPING (interp=exact, in check_interp_cost.py):
    treats each input pixel as a piecewise-constant cell and computes
    the exact geometric overlap between input and output cells, which
    by construction exactly preserves the integral (flux) of that
    piecewise-constant field -- not an approximation to get closer to,
    a provable identity. flux_ratio, not rms_error, is what measures
    this.

A scheme built to answer the first question is under no obligation to
also win the second, and vice versa: a sinc-family kernel's negative
side lobes are what let it resolve fine structure precisely at a
point, but those same negative weights can leave a small, same-signed
residual when many output pixels are later summed into a flux -- a
structural property of interpolation, not a bug, and not something a
better-tuned kernel removes. See docs/user/APP_REFERENCE.md's
"Choosing interp=" for the full discussion and measured numbers.

astropy's reproject_interp only supports order='nearest-neighbor'/
'bilinear'/'biquadratic'/'bicubic' (scipy spline orders 0-3) -- no
sinc/lanczos/gauss/somb/blockave equivalent exists. reproject_exact has
no order= at all (one fixed flux-conserving spherical-polygon-overlap
algorithm, a different method family, not tied to any specific AST
kernel). So the astropy comparison here is deliberately scoped to:
  - interp=nearest   vs reproject_interp(order='nearest-neighbor')
  - interp=linear    vs reproject_interp(order='bilinear')
  - flux conservation (extended-source case): our own best scheme vs
    reproject_exact
Every other AST scheme is benchmarked only against the analytic ground
truth (UR-12), not force-matched to an astropy order that isn't
actually equivalent.

interp=lanczos (the general a-parameter form) is excluded: it needs an
arbitrary interp_params=<support>,<a> with no natural default, unlike
lanczos2/lanczos3 (a fixed, no params needed) -- picking one `a` for
this benchmark would be an arbitrary choice dressed up as a result.

Standalone, on-demand tool (not wired into tests/run_tests.sh): this
reports numbers and a recommendation, it does not assert a pass/fail
tolerance -- see UR-8 for the separate value-level correctness
promotion that will reuse this same ground-truth machinery.

Usage: tests/check_interp_benchmark.py [--keep]
  --keep: keep the generated FITS fixtures/outputs under tests/output/
          (default: cleaned up after the report is written).
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_exact, reproject_interp

REPO_ROOT = Path(__file__).resolve().parent.parent
BIN = REPO_ROOT / "bin" / "reproject_cubes"
OUT_DIR = REPO_ROOT / "tests" / "output" / "interp_benchmark"

NX_IN = NY_IN = 64
CDELT_IN = 1.0 / 60.0  # 1 arcmin/pixel
CRVAL1, CRVAL2 = 180.0, 0.0
CRPIX_IN = (NX_IN + 1) / 2.0  # 32.5 -- pixel-grid centre, not a data point

# B2 (docs/dev/REPROJECT_EXACT_PLAN.md section 0): input and output are
# always the SAME system/projection here (only the pixel-grid geometry
# -- shift/rotation/scale -- differs between "case"s), so there is no
# cross-system conversion to get right, unlike check_target_projection.
# py's B1 -- CRVAL1/2 above are used as literal native coordinates in
# whichever system, needing no sofa_fk4-style correction. (lon_code,
# lat_code, RADESYS, EQUINOX), same convention as check_target_
# projection.py's own SYSTEMS dict.
SYSTEMS = {
    "icrs": ("RA--", "DEC-", None, None),
    "fk5": ("RA--", "DEC-", "FK5", 2000.0),
    "fk4": ("RA--", "DEC-", "FK4", 1950.0),
    "galactic": ("GLON", "GLAT", None, None),
}
PROJECTIONS = ("TAN", "SIN", "CAR", "AIT", "HPX")

# Every interp= this tool offers except lanczos2/lanczos3 (no params
# allowed) and the excluded general lanczos (see module header
# comment). gauss/sincgauss/blockave are NOT left at the omitted-
# params default (confirmed against AST's own mapping.c, 2026-10-13:
# gauss/sincgauss clamp an omitted FWHM to a floor of 0.1 pixels --
# far too narrow, the kernel underflows to zero at almost any
# sub-pixel offset; blockave's own box half-width is not auto-derived
# at all when 0 -- both produce all-NaN output with no warning, not a
# usable default, see docs/dev/REPROJECT_INTERP_OPTIONS_PLAN.md). Each
# gets an explicit, sensible interp_params instead -- FWHM/box
# half-width of a few pixels, matched to this benchmark's own ~1-3
# pixel source scale.
SCHEMES = [
    "nearest", "linear", "sinc", "sincsinc", "sinccos", "sincgauss",
    "blockave", "gauss", "somb", "sombcos", "lanczos2", "lanczos3",
]
SCHEME_PARAMS = {"sincgauss": "0,2.0", "gauss": "0,2.0", "blockave": "2"}

# What is run. Every scheme once with NO interp_params=, i.e. as a user
# gets it: since UR-19 gauss/sincgauss/blockave default to measured
# widths (FWHM 0.8 px, FWHM 8.0 px, half-width 1.0 px, ast_interp_mod.f90)
# instead of refusing. Then those three again at the widths above,
# chosen by hand for this benchmark before UR-19 existed and still used
# by check_interp_correctness.py's tolerances, labelled scheme@params so
# the two are never confused in the report.
DEFAULT_WIDTH_NOTE = {"gauss": "default FWHM 0.8 px", "sincgauss": "default FWHM 8.0 px",
                      "blockave": "default half-width 1.0 px"}
RUNS = [(s, s, None) for s in SCHEMES] + [
    (f"{s}@{p}", s, p) for s, p in SCHEME_PARAMS.items()]

# Scheme -> matched astropy reproject_interp order, only where a real
# equivalent exists (see module header comment).
ASTROPY_EQUIVALENT = {"nearest": "nearest-neighbor", "linear": "bilinear"}


def gaussian_sum(fx, fy, components):
    """components: list of (amp, cx, cy, sigma_x, sigma_y), all in
    input-pixel units. fx/fy may be numpy arrays (vectorised)."""
    val = np.zeros_like(np.asarray(fx, dtype=np.float64))
    for amp, cx, cy, sx, sy in components:
        val += amp * np.exp(-0.5 * (((fx - cx) / sx) ** 2 + ((fy - cy) / sy) ** 2))
    return val


def step_edge(fx, fy):
    """Sharp vertical edge through the input grid's own centre -- 1 on
    one side, 0 on the other, exact at any continuous fx (no smoothing
    in the ground-truth definition itself)."""
    return np.where(np.asarray(fx) > CRPIX_IN, 1.0, 0.0)


# Morphology name -> (value_fn(fx, fy), is_flux_source: bool, components-or-None)
EXTENDED_COMPONENTS = [
    (1.0, 24.0, 20.0, 4.0, 4.0),
    (0.6, 42.0, 44.0, 3.0, 3.0),
]
COMPACT_COMPONENTS = [
    (1.0, 24.0, 20.0, 0.9, 0.9),
    (0.6, 42.0, 44.0, 0.8, 0.8),
]
MORPHOLOGIES = {
    "extended": (lambda fx, fy: gaussian_sum(fx, fy, EXTENDED_COMPONENTS), True, EXTENDED_COMPONENTS),
    "compact": (lambda fx, fy: gaussian_sum(fx, fy, COMPACT_COMPONENTS), True, COMPACT_COMPONENTS),
    "step_edge": (lambda fx, fy: step_edge(fx, fy), False, None),
}


# HPX's own pixel-area convention (AST's definition, independent of
# UR-23's not-yet-implemented Nside-aligned-grid work): at the reference
# pixel, TAN/SIN/CAR/AIT's actual latitude cell size equals |CDELT2|
# directly, but HPX's equals |CDELT2| x 8/(3 pi) (equal-area cells;
# measured and documented already, docs/dev/BRANCH_STATUS_update-
# reproject.md UR-21: 3600.0000 x 3055.8867 arcsec for CDELT=1 deg,
# 3600*8/(3*pi) = 3055.77, matching). Longitude is unaffected. Without
# correcting for this, the SAME CDELT2 number gives a smaller actual
# HPX cell than every other projection, so a beam sized for the other
# four under-samples HPX specifically and UR-21 (correctly) refuses it
# in cases the other four projections accept -- found via B2's own
# extension to all 20 combos, 2026-10-24. Scaling CDELT2 by the inverse
# factor for HPX keeps the actual cell size equal to cdelt_mag for
# every projection.
_HPX_LAT_CORRECTION = 3.0 * np.pi / 8.0

# With the CDELT2 correction above, every projection's nominal cell
# size matches, but "compact"'s own beam (COMPACT_COMPONENTS' narrowest
# sigma) was originally sized to sit just barely above the aligned/
# rotated/upsampled cases' 2.0-cell UR-21 threshold for icrs/SIN alone
# (a knife-edge margin, not deliberate) -- small, genuine differences
# in AST's own per-projection cell-size measurement at the reference
# pixel (not captured by the simple nominal-CDELT correction above)
# put some of the other 19 combos just under that same threshold.
# Rather than replicate AST's own geometry in Python to predict this
# exactly per projection, "compact"'s declared beam gets a modest,
# uniform safety margin instead -- confirmed empirically (not just
# assumed) to clear all 20 combos' aligned/rotated/upsampled cases with
# comfortable margin, while downsampled (2x coarser) still correctly
# refuses for all 20, exactly as intended.
COMPACT_BEAM_MARGIN = 1.2


def cdelt2_mag(cdelt_mag, projection):
    return cdelt_mag * _HPX_LAT_CORRECTION if projection == "HPX" else cdelt_mag


def build_input_header(system="icrs", projection="SIN"):
    lon_code, lat_code, radesys, equinox = SYSTEMS[system]
    h = fits.Header()
    h["SIMPLE"] = True
    h["BITPIX"] = -32
    h["NAXIS"] = 2
    h["NAXIS1"] = NX_IN
    h["NAXIS2"] = NY_IN
    h["CTYPE1"] = f"{lon_code}-{projection}"
    h["CTYPE2"] = f"{lat_code}-{projection}"
    if radesys:
        h["RADESYS"] = radesys
    if equinox:
        h["EQUINOX"] = equinox
    h["CRVAL1"] = CRVAL1
    h["CRVAL2"] = CRVAL2
    h["CRPIX1"] = CRPIX_IN
    h["CRPIX2"] = CRPIX_IN
    h["CDELT1"] = -CDELT_IN
    h["CDELT2"] = cdelt2_mag(CDELT_IN, projection)
    h["BUNIT"] = "Jy/pixel"
    return h


def build_output_header(case, system="icrs", projection="SIN"):
    lon_code, lat_code, radesys, equinox = SYSTEMS[system]
    h = fits.Header()
    h["SIMPLE"] = True
    h["BITPIX"] = -32
    h["NAXIS"] = 2
    if case == "aligned":
        nx = ny = NX_IN
        cdelt = CDELT_IN
        crval1 = CRVAL1 + 0.5 * CDELT_IN  # half-pixel shift -- forces
        crval2 = CRVAL2 + 0.5 * CDELT_IN  # genuine sub-pixel resampling
        crpix = CRPIX_IN
        pc = (1.0, 0.0, 0.0, 1.0)
    elif case == "rotated":
        nx = ny = NX_IN
        cdelt = CDELT_IN
        crval1, crval2 = CRVAL1, CRVAL2
        crpix = CRPIX_IN
        theta = np.radians(30.0)
        pc = (np.cos(theta), -np.sin(theta), np.sin(theta), np.cos(theta))
    elif case == "upsampled":
        nx = ny = NX_IN * 2
        cdelt = CDELT_IN / 2.0
        crval1, crval2 = CRVAL1, CRVAL2
        crpix = (nx + 1) / 2.0
        pc = (1.0, 0.0, 0.0, 1.0)
    elif case == "downsampled":
        nx = ny = NX_IN // 2
        cdelt = CDELT_IN * 2.0
        crval1, crval2 = CRVAL1, CRVAL2
        crpix = (nx + 1) / 2.0
        pc = (1.0, 0.0, 0.0, 1.0)
    else:
        raise ValueError(case)
    h["NAXIS1"] = nx
    h["NAXIS2"] = ny
    h["CTYPE1"] = f"{lon_code}-{projection}"
    h["CTYPE2"] = f"{lat_code}-{projection}"
    if radesys:
        h["RADESYS"] = radesys
    if equinox:
        h["EQUINOX"] = equinox
    h["CRVAL1"] = crval1
    h["CRVAL2"] = crval2
    h["CRPIX1"] = crpix
    h["CRPIX2"] = crpix
    h["CDELT1"] = -cdelt
    h["CDELT2"] = cdelt2_mag(cdelt, projection)
    h["PC1_1"], h["PC1_2"], h["PC2_1"], h["PC2_2"] = pc
    return h


def sample_onto_grid(value_fn, nx, ny):
    """Sample a continuous ground-truth function onto an nx*ny FITS
    pixel grid (1-based pixel centres 1..nx, 1..ny) -- this is what any
    real discretized image is, the same discretization applied to both
    tools' own input."""
    xs = np.arange(1, nx + 1, dtype=np.float64)
    ys = np.arange(1, ny + 1, dtype=np.float64)
    fx, fy = np.meshgrid(xs, ys)
    return value_fn(fx, fy)


def ground_truth_on_output(value_fn, in_header, out_header):
    """Exact analytic value at each output pixel's own exact
    fractional INPUT-pixel coordinate -- never sampled from any
    resampled grid."""
    wcs_in = WCS(in_header)
    wcs_out = WCS(out_header)
    nx_out, ny_out = out_header["NAXIS1"], out_header["NAXIS2"]
    ox, oy = np.meshgrid(np.arange(1, nx_out + 1, dtype=np.float64),
                          np.arange(1, ny_out + 1, dtype=np.float64))
    ra, dec = wcs_out.wcs_pix2world(ox, oy, 1)
    fx, fy = wcs_in.wcs_world2pix(ra, dec, 1)
    return value_fn(fx, fy)


def run_reproject_cubes(infile, reffile, interp, interp_params=None):
    outfile = infile.with_name(infile.stem + "_REPROJ.FITS")
    if outfile.exists():
        outfile.unlink()
    args = [str(BIN), "mode=reference", f"reffile={reffile}", f"infiles={infile}",
            f"interp={interp}"]
    if interp_params is not None:
        args.append(f"interp_params={interp_params}")
    t0 = time.perf_counter()
    r = subprocess.run(args, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    if r.returncode != 0 or not outfile.exists():
        raise RuntimeError(f"reproject_cubes failed for interp={interp}: "
                            f"{r.stdout}\n{r.stderr}")
    data = fits.getdata(outfile).astype(np.float64)
    outfile.unlink()
    # stdout carries log_level=info's per-file stage-timing summary (the
    # exact_weights stage for interp=exact) and, for exact specifically,
    # its own overlap-weight memory lines -- returned, not discarded, so
    # a caller measuring cost (not just accuracy) does not need a second
    # run to get them.
    return data, elapsed, r.stdout


def rms_max_error(output_data, truth, valid_mask):
    diff = output_data[valid_mask] - truth[valid_mask]
    if diff.size == 0:
        return float("nan"), float("nan")
    return float(np.sqrt(np.mean(diff ** 2))), float(np.max(np.abs(diff)))


def total_flux(data, header, valid_mask):
    pixel_area_deg2 = abs(header["CDELT1"]) * abs(header["CDELT2"])
    return float(np.nansum(data[valid_mask])) * pixel_area_deg2


def analytic_total_flux(components):
    # integral of amp*exp(-0.5*((x-cx)/sx)^2-0.5*((y-cy)/sy)^2) over the
    # plane = amp*2*pi*sx*sy (pixel^2 units) -- convert to sky-area
    # units via the INPUT grid's own pixel area (deg^2/pixel^2), which
    # is the physically conserved quantity independent of any later
    # pixel-scale change.
    pixel_area_deg2 = CDELT_IN * CDELT_IN
    return sum(amp * 2 * np.pi * sx * sy for amp, _, _, sx, sy in components) * pixel_area_deg2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--keep", action="store_true")
    args = parser.parse_args()

    if not BIN.exists():
        print(f"ERROR: {BIN} not built -- run: make GPU=0 OMP=1 reproject_cubes")
        return 1

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # B2 (docs/dev/REPROJECT_EXACT_PLAN.md section 0): all 20 system x
    # projection combinations, not just icrs/SIN. Input and output
    # always share the same system/projection here (only the pixel-grid
    # geometry differs between "case"s), so no cross-system conversion
    # (and no sofa_fk4-style correction) is needed for the astropy
    # comparison below, unlike check_target_projection.py's B1.
    results = []
    for system in SYSTEMS:
        for projection in PROJECTIONS:
            in_header = build_input_header(system, projection)
            for morph_name, (value_fn, is_flux_source, components) in MORPHOLOGIES.items():
                infile = OUT_DIR / f"{system}_{projection}_{morph_name}.fits"
                input_data = sample_onto_grid(value_fn, NX_IN, NY_IN)
                # The input's beam is its narrowest feature, as for an observed
                # image. reproject_cubes refuses an output grid coarser than the
                # input whose cells that beam spans fewer than 2.0 times (UR-21),
                # so on the 2x downsampled case only the extended source runs;
                # the compact source and the step edge (no beam: an infinitely
                # sharp edge has no band limit) are refused there, and the
                # ranking below covers the cases that ran.
                morph_header = in_header.copy()
                beam_components = {"extended": EXTENDED_COMPONENTS,
                                    "compact": COMPACT_COMPONENTS}.get(morph_name)
                if beam_components:
                    sigma_px = min(min(c[3], c[4]) for c in beam_components)
                    margin = COMPACT_BEAM_MARGIN if morph_name == "compact" else 1.0
                    morph_header["BMAJ"] = (sigma_px * 2.0 * np.sqrt(2.0 * np.log(2.0)) *
                                             CDELT_IN * margin)
                    morph_header["BMIN"] = morph_header["BMAJ"]
                    morph_header["BPA"] = 0.0
                fits.PrimaryHDU(data=input_data.astype(np.float32), header=morph_header).writeto(
                    infile, overwrite=True)

                for case in ("aligned", "rotated", "upsampled", "downsampled"):
                    out_header = build_output_header(case, system, projection)
                    reffile = OUT_DIR / f"{system}_{projection}_{morph_name}_{case}_ref.fits"
                    nx_out, ny_out = out_header["NAXIS1"], out_header["NAXIS2"]
                    fits.PrimaryHDU(data=np.zeros((ny_out, nx_out), dtype=np.float32),
                                     header=out_header).writeto(reffile, overwrite=True)

                    truth = ground_truth_on_output(value_fn, in_header, out_header)
                    valid_truth_mask = np.isfinite(truth)

                    for scheme, run_scheme, run_params in RUNS:
                        try:
                            data, elapsed, _ = run_reproject_cubes(infile, reffile, run_scheme,
                                                                    run_params)
                        except RuntimeError as exc:
                            print(f"WARNING: {system}/{projection}/{morph_name}/{case}/"
                                  f"{scheme} failed: {exc}", file=sys.stderr)
                            continue
                        valid = valid_truth_mask & np.isfinite(data)
                        rms, worst = rms_max_error(data, truth, valid)
                        row = {"system": system, "projection": projection,
                               "morphology": morph_name, "case": case, "scheme": scheme,
                               "tool": "reproject_cubes", "rms_error": rms, "max_error": worst,
                               "seconds": elapsed, "valid_frac": float(np.mean(valid))}
                        if is_flux_source:
                            row["flux_ratio"] = (total_flux(data, out_header, valid) /
                                                  analytic_total_flux(components))
                        results.append(row)

                        astropy_order = ASTROPY_EQUIVALENT.get(scheme)
                        if astropy_order is not None:
                            t0 = time.perf_counter()
                            a_data, _ = reproject_interp((input_data, in_header), out_header,
                                                          order=astropy_order)
                            a_elapsed = time.perf_counter() - t0
                            a_valid = valid_truth_mask & np.isfinite(a_data)
                            a_rms, a_worst = rms_max_error(a_data, truth, a_valid)
                            a_row = {"system": system, "projection": projection,
                                     "morphology": morph_name, "case": case, "scheme": scheme,
                                     "tool": f"astropy:{astropy_order}", "rms_error": a_rms,
                                     "max_error": a_worst, "seconds": a_elapsed,
                                     "valid_frac": float(np.mean(a_valid))}
                            if is_flux_source:
                                a_row["flux_ratio"] = (total_flux(a_data, out_header, a_valid) /
                                                        analytic_total_flux(components))
                            results.append(a_row)

                    if is_flux_source:
                        t0 = time.perf_counter()
                        e_data, _ = reproject_exact((input_data, in_header), out_header)
                        e_elapsed = time.perf_counter() - t0
                        e_valid = valid_truth_mask & np.isfinite(e_data)
                        e_rms, e_worst = rms_max_error(e_data, truth, e_valid)
                        results.append({
                            "system": system, "projection": projection,
                            "morphology": morph_name, "case": case, "scheme": "(exact)",
                            "tool": "astropy:reproject_exact", "rms_error": e_rms,
                            "max_error": e_worst, "seconds": e_elapsed,
                            "valid_frac": float(np.mean(e_valid)),
                            "flux_ratio": total_flux(e_data, out_header, e_valid) /
                            analytic_total_flux(components),
                        })

    report_path = OUT_DIR / "report.json"
    report_path.write_text(json.dumps(results, indent=2))
    print_report(results)
    print(f"\nFull results: {report_path}")

    if not args.keep:
        for f in OUT_DIR.glob("*.fits"):
            f.unlink()
    return 0


def print_report(results):
    # B2: with all 20 system x projection combinations, a flat per-row
    # table runs to thousands of lines -- report.json already has every
    # row for anyone who needs that detail. What's useful to read here
    # instead is: (1) a per-combination sanity check, so a system or
    # projection that behaves differently from the rest stands out, and
    # (2) the same aggregate rankings as before, now pooled across all
    # 20 combinations rather than just icrs/SIN.
    print(f"{'system':<10} {'projection':<10} {'morph':<9} {'best scheme':<16} {'mean RMS':>12} "
          f"{'worst scheme':<16} {'mean RMS':>12}")
    print("-" * 90)
    for system in SYSTEMS:
        for projection in PROJECTIONS:
            for morph in ("extended", "compact"):
                by_scheme_combo: dict[str, list[float]] = {}
                for row in results:
                    if (row["system"] == system and row["projection"] == projection and
                            row["tool"] == "reproject_cubes" and row["morphology"] == morph):
                        by_scheme_combo.setdefault(row["scheme"], []).append(row["rms_error"])
                if not by_scheme_combo:
                    continue
                combo_ranking = sorted(by_scheme_combo.items(), key=lambda kv: np.mean(kv[1]))
                best_scheme, best_errs = combo_ranking[0]
                worst_scheme, worst_errs = combo_ranking[-1]
                print(f"{system:<10} {projection:<10} {morph:<9} {best_scheme:<16} "
                      f"{np.mean(best_errs):>12.3e} {worst_scheme:<16} {np.mean(worst_errs):>12.3e}")

    # Reported per morphology, never pooled: pooling extended (sigma
    # 3-4px, real curvature across a pixel) with compact (sigma 0.8-
    # 0.9px, far less curvature) folds two different-sized point-vs-
    # area-average gaps into one number, which can shift a scheme's
    # apparent ranking for reasons unrelated to its own interpolation
    # accuracy (docs/user/APP_REFERENCE.md's "Choosing interp=" has the
    # full explanation). Excludes step_edge throughout -- a sharp
    # discontinuity's Gibbs-ringing error is large for every scheme
    # except the degenerate nearest-neighbor case and would dominate/
    # distort either ranking.
    for morph in ("compact", "extended"):
        print(f"\n-- RMS ranking (interpolation accuracy): {morph} morphology only --")
        by_scheme: dict[str, list[float]] = {}
        for row in results:
            if row["tool"] == "reproject_cubes" and row["morphology"] == morph:
                by_scheme.setdefault(row["scheme"], []).append(row["rms_error"])
        ranking = sorted(by_scheme.items(), key=lambda kv: np.mean(kv[1]))
        for scheme, errs in ranking:
            secs = [r["seconds"] for r in results if r["scheme"] == scheme and
                    r["tool"] == "reproject_cubes" and r["morphology"] == morph]
            if "@" in scheme:
                at = f"  at interp_params={scheme.split('@', 1)[1]} (explicit)"
            elif scheme in DEFAULT_WIDTH_NOTE:
                at = f"  at its {DEFAULT_WIDTH_NOTE[scheme]} (no interp_params=)"
            else:
                at = ""
            print(f"  {scheme:<16} mean RMS {np.mean(errs):.4e}  "
                  f"mean seconds {np.mean(secs):.3f}{at}")
    print("  (test images are small (64x64) -- per-call overhead dominates; these"
          " timings are NOT a reliable cost comparison at production cube scale.)")
    # Whoever reads these rankings needs the parameter values in front
    # of them, not only in this file's own source: the three schemes
    # that take a kernel width rank where that width puts them, and a
    # different width is a different ranking.
    print("  (gauss/sincgauss/blockave appear twice: once with no interp_params=, at"
          " the default width a user gets, and once as name@params at a width chosen"
          " by hand for this benchmark. Their ranking depends on the width; the other"
          " nine schemes take none.)")

    print("\n-- Flux ratio (conservative-remapping accuracy) by morphology"
          " (compact/extended reported separately, not pooled -- a different question"
          " from the RMS ranking above: total integrated flux recovered, not pointwise"
          " reconstruction accuracy; `exact` conserves this by construction regardless"
          " of morphology) --")
    for morph in ("compact", "extended"):
        by_scheme_flux: dict[str, list[float]] = {}
        for row in results:
            if (row["tool"] == "reproject_cubes" and row["morphology"] == morph and
                    "flux_ratio" in row):
                by_scheme_flux.setdefault(row["scheme"], []).append(row["flux_ratio"])
        if not by_scheme_flux:
            continue
        print(f"  {morph}:")
        for scheme, ratios in sorted(by_scheme_flux.items(),
                                      key=lambda kv: abs(1.0 - np.mean(kv[1]))):
            print(f"    {scheme:<16} mean flux ratio {np.mean(ratios):.6f}")

    print("\n-- UR-13: our error / astropy error, matched schemes only (ratio near 1.0"
          " means comparable to astropy's own established algorithm on the identical"
          " case; no genuine astropy equivalent exists for the other AST schemes;"
          " pooled across all 20 system x projection combinations) --")
    for scheme, astropy_tool in (("nearest", "astropy:nearest-neighbor"),
                                  ("linear", "astropy:bilinear")):
        ratios = []
        for case in ("aligned", "rotated", "upsampled", "downsampled"):
            for morph in MORPHOLOGIES:
                theirs_rows = {(r["system"], r["projection"]): r["rms_error"]
                                for r in results if r["scheme"] == scheme and
                                r["tool"] == astropy_tool and r["case"] == case and
                                r["morphology"] == morph}
                for r in results:
                    if not (r["scheme"] == scheme and r["tool"] == "reproject_cubes" and
                             r["case"] == case and r["morphology"] == morph):
                        continue
                    theirs = theirs_rows.get((r["system"], r["projection"]))
                    if theirs is not None and theirs > 0:
                        ratios.append(r["rms_error"] / theirs)
        if ratios:
            print(f"  {scheme:<8} our_rms/astropy_rms: mean={np.mean(ratios):.4f} "
                  f"min={np.min(ratios):.4f} max={np.max(ratios):.4f}")

    print("\n-- UR-13: flux conservation (extended source) -- our best scheme vs"
          " reproject_exact, astropy's own dedicated flux-conserving algorithm;"
          " pooled across all 20 system x projection combinations --")
    best_scheme = ranking[0][0]
    for case in ("aligned", "rotated", "upsampled", "downsampled"):
        exact_by_combo = {(r["system"], r["projection"]): r["flux_ratio"] for r in results
                           if r["morphology"] == "extended" and r["case"] == case and
                           r["tool"] == "astropy:reproject_exact"}
        ours_by_combo = {(r["system"], r["projection"]): r["flux_ratio"] for r in results
                          if r["morphology"] == "extended" and r["case"] == case and
                          r["scheme"] == best_scheme and r["tool"] == "reproject_cubes"}
        combos = sorted(set(exact_by_combo) & set(ours_by_combo))
        if combos:
            exact_mean = np.mean([exact_by_combo[c] for c in combos])
            ours_mean = np.mean([ours_by_combo[c] for c in combos])
            print(f"  {case:<12} reproject_exact={exact_mean:.6f}  "
                  f"our_{best_scheme}={ours_mean:.6f}")


if __name__ == "__main__":
    sys.exit(main())
