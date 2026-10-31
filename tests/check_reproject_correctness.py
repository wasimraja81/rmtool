#!/usr/bin/env python3
"""
Value-level correctness checks for reproject_cubes (docs/dev/
REPROJECT_TEST_PLAN.md), regime: reference and input grids both well
sampled with respect to the beam. Reads the *_REPROJ.FITS outputs
already produced (by the calling test-suite section) from the fixtures
tests/make_reproject_truth_cubes.py writes, and reports [OK]/[FAIL] per
check in the same style as check_rm_peak.py. Thresholds below are the
ones proposed in REPROJECT_TEST_PLAN.md section 3, from the measurements
in REPROJECT_FLUX_NOTE.md -- both for AST__LINEAR, the only interpolation
scheme reproject_cubes offers today.
"""
from __future__ import annotations
import json
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits

DATA = Path(__file__).resolve().parent / "data"
TRUTH = json.loads((DATA / "reproject_truth.json").read_text())
NB = TRUTH["nb_g0"]
APERTURE_R = 15   # pixels; > 3*FWHM (largest FWHM here ~8.7 px)

ok = True


def report(flag_ok: bool, msg: str):
    global ok
    print(f"[{'OK' if flag_ok else 'FAIL'}] {msg}")
    if not flag_ok:
        ok = False


def load(name: str) -> np.ndarray:
    return fits.getdata(DATA / name).squeeze().astype(np.float64)


def check_ramp():
    """Check 1: I = a + b*x + c*y must reproduce exactly (tol 1e-4
    relative to the ramp's own range) away from the edges, for the
    sub-pixel-shift grid G1."""
    g0 = load("REPROJ_G0.RAMP.FITS")
    out = load("REPROJ_G1.RAMP_REPROJ.FITS")
    margin = 10
    interior = slice(margin, -margin)
    diff = np.abs(out[interior, interior] - g0[interior, interior])
    rng = g0.max() - g0.min()
    worst = float(np.nanmax(diff)) / rng
    report(worst <= 1e-4, f"ramp (G1 vs G0, interior): worst relative "
           f"error {worst:.2e} (tol 1e-4)")


def check_whole_pixel_shift():
    """Check 2: G5's own CRPIX is G0's shifted by a whole pixel, so a
    given SKY position sits at a different pixel INDEX on the two grids,
    but reprojecting G5 onto the G0 reference grid must reproduce G0's
    own truth array directly, position for position, once both are
    expressed on G0's grid -- that is the entire point of reprojection,
    and holds for any input grid sharing the same sky, not just a
    whole-pixel-shifted one. What makes G5 special, and worth its own
    check, is that the resample lands exactly on an input grid point, so
    linear interpolation should reproduce it with none of the peak loss
    a sub-pixel offset (G1, G3) causes."""
    g0 = load("REPROJ_G0.SOURCES.FITS")
    out = load("REPROJ_G5.SOURCES_REPROJ.FITS")
    margin = 10
    a = out[margin:-margin, margin:-margin]
    b = g0[margin:-margin, margin:-margin]
    diff = np.abs(a - b)
    worst = float(np.nanmax(diff))
    report(worst <= 1e-4, f"whole-pixel shift (G5 vs G0, interior): "
           f"worst absolute error {worst:.2e} (tol 1e-4)")


def _source_photometry(data: np.ndarray, x0: float, y0: float):
    xi, yi = int(round(x0)), int(round(y0))
    r = APERTURE_R
    sub = data[yi - r:yi + r, xi - r:xi + r]
    integral = float(np.nansum(sub)) / NB
    peak = float(np.nanmax(sub))
    # Flux-weighted centroid, in absolute pixel coordinates.
    yy, xx = np.mgrid[yi - r:yi + r, xi - r:xi + r].astype(np.float64)
    w = np.nan_to_num(sub, nan=0.0)
    w = np.clip(w, 0.0, None)
    wsum = w.sum()
    cx = float((w * xx).sum() / wsum)
    cy = float((w * yy).sum() / wsum)
    return integral, peak, cx, cy


def check_photometry_and_position(grid_name: str, integral_tol=0.005,
                                   peak_min=0.95):
    g0 = load("REPROJ_G0.SOURCES.FITS")
    out = load(f"REPROJ_{grid_name}.SOURCES_REPROJ.FITS")
    ratios, peaks, poserrs = [], [], []
    for src in TRUTH["sources"]:
        x0, y0 = src["x0"], src["y0"]
        int_truth, _, _, _ = _source_photometry(g0, x0, y0)
        int_out, peak_out, cx, cy = _source_photometry(out, x0, y0)
        ratios.append(int_out / int_truth)
        peaks.append(peak_out / src["amp"])
        poserrs.append(float(np.hypot(cx - x0, cy - y0)))
    ratios = np.array(ratios)
    peaks = np.array(peaks)
    poserrs = np.array(poserrs)

    report(bool(np.all(np.abs(ratios - 1.0) <= integral_tol)),
           f"{grid_name}: integral ratio per source within "
           f"{integral_tol*100:.1f}% (worst {np.max(np.abs(ratios-1.0))*100:.3f}%)")
    report(bool(np.all(peaks >= peak_min)),
           f"{grid_name}: peak ratio >= {peak_min} for every source "
           f"(mean {peaks.mean():.4f}, min {peaks.min():.4f})")
    report(bool(np.all(poserrs <= 0.02)),
           f"{grid_name}: centroid position error <= 0.02 pixel "
           f"(worst {poserrs.max():.4f})")


def check_nan_footprint():
    """Check 7 (regression guard, not an exactness claim -- the current
    scheme has no bad-value flag set, see REPROJECT_FLUX_NOTE.md section
    5): count NaN output pixels from an 8x8 NaN input block under G1's
    (0.3, 0.7) pixel sub-shift, and assert it stays near the measured
    baseline rather than growing without bound. 241 px measured here
    (not REPROJECT_FLUX_NOTE.md's own 81 px figure -- that was for a
    different, unspecified shift; margin growth from bilinear
    interpolation scales with how many input pixels each output pixel
    draws from, which depends on the sub-pixel offset, so the two
    numbers are not expected to match and neither is wrong)."""
    out = load("REPROJ_G1.FLATNAN_REPROJ.FITS")
    n_nan = int(np.isnan(out).sum())
    input_block = (2 * TRUTH["nan_block_half"]) ** 2
    baseline = 241
    report(n_nan <= int(baseline * 1.1),
           f"NaN footprint (G1, 8x8={input_block}px input block, "
           f"(0.3,0.7)px shift): {n_nan} output NaN pixels (regression "
           f"guard: measured baseline {baseline}, allow <= "
           f"{int(baseline*1.1)})")


def check_unit_handling():
    """Check 10: a flat field of 1.0 Jy/beam must stay 1.0 after
    reprojection through a small scale change (G3, x1.05), away from
    the edges -- BUNIT is never touched and values are not scaled by
    the pixel-area ratio."""
    out = load("REPROJ_G3.FLAT_REPROJ.FITS")
    margin = 10
    interior = out[margin:-margin, margin:-margin]
    worst = float(np.nanmax(np.abs(interior - 1.0)))
    report(worst <= 0.01, f"unit handling (G3, flat field, interior): "
           f"worst deviation from 1.0 is {worst:.4f} (tol 0.01)")


def check_nearest_ramp():
    """Tier 1 (interp=nearest, docs/dev/REPROJECT_INTERP_OPTIONS_PLAN.md
    section 5): AST__NEAREST must pick a single G1 input pixel per
    output pixel, not blend neighbours the way linear does. G1's own
    grid relates to G0 by a pure (dx_pix, dy_pix) pixel shift (no
    rotation/scale -- asserted below, since the formula here only holds
    for that case): the G0->G1 pixel map is x_G1 = x_G0 - dx_pix,
    y_G1 = y_G0 - dy_pix, so every output pixel's nearest-neighbour
    value differs from G0's own ramp truth by one constant offset,
    RAMP_B*(dx_pix - round(dx_pix)) + RAMP_C*(dy_pix - round(dy_pix)) --
    a closed-form prediction, not just "this differs from the
    linear-interp baseline." """
    g1 = TRUTH["grids"]["G1"]
    scale = g1["cdelt"] / TRUTH["pixel_deg"]
    assert g1["rotation_deg"] == 0.0 and abs(scale - 1.0) < 1e-12, (
        "check_nearest_ramp's closed-form offset assumes G1 is a pure "
        "shift (no rotation/scale) -- fixture changed, formula needs "
        "revisiting")
    dx_pix = TRUTH["crpix0"] - g1["crpix1"] * scale
    dy_pix = TRUTH["crpix0"] - g1["crpix2"] * scale
    ramp_b, ramp_c = TRUTH["ramp"]["b"], TRUTH["ramp"]["c"]
    const = (ramp_b * (dx_pix - round(dx_pix)) +
             ramp_c * (dy_pix - round(dy_pix)))

    g0 = load("REPROJ_G0.RAMP.FITS")
    out = load("REPROJ_G1.RAMP_REPROJ_NEAREST.FITS")
    margin = 10
    interior = slice(margin, -margin)
    expected = g0[interior, interior] + const
    worst = float(np.nanmax(np.abs(out[interior, interior] - expected)))
    report(worst <= 1e-4, f"interp=nearest (G1 ramp): worst absolute "
           f"error vs. closed-form prediction {worst:.2e} (tol 1e-4)")


def check_polarization_angle(grid_name: str, tol_deg=2.0):
    """Check 11: for a compact source with constant polarisation angle,
    0.5*atan2(U,Q) at the source peak must match the input angle -- AST's
    linear-interpolation weights are real and positive, so Q and U are
    scaled by the same local weights at any given output pixel and the
    angle is preserved."""
    q = load(f"REPROJ_{grid_name}.Q_REPROJ.FITS")
    u = load(f"REPROJ_{grid_name}.U_REPROJ.FITS")
    errs = []
    for src in TRUTH["sources"]:
        x0, y0 = src["x0"], src["y0"]
        xi, yi = int(round(x0)), int(round(y0))
        r = 3
        qs = np.nansum(q[yi - r:yi + r, xi - r:xi + r])
        us = np.nansum(u[yi - r:yi + r, xi - r:xi + r])
        ang = np.degrees(0.5 * np.arctan2(us, qs)) % 180.0
        expect = TRUTH["pol_angle_deg"] % 180.0
        err = min(abs(ang - expect), 180.0 - abs(ang - expect))
        errs.append(err)
    errs = np.array(errs)
    report(bool(np.all(errs <= tol_deg)),
           f"{grid_name}: polarisation angle preserved within "
           f"{tol_deg} deg (worst {errs.max():.3f} deg)")


def main():
    check_ramp()
    check_whole_pixel_shift()
    for gname in ("G1", "G2", "G3"):
        check_photometry_and_position(gname)
    check_nan_footprint()
    check_unit_handling()
    check_polarization_angle("G1")
    sys.exit(0 if ok else 1)


def main_interp():
    """Tier 1 interp=/interp_params= checks (section 63 in
    tests/run_tests.sh), run separately from main()'s fixed sequence
    above since they read extra *_REPROJ_NEAREST.FITS outputs that
    main()'s own section-62 caller never produces."""
    check_nearest_ramp()
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--interp":
        main_interp()
    else:
        main()
