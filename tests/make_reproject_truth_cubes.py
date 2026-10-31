#!/usr/bin/env python3
"""
Fixtures for reproject_cubes value-level correctness tests
(docs/dev/REPROJECT_TEST_PLAN.md), regime: both the reference grid and
every input grid are well sampled with respect to the beam (>=3 pixels
per FWHM) -- the only regime docs/dev/REPROJECT_FLUX_NOTE.md's
conclusion covers. The large-pixel-scale-factor regime (REPROJECT_FLUX_
NOTE.md section 8) is out of scope here on purpose.

Grid G0 is the reference. G1/G2/G3/G5 are the same sky content (sources,
a ramp, a flat field, a constant-angle Q/U pair, a NaN block), each
rendered directly on that grid's own pixel grid from the shared analytic
model (not by reprojecting G0 -- the point is an independently computed
truth), with a header relating it to G0 by a known linear pixel map:

  G1  shift (0.3, 0.7) pixel          -- sub-pixel alignment
  G2  rotate 30 degrees               -- rotation (circular sources only,
                                          see CIRC_FWHM comment below)
  G3  pixel scale x1.05 + shift       -- a small scale difference
  G5  shift by a whole pixel          -- exact-copy sanity check

The source beam matches docs/dev/REPROJECT_FLUX_NOTE.md's dancingghosts
measurements exactly (8.66 x 7.15 pixels on G0, at this fixture's own
0.001 deg/pixel scale, matching tests/make_test_cubes.py's convention),
so that note's already-measured thresholds (REPROJECT_TEST_PLAN.md
section 3) apply unchanged to this fixture.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from astropy.io import fits

OUTDIR = Path(__file__).resolve().parent / "data"

NX, NY = 160, 160
PIXEL_DEG = 0.001          # matches tests/make_test_cubes.py
CRVAL1, CRVAL2 = 180.0, 0.0
CRPIX0 = NX // 2 + 1        # 81, same for both axes (NX == NY)

# Elliptical beam matching REPROJECT_FLUX_NOTE.md section 2's
# dancingghosts measurement, expressed in sky degrees so every grid's
# own pixel-FWHM is correct for ITS OWN pixel scale (BMAJ/BMIN are sky
# quantities; pixel-FWHM = sky-FWHM / that grid's own CDELT).
BMAJ_DEG = 8.66 * PIXEL_DEG   # 8.66 pixels on G0 specifically
BMIN_DEG = 7.15 * PIXEL_DEG   # 7.15 pixels on G0 specifically

# Circular beam (geometric mean of the above) used ONLY for the rotation
# grid (G2): a circular Gaussian is rotation-invariant, so the source
# model needs no position-angle bookkeeping under a 30-degree pixel-grid
# rotation -- isolates the resampling behaviour under rotation from a
# separate (and separately testable) position-angle-rotation question.
CIRC_FWHM_DEG = float(np.sqrt(BMAJ_DEG * BMIN_DEG))

# Source layout: 3x3 grid, spacing far enough apart (>4 FWHM, the
# largest FWHM in pixels is ~8.7) that circular apertures of radius
# 3*FWHM (~26 px) never overlap, margin >=3*FWHM from the field edge.
# Fixed seed; sub-pixel jitter up to +-5 pixels around each nominal
# centre.
SEED = 20261010
_rng = np.random.default_rng(SEED)
_centres = [40.0, 80.0, 120.0]
SOURCES = []
for _iy, _cy in enumerate(_centres):
    for _ix, _cx in enumerate(_centres):
        _dx, _dy = _rng.uniform(-5.0, 5.0, size=2)
        SOURCES.append({
            "name": f"src_{_iy}_{_ix}",
            "x0": _cx + _dx,   # pixel position on G0, 0-indexed
            "y0": _cy + _dy,
            "amp": 1.0,
        })

# Ramp truth: I = A + B*x + C*y (x, y are 0-indexed G0 pixel coordinates)
RAMP_A, RAMP_B, RAMP_C = 0.3, 0.004, 0.006

# NaN block: an 8x8 region of the "flat" image, centred, for the
# footprint-margin check (REPROJECT_TEST_PLAN.md check 7).
NAN_BLOCK_HALF = 4   # 8x8 block = 2*NAN_BLOCK_HALF on a side

# Polarisation angle for the constant-angle Q/U pair (check 11), degrees.
POL_ANGLE_DEG = 35.0
POL_P0 = 1.0


def _base_header(nx: int, ny: int, crval1: float, crval2: float,
                  crpix1: float, crpix2: float, cdelt: float,
                  rotation_deg: float = 0.0) -> fits.Header:
    """4-axis header, same convention as tests/make_test_cubes.py's
    make_header (RA---SIN/DEC--SIN/FREQ/STOKES, Jy/beam), single
    channel/Stokes plane (this fixture does not feed rm_synthesis, so a
    realistic frequency axis is not needed -- only the sky WCS matters
    to reproject_cubes)."""
    h = fits.Header()
    h["SIMPLE"] = True
    h["BITPIX"] = -32
    h["NAXIS"] = 4
    h["NAXIS1"] = nx
    h["NAXIS2"] = ny
    h["NAXIS3"] = 1
    h["NAXIS4"] = 1
    h["CTYPE1"] = "RA---SIN"
    h["CTYPE2"] = "DEC--SIN"
    h["CTYPE3"] = "FREQ"
    h["CTYPE4"] = "STOKES"
    h["CRVAL1"] = crval1
    h["CRVAL2"] = crval2
    h["CRVAL3"] = 1.0e9
    h["CRVAL4"] = 1.0
    h["CRPIX1"] = crpix1
    h["CRPIX2"] = crpix2
    h["CRPIX3"] = 1.0
    h["CRPIX4"] = 1.0
    h["CDELT1"] = -cdelt
    h["CDELT2"] = cdelt
    h["CDELT3"] = 1.0e6
    h["CDELT4"] = 1.0
    h["CUNIT1"] = "deg"
    h["CUNIT2"] = "deg"
    h["CUNIT3"] = "Hz"
    h["BUNIT"] = "Jy/beam"
    h["OBJECT"] = "TEST_REPROJECT_TRUTH"
    h["EQUINOX"] = 2000.0
    if rotation_deg != 0.0:
        # CROTA2 convention (CDELT1/2 + CROTA2), matching what
        # copy_sky_rotation_matrix/copy_axis_keywords in
        # src/reproject_cubes.f90 read (it reads CROTA, falling back to
        # PC/CD only if present -- no PC/CD keyword is written here, so
        # the CROTA path is what gets exercised).
        h["CROTA2"] = rotation_deg
    h["HISTORY"] = "Synthetic fixture for reproject_cubes correctness tests"
    return h


def _grid_def(dx_pix=0.0, dy_pix=0.0, scale=1.0, rotation_deg=0.0):
    """Pixel grid parameters for a grid related to G0 by: rotate about
    G0's own reference pixel by rotation_deg, scale the pixel size by
    `scale`, then shift the reference pixel by (dx_pix, dy_pix) G0
    pixels. CRVAL stays fixed (all grids share the same sky origin;
    CRPIX/CDELT/CROTA carry the relationship) -- this keeps the pixel
    map between any two grids exactly linear and known, as
    REPROJECT_TEST_PLAN.md section 2 specifies."""
    cdelt = PIXEL_DEG * scale
    crpix = CRPIX0 / scale - dx_pix / scale, CRPIX0 / scale - dy_pix / scale
    return {"cdelt": cdelt, "crpix1": crpix[0], "crpix2": crpix[1],
            "rotation_deg": rotation_deg}


GRIDS = {
    "G0": _grid_def(),
    "G1": _grid_def(dx_pix=0.3, dy_pix=0.7),
    "G2": _grid_def(rotation_deg=30.0),
    "G3": _grid_def(dx_pix=0.3, dy_pix=0.7, scale=1.05),
    "G5": _grid_def(dx_pix=1.0, dy_pix=1.0),
}


def _world_of_g0_pixel(x0: float, y0: float) -> tuple:
    """RA/Dec (degrees) of a G0 pixel position (0-indexed), via G0's own
    tangent-plane (SIN) projection -- small field (160 pixels * 0.001 deg
    = 0.16 deg), so the flat tangent-plane approximation used here
    (consistent with reproject_cubes.f90's own stated scope: "a field
    under 0.15 degree where the linear mapping agrees with SIN", see
    REPROJECT_FLUX_NOTE.md section 9) is adequate for constructing a
    known truth, not a general-purpose SIN deprojection."""
    dx_deg = (x0 - (CRPIX0 - 1)) * (-PIXEL_DEG)
    dy_deg = (y0 - (CRPIX0 - 1)) * PIXEL_DEG
    ra = CRVAL1 + dx_deg / np.cos(np.radians(CRVAL2))
    dec = CRVAL2 + dy_deg
    return ra, dec


def _pixel_of_world(grid_name: str, ra: float, dec: float) -> tuple:
    """Inverse of _world_of_g0_pixel, generalised to any grid in GRIDS
    (same flat-tangent-plane approximation, same small-field scope).

    Matches the standard FITS CDELT+CROTA2 convention exactly (the CD
    matrix it implies: CD1_1=CDELT1*cos(CROTA2), CD1_2=-CDELT2*sin(CROTA2),
    CD2_1=CDELT1*sin(CROTA2), CD2_2=CDELT2*cos(CROTA2); xi=CD1_1*p1+CD1_2*p2,
    eta=CD2_1*p1+CD2_2*p2), inverted algebraically for (p1,p2) given
    (xi,eta). Earlier draft used a plain isotropic rotation of (px,py)
    instead, which is NOT equivalent once CDELT1 and CDELT2 have opposite
    signs (standard: CDELT1<0, CDELT2>0) -- that sign difference is a
    reflection, and composing a reflection with a rotation does not
    commute the way composing two rotations does. Caught by check
    G2 (the only rotated grid) failing hard (99% integral error, 18-pixel
    centroid error) while the non-rotated grids (G1/G3/G5) passed exactly
    -- the isotropic-rotation version was correct only in the degenerate
    rotation_deg=0 case, which is every grid except G2."""
    g = GRIDS[grid_name]
    xi = (ra - CRVAL1) * np.cos(np.radians(CRVAL2))   # dx_deg
    eta = dec - CRVAL2                                 # dy_deg
    theta = np.radians(g["rotation_deg"])
    cdelt1, cdelt2 = -g["cdelt"], g["cdelt"]
    p1 = (xi * np.cos(theta) + eta * np.sin(theta)) / cdelt1
    p2 = (-xi * np.sin(theta) + eta * np.cos(theta)) / cdelt2
    x = p1 + (g["crpix1"] - 1)
    y = p2 + (g["crpix2"] - 1)
    return x, y


def _render_sources(grid_name: str, fwhm_x_g0_pix: float,
                     fwhm_y_g0_pix: float, pa_deg_on_g0: float = 0.0
                     ) -> np.ndarray:
    """Render SOURCES on grid_name's own pixel grid. fwhm_*_g0_pix are
    the beam's FWHM in G0 pixels (converted internally to this grid's
    own pixel scale via BMAJ/BMIN degrees, so the same sky beam is
    represented correctly on every grid). pa_deg_on_g0 rotates with the
    grid for grids with a rotation (so an elliptical beam stays aligned
    with the pixel grid the same way G0's does) -- unused for the
    circular-beam case (fwhm_x_g0_pix == fwhm_y_g0_pix), since a
    circular Gaussian has no orientation to track."""
    g = GRIDS[grid_name]
    scale = g["cdelt"] / PIXEL_DEG
    fwhm_x = fwhm_x_g0_pix / scale
    fwhm_y = fwhm_y_g0_pix / scale
    pa = pa_deg_on_g0 + g["rotation_deg"]
    sig_x = fwhm_x / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    sig_y = fwhm_y / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    yy, xx = np.mgrid[0:NY, 0:NX].astype(np.float64)
    img = np.zeros((NY, NX), dtype=np.float32)
    theta = np.radians(pa)
    ct, st = np.cos(theta), np.sin(theta)
    for src in SOURCES:
        ra, dec = _world_of_g0_pixel(src["x0"], src["y0"])
        xc, yc = _pixel_of_world(grid_name, ra, dec)
        dx = xx - xc
        dy = yy - yc
        dxr = dx * ct + dy * st
        dyr = -dx * st + dy * ct
        img += src["amp"] * np.exp(-0.5 * ((dxr / sig_x) ** 2 +
                                            (dyr / sig_y) ** 2))
    return img


def _render_ramp(grid_name: str) -> np.ndarray:
    """I = RAMP_A + RAMP_B*x_g0 + RAMP_C*y_g0, evaluated at each of this
    grid's own pixels' G0-equivalent position -- a linear function of
    sky position stays linear (though with different coefficients) in
    any other grid's own pixel coordinates too, but it is simplest and
    equally valid to define and evaluate the truth directly in G0 pixel
    units throughout. Forward (pixel -> sky -> G0 pixel) counterpart of
    _pixel_of_world's inverse transform -- same CD-matrix convention,
    see that function's docstring for why the naive isotropic-rotation
    version this replaced was wrong for any rotated grid."""
    g = GRIDS[grid_name]
    yy, xx = np.mgrid[0:NY, 0:NX].astype(np.float64)
    p1 = xx - (g["crpix1"] - 1)
    p2 = yy - (g["crpix2"] - 1)
    theta = np.radians(g["rotation_deg"])
    cdelt1, cdelt2 = -g["cdelt"], g["cdelt"]
    xi = cdelt1 * np.cos(theta) * p1 - cdelt2 * np.sin(theta) * p2
    eta = cdelt1 * np.sin(theta) * p1 + cdelt2 * np.cos(theta) * p2
    x0 = -xi / PIXEL_DEG + (CRPIX0 - 1)
    y0 = eta / PIXEL_DEG + (CRPIX0 - 1)
    return (RAMP_A + RAMP_B * x0 + RAMP_C * y0).astype(np.float32)


def _write(path: Path, data: np.ndarray, grid_name: str, extra: dict | None = None):
    g = GRIDS[grid_name]
    hdr = _base_header(NX, NY, CRVAL1, CRVAL2, g["crpix1"], g["crpix2"],
                        g["cdelt"], g["rotation_deg"])
    hdr["BMAJ"] = BMAJ_DEG
    hdr["BMIN"] = BMIN_DEG
    hdr["BPA"] = 0.0
    if extra:
        for k, v in extra.items():
            hdr[k] = v
    out = data.reshape(1, 1, NY, NX)
    fits.PrimaryHDU(data=out, header=hdr).writeto(path, overwrite=True)


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    truth = {
        "pixel_deg": PIXEL_DEG, "nx": NX, "ny": NY, "crpix0": CRPIX0,
        "bmaj_deg": BMAJ_DEG, "bmin_deg": BMIN_DEG,
        "nb_g0": float(np.pi / (4.0 * np.log(2.0)) * BMAJ_DEG * BMIN_DEG /
                       (PIXEL_DEG * PIXEL_DEG)),
        "ramp": {"a": RAMP_A, "b": RAMP_B, "c": RAMP_C},
        "pol_angle_deg": POL_ANGLE_DEG, "pol_p0": POL_P0,
        "nan_block_half": NAN_BLOCK_HALF,
        "sources": [{"name": s["name"], "x0": s["x0"], "y0": s["y0"],
                     "amp": s["amp"]} for s in SOURCES],
        "grids": {},
    }

    for gname, g in GRIDS.items():
        truth["grids"][gname] = {k: float(v) for k, v in g.items()}
        is_rot = gname == "G2"
        fwhm_x = CIRC_FWHM_DEG / PIXEL_DEG if is_rot else BMAJ_DEG / PIXEL_DEG
        fwhm_y = CIRC_FWHM_DEG / PIXEL_DEG if is_rot else BMIN_DEG / PIXEL_DEG

        srcimg = _render_sources(gname, fwhm_x, fwhm_y)
        _write(OUTDIR / f"REPROJ_{gname}.SOURCES.FITS", srcimg, gname)

        ramp = _render_ramp(gname)
        _write(OUTDIR / f"REPROJ_{gname}.RAMP.FITS", ramp, gname)

        flat = np.ones((NY, NX), dtype=np.float32)
        flat_nan = flat.copy()
        c = NX // 2
        flat_nan[c - NAN_BLOCK_HALF:c + NAN_BLOCK_HALF,
                 c - NAN_BLOCK_HALF:c + NAN_BLOCK_HALF] = np.nan
        _write(OUTDIR / f"REPROJ_{gname}.FLAT.FITS", flat, gname)
        _write(OUTDIR / f"REPROJ_{gname}.FLATNAN.FITS", flat_nan, gname)

        # Constant polarisation angle: single compact source so the
        # angle check (test plan 11) probes the same resampling the
        # other checks do, not an extended/uniform field.
        pol_src = _render_sources(gname, fwhm_x, fwhm_y)
        chi = np.radians(POL_ANGLE_DEG)
        q = pol_src * POL_P0 * np.cos(2.0 * chi)
        u = pol_src * POL_P0 * np.sin(2.0 * chi)
        _write(OUTDIR / f"REPROJ_{gname}.Q.FITS", q.astype(np.float32), gname)
        _write(OUTDIR / f"REPROJ_{gname}.U.FITS", u.astype(np.float32), gname)

    (OUTDIR / "reproject_truth.json").write_text(json.dumps(truth, indent=2))
    print(f"Wrote reproject_cubes correctness fixtures to {OUTDIR}")


if __name__ == "__main__":
    main()
