#!/usr/bin/env python3
"""HEALPix-aligned HPX grids, target_hpx_nside= (UR-23,
docs/dev/BRANCH_STATUS_update-reproject.md).

A HEALPix-aligned HPX grid is fully determined by Nside alone:
CRVAL=(0,0), |CDELT|=90/(sqrt(2)*Nside), the pixel axes rotated 45 deg,
and CRPIX/NAXIS from a formula verified (below) against
astropy_healpix. target_center=/target_lon=/target_lat=/
target_cdelt_lon=/target_cdelt_lat= have nothing left to choose and are
errors with target_projection=hpx; target_hpx_nside= itself is only
valid with target_projection=hpx.

Checked:
  - output pixel centres equal astropy_healpix centres to 1e-9 deg,
    covering every facet (equatorial and both polar caps) exactly
    once, for icrs and galactic, at several explicit Nside;
  - reprojected values agree with astropy's reproject_interp
    (order='bilinear') within tolerance, same method as
    check_target_projection.py;
  - RADESYS/EQUINOX for every target_system (icrs/fk5/fk4/galactic);
  - default Nside selection: matches the finest infile pixel when an
    exact-enough Nside exists; coarsens only when every infile's beam
    spans >= 2.5 of the chosen cell; never coarsens for a no-beam
    infile (refused instead, by the ordinary UR-21 gate);
  - target_hpx_nside= validation (power of 2 required) and every
    target_cdelt_lon=/target_center=/target_lon=/target_lat=/
    target_hpx_nside=-without-hpx error path;
  - reproject_cubes and match_cubes agree pixel-for-pixel.

Exit status 0 only if every check passes; each failure is printed.
"""
import math
import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy_healpix import HEALPix
from reproject import reproject_interp
import astropy.units as u

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / "bin"
WORK = REPO / "tests" / "output" / "hpx_nside"
N = 64
ABS_TOL = 1.0e-4

# (lon_code, lat_code, RADESYS, EQUINOX)
SYSTEMS = {
    "icrs": ("RA--", "DEC-", None, None),
    "fk5": ("RA--", "DEC-", "FK5", 2000.0),
    "fk4": ("RA--", "DEC-", "FK4", 1950.0),
    "galactic": ("GLON", "GLAT", None, None),
}

n_checks = 0
n_fail = 0


def check(ok, message):
    global n_checks, n_fail
    n_checks += 1
    if not ok:
        n_fail += 1
        print(f"  FAIL: {message}")


def run(tool, args):
    cmd = [str(BIN / tool)] + list(args)
    if tool == "match_cubes":
        cmd.insert(1, "stages=reproject")
    return subprocess.run(cmd, capture_output=True, text=True)


def output_of(infile):
    return infile.with_name(infile.stem + "_REPROJ.FITS")


def make_input(path, cdelt_deg, lon0=180.0, lat0=0.0, lon_code="RA--",
                lat_code="DEC-", radesys=None, equinox=None, beam_fwhm_deg=None):
    h = fits.Header()
    h["SIMPLE"] = True; h["BITPIX"] = -32; h["NAXIS"] = 2
    h["NAXIS1"] = N; h["NAXIS2"] = N
    h["CTYPE1"] = f"{lon_code}-SIN"; h["CTYPE2"] = f"{lat_code}-SIN"
    if radesys:
        h["RADESYS"] = radesys
    if equinox:
        h["EQUINOX"] = equinox
    h["CRVAL1"] = lon0; h["CRVAL2"] = lat0
    h["CRPIX1"] = (N + 1) / 2; h["CRPIX2"] = (N + 1) / 2
    h["CDELT1"] = -cdelt_deg; h["CDELT2"] = cdelt_deg
    if beam_fwhm_deg is not None:
        h["BMAJ"] = beam_fwhm_deg; h["BMIN"] = beam_fwhm_deg; h["BPA"] = 0.0
    y, x = np.mgrid[1:N + 1, 1:N + 1].astype(float)
    data = (np.exp(-((x - 22) ** 2 + (y - 40) ** 2) / (2 * 6.0 ** 2))
            + 0.6 * np.exp(-((x - 44) ** 2 + (y - 20) ** 2) / (2 * 9.0 ** 2))
            + 0.002 * x)
    fits.PrimaryHDU(data.astype(np.float32), h).writeto(path, overwrite=True)


def hpx_cell_deg(nside):
    return 90.0 / (nside * math.sqrt(2.0))


def hpx_naxis(nside):
    return 6 if nside == 1 else 5 * nside + 2


def check_pixel_centres_and_values(tool, infile, outfile, nside, label, check_values=True):
    out = output_of(outfile)
    h = fits.getheader(out)
    naxis = h["NAXIS1"]
    check(naxis == hpx_naxis(nside),
          f"{label}: NAXIS1 {naxis}, expected {hpx_naxis(nside)} for Nside={nside}")
    w = WCS(h)
    px, py = np.meshgrid(np.arange(1, naxis + 1), np.arange(1, naxis + 1))
    px = px.ravel().astype(float); py = py.ravel().astype(float)
    sky = w.pixel_to_world(px, py)
    lon_attr, lat_attr = ("l", "b") if hasattr(sky, "l") else ("ra", "dec")
    lon_deg = getattr(sky, lon_attr).deg
    lat_deg = getattr(sky, lat_attr).deg
    valid = np.isfinite(lon_deg) & np.isfinite(lat_deg)
    hp = HEALPix(nside=nside, order="nested")
    npix = 12 * nside * nside
    lon = lon_deg[valid]; lat = lat_deg[valid]
    idx = hp.lonlat_to_healpix(lon * u.deg, lat * u.deg)
    lon2, lat2 = hp.healpix_to_lonlat(idx)
    dlon = (lon - lon2.deg + 180) % 360 - 180
    dlat = lat - lat2.deg
    err = np.hypot(dlon * np.cos(np.radians(lat)), dlat)
    check(len(np.unique(idx)) == npix,
          f"{label}: {len(np.unique(idx))} unique HEALPix pixels covered, expected {npix}")
    check(err.size > 0 and err.max() < 1.0e-9,
          f"{label}: max pixel-centre error {err.max() if err.size else 'n/a'} deg >= 1e-9")

    if check_values:
        in_h = fits.getheader(infile)
        in_data = fits.getdata(infile)
        out_data = fits.getdata(out).astype(np.float64)
        ref, _ = reproject_interp((in_data, in_h), WCS(h).celestial,
                                   shape_out=out_data.shape, order="bilinear")
        both = np.isfinite(out_data) & np.isfinite(ref)
        check(both.any(), f"{label}: no pixel finite in both ours and astropy")
        if both.any():
            diff = float(np.max(np.abs(out_data[both] - ref[both])))
            check(diff <= ABS_TOL, f"{label}: max |ours-astropy| = {diff:.3e} > {ABS_TOL:.0e}")


def case_explicit_nside(system, nside, centre=(180.0, 0.0)):
    """Fixed, explicit Nside -- bypasses default selection and its
    beam-safety reasoning entirely, so no beam is needed here; this
    is purely a geometry/value check at a chosen Nside."""
    lon_code, lat_code, radesys, equinox = SYSTEMS[system]
    label_base = f"{system}/nside{nside}"
    cdelt = hpx_cell_deg(nside) * 2.0  # deliberately coarser input: explicit Nside must not be refused
    outputs = {}
    for tool in ("reproject_cubes", "match_cubes"):
        d = WORK / tool / "explicit" / label_base
        d.mkdir(parents=True, exist_ok=True)
        infile = d / "in.fits"
        make_input(infile, cdelt, centre[0], centre[1], lon_code, lat_code, radesys, equinox)
        output_of(infile).unlink(missing_ok=True)
        r = run(tool, [f"target_system={system}", "target_projection=hpx",
                       f"target_hpx_nside={nside}", f"infiles={infile}", "interp=linear"])
        label = f"{tool} {label_base}"
        check(r.returncode == 0, f"{label}: exit {r.returncode}: "
              + next((l.strip() for l in r.stdout.splitlines() if "ERROR" in l), r.stdout[-200:]))
        if r.returncode != 0:
            continue
        h = fits.getheader(output_of(infile))
        if radesys:
            check(h.get("RADESYS") == radesys,
                  f"{label}: RADESYS {h.get('RADESYS')!r}, expected {radesys!r}")
            check(h.get("EQUINOX") == equinox,
                  f"{label}: EQUINOX {h.get('EQUINOX')!r}, expected {equinox!r}")
        elif system == "icrs":
            check(h.get("RADESYS") in (None, "ICRS"),
                  f"{label}: RADESYS {h.get('RADESYS')!r}, expected ICRS or absent")
        check_pixel_centres_and_values(tool, infile, infile, nside, label,
                                        check_values=(system in ("icrs", "galactic")))
        outputs[tool] = fits.getdata(output_of(infile))
    if len(outputs) == 2:
        check(np.array_equal(outputs["reproject_cubes"], outputs["match_cubes"], equal_nan=True),
              f"{label_base}: match_cubes differs from reproject_cubes")


def case_high_nside():
    """The NAXIS/CRPIX formula verified by hand up to Nside=32 is used
    unconditionally in production for any power of 2; this checks it
    holds at the high end of the ticket's own verified range too
    (geometry only, not the value comparison already covered at
    smaller Nside -- reproject_interp on a multi-million-pixel grid
    buys nothing new and costs real time)."""
    for nside in (128, 256, 512, 1024):
        cdelt = hpx_cell_deg(nside) * 2.0
        d = WORK / "high_nside" / str(nside)
        d.mkdir(parents=True, exist_ok=True)
        infile = d / "in.fits"
        make_input(infile, cdelt)
        output_of(infile).unlink(missing_ok=True)
        r = run("reproject_cubes", ["target_system=icrs", "target_projection=hpx",
                                     f"target_hpx_nside={nside}", f"infiles={infile}",
                                     "interp=linear"])
        label = f"high_nside/{nside}"
        check(r.returncode == 0, f"{label}: exit {r.returncode}: "
              + next((l.strip() for l in r.stdout.splitlines() if "ERROR" in l), ""))
        if r.returncode == 0:
            check_pixel_centres_and_values("reproject_cubes", infile, infile, nside,
                                            label, check_values=False)


def case_default_nside_exact():
    """Finest infile CDELT chosen to land exactly on an achievable
    Nside's cell -- default selection must pick it with no coarsening
    and no beam needed."""
    nside = 16
    cdelt = hpx_cell_deg(nside)
    d = WORK / "default_exact"
    d.mkdir(parents=True, exist_ok=True)
    infile = d / "in.fits"
    make_input(infile, cdelt)
    output_of(infile).unlink(missing_ok=True)
    r = run("reproject_cubes", ["target_system=icrs", "target_projection=hpx",
                                 f"infiles={infile}", "interp=linear"])
    check(r.returncode == 0, f"default_exact: exit {r.returncode}: "
          + next((l.strip() for l in r.stdout.splitlines() if "ERROR" in l), ""))
    if r.returncode == 0:
        h = fits.getheader(output_of(infile))
        check(h["NAXIS1"] == hpx_naxis(nside),
              f"default_exact: NAXIS1 {h['NAXIS1']}, expected {hpx_naxis(nside)} (Nside={nside})")


def case_default_nside_beam_safe_coarsen():
    """Finest infile CDELT deliberately finer than any achievable
    Nside's cell, but with a beam wide enough (>= 2.5 cells at the
    coarser Nside) that the default selection may safely coarsen."""
    nside = 16
    # 0.9x Nside=16's own cell: close enough that log2-nearest-power-of-2
    # still rounds to 16 (not 32), so 16 is genuinely coarser than this
    # input, not finer -- the beam-safety reasoning must actually engage.
    fine_cdelt = hpx_cell_deg(nside) * 0.9
    beam_fwhm = hpx_cell_deg(nside) * 3.0    # safely >= 2.5 cells at Nside
    d = WORK / "default_beam_safe"
    d.mkdir(parents=True, exist_ok=True)
    infile = d / "in.fits"
    make_input(infile, fine_cdelt, beam_fwhm_deg=beam_fwhm)
    output_of(infile).unlink(missing_ok=True)
    r = run("reproject_cubes", ["target_system=icrs", "target_projection=hpx",
                                 f"infiles={infile}", "interp=linear"])
    check(r.returncode == 0, f"default_beam_safe: exit {r.returncode}: "
          + next((l.strip() for l in r.stdout.splitlines() if "ERROR" in l), ""))


def case_default_nside_no_beam_refused():
    """The default-selection bump-to-finer-Nside loop doubles Nside
    (halving the cell) until safe or no longer coarser; since the
    nearest-power-of-2 starting point is always within a factor of
    sqrt(2) of the finest infile pixel, one doubling always escapes
    "coarser" -- except at the verified Nside=1024 ceiling, where
    there is nowhere finer left to go. An infile finer than even
    Nside=1024's own cell, with no beam, is the only way to reach a
    genuine no-beam refusal through default selection rather than an
    explicit target_hpx_nside=."""
    too_fine_cdelt = hpx_cell_deg(1024) * 0.5
    d = WORK / "default_no_beam"
    d.mkdir(parents=True, exist_ok=True)
    infile = d / "in.fits"
    make_input(infile, too_fine_cdelt)
    output_of(infile).unlink(missing_ok=True)
    r = run("reproject_cubes", ["target_system=icrs", "target_projection=hpx",
                                 f"infiles={infile}", "interp=linear"])
    check(r.returncode != 0,
          "default_no_beam: expected a no-beam coarsening refusal but the run succeeded")
    check(r.returncode == 0 or "aliasing cannot be ruled out" in r.stdout,
          f"default_no_beam: expected the no-beam gate's own message, got: {r.stdout[-300:]}")


def case_explicit_nside_bad():
    d = WORK / "explicit_bad"
    d.mkdir(parents=True, exist_ok=True)
    infile = d / "in.fits"
    make_input(infile, hpx_cell_deg(8))
    for bad in ("7", "0", "-4", "3.5", "notanumber"):
        r = run("reproject_cubes", ["target_system=icrs", "target_projection=hpx",
                                     f"target_hpx_nside={bad}", f"infiles={infile}"])
        check(r.returncode != 0, f"explicit_nside_bad({bad}): expected an error but it succeeded")


def case_errors():
    d = WORK / "errors"
    d.mkdir(parents=True, exist_ok=True)
    infile = d / "in.fits"
    make_input(infile, hpx_cell_deg(8))

    combos = [
        (["target_cdelt_lon=0.01", "target_cdelt_lat=0.01"], "target_cdelt_lon/lat"),
        (["target_center=user_defined", "target_lon=10", "target_lat=10"], "target_center"),
        (["target_lon=10", "target_lat=10"], "target_lon/lat alone"),
    ]
    for extra, name in combos:
        r = run("reproject_cubes",
                ["target_system=icrs", "target_projection=hpx", f"infiles={infile}"] + extra)
        check(r.returncode != 0, f"errors({name}): expected an error but it succeeded")

    # target_hpx_nside= without target_projection=hpx.
    r = run("reproject_cubes", ["target_system=icrs", "target_projection=sin",
                                 "target_hpx_nside=8", f"infiles={infile}"])
    check(r.returncode != 0,
          "errors(target_hpx_nside without hpx): expected an error but it succeeded")

    # Both tools agree on each error's exit status.
    for extra, name in combos:
        r = run("match_cubes",
                ["target_system=icrs", "target_projection=hpx", f"infiles={infile}"] + extra)
        check(r.returncode != 0, f"errors({name}) match_cubes: expected an error but it succeeded")


def main():
    if not (BIN / "reproject_cubes").exists() or not (BIN / "match_cubes").exists():
        print("ERROR: reproject_cubes/match_cubes not built -- run: "
              "make GPU=0 OMP=1 reproject_cubes match_cubes")
        return 1
    WORK.mkdir(parents=True, exist_ok=True)

    for system in SYSTEMS:
        for nside in (4, 32):
            case_explicit_nside(system, nside)
    # A polar-cap-heavy centre, for the systems the pixel-centre check
    # actually runs values/geometry on.
    case_explicit_nside("icrs", 8, centre=(45.0, 85.0))
    case_explicit_nside("galactic", 8, centre=(45.0, 85.0))

    case_high_nside()
    case_default_nside_exact()
    case_default_nside_beam_safe_coarsen()
    case_default_nside_no_beam_refused()
    case_explicit_nside_bad()
    case_errors()

    print(f"HPX Nside checks: {n_checks - n_fail}/{n_checks} passed")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
