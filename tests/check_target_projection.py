#!/usr/bin/env python3
"""target_projection= on the no-reffile path (UR-7/UR-8,
docs/dev/BRANCH_STATUS_update-reproject.md).

For each projection, target system, field centre and infile count, runs
reproject_cubes and match_cubes stages=reproject with no reffile and
checks the output against an independent reprojection: astropy's
reproject_interp(order='bilinear') of the same input onto the header
our tool wrote. UR-13 found our interp=linear equivalent to that
algorithm, so on a correct grid the two agree to float32 rounding,
while a wrong projection, orientation or footprint shows up as a large
difference.

Checked per run:
  - exit 0, and the output CTYPE carries the requested system and
    projection (or the inherited system, for projection alone);
  - values agree with astropy within ABS_TOL wherever both are finite;
  - no pixel is finite in ours where astropy has no value, i.e. the
    tool never invents data outside the input's coverage;
  - match_cubes is pixel-identical to reproject_cubes.
Plus the error paths: an unknown projection value, and projection alone
over infiles in different sky systems. Exit status 0 only if every
check passes; each failure is printed.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
from scipy.ndimage import map_coordinates

from sofa_fk4 import fk4_to_icrs

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / "bin"
WORK = REPO / "tests" / "output" / "target_projection"
N = 64
CDELT = 1.0 / 60.0          # 1 arcmin/pixel, a ~1 deg field
# Largest difference measured before this test existed, over the same
# matrix: 2.44e-6 on data of order 1. 1e-5 leaves headroom for compiler
# and platform rounding while staying far below any geometric error,
# which shows up at the level of the data itself (order 1e-1 or more).
ABS_TOL = 1.0e-5

# hpx excluded (UR-23, 2026-10-24): before UR-23, hpx accepted this
# matrix's arbitrary 60 arcsec CDELT like every other projection here,
# with no real HEALPix alignment. Now its cell size is quantised to
# Nside = 90/(sqrt(2)*CDELT), which this matrix's CDELT has no close
# power-of-2 match for within the verified Nside<=1024 range, and this
# matrix declares no beam -- so the no-beam safety gate correctly
# refuses rather than silently coarsening. hpx's own geometry/value
# correctness (pixel centres against astropy_healpix, values against
# astropy, default-Nside selection, beam-safety margin, every error
# path) is covered on its own, appropriate terms in
# tests/check_hpx_nside.py instead of being forced into this matrix.
PROJECTIONS = ("tan", "sin", "car", "ait")
# (lon_code, lat_code, RADESYS, EQUINOX) -- RADESYS/EQUINOX are None for
# icrs/galactic (neither needs a card: icrs is the default with none
# given, galactic has no equinox concept).
SYSTEMS = {
    "icrs": ("RA--", "DEC-", None, None),
    "fk5": ("RA--", "DEC-", "FK5", 2000.0),
    "fk4": ("RA--", "DEC-", "FK4", 1950.0),
    "galactic": ("GLON", "GLAT", None, None),
}
CENTRES = {"equator": (180.0, 0.0), "dec_m60": (80.0, -60.0)}

failures = []
n_checks = 0


def check(ok, message):
    global n_checks
    n_checks += 1
    if not ok:
        failures.append(message)


def make_input(path, ra0, dec0, lon_code="RA--", lat_code="DEC-", radesys=None, equinox=None):
    h = fits.Header()
    h["SIMPLE"] = True; h["BITPIX"] = -32; h["NAXIS"] = 2
    h["NAXIS1"] = N; h["NAXIS2"] = N
    h["CTYPE1"] = f"{lon_code}-SIN"; h["CTYPE2"] = f"{lat_code}-SIN"
    if radesys:
        h["RADESYS"] = radesys
    if equinox:
        h["EQUINOX"] = equinox
    h["CRVAL1"] = ra0; h["CRVAL2"] = dec0
    h["CRPIX1"] = (N + 1) / 2; h["CRPIX2"] = (N + 1) / 2
    h["CDELT1"] = -CDELT; h["CDELT2"] = CDELT
    y, x = np.mgrid[1:N + 1, 1:N + 1].astype(float)
    # smooth and asymmetric, so a flip or rotation cannot pass
    data = (np.exp(-((x - 22) ** 2 + (y - 40) ** 2) / (2 * 6.0 ** 2))
            + 0.6 * np.exp(-((x - 44) ** 2 + (y - 20) ** 2) / (2 * 9.0 ** 2))
            + 0.002 * x)
    fits.PrimaryHDU(data.astype(np.float32), h).writeto(path, overwrite=True)


def run(tool, args):
    cmd = [str(BIN / tool)] + list(args)
    if tool == "match_cubes":
        cmd.insert(1, "stages=reproject")
    return subprocess.run(cmd, capture_output=True, text=True)


def output_of(infile):
    return infile.with_name(infile.stem + "_REPROJ.FITS")


def output_pixel_to_input_pixel(infile, outfile, ix, iy):
    """Output pixel (ix, iy) -> input pixel (x, y), the one conversion
    both the reference-value comparison and the coverage check below
    need. When the output is FK4 and the input is genuinely a
    different system (every case this test exercises: ICRS input,
    target_system=fk4), the output sky position is converted to the
    input's frame via sofa_fk4 instead of astropy's own transform_to,
    which differs from AST by up to ~47 mas (REPROJECT_EXACT_PLAN.md
    sections 5d/5e/5f) -- enough to matter here because this test
    compares directly against astropy, unlike the UR-24 matrix tests'
    own arbiter. When input and output are *both* FK4 (projection
    alone, system inherited, on an FK4 infile -- the one case this
    test has where they match), no system conversion happens at all,
    so this is deliberately skipped: forcing it through ICRS and back
    would introduce a spurious round-trip, not avoid one. Every other
    system is untouched, exactly as before."""
    hin = fits.getheader(infile)
    hout = fits.getheader(outfile)
    win = WCS(hin).celestial
    wout = WCS(hout).celestial
    sky = wout.pixel_to_world(ix, iy)
    if hout.get("RADESYS") == "FK4" and hin.get("RADESYS") != "FK4":
        equinox_byear = hout.get("EQUINOX", 1950.0)
        ra_icrs, dec_icrs = fk4_to_icrs(sky.ra.deg, sky.dec.deg, equinox_byear)
        sky = SkyCoord(ra_icrs * u.deg, dec_icrs * u.deg, frame="icrs")
    return win.world_to_pixel(sky)


def compare_with_astropy(infile, outfile, label):
    ours = fits.getdata(outfile).astype(np.float64)
    ny, nx = ours.shape
    jj, ii = np.mgrid[0:ny, 0:nx]
    xi, yi = output_pixel_to_input_pixel(infile, outfile, ii, jj)
    if fits.getheader(outfile).get("RADESYS") == "FK4":
        data_in = fits.getdata(infile).astype(np.float64)
        ref = map_coordinates(data_in, [yi, xi], order=1, cval=np.nan, mode="constant")
    else:
        ref, _ = reproject_interp((fits.getdata(infile), WCS(fits.getheader(infile))),
                                  WCS(fits.getheader(outfile)).celestial,
                                  shape_out=ours.shape, order="bilinear")
    fo, fr = np.isfinite(ours), np.isfinite(ref)
    both = fo & fr
    check(both.any(), f"{label}: no pixel finite in both ours and astropy")
    if both.any():
        diff = float(np.max(np.abs(ours[both] - ref[both])))
        check(diff <= ABS_TOL, f"{label}: max |ours - astropy| = {diff:.3e} > {ABS_TOL:.0e}")
    # Coverage is the input's pixel AREA, [0.5, N+0.5] in FITS pixel
    # coordinates -- not astropy's finite mask, which stops at the last
    # pixel CENTRE. Where an output pixel centre lands exactly on the
    # input's outer pixel edge (it does when the output grid sits half a
    # pixel off the input's, as for sin onto the infile's own system),
    # AST returns the edge pixel's value and astropy returns NaN; that
    # point is inside the input's coverage, so it is not counted here.
    hin = fits.getheader(infile)
    jy, ix = np.nonzero(fo)
    px, py = output_pixel_to_input_pixel(infile, outfile, ix, jy)
    x1, y1 = np.asarray(px) + 1.0, np.asarray(py) + 1.0
    edge_tol = 1.0e-6
    outside = ((x1 < 0.5 - edge_tol) | (x1 > hin["NAXIS1"] + 0.5 + edge_tol) |
               (y1 < 0.5 - edge_tol) | (y1 > hin["NAXIS2"] + 0.5 + edge_tol) |
               ~np.isfinite(x1) | ~np.isfinite(y1))
    n_outside = int(np.sum(outside))
    check(n_outside == 0,
          f"{label}: {n_outside} finite output pixel(s) map outside the input's pixel area")


def one_case(proj, system, centre_name, n_files):
    ra0, dec0 = CENTRES[centre_name]
    lon_code, lat_code, radesys, equinox = SYSTEMS[system]
    case = f"{proj}/{system}/{centre_name}/{n_files}file"
    outputs = {}
    for tool in ("reproject_cubes", "match_cubes"):
        d = WORK / tool / proj / system / centre_name / f"{n_files}file"
        d.mkdir(parents=True, exist_ok=True)
        infiles = [d / "in1.fits"]
        make_input(infiles[0], ra0, dec0)
        if n_files == 2:
            # second field offset by ~0.4 deg, overlapping the first, so
            # the grid has to cover a union of footprints
            infiles.append(d / "in2.fits")
            make_input(infiles[1], ra0 + 0.4 / np.cos(np.radians(dec0)), dec0 + 0.2)
        for f in infiles:
            output_of(f).unlink(missing_ok=True)
        r = run(tool, [f"target_system={system}", f"target_projection={proj}",
                       "infiles=" + ",".join(str(f) for f in infiles), "interp=linear"])
        label = f"{tool} {case}"
        check(r.returncode == 0, f"{label}: exit {r.returncode}: "
              + next((l.strip() for l in r.stdout.splitlines() if "ERROR" in l), ""))
        if r.returncode != 0:
            continue
        for f in infiles:
            out = output_of(f)
            check(out.exists(), f"{label}: no output for {f.name}")
            if not out.exists():
                continue
            h = fits.getheader(out)
            want1 = f"{lon_code}-{proj.upper()}"
            want2 = f"{lat_code}-{proj.upper()}"
            check(h["CTYPE1"] == want1 and h["CTYPE2"] == want2,
                  f"{label}: CTYPE {h['CTYPE1']}/{h['CTYPE2']}, expected {want1}/{want2}")
            # Explicit, direct check -- UR-7/UR-8 (2026-10-24): the
            # no-reffile path silently wrote RADESYS=ICRS regardless of
            # the requested system, for every equatorial target_system
            # (CTYPE alone is the same string for ICRS/FK5/FK4, so it
            # could not catch this; only the value comparison below
            # caught it, indirectly, for FK4 specifically). This
            # asserts the header directly, so a regression here is
            # caught immediately rather than inferred from a value
            # mismatch that only shows up for FK4.
            if radesys:
                check(h.get("RADESYS") == radesys,
                      f"{label}: RADESYS {h.get('RADESYS')!r}, expected {radesys!r}")
                check(h.get("EQUINOX") == equinox,
                      f"{label}: EQUINOX {h.get('EQUINOX')!r}, expected {equinox!r}")
            elif system == "icrs":
                check(h.get("RADESYS") in (None, "ICRS"),
                      f"{label}: RADESYS {h.get('RADESYS')!r}, expected ICRS or absent")
            compare_with_astropy(f, out, f"{label} {f.name}")
            outputs.setdefault(f.name, {})[tool] = fits.getdata(out)
    for name, by_tool in outputs.items():
        if len(by_tool) == 2:
            check(np.array_equal(by_tool["reproject_cubes"], by_tool["match_cubes"],
                                 equal_nan=True),
                  f"{case} {name}: match_cubes differs from reproject_cubes")


def inheritance_and_errors():
    d = WORK / "inherit"
    d.mkdir(parents=True, exist_ok=True)
    for tool in ("reproject_cubes", "match_cubes"):
        # projection alone keeps the infile's own system
        for system, (lon_code, lat_code, radesys, equinox) in SYSTEMS.items():
            f = d / f"{tool}_{system}.fits"
            make_input(f, 80.0, -60.0, lon_code, lat_code, radesys, equinox)
            output_of(f).unlink(missing_ok=True)
            r = run(tool, ["target_projection=car", f"infiles={f}", "interp=linear"])
            label = f"{tool} projection alone on a {system} infile"
            check(r.returncode == 0, f"{label}: exit {r.returncode}")
            if r.returncode == 0 and output_of(f).exists():
                h = fits.getheader(output_of(f))
                check(h["CTYPE1"] == f"{lon_code}-CAR",
                      f"{label}: CTYPE1 {h['CTYPE1']}, expected {lon_code}-CAR (system kept)")
                if radesys:
                    check(h.get("RADESYS") == radesys,
                          f"{label}: RADESYS {h.get('RADESYS')!r}, expected {radesys!r} (system kept)")
                    check(h.get("EQUINOX") == equinox,
                          f"{label}: EQUINOX {h.get('EQUINOX')!r}, expected {equinox!r} (system kept)")
                elif system == "icrs":
                    check(h.get("RADESYS") in (None, "ICRS"),
                          f"{label}: RADESYS {h.get('RADESYS')!r}, expected ICRS or absent")
                compare_with_astropy(f, output_of(f), label)
        # unknown projection value
        f = d / f"{tool}_bad.fits"
        make_input(f, 180.0, 0.0)
        output_of(f).unlink(missing_ok=True)
        r = run(tool, ["target_system=icrs", "target_projection=mol", f"infiles={f}"])
        check(r.returncode != 0 and "target_projection must be" in r.stdout
              and not output_of(f).exists(),
              f"{tool}: target_projection=mol was not refused with a named error")
        # projection alone over infiles in different systems
        fa, fg = d / f"{tool}_mix_icrs.fits", d / f"{tool}_mix_gal.fits"
        make_input(fa, 80.0, -60.0)
        make_input(fg, 80.0, -60.0, "GLON", "GLAT")
        for f in (fa, fg):
            output_of(f).unlink(missing_ok=True)
        r = run(tool, ["target_projection=car", f"infiles={fa},{fg}"])
        check(r.returncode != 0 and "different" in r.stdout
              and not output_of(fa).exists() and not output_of(fg).exists(),
              f"{tool}: projection alone over icrs+galactic infiles was not refused")


def lattice_alignment():
    """UR-22: when system, projection and pixel scale all match the
    infile's, the no-reffile grid must keep the infile's size, pixel
    centres and header, for even and odd sizes alike. Before UR-22 an
    even-sized infile came back one pixel larger, its pixel centres half
    a pixel off (the grid's reference pixel centre sat on the infile's
    centroid, a pixel corner when the size is even). Values are compared
    to round-off: AST composes the two mappings without reducing them
    to an exact identity, which leaves offsets of order 1e-11 px."""
    d = WORK / "lattice"
    d.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(1)
    for tool in ("reproject_cubes", "match_cubes"):
        for n in (32, 33):
            for interp in ("linear", "lanczos3"):
                f = d / f"{tool}_{n}_{interp}.fits"
                h = fits.Header()
                h["SIMPLE"] = True; h["BITPIX"] = -32; h["NAXIS"] = 2
                h["NAXIS1"] = n; h["NAXIS2"] = n
                h["CTYPE1"] = "RA---SIN"; h["CTYPE2"] = "DEC--SIN"
                h["CRVAL1"] = 180.0; h["CRVAL2"] = 0.0
                h["CRPIX1"] = (n + 1) / 2; h["CRPIX2"] = (n + 1) / 2
                h["CDELT1"] = -0.001; h["CDELT2"] = 0.001
                data = rng.standard_normal((n, n)).astype(np.float32)
                fits.PrimaryHDU(data, h).writeto(f, overwrite=True)
                output_of(f).unlink(missing_ok=True)
                r = run(tool, ["target_system=icrs", f"infiles={f}", f"interp={interp}"])
                label = f"{tool} {n}x{n} {interp} identity grid"
                check(r.returncode == 0 and output_of(f).exists(), f"{label}: exit {r.returncode}")
                if r.returncode != 0 or not output_of(f).exists():
                    continue
                ho = fits.getheader(output_of(f))
                out = fits.getdata(output_of(f))
                check(out.shape == (n, n), f"{label}: output {out.shape[::-1]}, expected ({n}, {n})")
                for k in ("CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2"):
                    check(abs(ho[k] - h[k]) <= 1e-9, f"{label}: {k} {ho[k]} differs from the infile's {h[k]}")
                if out.shape == (n, n):
                    diff = float(np.max(np.abs(out.astype(float) - data.astype(float))))
                    check(diff <= 1e-9, f"{label}: max |output - infile| = {diff:.2e}, expected round-off")


def main():
    for tool in ("reproject_cubes", "match_cubes"):
        if not (BIN / tool).exists():
            print(f"ERROR: {BIN / tool} not built")
            return 1
    WORK.mkdir(parents=True, exist_ok=True)
    for proj in PROJECTIONS:
        for system in SYSTEMS:
            for centre_name in CENTRES:
                for n_files in (1, 2):
                    one_case(proj, system, centre_name, n_files)
    inheritance_and_errors()
    lattice_alignment()
    print(f"target_projection checks: {n_checks - len(failures)}/{n_checks} passed")
    for msg in failures:
        print(f"  FAIL: {msg}")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
