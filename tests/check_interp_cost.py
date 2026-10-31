#!/usr/bin/env python3
"""Accuracy versus cost: guidance for users (UR-24 section 7,
docs/dev/REPROJECT_EXACT_PLAN.md).

Compute cost is not an acceptance criterion and no target is set --
this is a standalone, on-demand report (like check_interp_benchmark.py,
not wired into tests/run_tests.sh), giving measured numbers for users
choosing between interp=exact and the interpolation kernels.

One representative source, in both of check_interp_benchmark.py's own
"extended" (sigma 3-4px) and "compact" (sigma 0.8-0.9px) morphologies,
per case -- not the full scheme sweep check_interp_benchmark.py already
covers. Methods: interp=exact, lanczos3 (the shipped default), linear,
sincgauss, and astropy's reproject_exact as an external reference.

Cases:
  - shifted: same pixel scale, sub-pixel shift only (check_interp_
    benchmark.py's own "aligned" case).
  - galactic: a genuine ICRS -> galactic FRAME conversion (not just a
    pixel-grid rotation -- check_interp_benchmark.py's own "rotated"
    case rotates the PC matrix within one system; this changes system).
    The ~1.4 mas AST/astropy ICRS<->galactic gap (docs/dev/
    AST_VS_ASTROPY_FRAMES.md) is used only to place the output's own
    CRVAL at the right galactic point for the test geometry -- 5-6
    orders of magnitude below anything this script measures, so it
    does not need SOFA-level care the way check_target_projection.py's
    own value comparisons do.
  - wide_field: output and input share one pixel scale (10 arcmin),
    so no resolution change is in play, but the field spans +/-20 deg
    -- far enough that SIN's own 1/cos(angle from reference pixel)
    growth is substantial at the edges, unlike every other case in
    this file's own small (~1 deg) fields.
  - coarser / finer: check_interp_benchmark.py's own "downsampled"/
    "upsampled" cases. compact's own beam does not clear the 2.0-cell
    sampling floor once coarser halves resolution, so compact/coarser
    is refused by design (see coarser_case's own comment) -- same
    behaviour check_interp_benchmark.py's compact morphology already
    has on its own downsampled case.

For each case/source/method: two different, established measurements,
not one scheme failing where another succeeds -- the distinction
geophysical regridding tools (ESMF/xESMF/CDO) draw between
"interpolation" and "conservative remapping" (bench.rms_max_error's
own module docstring in check_interp_benchmark.py has the full
version; see also docs/user/APP_REFERENCE.md's "Choosing interp="):

  - rms_error: INTERPOLATION accuracy -- how close each output pixel's
    own value is to the true analytic value AT THAT PIXEL'S OWN SKY
    COORDINATE, never a resampled grid, never an area average. What a
    Shannon-Whittaker-style kernel (lanczos3, sinc, etc.) is built to
    get right.
  - flux_ratio: CONSERVATIVE REMAPPING accuracy -- how close the total
    flux (value x pixel solid angle, summed over the field) is to the
    source's true analytic integral. What interp=exact is built to get
    right, by treating each input pixel as a piecewise-constant cell
    and computing exact geometric overlap -- a provable identity, not
    an approximation.

A scheme answering the first question is under no obligation to also
win the second: interpolation's negative kernel side lobes (what lets
it resolve fine structure precisely at a point) can leave a small,
same-signed residual once many output pixels are summed into a flux --
structural, not a bug. compact's own rms_error is the cleaner read of
the two sources (less curvature across one pixel than extended's, so
less of a gap between "the point value" and "the area average" to
begin with) -- see docs/user/APP_REFERENCE.md for the fuller
discussion and measured numbers across both.

Run time (interp=exact split into building the overlap weights and
applying them, parsed from reproject_cubes' own per-file stage-timing
log, log_level=info's default), and, for exact, the overlap weights'
own memory (parsed from its own printed estimate/actual-size lines).
"""
import platform
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_exact as astropy_reproject_exact
import astropy.units as u

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "tests"))
import check_interp_benchmark as bench  # noqa: E402

OUT_DIR = REPO / "tests" / "output" / "interp_cost"
METHODS = ("exact", "lanczos3", "linear", "sincgauss")
WIDE_CDELT = 10.0 / 60.0  # 10 arcmin/pixel
WIDE_NAXIS = 241  # +/- 20 deg at WIDE_CDELT


def parse_stage_seconds(stdout, stage):
    m = re.search(rf"\[timing\].*?\b{re.escape(stage)}\s+([0-9.]+)\s+([0-9.]+)", stdout)
    return float(m.group(1)) if m else None


def parse_weights_mib(stdout):
    m = re.search(r"overlap weights built, ([0-9.]+) MiB", stdout)
    return float(m.group(1)) if m else None


def shifted_case(source="extended"):
    in_header = bench.build_input_header("icrs", "SIN")
    out_header = bench.build_output_header("aligned", "icrs", "SIN")
    value_fn, _, components = bench.MORPHOLOGIES[source]
    return "shifted", in_header, out_header, value_fn, components, bench.ground_truth_on_output


def ground_truth_cross_system(value_fn, in_header, out_header, out_frame):
    """Like bench.ground_truth_on_output, but for a genuine system
    change: out_header's own world coordinates are in out_frame (e.g.
    "galactic"), not in_header's frame, so they are converted via
    astropy before being handed to in_header's own WCS. The ~1.4 mas
    AST/astropy ICRS<->galactic gap (docs/user/APP_REFERENCE.md) is
    5-6 orders of magnitude below anything this script measures."""
    wcs_in = WCS(in_header)
    wcs_out = WCS(out_header)
    nx_out, ny_out = out_header["NAXIS1"], out_header["NAXIS2"]
    ox, oy = np.meshgrid(np.arange(1, nx_out + 1, dtype=np.float64),
                          np.arange(1, ny_out + 1, dtype=np.float64))
    lon, lat = wcs_out.wcs_pix2world(ox, oy, 1)
    sky = SkyCoord(lon * u.deg, lat * u.deg, frame=out_frame).icrs
    fx, fy = wcs_in.wcs_world2pix(sky.ra.deg, sky.dec.deg, 1)
    return value_fn(fx, fy)


def galactic_case(source="extended"):
    in_header = bench.build_input_header("icrs", "SIN")
    out_header = bench.build_output_header("aligned", "galactic", "SIN")
    sky = SkyCoord(bench.CRVAL1 * u.deg, bench.CRVAL2 * u.deg, frame="icrs").galactic
    out_header["CRVAL1"] = sky.l.deg
    out_header["CRVAL2"] = sky.b.deg
    value_fn, _, components = bench.MORPHOLOGIES[source]
    truth_fn = lambda vf, ih, oh: ground_truth_cross_system(vf, ih, oh, "galactic")  # noqa: E731
    return "galactic", in_header, out_header, value_fn, components, truth_fn


def wide_field_case(source="extended"):
    in_header = fits.Header()
    in_header["SIMPLE"] = True; in_header["BITPIX"] = -32; in_header["NAXIS"] = 2
    in_header["NAXIS1"] = WIDE_NAXIS; in_header["NAXIS2"] = WIDE_NAXIS
    in_header["CTYPE1"] = "RA---SIN"; in_header["CTYPE2"] = "DEC--SIN"
    in_header["RADESYS"] = "ICRS"
    in_header["CRVAL1"], in_header["CRVAL2"] = bench.CRVAL1, bench.CRVAL2
    in_header["CRPIX1"] = in_header["CRPIX2"] = (WIDE_NAXIS + 1) / 2.0
    in_header["CDELT1"] = -WIDE_CDELT; in_header["CDELT2"] = WIDE_CDELT

    out_header = in_header.copy()
    # Sub-pixel shift only, same pixel scale, same wide extent -- the
    # one variable isolated here is distance from the reference pixel,
    # not resolution or system.
    out_header["CRVAL1"] = bench.CRVAL1 + 0.5 * WIDE_CDELT
    out_header["CRVAL2"] = bench.CRVAL2 + 0.5 * WIDE_CDELT

    centre_px = (WIDE_NAXIS + 1) / 2.0
    # One component at the reference pixel, one near the +/-20 deg edge
    # -- the edge one is where SIN's own growth away from the
    # reference pixel should show up as larger interpolation error.
    # Sigma matches bench.EXTENDED_COMPONENTS/COMPACT_COMPONENTS' own
    # convention (4.0px / 0.9,0.8px) so this case's two sources are the
    # same morphologies as every other case here, not a third one.
    sigma_by_source = {"extended": (4.0, 4.0), "compact": (0.9, 0.8)}
    sx0, sx1 = sigma_by_source[source]
    components = [
        (1.0, centre_px, centre_px, sx0, sx0),
        (0.6, centre_px + 0.42 * WIDE_NAXIS, centre_px, sx1, sx1),
    ]
    value_fn = lambda fx, fy: bench.gaussian_sum(fx, fy, components)  # noqa: E731
    return "wide_field", in_header, out_header, value_fn, components, bench.ground_truth_on_output


def coarser_case(source="extended"):
    in_header = bench.build_input_header("icrs", "SIN")
    out_header = bench.build_output_header("downsampled", "icrs", "SIN")
    value_fn, _, components = bench.MORPHOLOGIES[source]
    # Genuine 2x coarsening needs a declared beam (UR-21): the narrowest
    # component's own sigma, same convention check_interp_benchmark.py
    # uses for its own beam-bearing morphologies. compact's own
    # narrowest sigma does not clear the 2.0-cell sampling floor on a
    # 2x coarser grid (same as check_interp_benchmark.py's own
    # downsampled case) -- refused by design, not a bug; its rows here
    # will be empty except for the astropy baseline, which bypasses
    # this tool's own beam-sampling gate entirely.
    sigma_px = min(min(c[3], c[4]) for c in components)
    fwhm_deg = sigma_px * 2.0 * np.sqrt(2.0 * np.log(2.0)) * bench.CDELT_IN
    in_header["BMAJ"] = fwhm_deg; in_header["BMIN"] = fwhm_deg; in_header["BPA"] = 0.0
    return "coarser", in_header, out_header, value_fn, components, bench.ground_truth_on_output


def finer_case(source="extended"):
    in_header = bench.build_input_header("icrs", "SIN")
    out_header = bench.build_output_header("upsampled", "icrs", "SIN")
    value_fn, _, components = bench.MORPHOLOGIES[source]
    return "finer", in_header, out_header, value_fn, components, bench.ground_truth_on_output


def analytic_flux_deg2(components, cdelt_deg):
    return sum(amp * 2 * np.pi * sx * sy for amp, _, _, sx, sy in components) * cdelt_deg ** 2


def run_case(name, source, in_header, out_header, value_fn, components, cdelt_in_deg, truth_fn):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    nx_in, ny_in = in_header["NAXIS1"], in_header["NAXIS2"]
    nx_out, ny_out = out_header["NAXIS1"], out_header["NAXIS2"]
    input_data = bench.sample_onto_grid(value_fn, nx_in, ny_in)
    infile = OUT_DIR / f"{name}_{source}_in.fits"
    reffile = OUT_DIR / f"{name}_{source}_ref.fits"
    fits.PrimaryHDU(data=input_data.astype(np.float32), header=in_header).writeto(
        infile, overwrite=True)
    fits.PrimaryHDU(data=np.zeros((ny_out, nx_out), dtype=np.float32),
                     header=out_header).writeto(reffile, overwrite=True)

    truth = truth_fn(value_fn, in_header, out_header)
    valid_truth = np.isfinite(truth)
    analytic_flux = analytic_flux_deg2(components, cdelt_in_deg)

    rows = []
    for method in METHODS:
        params = bench.SCHEME_PARAMS.get(method)
        try:
            data, elapsed, stdout = bench.run_reproject_cubes(infile, reffile, method, params)
        except RuntimeError as exc:
            print(f"WARNING: {name}/{method} failed: {exc}", file=sys.stderr)
            continue
        valid = valid_truth & np.isfinite(data)
        rms, worst = bench.rms_max_error(data, truth, valid)
        flux = bench.total_flux(data, out_header, valid)
        row = {
            "case": name, "source": source, "method": method, "rms_error": rms,
            "max_error": worst, "flux_ratio": flux / analytic_flux, "seconds": elapsed,
            "weights_build_seconds": parse_stage_seconds(stdout, "exact_weights"),
            "weights_mib": parse_weights_mib(stdout),
        }
        rows.append(row)

    t0 = time.perf_counter()
    e_data, _ = astropy_reproject_exact((input_data, in_header), WCS(out_header).celestial,
                                         shape_out=(ny_out, nx_out))
    e_elapsed = time.perf_counter() - t0
    e_valid = valid_truth & np.isfinite(e_data)
    e_rms, e_worst = bench.rms_max_error(e_data, truth, e_valid)
    e_flux = bench.total_flux(e_data, out_header, e_valid)
    rows.append({
        "case": name, "source": source, "method": "astropy:reproject_exact",
        "rms_error": e_rms, "max_error": e_worst, "flux_ratio": e_flux / analytic_flux,
        "seconds": e_elapsed, "weights_build_seconds": None, "weights_mib": None,
    })

    infile.unlink(missing_ok=True)
    reffile.unlink(missing_ok=True)
    output_of = infile.with_name(infile.stem + "_REPROJ.FITS")
    output_of.unlink(missing_ok=True)
    return rows


def print_table(rows):
    header = (f"{'case':<10} {'source':<9} {'method':<24} {'rms_error':>11} {'max_error':>11} "
              f"{'flux_ratio':>10} {'seconds':>8} {'weights_build':>13} {'weights_MiB':>11}")
    print(header)
    print("-" * len(header))
    for r in rows:
        wb = f"{r['weights_build_seconds']:.3f}" if r["weights_build_seconds"] is not None else "-"
        wm = f"{r['weights_mib']:.2f}" if r["weights_mib"] is not None else "-"
        print(f"{r['case']:<10} {r['source']:<9} {r['method']:<24} {r['rms_error']:>11.3e} "
              f"{r['max_error']:>11.3e} {r['flux_ratio']:>10.6f} {r['seconds']:>8.3f} "
              f"{wb:>13} {wm:>11}")
    print("\n(rms_error measures INTERPOLATION accuracy -- value vs. the true analytic"
          " value AT THAT PIXEL'S OWN COORDINATE; flux_ratio measures CONSERVATIVE"
          " REMAPPING accuracy -- total integrated flux vs. the true analytic integral."
          " Different questions, not the same one with different winners -- see this"
          " file's own module docstring or docs/user/APP_REFERENCE.md's \"Choosing"
          " interp=\" for why a scheme can win one and lose the other. compact's"
          " rms_error is the cleaner interpolation-accuracy read of the two sources:"
          " extended's wider, more-curved source means more gap between \"the value at"
          " a point\" and \"the average over that point's pixel\" even for a"
          " hypothetically perfect reconstruction -- compare the two sources' rms_error"
          " for the same case/method to see it directly.)")


def main():
    if not bench.BIN.exists():
        print(f"ERROR: {bench.BIN} not built -- run: make GPU=0 OMP=1 reproject_cubes")
        return 1

    case_builders = [shifted_case, galactic_case, wide_field_case, coarser_case, finer_case]
    cdelt_by_case = {
        "shifted": bench.CDELT_IN, "galactic": bench.CDELT_IN, "wide_field": WIDE_CDELT,
        "coarser": bench.CDELT_IN, "finer": bench.CDELT_IN,
    }
    all_rows = []
    for build in case_builders:
        for source in ("extended", "compact"):
            name, in_header, out_header, value_fn, components, truth_fn = build(source)
            all_rows.extend(run_case(name, source, in_header, out_header, value_fn,
                                      components, cdelt_by_case[name], truth_fn))

    print_table(all_rows)
    print(f"\nMachine: {platform.platform()}, {platform.processor() or platform.machine()}")
    nthreads = subprocess.run(["nproc"], capture_output=True, text=True).stdout.strip()
    print(f"Threads available (nproc): {nthreads}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
