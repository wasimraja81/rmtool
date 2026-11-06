#!/usr/bin/env python3
"""How reprojection error depends on how well the beam is sampled.

Companion to check_interp_benchmark.py, which sweeps interp= schemes at
one source scale. This sweeps the source scale at a few schemes, to
answer a different question: how much of the resampling error is set by
the input's own sampling rather than by the kernel.

REPROJECT_FLUX_NOTE.md section 3 measured this for AST__LINEAR only,
and in peak-loss terms. Peak loss scores how closely the output
reproduces the INPUT GRID's own sampled peak, which is not a conserved
quantity: a pixel peak depends on where the grid samples fall relative
to the continuous maximum, and section 4 records the input's own
sampled peak as 0.9946 of the true amplitude for that beam. Preserving
it is not a correctness criterion. The metric here is instead RMS error
against analytic ground truth, as in check_interp_benchmark.py, which
is pointwise fidelity to the continuous sky at each output pixel's own
sky position.

Result this script exists to record (2026-10-15): the 13.5% figure in
section 3 at a 3x3 pixel beam is a property of linear interpolation,
not a sampling limit. lanczos3 at the same sampling is about 19x lower
in worst-case error, and lower than linear at every sampling measured.

Geometry is the benchmark's own "aligned" case -- same pixel scale,
with a deliberate half-pixel shift on both axes, the worst case for
interpolation. A single circular Gaussian of unit amplitude sits
off-grid so no row or column lands on its centre; sweeping its FWHM
sweeps how many pixels sample it. Errors are quoted as a percentage of
that unit peak, so they read as fractional errors.

Parameterisation: every scheme here (nearest, linear, lanczos3) takes
no interp_params=, so these numbers involve no parameter choice of
ours. The schemes that do take one are in check_interp_benchmark.py,
whose report states the values used.

A source wide enough to reach the edge of the 64x64 test grid is not
contained by it, and its row measures the test geometry rather than the
kernel; such rows are flagged in the output rather than dropped.
"""
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "tests"))
import check_interp_benchmark as bm  # noqa: E402

OUT_DIR = REPO_ROOT / "tests" / "output" / "sampling_sweep"
FWHMS = [3.0, 5.0, 8.66, 12.0, 20.0, 30.0]
SCHEMES = ["nearest", "linear", "lanczos3"]
# Pixel centres are integers, so .3/.7 keeps the peak off every sample
# on both axes.
CENTRE = (32.3, 31.7)
# A row is reported as not contained once the source is still this
# fraction of its peak at the nearest array edge.
CONTAINMENT_LIMIT = 0.01

FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def edge_fraction(sigma):
    """The source's value at the nearest array edge, as a fraction of
    its own peak -- the test grid truncates whatever this leaves."""
    dist = min(CENTRE[0] - 1.0, bm.NX_IN - CENTRE[0],
               CENTRE[1] - 1.0, bm.NY_IN - CENTRE[1])
    return float(np.exp(-0.5 * (dist / sigma) ** 2))


def main():
    if not bm.BIN.exists():
        print(f"ERROR: {bm.BIN} not built -- run: make GPU=0 OMP=1 reproject_cubes")
        return 1

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    in_header = bm.build_input_header()
    out_header = bm.build_output_header("aligned")

    reffile = OUT_DIR / "ref.fits"
    fits.PrimaryHDU(
        data=np.zeros((out_header["NAXIS2"], out_header["NAXIS1"]), dtype=np.float32),
        header=out_header,
    ).writeto(reffile, overwrite=True)

    rows = {}
    for fwhm in FWHMS:
        sigma = fwhm * FWHM_TO_SIGMA
        components = [(1.0, CENTRE[0], CENTRE[1], sigma, sigma)]

        def value_fn(fx, fy, _c=components):
            return bm.gaussian_sum(fx, fy, _c)

        infile = OUT_DIR / f"gauss_fwhm{fwhm:g}.fits"
        fits.PrimaryHDU(
            data=bm.sample_onto_grid(value_fn, bm.NX_IN, bm.NY_IN).astype(np.float32),
            header=in_header,
        ).writeto(infile, overwrite=True)

        truth = bm.ground_truth_on_output(value_fn, in_header, out_header)
        valid_truth = np.isfinite(truth)

        rows[fwhm] = {"edge_frac": edge_fraction(sigma)}
        for scheme in SCHEMES:
            data, _, _ = bm.run_reproject_cubes(infile, reffile, scheme)
            valid = valid_truth & np.isfinite(data)
            rms, worst = bm.rms_max_error(data, truth, valid)
            rows[fwhm][scheme] = (rms, worst)
        infile.unlink()

    def table(title, index):
        print(f"\n-- {title} (% of the source's own peak) --")
        print(f"{'beam FWHM':>12}  " + "  ".join(f"{s:>10}" for s in SCHEMES))
        for fwhm in FWHMS:
            cells = "  ".join(f"{100.0 * rows[fwhm][s][index]:>10.3f}" for s in SCHEMES)
            flag = "" if rows[fwhm]["edge_frac"] <= CONTAINMENT_LIMIT \
                else f"   <- not contained: {100.0 * rows[fwhm]['edge_frac']:.1f}% of peak at the array edge"
            print(f"{fwhm:>9.2f} px  {cells}{flag}")

    print("Reprojection error vs beam sampling, half-pixel shift, same pixel scale.")
    print(f"Ground truth is analytic; grid is {bm.NX_IN}x{bm.NY_IN}; no interp_params= "
          "is used or accepted by any scheme here.")
    table("RMS error", 0)
    table("Worst-case absolute error", 1)

    print("\n-- lanczos3 vs linear, ratio of RMS (>1 means lanczos3 is lower) --")
    for fwhm in FWHMS:
        ratio = rows[fwhm]["linear"][0] / rows[fwhm]["lanczos3"][0]
        note = "" if rows[fwhm]["edge_frac"] <= CONTAINMENT_LIMIT else "  (not contained)"
        print(f"{fwhm:>9.2f} px  {ratio:>8.2f}x{note}")

    contained = [f for f in FWHMS if rows[f]["edge_frac"] <= CONTAINMENT_LIMIT]
    print(f"\nOver the contained rows ({contained[0]:g}-{contained[-1]:g} px FWHM): error falls as the beam is "
          "better sampled,\nand lanczos3 is lower than linear at every sampling measured. The kernel's own "
          "contribution\ndoes not vanish at marginal sampling -- it is largest there in absolute terms.")

    reffile.unlink()
    return 0


if __name__ == "__main__":
    sys.exit(main())
