#!/usr/bin/env python3
"""Where UR-19's default kernel widths come from.

Three interp= schemes take a width AST cannot choose for itself:
gauss and sincgauss read a Gaussian FWHM in input-grid pixels,
blockave a box half-width. Omitting it used to be a hard error, which
left three of the twelve schemes unreachable unless the caller already
knew a number to pass. ast_interp_mod.f90 now supplies a default, and
this script is where those defaults are measured.

Selection rule: the width with the lowest WORST CASE over beam
samplings and output geometries, not the lowest average. A default
cannot know the caller's geometry, and the per-geometry optima
disagree -- gauss prefers a narrow width on an axis-aligned same-scale
grid and a wider one on a rotated or upsampled grid, by up to 8x. An
earlier pass that swept only the aligned case picked a width that was
8x off on a rotated grid, which is why every geometry is swept here.

Both worst-case RMS and worst-case maximum absolute error are
reported, and the chosen width has to hold on both. RMS alone is run
over the whole frame, so on a mostly-empty field it underweights the
source; the maximum catches a kernel that flattens the source while
leaving the empty regions right.

None of these widths scales with the beam. An earlier design derived
them from BMAJ and the pixel scale, on the argument that a Gaussian
kernel of FWHM K broadens a beam B to sqrt(B^2+K^2) so K should track
B. The measurement found no beam-proportional optimum to derive:
gauss and blockave degrade monotonically as the width grows, so their
best width is the narrowest that still resamples, and sincgauss
plateaus wide. The derivation was dropped.

Metric and geometry otherwise match check_interp_benchmark.py and
check_sampling_sweep.py: RMS against analytic ground truth, a single
off-grid Gaussian of unit amplitude, errors as a percentage of that
peak.
"""
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "tests"))
import check_interp_benchmark as bm  # noqa: E402

OUT_DIR = REPO_ROOT / "tests" / "output" / "interp_width_defaults"
# Contained in the 64x64 grid -- check_sampling_sweep.py's own
# containment test rules out anything wider.
BEAMS = [3.0, 8.66, 20.0]
CASES = ("aligned", "rotated", "upsampled", "downsampled")
# UR-24 B2: all 4 systems x 5 projections, not just icrs/SIN, so the
# "worst case" a default width is chosen against also covers any
# per-projection kernel-response difference, not just per-geometry.
SYSTEMS = tuple(bm.SYSTEMS)
PROJECTIONS = bm.PROJECTIONS
CENTRE = (32.3, 31.7)
FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))

# Widths swept per scheme, and the value ast_interp_mod.f90 ships.
# Keep these two in step: a change to the module's constants without a
# re-run here leaves the module asserting a measurement nobody made.
SWEEPS = [
    ("gauss", [0.5, 0.8, 1.5, 2.0, 3.0], 0.8),
    ("sincgauss", [2.0, 3.0, 4.0, 6.0, 8.0], 8.0),
    ("blockave", [0.5, 0.8, 1.0, 1.25, 2.0], 1.0),
]


def params_for(scheme, width):
    return f"{width:g}" if scheme == "blockave" else f"0,{width:g}"


def main():
    if not bm.BIN.exists():
        print(f"ERROR: {bm.BIN} not built -- run: make GPU=0 OMP=1 reproject_cubes")
        return 1
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rc = 0

    for scheme, widths, shipped in SWEEPS:
        rms_of = {w: [] for w in widths}
        max_of = {w: [] for w in widths}
        ref_rms, ref_max = [], []

        for system in SYSTEMS:
            for projection in PROJECTIONS:
                in_header = bm.build_input_header(system, projection)
                for beam in BEAMS:
                    sigma = beam * FWHM_TO_SIGMA
                    comps = [(1.0, CENTRE[0], CENTRE[1], sigma, sigma)]

                    def value_fn(fx, fy, _c=comps):
                        return bm.gaussian_sum(fx, fy, _c)

                    infile = OUT_DIR / f"{system}_{projection}_beam{beam:g}.fits"
                    # Each source is a beam of this FWHM, written to the header so
                    # reproject_cubes can check its sampling (UR-21): the 3 px beam
                    # on the 2x downsampled grid spans 1.5 cells and is refused,
                    # so the worst case below covers the runs that are allowed.
                    beam_header = in_header.copy()
                    beam_header["BMAJ"] = beam * bm.CDELT_IN
                    beam_header["BMIN"] = beam * bm.CDELT_IN
                    beam_header["BPA"] = 0.0
                    fits.PrimaryHDU(
                        data=bm.sample_onto_grid(value_fn, bm.NX_IN, bm.NY_IN).astype(np.float32),
                        header=beam_header).writeto(infile, overwrite=True)

                    for case in CASES:
                        out_header = bm.build_output_header(case, system, projection)
                        reffile = OUT_DIR / f"{system}_{projection}_ref_{case}.fits"
                        fits.PrimaryHDU(
                            data=np.zeros((out_header["NAXIS2"], out_header["NAXIS1"]),
                                          dtype=np.float32),
                            header=out_header).writeto(reffile, overwrite=True)
                        truth = bm.ground_truth_on_output(value_fn, in_header, out_header)
                        valid_truth = np.isfinite(truth)

                        for w in list(widths) + [None]:
                            scm = scheme if w is not None else "lanczos3"
                            prm = params_for(scheme, w) if w is not None else None
                            try:
                                data, _, _ = bm.run_reproject_cubes(infile, reffile, scm, prm)
                            except RuntimeError:
                                rms = worst = float("nan")
                            else:
                                valid = valid_truth & np.isfinite(data)
                                if float(np.mean(valid)) < 0.5:
                                    rms = worst = float("nan")
                                else:
                                    rms, worst = bm.rms_max_error(data, truth, valid)
                            if w is None:
                                ref_rms.append(rms)
                                ref_max.append(worst)
                            else:
                                rms_of[w].append(rms)
                                max_of[w].append(worst)
                        reffile.unlink()
                    infile.unlink()

        unit = "box half-width" if scheme == "blockave" else "Gaussian FWHM"
        print(f"\n{'='*76}\n{scheme}: {unit} in pixels, over "
              f"{len(BEAMS)} beams x {len(CASES)} geometries x "
              f"{len(SYSTEMS)} systems x {len(PROJECTIONS)} projections\n{'='*76}")
        print(f"{'width':>9}  {'worst RMS':>10}  {'worst max':>10}   (% of peak)")
        for w in widths:
            mark = "  <- shipped" if w == shipped else ""
            print(f"{w:>9.2f}  {100*np.nanmax(rms_of[w]):>10.3f}  "
                  f"{100*np.nanmax(max_of[w]):>10.3f}{mark}")
        print(f"{'lanczos3':>9}  {100*np.nanmax(ref_rms):>10.3f}  "
              f"{100*np.nanmax(ref_max):>10.3f}   (reference, takes no width)")

        # A tie counts as best. blockave's half-width is quantised by
        # AST -- 0.5 through 1.25 give the same output, down to an
        # identical blank footprint and peak -- so requiring the
        # numerically-first argmin would reject an equally good width
        # for no reason. Among tied widths the one to ship is the one
        # with margin from a failure: blockave goes all-NaN at 0.4, so
        # mid-plateau beats sitting next to the cliff.
        best_rms = float(np.nanmin([np.nanmax(rms_of[w]) for w in widths]))
        best_max = float(np.nanmin([np.nanmax(max_of[w]) for w in widths]))
        shipped_rms = float(np.nanmax(rms_of[shipped]))
        shipped_max = float(np.nanmax(max_of[shipped]))
        tied_rms = np.isclose(shipped_rms, best_rms, rtol=1.0e-9, atol=0.0)
        tied_max = np.isclose(shipped_max, best_max, rtol=1.0e-9, atol=0.0)
        print(f"  best worst-case RMS {100*best_rms:.3f}, shipped {100*shipped_rms:.3f}; "
              f"best worst-case maximum {100*best_max:.3f}, shipped {100*shipped_max:.3f}")
        if tied_rms and tied_max:
            print(f"  ast_interp_mod.f90 ships {shipped:g} -- best (or tied best) on both")
        else:
            print(f"  MISMATCH: {shipped:g} is not best or tied-best on both metrics -- "
                  f"re-check ast_interp_mod.f90's default for {scheme}")
            rc = 1

    print("\nShipped widths are in tests/../src/ast_interp_mod.f90 as "
          "default_gauss_fwhm_px,\ndefault_sincgauss_fwhm_px and "
          "default_blockave_halfwidth_px.")
    return rc


if __name__ == "__main__":
    sys.exit(main())
