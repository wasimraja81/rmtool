#!/usr/bin/env python3
"""
UR-8 (docs/dev/BRANCH_STATUS_update-reproject.md): value-level
correctness for every reproject_cubes interp= scheme except nearest
(already covered by check_nearest_ramp in check_reproject_correctness.
py) -- promotes sinc/sincsinc/sinccos/sincgauss/blockave/gauss/somb/
sombcos/lanczos2/lanczos3 and linear itself (previously checked only
as "bit-identical to the unnamed default", never under its own name)
from the "reaches AST, differs from linear" sanity level (tests/
run_tests.sh sections 63-66) to real correctness, against the same
analytic ground truth check_interp_benchmark.py already built (not
duplicated -- imported from it).

Tolerance per scheme is a regression guard, not an absolute spec,
same style as check_reproject_correctness.py's own NaN-footprint
check ("measured baseline X, allow <= Y"): 2x the mean RMS error
check_interp_benchmark.py measured for that scheme across the
"extended"/"compact" (realistic-source) morphologies and all 4
geometric cases, on 2026-10-13. The "step_edge" morphology is
deliberately excluded here -- its Gibbs-ringing error is large for
every non-trivial scheme and does not discriminate "is this scheme
implemented correctly" from "ringing near a sharp discontinuity is
inherent to any linear-phase kernel"; UR-12's own benchmark already
covers that comparison. A scheme regressing well past this margin
(e.g. silently falling back to nearest, or a broken kernel formula)
fails; ordinary floating-point/platform noise does not.

Exit 0 if every scheme's measured RMS stays within its own tolerance,
1 otherwise. Prints one PASS/FAIL line per scheme.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_interp_benchmark as bench  # noqa: E402

# Mean RMS measured by check_interp_benchmark.py (extended+compact
# morphologies, all 4 geometric cases), 2026-10-13, x2 safety margin.
TOLERANCE = {
    "linear": 3.1402e-03 * 2,
    "sinc": 6.2860e-03 * 2,
    "sincsinc": 3.2564e-03 * 2,
    "sinccos": 4.0601e-03 * 2,
    "sincgauss": 1.8012e-03 * 2,
    "blockave": 1.3539e-02 * 2,
    "gauss": 7.8667e-03 * 2,
    "somb": 1.8072e-03 * 2,
    "sombcos": 3.5972e-03 * 2,
    "lanczos2": 1.2835e-03 * 2,
    "lanczos3": 5.2296e-04 * 2,
}
CORRECTNESS_MORPHOLOGIES = ("extended", "compact")
CORRECTNESS_CASES = ("aligned", "rotated", "upsampled", "downsampled")
# UR-24 B2: all 4 systems x 5 projections, not just icrs/SIN.
SYSTEMS = tuple(bench.SYSTEMS)
PROJECTIONS = bench.PROJECTIONS

# Each input carries a beam equal to its own narrowest feature (the
# smallest Gaussian component's FWHM), as an observed image would.
# reproject_cubes refuses an output grid whose cells the beam spans
# fewer than 2.0 times (UR-21, docs/dev/BRANCH_STATUS_update-
# reproject.md): 2x downsampling leaves the extended source at about 3.5
# cells, which runs and is checked against truth, and the compact one at
# about 0.9, which aliases and must be refused rather than measured.
FWHM_PER_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))
COMPONENTS = {"extended": bench.EXTENDED_COMPONENTS, "compact": bench.COMPACT_COMPONENTS}


def beam_fwhm_deg(morph_name):
    sigma_px = min(min(c[3], c[4]) for c in COMPONENTS[morph_name])
    return sigma_px * FWHM_PER_SIGMA * bench.CDELT_IN


def expect_refusal(morph_name, case):
    if case != "downsampled":
        return False
    out_cell_deg = bench.CDELT_IN * 2.0
    return beam_fwhm_deg(morph_name) < 2.0 * out_cell_deg


def main():
    if not bench.BIN.exists():
        print(f"ERROR: {bench.BIN} not built -- run: make GPU=0 OMP=1 reproject_cubes")
        return 1

    bench.OUT_DIR.mkdir(parents=True, exist_ok=True)
    ok = True

    for system in SYSTEMS:
        for projection in PROJECTIONS:
            for morph_name in CORRECTNESS_MORPHOLOGIES:
                value_fn, _, _ = bench.MORPHOLOGIES[morph_name]
                in_header = bench.build_input_header(system, projection)
                in_header["BMAJ"] = beam_fwhm_deg(morph_name)
                in_header["BMIN"] = beam_fwhm_deg(morph_name)
                in_header["BPA"] = 0.0
                infile = bench.OUT_DIR / f"{system}_{projection}_{morph_name}_correctness.fits"
                input_data = bench.sample_onto_grid(value_fn, bench.NX_IN, bench.NY_IN)
                fits.PrimaryHDU(data=input_data.astype(np.float32), header=in_header).writeto(
                    infile, overwrite=True)

                for case in CORRECTNESS_CASES:
                    out_header = bench.build_output_header(case, system, projection)
                    reffile = (bench.OUT_DIR /
                               f"{system}_{projection}_{morph_name}_{case}_correctness_ref.fits")
                    nx_out, ny_out = out_header["NAXIS1"], out_header["NAXIS2"]
                    fits.PrimaryHDU(data=np.zeros((ny_out, nx_out), dtype=np.float32),
                                     header=out_header).writeto(reffile, overwrite=True)

                    truth = bench.ground_truth_on_output(value_fn, in_header, out_header)
                    valid_truth_mask = np.isfinite(truth)

                    for scheme, tolerance in TOLERANCE.items():
                        refused = expect_refusal(morph_name, case)
                        try:
                            data, _, _ = bench.run_reproject_cubes(
                                infile, reffile, scheme, bench.SCHEME_PARAMS.get(scheme))
                        except RuntimeError as exc:
                            if refused and "output cells; at least 2.0 are needed" in str(exc):
                                continue
                            print(f"[FAIL] {scheme}: {system}/{projection}/{morph_name}/"
                                  f"{case} run failed: {exc}")
                            ok = False
                            continue
                        if refused:
                            print(f"[FAIL] {scheme}: {system}/{projection}/{morph_name}/"
                                  f"{case} should have been refused"
                                  f" (beam under 2.0 output cells) but ran")
                            ok = False
                            continue
                        valid = valid_truth_mask & np.isfinite(data)
                        rms, _ = bench.rms_max_error(data, truth, valid)
                        if rms > tolerance:
                            print(f"[FAIL] {scheme}: {system}/{projection}/{morph_name}/"
                                  f"{case} RMS error {rms:.4e} exceeds tolerance "
                                  f"{tolerance:.4e}")
                            ok = False

                    reffile.unlink()
                infile.unlink()

    if ok:
        for scheme in TOLERANCE:
            print(f"[OK] {scheme}: within tolerance {TOLERANCE[scheme]:.4e} across "
                  f"{len(SYSTEMS)} systems x {len(PROJECTIONS)} projections x "
                  f"{CORRECTNESS_MORPHOLOGIES} x {CORRECTNESS_CASES}, except "
                  f"compact/downsampled, asserted refused as undersampled")
    else:
        print("[FAIL] one or more schemes exceeded their tolerance -- see above")

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
