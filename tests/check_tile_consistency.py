#!/usr/bin/env python3
"""Compare two rmclean_cubes runs (e.g. default tiling vs a forced small
tile_ra/tile_dec) for numerical agreement.

NOT byte-identical, deliberately: clean_complex's tiered peak refinement
(T3c) makes a floating-point threshold decision (fast-path vs escalate-to-
search) every CLEAN iteration, and gfortran's -O3 -march=native
auto-vectorization can reassociate that comparison's floating-point sum
differently depending on the runtime memory alignment of clean_complex's
own stack-allocated arguments -- which genuinely differs between tile
sizes (different automatic-array footprints in the enclosing call frames),
even though every value going INTO clean_complex is bit-identical (verified
directly, planning/RMCLEAN_INTEGRATION_PLAN.md ticket T4a). Confirmed this
is a real, pre-existing floating-point reassociation effect (not a tiling
logic bug): the SAME two tile configurations produce byte-identical output
when built at -O0 (no auto-vectorization). PHA at near-zero AMP is
mathematically ill-conditioned (atan2 of a near-zero complex number is
hypersensitive to tiny perturbations) and is therefore only compared where
AMP clears a floor.

DEFAULT amp_atol=5e-4 is anchored to the fixture's own noise floor, not
a round number: tests/data/truth.json's own noise_sigma=0.01 (per-
channel Q/U), 200 channels -> dirty/restored AMP's own median pixel
value is ~8e-4-1.5e-3 (measured directly, not assumed). The observed
worst-case diff across all 6 output cubes and every tile/thread
configuration tested is 1.7e-4 (CLEAN.AMP) -- 5e-4 gives ~3x margin
over that observed worst case while staying BELOW the noise floor, so
a genuine future regression corrupting noise-level pixels by more than
the fixture's own measurement noise would still be caught (the
original 1e-2, anchored to the CLEAN stopping threshold rather than
the noise floor, was ~12x looser than the noise floor itself and would
not have caught that).
"""
import sys
from pathlib import Path
import numpy as np
from astropy.io import fits


def compare_pair(base, other, suffix, amp_atol, pha_atol, amp_floor):
    amp_a = fits.getdata(f"{base}.{suffix.replace('PHA','AMP')}.RMCUBE.FITS")
    if suffix.endswith("AMP"):
        a = fits.getdata(f"{base}.{suffix}.RMCUBE.FITS")
        b = fits.getdata(f"{other}.{suffix}.RMCUBE.FITS")
        nan_a, nan_b = np.isnan(a), np.isnan(b)
        if not np.array_equal(nan_a, nan_b):
            return False, "NaN mask mismatch"
        diff = np.abs(a[~nan_a] - b[~nan_a])
        worst = float(np.max(diff)) if diff.size else 0.0
        return worst <= amp_atol, f"max|diff|={worst:.3e} (atol={amp_atol:.1e})"
    else:
        a = fits.getdata(f"{base}.{suffix}.RMCUBE.FITS")
        b = fits.getdata(f"{other}.{suffix}.RMCUBE.FITS")
        mask = (~np.isnan(a)) & (~np.isnan(b)) & (amp_a > amp_floor)
        wrapped = np.abs(np.angle(np.exp(1j * (a[mask] - b[mask]))))
        worst = float(np.max(wrapped)) if wrapped.size else 0.0
        return worst <= pha_atol, f"max|dphase| (AMP>{amp_floor})={worst:.3e} (atol={pha_atol:.1e})"


def main():
    if len(sys.argv) < 3:
        print("Usage: check_tile_consistency.py <base_outfile> <other_outfile> "
              "[amp_atol=5e-4] [pha_atol=0.05] [amp_floor=1e-2]")
        sys.exit(1)
    base, other = sys.argv[1], sys.argv[2]
    amp_atol = float(sys.argv[3]) if len(sys.argv) > 3 else 5e-4
    pha_atol = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05
    amp_floor = float(sys.argv[5]) if len(sys.argv) > 5 else 1e-2

    ok = True
    for suffix in ["CLEAN.AMP", "CLEAN.PHA", "RESID.AMP", "RESID.PHA",
                   "RESTORED.AMP", "RESTORED.PHA"]:
        passed, detail = compare_pair(base, other, suffix, amp_atol, pha_atol, amp_floor)
        flag = "OK" if passed else "FAIL"
        if not passed:
            ok = False
        print(f"[{flag}] {suffix}: {detail}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
