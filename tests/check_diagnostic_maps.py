#!/usr/bin/env python3
"""T19 Part B: cross-validate the 6 CLEAN diagnostic maps against the
6 real (already-independently-verified) RM cubes from the SAME run,
not just against themselves -- catches bugs a self-consistent-looking
map could still hide (e.g. reading the wrong pixel, an off-by-one in
the global-vs-per-block reset).
"""
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits


def main():
    if len(sys.argv) < 3:
        print("Usage: check_diagnostic_maps.py <outfile_base> <printed_mean_niter>")
        sys.exit(1)
    base, printed_mean_niter = Path(sys.argv[1]), float(sys.argv[2])
    ok = True

    niter = fits.getdata(f"{base}.NITER.MAP.FITS")
    stop_reason = fits.getdata(f"{base}.STOP_REASON.MAP.FITS")
    resid_peak = fits.getdata(f"{base}.RESID_PEAK.MAP.FITS")
    resid_rms = fits.getdata(f"{base}.RESID_RMS.MAP.FITS")
    n_components = fits.getdata(f"{base}.N_COMPONENTS.MAP.FITS")
    comp_rm_spread = fits.getdata(f"{base}.COMP_RM_SPREAD.MAP.FITS")

    if niter.dtype.byteorder in ("=", "<", ">") and niter.dtype.itemsize != 2:
        print(f"FAIL: NITER dtype expected 2-byte int, got {niter.dtype}")
        ok = False
    if stop_reason.dtype.itemsize != 1:
        print(f"FAIL: STOP_REASON dtype expected 1-byte, got {stop_reason.dtype}")
        ok = False
    if resid_peak.dtype.itemsize != 4 or resid_peak.dtype.kind != "f":
        print(f"FAIL: RESID_PEAK dtype expected float32, got {resid_peak.dtype}")
        ok = False

    # Cross-validate NITER against this SAME run's own printed global
    # mean -- independent of how the map itself was populated.
    computed_mean_niter = float(niter[niter > 0].mean())
    if abs(computed_mean_niter - printed_mean_niter) > 0.01:
        print(
            f"FAIL: NITER map mean ({computed_mean_niter:.4f}) != "
            f"printed global mean ({printed_mean_niter:.4f})"
        )
        ok = False

    # Cross-validate RESID_PEAK/N_COMPONENTS against the REAL RESID.AMP/
    # CLEAN.AMP output cubes from this same run, computed independently
    # here rather than trusting the same internal computation twice.
    resid_amp_cube = fits.getdata(f"{base}.RESID.AMP.RMCUBE.FITS")  # (nrm,ny,nx)
    clean_amp_cube = fits.getdata(f"{base}.CLEAN.AMP.RMCUBE.FITS")
    computed_resid_peak = np.nanmax(resid_amp_cube, axis=0)
    computed_n_components = np.sum(clean_amp_cube > 0, axis=0)

    cleaned_mask = niter > 0  # pixels this run actually CLEANed
    if not np.allclose(
        resid_peak[cleaned_mask], computed_resid_peak[cleaned_mask], rtol=1e-5
    ):
        print("FAIL: RESID_PEAK.MAP.FITS does not match max(RESID.AMP) per pixel")
        ok = False
    if not np.array_equal(
        n_components[cleaned_mask], computed_n_components[cleaned_mask]
    ):
        print(
            "FAIL: N_COMPONENTS.MAP.FITS does not match count(CLEAN.AMP>0) per pixel"
        )
        ok = False

    # STOP_REASON must be nonzero (1/2/3) exactly where CLEANed, and
    # exactly 0 elsewhere.
    if not np.array_equal(stop_reason[cleaned_mask] > 0, cleaned_mask[cleaned_mask]):
        print("FAIL: STOP_REASON is 0 (not-cleaned) at a pixel NITER says was cleaned")
        ok = False
    if np.any(stop_reason[~cleaned_mask] != 0):
        print("FAIL: STOP_REASON is nonzero at a pixel NITER says was NOT cleaned")
        ok = False

    # COMP_RM_SPREAD: NaN exactly where N_COMPONENTS<2, finite exactly
    # where N_COMPONENTS>=2 (the ticket's own explicit fill-value rule).
    spread_nan = np.isnan(comp_rm_spread)
    expect_nan = n_components < 2
    if not np.array_equal(spread_nan, expect_nan):
        print("FAIL: COMP_RM_SPREAD NaN pattern doesn't match N_COMPONENTS<2 exactly")
        ok = False
    if np.any(comp_rm_spread[~spread_nan] < 0):
        print("FAIL: COMP_RM_SPREAD has a negative value (should be a spread, >=0)")
        ok = False

    # RESID_RMS: sane range check (non-negative, finite everywhere
    # CLEANed).
    if np.any(resid_rms[cleaned_mask] < 0) or not np.all(
        np.isfinite(resid_rms[cleaned_mask])
    ):
        print("FAIL: RESID_RMS has a negative or non-finite value at a CLEANed pixel")
        ok = False

    # HISTORY cards on STOP_REASON's own header.
    hdr = fits.getheader(f"{base}.STOP_REASON.MAP.FITS")
    history = " ".join(str(h) for h in hdr.get("HISTORY", []))
    if "0=not cleaned" not in history or "3=stopped at auto_nsigma" not in history:
        print("FAIL: STOP_REASON.MAP.FITS header missing its own code-mapping HISTORY cards")
        ok = False

    if ok:
        print(
            f"OK: all 6 diagnostic maps cross-validated against "
            f"RESID.AMP/CLEAN.AMP ({cleaned_mask.sum()} CLEANed pixels, "
            f"{(n_components >= 2).sum()} with >=2 components)"
        )
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
