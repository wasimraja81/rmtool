#!/usr/bin/env python3
"""Check RM-synthesis AMP cube has peaks at expected RM values."""
import sys, json
from pathlib import Path
import numpy as np
from astropy.io import fits

def main():
    if len(sys.argv) < 3:
        print("Usage: check_rm_peak.py <amp_cube.fits> <truth.json>")
        sys.exit(1)
    amp_path, truth_path = Path(sys.argv[1]), Path(sys.argv[2])
    truth = json.loads(truth_path.read_text())
    with fits.open(amp_path) as hdul:
        data = hdul[0].data.squeeze()   # (nrm, ny, nx)
        hdr  = hdul[0].header
    crval3 = float(hdr["CRVAL3"])
    cdelt3 = float(hdr["CDELT3"])
    nrm    = data.shape[0]
    rm_axis = crval3 + np.arange(nrm) * cdelt3
    tol = 2.0 * abs(cdelt3)           # allow ±2 RM cells
    ok = True
    for src in truth["sources"]:
        x, y   = src["x"], src["y"]
        rm_exp = src["rm"]
        spectrum = data[:, y, x]
        peak_i   = int(np.argmax(spectrum))
        rm_found = float(rm_axis[peak_i])
        err      = abs(rm_found - rm_exp)
        flag = "OK" if err <= tol else "FAIL"
        if flag == "FAIL":
            ok = False
        print(f"[{flag}] {src['name']}: expected RM={rm_exp:+.1f}, "
              f"found RM={rm_found:+.2f} (err={err:.2f}, tol={tol:.2f})")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
