#!/usr/bin/env python3
"""
Compare two RM-synthesis output FITS cubes.

Modes:
  --exact   : assert bit-identical data arrays (CPU serial vs CPU OpenMP)
  --rtol N  : assert max relative difference ≤ N (GPU vs CPU reference)
  --peak-rtol N : assert max |a-b| / (peak of |a| along the RM axis, per
              pixel) ≤ N, and identical NaN patterns. Use this when the cube
              holds amplitudes that cancel to rounding noise (rem_mean>0):
              a pointwise relative difference is dominated by those values.
"""
import sys
import warnings
from pathlib import Path
import numpy as np
from astropy.io import fits

def load(path: str) -> np.ndarray:
    with fits.open(path) as hdul:
        return hdul[0].data.squeeze().astype(np.float64)

def main():
    args = sys.argv[1:]
    if len(args) < 2:
        print("Usage: compare_cubes.py <file_a> <file_b> [--exact | --rtol <val>]")
        sys.exit(1)

    file_a, file_b = args[0], args[1]
    mode = "exact"
    rtol = 1e-4
    i = 2
    while i < len(args):
        if args[i] == "--exact":
            mode = "exact"
        elif args[i] == "--rtol" and i + 1 < len(args):
            mode = "rtol"
            rtol = float(args[i + 1])
            i += 1
        elif args[i] == "--peak-rtol" and i + 1 < len(args):
            mode = "peak"
            rtol = float(args[i + 1])
            i += 1
        i += 1

    a = load(file_a)
    b = load(file_b)

    if a.shape != b.shape:
        print(f"[FAIL] Shape mismatch: {a.shape} vs {b.shape}")
        sys.exit(1)

    if mode == "exact":
        a32 = a.astype(np.float32)
        b32 = b.astype(np.float32)
        # equal_nan: two cubes blanked identically ARE bit-identical --
        # without it, NaN != NaN makes np.array_equal report a false
        # difference for any cube carrying masked pixels.
        if np.array_equal(a32, b32, equal_nan=True):
            print("[OK] Cubes are bit-identical (float32)")
            sys.exit(0)
        else:
            differ = (a32 != b32) & ~(np.isnan(a32) & np.isnan(b32))
            n_diff = int(np.sum(differ))
            max_diff = float(np.nanmax(np.abs(a32 - b32))) if n_diff else 0.0
            print(f"[FAIL] {n_diff} elements differ; max |diff| = {max_diff:.3e}")
            sys.exit(1)
    elif mode == "peak":
        if not np.array_equal(np.isnan(a), np.isnan(b)):
            print("[FAIL] NaN patterns differ")
            sys.exit(1)
        with warnings.catch_warnings():
            # a fully blanked pixel has no peak; its NaNs were matched above
            warnings.simplefilter("ignore", RuntimeWarning)
            peak = np.nanmax(np.abs(a), axis=0, keepdims=True) if a.ndim == 3 else np.nanmax(np.abs(a))
            ratio = np.abs(a - b) / np.where(peak > 0, peak, 1.0)
            worst = float(np.nanmax(ratio)) if not np.all(np.isnan(ratio)) else 0.0
        if worst <= rtol:
            print(f"[OK] max |diff| / pixel peak = {worst:.3e} ≤ peak-rtol={rtol:.3e}")
            sys.exit(0)
        else:
            print(f"[FAIL] max |diff| / pixel peak = {worst:.3e} > peak-rtol={rtol:.3e}")
            sys.exit(1)
    else:
        denom = np.where(np.abs(a) > 1e-30, np.abs(a), 1e-30)
        rel   = np.max(np.abs(a - b) / denom)
        if rel <= rtol:
            print(f"[OK] max relative diff = {rel:.3e} ≤ rtol={rtol:.3e}")
            sys.exit(0)
        else:
            print(f"[FAIL] max relative diff = {rel:.3e} > rtol={rtol:.3e}")
            sys.exit(1)

if __name__ == "__main__":
    main()
