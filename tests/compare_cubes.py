#!/usr/bin/env python3
"""
Compare two RM-synthesis output FITS cubes.

Modes:
  --exact   : assert bit-identical data arrays (CPU serial vs CPU OpenMP)
  --rtol N  : assert max relative difference ≤ N (GPU vs CPU reference)
"""
import sys
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
        i += 1

    a = load(file_a)
    b = load(file_b)

    if a.shape != b.shape:
        print(f"[FAIL] Shape mismatch: {a.shape} vs {b.shape}")
        sys.exit(1)

    if mode == "exact":
        a32 = a.astype(np.float32)
        b32 = b.astype(np.float32)
        if np.array_equal(a32, b32):
            print("[OK] Cubes are bit-identical (float32)")
            sys.exit(0)
        else:
            n_diff = int(np.sum(a32 != b32))
            max_diff = float(np.max(np.abs(a32 - b32)))
            print(f"[FAIL] {n_diff} elements differ; max |diff| = {max_diff:.3e}")
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
