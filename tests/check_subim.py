#!/usr/bin/env python3
"""
Check rm_synthesis subimage output against the full-resolution run.

  check_subim.py cube --full FULL.fits --sub SUB.fits --ra BLC,TRC,INC --dec BLC,TRC,INC
      SUB must equal FULL[:, dec_blc-1:dec_trc:dec_inc, ra_blc-1:ra_trc:ra_inc]
      bit for bit (NaN positions included), and its header must follow
          CRPIX_out = (CRPIX_full - blc) / inc + 1
          CDELT_out = CDELT_full * inc
      for both sky axes. BLC and TRC are 1-based and inclusive, as in the cfg.

  check_subim.py chan --freqfile sampled_freq.txt --blc B --trc T --inc I --f0 HZ --df HZ
      The sampled-frequency list must hold channels B, B+I, ... up to T
      (1-based), i.e. f0 + (channel-1)*df, and nothing else.

Exit status 0 if every check holds, 1 otherwise.
"""
import argparse
import sys

import numpy as np


def triple(text):
    blc, trc, inc = (int(v) for v in text.split(","))
    return blc, trc, inc


def check_cube(args):
    from astropy.io import fits

    full = fits.open(args.full)[0]
    sub = fits.open(args.sub)[0]
    a = full.data.squeeze()
    b = sub.data.squeeze()
    ra_blc, ra_trc, ra_inc = args.ra
    dec_blc, dec_trc, dec_inc = args.dec
    expect = a[:, dec_blc - 1:dec_trc:dec_inc, ra_blc - 1:ra_trc:ra_inc]
    ok = True
    if b.shape != expect.shape:
        print(f"[FAIL] shape {b.shape}, expected {expect.shape}")
        ok = False
    elif not np.array_equal(b, expect, equal_nan=True):
        n_diff = int(np.sum(b != expect))
        print(f"[FAIL] {n_diff} values are not those of the full-resolution cube")
        ok = False
    for axis, (blc, inc) in ((1, (ra_blc, ra_inc)), (2, (dec_blc, dec_inc))):
        crpix_full = full.header[f"CRPIX{axis}"]
        cdelt_full = full.header[f"CDELT{axis}"]
        crpix_want = (crpix_full - blc) / inc + 1
        cdelt_want = cdelt_full * inc
        crpix_got = sub.header[f"CRPIX{axis}"]
        cdelt_got = sub.header[f"CDELT{axis}"]
        if abs(crpix_got - crpix_want) > 1e-6:
            print(f"[FAIL] CRPIX{axis} = {crpix_got}, expected {crpix_want}")
            ok = False
        if abs(cdelt_got - cdelt_want) > 1e-9 * abs(cdelt_want):
            print(f"[FAIL] CDELT{axis} = {cdelt_got}, expected {cdelt_want}")
            ok = False
    if ok:
        print(f"[OK] {b.shape} cube equals the full-resolution slice; "
              f"CRPIX and CDELT follow blc and inc")
    return ok


def check_chan(args):
    freqs = []
    with open(args.freqfile) as handle:
        for line in handle:
            if line.strip() and not line.lstrip().startswith("#"):
                freqs.append(float(line.split()[0]))
    want = [args.f0 + (ch - 1) * args.df for ch in range(args.blc, args.trc + 1, args.inc)]
    if len(freqs) != len(want):
        print(f"[FAIL] {len(freqs)} channels listed, expected {len(want)}")
        return False
    if not np.allclose(freqs, want, rtol=0, atol=1.0):
        print(f"[FAIL] frequencies start {freqs[0]:.0f}, end {freqs[-1]:.0f}; "
              f"expected {want[0]:.0f} to {want[-1]:.0f}")
        return False
    print(f"[OK] {len(freqs)} channels, {freqs[0]:.0f} to {freqs[-1]:.0f} Hz")
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="mode", required=True)
    cube = sub.add_parser("cube")
    cube.add_argument("--full", required=True)
    cube.add_argument("--sub", required=True)
    cube.add_argument("--ra", type=triple, required=True)
    cube.add_argument("--dec", type=triple, required=True)
    chan = sub.add_parser("chan")
    chan.add_argument("--freqfile", required=True)
    chan.add_argument("--blc", type=int, required=True)
    chan.add_argument("--trc", type=int, required=True)
    chan.add_argument("--inc", type=int, required=True)
    chan.add_argument("--f0", type=float, required=True)
    chan.add_argument("--df", type=float, required=True)
    args = parser.parse_args()
    ok = check_cube(args) if args.mode == "cube" else check_chan(args)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
