#!/usr/bin/env python3
"""
Self-contained synthetic Q/U cube generator for RM-synthesis validation tests.
Does NOT require any real FITS reference file.

Two point-sources at known (pixel, RM) positions are injected:
  src_A: pixel (12, 10), RM = -5.0  rad/m^2
  src_B: pixel (22, 20), RM = +22.0 rad/m^2

The frequencies mimic a typical GMRT 550-750 MHz 200-channel setup.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from astropy.io import fits

OUTDIR = Path(__file__).resolve().parent / "data"

SOURCES = [
    {"name": "src_A", "x": 12, "y": 10, "rm": -5.0,  "p0": 1.0, "chi0_deg":  0.0},
    {"name": "src_B", "x": 22, "y": 20, "rm": 22.0, "p0": 0.8, "chi0_deg": 30.0},
]

# Bad channel test configuration
# BAD_CHANNEL_IDX: channel index (0-based) that will be NaN at certain pixels
# BAD_PIXELS: list of (x, y) coordinates where that channel is NaN
# FULLY_MASKED_PIXELS: list of (x, y) coordinates where ALL channels are NaN
BAD_CHANNEL_IDX = 50
BAD_PIXELS = [(12, 10), (12, 11), (13, 10)]  # src_A region with one channel bad
FULLY_MASKED_PIXELS = [(25, 25)]  # Completely masked pixel for NaN output test

NX, NY = 32, 32
N_CHAN  = 200
F_START = 550.0e6   # Hz
F_STEP  =   1.0e6   # Hz


def make_header(nx: int, ny: int, n_chan: int,
                f_start: float, f_step: float) -> fits.Header:
    h = fits.Header()
    h["SIMPLE"]  = True
    h["BITPIX"]  = -32
    h["NAXIS"]   = 4
    h["NAXIS1"]  = nx
    h["NAXIS2"]  = ny
    h["NAXIS3"]  = n_chan
    h["NAXIS4"]  = 1
    h["CTYPE1"]  = "RA---SIN"
    h["CTYPE2"]  = "DEC--SIN"
    h["CTYPE3"]  = "FREQ"
    h["CTYPE4"]  = "STOKES"
    h["CRVAL1"]  = 180.0
    h["CRVAL2"]  =   0.0
    h["CRVAL3"]  = f_start
    h["CRVAL4"]  =   1.0
    h["CRPIX1"]  = float(nx // 2 + 1)
    h["CRPIX2"]  = float(ny // 2 + 1)
    h["CRPIX3"]  = 1.0
    h["CRPIX4"]  = 1.0
    h["CDELT1"]  = -0.001         # deg
    h["CDELT2"]  =  0.001         # deg
    h["CDELT3"]  = f_step
    h["CDELT4"]  =  1.0
    h["CUNIT1"]  = "deg"
    h["CUNIT2"]  = "deg"
    h["CUNIT3"]  = "Hz"
    h["BUNIT"]   = "Jy/beam"
    h["OBJECT"]  = "TEST_RM_SYNTH"
    h["EQUINOX"] = 2000.0
    h["HISTORY"] = "Synthetic test cube for rm_synthesis validation"
    return h


def make_cubes_with_bad_channels(q_cube: np.ndarray, u_cube: np.ndarray,
                                 bad_chan_idx: int,
                                 bad_pixels: list[tuple[int, int]],
                                 fully_masked: list[tuple[int, int]]) -> tuple:
    """
    Create variants of Q/U cubes with injected bad channels (as marker values).
    
    Note: rm_synthesis.f checks for nullval=-999.0 to mask bad channels.
    We inject this marker value to simulate FITS NULL handling by CFITSIO.
    
    Args:
        q_cube, u_cube: original data arrays
        bad_chan_idx: channel index (0-based) to inject marker value
        bad_pixels: list of (x, y) where bad_chan_idx is marked bad
        fully_masked: list of (x, y) where ALL channels are marked bad
    
    Returns:
        (q_bad, u_bad): copies with IEEE NaN injected
    """
    q_bad = q_cube.copy()
    u_bad = u_cube.copy()
    
    # Use IEEE NaN - it is layout-independent and will be detected correctly
    # by the mask building code via NaN self-inequality check (x /= x)
    BAD_MARKER = np.nan
    
    # Inject NaN at one channel for specific pixels
    for x, y in bad_pixels:
        q_bad[0, bad_chan_idx, y, x] = BAD_MARKER
        u_bad[0, bad_chan_idx, y, x] = BAD_MARKER
    
    # Inject NaN at all channels for fully-masked pixels
    for x, y in fully_masked:
        q_bad[0, :, y, x] = BAD_MARKER
        u_bad[0, :, y, x] = BAD_MARKER
    
    return q_bad, u_bad


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    freq_hz   = F_START + np.arange(N_CHAN, dtype=np.float64) * F_STEP
    c         = 299792458.0
    lambda_sq = (c / freq_hz) ** 2

    q_cube = np.zeros((1, N_CHAN, NY, NX), dtype=np.float32)
    u_cube = np.zeros_like(q_cube)

    rng = np.random.default_rng(20260624)

    for src in SOURCES:
        chi0 = np.deg2rad(src["chi0_deg"])
        for ci, lsq in enumerate(lambda_sq):
            angle = 2.0 * (chi0 + src["rm"] * lsq)
            q_cube[0, ci, src["y"], src["x"]] += src["p0"] * np.cos(angle)
            u_cube[0, ci, src["y"], src["x"]] += src["p0"] * np.sin(angle)

    # Low-level Gaussian noise
    noise = 0.01
    q_cube += rng.normal(0.0, noise, q_cube.shape).astype(np.float32)
    u_cube += rng.normal(0.0, noise, u_cube.shape).astype(np.float32)

    hdr = make_header(NX, NY, N_CHAN, F_START, F_STEP)
    
    # Write normal cubes (without bad channels)
    q_path = OUTDIR / "TEST.Q.FITSCUBE"
    u_path = OUTDIR / "TEST.U.FITSCUBE"
    fits.PrimaryHDU(data=q_cube, header=hdr).writeto(q_path, overwrite=True)
    fits.PrimaryHDU(data=u_cube, header=hdr).writeto(u_path, overwrite=True)

    # Write cubes with bad channels (NaN at specific pixels)
    q_bad, u_bad = make_cubes_with_bad_channels(q_cube, u_cube,
                                                BAD_CHANNEL_IDX,
                                                BAD_PIXELS,
                                                FULLY_MASKED_PIXELS)
    q_path_bad = OUTDIR / "TEST_BADCHAN.Q.FITSCUBE"
    u_path_bad = OUTDIR / "TEST_BADCHAN.U.FITSCUBE"
    fits.PrimaryHDU(data=q_bad, header=hdr).writeto(q_path_bad, overwrite=True)
    fits.PrimaryHDU(data=u_bad, header=hdr).writeto(u_path_bad, overwrite=True)

    manifest = {
        "nx": NX, "ny": NY, "n_chan": N_CHAN,
        "freq_start_hz": F_START, "freq_step_hz": F_STEP,
        "sources": SOURCES,
        "noise_sigma": noise,
        "files": {"Q": str(q_path), "U": str(u_path)},
        "bad_channel_test": {
            "files_badchan": {"Q": str(q_path_bad), "U": str(u_path_bad)},
            "bad_channel_idx": BAD_CHANNEL_IDX,
            "bad_pixels": BAD_PIXELS,
            "fully_masked_pixels": FULLY_MASKED_PIXELS,
        },
    }
    (OUTDIR / "truth.json").write_text(json.dumps(manifest, indent=2))
    print(f"Wrote {q_path}")
    print(f"Wrote {u_path}")
    print(f"Wrote {q_path_bad}")
    print(f"Wrote {u_path_bad}")
    print(f"Wrote {OUTDIR / 'truth.json'}")


if __name__ == "__main__":
    main()
