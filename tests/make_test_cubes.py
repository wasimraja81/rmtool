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

# Multi-band tomography (T1 ticket, planning/MULTI_BAND_TOMOGRAPHY_PLAN.md):
# a second synthetic band, same RA/Dec geometry (NX/NY/CRVAL/CRPIX/CDELT)
# and the same injected SOURCES as the primary band above, but a distinct
# frequency range -- for exercising the new N-band geometry-validation
# accept path (T1's own scope stops at "geometry validated, frequency
# merge not yet implemented"; it does not need scientifically-motivated
# band parameters the way the Sec 10 thesis-derived scenario does).
BAND2_N_CHAN  = 150
BAND2_F_START = 800.0e6   # Hz
BAND2_F_STEP  =   1.0e6   # Hz


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


def make_source_cubes(freq_hz: np.ndarray, rng: np.random.Generator,
                       noise: float = 0.01) -> tuple:
    """
    Build Q/U cubes for the given frequency axis, injecting SOURCES at
    their fixed pixel/RM positions -- shared by the primary band and the
    multi-band-tomography second band (same sources, different frequency
    axis), so both bands' fixtures come from one code path.
    """
    n_chan = freq_hz.size
    c = 299792458.0
    lambda_sq = (c / freq_hz) ** 2

    q_cube = np.zeros((1, n_chan, NY, NX), dtype=np.float32)
    u_cube = np.zeros_like(q_cube)

    for src in SOURCES:
        chi0 = np.deg2rad(src["chi0_deg"])
        for ci, lsq in enumerate(lambda_sq):
            angle = 2.0 * (chi0 + src["rm"] * lsq)
            q_cube[0, ci, src["y"], src["x"]] += src["p0"] * np.cos(angle)
            u_cube[0, ci, src["y"], src["x"]] += src["p0"] * np.sin(angle)

    q_cube += rng.normal(0.0, noise, q_cube.shape).astype(np.float32)
    u_cube += rng.normal(0.0, noise, u_cube.shape).astype(np.float32)
    return q_cube, u_cube


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

    # Primary band. Refactored to call the shared make_source_cubes()
    # helper below, but the sequence of numpy calls (zeros, deterministic
    # source injection, then one rng.normal draw per cube from a
    # freshly-seeded generator) is unchanged from before that refactor --
    # this cube's bytes must stay bit-identical, since every existing
    # test's baseline output was generated from it.
    freq_hz = F_START + np.arange(N_CHAN, dtype=np.float64) * F_STEP
    rng = np.random.default_rng(20260624)
    noise = 0.01
    q_cube, u_cube = make_source_cubes(freq_hz, rng, noise)

    hdr = make_header(NX, NY, N_CHAN, F_START, F_STEP)

    # Write normal cubes (without bad channels)
    q_path = OUTDIR / "TEST.Q.FITSCUBE"
    u_path = OUTDIR / "TEST.U.FITSCUBE"
    fits.PrimaryHDU(data=q_cube, header=hdr).writeto(q_path, overwrite=True)
    fits.PrimaryHDU(data=u_cube, header=hdr).writeto(u_path, overwrite=True)

    # Multi-band tomography (T5 ticket, planning/MULTI_BAND_TOMOGRAPHY_PLAN.md):
    # split the primary band's own already-built arrays (not regenerated)
    # into two CONTIGUOUS halves, channels 1..100 and 101..200 -- for the
    # split-band identity test (contiguous multi-band split must reproduce
    # the undivided cube's output exactly). Slicing the same array
    # guarantees pixel-for-pixel correspondence to TEST.Q/U.FITSCUBE by
    # construction; a re-synthesized pair of half-cubes would not give
    # that guarantee even with the "same" RNG/signal parameters.
    n_split = N_CHAN // 2
    q_split_lo, u_split_lo = q_cube[:, :n_split], u_cube[:, :n_split]
    q_split_hi, u_split_hi = q_cube[:, n_split:], u_cube[:, n_split:]
    hdr_split_lo = make_header(NX, NY, n_split, F_START, F_STEP)
    hdr_split_hi = make_header(NX, NY, N_CHAN - n_split,
                                F_START + n_split * F_STEP, F_STEP)
    q_path_split_lo = OUTDIR / "TEST_SPLIT_LO.Q.FITSCUBE"
    u_path_split_lo = OUTDIR / "TEST_SPLIT_LO.U.FITSCUBE"
    q_path_split_hi = OUTDIR / "TEST_SPLIT_HI.Q.FITSCUBE"
    u_path_split_hi = OUTDIR / "TEST_SPLIT_HI.U.FITSCUBE"
    fits.PrimaryHDU(data=q_split_lo, header=hdr_split_lo).writeto(
        q_path_split_lo, overwrite=True)
    fits.PrimaryHDU(data=u_split_lo, header=hdr_split_lo).writeto(
        u_path_split_lo, overwrite=True)
    fits.PrimaryHDU(data=q_split_hi, header=hdr_split_hi).writeto(
        q_path_split_hi, overwrite=True)
    fits.PrimaryHDU(data=u_split_hi, header=hdr_split_hi).writeto(
        u_path_split_hi, overwrite=True)

    # Write cubes with bad channels (NaN at specific pixels)
    q_bad, u_bad = make_cubes_with_bad_channels(q_cube, u_cube,
                                                BAD_CHANNEL_IDX,
                                                BAD_PIXELS,
                                                FULLY_MASKED_PIXELS)
    q_path_bad = OUTDIR / "TEST_BADCHAN.Q.FITSCUBE"
    u_path_bad = OUTDIR / "TEST_BADCHAN.U.FITSCUBE"
    fits.PrimaryHDU(data=q_bad, header=hdr).writeto(q_path_bad, overwrite=True)
    fits.PrimaryHDU(data=u_bad, header=hdr).writeto(u_path_bad, overwrite=True)

    # Multi-band tomography (T1 ticket): a second band, same RA/Dec
    # geometry, distinct frequency range, same injected sources -- for
    # exercising the new N-band geometry-validation accept path. Drawn
    # from the same rng, continuing its stream *after* the primary band's
    # draws above, so the primary band's bytes are unaffected by this
    # addition (rng state already consumed for band 1 by this point).
    freq_hz_band2 = BAND2_F_START + np.arange(BAND2_N_CHAN, dtype=np.float64) * BAND2_F_STEP
    q_cube2, u_cube2 = make_source_cubes(freq_hz_band2, rng, noise)
    hdr2 = make_header(NX, NY, BAND2_N_CHAN, BAND2_F_START, BAND2_F_STEP)
    q_path2 = OUTDIR / "TEST_BAND2.Q.FITSCUBE"
    u_path2 = OUTDIR / "TEST_BAND2.U.FITSCUBE"
    fits.PrimaryHDU(data=q_cube2, header=hdr2).writeto(q_path2, overwrite=True)
    fits.PrimaryHDU(data=u_cube2, header=hdr2).writeto(u_path2, overwrite=True)

    # Deliberately geometry-mismatched variant of band 2 (RA reference
    # value shifted) -- for exercising the new N-band geometry-validation
    # loud-refuse path. Same data, header-only difference; only the Q file
    # is needed since the mismatch check runs per-file.
    hdr2_mismatch = make_header(NX, NY, BAND2_N_CHAN, BAND2_F_START, BAND2_F_STEP)
    hdr2_mismatch["CRVAL1"] = hdr2_mismatch["CRVAL1"] + 5.0
    q_path2_mismatch = OUTDIR / "TEST_BAND2_MISMATCH.Q.FITSCUBE"
    fits.PrimaryHDU(data=q_cube2, header=hdr2_mismatch).writeto(
        q_path2_mismatch, overwrite=True)

    manifest = {
        "nx": NX, "ny": NY, "n_chan": N_CHAN,
        "freq_start_hz": F_START, "freq_step_hz": F_STEP,
        "sources": SOURCES,
        "noise_sigma": noise,
        "files": {"Q": str(q_path), "U": str(u_path)},
        "multiband_test": {
            "band2_n_chan": BAND2_N_CHAN,
            "band2_freq_start_hz": BAND2_F_START,
            "band2_freq_step_hz": BAND2_F_STEP,
            "files_band2": {"Q": str(q_path2), "U": str(u_path2)},
            "files_band2_mismatch": {"Q": str(q_path2_mismatch)},
        },
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
    print(f"Wrote {q_path_split_lo}")
    print(f"Wrote {u_path_split_lo}")
    print(f"Wrote {q_path_split_hi}")
    print(f"Wrote {u_path_split_hi}")
    print(f"Wrote {q_path_bad}")
    print(f"Wrote {u_path_bad}")
    print(f"Wrote {q_path2}")
    print(f"Wrote {u_path2}")
    print(f"Wrote {q_path2_mismatch}")
    print(f"Wrote {OUTDIR / 'truth.json'}")


if __name__ == "__main__":
    main()
