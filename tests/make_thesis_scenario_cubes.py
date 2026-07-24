#!/usr/bin/env python3
"""
Synthetic P-band / L-band Q/U cube generator reproducing the multi-band
tomography scenario in planning/MULTI_BAND_TOMOGRAPHY_PLAN.md Sec 10,
grounded in Raja (2014) "Faraday Slicing Polarized Radio Sources",
Chapter 6 Table 6.1/6.2.

Bands (Table 6.1, exact):
  P: nu_c=300 MHz, dnu=30 MHz   (edges 285-315 MHz)
  L: nu_c=1200 MHz, dnu=120 MHz (edges 1140-1260 MHz)

Sky model, all at one pixel (32x32 image, source at the centre):
  1. Point source (Faraday-thin): RM=-100 rad/m^2, amplitude 15 Jy/(rad/m^2)
     -- reproduced exactly from Table 6.2.
  2. Faraday-thick top-hat: RM=100..130 rad/m^2 (centre 115, thickness 30),
     integrated amplitude 5 Jy/(rad/m^2) -- reproduced exactly from
     Table 6.2 ("modeled as a top-hat function along RM"). Approximated
     here as a fine sum of many closely-spaced thin components spanning
     the top-hat (slice spacing << the achievable RM resolution in this
     scenario, so the discretisation error is negligible relative to what
     is being measured).
  3. F2/F3 -- an addition beyond the thesis's own demo (plan Sec 10,
     confirmed with user 2026-07-20/21): a close Faraday-thin pair at
     RM=250 and RM=290 (separation 40 rad/m^2), amplitude 8 Jy/(rad/m^2)
     each, chosen using the thesis's own real numbers (Table 6.1) to be
     resolved by P alone/combined but blended at L alone.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from astropy.io import fits

OUTDIR = Path(__file__).resolve().parent / "data"

NX, NY = 32, 32
SRC_X, SRC_Y = 16, 16  # centre pixel

C_LIGHT = 299792458.0

# --- Bands (Table 6.1, exact) ---
BAND_P = {"name": "P", "nu_c": 300.0e6, "dnu": 30.0e6, "n_chan": 61}
BAND_L = {"name": "L", "nu_c": 1200.0e6, "dnu": 120.0e6, "n_chan": 121}

# --- Sky model (Table 6.2 + Sec 10's F2/F3 addition) ---
POINT_RM = -100.0
POINT_AMP = 15.0

THICK_RM_LO = 100.0
THICK_RM_HI = 130.0
THICK_AMP_TOTAL = 5.0
THICK_N_SLICES = 61  # spacing = 30/60 = 0.5 rad/m^2, well below any
                      # achievable delta-RM in this scenario (14.7-250.1)

F2_RM = 250.0
F3_RM = 290.0
F23_AMP = 8.0

CHI0_DEG = 0.0


def make_header(nx: int, ny: int, n_chan: int, nu_c: float, dnu: float) -> fits.Header:
    f_step = dnu / n_chan
    f_start = nu_c - dnu / 2.0 + f_step / 2.0  # first channel's centre freq
    h = fits.Header()
    h["SIMPLE"] = True
    h["BITPIX"] = -32
    h["NAXIS"] = 4
    h["NAXIS1"] = nx
    h["NAXIS2"] = ny
    h["NAXIS3"] = n_chan
    h["NAXIS4"] = 1
    h["CTYPE1"] = "RA---SIN"
    h["CTYPE2"] = "DEC--SIN"
    h["CTYPE3"] = "FREQ"
    h["CTYPE4"] = "STOKES"
    h["CRVAL1"] = 180.0
    h["CRVAL2"] = 0.0
    h["CRVAL3"] = f_start
    h["CRVAL4"] = 1.0
    h["CRPIX1"] = float(nx // 2 + 1)
    h["CRPIX2"] = float(ny // 2 + 1)
    h["CRPIX3"] = 1.0
    h["CRPIX4"] = 1.0
    h["CDELT1"] = -0.001
    h["CDELT2"] = 0.001
    h["CDELT3"] = f_step
    h["CDELT4"] = 1.0
    h["CUNIT1"] = "deg"
    h["CUNIT2"] = "deg"
    h["CUNIT3"] = "Hz"
    h["BUNIT"] = "Jy/beam"
    h["OBJECT"] = "MULTIBAND_THESIS_SCENARIO"
    h["EQUINOX"] = 2000.0
    h["HISTORY"] = "Sec 10 thesis-scenario cube (Raja 2014, Table 6.1/6.2)"
    return h


def components():
    """List of (rm, amplitude) point-like contributions approximating the
    full sky model (point source + thick top-hat discretised into thin
    slices + F2/F3)."""
    comps = [(POINT_RM, POINT_AMP)]
    thick_slice_amp = THICK_AMP_TOTAL / THICK_N_SLICES
    for rm in np.linspace(THICK_RM_LO, THICK_RM_HI, THICK_N_SLICES):
        comps.append((float(rm), thick_slice_amp))
    comps.append((F2_RM, F23_AMP))
    comps.append((F3_RM, F23_AMP))
    return comps


def make_band_cubes(band: dict, rng: np.random.Generator, noise: float = 0.005):
    n_chan = band["n_chan"]
    f_step = band["dnu"] / n_chan
    f_start = band["nu_c"] - band["dnu"] / 2.0 + f_step / 2.0
    freq_hz = f_start + np.arange(n_chan, dtype=np.float64) * f_step
    lambda_sq = (C_LIGHT / freq_hz) ** 2

    q_cube = np.zeros((1, n_chan, NY, NX), dtype=np.float32)
    u_cube = np.zeros_like(q_cube)

    chi0 = np.deg2rad(CHI0_DEG)
    for rm, amp in components():
        angle = 2.0 * (chi0 + rm * lambda_sq)
        q_cube[0, :, SRC_Y, SRC_X] += (amp * np.cos(angle)).astype(np.float32)
        u_cube[0, :, SRC_Y, SRC_X] += (amp * np.sin(angle)).astype(np.float32)

    q_cube += rng.normal(0.0, noise, q_cube.shape).astype(np.float32)
    u_cube += rng.normal(0.0, noise, u_cube.shape).astype(np.float32)
    return q_cube, u_cube, freq_hz


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20260721)

    written = {}
    for band in (BAND_P, BAND_L):
        q_cube, u_cube, freq_hz = make_band_cubes(band, rng)
        hdr = make_header(NX, NY, band["n_chan"], band["nu_c"], band["dnu"])
        q_path = OUTDIR / f"THESIS_{band['name']}.Q.FITSCUBE"
        u_path = OUTDIR / f"THESIS_{band['name']}.U.FITSCUBE"
        fits.PrimaryHDU(data=q_cube, header=hdr).writeto(q_path, overwrite=True)
        fits.PrimaryHDU(data=u_cube, header=hdr).writeto(u_path, overwrite=True)
        written[band["name"]] = {
            "Q": str(q_path), "U": str(u_path),
            "n_chan": band["n_chan"], "nu_c": band["nu_c"], "dnu": band["dnu"],
        }
        print(f"Wrote {q_path}")
        print(f"Wrote {u_path}")

    manifest = {
        "src_x": SRC_X, "src_y": SRC_Y,
        "point_rm": POINT_RM, "point_amp": POINT_AMP,
        "thick_rm_lo": THICK_RM_LO, "thick_rm_hi": THICK_RM_HI,
        "thick_amp_total": THICK_AMP_TOTAL,
        "f2_rm": F2_RM, "f3_rm": F3_RM, "f23_amp": F23_AMP,
        "bands": written,
    }
    truth_path = OUTDIR / "thesis_scenario_truth.json"
    truth_path.write_text(json.dumps(manifest, indent=2))
    print(f"Wrote {truth_path}")


if __name__ == "__main__":
    main()
