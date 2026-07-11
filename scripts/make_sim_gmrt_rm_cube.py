#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from astropy.io import fits


REPO = Path(__file__).resolve().parents[1]
REF_Q = Path(
    "/home/wasim/softwares/CURR_DEVEL/fitsio_utils/myfitsio.1.0/DATA/"
    "MY_CASA.QCL001.1.FITSCUBE"
)
OUTDIR = REPO / "scratch" / "sim_gmrt_rm"


def frequencies_from_header(header: fits.Header) -> np.ndarray:
    nchan = int(header["NAXIS3"])
    crval3 = float(header["CRVAL3"])
    crpix3 = float(header.get("CRPIX3", 1.0))
    cdelt3 = float(header["CDELT3"])
    chan = np.arange(1, nchan + 1, dtype=np.float64)
    return crval3 + (chan - crpix3) * cdelt3


def gaussian_stamp(nx: int, ny: int, x0: float, y0: float, fwhm_pix: float) -> np.ndarray:
    sigma = fwhm_pix / 2.354820045
    y, x = np.mgrid[0:ny, 0:nx]
    rr2 = (x - x0) ** 2 + (y - y0) ** 2
    return np.exp(-0.5 * rr2 / sigma**2).astype(np.float32)


def build_header(ref_header: fits.Header, nx: int, ny: int) -> fits.Header:
    header = fits.Header()
    header["SIMPLE"] = True
    header["BITPIX"] = -32
    header["NAXIS"] = 4
    header["NAXIS1"] = nx
    header["NAXIS2"] = ny
    header["NAXIS3"] = int(ref_header["NAXIS3"])
    header["NAXIS4"] = 1

    for key in [
        "CTYPE1",
        "CTYPE2",
        "CTYPE3",
        "CTYPE4",
        "CUNIT1",
        "CUNIT2",
        "CUNIT3",
        "CUNIT4",
        "CRVAL1",
        "CRVAL2",
        "CRVAL3",
        "CRVAL4",
        "CDELT1",
        "CDELT2",
        "CDELT3",
        "CDELT4",
        "EPOCH",
        "OBJECT",
        "TELESCOP",
        "OBSERVER",
        "BUNIT",
    ]:
        if key in ref_header:
            header[key] = ref_header[key]

    header["CRPIX1"] = float((nx + 1) / 2.0)
    header["CRPIX2"] = float((ny + 1) / 2.0)
    header["CRPIX3"] = float(ref_header.get("CRPIX3", 1.0))
    if "CRPIX4" in ref_header:
        header["CRPIX4"] = float(ref_header["CRPIX4"])
    else:
        header["CRPIX4"] = 1.0

    header["OBJECT"] = "SIM_GMRT_RM_TEST"
    header["BUNIT"] = "Jy/beam"
    header["HISTORY"] = "Synthetic Q/U/I cubes for RM-synthesis validation"
    header["HISTORY"] = "Two polarized point sources injected at known RMs"
    return header


def main() -> None:
    nx = 64
    ny = 64
    noise_sigma = 0.02
    seed = 20260624

    OUTDIR.mkdir(parents=True, exist_ok=True)
    ref_header = fits.getheader(REF_Q)
    freq_hz = frequencies_from_header(ref_header)
    lambda_sq = (299792458.0 / freq_hz) ** 2

    q_cube = np.zeros((1, len(freq_hz), ny, nx), dtype=np.float32)
    u_cube = np.zeros_like(q_cube)
    i_cube = np.zeros_like(q_cube)

    sources = [
        {
            "name": "src_rm_minus_5",
            "x": 18.0,
            "y": 21.0,
            "fwhm_pix": 2.4,
            "p0": 0.85,
            "i0": 4.0,
            "rm": -5.0,
            "chi0_deg": 12.0,
            "spec_index": -0.65,
        },
        {
            "name": "src_rm_plus_22",
            "x": 45.0,
            "y": 39.0,
            "fwhm_pix": 3.1,
            "p0": 0.62,
            "i0": 3.1,
            "rm": 22.0,
            "chi0_deg": -28.0,
            "spec_index": -0.9,
        },
    ]

    freq0 = float(freq_hz[len(freq_hz) // 2])
    rng = np.random.default_rng(seed)

    for source in sources:
        stamp = gaussian_stamp(nx, ny, source["x"], source["y"], source["fwhm_pix"])
        chi0 = np.deg2rad(source["chi0_deg"])
        spec_scale = (freq_hz / freq0) ** source["spec_index"]
        for chan, lsq in enumerate(lambda_sq):
            pol_amp = source["p0"] * spec_scale[chan]
            angle = 2.0 * (chi0 + source["rm"] * lsq)
            q_cube[0, chan] += (pol_amp * np.cos(angle) * stamp).astype(np.float32)
            u_cube[0, chan] += (pol_amp * np.sin(angle) * stamp).astype(np.float32)
            i_cube[0, chan] += (source["i0"] * spec_scale[chan] * stamp).astype(np.float32)

    # Add low-level diffuse-like correlated structure plus thermal noise.
    for chan in range(len(freq_hz)):
        background = rng.normal(0.0, noise_sigma * 0.35, size=(ny, nx)).astype(np.float32)
        smooth = (
            background
            + np.roll(background, 1, axis=0)
            + np.roll(background, -1, axis=0)
            + np.roll(background, 1, axis=1)
            + np.roll(background, -1, axis=1)
        ) / 5.0
        q_cube[0, chan] += smooth
        u_cube[0, chan] += 0.9 * smooth
        q_cube[0, chan] += rng.normal(0.0, noise_sigma, size=(ny, nx)).astype(np.float32)
        u_cube[0, chan] += rng.normal(0.0, noise_sigma, size=(ny, nx)).astype(np.float32)
        i_cube[0, chan] += rng.normal(0.0, noise_sigma * 0.5, size=(ny, nx)).astype(np.float32)

    header = build_header(ref_header, nx, ny)
    q_path = OUTDIR / "SIM_GMRT.Q.FITSCUBE"
    u_path = OUTDIR / "SIM_GMRT.U.FITSCUBE"
    i_path = OUTDIR / "SIM_GMRT.I.FITSCUBE"

    fits.PrimaryHDU(data=q_cube, header=header).writeto(q_path, overwrite=True)
    fits.PrimaryHDU(data=u_cube, header=header).writeto(u_path, overwrite=True)
    fits.PrimaryHDU(data=i_cube, header=header).writeto(i_path, overwrite=True)

    manifest = {
        "seed": seed,
        "noise_sigma_qu": noise_sigma,
        "shape_fits": [1, len(freq_hz), ny, nx],
        "shape_fortran_interpretation": [nx, ny, len(freq_hz), 1],
        "reference_q_header": str(REF_Q),
        "freq_hz_first": float(freq_hz[0]),
        "freq_hz_last": float(freq_hz[-1]),
        "freq_hz_step": float(freq_hz[1] - freq_hz[0]),
        "sources": sources,
        "files": {
            "Q": str(q_path),
            "U": str(u_path),
            "I": str(i_path),
        },
    }
    (OUTDIR / "SIM_GMRT.truth.json").write_text(json.dumps(manifest, indent=2))

    print("Wrote synthetic cubes:")
    print(q_path)
    print(u_path)
    print(i_path)
    print("Truth manifest:")
    print(OUTDIR / "SIM_GMRT.truth.json")


if __name__ == "__main__":
    main()