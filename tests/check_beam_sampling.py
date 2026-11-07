#!/usr/bin/env python3
"""Beam sampling on the output grid (UR-21,
docs/dev/BRANCH_STATUS_update-reproject.md).

Reprojection interpolates; it reproduces what the output grid can
sample and aliases what it cannot. Both reprojecting tools therefore
refuse a grid COARSER than an infile's own whose cells the beam's minor
axis spans fewer than 2.0 times, warn between 2.0 and 2.5 (barely
Nyquist), and never enlarge a beam themselves. On a grid no coarser
than the infile's own nothing new can alias: an infile that already
undersamples its beam is warned about, and one with no beam information
is accepted. On a coarser grid an infile with no beam is refused.

Checked here, each against the tool's exit status, its message and
whether an output was written:
  reproject_cubes, no reffile and with a reffile: pass, warn, refuse,
  the printed advice, beams from a header and from a per-channel BEAMS
  table (narrowest channel used), no beam on an equal and a coarser
  grid, no beam on an equal grid whose reference point lies away from
  the data (accepted; see below), an already-undersampled infile on its
  own grid;
  match_cubes: stages=reproject and stages=both order=reproject_convolve
  refuse at the native beam; order=convolve_reproject refuses a derived
  or explicit common beam under 2.0 cells, warns between 2.0 and 2.5,
  and runs the beam reproject_cubes' advice names, writing that beam to
  the output header.
"""
import re
import shlex
import subprocess
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / "bin"
WORK = REPO / "tests" / "output" / "beam_sampling"
CELL_DEG = 0.001                     # 3.6 arcsec input cells
N = 64

failures, n_checks = [], 0


def check(ok, msg):
    global n_checks
    n_checks += 1
    if not ok:
        failures.append(msg)


def cube(path, beam_arcsec=None, beams_table=None, nchan=1, cell_deg=CELL_DEG):
    """beam_arcsec: header BMAJ=BMIN; beams_table: per-channel BMIN list
    (BMAJ = BMIN), written as a CASA BEAMS table with CASAMBM=T."""
    h = fits.Header()
    h["SIMPLE"] = True; h["BITPIX"] = -32
    h["NAXIS"] = 3 if nchan > 1 else 2
    h["NAXIS1"] = N; h["NAXIS2"] = N
    h["CTYPE1"] = "RA---SIN"; h["CTYPE2"] = "DEC--SIN"
    h["CRVAL1"] = 180.0; h["CRVAL2"] = -30.0
    h["CRPIX1"] = (N + 1) / 2; h["CRPIX2"] = (N + 1) / 2
    h["CDELT1"] = -cell_deg; h["CDELT2"] = cell_deg
    if nchan > 1:
        h["NAXIS3"] = nchan; h["CTYPE3"] = "FREQ"; h["CRVAL3"] = 1.0e9
        h["CRPIX3"] = 1.0; h["CDELT3"] = 1.0e6
    if beam_arcsec is not None:
        h["BMAJ"] = beam_arcsec / 3600; h["BMIN"] = beam_arcsec / 3600; h["BPA"] = 0.0
    shape = (nchan, N, N) if nchan > 1 else (N, N)
    hdus = [fits.PrimaryHDU(np.ones(shape, np.float32), h)]
    if beams_table is not None:
        hdus[0].header["CASAMBM"] = True
        b = np.array(beams_table, dtype=np.float32)
        cols = [fits.Column("BMAJ", "E", "arcsec", array=b),
                fits.Column("BMIN", "E", "arcsec", array=b),
                fits.Column("BPA", "E", "deg", array=np.zeros_like(b)),
                fits.Column("CHAN", "J", array=np.arange(len(b), dtype=np.int32)),
                fits.Column("POL", "J", array=np.zeros(len(b), dtype=np.int32))]
        t = fits.BinTableHDU.from_columns(cols); t.name = "BEAMS"
        hdus.append(t)
    fits.HDUList(hdus).writeto(path, overwrite=True)


def run(args):
    r = subprocess.run([str(a) for a in args], capture_output=True, text=True, cwd=WORK)
    return r.returncode, r.stdout + r.stderr


def outputs(stem):
    return sorted(p.name for p in WORK.glob(stem + "_*.FITS"))


def clean(stem):
    for p in WORK.glob(stem + "_*.FITS"):
        p.unlink()


def expect(label, args, stem, exit_ok, written, must=(), must_not=()):
    clean(stem)
    code, out = run(args)
    check((code == 0) == exit_ok, f"{label}: exit {code}, expected {'0' if exit_ok else 'non-zero'}")
    check(bool(outputs(stem)) == written, f"{label}: outputs {outputs(stem)}, expected {'some' if written else 'none'}")
    for s in must:
        check(s in out, f"{label}: message lacks {s!r}")
    for s in must_not:
        check(s not in out, f"{label}: message unexpectedly has {s!r}")
    return out


def cdelt(c):
    return [f"target_cdelt_lon={c}", f"target_cdelt_lat={c}"]


def main():
    for tool in ("reproject_cubes", "match_cubes"):
        if not (BIN / tool).exists():
            print(f"ERROR: {BIN / tool} not built"); return 1
    WORK.mkdir(parents=True, exist_ok=True)
    R, M = BIN / "reproject_cubes", BIN / "match_cubes"
    NOGRID = ["target_system=icrs"]

    # --- reproject_cubes, beam in the header (10" on 3.6" cells) ---
    cube(WORK / "b10.fits", beam_arcsec=10.0)
    expect("10\" beam, same cells", [R, *NOGRID, "infiles=b10.fits"], "b10", True, True,
           must_not=("barely Nyquist", "ERROR"))
    expect("10\" beam, 4.68\" cells (2.14)", [R, *NOGRID, *cdelt(0.0013), "infiles=b10.fits"],
           "b10", True, True, must=("barely Nyquist",))
    out = expect("10\" beam, 5.4\" cells (1.85)", [R, *NOGRID, *cdelt(0.0015), "infiles=b10.fits"],
                 "b10", False, False,
                 must=("spans only 1.85 output cells; at least 2.0 are needed",
                       "match_cubes stages=both order=convolve_reproject",
                       "convolve_cubes infiles=b10.fits", "beamfiles="))
    # the advice must run: take its match_cubes line, add beamfiles (this
    # cube has a header beam, no BEAMS table, as the advice says)
    adv = next((l.strip() for l in out.splitlines() if l.strip().startswith("match_cubes stages=both")), "")
    check(bool(adv), "advice: no match_cubes command printed")
    if adv:
        cube(WORK / "b10c.fits", beam_arcsec=10.0, nchan=3)
        (WORK / "b10c_beams.txt").write_text("".join(f"{i} 10.0 10.0 0.0\n" for i in (1, 2, 3)))
        words = [w.replace("infiles=b10.fits", "infiles=b10c.fits") for w in shlex.split(adv)[1:]]
        clean("b10c")
        code, aout = run([M, *words, "beamfiles=b10c_beams.txt"])
        check(code == 0 and outputs("b10c"), f"advice command did not run: exit {code}")
        tb = float(re.search(r"target_bmaj=([0-9.]+)", adv).group(1))
        if outputs("b10c"):
            h = fits.getheader(WORK / outputs("b10c")[0])
            check(abs(h["BMAJ"] * 3600 - tb) < 1e-3,
                  f"advice command: output BMAJ {h['BMAJ']*3600:.3f}\", expected {tb}")
            check("barely Nyquist" not in aout, "advice command: still warned barely Nyquist")

    # --- beams from a BEAMS table: the narrowest channel decides ---
    cube(WORK / "tbl.fits", beams_table=[10.0, 6.0, 10.0], nchan=3)
    expect("BEAMS table, narrowest 6\", 4.68\" cells (1.28)",
           [R, *NOGRID, *cdelt(0.0013), "infiles=tbl.fits"], "tbl", False, False,
           must=("BMIN=6.000",))

    # --- no beam information ---
    cube(WORK / "nb.fits")
    expect("no beam, same cells", [R, *NOGRID, "infiles=nb.fits"], "nb", True, True,
           must=("has no beam information; not checked",))
    expect("no beam, coarser cells", [R, *NOGRID, *cdelt(0.002), "infiles=nb.fits"], "nb", False, False,
           must=("has no beam information", "aliasing cannot be ruled out"))
    # Pixel sizes are compared at each grid's reference pixel, where the
    # projection adds no distortion. With the output's reference point
    # 0.46 deg from the data, its SIN pixels there cover
    # sec(0.46 deg) - 1 = 3.3e-5 more sky than nominal; that is projection geometry, not a
    # coarser grid, so equal pixel sizes must still be accepted.
    away = ["target_center=user_defined", "target_lon=180.35", "target_lat=-30.35"]
    for tool, extra in ((R, []), (M, ["stages=reproject"])):
        expect(f"{tool.name}: no beam, same cells, reference point 0.46 deg away",
               [tool, *extra, *NOGRID, *away, "infiles=nb.fits"], "nb", True, True,
               must=("has no beam information; not checked",), must_not=("ERROR",))
    # A coarser pixel size set by the user is still refused, however small.
    expect("no beam, cells 1% coarser", [R, *NOGRID, *cdelt(0.00101), "infiles=nb.fits"], "nb", False, False,
           must=("aliasing cannot be ruled out",))

    # --- an infile already undersampled on its own grid: warned, not refused ---
    cube(WORK / "b5.fits", beam_arcsec=5.0)          # 1.39 of its own cells
    expect("5\" beam on its own 3.6\" cells", [R, *NOGRID, "infiles=b5.fits"], "b5", True, True,
           must=("already undersampled",), must_not=("at least 2.0 are needed",))

    # --- the reffile path ---
    cube(WORK / "coarse_ref.fits", cell_deg=0.002)
    cube(WORK / "b10r.fits", beam_arcsec=10.0)
    expect("reffile with 7.2\" cells", [R, "mode=reference", "reffile=coarse_ref.fits", "infiles=b10r.fits"],
           "b10r", False, False, must=("at least 2.0 are needed", "footprint_mode=reference reffile=coarse_ref.fits"))

    # --- match_cubes ---
    cube(WORK / "m.fits", beam_arcsec=10.0, nchan=3)
    (WORK / "m_beams.txt").write_text("".join(f"{i} 10.0 10.0 0.0\n" for i in (1, 2, 3)))
    G = [*NOGRID, *cdelt(0.002), "infiles=m.fits"]
    expect("match stages=reproject, coarse", [M, "stages=reproject", *G], "m", False, False,
           must=("at least 2.0 are needed",))
    expect("match order=reproject_convolve, coarse", [M, "stages=both", "order=reproject_convolve", *G,
           "beamfiles=m_beams.txt"], "m", False, False, must=("at least 2.0 are needed",))
    expect("match derived common beam too narrow", [M, "stages=both", *G, "beamfiles=m_beams.txt"],
           "m", False, False, must=("the derived common beam", "at least 2.0 are needed"))
    expect("match explicit target 12\" (1.67 cells)", [M, "stages=both", *G, "beamfiles=m_beams.txt",
           "target_bmaj=12", "target_bmin=12", "target_bpa=0"], "m", False, False,
           must=("target_bmin=12.0000", "at least 2.0 are needed"))
    expect("match explicit target 16\" (2.22 cells)", [M, "stages=both", *G, "beamfiles=m_beams.txt",
           "target_bmaj=16", "target_bmin=16", "target_bpa=0"], "m", True, True, must=("barely Nyquist",))
    expect("match explicit target 18.01\" (2.50 cells)", [M, "stages=both", *G, "beamfiles=m_beams.txt",
           "target_bmaj=18.01", "target_bmin=18.01", "target_bpa=0"], "m", True, True,
           must_not=("barely Nyquist", "ERROR"))

    print(f"beam sampling checks: {n_checks - len(failures)}/{n_checks} passed")
    for f in failures:
        print(f"  FAIL: {f}")
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.exit(main())
