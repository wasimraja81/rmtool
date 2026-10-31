#!/usr/bin/env python3
"""Where Starlink AST and astropy put the same sky point in another sky
system (galactic, FK5 or FK4), and how far apart.

reproject_cubes and match_cubes do every sky conversion with AST; astropy
(and its reproject package) use their own. This script converts the same
ICRS positions to the chosen system with both -- AST's frame built from
FITS header cards, as the tools build theirs -- and plots the difference
as a vector field: each arrow points from astropy's position to AST's, its length and
colour give the size of the offset. A second panel shows what is left
after removing the single small rotation that best fits the offsets: if
only round-off remains, the difference between the two libraries is that
rotation. It makes no claim about which library is more correct.

Cases: the whole sky; a 60x40 deg field centred on RA 0, Dec 0; the
40x40 deg cap around the south celestial pole.

Usage: compare_ast_astropy_frames.py [--system galactic|fk5|fk4] [--outdir DIR]
Builds scripts/ast_icrs_to_frame.f90 (needs gfortran, gcc and the AST
libraries, as for the rest of rmtool) into DIR on first use.
"""
import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, FK4, FK5
import astropy.units as u

REPO = Path(__file__).resolve().parent.parent
MAS = 3.6e6  # milliarcseconds per degree


def build_ast_tool(outdir):
    exe = outdir / "ast_icrs_to_frame"
    if exe.exists():
        return exe
    stub = outdir / "ast_grf_stub.o"
    subprocess.run(["gcc", "-c", str(REPO / "src" / "ast_grf_stub.c"), "-o", str(stub)], check=True)
    subprocess.run(["gfortran", "-std=gnu", "-fallow-argument-mismatch", "-ffree-line-length-none", "-O2",
                    "-o", str(exe), str(REPO / "scripts" / "ast_icrs_to_frame.f90"), str(stub),
                    "-lstarlink_ast", "-lstarlink_ast_err", "-lstarlink_ast_grf3d"], check=True)
    return exe


def ast_convert(exe, ra, dec, outdir, system):
    pts = outdir / "points_icrs.txt"
    np.savetxt(pts, np.column_stack([ra, dec]), fmt="%.17g")
    r = subprocess.run([str(exe), str(pts), system], capture_output=True, text=True, check=True)
    lb = np.loadtxt(r.stdout.splitlines())
    return lb[:, 0], lb[:, 1]


def unit(lon, lat):
    lo, la = np.radians(lon), np.radians(lat)
    return np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)


def offsets(l_ref, b_ref, l_other, b_other):
    """Offset of (l_other, b_other) from (l_ref, b_ref), in mas, along
    increasing longitude (dl cos b) and latitude (db)."""
    dl = ((l_other - l_ref + 180.0) % 360.0) - 180.0
    return dl * np.cos(np.radians(b_ref)) * MAS, (b_other - b_ref) * MAS


def fit_rotation(p, q):
    """Small rotation w (radians, as a vector) best fitting q - p = w x p."""
    d = (q - p).reshape(-1)
    # w x p as a linear function of w: rows for each point's x, y, z
    A = np.zeros((p.shape[0] * 3, 3))
    A[0::3] = np.column_stack([np.zeros(len(p)), p[:, 2], -p[:, 1]])
    A[1::3] = np.column_stack([-p[:, 2], np.zeros(len(p)), p[:, 0]])
    A[2::3] = np.column_stack([p[:, 1], -p[:, 0], np.zeros(len(p))])
    w, *_ = np.linalg.lstsq(A, d, rcond=None)
    return w


def rotate_small(p, w):
    q = p + np.cross(w, p)
    return q / np.linalg.norm(q, axis=1)[:, None]


def to_lonlat(v):
    return np.degrees(np.arctan2(v[:, 1], v[:, 0])) % 360.0, np.degrees(np.arcsin(np.clip(v[:, 2], -1, 1)))


FRAMES = {"galactic": "galactic", "fk5": FK5(equinox="J2000"), "fk4": FK4(equinox="B1950")}
LABELS = {"galactic": ("l", "b"), "fk5": ("RA", "Dec"), "fk4": ("RA", "Dec")}


def case(name, ra, dec, exe, outdir, system, w_all=None):
    l_ast, b_ast = ast_convert(exe, ra, dec, outdir, system)
    g = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs").transform_to(FRAMES[system])
    l_py, b_py = g.spherical.lon.deg, g.spherical.lat.deg
    ex, ey = offsets(l_py, b_py, l_ast, b_ast)
    mag = np.hypot(ex, ey)
    p, q = unit(l_py, b_py), unit(l_ast, b_ast)
    w = fit_rotation(p, q) if w_all is None else w_all
    lr, br = to_lonlat(rotate_small(p, w))
    rx, ry = offsets(lr, br, l_ast, b_ast)
    rmag = np.hypot(rx, ry)

    fig, axs = plt.subplots(1, 2, figsize=(14, 5.5))
    for ax, (vx, vy, m, title) in zip(axs, ((ex, ey, mag, "AST minus astropy"),
                                            (rx, ry, rmag, "after removing the best-fit rotation"))):
        lplot = ((l_py + 180.0) % 360.0) - 180.0
        qv = ax.quiver(lplot, b_py, vx, vy, m, cmap="viridis", angles="xy", pivot="middle")
        cb = fig.colorbar(qv, ax=ax)
        cb.set_label("offset (mas)")
        ax.set_xlabel(f"{system} {LABELS[system][0]} (deg, -180..180)")
        ax.set_ylabel(f"{system} {LABELS[system][1]} (deg)")
        ax.set_title(f"{name}, ICRS to {system}: {title}\nmax {m.max():.3g} mas, median {np.median(m):.3g} mas", fontsize=10)
        ax.invert_xaxis()
    fig.text(0.5, 0.005, f"Arrow: from astropy's {system} position to AST's; length and colour: size of the "
             "offset (arrow lengths are scaled per panel).", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    png = outdir / f"ast_vs_astropy_{system}_{name.replace(' ', '_').replace(',', '')}.png"
    fig.savefig(png, dpi=110)
    plt.close(fig)
    return dict(name=name, n=len(ra), max=mag.max(), min=mag.min(), median=float(np.median(mag)),
                rmax=rmag.max(), rmedian=float(np.median(rmag)), w=w, png=png)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--outdir", default="ast_vs_astropy_output")
    ap.add_argument("--system", default="galactic", choices=sorted(FRAMES))
    args = ap.parse_args()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    exe = build_ast_tool(outdir)

    ra_s, dec_s = np.meshgrid(np.arange(0.0, 360.0, 6.0), np.arange(-87.0, 88.0, 6.0))
    ra_a, dec_a = np.meshgrid(np.arange(-30.0, 30.01, 2.0) % 360.0, np.arange(-20.0, 20.01, 2.0))
    ra_b, dec_b = np.meshgrid(np.arange(0.0, 360.0, 15.0), np.arange(-89.0, -69.9, 1.0))
    results = []
    sky = case("whole sky", ra_s.ravel(), dec_s.ravel(), exe, outdir, args.system)
    results.append(sky)
    # the field cases use the whole-sky rotation, so their residual panels
    # test whether that one rotation describes them too
    results.append(case("60x40 deg field at RA 0, Dec 0", ra_a.ravel(), dec_a.ravel(), exe, outdir, args.system, sky["w"]))
    results.append(case("south celestial pole cap", ra_b.ravel(), dec_b.ravel(), exe, outdir, args.system, sky["w"]))

    w = sky["w"]
    ang = np.linalg.norm(w)
    axis_l, axis_b = to_lonlat((w / ang)[None, :])
    print(f"AST minus astropy, ICRS -> {args.system}")
    print(f"best-fit rotation (whole sky): {np.degrees(ang) * MAS:.4g} mas about the axis at "
          f"lon = {axis_l[0]:.3f}, lat = {axis_b[0]:.3f} deg ({args.system})")
    print(f"{'case':<34} {'points':>6} {'offset min':>11} {'median':>9} {'max':>9}   "
          f"{'after rotation median':>21} {'max':>9}   (mas)")
    for r in results:
        print(f"{r['name']:<34} {r['n']:>6} {r['min']:>11.4g} {r['median']:>9.4g} {r['max']:>9.4g}   "
              f"{r['rmedian']:>21.3g} {r['rmax']:>9.3g}")
    for r in results:
        print(f"plot: {r['png']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
