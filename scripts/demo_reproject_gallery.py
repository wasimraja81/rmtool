#!/usr/bin/env python3
"""Generates the reprojection gallery plots in docs/user/EXAMPLES.md.

Two real, public datasets (committed as compact float32 cutouts in
docs/user/demo_data/, fetched via CDS's hips2fits cutout service --
see that directory's own small header comment for the exact fetch
commands, kept there so the images can be regenerated at a different
size/FOV without re-deriving them):

  - Haslam 408 MHz all-sky map (Remazeilles et al. 2014 reprocessing,
    via the CDS/P/Haslam408 HiPS survey), in galactic coordinates,
    AIT, whole sky. Used to show target_projection=hpx: the bundled
    cutout is deliberately in galactic/AIT, NOT already HPX, so the
    HPX output here is this script's own reproject_cubes run, not a
    pre-made HPX image.
  - MeerKAT 1.28 GHz Galactic Centre mosaic (Heywood et al. 2022,
    initial data release, via the CDS/P/MeerKAT/Galactic-Centre-
    1284MHz-StokesI HiPS survey), a 1.5 deg cutout around Sgr A* in
    galactic coordinates. Used to show target_system=icrs: the
    galactic plane's own tilt relative to RA/Dec makes the frame
    conversion visually obvious.

Neither reprojection needs a declared beam: both cases either keep or
increase resolution (never genuinely coarser than the input), so
UR-21's beam-sampling gate has nothing to refuse.
"""
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import wcs_to_celestial_frame

REPO = Path(__file__).resolve().parent.parent
BIN = REPO / "bin" / "reproject_cubes"
DATA = REPO / "docs" / "user" / "demo_data"
OUT_IMAGES = REPO / "docs" / "user" / "images"
WORK = REPO / "scratch" / "demo_gallery"

HASLAM = DATA / "haslam_408MHz_galactic.fits"
MEERKAT = DATA / "meerkat_gc_1284MHz_galactic.fits"


def run(args):
    r = subprocess.run([str(BIN)] + args, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"reproject_cubes failed: {' '.join(args)}\n{r.stdout}\n{r.stderr}")
    print(r.stdout)


def output_of(infile):
    return infile.with_name(infile.stem + "_REPROJ.FITS")


def label_graticule(ax, w, naxis1, naxis2, lon_vals, lat_vals, min_sep_frac=0.06,
                     meridians_from_centre=False):
    """Text-label a GLON/GLAT graticule directly on the lines themselves.

    WCSAxes finds border tick values by walking the plot's rectangular
    frame for coordinate crossings; that fails (too few, or clustered,
    ticks) whenever the valid-sky region doesn't meet the frame as a
    clean continuous sweep -- an HPX facet seam, or a whole-sky
    ellipse's frame only touching at isolated points. Sampling each
    line directly via this WCS's own world2pix, and placing one label
    per line wherever it lands that is far enough from every label
    already placed, works regardless of the frame's own shape.

    meridians_from_centre: for a whole-sky projection where every
    meridian (lon line, swept over latitude) converges to a single
    point at each pole, searching outward from that pole crams every
    meridian's label together there -- searching from the sweep's own
    centre instead finds each one where it is most spread out from its
    neighbours (the equator, for a whole-sky ellipse). Leave False for
    a projection (e.g. HPX) with no such single-point pole convergence
    within the plotted facets, where this does the opposite -- clumps
    labels toward the grid's own centre instead of spreading them
    across facets.
    """
    placed = []  # pixel positions already used, to space labels apart
    min_sep = min_sep_frac * min(naxis1, naxis2)
    frame = wcs_to_celestial_frame(w)

    def label_line(world_fixed_axis, fixed_val, sweep, fmt):
        sweep = np.asarray(sweep)
        # SkyCoord + world_to_pixel takes (lon, lat) positionally by the
        # frame's own convention and is therefore axis-order-safe --
        # w.wcs_world2pix(a, b, ...) instead assumes a is whichever axis
        # is CTYPE1, which for this grid's own HPX output is GLAT, not
        # GLON (w.wcs.lat == 0 here), the opposite of every other
        # projection this script plots. Passing (a, b) in the wrong
        # order there doesn't error -- it silently solves for a
        # different point, which is what produced wrong, inconsistent
        # positions before this was caught.
        if world_fixed_axis == "lon":
            sc = SkyCoord(np.full_like(sweep, fixed_val) * u.deg, sweep * u.deg, frame=frame)
        else:
            sc = SkyCoord(sweep * u.deg, np.full_like(sweep, fixed_val) * u.deg, frame=frame)
        px, py = w.world_to_pixel(sc)
        px, py = np.atleast_1d(px), np.atleast_1d(py)
        valid = np.isfinite(px) & np.isfinite(py) & (px >= 0) & (px < naxis1) & (py >= 0) & (py < naxis2)
        idx = np.flatnonzero(valid)
        if idx.size == 0:
            return
        if world_fixed_axis == "lon" and meridians_from_centre:
            idx = idx[np.argsort(np.abs(sweep[idx] - np.median(sweep)))]
        for i in idx:
            p = (px[i], py[i])
            if all((p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2 >= min_sep ** 2 for q in placed):
                placed.append(p)
                ax.text(p[0], p[1], fmt, color="white", fontsize=7,
                        ha="center", va="center",
                        bbox=dict(boxstyle="round,pad=0.1", fc="black", alpha=0.6, ec="none"))
                return

    for lon in lon_vals:
        label_line("lon", lon, np.linspace(-89, 89, 400), f"{lon:g}°")
    for lat in lat_vals:
        label_line("lat", lat, np.linspace(0, 359.9, 800), f"{lat:g}°")


def label_points(ax, w, naxis1, naxis2, points):
    """Label an exact (lon, lat) sky point each, not a searched line.

    Unlike label_graticule's line search, this places every label at
    the one point requested -- for putting a whole value sequence
    (e.g. every GLON along GLAT=0) on one single, consistent, visible
    line, the way the input plot's own equator sequence reads, rather
    than wherever a generic search happens to land.
    """
    frame = wcs_to_celestial_frame(w)
    for lon, lat, fmt in points:
        # See label_graticule's own comment -- SkyCoord's positional
        # (lon, lat) is resolved via the frame itself, not by assuming
        # which WCS axis is which, so this stays correct regardless of
        # this file's own CTYPE1/CTYPE2 axis order.
        sc = SkyCoord(lon * u.deg, lat * u.deg, frame=frame)
        px, py = w.world_to_pixel(sc)
        px, py = float(px), float(py)
        if not (np.isfinite(px) and np.isfinite(py)):
            continue
        if not (0 <= px < naxis1 and 0 <= py < naxis2):
            continue
        ax.text(px, py, fmt, color="white", fontsize=7, ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.1", fc="black", alpha=0.6, ec="none"))


def savefig_sqrt(fits_path, title, path, cmap="inferno", pct=99.5,
                  lon_labels=None, lat_labels=None):
    h = fits.getheader(fits_path)
    data = fits.getdata(fits_path)
    w = WCS(h)
    vmax = np.nanpercentile(data, pct)
    vmax = vmax if np.isfinite(vmax) and vmax > 0 else 1.0

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection=w)
    ax.imshow(np.sqrt(np.clip(data, 0, None)), origin="lower", cmap=cmap,
              vmax=np.sqrt(vmax))
    ax.set_title(title)
    ax.coords[0].set_axislabel(w.wcs.ctype[0])
    ax.coords[1].set_axislabel(w.wcs.ctype[1])
    if lon_labels is not None:
        # A whole-sky ellipse's frame only touches its own bounding box
        # at isolated points (e.g. one single GLAT value on the left/
        # right edges), so WCSAxes' border-tick search badly under-
        # labels it -- same underlying issue as HPX's facet seams, fixed
        # the same way. Left as automatic (None) for a small pointed
        # field, where the default border-tick locator already behaves
        # well (see the MeerKAT plots).
        #
        # The grid lines drawn below must use this SAME spacing as the
        # label values passed in -- grid() left on its own default
        # locator draws lines at whatever spacing it picks, which does
        # not match label_graticule's explicit values, so the labels
        # would sit on, or near, the wrong lines.
        ax.coords[0].set_ticks(spacing=(lon_labels[1] - lon_labels[0]) * u.deg)
        ax.coords[1].set_ticks(spacing=(lat_labels[1] - lat_labels[0]) * u.deg)
        ax.coords[0].set_ticklabel_visible(False)
        ax.coords[1].set_ticklabel_visible(False)
        naxis1, naxis2 = h["NAXIS1"], h["NAXIS2"]
        label_graticule(ax, w, naxis1, naxis2, lon_labels, lat_labels,
                         meridians_from_centre=True)
        # set_axislabel() above never actually draws once its axis has
        # no (border) ticks -- astropy ties the two together regardless
        # of why the ticks are off. Figure-level text is not bound by
        # that and still names the axes correctly here, since (unlike
        # the HPX case) this plot's x/y pixel axes do align with
        # GLON/GLAT -- no rotation to misrepresent.
        plt.tight_layout()
        fig.text(0.5, 0.02, w.wcs.ctype[0], ha="center")
        fig.text(0.02, 0.5, w.wcs.ctype[1], va="center", rotation="vertical")
    ax.coords.grid(color="cyan", alpha=0.3, ls="solid")
    plt.tight_layout()
    plt.savefig(path, dpi=110)
    plt.close()


def hpx_plot_with_graticule(out_fits, title, path):
    h = fits.getheader(out_fits)
    data = fits.getdata(out_fits)
    w = WCS(h)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection=w)
    vmax = np.nanpercentile(data, 99.5)
    ax.imshow(np.log10(np.clip(data, 1e-3, None)), origin="lower", cmap="inferno",
              vmax=np.log10(vmax) if vmax > 0 else None)

    # astropy's WCSAxes draws the GLON/GLAT grid through this grid's own
    # WCS -- this is what makes HPX's 45 deg pixel-axis rotation visible
    # as a diamond lattice against the plot's own pixel rows/columns,
    # not the image alone.
    ax.coords.grid(color="cyan", alpha=0.6, ls="solid")
    ax.coords[0].set_ticks(spacing=30 * u.deg)
    ax.coords[1].set_ticks(spacing=30 * u.deg)
    ax.coords[0].set_ticklabel_visible(False)
    ax.coords[1].set_ticklabel_visible(False)
    ax.set_title(title)

    naxis1, naxis2 = h["NAXIS1"], h["NAXIS2"]
    # GLON along GLAT=0 (the bright galactic plane -- one line, one
    # consistent sequence, same presentation as the input plot's own
    # equator) and GLAT along a single separate reference meridian
    # (GLON=90), not scattered across whichever line a generic search
    # happens to land on -- the two sequences would otherwise interleave
    # on top of each other near the plot's own centre.
    label_points(ax, w, naxis1, naxis2,
                 [(lon, 0, f"{lon:g}°") for lon in np.arange(0, 360, 30)])
    label_points(ax, w, naxis1, naxis2,
                 [(90, lat, f"{lat:g}°") for lat in np.arange(-60, 61, 30) if lat != 0])

    # The pixel x/y axes here are NOT GLON/GLAT -- that is the 45 deg
    # rotation this plot exists to show, so labelling them as such would
    # misstate the one thing the plot demonstrates. GLON/GLAT values are
    # read from the gridline labels above instead; "pixel" here (figure-
    # level text -- WCSAxes hides its own plain x/y axis outright, so
    # ax.set_xlabel() alone never renders) names what the frame border
    # actually is.
    plt.tight_layout()
    fig.text(0.5, 0.01, "pixel", ha="center")
    fig.text(0.01, 0.5, "pixel", va="center", rotation="vertical")
    plt.savefig(path, dpi=110)
    plt.close()


def main():
    if not BIN.exists():
        print(f"ERROR: {BIN} not built -- run: make GPU=0 OMP=1 reproject_cubes")
        return 1
    WORK.mkdir(parents=True, exist_ok=True)
    OUT_IMAGES.mkdir(parents=True, exist_ok=True)

    # --- Haslam: galactic -> HPX ---
    haslam_work = WORK / HASLAM.name
    haslam_work.write_bytes(HASLAM.read_bytes())
    output_of(haslam_work).unlink(missing_ok=True)
    run(["target_projection=hpx", f"infiles={haslam_work}", "interp=linear"])
    savefig_sqrt(haslam_work,
                 "Input: Haslam 408 MHz, galactic coordinates (AIT)",
                 OUT_IMAGES / "demo_haslam_input_galactic.png",
                 lon_labels=np.arange(0, 360, 60), lat_labels=np.arange(-60, 61, 30))
    hpx_plot_with_graticule(output_of(haslam_work),
                             "reproject_cubes target_projection=hpx -- graticule every "
                             "30 deg GLON / 30 deg GLAT",
                             OUT_IMAGES / "demo_haslam_output_hpx.png")

    # --- MeerKAT: galactic -> ICRS ---
    meerkat_work = WORK / MEERKAT.name
    meerkat_work.write_bytes(MEERKAT.read_bytes())
    output_of(meerkat_work).unlink(missing_ok=True)
    run(["target_system=icrs", "target_projection=sin", f"infiles={meerkat_work}",
         "interp=linear"])
    savefig_sqrt(meerkat_work,
                 "Input: MeerKAT 1.28 GHz Galactic Centre, galactic coordinates (SIN)",
                 OUT_IMAGES / "demo_meerkat_input_galactic.png")
    savefig_sqrt(output_of(meerkat_work),
                 "reproject_cubes target_system=icrs target_projection=sin",
                 OUT_IMAGES / "demo_meerkat_output_icrs.png")

    print(f"\nPlots written to {OUT_IMAGES}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
