#!/usr/bin/env python3
"""Animate a FITS cube's 3rd axis as a sequence of 2D face-on frames.

Steps through a cube's 3rd axis (Faraday depth, frequency, velocity --
whatever it represents) one plane at a time, rendering each plane as a
plain 2D image (matplotlib) with a consistent color scale across all
frames, and encodes the sequence into a video -- like scrolling through
slices in a medical-imaging viewer, not a 3D perspective render. General
purpose: not specific to any one cube's axis meaning, though the
defaults here point at rmtool's own RM-CLEAN RESTORED.AMP cube.

Frames are rendered once to PNG (sized to the crop's own aspect ratio,
so there's no letterboxing to crop out afterward), then encoded via
ffmpeg into a single size-optimized, palette-encoded GIF -- small
enough to commit and embed inline in README.md, which is the only
output this tool needs to produce.

Style:
- Spatial pixel bounds (--ra-min/--ra-max/--dec-min/--dec-max) are
  1-indexed, FITS/DS9 convention, inclusive
- Axis-3 bounds (--axis3-min/--axis3-max) are physical WCS units
  (e.g. rad/m^2 for a Faraday-depth cube), converted to pixel planes
  via the cube's own WCS at runtime -- never hardcoded plane indices.
  Omit both to animate every plane in the cube.
- `--plane-stride` skips planes at render time (not just at encode
  time), so a lighter preview doesn't waste time rendering frames it's
  about to discard.

Example:
  ~/venv/rmtool/bin/python scripts/animate_fits_cube.py --test-frames 5
  ~/venv/rmtool/bin/python scripts/animate_fits_cube.py
"""

from __future__ import annotations

import argparse
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

DEFAULT_FITS = (
    "/scratch/tmp/pipeline_multiband_e2e/"
    "multiband_wallaby_emu_cleaned.RESTORED.AMP.RMCUBE.FITS"
)


@dataclass
class Crop:
    ra_min: int
    ra_max: int
    dec_min: int
    dec_max: int


def compute_axis3_plane_slice(header: fits.Header, lo: float, hi: float) -> slice:
    """Axis-3 physical bounds -> 0-indexed numpy slice, via the cube's own WCS."""
    w3 = WCS(header).sub([3])
    p_lo, p_hi = w3.wcs_world2pix([lo, hi], 1)[0]
    k_min = math.ceil(min(p_lo, p_hi))
    k_max = math.floor(max(p_lo, p_hi))
    if k_max < k_min:
        raise ValueError(f"Axis-3 range [{lo}, {hi}] contains no whole planes in this cube")
    return slice(k_min - 1, k_max)


def block_mean_2d(arr: np.ndarray, factor: int) -> np.ndarray:
    """Block-average downsample the last two axes of a (plane, dec, ra) array."""
    if factor <= 1:
        return arr
    n, ny, nx = arr.shape
    ny_t, nx_t = (ny // factor) * factor, (nx // factor) * factor
    arr = arr[:, :ny_t, :nx_t]
    arr = arr.reshape(n, ny_t // factor, factor, nx_t // factor, factor)
    return arr.mean(axis=(2, 4))


def load_cube_slice(
    fits_path: str, crop: Crop, axis3_min: float | None, axis3_max: float | None
) -> tuple[np.ndarray, fits.Header, slice]:
    """Read only the requested (axis3, dec, ra) subregion via a memmapped slice."""
    with fits.open(fits_path, memmap=True) as hdul:
        header = hdul[0].header
        if axis3_min is None and axis3_max is None:
            axis3_slice = slice(0, header["NAXIS3"])
        else:
            axis3_slice = compute_axis3_plane_slice(header, axis3_min, axis3_max)
        ra_slice = slice(crop.ra_min - 1, crop.ra_max)
        dec_slice = slice(crop.dec_min - 1, crop.dec_max)
        data = np.asarray(hdul[0].data[axis3_slice, dec_slice, ra_slice])
    return data, header, axis3_slice


def render_frames(
    cube: np.ndarray,
    axis3_vals: np.ndarray,
    args: argparse.Namespace,
    frames_dir: Path,
    plane_indices: np.ndarray,
) -> None:
    nz = cube[cube > 0]
    vmin, vmax = np.percentile(nz, [args.vmin_percentile, args.vmax_percentile])

    frames_dir.mkdir(parents=True, exist_ok=True)
    # Figure sized to the crop's own aspect ratio -- matplotlib's default
    # equal-aspect imshow would otherwise letterbox/pillarbox a mismatched
    # figure size, wasting encoded pixels (and GIF palette budget) on flat
    # background bars that then have to be cropped out after the fact.
    n_dec, n_ra = cube.shape[1], cube.shape[2]
    height = 1024
    width = max(1, round(height * n_ra / n_dec))
    dpi = 100
    fig = plt.figure(figsize=(width / dpi, height / dpi), dpi=dpi, facecolor=args.background_color)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_facecolor(args.background_color)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    n_total = cube.shape[0]
    n_frames = len(plane_indices)
    for i, plane in enumerate(plane_indices):
        ax.clear()
        ax.set_xticks([])
        ax.set_yticks([])
        ax.imshow(cube[plane], origin="lower", cmap=args.colormap, vmin=vmin, vmax=vmax)
        ax.text(
            0.02,
            0.97,
            args.label_fmt.format(value=axis3_vals[plane], index=plane + 1, total=n_total),
            transform=ax.transAxes,
            color="white",
            fontsize=14,
            va="top",
            ha="left",
            fontfamily="monospace",
        )
        fig.savefig(frames_dir / f"frame_{i:04d}.png", facecolor=fig.get_facecolor())
        if (i + 1) % 30 == 0 or i == n_frames - 1:
            print(f"[render_frames] {i + 1}/{n_frames}", file=sys.stderr)

    plt.close(fig)


def encode_gif(frames_dir: Path, args: argparse.Namespace) -> None:
    out_gif = Path(args.out_gif)
    out_gif.parent.mkdir(parents=True, exist_ok=True)

    frame_pattern = str(frames_dir / "frame_%04d.png")
    filter_base = f"fps={args.gif_fps},scale={args.gif_width}:-1:flags=lanczos"

    palette = frames_dir.parent / "palette.png"
    subprocess.run(
        [
            "ffmpeg", "-y",
            "-framerate", str(args.gif_fps),
            "-i", frame_pattern,
            "-vf", f"{filter_base},palettegen",
            str(palette),
        ],
        check=True,
    )
    subprocess.run(
        [
            "ffmpeg", "-y",
            "-framerate", str(args.gif_fps),
            "-i", frame_pattern,
            "-i", str(palette),
            "-filter_complex",
            f"{filter_base}[x];[x][1:v]paletteuse=dither=bayer:bayer_scale=3",
            str(out_gif),
        ],
        check=True,
    )
    print(f"[encode_gif] wrote {out_gif}", file=sys.stderr)

    if not args.keep_frames:
        shutil.rmtree(frames_dir, ignore_errors=True)
        palette.unlink(missing_ok=True)


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Animate a FITS cube's 3rd axis as a sequence of 2D face-on frames"
    )
    p.add_argument("--fits", default=DEFAULT_FITS, help="Path to the FITS cube")
    p.add_argument("--ra-min", type=int, default=7000, help="1-indexed FITS RA pixel, inclusive")
    p.add_argument("--ra-max", type=int, default=10300)
    p.add_argument("--dec-min", type=int, default=2000, help="1-indexed FITS Dec pixel, inclusive")
    p.add_argument("--dec-max", type=int, default=6000)
    p.add_argument(
        "--axis3-min",
        type=float,
        default=None,
        help="Physical units (e.g. rad/m^2); omit with --axis3-max to animate every plane",
    )
    p.add_argument("--axis3-max", type=float, default=None)
    p.add_argument("--spatial-downsample", type=int, default=2, help="Block-average RA/Dec")
    p.add_argument(
        "--plane-stride",
        type=int,
        default=2,
        help="Animate every Nth plane (keeps the GIF small; skipped planes are never rendered)",
    )
    p.add_argument("--vmin-percentile", type=float, default=50.0)
    p.add_argument("--vmax-percentile", type=float, default=99.5)
    p.add_argument("--colormap", default="inferno")
    p.add_argument(
        "--label-fmt",
        default="RM = {value:+.1f} rad/m^2   ({index}/{total})",
        help="Per-frame text label; {value}/{index}/{total} available",
    )
    p.add_argument("--gif-fps", type=int, default=5)
    p.add_argument("--gif-width", type=int, default=360)
    p.add_argument("--background-color", default="#0b0b14")
    p.add_argument("--frames-dir", default="scratch/fits_cube_anim/frames")
    p.add_argument("--out-gif", default="docs/user/images/rmclean_channelmap_demo.gif")
    p.add_argument("--keep-frames", action="store_true")
    p.add_argument(
        "--test-frames",
        type=int,
        default=None,
        help="Render only the first N (strided) planes and skip encoding (smoke test)",
    )
    return p


def main() -> int:
    args = build_argparser().parse_args()
    crop = Crop(args.ra_min, args.ra_max, args.dec_min, args.dec_max)

    print(f"[main] loading from {args.fits}", file=sys.stderr)
    cube, header, axis3_slice = load_cube_slice(args.fits, crop, args.axis3_min, args.axis3_max)
    print(f"[main] loaded shape (axis3,dec,ra)={cube.shape}", file=sys.stderr)

    cube = block_mean_2d(cube, args.spatial_downsample)
    print(f"[main] after spatial downsample: {cube.shape}", file=sys.stderr)

    w3 = WCS(header).sub([3])
    plane_idx = np.arange(axis3_slice.start + 1, axis3_slice.start + cube.shape[0] + 1, dtype=float)
    axis3_vals = w3.wcs_pix2world(plane_idx, 1)[0]

    plane_indices = np.arange(0, cube.shape[0], args.plane_stride)
    if args.test_frames is not None:
        plane_indices = plane_indices[: args.test_frames]
    frames_dir = Path(args.frames_dir)
    render_frames(cube, axis3_vals, args, frames_dir, plane_indices)

    if args.test_frames is not None:
        print(f"[main] smoke test complete: {len(plane_indices)} frames in {frames_dir}", file=sys.stderr)
        return 0

    encode_gif(frames_dir, args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
