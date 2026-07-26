#!/usr/bin/env python3
"""
Plot tests/thesis_scenario_rmclean.f90's own CSV output in the same
panel style as Raja (2014) "Faraday Slicing Polarized Radio Sources"
Chapter 6 Figures 6.1/6.2/6.3: dirty Re/Im/Amp + phase (left/top),
cleaned (restored) amplitude + phase (right/top), the true input model
and the lambda-squared coverage (bottom two rows, same on both sides --
matching the thesis's own "bottom two rows have redundancy" layout,
Figure 6-1's caption).

Usage:
  plot_thesis_scenario_rmclean.py <csv_dir> <out_png_dir>

Reads <csv_dir>/<slug>_profile.csv and <slug>_lsq.csv for each of the
three thesis cases (fig6_1_Lalone, fig6_2_Palone, fig6_3_PLcombined) and
writes <out_png_dir>/<slug>.png.
"""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASES = [
    ("fig6_1_Lalone", "Figure 6-1 analogue: L-band alone (poor RM-resolution)"),
    ("fig6_2_Palone", "Figure 6-2 analogue: P-band alone (poor sensitivity to extended structure)"),
    ("fig6_3_PLcombined", "Figure 6-3 analogue: P-band + L-band combined"),
]

POINT_RM, POINT_AMP = -100.0, 15.0
THICK_LO, THICK_HI, THICK_AMP = 100.0, 130.0, 5.0
NITER, GAIN = 1000, 0.100


def load_profile(csv_dir: Path, slug: str):
    path = csv_dir / f"{slug}_profile.csv"
    data = np.genfromtxt(path, delimiter=",", names=True)
    return data


def load_lsq(csv_dir: Path, slug: str):
    path = csv_dir / f"{slug}_lsq.csv"
    data = np.genfromtxt(path, delimiter=",", names=True)
    return data["lambda_sq"], data["band"]


def plot_case(csv_dir: Path, out_dir: Path, slug: str, title: str):
    prof = load_profile(csv_dir, slug)
    lsq, band = load_lsq(csv_dir, slug)
    rm = prof["RM"]

    fig, axes = plt.subplots(4, 2, figsize=(11, 13))
    fig.suptitle(title, fontsize=13)

    # Row 1: dirty Re/Im/Amp (left), cleaned/restored amplitude (right)
    ax = axes[0, 0]
    ax.plot(rm, prof["dirty_re"], color="tab:red", lw=0.7, label="Re")
    ax.plot(rm, prof["dirty_im"], color="tab:green", lw=0.7, label="Im")
    ax.plot(rm, prof["dirty_amp"], color="black", lw=1.0, label="Amp")
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_title("Dirty RM-Profiles")
    ax.set_ylabel("Real, Imag & Amp")
    ax.legend(fontsize=7, loc="upper right")

    ax = axes[0, 1]
    ax.plot(rm, prof["clean_amp"], color="black", lw=1.0)
    ax.set_title("Cleaned RM profile")
    ax.set_ylabel("Restored Amplitude")
    ax.text(0.98, 0.95, f"niter: {NITER}\nloop-gain: {GAIN:.3f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8)

    # Row 2: dirty phase (left), restored phase (right)
    ax = axes[1, 0]
    ax.plot(rm, prof["dirty_phase_deg"], ".", color="black", ms=1.5)
    ax.set_ylabel("Phase (degrees)")
    ax.set_ylim(-180, 180)

    ax = axes[1, 1]
    ax.plot(rm, prof["clean_phase_deg"], ".", color="black", ms=1.5)
    ax.set_ylabel("Restored Phase (deg)")
    ax.set_ylim(-180, 180)

    # Rows 3-4: input model and lambda^2 coverage, shown on both columns
    # (matching the thesis's own redundant bottom-two-rows layout).
    for col in range(2):
        ax = axes[2, col]
        ax.plot(rm, prof["input_amp"], color="black", lw=1.2)
        ax.set_ylabel("Input Amplitude")
        ax.set_xlabel("Rotation Measure (rad/m**2)")
        ax.text(0.02, 0.95, f"AMP: {POINT_AMP:.1f}\nPHA: 0.0\nRM: {POINT_RM:.0f} rad/m**2",
                transform=ax.transAxes, ha="left", va="top", fontsize=7)
        ax.text(0.98, 0.95,
                f"Faraday thick component\nFaraday thickness: {THICK_HI-THICK_LO:.1f}\n"
                f"Amplitude: {THICK_AMP:.1f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=7, color="tab:orange")

        ax = axes[3, col]
        lsq_min, lsq_max = float(lsq.min()), float(lsq.max())
        # Draw one red segment PER BAND, not a single bar spanning
        # min-to-max: for multi-band data (e.g. P+L combined) a single
        # bar would misleadingly imply continuous coverage across the
        # large GAP between bands. Band membership comes straight from
        # the Fortran side's own band_offset/band_nz (written as the
        # "band" column), not inferred from lambda^2-spacing statistics
        # -- a spacing-based gap threshold is fragile because one band's
        # own channel-to-channel spacing (e.g. P-band's, coarser near
        # the low-frequency end) can be orders of magnitude larger than
        # another band's (e.g. L-band's), so no single global threshold
        # can separate a real inter-band gap from a coarser band's own
        # ordinary channel spacing.
        for b in np.unique(band):
            b_lsq = lsq[band == b]
            ax.hlines(1.0, b_lsq.min(), b_lsq.max(), color="red", lw=3)
        ax.set_xlim(left=0.0)
        ax.set_ylim(0.5, 1.5)
        ax.set_yticks([])
        ax.set_xlabel("Wavelength-squared (m**2)")
        span = lsq_max - lsq_min
        max_rm_scale = np.pi / span if span > 0 else float("nan")
        ax.text(0.02, 0.9, f"min Lsq: {lsq_min:.3f}", transform=ax.transAxes,
                ha="left", va="top", fontsize=7)
        ax.text(0.98, 0.9, f"max RM-scale: {max_rm_scale:.3f}", transform=ax.transAxes,
                ha="right", va="top", fontsize=7)
        ax.text(0.02, 0.1, f"Lsq_span: {span:.3f}", transform=ax.transAxes,
                ha="left", va="bottom", fontsize=7)

    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out_path = out_dir / f"{slug}.png"
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    return out_path


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    csv_dir = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)
    for slug, title in CASES:
        out_path = plot_case(csv_dir, out_dir, slug, title)
        print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
