#!/usr/bin/env python3
"""
Check the Sec 10 thesis-grounded multi-band scenario
(planning/MULTI_BAND_TOMOGRAPHY_PLAN.md Sec 10; Raja 2014 Table 6.1/6.2):
point source + Faraday-thick top-hat + F2/F3 pair, recovered from P-band
alone, L-band alone, and the P+L combined synthesis.

Usage:
  check_thesis_scenario.py <amp_P.fits> <amp_L.fits> <amp_PL.fits> <truth.json>
"""
import sys, json
from pathlib import Path
import numpy as np
from astropy.io import fits


def load_spectrum(amp_path, x, y):
    with fits.open(amp_path) as hdul:
        data = hdul[0].data.squeeze()  # (nrm, ny, nx)
        hdr = hdul[0].header
    crval3 = float(hdr["CRVAL3"])
    cdelt3 = float(hdr["CDELT3"])
    nrm = data.shape[0]
    rm_axis = crval3 + np.arange(nrm) * cdelt3
    return rm_axis, data[:, y, x]


def window(rm_axis, spectrum, lo, hi):
    mask = (rm_axis >= lo) & (rm_axis <= hi)
    return rm_axis[mask], spectrum[mask]


def peak_rm(rm_axis, spectrum):
    i = int(np.argmax(spectrum))
    return float(rm_axis[i]), float(spectrum[i])


def local_stat(rm_axis, spectrum, rm_centre, half_width, stat=np.max):
    mask = (rm_axis >= rm_centre - half_width) & (rm_axis <= rm_centre + half_width)
    vals = spectrum[mask]
    return float(stat(vals)) if len(vals) else 0.0


def targeted_resolved(rm_axis, spectrum, rm_a, rm_b, half_width=4.0, dip_frac=0.7):
    """Resolved-vs-blended check targeted at two *known* expected peak
    positions (rather than blind local-maxima counting, which is thrown
    off by unrelated dirty-beam sidelobe ripple elsewhere in the window --
    see the P-alone/PL false-positive-peak note in Sec 10's Implementation
    notes). 'Resolved' means the local max near the midpoint is clearly
    lower than the local maxima at rm_a and rm_b themselves."""
    peak_a = local_stat(rm_axis, spectrum, rm_a, half_width)
    peak_b = local_stat(rm_axis, spectrum, rm_b, half_width)
    mid = 0.5 * (rm_a + rm_b)
    dip = local_stat(rm_axis, spectrum, mid, half_width)
    resolved = dip < dip_frac * min(peak_a, peak_b)
    return resolved, peak_a, peak_b, dip


def main():
    if len(sys.argv) < 5:
        print("Usage: check_thesis_scenario.py <amp_P> <amp_L> <amp_PL> <truth.json>")
        sys.exit(1)
    amp_p, amp_l, amp_pl, truth_path = sys.argv[1:5]
    truth = json.loads(Path(truth_path).read_text())
    x, y = truth["src_x"], truth["src_y"]
    point_rm = truth["point_rm"]
    thick_lo, thick_hi = truth["thick_rm_lo"], truth["thick_rm_hi"]
    f2_rm, f3_rm = truth["f2_rm"], truth["f3_rm"]

    ok = True

    def check(cond, label):
        nonlocal ok
        print(f"[{'OK' if cond else 'FAIL'}] {label}")
        if not cond:
            ok = False

    rm_p, spec_p = load_spectrum(amp_p, x, y)
    rm_l, spec_l = load_spectrum(amp_l, x, y)
    rm_pl, spec_pl = load_spectrum(amp_pl, x, y)

    # --- Point source: recovered accurately at P alone and P+L combined.
    # (Not asserted at L alone -- delta_RM=250 there is coarser than the
    # thick component's own span, so the point source and thick component
    # blend into one feature; that blending is itself checked below.)
    for label, rm_axis, spec in (("P-alone", rm_p, spec_p), ("P+L", rm_pl, spec_pl)):
        pk_rm, _ = peak_rm(rm_axis, spec)
        cdelt = abs(rm_axis[1] - rm_axis[0])
        tol = max(2.0 * cdelt, 5.0)
        check(abs(pk_rm - point_rm) <= tol,
              f"{label}: point source recovered near RM={point_rm:+.1f} "
              f"(found {pk_rm:+.2f}, tol={tol:.2f})")

    # --- Thick top-hat component: washed out at P alone, revealed
    # (substantially more recovered signal) when combined with L.
    _, thick_win_p = window(rm_p, spec_p, thick_lo, thick_hi)
    _, thick_win_pl = window(rm_pl, spec_pl, thick_lo, thick_hi)
    thick_peak_p = float(np.max(thick_win_p)) if len(thick_win_p) else 0.0
    thick_peak_pl = float(np.max(thick_win_pl)) if len(thick_win_pl) else 0.0
    check(thick_peak_pl > 2.0 * thick_peak_p,
          f"Thick component revealed by combining bands: peak in "
          f"[{thick_lo:.0f},{thick_hi:.0f}] rad/m^2 is {thick_peak_pl:.3f} "
          f"(P+L) vs {thick_peak_p:.3f} (P alone), expected P+L > 2x P-alone")

    # --- F2/F3 close pair: unresolved at L alone, resolved at P alone --
    # both checked with a dip-vs-peaks comparison targeted at the known
    # expected RM positions (robust against unrelated dirty-beam sidelobe
    # ripple elsewhere in the window, unlike blind peak-counting).
    resolved_p, pa, pb, dip = targeted_resolved(rm_p, spec_p, f2_rm, f3_rm)
    check(resolved_p,
          f"F2/F3 resolved at P-alone (peaks {pa:.2f}/{pb:.2f}, dip {dip:.2f})")
    resolved_l, pa, pb, dip = targeted_resolved(rm_l, spec_l, f2_rm, f3_rm)
    check(not resolved_l,
          f"F2/F3 blended at L-alone (peaks {pa:.2f}/{pb:.2f}, dip {dip:.2f})")

    # --- F2/F3 at P+L combined: NOT asserted as cleanly "resolved" here.
    # This codebase has no RM-CLEAN deconvolution (grep confirms no
    # RMCLEAN/rm-clean anywhere in src/) -- the raw dirty spectrum from
    # combining two widely-separated bands (P+L) carries fine-period
    # sidelobe ringing from the large total lambda^2 span (thesis Sec
    # 2.5's "dirty RM response... side-lobe levels" discussion), whose
    # amplitude at 2 rad/m^2 sampling can rival the true F2/F3 dip even
    # though the underlying resolution (delta_RM=14.7) is fine enough in
    # principle. A simple dip-vs-peak comparison on the un-cleaned
    # spectrum is not a reliable test of this specific claim -- this is a
    # genuine, documented limitation (see Sec 10 addendum / this ticket's
    # Evidence section), not something this check papers over. What *is*
    # checked: both expected positions still carry a real, elevated
    # signal (confirming the merge computed something meaningful there,
    # not a garbage/near-zero DFT result).
    peak_a = local_stat(rm_pl, spec_pl, f2_rm, 4.0)
    peak_b = local_stat(rm_pl, spec_pl, f3_rm, 4.0)
    # Absolute floor, not relative to the window median: the combined
    # spectrum's ringing (see comment above) inflates the local median
    # enough that a relative comparison isn't meaningful here. Both true
    # components have injected amplitude 8; a real detection should sit
    # well above the injected per-channel noise level (sigma=0.005).
    check(peak_a > 2.0 and peak_b > 2.0,
          f"F2/F3 both still show elevated signal at P+L combined "
          f"(peak_a={peak_a:.2f}, peak_b={peak_b:.2f}, expected >2.0) -- "
          f"NOT a resolved-vs-blended claim, see comment above")

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
