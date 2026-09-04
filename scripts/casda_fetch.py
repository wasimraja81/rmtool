#!/home/wasim/venv/rmtool/bin/python3
"""Check CASDA for calibrated POSSUM Stokes Q/U cubes usable by this
project's RM-synthesis tomography, for a given sky position -- and, once
confirmed, fetch them.

Requires astroquery (not in this project's core requirements -- install
once into the project venv: ~/venv/rmtool/bin/pip install astroquery).

Password handling: astroquery's login() takes no password argument at all
-- there is no way to pass one via CLI arg or env var. The first run below
prompts for it interactively and saves it to the OS keyring
(store_password=True); every later run reads it back from the keyring
silently, with no prompt, which is what makes an unattended pipeline run
possible after that first interactive use.

Usage:

    python3 scripts/casda_fetch.py --target dancingghosts --username you@example.org
    python3 scripts/casda_fetch.py --target "Kes 27" --username you@example.org
    python3 scripts/casda_fetch.py --ra 323.56667 --dec -53.61528 --radius 8 \
        --username you@example.org

--target takes either of two kinds of value, checked in this order:
(1) one of this project's own presets -- msh15-56, kes27, cena,
dancingghosts -- objects investigated for this project, see
docs/dev/CASDA_FETCH_PLAN.md section 4 for citations and how each was
found; or (2) if not one of those, any other name resolvable via astropy
SkyCoord.from_name (SIMBAD/NED-style community naming, via Sesame/CDS).
Either way, the resolved coordinates are looked up live at run time, not
hardcoded here, so a stale/wrong coordinate can't silently creep in --
and the script prints which of the two paths it used, so that's never
ambiguous after the fact.

The one question this tool answers: for this position, how many
observations exist, at how many distinct frequency bands, with a
complete Stokes Q and U pair under dataproduct_subtype=cont.restored.3d
-- confirmed against Kes 27's own CASDA data to be the calibrated,
per-channel, science-ready cube product, as opposed to cont.cleanmodel/residual/
weight.3d (intermediate CLEAN products) or the .t0/.t1 Taylor-term
variants (frequency-collapsed MFS images, not usable for RM synthesis).
2+ such bands means multi-band tomography is possible; exactly 1 means
single-band only; 0 means this position isn't usable yet with what CASDA
currently has public.

If nothing matches, this tool falls back to showing what collections and
subtypes DO exist nearby, so "nothing found" is diagnosable rather than a
dead end.

Add --fetch to actually stage and download the Q/U cubes from every
band found complete:

    python3 scripts/casda_fetch.py --target dancingghosts --username you@example.org \
        --fetch --outdir ./data
"""

from __future__ import annotations

import argparse
import os
import re
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.casda import Casda

DEFAULT_PROJECT = "POSSUM"
DEFAULT_SUBTYPE = "cont.restored.3d"

PRESET_TARGETS = {
    "msh15-56": "MSH 15-56",
    "kes27": "Kes 27",
    "cena": "Centaurus A",
    "dancingghosts": "PKS 2130-538",
}

SPEED_OF_LIGHT_M_PER_S = 299792458.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--username", required=True,
        help="OPAL account email (password is read from the OS keyring -- "
             "see the password-handling note in this script's own module "
             "docstring)",
    )
    target = parser.add_mutually_exclusive_group(required=True)
    target.add_argument(
        "--target",
        help="Object to search around -- checked against this project's "
             "own preset list first (see module docstring), then falls "
             "back to community name resolution (SkyCoord.from_name) if "
             "not found there. Prints which path was used.",
    )
    target.add_argument(
        "--ra", type=float,
        help="RA in degrees (use together with --dec)",
    )
    parser.add_argument(
        "--dec", type=float, default=None,
        help="Dec in degrees (required if --ra is given)",
    )
    parser.add_argument(
        "--radius", type=float, default=15.0,
        help="Search/cutout radius in arcmin (default: 15)",
    )
    parser.add_argument(
        "--project", default=DEFAULT_PROJECT,
        help=f"Match rows whose obs_collection contains this substring "
             f"(default: {DEFAULT_PROJECT}). Give an empty string to see "
             f"every collection instead of restricting to one.",
    )
    parser.add_argument(
        "--subtype", default=DEFAULT_SUBTYPE,
        help=f"Match rows whose dataproduct_subtype contains this "
             f"substring (default: {DEFAULT_SUBTYPE}, confirmed to be the "
             f"calibrated per-channel science cube -- see module docstring). "
             f"Give an empty string to see every subtype instead.",
    )
    parser.add_argument(
        "--cutout-radius", type=float, default=None,
        help="Radius in arcmin for the downloaded cutout (with "
             "--fetch) -- independent of --radius, which only controls "
             "the search. Defaults to --radius's value if not given.",
    )
    parser.add_argument(
        "--fetch", action="store_true",
        help="Actually request and download the Q/U cubes from every "
             "band found complete (default: report only, fetch nothing)",
    )
    parser.add_argument(
        "--outdir", default=".",
        help="Directory to save downloaded cutouts into (with --fetch)",
    )
    args = parser.parse_args()
    if args.ra is not None and args.dec is None:
        parser.error("--ra given without --dec")
    return args


def resolve_target(args: argparse.Namespace) -> SkyCoord:
    if args.target:
        if args.target in PRESET_TARGETS:
            name = PRESET_TARGETS[args.target]
            print(f"'{args.target}' matched this project's own preset "
                  f"list -> resolving '{name}' (internal lookup)...")
            return SkyCoord.from_name(name)
        print(f"'{args.target}' is not one of this project's own presets "
              f"({', '.join(sorted(PRESET_TARGETS))}) -- resolving via "
              f"community name lookup (SkyCoord.from_name, Sesame/CDS)...")
        return SkyCoord.from_name(args.target)
    return SkyCoord(ra=args.ra * u.deg, dec=args.dec * u.deg)


def freq_range_mhz(row) -> Optional[tuple]:
    """(freq_lo_mhz, freq_hi_mhz) from ObsCore's em_min/em_max (metres,
    wavelength -- so em_max, the LONGER wavelength, is the LOWER
    frequency edge). None if either value is missing/masked, rather than
    letting float() silently turn a masked value into nan."""
    em_min, em_max = row["em_min"], row["em_max"]
    if np.ma.is_masked(em_min) or np.ma.is_masked(em_max):
        return None
    if em_min is None or em_max is None:
        return None
    return (
        SPEED_OF_LIGHT_M_PER_S / float(em_max) / 1e6,
        SPEED_OF_LIGHT_M_PER_S / float(em_min) / 1e6,
    )


def sbid_dirname(obs_id: str) -> str:
    """'ASKAP-74876' -> 'SB74876', matching the SB<n> convention already
    embedded in the filenames themselves. Falls back to the raw obs_id
    unchanged for anything that doesn't match the ASKAP-<digits> pattern,
    rather than guessing -- seen so far only for our target subtype/pol
    rows, not the project-code-style obs_ids (e.g. 'AS103') that show up
    for other, unfetched subtypes."""
    m = re.match(r"^ASKAP-(\d+)$", obs_id)
    return f"SB{m.group(1)}" if m else obs_id


@dataclass
class Observation:
    obs_id: str
    freq: Optional[tuple]
    pol_counts: dict = field(default_factory=dict)

    @property
    def has_qu(self) -> bool:
        return (any("Q" in p for p in self.pol_counts) and
                any("U" in p for p in self.pol_counts))


def summarise_observations(table) -> list:
    """One Observation per distinct obs_id in table, with its frequency
    range and a count of each pol_states value seen."""
    by_obs: dict = {}
    for row in table:
        obs_id = str(row["obs_id"])
        obs = by_obs.setdefault(obs_id, Observation(obs_id, None))
        freq = freq_range_mhz(row)
        if freq is not None:
            obs.freq = freq
        pol = str(row["pol_states"])
        obs.pol_counts[pol] = obs.pol_counts.get(pol, 0) + 1
    return sorted(by_obs.values(), key=lambda o: o.freq or (0.0, 0.0))


def print_observation_table(observations: list) -> None:
    print(f"\n{'obs_id':<16}{'freq range (MHz)':<20}"
          f"{'pol_states seen (count)':<45}{'Q & U?'}")
    for obs in observations:
        freq_str = f"{obs.freq[0]:.1f}-{obs.freq[1]:.1f}" if obs.freq else "unknown"
        pol_str = ", ".join(f"{p}:{n}" for p, n in sorted(obs.pol_counts.items()))
        print(f"{obs.obs_id:<16}{freq_str:<20}{pol_str:<45}"
              f"{'yes' if obs.has_qu else 'no'}")


def print_verdict(observations: list) -> list:
    """Prints the bottom-line answer and returns the list of Observations
    that are actually usable (Q+U complete, known frequency)."""
    usable = [o for o in observations if o.has_qu and o.freq is not None]
    bands = sorted(set(o.freq for o in usable))
    print(f"\n{len(usable)} of {len(observations)} observation(s) have "
          f"both Q and U present.")
    if len(bands) >= 2:
        print(f"{len(bands)} distinct frequency band(s) with a complete "
              f"Q/U pair -- multi-band tomography looks possible:")
        for lo, hi in bands:
            print(f"  {lo:.1f}-{hi:.1f} MHz")
    elif len(bands) == 1:
        lo, hi = bands[0]
        print(f"Only 1 distinct band with a complete Q/U pair "
              f"({lo:.1f}-{hi:.1f} MHz) -- single-band synthesis only "
              f"for this position with currently public data.")
    else:
        print("No band with both Q and U found -- this position isn't "
              "usable yet with what CASDA currently has public.")
    return usable


def print_fallback_diagnostics(table) -> None:
    print("\nNothing matched --project/--subtype. Here's what IS public "
          "near this position, to diagnose why:")
    for column in ("obs_collection", "dataproduct_subtype"):
        values = sorted(set(str(v) for v in table[column]))
        print(f"\nDistinct {column} values ({len(values)}):")
        for v in values:
            n = sum(1 for row in table[column] if str(row) == v)
            print(f"  {v!r}: {n} row(s)")


def main() -> None:
    args = parse_args()
    coord = resolve_target(args)
    print(
        f"Target position: {coord.to_string('hmsdms')} "
        f"({coord.ra.deg:.5f}, {coord.dec.deg:.5f} deg)"
    )

    casda = Casda()
    casda.login(username=args.username, store_password=True)

    print(f"Querying CASDA within {args.radius} arcmin...")
    result = Casda.query_region(coord, radius=args.radius * u.arcmin)
    public = Casda.filter_out_unreleased(result) if len(result) else result
    print(f"{len(result)} product(s) found, {len(public)} publicly released.")
    if len(public) == 0:
        return

    subset = public
    if args.project:
        mask = [args.project in str(v) for v in subset["obs_collection"]]
        subset = subset[mask]
    if args.subtype:
        mask = [args.subtype in str(v) for v in subset["dataproduct_subtype"]]
        subset = subset[mask]

    print(f"\n{len(subset)} row(s) match --project={args.project!r} "
          f"--subtype={args.subtype!r}.")
    if len(subset) == 0:
        print_fallback_diagnostics(public)
        return

    observations = summarise_observations(subset)
    print_observation_table(observations)
    usable = print_verdict(observations)

    if not args.fetch:
        print("\n(discovery only -- pass --fetch to download)")
        return

    if not usable:
        print("\nNothing usable to fetch.")
        return

    cutout_radius = args.cutout_radius if args.cutout_radius is not None else args.radius
    print(f"\nFetching {len(usable)} band(s) into {args.outdir}/<SBID>/, "
          f"cutout radius={cutout_radius} arcmin...")

    # One subdirectory per observation (= one SBID = one frequency band),
    # not one flat dump -- download_files() otherwise names files only by
    # their own FITS filename, which doesn't group by band on its own.
    for obs in usable:
        mask = [str(v) == obs.obs_id for v in subset["obs_id"]]
        rows = subset[mask]
        mask_qu = [str(v) in ("/Q/", "/U/") for v in rows["pol_states"]]
        rows = rows[mask_qu]

        obs_dir = os.path.join(args.outdir, sbid_dirname(obs.obs_id))
        os.makedirs(obs_dir, exist_ok=True)
        freq_str = f"{obs.freq[0]:.1f}-{obs.freq[1]:.1f} MHz" if obs.freq else "unknown freq"
        print(f"\n{obs.obs_id} ({freq_str}), {len(rows)} file(s) -> {obs_dir}/")
        urls = casda.cutout(rows, coordinates=coord, radius=cutout_radius * u.arcmin)
        files = casda.download_files(urls, savedir=obs_dir)
        for path in files:
            print(f"  wrote {path}")

    print(
        "\nNEXT STEP -- do not skip: these are raw archival cubes, not yet "
        "resolution-matched. Each band's own channels almost certainly have "
        "a channel-varying native beam, and different bands will have "
        "different resolutions from each other regardless. Run "
        "convolve_cubes across ALL of these files together (all bands, "
        "both Q and U) before rm_synthesis -- skipping this will not fail "
        "the run outright (resolution matching is a warning, not a hard "
        "requirement -- see APP_REFERENCE.md's reference_band callout) but "
        "will silently weaken the result. See QUICKSTART.md section 4a."
    )


if __name__ == "__main__":
    main()
