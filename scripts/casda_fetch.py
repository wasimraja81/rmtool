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

Add --fetch to stage and download the I/Q/U/V cubes (one per
Stokes parameter found) from every band found Q/U-complete:

    python3 scripts/casda_fetch.py --target dancingghosts --username you@example.org \
        --fetch --outdir ./data

If a fetched Q/U cube doesn't already carry a per-channel BEAMS table of
its own (some SBIDs' headers have one, older/pilot-era ones may not --
checked directly per file, not guessed from vintage), --fetch also pulls
the SBID's diagnostics evaluation file from CASDA and curates the one
beamlog it contains into the ASCII format convolve_cubes' beamfiles=
expects, so there's always a target beam to convolve to even when the
FITS file alone doesn't carry one.
"""

from __future__ import annotations

import argparse
import os
import re
import tarfile
from dataclasses import dataclass, field
from typing import Optional
from xml.etree import ElementTree as ET

import numpy as np
import pyvo
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
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
        help="Actually request and download the I/Q/U/V cubes from every "
             "band found Q/U-complete (default: report only, fetch nothing)",
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


STOKES_AXIS_MAP = {1: "I", 2: "Q", 3: "U", 4: "V"}


def stokes_letter_from_header(path: str) -> str:
    """Stokes parameter of a downloaded FITS file, from its own WCS
    STOKES axis (CTYPE<n>='STOKES', CRVAL<n> -> 1=I/2=Q/3=U/4=V). CASDA's
    cutout service renames every output to its own generic
    cutout-<job>-imagecube-<id>.fits, discarding the original archive
    filename entirely (confirmed by inspecting downloaded files), so
    stokes_letter_from_filename's filename matching only works
    pre-fetch, on the query result's own filename column, never on a
    downloaded file -- this is the post-fetch equivalent."""
    header = fits.getheader(path)
    for axis in (3, 4, 1, 2):
        if str(header.get(f"CTYPE{axis}", "")).strip().upper() == "STOKES":
            crval = header.get(f"CRVAL{axis}")
            if crval is not None:
                return STOKES_AXIS_MAP.get(int(round(crval)), "UNKNOWN")
            break
    return "UNKNOWN"


def has_beams_table(path: str) -> bool:
    """Whether this downloaded FITS file already carries a CASA-style
    per-channel BEAMS binary table extension -- the ground truth for
    whether convolve_cubes' beamfiles=auto can read a usable beam
    straight from the file, as opposed to the CASAMBM header keyword,
    which is not a reliable proxy for it: confirmed against a
    SB74876 Q/U pair that the .conv cube's BEAMS table (when present)
    matches the observatory's own beamlog exactly on all 288 channels,
    but older/pilot-era headers can omit the CASAMBM keyword even when
    the underlying convolution was done, so its absence alone
    proves nothing either way -- the extension itself must be checked
    directly."""
    with fits.open(path) as hdul:
        return any(hdu.name == "BEAMS" for hdu in hdul)


def sbid_dirname(obs_id: str) -> str:
    """'ASKAP-74876' -> 'SB74876', matching the SB<n> convention already
    embedded in the filenames themselves. Falls back to the raw obs_id
    unchanged for anything that doesn't match the ASKAP-<digits> pattern,
    rather than guessing -- seen so far only for our target subtype/pol
    rows, not the project-code-style obs_ids (e.g. 'AS103') that show up
    for other, unfetched subtypes."""
    m = re.match(r"^ASKAP-(\d+)$", obs_id)
    return f"SB{m.group(1)}" if m else obs_id


def sbid_number(obs_id: str) -> Optional[int]:
    """Bare SBID integer from an obs_id like 'ASKAP-74876' -- None if it
    doesn't match, same convention/caveat as sbid_dirname."""
    m = re.match(r"^ASKAP-(\d+)$", obs_id)
    return int(m.group(1)) if m else None


STOKES_LETTERS = ("i", "q", "u", "v")


def stokes_letter_from_filename(filename: str) -> Optional[str]:
    """Lowercase i/q/u/v from the known archive filename convention
    (image.restored.<i|q|u|v>.<field>.SB<n>.contcube[.conv].fits), or
    None if it doesn't match -- only meaningful pre-fetch, on the query
    result's own filename column (see identify_from_header for why this
    can't be used on a downloaded file)."""
    m = re.search(r"\.restored\.([iquv])\.", filename, re.IGNORECASE)
    return m.group(1).lower() if m else None


def select_iquv_mask(rows) -> list:
    """Boolean mask selecting at most one file per Stokes parameter (I,
    Q, U, V) from one observation's own cont.restored.3d rows, using the
    known filename convention rather than pol_states alone -- pol_states
    alone can match more than one file per Stokes parameter for the same
    observation (confirmed against SB74876, a single-field observation
    that still returned 4 pol_states=/Q/,/U/ rows instead of the 2
    expected -- a raw/.conv duplicate pair for each of Q and U, not two
    separate fields).

    Q and U: .conv. is a HARD, ALL-OR-NOTHING requirement for the WHOLE
    observation, not a per-letter preference -- a plain (non-.conv) file
    hasn't had its 36 PAF beams matched to each other for that channel,
    so using it would silently combine beam-wise-mismatched data into
    the tomography input. If EITHER Q or U lacks a .conv. variant,
    NOTHING is selected for this observation at all -- not even I/V --
    rather than partially fetching whatever does qualify, since an
    observation missing .conv. Q/U isn't usable for tomography regardless
    of what else is available for it. The caller is expected to have
    already excluded such observations via Observation.has_conv_qu (see
    print_verdict), so hitting this path at all indicates a mismatch
    between that filtering and what's in `rows`; it's still
    enforced here directly rather than trusted to the caller, since a
    partial fetch would be actively misleading, not just redundant.

    I and V: .conv. is a soft preference (they're for source verification/
    noise-floor checks per the user's own stated use, not tomography
    itself) -- falls back to the plain file if no .conv. variant exists,
    but prints a loud warning, since Q/U being .conv while I/V aren't is
    an anomaly worth flagging rather than accepting quietly. This only
    ever applies once Q and U have already both been confirmed .conv,
    per the gate above."""
    by_letter: dict = {}
    for i, filename in enumerate(rows["filename"]):
        letter = stokes_letter_from_filename(str(filename))
        if letter:
            by_letter.setdefault(letter, []).append(i)

    def has_conv(letter: str) -> bool:
        return any(".conv." in str(rows["filename"][i])
                    for i in by_letter.get(letter, []))

    if not (has_conv("q") and has_conv("u")):
        print(f"    ERROR: Q and/or U lack a .conv. variant -- this whole "
              f"observation is not usable for tomography. Fetching NOTHING "
              f"for it, including I/V (this observation should not have "
              f"reached fetch at all -- check has_conv_qu filtering "
              f"upstream).")
        return [False] * len(rows)

    keep = set()
    for letter in STOKES_LETTERS:
        idxs = by_letter.get(letter, [])
        if not idxs:
            print(f"    Stokes {letter.upper()}: no file found for this observation")
            continue
        conv_idxs = [i for i in idxs if ".conv." in str(rows["filename"][i])]
        if conv_idxs:
            chosen = conv_idxs[0]
            print(f"    Stokes {letter.upper()}: {len(idxs)} file(s) matched, "
                  f"using .conv. variant: {rows['filename'][chosen]}")
            keep.add(chosen)
        else:
            # Only I/V can reach here -- Q/U conv-presence was already
            # confirmed by the gate above.
            chosen = idxs[0]
            print(f"    WARNING: Stokes {letter.upper()} has no .conv. "
                  f"variant even though Q/U do -- anomaly. Falling back to "
                  f"the plain (non-.conv) file: {rows['filename'][chosen]}")
            keep.add(chosen)
    return [i in keep for i in range(len(rows))]


DATALINK_NS = {"vot": "http://www.ivoa.net/xml/VOTable/v1.3"}


def find_diagnostics_datalink_url(sbid: int) -> Optional[str]:
    """The DataLink access_url for this SBID's diagnostics evaluation
    file -- the tarball with per-beam ASCII beamlog files inside, whose
    Target BMAJ/BMIN/BPA columns record the per-channel common beam
    this SBID's .conv cubes were beam-wise matched to. This table is
    NOT part of ivoa.obscore (confirmed: query_region only ever returns
    dataproduct_type cube/visibility for a position search) -- it was
    found only by listing every table CASDA's TAP service exposes.
    Public metadata, no login needed for this query; only the
    file download requires authentication."""
    tap = pyvo.dal.TAPService("https://casda.csiro.au/casda_vo_tools/tap")
    query = (
        "SELECT filename, access_url FROM casda.observation_evaluation_file "
        f"WHERE sbid = {int(sbid)} AND filename LIKE 'diagnostics%'"
    )
    result = tap.search(query)
    if len(result) == 0:
        return None
    return str(result[0]["access_url"])


def resolve_datalink_sync_url(datalink_url: str, auth) -> Optional[str]:
    """Follow one DataLink XML indirection to the sync?id=...
    download URL -- both ObsCore's and the evaluation-file table's own
    access_url values point at a DataLink XML document describing the
    file, not the file itself (confirmed directly: fetching an
    access_url as a tar fails with 'invalid header'). Picks the row
    with semantics='#this' and a plain sync endpoint, not the
    Pawsey-internal one (only reachable from inside Pawsey's network)."""
    resp = requests.get(datalink_url, auth=auth, timeout=60)
    resp.raise_for_status()
    root = ET.fromstring(resp.content)
    fields = [f.get("name") for f in root.findall(".//vot:FIELD", DATALINK_NS)]
    for tr in root.findall(".//vot:TR", DATALINK_NS):
        values = [td.text or "" for td in tr.findall("vot:TD", DATALINK_NS)]
        row = dict(zip(fields, values))
        url = row.get("access_url", "")
        if row.get("semantics") == "#this" and url and "/sync/pawsey" not in url:
            return url
    return None


def fetch_beamlog_text(sync_url: str, sbid: int) -> Optional[str]:
    """Stream the diagnostics tarball and return the text content of
    the first Stokes-Q .conv beamlog found, whichever of the 36 PAF
    beam numbers it happens to be -- confirmed against SB74876
    data that Target BMAJ/BMIN/BPA is identical across all 36 beams and
    across Stokes parameters (the restoring beam comes from the (u,v)
    coverage/imaging weights per channel, the same for every Stokes
    combination), so any one match is as good as any other. CASDA's
    sync download endpoint does not support HTTP Range requests, and
    its SODA service for evaluation files declares no sub-file
    parameter in its own DataLink inputParams (confirmed directly) --
    there is no way to ask CASDA for just this one tar member, so this
    streams from the start of the tar and stops as soon as a match is
    read, without ever writing the tar to disk or reading past the
    match. Returns None if the tar is exhausted with no match.

    No `auth` here, unlike resolve_datalink_sync_url's own fetch -- the
    id= token in sync_url is itself a scoped, pre-authorized credential;
    confirmed directly that adding HTTP Basic Auth on top of it makes
    CASDA reject the request with 401 (application/json error body)
    where the same URL with no auth at all returns 200 application/x-tar."""
    pattern = re.compile(rf"beamlog\.image\.restored\.q\..*SB{sbid}.*\.beam\d+\.conv\.txt$")
    resp = requests.get(sync_url, stream=True, timeout=600)
    resp.raise_for_status()
    try:
        with tarfile.open(fileobj=resp.raw, mode="r|*") as tar:
            for member in tar:
                if pattern.search(member.name):
                    fh = tar.extractfile(member)
                    if fh is not None:
                        return fh.read().decode("utf-8", errors="replace")
    finally:
        resp.close()
    return None


def curate_beamlog(raw_text: str, out_path: str) -> int:
    """Reformat a CASDA beamlog's 0-indexed 'Target BMAJ/BMIN/BPA'
    columns (the per-channel beam this SBID's .conv cubes were already
    matched to, beam-wise across the 36 PAF beams) into convolve_cubes'
    own beamfiles= ASCII format: 1-indexed "channel bmaj_arcsec
    bmin_arcsec bpa_deg" (src/convolve_cubes.f90's read_beams_ascii,
    see also cfg/example_beamLog.txt). Returns the number of channel
    rows written."""
    n = 0
    with open(out_path, "w") as out:
        out.write(
            "# Target BMAJ/BMIN/BPA curated from a CASDA diagnostics "
            "beamlog by scripts/casda_fetch.py.\n"
            "# channel bmaj_arcsec bmin_arcsec bpa_deg (1-indexed)\n"
        )
        for line in raw_text.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            ch0, bmaj, bmin, bpa = int(parts[0]), parts[1], parts[2], parts[3]
            out.write(f"{ch0 + 1} {bmaj} {bmin} {bpa}\n")
            n += 1
    return n


def fetch_and_curate_beamlog(obs_id: str, auth, obs_dir: str) -> Optional[str]:
    """Find, stream, and curate the one beamlog this observation needs,
    for when its .conv Q/U cubes don't already carry a usable BEAMS
    table of their own (see has_beams_table). Prints progress
    throughout since streaming a multi-GB tar can take minutes; never
    raises -- returns the curated file's path, or None (with a printed
    reason) if anything along the way didn't turn up what was needed,
    so one observation's beamlog trouble doesn't take down the rest of
    a multi-observation fetch."""
    sbid = sbid_number(obs_id)
    if sbid is None:
        print(f"    could not parse an SBID number out of {obs_id!r} -- "
              f"skipping beamlog fetch")
        return None
    try:
        print(f"    looking up {obs_id}'s diagnostics evaluation file...")
        datalink_url = find_diagnostics_datalink_url(sbid)
        if datalink_url is None:
            print(f"    no diagnostics evaluation file published for "
                  f"SB{sbid} -- cannot curate a beamlog for it")
            return None
        sync_url = resolve_datalink_sync_url(datalink_url, auth)
        if sync_url is None:
            print(f"    could not resolve a download URL from the "
                  f"diagnostics file's DataLink response -- cannot "
                  f"curate a beamlog")
            return None
        print(f"    streaming the diagnostics tarball for its beamlog "
              f"(this can take a few minutes -- stops as soon as a "
              f"match is found)...")
        raw_text = fetch_beamlog_text(sync_url, sbid)
        if raw_text is None:
            print(f"    no matching beamlog found in the diagnostics "
                  f"tarball for SB{sbid}")
            return None
    except Exception as exc:
        print(f"    FAILED to fetch beamlog: {type(exc).__name__}: {exc}")
        return None

    out_path = os.path.join(obs_dir, "beamlog_target.txt")
    n = curate_beamlog(raw_text, out_path)
    print(f"    wrote {out_path} ({n} channels) -- pass this as "
          f"convolve_cubes' beamfiles= for this observation's Q/U cubes")
    return out_path


@dataclass
class Observation:
    obs_id: str
    freq: Optional[tuple]
    pol_counts: dict = field(default_factory=dict)
    conv_letters: set = field(default_factory=set)

    @property
    def has_qu(self) -> bool:
        return (any("Q" in p for p in self.pol_counts) and
                any("U" in p for p in self.pol_counts))

    @property
    def has_conv_qu(self) -> bool:
        """Whether a .conv. (common-beam) variant is available for BOTH
        Q and U specifically -- known pre-fetch, from the filename alone,
        no download needed. This is what --fetch will select
        for Q/U (select_iquv_mask always prefers .conv. when present),
        so it's a more accurate preview of fetch-time behaviour than
        has_qu alone."""
        return "q" in self.conv_letters and "u" in self.conv_letters


def summarise_observations(table) -> list:
    """One Observation per distinct obs_id in table, with its frequency
    range, a count of each pol_states value seen, and which Stokes
    letters have a .conv. filename variant available."""
    by_obs: dict = {}
    for row in table:
        obs_id = str(row["obs_id"])
        obs = by_obs.setdefault(obs_id, Observation(obs_id, None))
        freq = freq_range_mhz(row)
        if freq is not None:
            obs.freq = freq
        pol = str(row["pol_states"])
        obs.pol_counts[pol] = obs.pol_counts.get(pol, 0) + 1
        filename = str(row["filename"])
        letter = stokes_letter_from_filename(filename)
        if letter and ".conv." in filename:
            obs.conv_letters.add(letter)
    return sorted(by_obs.values(), key=lambda o: o.freq or (0.0, 0.0))


def print_observation_table(observations: list) -> None:
    print(f"\n{'obs_id':<16}{'freq range (MHz)':<20}"
          f"{'pol_states seen (count)':<45}{'Q & U?':<8}{'.conv Q&U?'}")
    for obs in observations:
        freq_str = f"{obs.freq[0]:.1f}-{obs.freq[1]:.1f}" if obs.freq else "unknown"
        pol_str = ", ".join(f"{p}:{n}" for p, n in sorted(obs.pol_counts.items()))
        print(f"{obs.obs_id:<16}{freq_str:<20}{pol_str:<45}"
              f"{'yes' if obs.has_qu else 'no':<8}"
              f"{'yes' if obs.has_conv_qu else 'no'}")


def print_verdict(observations: list) -> list:
    """Prints the bottom-line answer and returns the list of Observations
    usable for tomography. Gates on has_conv_qu, NOT
    has_qu -- a hard requirement, not just Q+U present in any form. A
    plain (non-.conv) Q/U pair hasn't had its 36 PAF beams matched to
    each other for that channel; using it would silently combine
    beam-wise-mismatched data, so it does not count as usable here even
    though the observation table above still shows it under "Q & U?" for
    diagnostic visibility (compare that column against ".conv Q&U?" to
    see exactly why an observation was excluded)."""
    usable = [o for o in observations if o.has_conv_qu and o.freq is not None]
    bands = sorted(set(o.freq for o in usable))
    print(f"\n{len(usable)} of {len(observations)} observation(s) have "
          f"both Q and U present AS .conv (beam-wise matched).")
    if len(bands) >= 2:
        print(f"{len(bands)} distinct frequency band(s) with a complete "
              f".conv Q/U pair -- multi-band tomography looks possible:")
        for lo, hi in bands:
            print(f"  {lo:.1f}-{hi:.1f} MHz")
    elif len(bands) == 1:
        lo, hi = bands[0]
        print(f"Only 1 distinct band with a complete .conv Q/U pair "
              f"({lo:.1f}-{hi:.1f} MHz) -- single-band synthesis only "
              f"for this position with currently public data.")
    else:
        print("No band with both Q and U present as .conv -- this "
              "position isn't usable yet with what CASDA currently has "
              "public.")
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
    distinct_bands = len(set(o.freq for o in usable))
    print(f"\nFetching {len(usable)} observation(s) across {distinct_bands} "
          f"distinct band(s) into {args.outdir}/<SBID>/, cutout "
          f"radius={cutout_radius} arcmin...")

    # One subdirectory per observation (one SBID) -- NOT necessarily one
    # per distinct band, since more than one observation/SBID can share
    # the same frequency band (observed directly: 3 Q+U-complete
    # observations for PKS 2130-538, all at the same 799.5-1087.5 MHz
    # band).
    #
    # All of an observation's (up to 4) selected files are requested in
    # ONE batched cutout() call, not one call per file -- CASDA processes
    # them together server-side, which matters for wall time (a single
    # ~1GB cutout has been observed to take upwards of a minute; doing
    # that 4 separate times, with 4 separate job-submit-and-poll cycles,
    # would multiply it for no benefit once file identification no longer
    # depends on request order -- see below). Each OBSERVATION is still
    # fetched independently of the others, and failures don't abort the
    # rest -- a mid-download network drop on one SBID (observed directly:
    # requests.exceptions.ConnectionError partway through a ~1GB transfer)
    # would otherwise take an unhandled exception all the way up and kill
    # the whole run, discarding every observation not yet reached even
    # though they're independent CASDA jobs. The trade-off from batching
    # within one observation: a network drop partway through now risks
    # that observation's whole batch rather than just one file -- accepted
    # for the wall-time win, since select_iquv_mask has already pared the
    # batch down to at most 4 well-defined files, not a large redundant
    # set.
    #
    # Identifying which downloaded file is which is NOT done from
    # cutout()'s response order (not a documented guarantee to match the
    # input row order under async job completion) and NOT from the
    # downloaded filename (CASDA's cutout service renames every output to
    # its own generic cutout-<job>-imagecube-<id>.fits, discarding the
    # original name entirely -- confirmed directly). Instead, each
    # downloaded file's own header is read after the fact
    # (stokes_letter_from_header) -- authoritative, not inferred.
    obs_results = []
    for obs in usable:
        mask = [str(v) == obs.obs_id for v in subset["obs_id"]]
        rows = subset[mask]
        mask_pol = [str(v) in ("/I/", "/Q/", "/U/", "/V/") for v in rows["pol_states"]]
        rows = rows[mask_pol]

        obs_dir = os.path.join(args.outdir, sbid_dirname(obs.obs_id))
        os.makedirs(obs_dir, exist_ok=True)
        freq_str = f"{obs.freq[0]:.1f}-{obs.freq[1]:.1f} MHz" if obs.freq else "unknown freq"
        print(f"\n{obs.obs_id} ({freq_str}) -> {obs_dir}/")
        rows = rows[select_iquv_mask(rows)]
        print(f"  {len(rows)} file(s) selected (one per Stokes parameter found, "
              f".conv preferred over plain restored where both exist)")

        if len(rows) == 0:
            print(f"  SKIPPED: Q and/or U lack a .conv. variant -- nothing "
                  f"fetched for this observation (see ERROR above)")
            obs_results.append((obs.obs_id, "SKIPPED"))
            continue

        # Known pre-fetch, from the archive filename (same signal
        # select_iquv_mask already used to pick these rows) -- not
        # re-derived from the download, since the download's own header
        # signal (CASAMBM) is not a reliable proxy for it (see
        # has_beams_table).
        conv_by_letter = {
            stokes_letter_from_filename(str(fn)).upper(): ".conv." in str(fn)
            for fn in rows["filename"]
        }

        n_ok = 0
        n_total = len(rows)
        written: dict = {}
        try:
            urls = casda.cutout(rows, coordinates=coord, radius=cutout_radius * u.arcmin)
            downloaded = casda.download_files(urls, savedir=obs_dir)
            fits_paths = [p for p in downloaded if not p.endswith(".checksum")]
            for path in fits_paths:
                letter = stokes_letter_from_header(path)
                tag = f"{letter}_conv" if conv_by_letter.get(letter) else letter
                new_fits = os.path.join(obs_dir, f"{tag}.fits")
                os.replace(path, new_fits)
                print(f"  wrote {new_fits}")
                checksum_src = path + ".checksum"
                if os.path.exists(checksum_src):
                    os.replace(checksum_src, new_fits + ".checksum")
                written[letter] = new_fits
                n_ok += 1
        except Exception as exc:
            print(f"  FAILED: {type(exc).__name__}: {exc}")
            print(f"  (likely a transient network/server issue -- other "
                  f"observations are unaffected; retry {obs.obs_id} alone "
                  f"later if needed)")
        if n_ok == n_total:
            status = "OK"
        elif n_ok == 0:
            status = "FAILED"
        else:
            status = f"PARTIAL ({n_ok}/{n_total} files)"
        obs_results.append((obs.obs_id, status))

        if "Q" in written and "U" in written:
            if has_beams_table(written["Q"]) and has_beams_table(written["U"]):
                print(f"  Q and U both already carry a per-channel BEAMS "
                      f"table -- convolve_cubes' beamfiles=auto will read "
                      f"it directly, no separate beamlog needed.")
            else:
                print(f"  Q and/or U's .conv cube has no BEAMS table of "
                      f"its own -- fetching the observatory's beamlog to "
                      f"curate one:")
                fetch_and_curate_beamlog(obs.obs_id, casda._auth, obs_dir)

    print("\nFetch summary:")
    for obs_id, status in obs_results:
        print(f"  {obs_id}: {status}")

    print(
        "\nNEXT STEP -- do not skip: a band whose Q/U came back tagged "
        "_conv already has one common beam within that band, but "
        "different bands still very likely have different beams from "
        "each other; a band that fell back to the plain (non-_conv) "
        "file still has a channel-varying native beam within itself "
        "too. Either way, run convolve_cubes across the Q and U cubes "
        "together (all bands) before rm_synthesis -- the I/V cubes are "
        "for visual source-verification and noise-floor checks and "
        "don't need the same treatment. If an SBID directory has its own "
        "beamlog_target.txt, that Q/U cube had no BEAMS table of its own "
        "-- pass beamfiles=<path to it> for that pair; otherwise "
        "beamfiles=auto (the default) already reads the beam straight "
        "from the FITS file. Skipping this will not fail the run "
        "outright (resolution matching is a warning, not a hard "
        "requirement -- see APP_REFERENCE.md's reference_band callout) but "
        "will silently weaken the result. See QUICKSTART.md section 4a."
    )


if __name__ == "__main__":
    main()
