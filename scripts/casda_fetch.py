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

A position with observation history routinely has more than one SBID at
the same frequency band (confirmed: PKS 2130-538 has three separate,
fully usable SBIDs at 799.5-1087.5 MHz alone) and sometimes more than one
reprocessed version of the same SBID's cube (confirmed: a ".vN."
filename tag) -- both need a decision, not a guess. --run-mode controls
how:

    --run-mode=dry (default): report only, fetch nothing. Numbers every
        distinct frequency band and lists every candidate SBID per band
        (latest observation first) plus any per-file version ambiguity --
        enough detail to build a --select list without guessing.

    --run-mode=auto: fetch + automatically use every band found, the
        latest-observed SBID within each, and the latest reprocessed
        version of each file -- for when you just want the most complete
        result CASDA currently supports, no manual choices.

    --run-mode=select --select="<list>": fetch + exactly what you name.
        <list> is a comma-separated "band[:sbid[:version]]" list (band
        numbers from a prior dry run) -- e.g. --select "1,2:51797,3:74876:v2"
        processes bands 1/2/3, using SBID 51797 for band 2 and, for band
        3's SBID 74876, requesting version 2 of whichever Stokes file(s)
        have one (Q and U are reprocessed independently -- a Stokes
        parameter without its own v2 falls back to its own latest
        version, not treated as an error); omitting sbid/version
        defaults that axis to latest. A single band number alone is single-band
        mode.

    python3 scripts/casda_fetch.py --target dancingghosts --username you@example.org
    python3 scripts/casda_fetch.py --target dancingghosts --username you@example.org \
        --run-mode=auto --outdir ./data
    python3 scripts/casda_fetch.py --target dancingghosts --username you@example.org \
        --run-mode=select --select "1,2:51797" --outdir ./data

If a fetched Q/U cube doesn't already carry a per-channel BEAMS table of
its own (some SBIDs' headers have one, older/pilot-era ones may not --
checked directly per file, not guessed from vintage), fetching also pulls
the SBID's diagnostics evaluation file from CASDA and curates the one
beamlog it contains into the ASCII format convolve_cubes' beamfiles=
expects, so there's always a target beam to convolve to even when the
FITS file alone doesn't carry one. A ready-to-run pipeline cfg (chaining
match_cubes into rm_synthesis) is generated from whatever was fetched.
"""

from __future__ import annotations

import argparse
import math
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
from astropy.time import Time
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
        "--run-mode", choices=("dry", "auto", "select"), default="dry",
        help="dry (default): report only, fetch nothing. auto: fetch "
             "using every band found, the latest-observed SBID per band, "
             "and the latest reprocessed version of each file. select: "
             "fetch exactly what --select names. See module docstring "
             "for the full picture.",
    )
    parser.add_argument(
        "--select",
        help="Required with --run-mode=select: comma-separated "
             "'band[:sbid[:version]]' list, e.g. '1,2:51797,3:74876:v2' "
             "-- band numbers come from a prior --run-mode=dry report. "
             "Omitting sbid/version defaults that axis to latest.",
    )
    parser.add_argument(
        "--outdir", default=".",
        help="Directory to save downloaded cutouts into (with "
             "--run-mode=auto or =select)",
    )
    args = parser.parse_args()
    if args.ra is not None and args.dec is None:
        parser.error("--ra given without --dec")
    if args.run_mode == "select" and not args.select:
        parser.error("--run-mode=select requires --select=<list>")
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


def version_number_from_filename(filename: str) -> int:
    """Reprocessing version tag from the known .vN. filename convention
    (confirmed: image.restored.i.SB10083.contcube.v2.fits) -- 1 if
    no tag is present (the implicit baseline version). Matches the tag
    anywhere in the filename rather than assuming a fixed position
    relative to .conv. -- only a non-.conv example has been confirmed
    so far, and that same SBID's .conv sibling had no version tag at
    all, so a version tag is not guaranteed to apply uniformly across
    one SBID's whole file set (see docs/dev/CASDA_FETCH_PLAN.md T1a)."""
    m = re.search(r"\.v(\d+)\.", filename, re.IGNORECASE)
    return int(m.group(1)) if m else 1


def select_iquv_mask(rows, version: Optional[int] = None) -> list:
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
    the `usable` filter in main()), so hitting this path at all
    indicates a mismatch between that filtering and what's in `rows`;
    it's still enforced here directly rather than trusted to the
    caller, since a partial fetch would be actively misleading, not
    just redundant.

    I and V: .conv. is a soft preference (they're for source verification/
    noise-floor checks per the user's own stated use, not tomography
    itself) -- falls back to the plain file if no .conv. variant exists,
    but prints a loud warning, since Q/U being .conv while I/V aren't is
    an anomaly worth flagging rather than accepting quietly. This only
    ever applies once Q and U have already both been confirmed .conv,
    per the gate above.

    version: reprocessing version to request (see
    version_number_from_filename) -- None (the default) picks the
    highest version present per (Stokes letter, conv-status) group.
    Resolved INDEPENDENTLY per Stokes letter, never shared between Q
    and U the way .conv is: Q and U are reprocessed independently (a
    fix that touches Q's calibration need not touch U's, and a version
    tag is confirmed not to apply uniformly across one SBID's whole
    file set -- see version_number_from_filename), so requiring a
    version requested for one to also be present on the other would be
    wrong, not just strict. A letter missing the requested version
    falls back to its own highest available version -- expected
    behaviour, not an anomaly, so this is not a warning the way the
    conv-status fallback below is."""
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

    def pick_version(idxs: list, letter: str) -> int:
        versions = {i: version_number_from_filename(str(rows["filename"][i]))
                    for i in idxs}
        if version is not None:
            matches = [i for i in idxs if versions[i] == version]
            if matches:
                return matches[0]
            available = sorted(set(versions.values()))
            print(f"    Stokes {letter.upper()}: requested v{version} not "
                  f"available here (has: {available}) -- using the highest "
                  f"available version instead.")
        return max(idxs, key=lambda i: versions[i])

    keep = set()
    for letter in STOKES_LETTERS:
        idxs = by_letter.get(letter, [])
        if not idxs:
            print(f"    Stokes {letter.upper()}: no file found for this observation")
            continue
        conv_idxs = [i for i in idxs if ".conv." in str(rows["filename"][i])]
        if conv_idxs:
            chosen = pick_version(conv_idxs, letter)
            v = version_number_from_filename(str(rows["filename"][chosen]))
            print(f"    Stokes {letter.upper()}: {len(idxs)} file(s) matched, "
                  f"using .conv. variant (v{v}): {rows['filename'][chosen]}")
            keep.add(chosen)
        else:
            # Only I/V can reach here -- Q/U conv-presence was already
            # confirmed by the gate above.
            chosen = pick_version(idxs, letter)
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
    conv_letters: set = field(default_factory=set)
    t_min: Optional[float] = None
    file_versions: dict = field(default_factory=dict)
    s_ra: Optional[float] = None
    s_dec: Optional[float] = None
    s_fov: Optional[float] = None

    @property
    def has_conv_qu(self) -> bool:
        """Whether a .conv. (common-beam) variant is available for BOTH
        Q and U specifically -- known pre-fetch, from the filename alone,
        no download needed. This is what fetching will select for Q/U
        (select_iquv_mask always prefers .conv. when present, and
        highest version -- see file_versions)."""
        return "q" in self.conv_letters and "u" in self.conv_letters


def summarise_observations(table) -> list:
    """One Observation per distinct obs_id in table: its frequency
    range, which Stokes letters have a .conv. filename variant
    available, its observation start time (t_min, MJD -- an ObsCore
    field, confirmed present and populated live, not needing a
    download), every reprocessing version seen per (Stokes letter,
    conv-status) file group (see version_number_from_filename), and its
    field footprint (s_ra/s_dec/s_fov, degrees -- also ObsCore fields,
    describe the pointing/mosaic as a whole, so read once per obs_id
    rather than per row; see observation_covers_target for why this is
    needed at all)."""
    by_obs: dict = {}
    for row in table:
        obs_id = str(row["obs_id"])
        obs = by_obs.setdefault(obs_id, Observation(obs_id, None))
        freq = freq_range_mhz(row)
        if freq is not None:
            obs.freq = freq
        if not np.ma.is_masked(row["t_min"]):
            t_min = float(row["t_min"])
            obs.t_min = t_min if obs.t_min is None else min(obs.t_min, t_min)
        if obs.s_fov is None and not np.ma.is_masked(row["s_fov"]):
            obs.s_ra = float(row["s_ra"])
            obs.s_dec = float(row["s_dec"])
            obs.s_fov = float(row["s_fov"])
        filename = str(row["filename"])
        letter = stokes_letter_from_filename(filename)
        if letter:
            is_conv = ".conv." in filename
            if is_conv:
                obs.conv_letters.add(letter)
            key = (letter, is_conv)
            obs.file_versions.setdefault(key, set()).add(
                version_number_from_filename(filename))
    return sorted(by_obs.values(), key=lambda o: o.freq or (0.0, 0.0))


def observation_covers_target(obs: Observation, coord: SkyCoord,
                               cutout_radius_arcmin: float) -> bool:
    """Whether this observation's field plausibly covers a cutout of
    the given radius around coord -- a conservative circular
    approximation of the field's own (possibly irregular, polygon)
    footprint using its own s_ra/s_dec/s_fov, since CASDA's
    query_region cone search only guarantees the field is near the
    search circle, not that the target position (let alone a cutout
    radius around it) is within the field's coverage. Confirmed: a
    field's s_fov can be tens of degrees across (EMU_2132-51's own is
    ~54.7 deg) with the target well off-center, so an edge-of-search-
    radius match could plausibly miss the target entirely. Uses the
    field's own nominal FOV circle rather than the polygon in
    ObsCore's own s_region -- this is a sanity check, not a precise
    containment test."""
    if obs.s_fov is None or obs.s_ra is None or obs.s_dec is None:
        return False
    field_center = SkyCoord(ra=obs.s_ra * u.deg, dec=obs.s_dec * u.deg)
    sep_deg = float(coord.separation(field_center).deg)
    return bool(sep_deg + (cutout_radius_arcmin / 60.0) < obs.s_fov / 2.0)


def filter_target_coverage(observations: list, coord: SkyCoord,
                            cutout_radius_arcmin: float) -> list:
    """Excludes any observation whose own field footprint doesn't
    plausibly cover a cutout_radius_arcmin cutout around coord (see
    observation_covers_target) -- printed, not silent, since this can
    otherwise silently remove the only candidate for a band, or make
    "latest wins" pick between two SBIDs where the latest one doesn't
    even reach the target (see docs/dev/CASDA_FETCH_PLAN.md T1a)."""
    kept = []
    for obs in observations:
        if observation_covers_target(obs, coord, cutout_radius_arcmin):
            kept.append(obs)
        else:
            sep = "unknown (s_fov not populated)"
            if obs.s_ra is not None and obs.s_dec is not None:
                field_center = SkyCoord(ra=obs.s_ra * u.deg, dec=obs.s_dec * u.deg)
                sep = f"{coord.separation(field_center).deg:.2f} deg from field centre, field s_fov={obs.s_fov:.1f} deg"
            print(f"  excluding {obs.obs_id}: field footprint doesn't "
                  f"plausibly cover a {cutout_radius_arcmin:.1f} arcmin "
                  f"cutout at the target ({sep})")
    return kept


@dataclass
class Band:
    number: int
    freq: tuple
    observations: list  # Observation, sorted latest t_min first


def group_into_bands(usable: list) -> list:
    """Groups usable (Q/U-.conv-complete) observations by exact
    frequency-range match, numbered 1..N ascending by frequency for
    stable reference in --select -- confirmed that more than one
    SBID can share a band (PKS 2130-538: three SBIDs at
    799.5-1087.5 MHz). Each band's own observations are sorted
    latest-t_min-first, so index 0 is always the --run-mode=auto pick.

    Band numbering is only stable between two runs if CASDA's public
    holdings for this position haven't changed in between -- a --select
    built from an older dry run could silently refer to a different
    band if new data appeared meanwhile. Not otherwise guarded against;
    each report shows the frequency range next to the number so a
    mismatch is at least visible, and resolve_selection re-derives
    numbering fresh from the CURRENT query every run, never trusting a
    stale mapping."""
    by_freq: dict = {}
    for obs in usable:
        by_freq.setdefault(obs.freq, []).append(obs)
    bands = []
    for number, freq in enumerate(sorted(by_freq), start=1):
        obs_list = sorted(
            by_freq[freq],
            key=lambda o: o.t_min if o.t_min is not None else -1.0,
            reverse=True,
        )
        bands.append(Band(number=number, freq=freq, observations=obs_list))
    return bands


def mjd_to_date(mjd: float) -> str:
    return Time(mjd, format="mjd").datetime.strftime("%Y-%m-%d")


def print_band_report(bands: list, n_observations_total: int) -> None:
    """The dry-run (and pre-fetch confirmation) report: one entry per
    distinct frequency band, every candidate SBID within it
    (latest-observed first), and any per-file version ambiguity --
    enough detail to build a --select list without guessing (see
    docs/dev/CASDA_FETCH_PLAN.md T1a)."""
    n_usable = sum(len(b.observations) for b in bands)
    print(f"\n{n_usable} of {n_observations_total} observation(s) have "
          f"both Q and U present AS .conv (beam-wise matched), across "
          f"{len(bands)} distinct frequency band(s).")
    if not bands:
        print("No band with both Q and U present as .conv -- this "
              "position isn't usable yet with what CASDA currently has "
              "public.")
        return
    if len(bands) >= 2:
        print("Multi-band tomography looks possible.")

    for band in bands:
        lo, hi = band.freq
        print(f"\nBand {band.number}: {lo:.1f}-{hi:.1f} MHz "
              f"({len(band.observations)} candidate SBID(s))")
        for i, obs in enumerate(band.observations):
            sbid = sbid_number(obs.obs_id)
            date_str = mjd_to_date(obs.t_min) if obs.t_min is not None else "unknown date"
            default_tag = " -- DEFAULT (latest observation)" if i == 0 else ""
            print(f"  SBID {sbid} (observed {date_str}){default_tag}")
            for (letter, is_conv), versions in sorted(obs.file_versions.items()):
                if len(versions) > 1:
                    conv_str = ".conv" if is_conv else "(plain)"
                    version_str = ", ".join(f"v{v}" for v in sorted(versions))
                    print(f"      Stokes {letter.upper()} {conv_str} has "
                          f"{len(versions)} versions: {version_str} -- "
                          f"v{max(versions)} used by default")

    print(
        "\nTo fetch: --run-mode=auto (every band, latest SBID/version "
        "throughout) or --run-mode=select --select=\"<list>\", <list> a "
        "comma-separated band[:sbid[:version]] list using the band "
        "numbers above, e.g. --select \"1,2:51797,3:74876:v2\" -- "
        "omitting sbid/version defaults that axis to latest; a single "
        "band number alone is single-band mode."
    )


def parse_select(select_str: str) -> dict:
    """band_number -> (sbid_number_or_None, version_or_None) from a
    comma-separated 'band[:sbid[:version]]' list (see --help/module
    docstring). An omitted sbid or version means "latest" for that
    axis. Raises ValueError with a human-readable reason on malformed
    input -- the caller is expected to print it and stop, not raise
    past the user."""
    result = {}
    for entry in select_str.split(","):
        entry = entry.strip()
        if not entry:
            continue
        parts = entry.split(":")
        try:
            band_num = int(parts[0])
        except ValueError:
            raise ValueError(f"invalid band number in --select entry {entry!r}")
        sbid = None
        if len(parts) > 1 and parts[1]:
            try:
                sbid = int(parts[1])
            except ValueError:
                raise ValueError(f"invalid SBID in --select entry {entry!r}")
        version = None
        if len(parts) > 2 and parts[2] and parts[2].lower() != "latest":
            m = re.match(r"^v?(\d+)$", parts[2], re.IGNORECASE)
            if not m:
                raise ValueError(f"invalid version in --select entry {entry!r}")
            version = int(m.group(1))
        result[band_num] = (sbid, version)
    return result


def resolve_selection(bands: list, run_mode: str, select_str: Optional[str]) -> Optional[list]:
    """[(Observation, version_or_None), ...] to fetch, or None (with a
    printed reason already shown) if the request can't be resolved.

    --run-mode=auto: every band, each band's own latest-observed SBID
    (index 0, see group_into_bands), version=None (select_iquv_mask's
    own latest-version default) -- "fix ALL cfg params to auto" per the
    user's own framing of this mode.

    --run-mode=select: exactly what --select names, with band/SBID all
    validated against what group_into_bands found for THIS
    run's own query (never a stale/remembered mapping from an earlier
    dry run -- see group_into_bands' own caveat on this)."""
    if run_mode == "auto":
        return [(band.observations[0], None) for band in bands]

    by_number = {band.number: band for band in bands}
    try:
        selection = parse_select(select_str)
    except ValueError as exc:
        print(f"\nERROR: {exc}")
        return None
    if not selection:
        print("\nERROR: --select is empty")
        return None

    chosen = []
    for band_num, (sbid, version) in selection.items():
        band = by_number.get(band_num)
        if band is None:
            available = sorted(by_number)
            print(f"\nERROR: no band numbered {band_num} in this run's "
                  f"discovery result -- available band numbers: "
                  f"{available} (re-run with --run-mode=dry first if "
                  f"unsure of the current numbering; it's derived fresh "
                  f"from CASDA's current holdings every run, not cached)")
            return None
        if sbid is None:
            obs = band.observations[0]
        else:
            obs = next((o for o in band.observations
                        if sbid_number(o.obs_id) == sbid), None)
            if obs is None:
                available = [sbid_number(o.obs_id) for o in band.observations]
                print(f"\nERROR: SBID {sbid} is not a candidate for band "
                      f"{band_num} -- available: {available}")
                return None
        chosen.append((obs, version))
    return chosen


def check_band_overlaps(selection: list) -> bool:
    """True if no two observations about to be fetched together
    overlap in frequency -- catches both a spectral overlap
    between two truly distinct bands and the same band mistakenly
    split into two by floating-point noise in em_min/em_max between
    SBIDs (near-duplicate ranges trivially overlap almost entirely
    either way, so one check covers both causes). Prints the
    overlapping pair and by how much, and returns False, if any
    overlap is found -- the caller must not proceed to fetch, since
    combining overlapping channels into one multi-band rm_synthesis
    run would double-count that frequency range in the combined RMSF/
    sensitivity, silently biasing the result. Resolving this (picking
    only one of the overlapping bands via --select) is left to the
    user -- this only detects and reports it."""
    ok = True
    for i in range(len(selection)):
        for j in range(i + 1, len(selection)):
            obs1, _ = selection[i]
            obs2, _ = selection[j]
            lo1, hi1 = obs1.freq
            lo2, hi2 = obs2.freq
            if lo1 < hi2 and lo2 < hi1:
                overlap = min(hi1, hi2) - max(lo1, lo2)
                print(f"\nERROR: {obs1.obs_id} ({lo1:.3f}-{hi1:.3f} MHz) and "
                      f"{obs2.obs_id} ({lo2:.3f}-{hi2:.3f} MHz) overlap in "
                      f"frequency by {overlap:.3f} MHz -- combining both "
                      f"into one multi-band run would double-count that "
                      f"range in the combined RMSF/sensitivity. Choose only "
                      f"one of these two (--select), or otherwise resolve "
                      f"this yourself before fetching.")
                ok = False
    return ok


def print_fallback_diagnostics(table) -> None:
    print("\nNothing matched --project/--subtype. Here's what IS public "
          "near this position, to diagnose why:")
    for column in ("obs_collection", "dataproduct_subtype"):
        values = sorted(set(str(v) for v in table[column]))
        print(f"\nDistinct {column} values ({len(values)}):")
        for v in values:
            n = sum(1 for row in table[column] if str(row) == v)
            print(f"  {v!r}: {n} row(s)")


def sanitize_run_name(s: str) -> str:
    """Filesystem- and cfg-value-safe run name from a target string or
    coordinate pair -- alphanumerics kept, everything else collapsed to
    a single underscore."""
    s = re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")
    return s or "rmtool_fetch"


def stage_symlink(real_path: str, staging_dir: str, link_name: str) -> str:
    """Symlink real_path into staging_dir under link_name. Needed because
    scripts/run_pipeline.sh's own pipeline cfg takes match_infileQ=/
    match_infileU= as BARE filenames under one shared match_input_path=,
    and its own symlink step (which relocates match_cubes' writes onto
    outdir) does not create intermediate directories -- confirmed
    directly from that script -- so a per-SBID subdirectory can't appear
    inside the filename the way this tool's own <outdir>/SB<n>/ fetch
    layout would otherwise suggest. A pre-existing symlink pointing at
    the same file is reused across calls (harmless); one pointing
    elsewhere is a conflict -- raises rather than silently
    repointing, mirroring run_pipeline.sh's own policy for exactly this
    situation."""
    link_path = os.path.join(staging_dir, link_name)
    real_path = os.path.abspath(real_path)
    if os.path.islink(link_path) or os.path.exists(link_path):
        if os.path.realpath(link_path) != os.path.realpath(real_path):
            raise RuntimeError(
                f"{link_path} already points elsewhere -- remove it first "
                f"if this is intentional")
    else:
        os.symlink(real_path, link_path)
    return link_path


RM_FAC = 3.14159265358979  # rm_synthesis' own fac= default; matches
# every hand-written cfg in this repo -- see cfg/rmsynth.cfg.
GALACTIC_BEG_RM = -1000.0
GALACTIC_END_RM = 1000.0


def band_lambda_sq_span(fits_path: str) -> float:
    """dlam2_band -- the half-channel-edge-extended lambda^2 span (m^2)
    of this cube's own FREQ axis, mirroring src/rm_synthesis.f90's own
    multi-band RM-range diagnostic exactly (its lam2_lo_edge/
    lam2_hi_edge/dlam2_band computation) so the generated nrm is
    consistent with what that diagnostic would itself report for this
    data. Assumes the full channel range (no chan_blc/chan_inc
    subsetting -- matches this tool's own generated cfg) and Hz-scale
    CRVAL/CDELT (confirmed directly against a downloaded cube's
    header: CTYPE4='FREQ', CRVAL4~8.0e8, CDELT4=1.0e6)."""
    header = fits.getheader(fits_path)
    freq_axis = None
    for axis in (1, 2, 3, 4):
        if str(header.get(f"CTYPE{axis}", "")).strip().upper().startswith("FREQ"):
            freq_axis = axis
            break
    if freq_axis is None:
        raise ValueError(f"{fits_path}: no FREQ axis in header")
    crval = float(header[f"CRVAL{freq_axis}"])
    cdelt = float(header[f"CDELT{freq_axis}"])
    crpix = float(header[f"CRPIX{freq_axis}"])
    naxis = int(header[f"NAXIS{freq_axis}"])
    z1 = crval - (crpix - 1.0) * cdelt
    zn = z1 + (naxis - 1.0) * cdelt
    lo_freq = min(z1, zn) - 0.5 * abs(cdelt)
    hi_freq = max(z1, zn) + 0.5 * abs(cdelt)
    lam2_lo_edge = (SPEED_OF_LIGHT_M_PER_S / lo_freq) ** 2
    lam2_hi_edge = (SPEED_OF_LIGHT_M_PER_S / hi_freq) ** 2
    return lam2_lo_edge - lam2_hi_edge


def compute_rm_range(q_paths: list) -> tuple:
    """(beg_rm, end_rm, nrm) for the generated rm_synthesis template.

    use_auto_rm_range=1 is hard-rejected by rm_synthesis for any
    multi-band run (src/rm_synthesis.f90:1002-1008) -- the auto-range
    heuristic infers the dataset's min/max frequency from a merged
    lambda^2 array's own endpoints, which only holds for one
    internally-monotonic band. This project uses a fixed Galactic-scale
    range instead, uniformly, even for single-band fetches, for one
    consistent generated-cfg policy rather than branching on band count
    (confirmed: not wanting to spend time processing RMs outside the
    range of interest either way).

    nrm is the number of RM samples across that range at this data's
    own combined delta_RM (RMSF FWHM) resolution -- mirroring
    src/rm_synthesis.f90's own multi-band diagnostic formula exactly
    (combined delta_RM = fac / sum(per-band lambda^2 span)) so the
    generated value is consistent with what that diagnostic would
    itself report for this data, rather than a value derived in
    isolation. ofac (still 4 in the generated template) provides the
    oversampling on top of this nrm, exactly as it already does for
    every hand-written cfg in this repo (nrm_out = nrm * ofac) -- this
    is not a second, independent oversampling factor."""
    total_dlam2 = sum(band_lambda_sq_span(p) for p in q_paths)
    combined_delta_rm = RM_FAC / total_dlam2
    nrm = math.ceil((GALACTIC_END_RM - GALACTIC_BEG_RM) / combined_delta_rm)
    return GALACTIC_BEG_RM, GALACTIC_END_RM, max(nrm, 1)


def write_rmsynth_template(path: str, beg_rm: float, end_rm: float, nrm: int) -> None:
    body = f"""\
# rm_synthesis config auto-generated by scripts/casda_fetch.py.
# path=/infileQ=/infileU=/outfile= below are placeholders -- run_pipeline.sh
# overwrites all four automatically from the match stage's own manifest
# when this is run together with the match stage (stages=match,rmsynth in
# the pipeline cfg that references this template). Only fill these in by
# hand if running rm_synthesis on this template directly, after a
# separate match_cubes run of your own.
path = .
infileQ = placeholder.fits
infileU = placeholder.fits
outfile = placeholder

subim = n
rem_mean = 0

# Bias correction is a documented no-op in the current rm_synthesis (see
# docs/user/APP_REFERENCE.md) -- these four are required regardless of
# remove_qu_bias, always 0.0 while it's off.
remove_qu_bias = n
resiQ = 0.0
slopeQ = 0.0
resiU = 0.0
slopeU = 0.0

ofac = 4
fac = {RM_FAC}

# use_auto_rm_range=1 is rejected outright by rm_synthesis for
# multi-band runs (src/rm_synthesis.f90:1002-1008) -- this project uses
# a fixed Galactic-scale range and a computed nrm instead, uniformly
# (see docs/dev/CASDA_FETCH_PLAN.md T1b). nrm below is this fetch's own
# combined delta_RM (RMSF FWHM) sampled across the range below; ofac
# above provides the oversampling on top of it.
use_auto_rm_range = 0
beg_rm = {beg_rm}
end_rm = {end_rm}
nrm = {nrm}

output_mode = ap
ap_angle_mode = phase

write_mask_output = y
write_nvalid_output = y
cubestat = n
use_gpu = n
mem_frac_ram = 0.25
log_level = info
"""
    with open(path, "w") as f:
        f.write(body)


def generate_pipeline_cfg(outdir: str, run_name: str, records: list) -> Optional[str]:
    """records: one (sbid_dir, q_path, u_path, beamfile_or_None) tuple
    per fetched observation -- already resolved to exactly one SBID per
    band by resolve_selection before this is called, never "everything
    usable" (see docs/dev/CASDA_FETCH_PLAN.md T1a on why that would be
    wrong: separate SBIDs at the same band are re-observations/
    reprocessings, not independent bands). Writes a flattened staging
    directory of symlinks, an rm_synthesis template (RM range/nrm from
    compute_rm_range, using every record's own Q cube header), and a
    top-level scripts/run_pipeline.sh-compatible pipeline cfg chaining
    match->rmsynth, so `./scripts/run_pipeline.sh <returned path>` runs
    end to end with no manual editing. A single fetched observation uses
    match_args' stages=convolve (nothing to reproject -- one
    observation's own Q/U already share one native grid); two or more
    use stages=both (separate SBIDs/bands may not share a grid).

    rmclean is deliberately left out of the generated stages= -- its
    CLEAN stopping criteria are a choice for the user to make (see
    docs/user/EXAMPLES.md's own dedicated section on this), not
    something to default silently.

    Returns the pipeline cfg's path, or None if `records` is empty
    (nothing fetched to build a pipeline around)."""
    if not records:
        return None

    staging_dir = os.path.join(outdir, "pipeline_staging")
    os.makedirs(staging_dir, exist_ok=True)

    q_names, u_names, beamfiles = [], [], []
    for sbid_dir, q_path, u_path, beamfile in records:
        q_link = stage_symlink(q_path, staging_dir, f"{sbid_dir}_Q.fits")
        u_link = stage_symlink(u_path, staging_dir, f"{sbid_dir}_U.fits")
        q_names.append(os.path.basename(q_link))
        u_names.append(os.path.basename(u_link))
        beamfiles.append(os.path.abspath(beamfile) if beamfile else "auto")

    match_stages = "convolve" if len(records) == 1 else "both"
    match_args_parts = [f"stages={match_stages}", "order=convolve_reproject"]
    if match_stages == "both":
        match_args_parts.append("footprint_mode=intersection")
        first_q_abspath = os.path.abspath(os.path.join(staging_dir, q_names[0]))
        match_args_parts.append(f"reffile={first_q_abspath}")
    # Q-block then U-block, matching infiles='s own order (all Q bands,
    # then all U bands -- see scripts/run_pipeline.sh's match stage) --
    # same value in both blocks per band, since Q and U share one
    # restoring beam (confirmed, see docs/dev/CASDA_FETCH_PLAN.md).
    beamfiles_list = beamfiles + beamfiles
    match_args_parts.append(f"beamfiles={','.join(beamfiles_list)}")
    match_args_parts.append("mem_frac_ram=0.25")
    match_args = " ".join(match_args_parts)

    beg_rm, end_rm, nrm = compute_rm_range([q_path for _, q_path, _, _ in records])
    rmsynth_template_path = os.path.join(outdir, f"{run_name}_rmsynth.cfg")
    write_rmsynth_template(rmsynth_template_path, beg_rm, end_rm, nrm)

    pipeline_cfg_path = os.path.join(outdir, f"{run_name}_pipeline.cfg")
    with open(pipeline_cfg_path, "w") as f:
        f.write(
            f"# Pipeline config auto-generated by scripts/casda_fetch.py -- run with:\n"
            f"#   ./scripts/run_pipeline.sh {pipeline_cfg_path}\n"
            f"# rmclean is not included by default -- its CLEAN stopping criteria are\n"
            f"# a choice for you to make (see docs/user/EXAMPLES.md); add\n"
            f"# 'rmclean' to stages= and an rmclean_cfg_template= once you've looked\n"
            f"# at the rmsynth output and decided on a stopping criterion.\n"
            f"stages = match,rmsynth\n"
            f"outdir = {outdir}\n"
            f"run_name = {run_name}\n"
            f"\n"
            f"match_input_path = {staging_dir}\n"
            f"match_infileQ = {','.join(q_names)}\n"
            f"match_infileU = {','.join(u_names)}\n"
            f"match_args = {match_args}\n"
            f"\n"
            f"rmsynth_cfg_template = {rmsynth_template_path}\n"
            f"rmsynth_backend = auto\n"
            f"rmsynth_omp_threads = 6\n"
        )
    return pipeline_cfg_path


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

    cutout_radius = args.cutout_radius if args.cutout_radius is not None else args.radius

    observations = summarise_observations(subset)
    usable = [o for o in observations if o.has_conv_qu and o.freq is not None]
    usable = filter_target_coverage(usable, coord, cutout_radius)
    bands = group_into_bands(usable)
    print_band_report(bands, len(observations))

    if args.run_mode == "dry":
        print("\n(dry run -- pass --run-mode=auto or --run-mode=select "
              "--select=... to fetch)")
        return

    if not bands:
        print("\nNothing usable to fetch.")
        return

    selection = resolve_selection(bands, args.run_mode, args.select)
    if selection is None:
        return

    if not check_band_overlaps(selection):
        return

    print(f"\nFetching {len(selection)} observation(s) into "
          f"{args.outdir}/<SBID>/, cutout radius={cutout_radius} arcmin...")

    # One subdirectory per observation (one SBID) -- NOT necessarily one
    # per distinct band; resolve_selection has already reduced this to
    # exactly one SBID per chosen band (see docs/dev/CASDA_FETCH_PLAN.md
    # T1a on why fetching every Q/U-complete SBID at a shared band would
    # be wrong -- they're re-observations/reprocessings, not independent
    # bands).
    #
    # All of an observation's (up to 4) selected files are requested in
    # ONE batched cutout() call, not one call per file -- CASDA processes
    # them together server-side, which matters for wall time (a single
    # ~1GB cutout has been observed to take upwards of a minute; doing
    # that 4 separate times, with 4 separate job-submit-and-poll cycles,
    # would multiply it for no benefit once file identification no longer
    # depends on request order -- see below). Each observation is still
    # ATTEMPTED independently of the others, so one network drop doesn't
    # stop the script from at least trying the rest and reporting on all
    # of them -- but the run as a whole only produces a pipeline cfg if
    # EVERY requested observation fully succeeds (see the all-succeeded
    # check after this loop): once --run-mode has resolved a specific
    # set of bands as the user's intent, silently handing back a cfg for
    # fewer bands than that would misrepresent what was achieved, which
    # is worse than refusing outright and saying why.
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
    pipeline_records = []
    all_succeeded = True
    for obs, version in selection:
        mask = [str(v) == obs.obs_id for v in subset["obs_id"]]
        rows = subset[mask]
        mask_pol = [str(v) in ("/I/", "/Q/", "/U/", "/V/") for v in rows["pol_states"]]
        rows = rows[mask_pol]

        obs_dir = os.path.join(args.outdir, sbid_dirname(obs.obs_id))
        os.makedirs(obs_dir, exist_ok=True)
        freq_str = f"{obs.freq[0]:.1f}-{obs.freq[1]:.1f} MHz" if obs.freq else "unknown freq"
        print(f"\n{obs.obs_id} ({freq_str}) -> {obs_dir}/")
        rows = rows[select_iquv_mask(rows, version=version)]
        print(f"  {len(rows)} file(s) selected (one per Stokes parameter found, "
              f".conv preferred over plain restored where both exist)")

        if len(rows) == 0:
            print(f"  SKIPPED: Q and/or U lack a .conv. variant -- nothing "
                  f"fetched for this observation (see ERROR above)")
            obs_results.append((obs.obs_id, "SKIPPED"))
            all_succeeded = False
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
            print(f"  (likely a transient network/server issue -- the rest "
                  f"of this run will still be attempted, but since this "
                  f"observation didn't fully succeed, no pipeline config "
                  f"will be generated at the end -- see the summary below)")
        if n_ok == n_total:
            status = "OK"
        elif n_ok == 0:
            status = "FAILED"
        else:
            status = f"PARTIAL ({n_ok}/{n_total} files)"

        if "Q" in written and "U" in written:
            beamfile = None
            if has_beams_table(written["Q"]) and has_beams_table(written["U"]):
                print(f"  Q and U both already carry a per-channel BEAMS "
                      f"table -- convolve_cubes' beamfiles=auto will read "
                      f"it directly, no separate beamlog needed.")
            else:
                print(f"  Q and/or U's .conv cube has no BEAMS table of "
                      f"its own -- fetching the observatory's beamlog to "
                      f"curate one:")
                beamfile = fetch_and_curate_beamlog(obs.obs_id, casda._auth, obs_dir)
                if beamfile is None:
                    status = "FAILED (beamlog)"
                    print(f"  {obs.obs_id}'s .conv cube has no BEAMS table "
                          f"and no beamlog could be curated for it -- "
                          f"convolve_cubes' beamfiles=auto would fail on "
                          f"this file, so this observation cannot be used "
                          f"for tomography without manual intervention.")
            pipeline_records.append(
                (sbid_dirname(obs.obs_id), written["Q"], written["U"], beamfile))

        if status != "OK":
            all_succeeded = False
        obs_results.append((obs.obs_id, status))

    print("\nFetch summary:")
    for obs_id, status in obs_results:
        print(f"  {obs_id}: {status}")

    if not all_succeeded:
        print(
            "\nABORTED: not every requested observation fully succeeded "
            "(see Fetch summary above) -- no pipeline config generated, "
            "since it would silently deliver less than what was asked "
            "for. Re-run --run-mode=select with just the observation(s) "
            "that succeeded, or retry once the failing one(s) are "
            "resolved."
        )
        return

    run_name = sanitize_run_name(
        args.target if args.target else f"ra{coord.ra.deg:.3f}_dec{coord.dec.deg:.3f}")
    pipeline_cfg_path = generate_pipeline_cfg(args.outdir, run_name, pipeline_records)
    print(
        f"\nPipeline config generated from all {len(pipeline_records)} "
        f"requested observation(s): {pipeline_cfg_path}\n"
        f"Run it with: ./scripts/run_pipeline.sh {pipeline_cfg_path}\n"
        f"This chains match_cubes (resolution/grid matching -- "
        f"beamfiles= already set per band, curated or auto as printed "
        f"above) into rm_synthesis. It does NOT include rmclean_cubes "
        f"or use the I/V cubes (source-verification/noise-floor only, "
        f"not part of tomography) -- see the generated cfg's own "
        f"header comments to extend it."
    )


if __name__ == "__main__":
    main()
