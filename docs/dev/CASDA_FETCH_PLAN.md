# Fetching calibrated Q/U multiband cubes from public archives — Plan

**Status: T0 (discovery, fetch, beam curation) working end-to-end against
CASDA data. T1 (pipeline wiring) not started.**

## 1. Motivation

For rmtool to be widely usable, it should be easy for a new user to try it
against calibrated Stokes Q/U multiband data, not just the synthetic test
fixtures under `tests/data/`. This plan covers scripted fetching of small
cutouts from public radio astronomy archives, ideally wired into the
e2e pipeline (`scripts/run_pipeline.sh`) as a `fetch` stage ahead of
`match`.

## 2. Scope decision (confirmed with user, 2026-09-04)

**CASDA only, for now.** CASDA (CSIRO ASKAP Science Data Archive) hosts
every ASKAP survey product (WALLABY, EMU, POSSUM, RACS, ...) from one
telescope with one internally self-consistent calibration pipeline.
Mixing in a second archive (e.g. VLA/VLASS, a different telescope with its
own independent calibration chain) was explicitly ruled out for this first
pass — the user's reasoning: cross-archive combination risks flux-scale/
calibration-convention mismatches that don't exist when staying within one
archive's own products.

Within CASDA, the relevant survey is **POSSUM** (project code **AS203**),
the ASKAP polarization survey riding commensally on WALLABY/EMU
observations — WALLABY (HI line) and EMU (plain continuum) do not
themselves produce Stokes Q/U products; POSSUM does, in the same fields.
(An earlier, pilot-era project code AS103 was checked and set aside —
POSSUM's own data page describes that code's CASDA holdings as Stokes I/V
only; AS203 is the full-survey code with public I/Q/U/V cubes.)

## 3. Access mechanism

CASDA implements IVOA-standard **TAP** (query/discovery) + **SODA**
(server-side cutouts) + **Datalink** — confirmed directly from CSIRO's own
`casda-samples/cutouts.py` sample script and the CASDA VO Tools repo
(`csiro-rds/casda_vo_tools`). `astroquery.casda` is a maintained Python
wrapper around exactly this, and is what this plan builds on rather than
hand-rolling TAP/ADQL queries or scraping the CASDA web UI. (`pyvo` was
also considered as a more archive-agnostic alternative to `astroquery`'s
per-archive submodules — more general, but `astroquery.casda`
already covers everything CASDA-specific needed here, so there's no
concrete win from swapping yet; worth revisiting only if a second
VO-compliant archive is ever added to this plan's scope.)

Confirmed API surface (`astroquery.casda.CasdaClass`):
- `Casda().login(*, username=None, store_password=False,
  reenter_password=False)` — **no password parameter at all.** The only
  way to authenticate non-interactively (needed for an unattended pipeline
  run) is a one-time interactive setup that stores the password in the OS
  keyring (`store_password=True`), after which every later
  `login(username=...)` call reads it back silently. There is no way to
  pass a password via CLI arg or env var through this API — by design, so
  a password never ends up in a config file, shell history, or process
  list.
- `Casda.query_region(coordinates, radius=..., ...)` — class-level (not
  instance) call; returns an `ivoa.obscore`-conformant table.
- `Casda.filter_out_unreleased(result)` — also class-level; must be called
  explicitly or the result can include proprietary/embargoed rows the
  account can't download.
- `casda.cutout(table, *, coordinates=None, radius=..., height=None,
  width=None, band=None, channel=None)` — instance call; spatial cutout via
  `coordinates`+`radius` (or `height`/`width`), spectral cutout via `band=`
  (frequency/wavelength range) or `channel=` (channel-index range). Returns
  staged URLs, not files.
- `casda.download_files(urls, savedir=...)` — instance call, the download
  step.

**Credentials:** an OPAL account (ATNF's Online Proposal Access & Links,
free self-registration at `opal.atnf.csiro.au` — email/name/affiliation/
password, no institutional access needed) is what `login()` is built
around and what CASDA's own docs recommend "for general use." A Pawsey HPC
account is a documented alternative but requires an existing Pawsey
allocation — too high a bar for "widely usable," so OPAL is this plan's
default/documented path.

## 4. Candidate target objects

Cited (not guessed) objects investigated for this project. Each was
checked against live CASDA query results, not just literature — a
published detection does not necessarily mean the underlying cube is
public under project AS203 today (see MSH 15-56 below, the concrete case
where this distinction mattered).

- **G326.3−1.8 (MSH 15−56)** — supernova remnant with a pulsar wind nebula,
  cited in the EMU/POSSUM Galactic Pilot Field catalog (arXiv:2307.01948)
  for "clear polarization... visible in both the polarized intensity and
  rotation measure maps." **Checked against live CASDA data: no
  `cont.restored.3d` POSSUM cube is currently public for this specific
  field** (POSSUM-collection rows exist nearby, but none of the confirmed
  science-cube subtype) — a gap between "used in a published catalog" and
  "public in the archive today," not a tool bug. Set aside; may become
  usable as POSSUM's ongoing survey releases more fields.
- **G327.4+0.4 (Kes 27)** — same catalog, polarization detected along the
  NE edge. **Confirmed working**: `cont.restored.3d` Q/U cubes exist for
  two fields near this position (`EMU_1554-55_band1`/`SB43773` and
  `EMU_1605-51`/`SB74627`), in the exact filename convention
  (`image.restored.{q,u}.<field>.SB<n>.contcube[.conv].fits`) already used
  by this project's own WALLABY+EMU e2e cfg. Frequency-band completeness
  not yet re-checked with the current script version.
- **Centaurus A** (southern lobe) — ASKAP polarization data exists and
  produced, per the discovery paper, "amongst the most detailed
  [polarization/RM maps] ever made for radio lobes." Not yet re-checked
  with the current script version whether this is a project-AS203/POSSUM
  product specifically, as opposed to an earlier ASKAP commissioning/
  early-science project.
- **PKS 2130-538, "Dancing Ghosts"** (preset `dancingghosts`) — a bent-tail
  pair of interacting galaxies in Abell 3785, ~7.1 arcmin total extent
  (smaller than Kes 27; much smaller than the ~43 arcmin NGC 2663 or
  ~20 arcmin Kes 27 alternatives considered and rejected for being too
  large for a demo cutout). Sourced from the EMU survey's own curated
  "interesting objects" list (emu-survey.org), not a literature search.
  Position independently confirmed via a SIMBAD cone search (exact match,
  5 decimal places, to `PKS 2130-53`) rather than trusted from
  `SkyCoord.from_name` alone. **Current strongest candidate** — small,
  bent-tail morphology, ASKAP-detected polarization expected for this
  source class. Frequency-band completeness check in progress.

All names are resolved live via `SkyCoord.from_name()` at run time
(Sesame/CDS lookup), not hardcoded as RA/Dec in this doc or in the script —
avoids a transcription error silently pointing the tool at the wrong patch
of sky.

## T0 — `scripts/casda_fetch.py`: is a position usable for tomography?

**Motivation:** the tool's job is to answer one question directly for a
given sky position — how many observations exist, at how many distinct
frequency bands, with a complete Stokes Q/U pair under the subtype
needed for RM synthesis. Two facts, confirmed against CASDA
data (not assumed), make this answerable:
- **`dataproduct_subtype=cont.restored.3d`** is the calibrated,
  per-channel, science-ready cube product — confirmed directly from Kes 27's own
  working files. `cont.cleanmodel/residual/weight.3d` are intermediate
  CLEAN byproducts, not usable directly; `cont.*.t0`/`.t1` are
  frequency-collapsed Taylor-term MFS images, not per-channel data RM
  synthesis needs at all.
- **`pol_states`** (an ObsCore field, e.g. `/Q/`, `/U/`, `/I/`, `/V/` —
  one value per row, one row per Stokes parameter's own file) is what
  distinguishes Stokes parameters — `dataproduct_subtype` alone
  does not.

**Design** (rewritten from an earlier, more exploratory version after it
accumulated ad hoc truncation fixes rather than a clear default view):
resolve a target from a single `--target` flag (checked against this
project's own preset list first, falling back to community name
resolution via `SkyCoord.from_name` if not found there — `--name` as a
separate flag was consolidated away, since it and a non-preset `--target`
value did exactly the same lookup) or explicit `--ra`/`--dec`, query
CASDA within `--radius` (default 15 arcmin), filter to public-only, then to
`--project` (default `POSSUM`, matched against `obs_collection`) and
`--subtype` (default `cont.restored.3d`) — both overridable, both
confirmed-correct defaults rather than placeholders. Groups the result by
`obs_id` (one scheduling block = one observation = one frequency
setup) into a single table: frequency range, `pol_states` seen with
counts, and whether that observation has both Q and U. The bottom-line
verdict counts how many distinct frequency bands are Q+U-complete: 2+
means multi-band tomography is possible, exactly 1 means single-band
only, 0 means this position isn't usable yet. If nothing matches the
project/subtype filters at all, falls back to a plain `obs_collection`/
`dataproduct_subtype` breakdown of what IS public nearby, so "nothing
found" is diagnosable rather than a dead end.

**Two bugs found and fixed during earlier development, both about
hidden truncation** (worth recording since they're easy to reintroduce):
astropy's `Table.pprint()` truncates to an arbitrary head/tail slice by
default, keyed off *terminal screen height* unless `max_lines=-1` is
passed explicitly — and a first attempt at working around this by
grouping preview rows per `dataproduct_subtype` still capped at a fixed
number per group, which did nothing once a single `--subtype` filter had
already narrowed the result to one group. The rewritten design avoids
both by building the observation table from plain Python dicts (no
`pprint` at all) and by making `--subtype`/`--project` defaults do the
narrowing up front, rather than asking a user to iteratively guess a
filter against a truncated view.

**Status:** rewritten script in place. Confirmed against CASDA data
for Kes 27 (multi-band: bands at ~800-1088 MHz and ~1296-1440 MHz both
present, completeness not yet re-verified with this version) and PKS
2130-538 (position confirmed; `--fetch` run end-to-end, see below).
MSH 15-56 confirmed NOT currently usable (see §4). Centaurus A not yet
re-checked with this version.

### Fetching (`--fetch`): what gets selected and downloaded

**Q and U's `.conv.` variant is an all-or-nothing requirement for the
whole observation, not a per-file preference.** `.conv.` in the archive
filename marks that the observatory has already beam-wise matched the
36 ASKAP PAF (Phased Array Feed) beams to each other for that channel —
a different axis from the channel-to-channel matching this project's own
`convolve_cubes` does (see the next section). A plain (non-`.conv.`) Q or
U file hasn't had that beam-wise matching done, so using it would
silently combine spatially-mismatched data into the tomography input.
Confirmed directly against SB74876 that `pol_states` alone is not enough
to select the right file: it matched 4 rows for Q (a raw/`.conv.`
duplicate pair), not the 2 expected, for what is a single-field
observation. `select_iquv_mask()` in `scripts/casda_fetch.py` resolves
this from the filename directly (`stokes_letter_from_filename`,
`.conv.` substring check) — **if EITHER Q or U lacks a `.conv.` variant,
NOTHING is fetched for that whole observation**, not even I/V, since an
observation missing `.conv.` Q/U isn't usable for tomography regardless
of what else is available. `print_verdict()` gates its own usable-count
on this (`Observation.has_conv_qu`), so a `--fetch` run should never
reach this all-or-nothing check with a failing observation in
practice — it's enforced again at fetch time regardless, since a partial
fetch would be actively misleading. An observation that hits this is
reported as `SKIPPED` in the final fetch summary, distinct from
`OK`/`PARTIAL`/`FAILED`.

**All four Stokes parameters are fetched, not just Q/U:** I helps a user
visually confirm there's signal at the target position; V is a
practical noise-floor/data-quality check (best-achievable noise, no
astrophysical V signal expected). I/V prefer `.conv.` too when available
but fall back to the plain file with a loud warning if Q/U have `.conv.`
and I/V anomalously don't — this only matters for I/V, since Q/U are
hard-gated above.

**Batching and resilience:** one observation's up to 4 selected files are
requested in a single batched `cutout()` call (not one call per file —
each cutout job has been observed to take upwards of a minute, and
CASDA's own async job-submit-and-poll cycle would otherwise multiply
that for no benefit). Each *observation* is still fetched independently
of the others; a mid-download network drop on one SBID (observed
directly: `requests.exceptions.ConnectionError` partway through a ~1GB
transfer) does not abort observations not yet reached. Because CASDA's
cutout service renames every downloaded file to its own generic
`cutout-<job>-imagecube-<id>.fits` (discarding the original archive
filename entirely — confirmed directly), which of the batch's downloaded
files is which Stokes parameter is recovered from each file's own header
(`stokes_letter_from_header`: the standard FITS WCS STOKES axis,
`CTYPE<n>='STOKES'`/`CRVAL<n>` → 1=I/2=Q/3=U/4=V) — reliable, unlike
conv-status (next section), which is instead carried forward from what
was already known pre-fetch about that Stokes letter's selected file,
not re-derived from the download.

### Beam matching: curating a target beamlog when the FITS file has none

**Two distinct axes of "matching a beam," easy to conflate:**
1. **Beam-wise** (within one channel, across the 36 PAF beams) — what
   `.conv.` in the archive filename marks as already done by the
   observatory.
2. **Channel-wise** (across the band, for one final position) — each
   channel can still have its own, slightly different achieved beam even
   after step 1; RM synthesis needs ONE common beam across all channels
   of a band. This is what this project's own `convolve_cubes`/
   `match_cubes` do, and is NOT satisfied by `.conv.` alone.

`convolve_cubes` needs to know each channel's *current* beam to convolve
it to the common target. Its default (`beamfiles=auto`) reads that
straight from the FITS file's own CASA-style `BEAMS` binary table
extension. **Whether that table is present is a per-file fact, not a
reliable function of "pilot-era vs not," and must be checked on the
downloaded file directly (`has_beams_table` — checks for the extension
itself, not the `CASAMBM` header keyword).** Confirmed against
SB74876 data (a comparatively recent SBID): its `.conv.` Q cube header
carries `CASAMBM=T` and a complete 288-row `BEAMS` table that matches
the observatory's own beamlog exactly, channel for channel, including
which 2 channels are flagged bad (the FITS table's own placeholder
value, `~1.18e-38`, on channels 160/177, matching the beamlog's own gaps
at those same channels). The plain/raw Q cube has neither the keyword
nor the extension. Older/pilot-era SBIDs are expected to lack the
`BEAMS` extension even on their `.conv.` cube (per the user's own
domain knowledge of when ASKAPsoft started writing it) — `CASAMBM`'s
absence on its own does not distinguish "pilot-era, genuinely no
per-channel beam info" from "just an unset keyword," so the extension
itself is what's checked, not the keyword.

**When a fetched `.conv.` Q or U cube has no `BEAMS` table, the target
beam it needs instead lives in an ASCII "beamlog" file inside a
separate CASDA product: the SBID's diagnostics *evaluation file*
(`diagnostics-SB<n>.tar`).** This is not part of `ivoa.obscore` (the
table `query_region()` queries — confirmed both from CASDA's own user
guide and directly: a position search only ever returns
`dataproduct_type` `cube`/`visibility`) but a wholly separate TAP table,
`casda.observation_evaluation_file` (columns `sbid`, `project_code`,
`format`, `filename`, `filesize`, `access_url`, `released_date`,
`observation_id`, `project_id`), found only by listing every table
CASDA's TAP service exposes (`pyvo.dal.TAPService(...).tables`) rather
than trusting a documentation summary that called the CASDA web UI the
only way to reach it.

Both this table's `access_url` and ObsCore's own `access_url` point at a
**DataLink XML document describing the file, not the file itself**
(confirmed directly: fetching one as a tar fails with "invalid header").
`resolve_datalink_sync_url()` follows this one indirection to the
`sync?id=...` download URL — the row with `semantics="#this"` and a
plain (not `/sync/pawsey`, Pawsey-internal-network-only) endpoint.

**No surgical single-file fetch is possible — confirmed, not assumed.**
CASDA's sync download endpoint does not support HTTP Range requests
(confirmed: `Accept-Ranges: None`, a `Range` header is silently ignored).
The SODA async service CASDA offers for evaluation files declares only
one `inputParam` (`ID`) in its own DataLink XML — no `POS`/`BAND`-style
sub-file parameter like the image-cutout SODA service has — so there is
no server-side way to ask for just the one beamlog inside a multi-GB
tar. `fetch_beamlog_text()` instead streams the tar
(`tarfile.open(fileobj=resp.raw, mode='r|*')`) from the start and stops
as soon as a matching member is read, which avoids a full download
whenever the match isn't near the very end, but is still a
possibly multi-minute cost — a rejected alternative (estimating how many
bytes precede the beamlog by inspecting one tar's layout) does not
generalize, since different SBIDs' evaluation tars lay out their
contents differently.

Only **one** beamlog needs fetching per observation, not one per beam or
per Stokes parameter: confirmed by diffing all 36 beams' `beamlog`
files for Stokes Q channel-for-channel (all 288 channels byte-identical
across every beam — this is the direct meaning of "beam-wise matched"),
and separately confirmed that Stokes U's own `BEAMS` table (from a
downloaded cube) matches Stokes Q's beamlog exactly too, with 0
mismatches across all 288 channels — the restoring beam comes from the
(u,v) coverage and imaging weights for a channel, not from which Stokes
combination forms the image, so this is not specific to this one SBID.
`fetch_beamlog_text()` therefore always targets a Stokes-Q `.conv.`
beamlog, whichever of the 36 beam numbers is encountered first while
streaming.

**Auth gotcha, confirmed directly (cost one failed live test to catch):**
the resolved `sync?id=...` URL is itself a scoped, pre-authorized
credential — passing HTTP Basic Auth on top of it (the same
`casda._auth` tuple needed to fetch the DataLink XML one step earlier)
makes CASDA reject the request with `401` (`application/json` error
body), where the identical URL with *no* auth at all returns `200`
(`application/x-tar`). `fetch_beamlog_text()` takes no `auth` parameter
for exactly this reason; `resolve_datalink_sync_url()`, one step
earlier, still needs it.

The found beamlog's 0-indexed `"Target BMAJ" "Target BMIN" "Target
BPA"` columns (the per-channel beam the observatory already convolved
this SBID's `.conv.` cubes to, beam-wise) are reformatted by
`curate_beamlog()` into `convolve_cubes`' own `beamfiles=` ASCII
format — 1-indexed `channel bmaj_arcsec bmin_arcsec bpa_deg`
(`src/convolve_cubes.f90`'s `read_beams_ascii`, `cfg/example_beamLog.txt`)
— and written to `<obs_dir>/beamlog_target.txt`. No FITS-injection route
(appending a `BEAMS` table + `CASAMBM=T` back into the downloaded cube)
was considered further once the ASCII curation path was confirmed
sufficient for `convolve_cubes`' own `beamfiles=` mechanism.

**T1 (wiring `--fetch` output into `scripts/run_pipeline.sh` as a
`fetch` stage) not yet started.**
