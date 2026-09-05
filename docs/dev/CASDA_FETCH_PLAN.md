# Fetching calibrated Q/U multiband cubes from public archives — Plan

**Status: T0 (discovery, fetch, beam curation), T1 (pipeline-cfg
generation), T1a (band/SBID/version selection), and T1b (RM-range
policy) implemented, not yet committed.**

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
absence on its own does not distinguish "pilot-era, no
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

### T1 — generating a ready-to-run pipeline cfg from `--fetch` output

**In progress, not yet reviewed/committed.** Rather than a `fetch` stage
built into `scripts/run_pipeline.sh` itself, `--fetch` now also calls
`generate_pipeline_cfg()`, which writes a flattened `pipeline_staging/`
symlink directory (needed because `run_pipeline.sh`'s own
`match_infileQ=`/`match_infileU=` are bare filenames under one shared
`match_input_path=`, and its symlink step doesn't create intermediate
directories — a per-SBID subdirectory in the filename, matching this
tool's own `<outdir>/SB<n>/` layout, isn't an option), an `rm_synthesis`
cfg template, and a top-level `run_pipeline.sh`-compatible pipeline cfg
chaining `match->rmsynth` — so `./scripts/run_pipeline.sh
<generated>.cfg` runs end to end with no manual editing. `rmclean` is
deliberately excluded from the generated `stages=`: its CLEAN stopping
criteria are a choice for the user to make (see
`docs/user/EXAMPLES.md`'s own dedicated section on this), not something
to default silently.

**Correction (2026-09-05): the template's `use_auto_rm_range=1` default
was wrong, not just a stylistic choice to revisit.** `use_auto_rm_range=1`
is hard-rejected by `rm_synthesis` itself for any multi-band run
(`src/rm_synthesis.f90:1002-1008`, `'ERROR: use_auto_rm_range=1 is not
supported for multi-band runs'`) — a generated multi-band cfg using it
would fail outright, not just produce a suboptimal range. The reason is
load-bearing, not just conservative: the auto-range heuristic
(`src/rm_synthesis_mod.f90:774-799`) infers the dataset's min/max
frequency from a merged λ² array's own endpoints, which only holds for
one internally-monotonic band — concatenating bands in any order but
frequency order would silently corrupt it. See "RM range and `nrm`"
below for the corrected default.

`match_args`' own `stages=` is chosen from how many observations were
fetched: a single observation uses `stages=convolve` (nothing
to reproject — one observation's own Q/U already share one native
grid); two or more use `stages=both` with `footprint_mode=intersection`
(separate SBIDs are separate scheduling blocks and may not share a grid
even at the same frequency band, not just across distinct bands).

**Open, scoped but not yet implemented (user note, 2026-09-05):**
CASDA's public holdings for a given position can change over time —
multi-band data may simply not exist yet even when the position is
otherwise usable. The generator should not just silently build whatever
`stages=`/band combination happens to fall out of one fetch run's
result. Needed: let the user choose, after seeing what was
fetched, whether to (a) proceed single-band if multi-band isn't
available, (b) proceed with a specific subset of bands even when more
are available (e.g. run single-band first before committing to full
multi-band tomography), or (c) hold off entirely. Current behaviour
(this session) is fully automatic over whatever `usable` contains —
correct as a mechanism (single vs multi cfg shape does need
to track the runtime fetch result, confirmed by the user), but the
*selection* of which fetched observations feed the generator should
become a user decision, not an automatic "use everything usable."
`generate_pipeline_cfg(outdir, run_name, records)`'s `records` parameter
is already a plain, explicit list (not derived internally from global
state), so filtering it by user choice before the call is the natural
extension point — no rework of the generator itself anticipated.

**This turned out to need its own analysis pass, confirmed via a live
CASDA query, not just a design opinion** — see below.

### T1a — SBID and version selection (confirmed via a live CASDA query)

**The "more than one SBID per band" case is not an edge case — it's the
normal case for a target with observation history.** Querying
PKS 2130-538 live (`cont.restored.3d`, all 12 distinct `obs_id`s CASDA
currently has for this position) shows the 799.5-1087.5 MHz band alone
has **three** separate, fully Q/U-`.conv`-complete SBIDs: `ASKAP-45811`
(t_min MJD 59907.26), `ASKAP-51797` (MJD 60157.55), `ASKAP-74876`
(MJD 60874.58) — three re-observations/reprocessings of what tomography
needs to treat as ONE band, not three independent ones.
`generate_pipeline_cfg()` as committed today would treat all three as
separate "observations" and feed all three into
`match_cubes`/`rm_synthesis` as if they were three distinct bands —
silently wrong, confirmed against this example, not hypothetical.

**Two independent selection problems, both confirmed:**

1. **Which SBID, when 2+ share a band.** `t_min`/`t_max` (MJD) are
   standard ObsCore columns, confirmed present and populated for every
   row queried above — available from the discovery query alone, no
   download needed, so a "latest observation" choice is fully knowable
   in dry-run mode. User's proposed default: pick the SBID with the
   latest `t_min` per band.
2. **Which file version, when the SAME SBID has more than one cube.**
   Confirmed (not hypothetical): `SB10083` (same field) has both
   `image.restored.i.SB10083.contcube.fits` and
   `image.restored.i.SB10083.contcube.v2.fits` alongside it — a `.vN.`
   reprocessing tag. Notably, the tag is NOT applied uniformly across
   that SBID's whole file set in this example (`v2` exists for the
   plain Stokes-I file but not for its `.conv` counterpart, nor for
   Stokes V at all) — so version selection has to be resolved per
   (SBID, Stokes letter, conv-status) file group, not once per SBID.
   User's proposed default: highest `vN` wins per file group; an
   unversioned filename is the implicit baseline (`v1`).

**Implemented (2026-09-05), confirmed with the user before building:**
`--fetch` is retired in favor of `--run-mode` (`dry` default, `auto`,
`select`):

- `group_into_bands()` groups usable observations by exact frequency
  match, numbered 1..N ascending, each band's own SBIDs sorted
  latest-`t_min`-first (index 0 = the `auto` pick). `print_band_report()`
  (the `--run-mode=dry` output, always printed regardless of mode) lists
  every band, every candidate SBID with its observation date, and any
  (Stokes letter, conv-status) file group with more than one `.vN.`
  version present, flagging the default in both cases.
- `--run-mode=auto`: every band, each one's latest-`t_min` SBID, every
  file's latest version -- "fix ALL cfg params to auto" per the user's
  own framing.
- `--run-mode=select --select="<list>"`: `<list>` is a comma-separated
  `band[:sbid[:version]]` list (band numbers from a prior dry run);
  omitting `sbid`/`version` defaults that axis to latest. Resolved by
  `resolve_selection()`, validated against band/SBID numbers derived
  fresh from THIS run's own query, never a cached mapping (dry-run band
  numbers can drift between runs if CASDA's holdings change in between
  -- not otherwise guarded against beyond showing the frequency range
  next to every number so a mismatch is visible).
- Version selection threads into `select_iquv_mask(rows, version=...)`,
  resolved **independently per Stokes letter, never shared between Q
  and U the way `.conv` is** (correction, 2026-09-05: an earlier version
  of this hard-aborted the whole observation if a requested version
  wasn't present on both Q and U, which the user caught as wrong -- Q
  and U are recalibrated independently, e.g. a fix that only touches Q
  need not touch U, so requiring the same version number on both is not
  a legitimate requirement the way `.conv` is). A letter
  missing the requested version falls back to its own highest available
  version, no warning -- expected, not an anomaly. `.conv` presence
  remains the only cross-Stokes hard, all-or-nothing gate.
- This is scoped to *which observation feeds the pipeline cfg*, not the
  raw per-observation fetch mechanics, which are otherwise unchanged.

### T1b — RM range and `nrm`: implemented (2026-09-05)

Per the correction above, `use_auto_rm_range=1` cannot be the generated
default for multi-band. Implemented, checked against this codebase's own
existing (informational-only) multi-band diagnostic rather than derived
from a textbook formula in isolation (`src/rm_synthesis.f90:2264-2278`),
confirmed with the user before building:

- `beg_rm = -1000.0`, `end_rm = 1000.0` (a Galactic-scale range, fixed
  rather than data-derived) -- applied **uniformly, including
  single-band fetches**, one consistent generated-cfg policy regardless
  of band count (the user's own reasoning: not wanting to spend time
  processing RM values outside the range of interest either way, even
  though `use_auto_rm_range=1` is only *hard-rejected* for `nbands>1`).
- `nrm` computed by `compute_rm_range()` from the **combined δRM** (RMSF
  FWHM) the fetched bands' own λ² coverage dictates, reading each
  record's own downloaded FITS header (`band_lambda_sq_span()`: FREQ
  axis `CRVAL`/`CDELT`/`CRPIX`/`NAXIS`, confirmed Hz-scale directly
  against a downloaded cube's header) -- not an ObsCore column, which
  only carries the total span, not the per-channel resolution needed
  here. Uses the SAME formula the existing multi-band diagnostic already
  prints as informational output (not a new derivation): per band,
  `dlam2_band = lam2_lo_edge - lam2_hi_edge` (half-channel-edge-extended
  λ² span); combined `delta_RM = cfg.fac / sum(dlam2_band over all
  bands)` (`src/rm_synthesis.f90:2264-2278`, `cfg.fac` = π by default).
  `nrm = ceil((end_rm - beg_rm) / combined_delta_RM)` -- matches this
  project's own multi-band cfg
  (`cfg/rmsynth-multiband-wallaby-emu.cfg`: `nrm=29` for a ±500 range
  with δRM=35.70 rad/m², and 1000/35.70≈28.0, consistent with one nrm
  sample per resolution element before oversampling). `ofac` stays at
  its existing default (`4`) as the oversampling step on top of `nrm`
  (`nrm_out = nrm * ofac`, already how every cfg in this repo uses it)
  -- confirmed with the user that "oversample δRM by some factor" reads
  as `ofac`'s existing role, not a second independent factor.
- Sanity-checked against a downloaded SB74876 Q cube (288 channels,
  1 MHz/chan, ~800-1088 MHz): δRM≈48.6 rad/m² -> `nrm=42` for the
  ±1000 range -- larger (coarser) than the 2-band WALLABY+EMU example's
  35.70 rad/m², consistent with a single, narrower band having less
  total λ² coverage than a combined multi-band span.

**Still open:** `rmclean_cubes` isn't wired into `--run-mode=auto`/
`select` at all yet (deliberately -- see the `generate_pipeline_cfg()`
docstring); the user's SBID-selection framework ("quit and let user
specify... latest processing could be the default") should extend to it
the same way once it's added, rather than being designed separately.

### T1c — three corner cases found by direct review, implemented (2026-09-05)

Asked to precisely enumerate corner cases in fetching a small cutout
set for tomography and say whether the implementation was immune to
each. Three were not, and are now fixed:

1. **Two SBIDs sharing a band could be different sky fields (adjacent/
   overlapping mosaic tiles), not reprocessings of the same field** --
   `Observation`/`Band` never captured or compared field identity, only
   frequency + observation time. Resolved by checking target
   encompassment directly rather than field identity: CASDA's
   `query_region` cone search only guarantees a field is *near* the
   search circle, not that the target position is within its
   coverage -- confirmed live (SB74876/EMU_2132-51 has `s_fov=54.7 deg`,
   a wide mosaic, with the target 2.5 deg off its own field centre --
   comfortably covered here, but a field that only clips the search
   radius at its edge would not be). `observation_covers_target()`
   checks `separation(target, field_centre) + cutout_radius < s_fov/2`
   using ObsCore's own `s_ra`/`s_dec`/`s_fov` (a circular approximation
   of the field's possibly-irregular polygon footprint in
   `s_region` -- a sanity check, not a precise containment test).
   `filter_target_coverage()` excludes (with a printed reason) any
   candidate that fails this, before band grouping ever sees it. Once
   every remaining candidate reaches the target, "latest
   wins" is the right default, per the user's own reasoning.
2. **Two legitimately different bands could overlap in frequency**, and
   separately, **the same band could be split into two by
   floating-point noise** in `em_min`/`em_max` between SBIDs (`Band`
   grouping keys on exact tuple equality). Both would silently
   double-count the shared frequency range in a combined multi-band
   run, biasing the RMSF/sensitivity. One check handles both (a
   near-duplicate range trivially overlaps almost entirely too):
   `check_band_overlaps()` runs on the final resolved selection, right
   before fetching, and quits with a message naming the overlapping
   pair and by how much if any two selected observations overlap.
   Resolving it (picking only one) is left to the user.
3. **A partial multi-band fetch (one requested observation fails or
   its beamlog can't be curated) used to still generate a pipeline cfg
   for whatever succeeded**, silently downgrading a stated multi-band
   intent (e.g. from `--select "1,2"`) to single-band. Fixed:
   `main()`'s fetch loop now tracks `all_succeeded` across every
   requested observation (a `SKIPPED`/`FAILED`/`PARTIAL` cutout status,
   or a beamlog curation failure when one was needed, each
   set it `False`) and, after attempting every requested observation
   (not stopping at the first failure -- confirmed with the user:
   complete diagnostics in one pass are worth the extra time on a run
   that's going to be rejected anyway), refuses to generate a pipeline
   cfg at all if anything short of full success occurred, printing
   which observation(s) didn't make it and why instead. This also
   closes an independently-existing bug: a failed beamlog
   curation used to silently fall back to `beamfiles=auto` in the
   generated cfg even when a `BEAMS` table was already confirmed
   absent -- guaranteed to fail at `match_cubes` runtime rather than
   inform the user immediately.

**One implementation-level fix along the way**: `observation_covers_target()`
originally returned `numpy.bool_` (from `SkyCoord.separation()`), not
Python's native `bool`, despite its own declared `-> bool` signature --
harmless under a plain `if`, but wrong under a strict `is True`/`is
False` check, which is how the mismatch was caught. Wrapped the return
in `bool(...)`.

The two accepted-as-is items from this same review: no validation that
a fetched cutout has usable signal (Stokes I is exactly why it gets
fetched -- the user inspects it themselves), and `rmclean_cubes` still
not wired into `--run-mode` (tracked above, unchanged by this pass).

### T1d — live-run fixes (2026-09-05): a starter rmclean cfg, and absolute paths

Two issues surfaced from the first live `--run-mode=auto` run against
`dancingghosts` (SBID 74876, single band), both fixed:

1. **A starter `rmclean_cubes` cfg is now generated alongside the other
   two, even though it isn't chained into `stages=` by default** --
   `write_rmclean_template()`, using a single-band structure (deliberately
   omits `min_valid_chan_frac`,
   which only matters for a multi-band merge's spatially-varying
   combined-band coverage, per `cfg/rmclean-multiband-wallaby-emu.cfg`)
   with two user-chosen values: `abs_flux_floor=20uJy`, `nwriters=2` --
   both flagged in the generated file's own comment as starting points
   to re-check against this run's own dirty cube, not copied
   measurements from either existing project cfg. The pipeline cfg's
   own header comment and `casda_fetch.py`'s final terminal message
   both now name this file's exact path directly, so extending to
   CLEAN means editing that file and adding two lines to the pipeline
   cfg, not finding a template to copy from `cfg/` first.
2. **A relative `--outdir` broke `scripts/run_pipeline.sh`'s own
   match-stage symlinking.** That script builds its own second symlink
   layer (`match_input_symlinks/`) via `ln -s "${f}" "${link}"`, where
   `${f}` comes from this tool's own `match_input_path=`. A relative
   `${f}` makes `ln -s` write a relative target, which POSIX resolves
   against the *new* symlink's own directory (`match_input_symlinks/`),
   not any directory `casda_fetch.py` or the user had in mind --
   producing `ERROR: cannot open FITS file` once `match_cubes` tried to
   read through it, confirmed directly against the user's own live run
   output. Independent of `stage_symlink()`'s own symlinks, which
   already used an absolute target. Fixed at the source: `args.outdir`
   is absolutized (`os.path.abspath`) once, immediately after
   `parse_args()`, so every path derived from it downstream (`obs_dir`,
   `pipeline_staging/`, and all three generated cfgs' own path-valued
   keys) is absolute regardless of which directory `casda_fetch.py` is
   invoked from.

### T1e — re-fetch avoidance and UX fixes from a second live run (2026-09-05)

1. **`--if-exists=check|refetch|reuse` (default `check`)** -- the fetch
   loop previously had no way to know a re-run's expected output files
   (predicted from `conv_by_letter`, already known pre-fetch -- the
   only checkable signal, since CASDA's cutout service renames every
   download and a version tag isn't recoverable from a downloaded file
   either, see `version_number_from_filename`) already existed on disk,
   so re-running after fixing the path bug above would have
   unconditionally resubmitted a CASDA cutout job and re-downloaded
   everything. Presence-only, not correctness -- deliberately not
   checked against version (per the user's own framing: "we do not
   technically need to check for correctness"). `check` (default)
   reports what's already there and aborts rather than silently guess
   either way; `reuse` skips the download and treats the existing files
   as this run's own output; `refetch` is today's unconditional
   behaviour, unchanged.
2. **A "stale symlink" error surfaced from `run_pipeline.sh` itself**
   on the corrected second attempt -- not a new bug: `match_input_symlinks/`
   is a directory `run_pipeline.sh` owns and manages entirely on its
   own (unrelated to this tool's own `pipeline_staging/`), and it still
   had a symlink left over from the FIRST (pre-fix) attempt, whose
   target canonicalizes to a dangling path under the old relative-path
   scheme. `run_pipeline.sh`'s own guard against silently repointing an
   existing symlink is doing exactly what it's supposed to; the fix is
   just removing that stale, fully-regeneratable directory before
   re-running -- not a `casda_fetch.py` change.
3. **A progress message added before the `casda.cutout()` call**, which
   stages the request server-side and polls until the job completes --
   previously silent for however long that takes (confirmed: upwards of
   a minute per observation), which could read as a hang.
4. **`scripts/run_pipeline.sh`'s own stale-symlink error message was
   ambiguous about which of its two paths to delete** ("stale symlink
   at X does not point at Y -- remove it first" -- "it" reads as
   referring to either). Reworded to name both paths explicitly next to
   the instruction: remove `${link}` (`run_pipeline.sh`'s own, gets
   recreated automatically), never `${f}` (the input data). Not a
   `casda_fetch.py` file, but the same live-run session surfaced it.

### T1f — two more `run_pipeline.sh`/rm_synthesis fixes from continued live runs (2026-09-05)

1. **`cubestat` in the generated rmsynth template changed from `n` to
   `y`.** The user asked why their run produced no PEAK/RM_PEAK/
   ANG_PEAK/SNR maps -- confirmed against source
   (`src/rm_synthesis_mod.f90:1564-1679`, gated only by `cfg%cubestat`,
   `src/rm_synthesis.f90:4834`) that `cubestat=n` was the sole cause,
   and that this was the generator's own unprompted default, never
   confirmed with the user the way the RM-range/rmclean values were.
   Peak/SNR maps are exactly the kind of first-look diagnostic this
   tool's own users want, so the default flips to `y`.
2. **`run_pipeline.sh` now clears only symlinks (not regular files) from
   `match_input_symlinks/` at the start of its own match stage**, so a
   symlink left over from an earlier, separate invocation (confirmed
   case: one built from a relative `match_input_path` canonicalizes
   to a stale, dangling target once a later run supplies an absolute
   one, tripping the "stale symlink" guard even though nothing is
   wrong) no longer blocks a clean retry. Deliberately scoped
   to symlinks only (`find ... -type l -delete`, verified empirically
   to leave regular files untouched) -- NOT a blanket `rm -rf` of that
   directory, since `convolve_cubes`/`match_cubes` write their own
   output (`<symlink_name><outsuffix>`) alongside these same symlinks,
   and that output is computed data, not ephemeral infrastructure.
   User's own explicit call: auto-clean the symlinks, but leave the
   output-already-exists case exactly as it was (an informative error,
   manual removal) -- "we have provided enough info on what to do to
   the user."

### T1g — `--if-exists=clean`: a fresh start, opt-in (2026-09-05)

The user's own request: a way to fetch-and-run completely from scratch,
overwriting existing links/outputs, without hunting through several
directories by hand each time (fetch output, `pipeline_staging/`, the 3
generated cfgs, and everything `scripts/run_pipeline.sh` itself would
have produced against those cfgs). Added as a 4th `--if-exists` choice
rather than a separate flag, since it's the same underlying question
("what to do about existing state") just answered more aggressively.

`clean_pipeline_artifacts(outdir, run_name)` removes, by exact/
predictable name only -- never a blanket wipe of `outdir` itself:
`pipeline_staging/`, `match_input_symlinks/` (here, unlike
`run_pipeline.sh`'s own default per-run symlink-only cleanup above,
the *whole* directory including any `_CONV.FITS` leftover -- this is
an explicit, opt-in "start over" request, not `run_pipeline.sh`'s own
quiet default behaviour, so full thoroughness is the right call here
and doesn't contradict that earlier, more conservative decision), the
3 generated `*_pipeline.cfg`/`*_rmsynth.cfg`/`*_rmclean.cfg` files, and
a glob on `<run_name>.*` (catches rm_synthesis' own AMP/PHA/MASK/
NVALID/PEAK/RM_PEAK/ANG_PEAK/SNR outputs and `run_pipeline.sh`'s own
`<run_name>.provenance/` directory in one pattern, since all of them
share that one dot-prefixed naming convention -- confirmed directly
from `run_pipeline.sh`'s own `RMSYNTH_OUTFILE="${OUTDIR}/${RUN_NAME}"`
and its printed `"[pipeline] Provenance: <run_name>.provenance"` line
-- rather than hardcoding every exact suffix, so this stays correct
regardless of `output_mode`/`cubestat`/etc.). Each selected
observation's own `<outdir>/SB<n>/` directory is removed separately,
per observation, inside the fetch loop itself (not this function's
job), since only `clean_pipeline_artifacts` needs `run_name` and only
the fetch loop knows which specific SBIDs are about to be re-fetched.

Verified against a fixture mirroring a full run's own artifact set
(all of the above, plus an unrelated file and a `dancingghosts2_*`
cfg sharing a name prefix): removed exactly the 10 expected items,
left the unrelated file and the prefix-sharing-but-distinct cfg
untouched, confirming the glob doesn't over-match.

### T1h — `--chain-rmclean`: opt-in, with an explicit caveat (2026-09-05)

The user's own request, after discussing the WALLABY+EMU/single-band
rmclean cfgs' shared three-way stopping strategy (`auto_nsigma=1.0`
per-pixel adaptive threshold, plus `abs_flux_floor` as a fixed backstop
for pixels where that estimate is biased, plus `niter` as a hard cap --
see T1d above): "I think we should allow --chain-rmclean for a user to
specify. We should print out what the lacunaes are of that chaining
without looking at dirty cubes, and that the numbers are based on
typical EMU fields." The "EMU fields" framing is confirmed accurate:
the user separately confirmed that single-band dataset (the source of
the `abs_flux_floor=20uJy`/`niter=500`/`gain=0.1` values baked into
`write_rmclean_template`) is EMU-only band, even though that isn't
documented anywhere in its own cfg files (checked its own `.fullim`/
`.e2e` cfgs directly -- see this project's own memory notes for this
fact's provenance).

`generate_pipeline_cfg()` takes a new `chain_rmclean: bool` parameter
(default `False`, so every existing call site keeps today's behaviour
unless it opts in). When set: `stages = match,rmsynth,rmclean` and
`rmclean_cfg_template = <path>` are emitted instead of the default
`stages = match,rmsynth` with no rmclean wiring; the generated pipeline
cfg's own header comment is replaced (not appended to) with a caveat
enumerating exactly what hasn't been checked before CLEAN runs:

- this run's own dirty cube (`<run_name>.PEAK.MAP.FITS`/
  `.SNR.MAP.FITS`) has not been inspected to judge whether the
  criteria in the generated rmclean cfg suit it, since CLEAN starts
  immediately after rmsynth in the same invocation;
- `abs_flux_floor`/`niter`/`gain` are carried over from another
  EMU-band field's own measured dirty-cube noise properties, not this
  run's own data;
- `auto_nsigma=1.0`'s per-pixel estimate is confirmed biased for some
  pixels (per that same cfg's own validation notes), so how well the
  fixed floor suits this field matters most exactly where the
  per-pixel estimate is least dependable.

`main()`'s own terminal summary at the end of a successful fetch prints
the same caveat when `--chain-rmclean` was given, and points to
`--chain-rmclean` itself in the closing suggestion when it wasn't.

Verified with a synthetic two-branch test (`chain_rmclean=False`/
`True` against the same fixture records): confirmed the default branch
has no active (non-comment) `rmclean_cfg_template` line and keeps
`stages = match,rmsynth`, and the chained branch emits exactly
`stages = match,rmsynth,rmclean` plus the correct
`rmclean_cfg_template=` path, both by checking active (non-`#`) lines
specifically -- the default branch's own header comment mentions
`rmclean_cfg_template =` as instructional text, so a plain substring
check would have been a false negative for the "not wired in" case.
