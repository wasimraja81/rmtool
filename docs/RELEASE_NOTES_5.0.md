# Release Notes 5.0

Status: **in preparation on `multi-band-tomography`, not yet merged to
`develop`/`main`, not yet tagged.** These notes are written anticipating
that merge — see `CHANGELOG.md`'s own `[5.0]` entry and
`planning/MULTI_BAND_TOMOGRAPHY_PLAN.md` (tickets T0-T12) for the living
record this document summarizes.

## Summary

5.0 is a multi-band Faraday tomography milestone — the largest body of
work this project has shipped in one release. Real spectro-polarimetric
surveys routinely need to combine several frequency bands (different
receivers, different pointings, different epochs) into one RM synthesis
analysis, and doing that correctly requires three things that didn't
exist before this release: `rm_synthesis` itself needs to understand
"more than one input file", the bands need to actually be on the same
sky grid and the same angular resolution before combining them makes
scientific sense, and the output files need to honestly say what
resolution (or lack of one) went into them. This release adds all three:
multi-band support built directly into `rm_synthesis` (comma-list config
schema, frequency merge, full tiling/GPU/IO-parallelism support), two new
standalone preprocessing tools (`reproject_cubes` for sky-grid alignment,
`convolve_cubes` for resolution matching), and beam-metadata propagation
fixes to `rm_synthesis`'s own output headers.

Every ticket in this release was held to the same bar: for the
`nbands=1` (single-band) case, output must remain bit-identical to
before this branch existed, checked after every single change, with no
exceptions across the whole release. For genuinely new multi-band
behaviour, the bar was bit-identical against an independently
constructed reference wherever one could be built (e.g. a contiguous
channel split of an existing single-band cube must reproduce that same
cube's own result exactly), not merely "produces a plausible-looking
number."

## Highlights

### Multi-band RM synthesis, built into `rm_synthesis` itself

`infileQ`/`infileU`/`resiQ`/`slopeQ`/`resiU`/`slopeU`/`infileI`/`path_I`/
`badchan_file`/`chan_blc`/`chan_trc`/`chan_inc` all now accept
comma-separated lists — the number of bands is however many entries you
list, no separate `nbands` key, and a config with no commas anywhere is
completely unaffected. Every band's own geometry (RA/Dec WCS, NAXIS,
frequency axis) is validated against a `reference_band` before any
compute starts, loudly refusing on mismatch rather than silently
combining data that doesn't line up.

Once validated, every band's channels are concatenated into one merged
frequency/λ² spectrum and run through the same RM-synthesis compute path
a single-band cube already used — no separate multi-band numerical
kernel, no deduplication of overlapping channels, no sort requirement
(the underlying DFT sum doesn't care what order its terms arrive in).
This was verified end to end against a real thesis-published scenario
(P-band 300/30 MHz, L-band 1200/120 MHz, reproduced from Table 6.1): a
point source recovered accurately alone and combined, a Faraday-thick
component whose peak amplitude in its own RM window is revealed ~9x
larger when P and L are combined than at P alone (the "washout at one
band, revealed by combining bands" physics this whole effort exists to
demonstrate), and two close RM components resolved at the higher-resolution
band, blended at the lower-resolution one, both still visible when
combined.

What used to be single-band-only restrictions are now fully supported for
multi-band: multi-tile runs (not just images small enough for one RAM
tile), per-band channel sub-range selection (reject bad edge channels or
hand-pick a good range independently per band), per-band bad-channel
files, GPU offload (verified on real GPU hardware, including the VRAM
staging path), and `io_read_threads`/`io_overlap` I/O parallelism. Every
one of these had an initial "not yet implemented for multi-band" stop
when multi-band support first landed — each was investigated on its own
ticket (frequently by direct code inspection *before* writing any change,
on the discipline of "is this restriction actually necessary, or just
conservative"), and in most cases the underlying code turned out to
already be band-count-agnostic, needing only the blocking check removed
and a real verification run to prove it.

Since `use_auto_rm_range=1`'s existing heuristic isn't safe across
bands, a new diagnostic logs `δRM`, `max RM scale` (per band and
combined), and per-band `ΔRM` — informational only, doesn't auto-select
`beg_rm`/`end_rm`/`nrm` — verified against the same published Table 6.1
to within about 1%.

### `reproject_cubes` — cross-band sky-grid alignment

A new standalone tool (own binary, `make reproject_cubes`, independent
of the main `rm_synthesis` build) that reprojects two or more FITS cubes
onto one common sky grid, using Starlink AST for WCS handling and
`astResampleR` for the actual resampling. Three footprint modes
(`intersection`/`union`/`reference`); full header propagation, including
`CROTA`/`PCi_j`/`CDi_j` sky rotation (needed for real ASKAP data, not a
hypothetical case); `mem_frac_ram`-budgeted block I/O and OpenMP
parallelism across planes, mirroring `rm_synthesis`'s own tile-planning
concept.

### `convolve_cubes` — cross-band resolution matching

A second new standalone tool (own binary, `make convolve_cubes`)
convolving cubes — across one or several input files at once — to a
single common angular resolution, built on two new pure-computation
modules:

- **`gaussft_mod`** (`src/gaussft.f90`): FFT-domain deconvolve-then-
  reconvolve between elliptical-Gaussian PSFs, thread-safe for OpenMP via
  a plan-once/execute-many split (FFTW's own planner isn't thread-safe;
  its "new-array execute" form is, verified directly — 16 threads sharing
  one plan reproduces a serial run bit-for-bit).
- **`commonbeam_mod`** (`src/commonbeam.f90`): finds the smallest beam
  every one of N per-channel PSFs can be deconvolved from (real ASKAP
  per-channel position angle varies by more than 90 degrees across a
  band, so "just take the largest beam" is not generally correct).
  Follows the same approach CASA and the `radio_beam` Python package
  use, verified directly against `radio_beam` 0.3.9 on a real 286-channel
  ASKAP beam table.

`convolve_cubes` reads per-channel beams from a CASA-style `BEAMS`
binary table (auto-detected) or a portable ASCII/CSV beam log — see
`cfg/example_beamLog.txt`/`.csv` for ready-to-adapt examples, so the
format doesn't need to be reverse-engineered from source. A channel is
treated as bad (its output plane written all-NaN, not convolved, and
automatically excluded by `rm_synthesis`'s own NaN detection downstream)
if it's missing from the beam source entirely, or listed with BMAJ or
BMIN equal to 0.

### `rm_synthesis` beam-metadata propagation

Previously, none of `rm_synthesis`'s 8 output products (AMP/PHA cubes,
mask, nvalid, and — when `cubestat=y` — the peak/rmpeak/angpeak/snr maps)
carried any beam information at all. Now `BMAJ`/`BMIN`/`BPA` propagate
from the input Q cube to every one of them. If the input still has
`CASAMBM=T` (a genuinely per-channel-varying beam that hasn't been run
through `convolve_cubes` yet), that propagated scalar means nothing on
its own — so the flux-derived outputs also get `CASAMBM=T` plus the
input's own real per-channel `BEAMS` table attached as an extension,
plus a `HISTORY` note explaining why, so a user who notices the
unexpected extension is prompted to ask exactly the right question
("have we processed this correctly?"). The mask and valid-channel-count
outputs deliberately don't get this — they're validity bookkeeping, not
flux data, and a beam extension there would only be confusing. In
multi-band mode, every band's own beam metadata is now cross-checked
against the reference band's, with a runtime warning (not a hard error)
on any mismatch.

## Validation

- All 4 build flavours (`scratch/make_all.sh`) clean throughout; full
  `tests/run_tests.sh` grew from 28/28 at the start of this branch to
  49/49 by the end, re-run clean after every change.
- The `nbands=1` bit-identical sweep (140/140 FITS outputs against a
  frozen pre-branch baseline) was checked after every single ticket in
  this release with no exceptions — the primary correctness instrument
  for a branch whose own design plan explicitly gave up "bit-identical
  by construction" for the genuinely new multi-band code paths.
- Multi-band-specific correctness gates, each verified against a
  purpose-built reference rather than "runs without crashing": a
  contiguous channel split of an existing single-band cube reproduces
  that cube's own undivided result bit-identically; a multi-tile run of
  multi-band data reproduces the single-tile result of the same data
  bit-identically; per-band channel sub-range selection composed with
  the channel-split test reproduces the same known-good reference;
  per-band bad-channel flagging of a given raw channel reproduces
  flagging that same channel via the single-band mechanism.
- `reproject_cubes`, `gaussft_mod`, `commonbeam_mod`, and `convolve_cubes`
  each independently verified against real ASKAP data, independently
  computed ground truth (Python/astropy, `radio_beam`), or both — see
  `planning/MULTI_BAND_TOMOGRAPHY_PLAN.md` tickets T10-T11 for the full
  evidence trail.
- `rm_synthesis` beam-metadata propagation verified by injection (real
  BMAJ/BMIN/BPA, `CASAMBM=T`/`BEAMS`, mismatched multi-band beams) and by
  a real end-to-end run: a `convolve_cubes`-produced NaN bad-channel
  plane fed into `rm_synthesis` with no `badchan_file` at all was
  correctly excluded automatically via `rm_synthesis`'s existing
  (default-on) NaN detection.

## Compatibility and behaviour notes

- Every existing single-band config file is unaffected — the comma-list
  schema is purely additive (no commas anywhere behaves exactly as
  before), and every new multi-band code path is gated behind "more than
  one band" checks that are dead code for `nbands=1`.
- `rm_synthesis`'s Q/U/I/mask input cubes are now opened `READONLY`
  instead of `READWRITE` (a fix, not a new feature — see `CHANGELOG.md`).
  No behavioural change for any correct usage; closes a real, if latent,
  risk to input data.
- Two new standalone binaries (`reproject_cubes`, `convolve_cubes`) with
  their own extra build dependencies (Starlink AST, FFTW3) — the main
  `rm_synthesis` build and its dependencies (gfortran, CFITSIO) are
  unaffected; see `BUILD.md`.
- Output header changes: `BMAJ`/`BMIN`/`BPA` (and, when applicable,
  `CASAMBM`/a `BEAMS` extension/`HISTORY` cards) may now appear on
  `rm_synthesis` outputs where they never did before, if the input cube
  itself carries this metadata. No existing header keyword was removed,
  renamed, or changed in meaning.
- No cfg keys were removed or renamed. New per-band keys
  (`chan_blc`/`chan_trc`/`chan_inc`) are optional; `badchan_file` remains
  required per band (same requirement it already had for single-band).

## What shipped in this release

- Multi-band `rm_synthesis`: comma-list config schema, N-band geometry
  validation, frequency/λ² merge, RM-range/resolution diagnostic,
  multi-tile support, per-band channel sub-range selection, per-band
  bad-channel files, GPU offload, `io_read_threads`/`io_overlap` support.
- `src/reproject_cubes.f90`, `src/gaussft.f90`, `src/commonbeam.f90`,
  `src/convolve_cubes.f90` — new standalone tools and modules.
- `cfg/example_beamLog.txt`/`.csv` — ready-to-adapt beam-log examples.
- `rm_synthesis` beam-metadata propagation and input-safety fixes.
- Full documentation pass: README, `docs/ARCHITECTURE.md`, `BUILD.md`,
  `QUICKSTART.md`, `cfg/CONFIG_README.md` all updated to cover the new
  toolchain (previously undocumented anywhere, including the earlier
  `reproject_cubes` work).

## What's next (beyond this release)

- A full run of the preprocessing toolchain against the complete 23GB
  real ASKAP cube this work targets — only cutouts and synthetic data
  verified so far.
- Merging `multi-band-tomography` into `develop`, then `main`, and
  tagging `5.0` there — not yet done as of this document.
