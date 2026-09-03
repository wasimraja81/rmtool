# App Reference

Complete, current parameter reference for every command-line tool in
rmtool: what each one does, exactly how to invoke it, every parameter
it accepts (name, default, required/optional, meaning), the output
files it produces, and one runnable example.

This document is the single source of truth for "what does parameter X
do" across all 5 tools. `README.md` and `QUICKSTART.md` give shorter,
task-oriented walkthroughs and link here for the full detail;
[docs/user/EXAMPLES.md](EXAMPLES.md) walks through choosing and combining
these parameters for real-world scenarios (multi-band, RM-CLEAN
stopping criteria, memory tuning, and more). `--help` on any tool
prints the same information this document does, generated from the
same source. If the two ever disagree, trust `--help` (it's generated
from the running binary) and please report the doc as stale.

## Contents

1. [rm_synthesis](#1-rm_synthesis)
2. [reproject_cubes](#2-reproject_cubes)
3. [convolve_cubes](#3-convolve_cubes)
4. [match_cubes](#4-match_cubes)
5. [rmclean_cubes](#5-rmclean_cubes)

---

## Shared Parameters

**`mem_frac_ram`** (all 5 tools): a fraction of the machine's **total**
system RAM.

---

## 1. `rm_synthesis`

The main RM-synthesis tool. Given Stokes Q and U spectral FITS
cubes (optionally several bands at once — "multi-band tomography"), it
performs pixel-by-pixel Faraday rotation-measure synthesis (a
Fourier-like transform from λ² to Faraday depth φ), producing RM cubes
of linear polarized intensity and polarization angle (or real/
imaginary) as a function of RA, Dec, RM. It supports bad-channel
removal, subimage/sub-channel extraction, automatic or manual RM-grid 
selection, mask/peak-map diagnostics, GPU offload, and a memory-budgeted, multi-threaded, tileable I/O pipeline.

### Invocation

```bash
bin/rm_synthesis <cfgfile> 
```

Purely positional — **not** `key=value` on the command line. `<cfgfile>`
is a required path to a `KEY=VALUE` text config file (one key per line,
`#` for comments; strict parser — duplicate keys, unknown keys, and
unparsable values are all hard errors). See 
[cfg/rmsynth.cfg](../../cfg/rmsynth.cfg) for an annotated template.

### Config keys

**Input / output paths (required):**

| Key | Default | Meaning |
|---|---|---|
| `path` | — | Directory containing the input Q/U cubes. |
| `infileQ` | — | Stokes-Q input filename, relative to `path`. Comma-separated for multi-band (e.g. `bandA.fits,bandB.fits`); band count is derived from this list's length. |
| `infileU` | — | Stokes-U input filename; must have the same band count as `infileQ`. |
| `reference_band` | `1` | Which band (1-indexed) is treated as the primary band in multi-band mode — see callout below. |
| `outfile` | — | Output basename; cube-type suffixes are appended automatically. |

> **`reference_band`:** `rm_synthesis` expects matched
> input across bands -- same geometry, and same resolution across the
> full spectral axis. Geometry (sky-axis WCS) is checked strictly: a
> mismatch between any band and the reference band, or between either
> band's own Q and U cubes, is a hard error -- pixel scale, reference
> position, CTYPE, pixel-grid rotation, and the celestial reference
> frame (RADESYS/EQUINOX/LONPOLE/LATPOLE) are checked. Matching resolution across
> frequency channels -- within a band, and across bands -- is the
> user's own responsibility. `rm_synthesis` only warns on a mismatch
> triggered when `CASAMBM=T` within a band (see "Beam metadata" under Output files), or when BMAJ/BMIN/BPA across bands disagrees. 
> `reference_band` picks whose beam values get written to the output
> headers regardless of whether that check passed. 
> Note: If the beams don't match, that propagated metadata is 
> correspondingly less useful as a description of the run -- and more so in the 
> multi-band mode, where only the reference band's own 
> per-channel beam data is ever captured in the output at all; a
> mismatching non-reference band's own beams are never recorded
> anywhere in the output, only flagged in the run log.
> `infileQ`/`infileU`/`infileI`, `badchan_file`, and the channel-window
> limits use each band's own value regardless of `reference_band`.

**Bad-channel handling (optional):**

| Key | Default | Meaning |
|---|---|---|
| `badchan_file` (alias `global_badchan_file`) | omitted (no removal) | One channel index per line. Omit this key entirely for no bad-channel removal. When given, comma-separated per-band list for multi-band — every band needs its own entry: a real file path, or the literal value `none`. A blank entry is a hard parse-time error once this key is given at all. |

**Subimage extraction (`subim` required; the rest optional, default = full cube):**

| Key | Default | Meaning |
|---|---|---|
| `subim` | — | `y`: only process the ranges below; `n`: full cube (ranges ignored). |
| `subim_ra_blc` / `subim_ra_trc` / `subim_ra_inc` | `1` / `0` (=max) / `1` | RA first pixel (≥1) / last pixel (0=max) / step (≥1). |
| `subim_dec_blc` / `subim_dec_trc` / `subim_dec_inc` | `1` / `0` (=max) / `1` | Same, for Dec. |
| `subim_chan_blc` / `subim_chan_trc` / `subim_chan_inc` | `0` (=first) / `0` (=max) / `1` | Same, for channel; per-band comma list allowed. |
| `subim_parfile` | `subimage.par` | **Deprecated, no-op.** Accepted by the parser for backward compatibility but never read — superseded entirely by the `subim_*` KEY=VALUE keys above. Setting it has no effect. |

**Q/U processing & bias correction: TODO (not implemented)**

> `resiQ`/`slopeQ`/`resiU`/`slopeU` are parsed,
> required, and stored, but never combined with the I-cube data in a
> bias-correction computation anywhere in
> `rm_synthesis`/`rm_synthesis_mod` — the I-cube is read into memory
> when `remove_qu_bias=y` but never applied to Q/U afterward. Multi-band
> runs (more than one band) with `remove_qu_bias=y` hit an explicit
> hard stop (`ERROR: multi-band Q/U bias correction is not yet
> implemented.`). Single-band `remove_qu_bias=y` has no such stop but
> also performs no correction — the feature is a no-op today regardless
> of `remove_qu_bias`'s value. This is however left as a TODO that will allow 
> for example, leakage correction, as well as tomography on suitably normalised 
> polarised fractions. 

| Key | Default | Meaning |
|---|---|---|
| `rem_mean` | — (required) | `0`/`1`: subtract per-pixel mean from Q/U before synthesis. |
| `remove_qu_bias` | — (required) | `y`/`n`: apply I-cube-based bias correction (needs `path_I`/`infileI`). |
| `resiQ` / `slopeQ` / `resiU` / `slopeU` | `0.0` each (required) | Bias-correction residual/slope terms. Always required (a hard parse error if any is missing, regardless of `remove_qu_bias`) — but the *values* only matter when `remove_qu_bias=y`; leave at `0.0` otherwise. |
| `path_I` | = `path` | Directory of the I-cube; required if `remove_qu_bias=y`. |
| `infileI` | — | Stokes-I input filename; required if `remove_qu_bias=y`; per-band comma list for multi-band. |

**RM sampling (`ofac`/`fac`/`use_auto_rm_range` required):**

| Key | Default | Meaning |
|---|---|---|
| `ofac` | `4` | Oversampling factor (≥1); `nrm_out = nrm * ofac`. |
| `fac` | `3.14159265358979` (π) | λ² conversion constant — leave as-is. |
| `use_auto_rm_range` | — | `0`: use manual `beg_rm`/`end_rm`/`nrm` grid; `1`: derive the RM range from the data (those three become optional overrides). |
| `beg_rm` | `-50.0` | RM grid start (rad/m²); required if `use_auto_rm_range=0`. |
| `end_rm` (alias `max_rm`) | `50.0` | RM grid end; must be > `beg_rm`; required if `use_auto_rm_range=0`. |
| `nrm` (alias `nrm_out`) | `100` | RM samples before oversampling (≥1); required if `use_auto_rm_range=0`. |

**Output format (optional):**

| Key | Default | Meaning |
|---|---|---|
| `output_mode` | `ap` | `ap`: write amplitude+phase; `ri`: write real+imaginary instead. |
| `ap_angle_mode` | `phase` | (`output_mode=ap` only) `phase`: PHA cube is `arg(F)`; `pol`: PHA cube is `0.5*arg(F)` (polarization angle), written as `.POLA.RMCUBE.FITS`. |
| `lsq_ref_mode` | `zero` | `zero\|mid\|centroid\|min\|max\|fixed` — the λ² phase reference the dirty cube is built at. `zero` preserves this project's historical convention. The value actually used is recorded in the output's `LSQREF` FITS header keyword. |
| `lsq_ref_fixed_value` | `0.0` | Required if `lsq_ref_mode=fixed`. |

**Masking & optional outputs:**

| Key | Default | Meaning |
|---|---|---|
| `mask_cube_file` | empty (= `<outfile>.MASK.CUBE.FITS`) | Override the mask-cube output path. |
| `mask_input_cube_file` | empty (= none) | Externally supplied FITS mask applied on top of NaN/Inf and bad-channel masking. |
| `mask_trust_mode` | `safe` | `safe`: tolerate minor mask/data shape mismatches; `strict`: reject on any mismatch. |
| `write_mask_output` | `y` | Write the per-channel MASK cube. |
| `write_nvalid_output` | `y` | Write the per-pixel NVALID map. |

**Cubestat / peak maps:**

| Key | Default | Meaning |
|---|---|---|
| `cubestat` | `n` | `y`: also write PEAK/RM_PEAK/ANG_PEAK/SNR 2-D maps. |

**GPU:**

| Key | Default | Meaning |
|---|---|---|
| `use_gpu` (alias `use_gpus`) | `n` | Request GPU offload; on a CPU-only binary, prints a warning and falls back to CPU. |
| `gpu_vram_mib` | `0` (auto-detect) | Override detected VRAM size (MiB) for sub-block planning. |
| `mem_frac_vram` | `0.70` | Fraction (0,0.95] of VRAM budgeted per compute sub-block. |

**Tile memory planning & I/O parallelism:**

| Key | Default | Meaning |
|---|---|---|
| `tile_auto` | `y` | Auto-size tiles from `mem_frac_ram` (recommended). |
| `tile_ra` / `tile_dec` | `0`/`0` (auto) | Manual tile-size override; ignored when `tile_auto=y`. |
| `mem_frac_ram` | `0.25` | Fraction (0,0.95] of system RAM budgeted per tile (see [Shared Parameters](#shared-parameters) above). |
| `io_read_threads` | `1` | N independent read-only FITS handles per input cube. |
| `nwriters` | `1` | N-way parallel AMP/PHA writes via raw stream I/O. |
| `io_overlap` | `n` | `y`: overlap tile N's write with tile N+1's read/compute on a background thread. |

**Logging & timing:**

| Key | Default | Meaning |
|---|---|---|
| `log_level` | `info` | `error\|warn\|info\|debug`. |
| `timing_enabled` | `n` | Print the Timing summary / Macro breakdown. |
| `timing_tile_enabled` | `n` | Include tile-level stage timers. |
| `timing_io_enabled` | `n` | Include I/O stage timers. |
| `log_output_file` | empty (stdout) | Non-empty: append to this file instead. |
| `timing_csv_file` | empty (no CSV) | Non-empty: append one benchmark CSV row per run. |

**Misc:**

| Key | Default | Meaning |
|---|---|---|
| `dry_run` | `n` | `y`: read headers, run the tile planner, write `tile_autotune.cfg`/`runtime_estimate.txt`, then exit — no pixel data read, no output written. |

Boolean keys accept `1`/`y`/`Y`/`t`/`T` (first non-blank character) as
true; anything else (including blank) is false.

### Output files

- `output_mode=ap` + `ap_angle_mode=phase` (default): `<outfile>.AMP.RMCUBE.FITS`, `<outfile>.PHA.RMCUBE.FITS`
- `output_mode=ap` + `ap_angle_mode=pol`: `<outfile>.AMP.RMCUBE.FITS`, `<outfile>.POLA.RMCUBE.FITS`
- `output_mode=ri`: `<outfile>.REAL.RMCUBE.FITS`, `<outfile>.IMAG.RMCUBE.FITS`
- Always (unless disabled): `<outfile>.MASK.CUBE.FITS`, `<outfile>.NVALID.MAP.FITS`
- If `cubestat=y`: `<outfile>.PEAK.MAP.FITS`, `<outfile>.RM_PEAK.MAP.FITS`, `<outfile>.ANG_PEAK.MAP.FITS`, `<outfile>.SNR.MAP.FITS`

**Beam metadata:** if the reference band's own Q cube has `CASAMBM=T`
(a per-channel CASA multi-beam cube, e.g. not yet run through
`convolve_cubes`), the AMP/PHA cubes and PEAK/RMPEAK/ANGPEAK/SNR maps
get that same `CASAMBM=T` plus the reference band's per-channel
`BEAMS` table — this is the one to trust; ignore the single
`BMAJ`/`BMIN`/`BPA` these outputs also carry, since with a varying
beam it's just the reference band's nominal value, not a common
resolution. Only the reference band's own table is ever attached this
way — a non-reference band's own channel-to-channel variation is only
flagged in the run log, never captured in the output `BEAMS` table. Otherwise (no `CASAMBM`), those same outputs simply
carry the input's single `BMAJ`/`BMIN`/`BPA` unchanged. (MASK/NVALID
are validity bookkeeping, not flux data, and always just get the plain
scalar.) A mismatch between bands' own beam metadata is a runtime
warning, not a hard error — RM synthesis itself doesn't depend on beam
metadata for correctness.

### Example

```bash
make OMP=1 GPU=0
bin/rm_synthesis_release_cpu_omp cfg/rmsynth.cfg
```

[cfg/rmsynth.cfg](../../cfg/rmsynth.cfg) is a full annotated template
with every key shown in context, grouped by section (required vs.
optional, with defaults). See also [cfg/CONFIG_README.md](../../cfg/CONFIG_README.md)
for the parser's own rules (required keys, validation, output-filename
conventions).

---

## 2. `reproject_cubes`

Reprojects two or more FITS cubes onto one common sky grid using
Starlink AST for WCS handling and `astResampleR` for resampling —
needed before a multi-band `rm_synthesis` run whose input bands don't
already share the same sky grid. Auto-detects the 2 sky axes regardless
of axis order/adjacency, and loops generically over every combination
of non-sky axes (channel, Stokes, etc). Processes the cube in OpenMP-
parallelised, memory-budgeted blocks of planes, with optional
background-thread write overlap.

### Invocation

```bash
bin/reproject_cubes mode=<intersection|union|reference> reffile=<reference_file> infiles=<file>[,<file>...]
    [mem_frac_ram=<fraction>] [io_overlap=y|n] [nwriters=<n>]
    [log_level=<level>] [timing_enabled=y|n] [log_output_file=<path>] [dry_run=y|n]
bin/reproject_cubes --config <cfgfile>
bin/reproject_cubes --config <cfgfile> mode=<...> [reffile=<...>] [infiles=<...>]
bin/reproject_cubes --help | -h
```

No positional arguments — every value is `key=value` (no spaces around
`=`), given on the command line and/or via `--config <cfgfile>` (same
`key=value` file syntax, `#`/`;` comments). If both are given, each CLI
`key=value` overrides only that same key from the config file.

### Parameters

| Key | Default | Required | Meaning |
|---|---|---|---|
| `mode` | — | yes | `intersection`: output grid shrinks to the overlap of all inputs with the reference. `union`: output grid grows to cover all inputs and the reference. `reference`: output grid is the reference file's own extent. Zero overlap between any input and the running grid is always a hard failure regardless of mode. |
| `reffile` | — | yes | Reference FITS cube whose WCS anchors the output grid. |
| `infiles` | — | yes | Comma-separated list of 1–50 input FITS cube paths. |
| `mem_frac_ram` | `0.25` | no | Fraction (0,0.95] of system RAM budgeted for one read/resample/write block of planes (see [Shared Parameters](#shared-parameters) above). Threads parallelize only *within* one block, never across blocks — too small a value both increases CFITSIO call count and discards most of the OpenMP speedup (a WARNING is printed if this looks likely). |
| `io_overlap` | `n` | no | `y`: write each block on a background thread, overlapped with the *next* block's read+resample. Only one background write is ever in flight. |
| `nwriters` | `1` | no | `>1`: split each block's own write across this many disjoint writer threads (clamped to `[1, OMP_NUM_THREADS]`) instead of one. Only useful alongside `io_overlap=y` on a disk that benefits from concurrent writes (e.g. NVMe) — see `dry_run` below for a per-machine suggestion. |
| `log_level` | `info` | no | `error\|warn\|info\|debug`. At `info` and above, a per-file stage-timing summary (seconds and % of that file's own total, per stage — e.g. `block_read`, `block_convolve`/`reproject_compute`, `block_write`) is always printed after each file finishes, regardless of `timing_enabled`. `debug` additionally logs one `tile_thread` line per block per worker thread (read/process/write start and end times) — fine-grained enough to drive `scripts/plot_tile_async_swimlane.py`'s per-thread timing plot, useful for seeing exactly where threads stall or serialize. |
| `timing_enabled` | `n` | no | Print a whole-run (all files combined) stage timing summary at exit, in addition to the always-on per-file summary described under `log_level` above. |
| `log_output_file` | empty (stdout) | no | Non-empty: append log/timing output to this file instead. |
| `dry_run` | `n` | no | `y`: check the first `infiles=` entry's own target disk (rotational vs. SSD/NVMe, via `/sys/block/.../queue/rotational`), print a suggested `io_overlap`/`nwriters` starting point for that disk type, write it to `reproject_cubes_dryrun.cfg`, then exit — no pixel data read, no output written. |

### Output files

One reprojected cube per input, named `<infile-without-extension>_REPROJ.FITS` (not configurable). A pre-existing file at that path aborts the run before anything is touched.

### Example

```bash
make reproject_cubes
bin/reproject_cubes mode=intersection reffile=ref.fits infiles=a.fits,b.fits
```

---

## 3. `convolve_cubes`

Convolves all channels, across all input files, to one common angular
resolution — needed whenever channels don't already share one
resolution, whether that's because of combining bands observed at
different resolutions, or because a single band's own native beam
varies channel-to-channel with frequency (the normal case for any cube
with channel-wise varying beams). Reads
per-channel beams from a CASA-style `BEAMS` binary table
(`beamfiles=auto`) or a portable ASCII/CSV beam log (see
`cfg/example_beamLog.txt`/`.csv`). A channel missing
from the beam log, or with `BMAJ`/`BMIN`=0, is written as an all-NaN
plane rather than convolved (and is automatically excluded by
`rm_synthesis`'s own NaN detection later). Each input must have exactly
2 sky axes plus one non-degenerate FREQ axis; run separate Stokes
slices as separate `infiles`.

### Invocation

```bash
bin/convolve_cubes infiles=<file1>[,<file2>...] [outsuffix=<suffix>]
    [beamfiles=<spec1>[,<spec2>...]] [badchan_file=<file1>[,<file2>...]]
    [target_bmaj=<arcsec> target_bmin=<arcsec> target_bpa=<deg>]
    [max_common_bmaj=<arcsec>] [mem_frac_ram=<fraction>]
    [npts=<n>] [khachiyan_tol=<tol>] [io_overlap=y|n] [nwriters=<n>]
    [log_level=<level>] [timing_enabled=y|n] [log_output_file=<path>] [dry_run=y|n]
bin/convolve_cubes --config <cfgfile>
bin/convolve_cubes --help | -h
```

Same `key=value`-only, no-positional-args convention as `reproject_cubes`.

### Parameters

| Key | Default | Required | Meaning |
|---|---|---|---|
| `infiles` | — | yes | 1–50 comma-separated FITS cube paths. |
| `outsuffix` | `_CONV.FITS` | no | Appended to each infile's own path (trailing `.fits`/`.FITS` stripped first). |
| `beamfiles` | `auto` for every input | no | Comma list, one entry per infile, in order. Each entry is either the literal word `auto` (read that infile's own CASA BEAMS table) or a path to an ASCII/CSV beam-log file (`channel bmaj_arcsec bmin_arcsec bpa_deg`, 1-indexed, `#`-comments allowed). If given, must list exactly as many entries as `infiles`. |
| `badchan_file` | none | no | One entry per `infiles` entry, same order, comma-separated — the list must have exactly as many entries as `infiles`, including a position for every file with nothing to list: give that position the literal value `none` — do not leave it out or leave it blank. E.g. with `infiles=a.fits,b.fits`, `badchan_file=none,file2.txt` gives `a.fits` no manual list and `b.fits` the second position (`file2.txt`); `badchan_file=,file2.txt` (blank first position) is a hard error once `badchan_file=` is given at all. Omitting `badchan_file=` entirely is unaffected — still valid, still means no list for any infile. Each file is one channel index per line, 1-indexed — same convention as `rm_synthesis`'s `global_badchan_file`. Independent of, and in addition to, automatic bad-channel detection from each infile's own CASA BEAMS table (a degenerate near-zero beam entry) — use this for channels known bad for reasons a beam table alone wouldn't capture (e.g. RFI). Either way, a bad channel is written as an all-NaN plane, not convolved. |
| `target_bmaj` / `target_bmin` / `target_bpa` | none (auto-derive) | no | Explicit target beam (arcsec/arcsec/deg) — give all three together to skip automatic common-beam derivation entirely. |
| `max_common_bmaj` | none | no | If the AUTO-derived common beam's BMAJ exceeds this (arcsec), refuse to proceed. Ignored when an explicit target is given. |
| `mem_frac_ram` | `0.25` | no | Fraction (0,0.95] of system RAM budgeted for one read/convolve/write block of planes (see [Shared Parameters](#shared-parameters) above). |
| `npts` | `2000` | no | Boundary points sampled per beam (≥12), passed to `commonbeam_mod`'s `find_common_beam`. |
| `khachiyan_tol` | `1.0e-5` | no | Khachiyan-algorithm convergence tolerance for the common-beam fit. |
| `io_overlap` | `n` | no | `y`: write each block on a background thread, overlapped with the *next* block's read+convolve. |
| `nwriters` | `1` | no | `>1`: split each block's own write across this many disjoint writer threads (clamped to `[1, OMP_NUM_THREADS]`) instead of one. Only useful alongside `io_overlap=y` on a disk that benefits from concurrent writes (e.g. NVMe) — see `dry_run` below for a per-machine suggestion. |
| `log_level` | `info` | no | `error\|warn\|info\|debug`. At `info` and above, a per-file stage-timing summary (seconds and % of that file's own total, per stage — e.g. `block_read`, `convolve_compute`, `block_write`) is always printed after each file finishes, regardless of `timing_enabled`. `debug` additionally logs one `tile_thread` line per block per worker thread (read/convolve/write start and end times) — fine-grained enough to drive `scripts/plot_tile_async_swimlane.py`'s per-thread timing plot, useful for seeing exactly where threads stall or serialize. |
| `timing_enabled` | `n` | no | Print a whole-run (all files combined) stage timing summary at exit, in addition to the always-on per-file summary described under `log_level` above. |
| `log_output_file` | empty (stdout) | no | Non-empty: append log/timing output to this file instead. |
| `dry_run` | `n` | no | `y`: check the first `infiles=` entry's own target disk (rotational vs. SSD/NVMe, via `/sys/block/.../queue/rotational`), print a suggested `io_overlap`/`nwriters` starting point for that disk type, write it to `convolve_cubes_dryrun.cfg`, then exit — no pixel data read, no output written. |

### Output files

One convolved cube per input, `<infile-without-extension><outsuffix>` (default `_CONV.FITS`), carrying `CASAMBM=T` plus a freshly synthesized `BEAMS` table (one row per channel — the common target beam for every channel actually convolved, and the same degenerate sentinel CASA itself uses for a bad/skipped channel).

### Example

```bash
make convolve_cubes
bin/convolve_cubes infiles=bandA.fits,bandB.fits mem_frac_ram=0.25
```

---

## 4. `match_cubes`

Consolidates `reproject_cubes` and `convolve_cubes` into one tool that
can run either stage alone, or chain both **through memory** with no
intermediate FITS file written — a real saving for 200GB+ cubes.
Default chain order is convolve-then-reproject (low-pass-filtering
before resampling avoids baking interpolation/aliasing error into a
marginally-sampled native PSF, and is usually cheaper too).
`stages=convolve`/`both` inherit `convolve_cubes`'s stricter axis scope
(exactly 2 sky axes + 1 FREQ axis); `stages=reproject` alone keeps the
fully general N-dimensional handling. Includes a "skip-if-already-
matched" check (tight tolerances only, never a fuzzy match) and an
optional machine-readable manifest of which files were skipped vs.
processed.

### Invocation

```bash
bin/match_cubes stages=reproject|convolve|both
    [order=convolve_reproject|reproject_convolve]
    infiles=<file1>[,<file2>...]
    [footprint_mode=intersection|union|reference] [reffile=<file>]
    [beamfiles=<spec1>[,<spec2>...]] [badchan_file=<file1>[,<file2>...]]
    [target_bmaj=<arcsec> target_bmin=<arcsec> target_bpa=<deg>]
    [max_common_bmaj=<arcsec>] [mem_frac_ram=<fraction>] [outsuffix=<suffix>]
    [npts=<n>] [khachiyan_tol=<tol>] [manifest=<path>] [io_overlap=y|n] [nwriters=<n>]
    [log_level=<level>] [timing_enabled=y|n] [log_output_file=<path>] [dry_run=y|n]
bin/match_cubes --config <cfgfile>
bin/match_cubes --help | -h
```

Same `key=value`-only convention as the other two.

### Parameters

| Key | Default | Required | Meaning |
|---|---|---|---|
| `stages` | — | yes | `reproject \| convolve \| both`. |
| `order` | `convolve_reproject` | no | `convolve_reproject \| reproject_convolve`; only meaningful when `stages=both`. |
| `infiles` | — | yes | 1–50 comma-separated FITS cube paths. |
| `footprint_mode` | — | if `stages` includes `reproject` | `intersection \| union \| reference` — same semantics as `reproject_cubes`'s `mode`. |
| `reffile` | — | if `stages` includes `reproject` | Reference cube path. |
| `beamfiles` | `auto` per input | no | Used when `stages` includes `convolve` — same semantics as `convolve_cubes`. |
| `badchan_file` | none | no | Used when `stages` includes `convolve` — same semantics as `convolve_cubes`. |
| `target_bmaj` / `target_bmin` / `target_bpa` | none | no | Used when `stages` includes `convolve` — same semantics as `convolve_cubes`. |
| `max_common_bmaj` | none | no | Used when `stages` includes `convolve` — same semantics as `convolve_cubes`. |
| `mem_frac_ram` | `0.25` | no | Fraction (0,0.95] of RAM budgeted for one read/process/write block (see [Shared Parameters](#shared-parameters) above). |
| `outsuffix` | `_REPROJ.FITS` / `_CONV.FITS` / `_MATCHED.FITS` | no | Default depends on `stages` if not given explicitly. |
| `npts` | `2000` | no | Used when `stages` includes `convolve` — same as `convolve_cubes`. |
| `khachiyan_tol` | `1.0e-5` | no | Used when `stages` includes `convolve` — same as `convolve_cubes`. |
| `manifest` | none | no | Path to write a tab-separated `<infile> SKIPPED\|PROCESSED <effective_path>` record, one line per input. Aborts if the manifest path already exists. |
| `io_overlap` | `n` | no | `y`: write each block on a background thread, overlapped with the *next* block's read+process. |
| `nwriters` | `1` | no | `>1`: split each block's own write across this many disjoint writer threads (clamped to `[1, OMP_NUM_THREADS]`) instead of one. Only useful alongside `io_overlap=y` on a disk that benefits from concurrent writes (e.g. NVMe) — see `dry_run` below for a per-machine suggestion. |
| `log_level` | `info` | no | `error\|warn\|info\|debug`. At `info` and above, a per-file stage-timing summary (seconds and % of that file's own total, per stage — e.g. `block_read`, `convolve_compute`, `reproject_compute`, `block_write`) is always printed after each file finishes, regardless of `timing_enabled`. `debug` additionally logs one `tile_thread` line per block per worker thread (read/process/write start and end times) — fine-grained enough to drive `scripts/plot_tile_async_swimlane.py`'s per-thread timing plot, useful for seeing exactly where threads stall or serialize. |
| `timing_enabled` | `n` | no | Print a whole-run (all files combined) stage timing summary at exit, in addition to the always-on per-file summary described under `log_level` above. |
| `log_output_file` | empty (stdout) | no | Non-empty: append log/timing output to this file instead. |
| `dry_run` | `n` | no | `y`: check the first `infiles=` entry's own target disk (rotational vs. SSD/NVMe, via `/sys/block/.../queue/rotational`), print a suggested `io_overlap`/`nwriters` starting point for that disk type, write it to `match_cubes_dryrun.cfg`, then exit — no pixel data read, no output written. |

### Reference-band advisory

When `stages` includes `reproject` and `infiles=` spans more than one
distinct sky pointing/pixel scale, `match_cubes` checks whether
`reffile=` was chosen from the group of inputs with the most total
channels, and prints an `ADVISORY:` message if not. Matching a file
onto its own kind of pointing (the reference's own group) is fast;
matching it onto a different one is much slower — so the
group with the most channels should generally be the one that defines
the output grid, since that's the group that gets the cheap treatment.
This is advice only: `reffile=` is never changed automatically, since
it also determines the output grid extent, a choice this program
leaves to you.

### Output files

One output cube per input (unless skipped — see below), `<infile-without-extension><outsuffix>`. `manifest=<path>` additionally writes the skip/process record described above.

### Example

```bash
make match_cubes
bin/match_cubes stages=both order=convolve_reproject \
  footprint_mode=intersection reffile=ref.fits infiles=bandA.fits,bandB.fits \
  mem_frac_ram=0.25
```

---

## 5. `rmclean_cubes`

Standalone RM-CLEAN (Högbom-style deconvolution), driving `rmclean_mod`
(`src/rmclean.f90`, pure computation, no FITS I/O of its own) against a
real *dirty* AMP/PHA cube pair that `rm_synthesis` itself wrote, plus
that run's own `.MASK.CUBE.FITS`/`CHANFREQ` table (which records true
per-channel frequency, needed to build each pixel's own RMSF from its
own valid-channel subset). Writes 6 output cubes
(`CLEAN`/`RESID`/`RESTORED` × `AMP`/`PHA`). Refuses to run (Gate 0) if
the cube's existing RM grid doesn't resolve the RMSF FWHM at
`min_samples_per_fwhm` samples — this tool cannot resample the RM axis,
only validate it. CLEAN uses a tiered peak-refinement strategy (a cheap
matched-filter fit, escalating to a full local search only when misfit
exceeds `refine_nsigma`×noise), parallelised one OpenMP thread per
pixel, with a mask-pattern cache sharing one RMSF table across every
pixel with the same valid-channel pattern. Fully memory-budgeted and
tileable exactly like `rm_synthesis` — both the AMP/PHA float cubes
*and* the mask cube are read in RAM-budgeted tiles (not held whole in
memory), so a cube far larger than available RAM CLEANs in bounded
memory regardless of machine size. GPU support is not yet implemented
for this tool.

> **If you're migrating an old cfg**: `threshold=`, `threshold_snr=`,
> `noise_percentile=`, `noise_nlos=`, and `noise_seed=` have all been
> **removed** (not aliased) — a cfg using any of them will now fail to
> parse. Replace with `abs_flux_floor=`/`auto_nsigma=` below.

### Invocation

```bash
bin/rmclean_cubes ampfile=<f> phafile=<f> maskfile=<f> outfile=<base>
    [abs_flux_floor=<v>] [auto_nsigma=<n>]
    [niter=<n>] [gain=<g>] [min_samples_per_fwhm=<f>] [refine_nsigma=<f>] [table_oversample=<n>] [restore_fwhm=<v>]
    [trace_ix=<n>] [trace_iy=<n>] [log_every=<n>]
    [lsq_ref_report_mode=intrinsic|centroid|min|max|mid|fixed] [lsq_ref_report_value=<v>]
    [lsq_ref_compute_mode=native|zero|centroid|min|max|fixed] [lsq_ref_compute_value=<v>]
    [mask_pattern_cache_max=<n>] [cache_eviction_policy=belady|hitcount] [min_valid_chan_frac=<f>]
    [mem_frac_ram=<f>] [tile_ra=<n>] [tile_dec=<n>] [tile_auto=y|n]
    [io_read_threads=<n>] [nwriters=<n>] [io_overlap=y|n]
    [log_level=<level>] [timing_enabled=y|n] [log_output_file=<path>]
    [write_clean_diagnostics=y|n]
bin/rmclean_cubes --config <cfgfile>
bin/rmclean_cubes --help | -h
```

Same `key=value`-only convention as the other 3 modernized tools.

### Parameters

**Required (checked as a group — missing any one is a hard error):**

| Key | Meaning |
|---|---|
| `ampfile` | The dirty `AMP.RMCUBE.FITS` that `rm_synthesis` wrote. |
| `phafile` | The matching dirty `PHA.RMCUBE.FITS`. |
| `maskfile` | The matching `MASK.CUBE.FITS` (with its own `CHANFREQ` table). |
| `outfile` | Base name for the 6 output cubes. |

**CLEAN stopping criteria** — three independent, unambiguously-named
conditions, checked every iteration in this order. `niter` is always
the hard backstop; the other two are each independently optional and
**may be combined** (whichever fires first wins for that pixel);
neither given is valid too — `niter` alone then governs, and a NOTE is
printed saying so.

| Key | Default | Meaning |
|---|---|---|
| `niter` | `500` | Hard maximum Hogbom CLEAN iterations. |
| `gain` | `0.1` | Hogbom CLEAN loop gain. |
| `abs_flux_floor` | none (opt-in) | Stop the instant a pixel's own peak amplitude drops to/below this literal fixed flux value. Accepts a bare number (native AMP-cube units) or a value suffixed with `Jy`/`mJy`/`uJy`, **no space** (e.g. `abs_flux_floor=10mJy`), converted against the AMP cube's own `BUNIT`. No noise/baseline adjustment — a pure "brightest remaining feature below X" comparison. |
| `auto_nsigma` | none (opt-in) | Stop the instant a pixel's own peak amplitude drops to/below `auto_nsigma` × that *same pixel's own* noise sigma. The sigma is estimated **once per pixel**, from that pixel's own full dirty amplitude spectrum's interquartile range, converted to sigma via the Rayleigh-distribution analytic relation (dirty amplitude is Rayleigh-, not Gaussian-, distributed). |

**RM-resolution / peak-refinement tuning:**

These four rarely need touching — Gate 0 and the tiered peak-refinement
strategy work well at their defaults. They exist for tuning, not
day-to-day use.

| Key | Default | Meaning |
|---|---|---|
| `min_samples_per_fwhm` | `2.0` | Gate 0's RM-resolution check: `\|CDELT3\| ≤ fwhm/min_samples_per_fwhm`. `\|CDELT3\|` (the RM grid spacing) is already fixed by whatever `rm_synthesis` run produced the input cube — this tool never resamples the RM axis, it only validates it, and `stop`s immediately if the check fails. `min_samples_per_fwhm` doesn't touch the grid; it only sets how strict that check is. Raising it demands a finer grid (more likely to fail on a given cube); lowering it accepts a coarser one (more likely to pass), down to a hard floor of `1.0` (the parser rejects anything less — below 1 sample per FWHM the grid can't resolve the RMSF's main lobe at all). So if Gate 0 fails at the default `2.0`, the fix is either to re-run `rm_synthesis` with a finer `cdelt3`, or to explicitly set `min_samples_per_fwhm=` lower here (as low as `1.0`) to accept the grid you already have. |
| `refine_nsigma` | `3.0` | Controls the two-tier peak-refinement strategy CLEAN runs every iteration: a cheap fixed-location matched-filter fit runs first, and is accepted as-is if its leftover misfit is within `refine_nsigma`× that iteration's noise estimate; only when the misfit exceeds this threshold does the slower full local search run instead. Raising it accepts the cheap fit more often (faster, slightly less precise peak location); lowering it escalates to the full search more often (slower, more precise). |
| `table_oversample` | `20` | Fineness of the per-pattern RMSF lookup table (see *Mask-pattern cache* below). Grid spacing = native `\|CDELT3\|` / `table_oversample`. Higher values interpolate the RMSF more accurately between grid points, at the cost of a larger one-time table build. This is a floor, not a fixed value: `rmclean_cubes` raises it automatically for a run (printing a `NOTE`) if the chosen `lsq_ref_compute` would otherwise undersample the table, and never lowers it below what you set. Raise it manually for extra interpolation accuracy on high-S/N data. |
| `restore_fwhm` | derived from the data | The restoring beam FWHM (rad/m²) used to convolve CLEAN components into the RESTORED cube. By default this is computed automatically from the RMSF itself (`compute_rmsf_fwhm_multiband`); set this to force a specific value instead, e.g. to match a restoring beam used elsewhere for comparison. |

**Debug/tracing:**

| Key | Default | Meaning |
|---|---|---|
| `trace_ix` / `trace_iy` | `0` (disabled) | 1-indexed GLOBAL pixel coordinates to log a per-iteration CLEAN trend for (peak_val/rms_val/cumulative cleaned flux, exact stop_reason). One pixel only — a debugging aid, not for production runs. |
| `log_every` | `50` | Throttles the trace log lines to every Nth iteration; iteration 1 and the final iteration are always logged regardless. |

**Phase-reference (`lsq_ref`) controls:**

| Key | Default | Meaning |
|---|---|---|
| `lsq_ref_report_mode` | `intrinsic` (=0.0) | `intrinsic\|centroid\|min\|max\|mid\|fixed` — where to report the derotated chi0/restored phase (a safe, independent post-processing choice). |
| `lsq_ref_report_value` | — | Required only if `lsq_ref_report_mode=fixed`. |
| `lsq_ref_compute_mode` | `mid` (recommended) | `native\|zero\|centroid\|min\|max\|fixed` — the reference this program's own RMSF table/CLEAN computation uses internally. `mid` minimizes per-pixel escalation-search cost and does not change the RM grid size. |
| `lsq_ref_compute_value` | — | Used only with `lsq_ref_compute_mode=fixed`. |

> The `LSQREF` keyword on every CLEAN/RESID/RESTORED output records
> `lsq_ref_report` — the reference the phase actually on disk is
> derotated to — not `lsq_ref_compute` (an internal-only value that
> never appears in any output) and not the input dirty cube's own
> `LSQREF` (`lsq_ref_native`, which can differ from `lsq_ref_report`
> whenever either tool is set away from its own default).

**Mask-pattern cache:**

> The RMSF table (`table_oversample` above) is a single, fixed,
> unit-strength curve — the RMSF's own value at each of many
> finely-spaced offsets — built once per distinct valid-channel pattern
> and never shifted, rescaled, or otherwise modified for the rest of
> the run. Peak location and amplitude are always found beforehand (see
> *RM-resolution / peak-refinement tuning* above) and never come from
> this table. Each CLEAN iteration, for each of the `nrm` native output
> points independently, what moves is the *query* — that point's own
> distance to the found peak — not the table: the two table entries
> bracketing that distance are
> blended, rotated by the peak's own phase, scaled by its own measured
> flux, and subtracted, one output point at a time. Building the whole
> per-point beam this way, `nrm` separate small lookups against one
> unmoving table, is what avoids ever having to shift a discretely
> sampled RMSF onto a fractional peak location as a single operation.

`belady` is L. A. Belady's offline-optimal page-replacement algorithm
(also called OPT/MIN) — L. A. Belady, "A Study of Replacement
Algorithms for a Virtual-Storage Computer," *IBM Systems Journal*,
vol. 5, no. 2, pp. 78–101, 1966. The admission-control refinement and
the whole-mask-cube pre-scan that supplies its lookahead are this
project's own addition on top of that algorithm.

| Key | Default | Meaning |
|---|---|---|
| `mask_pattern_cache_max` | `4096` | Pixels sharing the same valid-channel mask pattern share one RMSF table: it's built once, on that pattern's first encounter during the real per-tile CLEAN pass (not during the lookahead pre-scan below, which only records which patterns exist and where — it never builds a real table itself), and reused by every later pixel with the same pattern. Past this many distinct patterns resident at once, additional new patterns fall back to a one-off table built and discarded per pixel (a performance safety valve, not a correctness issue). |
| `cache_eviction_policy` | `belady` | `belady\|hitcount` — which cached pattern to evict once `mask_pattern_cache_max` is reached, rather than simply refusing new entries. `belady`: offline-optimal, informed by a one-time whole-mask-cube pre-scan (Pass 0) that evicts whichever cached pattern isn't needed again for the longest time (or ever) — no tuning parameter, provably the best possible for a given cache size. `hitcount`: a plain least-frequently-used fallback needing no pre-scan — cheaper to start, but not competitive with `belady`; an opt-out, not the recommended choice. Pass 0 has its own safety valve: if distinct patterns exceed 10% of total image pixels (too little reuse for a lookahead cache to be worth its cost), it aborts with a warning and falls back to `hitcount` for that run. |
| `min_valid_chan_frac` | `0.0` (off) | Pixels whose own valid-channel fraction (summed across ALL bands, out of the full channel count) falls below this are never CLEANed — output NaN across every RM-bin instead. `0.0` means every pixel with at least 1 valid channel is CLEANed (original behaviour); e.g. `0.7` skips any pixel with less than 70% of channels valid. Sparse-coverage pixels have a coarser/noisier RMSF regardless of compute spent on them, so skipping them saves real time — the skipped pixels are also never registered in Pass 0's pattern registry or the runtime cache. |

When `cache_eviction_policy=belady`, Pass 0 finishes by printing a
**build-time advisory** (to both stdout and `log_output_file` if set):
a per-block, per-total estimate of table-*generation* time only (never
CLEAN time, which depends on the sky and can't be predicted from the
mask alone), derived from a full Belady-cache simulation plus a small
per-block timing sample — not just a naive "first occurrence" count.
Informational only; it does not affect the run. It also writes
`<outfile>.advisory.csv` (one row per block) — fine-grained enough to
drive `scripts/plot_rmclean_advisory.py`'s predicted-vs-actual plot,
useful for checking a run's real progress against what Pass 0 expected
(`--csv` is auto-detected from a `log_level=debug` run log via
`--actual-log` alone, no need to locate the CSV by hand).

**Memory/tiling (identical scheme/defaults to `rm_synthesis`, now covering the mask cube too):**

| Key | Default | Meaning |
|---|---|---|
| `mem_frac_ram` | `0.25` | Fraction (0,0.95] of system RAM budgeted for one tile's own read+compute+write working set — now including the mask's own footprint (see [Shared Parameters](#shared-parameters) above). |
| `tile_ra` / `tile_dec` | `0`/`0` (auto) | Manual tile-size override in pixels; ignored unless `tile_auto=n`. |
| `tile_auto` | `y` | `y`: derive `tile_ra`/`tile_dec` from `mem_frac_ram`; `n`: use the manual override (still clamped to image size). |

**I/O parallelism (identical scheme to `rm_synthesis`):**

| Key | Default | Meaning |
|---|---|---|
| `io_read_threads` | `1` | Parallel chunked reads of each tile's own AMP/PHA/mask depth range. |
| `nwriters` | `1` | Parallel chunked writes of each tile's 6 output cubes; `>1` bypasses CFITSIO for pixel writes (raw stream writes at computed byte offsets). |
| `io_overlap` | `n` | `y`: write each tile's 6 output cubes on a background thread while the next tile's read/compute proceeds. |

**Logging & timing:**

| Key | Default | Meaning |
|---|---|---|
| `log_level` | `info` | `error\|warn\|info\|debug`. Any unrecognized value silently falls back to `info` rather than erroring — true of every tool in this package, not specific to `rmclean_cubes`. |
| `timing_enabled` | `n` | Print a stage timing summary. |
| `log_output_file` | empty (stdout) | Non-empty: append to this file instead. |

**Per-pixel diagnostic maps:**

| Key | Default | Meaning |
|---|---|---|
| `write_clean_diagnostics` | `y` | Write 7 additional per-pixel 2D diagnostic maps at run end (see Output files below). Negligible extra compute (all 7 are computed from values CLEAN already produces); the memory cost is a small fixed amount for the whole image regardless of tile size (~21 bytes/pixel), already accounted for in the `mem_frac_ram` tile-size budget. Set `n` to skip them. |

### Output files

`<outfile>.CLEAN.AMP.RMCUBE.FITS`, `<outfile>.CLEAN.PHA.RMCUBE.FITS`,
`<outfile>.RESID.AMP.RMCUBE.FITS`, `<outfile>.RESID.PHA.RMCUBE.FITS`,
`<outfile>.RESTORED.AMP.RMCUBE.FITS`, `<outfile>.RESTORED.PHA.RMCUBE.FITS`.

If `write_clean_diagnostics=y` (the default), also 7 per-pixel 2D maps:
`<outfile>.NITER.MAP.FITS`, `<outfile>.STOP_REASON.MAP.FITS`,
`<outfile>.RESID_PEAK.MAP.FITS`, `<outfile>.RESID_RMS.MAP.FITS`,
`<outfile>.N_COMPONENTS.MAP.FITS`, `<outfile>.COMP_RM_SPREAD.MAP.FITS`,
`<outfile>.TOTAL_POL_FLUX.MAP.FITS`.

### Example

```bash
make rmclean_cubes
bin/rmclean_cubes ampfile=out.AMP.RMCUBE.FITS phafile=out.PHA.RMCUBE.FITS \
  maskfile=out.MASK.CUBE.FITS outfile=out_cleaned \
  abs_flux_floor=20uJy auto_nsigma=1.0 niter=500 gain=0.1
```

Or via a config file — see `cfg/rmclean-example.cfg` for a fully
annotated template covering every key above with worked-through
reasoning for each default.
