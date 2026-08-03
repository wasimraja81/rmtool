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

## 1. `rm_synthesis`

The main RM-synthesis tool. Given Stokes Q and U spectral-line FITS
cubes (optionally several bands at once — "multi-band tomography"), it
performs pixel-by-pixel Faraday rotation-measure synthesis (a
Fourier-like transform from λ² to Faraday depth φ), producing RM cubes
of linear polarized intensity and polarization angle (or real/
imaginary) as a function of RA, Dec, RM. It supports bad-channel
removal, Q/U mean-subtraction and I-cube bias correction, subimage/
sub-channel extraction, automatic or manual RM-grid selection, mask/
peak-map diagnostics, GPU offload, and a memory-budgeted, multithreaded,
tileable I/O pipeline.

### Invocation

```bash
bin/rm_synthesis <cfgfile> [<addreq>]
```

Purely positional — **not** `key=value` on the command line. `<cfgfile>`
is a required path to a `KEY=VALUE` text config file (one key per line,
`#` for comments; strict parser — duplicate keys, unknown keys, and
unparsable values are all hard errors). `<addreq>` is accepted but not
currently used by anything downstream; omit it. There is no `--help`
flag for this tool — running it with no arguments prints a two-line
usage pointer, not a full option list (that's what this document and
the README's own "Configuration" section are for).

### Config keys

**Input / output paths (required):**

| Key | Default | Meaning |
|---|---|---|
| `path` | — | Directory containing the input Q/U cubes. |
| `infileQ` | — | Stokes-Q input filename, relative to `path`. Comma-separated for multi-band (e.g. `bandA.fits,bandB.fits`); band count is derived from this list's length. |
| `infileU` | — | Stokes-U input filename; must have the same band count as `infileQ`. |
| `reference_band` | `1` | Which band (1-indexed) populates single-band-era scalar fields (see below) when running multi-band. |
| `outfile` | — | Output basename; cube-type suffixes are appended automatically. |

**Bad-channel handling (required):**

| Key | Default | Meaning |
|---|---|---|
| `remove_badchan` | — | `y`/`n`: drop channels listed in `badchan_file`. |
| `badchan_file` (alias `global_badchan_file`) | — | One channel index per line. Must be present even when `remove_badchan=n` (the file need not exist in that case). Comma-separated per-band list for multi-band. |

**Subimage extraction (`subim` required; the rest optional, default = full cube):**

| Key | Default | Meaning |
|---|---|---|
| `subim` | — | `y`: only process the ranges below; `n`: full cube (ranges ignored). |
| `subim_ra_blc` / `subim_ra_trc` / `subim_ra_inc` | `1` / `0` (=max) / `1` | RA first pixel (≥1) / last pixel (0=max) / step (≥1). |
| `subim_dec_blc` / `subim_dec_trc` / `subim_dec_inc` | `1` / `0` (=max) / `1` | Same, for Dec. |
| `subim_chan_blc` / `subim_chan_trc` / `subim_chan_inc` | `0` (=first) / `0` (=max) / `1` | Same, for channel; per-band comma list allowed. |

**Q/U processing & bias correction:**

| Key | Default | Meaning |
|---|---|---|
| `rem_mean` | — (required) | `0`/`1`: subtract per-pixel mean from Q/U before synthesis. |
| `remove_qu_bias` | — (required) | `y`/`n`: apply I-cube-based bias correction (needs `path_I`/`infileI`). |
| `resiQ` / `slopeQ` / `resiU` / `slopeU` | `0.0` each | Bias-correction residual/slope terms; only meaningful when `remove_qu_bias=y`. |
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
| `mem_frac_ram` | `0.25` | Fraction (0,0.95] of system RAM budgeted per tile. |
| `io_read_threads` | `1` | N independent read-only FITS handles per input cube. |
| `io_write_threads` | `1` | N-way parallel AMP/PHA writes via raw stream I/O. |
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

### Example

```bash
make OMP=1 GPU=0
bin/rm_synthesis_release_cpu_omp cfg/rmsynth.cfg
```

See `README.md`'s "Configuration" section for a full annotated example
cfg file with every key shown in context.

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
| `mem_frac_ram` | `0.25` | no | Fraction (0,0.95] of system RAM budgeted for one read/resample/write block of planes. Threads parallelize only *within* one block, never across blocks — too small a value both increases CFITSIO call count and discards most of the OpenMP speedup (a WARNING is printed if this looks likely). |
| `io_overlap` | `n` | no | `y`: write each block on a background thread, overlapped with the *next* block's read+resample. Only one background write is ever in flight. |
| `log_level` | `info` | no | `error\|warn\|info\|debug`. |
| `timing_enabled` | `n` | no | Print a stage timing summary. |
| `log_output_file` | empty (stdout) | no | Non-empty: append log/timing output to this file instead. |

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
resolution — needed before a multi-band `rm_synthesis` run whose bands
were observed at different resolutions. Reads per-channel beams from a
CASA-style `BEAMS` binary table (`beamfiles=auto`) or a portable ASCII/
CSV beam log (see `cfg/example_beamLog.txt`/`.csv`). A channel missing
from the beam log, or with `BMAJ`/`BMIN`=0, is written as an all-NaN
plane rather than convolved (and is automatically excluded by
`rm_synthesis`'s own NaN detection later). Each input must have exactly
2 sky axes plus one non-degenerate FREQ axis; run separate Stokes
slices as separate `infiles`.

### Invocation

```bash
bin/convolve_cubes infiles=<file1>[,<file2>...] [outsuffix=<suffix>]
    [beamfiles=<spec1>[,<spec2>...]] [badchan_file=<file>]
    [target_bmaj=<arcsec> target_bmin=<arcsec> target_bpa=<deg>]
    [max_common_bmaj=<arcsec>] [mem_frac_ram=<fraction>]
    [npts=<n>] [khachiyan_tol=<tol>] [io_overlap=y|n]
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
| `badchan_file` | none | no | One channel index per line, 1-indexed — same convention as `rm_synthesis`'s `global_badchan_file`. |
| `target_bmaj` / `target_bmin` / `target_bpa` | none (auto-derive) | no | Explicit target beam (arcsec/arcsec/deg) — give all three together to skip automatic common-beam derivation entirely. |
| `max_common_bmaj` | none | no | If the AUTO-derived common beam's BMAJ exceeds this (arcsec), refuse to proceed. Ignored when an explicit target is given. |
| `mem_frac_ram` | `0.25` | no | Fraction (0,0.95] of system RAM budgeted for one read/convolve/write block of planes. |
| `npts` | `2000` | no | Boundary points sampled per beam (≥12), passed to `commonbeam_mod`'s `find_common_beam`. |
| `khachiyan_tol` | `1.0e-5` | no | Khachiyan-algorithm convergence tolerance for the common-beam fit. |
| `io_overlap` | `n` | no | `y`: write each block on a background thread, overlapped with the *next* block's read+convolve. |
| `log_level` | `info` | no | `error\|warn\|info\|debug`. |
| `timing_enabled` | `n` | no | Print a stage timing summary. |
| `log_output_file` | empty (stdout) | no | Non-empty: append log/timing output to this file instead. |

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
    [beamfiles=<spec1>[,<spec2>...]] [badchan_file=<file>]
    [target_bmaj=<arcsec> target_bmin=<arcsec> target_bpa=<deg>]
    [max_common_bmaj=<arcsec>] [mem_frac_ram=<fraction>] [outsuffix=<suffix>]
    [npts=<n>] [khachiyan_tol=<tol>] [manifest=<path>] [io_overlap=y|n]
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
| `mem_frac_ram` | `0.25` | no | Fraction (0,0.95] of RAM budgeted for one read/process/write block. |
| `outsuffix` | `_REPROJ.FITS` / `_CONV.FITS` / `_MATCHED.FITS` | no | Default depends on `stages` if not given explicitly. |
| `npts` | `2000` | no | Used when `stages` includes `convolve` — same as `convolve_cubes`. |
| `khachiyan_tol` | `1.0e-5` | no | Used when `stages` includes `convolve` — same as `convolve_cubes`. |
| `manifest` | none | no | Path to write a tab-separated `<infile> SKIPPED\|PROCESSED <effective_path>` record, one line per input. Aborts if the manifest path already exists. |
| `io_overlap` | `n` | no | `y`: write each block on a background thread, overlapped with the *next* block's read+process. |
| `log_level` | `info` | no | `error\|warn\|info\|debug`. |
| `timing_enabled` | `n` | no | Print a stage timing summary. |
| `log_output_file` | empty (stdout) | no | Non-empty: append log/timing output to this file instead. |

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
    [mask_pattern_cache_max=<n>] [mem_frac_ram=<f>] [tile_ra=<n>] [tile_dec=<n>] [tile_auto=y|n]
    [io_read_threads=<n>] [io_write_threads=<n>] [io_overlap=y|n]
    [log_level=<level>] [timing_enabled=y|n] [log_output_file=<path>]
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

| Key | Default | Meaning |
|---|---|---|
| `min_samples_per_fwhm` | `2.0` (floor `1.0`) | Gate 0's RM-resolution criterion: the cube's `\|CDELT3\|` must be ≤ `fwhm/min_samples_per_fwhm`. |
| `refine_nsigma` | `3.0` | Escalation threshold for tiered peak refinement — the cheap fixed-location fit is accepted when its leftover misfit is within `refine_nsigma`×the per-iteration noise estimate; beyond that a full local search runs. |
| `table_oversample` | `20` | Interpolation-table fineness for the RMSF offset table. |
| `restore_fwhm` | derived from the data | Override the restoring beam FWHM (rad/m²). |

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

**Mask-pattern cache:**

| Key | Default | Meaning |
|---|---|---|
| `mask_pattern_cache_max` | `4096` | Pixels sharing the same valid-channel mask pattern share one RMSF table, built once during an incremental per-tile pre-scan; past this many distinct patterns, additional patterns fall back to a one-off table per pixel (a performance safety valve, not a correctness issue). |

**Memory/tiling (identical scheme/defaults to `rm_synthesis`, now covering the mask cube too — see `docs/user/PARALLELISM.md`):**

| Key | Default | Meaning |
|---|---|---|
| `mem_frac_ram` | `0.25` | Fraction (0,0.95] of system RAM budgeted for one tile's own read+compute+write working set — now including the mask's own footprint. |
| `tile_ra` / `tile_dec` | `0`/`0` (auto) | Manual tile-size override in pixels; ignored unless `tile_auto=n`. |
| `tile_auto` | `y` | `y`: derive `tile_ra`/`tile_dec` from `mem_frac_ram`; `n`: use the manual override (still clamped to image size). |

**I/O parallelism (identical scheme to `rm_synthesis`):**

| Key | Default | Meaning |
|---|---|---|
| `io_read_threads` | `1` | Parallel chunked reads of each tile's own AMP/PHA/mask depth range. |
| `io_write_threads` | `1` | Parallel chunked writes of each tile's 6 output cubes; `>1` bypasses CFITSIO for pixel writes (raw stream writes at computed byte offsets). |
| `io_overlap` | `n` | `y`: write each tile's 6 output cubes on a background thread while the next tile's read/compute proceeds. |

**Logging & timing:**

| Key | Default | Meaning |
|---|---|---|
| `log_level` | `info` | `error\|warn\|info\|debug`. Any unrecognized value silently falls back to `info` rather than erroring — true of every tool in this package, not specific to `rmclean_cubes`. |
| `timing_enabled` | `n` | Print a stage timing summary. |
| `log_output_file` | empty (stdout) | Non-empty: append to this file instead. |

### Output files

`<outfile>.CLEAN.AMP.RMCUBE.FITS`, `<outfile>.CLEAN.PHA.RMCUBE.FITS`,
`<outfile>.RESID.AMP.RMCUBE.FITS`, `<outfile>.RESID.PHA.RMCUBE.FITS`,
`<outfile>.RESTORED.AMP.RMCUBE.FITS`, `<outfile>.RESTORED.PHA.RMCUBE.FITS`.

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
