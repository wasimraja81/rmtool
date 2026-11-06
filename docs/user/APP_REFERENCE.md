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
prints a summary of the same parameters. The two are written and
maintained separately, so if they ever disagree, what the tool does is
what counts; please report the difference.

Every invocation below is shown as run from `$RMTOOL_HOME` -- the
directory you cloned this repository into (see
[QUICKSTART.md](../../QUICKSTART.md#where-things-go) for the full
convention).

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

**`badchan_file`** (`rm_synthesis`, `convolve_cubes`, `match_cubes`): one
manually-specified bad channel index per line, 1-indexed — these
channels are rejected the same way in all three tools. A multi-input
list is comma-separated, one entry per input (band, for `rm_synthesis`;
file, for `convolve_cubes`/`match_cubes`) — give an input with nothing
to list the literal value `none`; a blank entry is a hard parse-time
error once the key is given at all. Omitting the key entirely means no
manual bad-channel removal for anything.

**Thread/device environment variables** (all tools, set in the shell,
not the cfg file):

| Variable | Recommended value | Effect |
|---|---|---|
| `OMP_NUM_THREADS` | your own choice (e.g. 6-8) | CPU thread count |
| `OMP_PROC_BIND` | `close` | Thread affinity (CPU mode) |
| `OMP_PLACES` | `cores` | Thread placement (CPU mode) |
| `OMP_TARGET_OFFLOAD` | `MANDATORY` | GPU mode: `MANDATORY` aborts if no GPU is found (with `default`, a missing GPU or driver gives a CPU run and no error); set `DISABLED` for host-fallback testing |
| `OMP_DEFAULT_DEVICE` | `0` | GPU device index (multi-GPU systems) |

`scripts/run_pipeline.sh` sets thread count from its own
`rmsynth_omp_threads=`/cfg key instead, so these only need setting by
hand when invoking a binary directly.

**FFTW wisdom** (`convolve_cubes`, and `match_cubes` when it convolves):
the FFT plans each run measures are saved to `wisdom/fftw_wisdom_sp.dat`
under the current directory, and later runs of the same size reuse them.
Set `RMTOOL_WISDOM_DIR` to keep them somewhere else, e.g. when the
current directory is not writable (the Docker image sets it to
`/tmp/rmtool_wisdom`). If the file cannot be saved, the tool prints
`WARNING: could not save FFTW wisdom ...` and carries on, but every later
run then re-measures its plans, so two runs on the same input can differ
in the last bits.

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
>
> **Note:** *If the beams don't match, that propagated metadata is
> correspondingly less useful as a description of the run -- and more
> so in multi-band mode, where only the reference band's own
> per-channel beam data is ever captured in the output at all; a
> mismatching non-reference band's own beams are never recorded
> anywhere in the output, only flagged in the run log.*
>
> `infileQ`/`infileU`/`infileI`, `badchan_file`, and the channel-window
> limits use each band's own value regardless of `reference_band`.

**Bad-channel handling (optional):**

| Key | Default | Meaning |
|---|---|---|
| `badchan_file` (alias `global_badchan_file`) | omitted (no removal) | See [Shared Parameters](#shared-parameters) above. |

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
| `path_I` | — | Directory of the I-cube; required if `remove_qu_bias=y`. |
| `infileI` | — | Stokes-I input filename; required if `remove_qu_bias=y`; per-band comma list for multi-band. |

**RM sampling (`ofac`/`fac`/`use_auto_rm_range` required):**

| Key | Default | Meaning |
|---|---|---|
| `ofac` | — (required) | Oversampling factor (≥1); `nrm_out = nrm * ofac`. `4` is a common choice. |
| `fac` | — (required) | λ² conversion constant; set it to `3.14159265358979` (π). |
| `use_auto_rm_range` | — (required) | `0`: use manual `beg_rm`/`end_rm`/`nrm` grid; `1`: derive the RM range from the data (those three become optional overrides). |
| `beg_rm` | `-50.0` | RM grid start (rad/m²); required if `use_auto_rm_range=0`. |
| `end_rm` (alias `max_rm`) | `50.0` | RM grid end; must be > `beg_rm`; required if `use_auto_rm_range=0`. |
| `nrm` (alias `nrm_out`) | `100` | RM samples before oversampling (≥1); required if `use_auto_rm_range=0`. |

**Output format (optional):**

| Key | Default | Meaning |
|---|---|---|
| `output_mode` | `ap` | `ap`: write amplitude+phase; `ri`: write real+imaginary instead. |
| `ap_angle_mode` | `phase` | (`output_mode=ap` only) `phase`: PHA cube is `arg(F)`; `pol`: PHA cube is `0.5*arg(F)` (polarization angle), written as `.POLA.RMCUBE.FITS`. |
| `lsq_ref_mode` | `zero` | `zero\|mid\|centroid\|min\|max\|fixed` — the λ² phase reference the dirty cube is built at. `zero` preserves this project's historical convention. The value used is recorded in the output's `LSQREF` FITS header keyword. |
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
| `gpu_vram_mib` | `0` (not set) | VRAM size (MiB) used for sub-block planning. There is no device query: when `0`, the environment variable `GPU_MEM_MIB` is used if set, otherwise 4096 MiB is assumed and a warning is printed. Set it to your card's memory. |
| `mem_frac_vram` | `0.70` | Fraction (0,0.95] of VRAM budgeted in total for the sub-blocks in flight; with two in flight, each gets half. |

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

**Beam metadata:** AMP/PHA, and the PEAK/RM_PEAK/ANG_PEAK/SNR maps when
`cubestat=y`, carry the reference band's `BMAJ`/`BMIN`/`BPA` — plus its
`CASAMBM`/`BEAMS` table, if the reference band's own Q cube has one
(see `reference_band` above for what that table does and doesn't cover
in multi-band mode). MASK/NVALID always carry just the plain scalar,
regardless. The usefulness of this from a science perspective is only when the beams across the spectral axis match. If they do not match, the code does not forbid the rm-synthesis hoping these may be useful for certain debugging. 

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

Reprojects one or more FITS cubes onto one common sky grid using
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
    [interp=<scheme>] [interp_params=<p1>[,<p2>]]
    [log_level=<level>] [timing_enabled=y|n] [log_output_file=<path>] [dry_run=y|n]
    [min_free_pct_after_write=<percent>]
bin/reproject_cubes [target_system=<icrs|fk5|fk4|galactic>] [target_projection=<tan|sin|car|ait|hpx>]
    infiles=<file>[,<file>...]
    [target_center=<infiles_centroid|user_defined>] [target_lon=<lon> target_lat=<lat>]
    [target_cdelt_lon=<deg> target_cdelt_lat=<deg>]
    [mem_frac_ram=<fraction>] [io_overlap=y|n] [nwriters=<n>]
    [interp=<scheme>] [interp_params=<p1>[,<p2>]]
    [log_level=<level>] [timing_enabled=y|n] [log_output_file=<path>] [dry_run=y|n]
    [min_free_pct_after_write=<percent>]
bin/reproject_cubes --config <cfgfile>
bin/reproject_cubes --help | -h
```

No positional arguments — every value is `key=value` (no spaces around
`=`), given on the command line and/or via `--config <cfgfile>` (same
`key=value` file syntax, `#`/`;` comments). If both are given, each CLI
`key=value` overrides only that same key from the config file.

`infiles=` is always required, and exactly one of two things besides
it: `mode=`+`reffile=` together (align onto an existing reference
file), or `target_system=` and/or `target_projection=` with no
`reffile=` at all (build a fresh output grid from the infiles' own
combined footprint — useful when there's no natural "reference band"
to pick, e.g. a single file you want re-expressed in a different sky
system or projection). Giving `reffile=`
together with any `target_*` key is a hard, named error, not a silent
ignore.

### Parameters

| Key | Default | Required | Meaning |
|---|---|---|---|
| `mode` | — | yes, with `reffile=` | `intersection`: output grid shrinks to the overlap of all inputs with the reference. `union`: output grid grows to cover all inputs and the reference. `reference`: output grid is the reference file's own extent. Zero overlap between any input and the output grid is always a hard failure, in every mode. |
| `reffile` | — | yes, with `mode=` | Reference FITS cube whose WCS anchors the output grid. Mutually exclusive with every `target_*` key below. |
| `infiles` | — | yes | Comma-separated list of 1–50 input FITS cube paths. |
| `target_system` | the infiles' own system | with no `reffile=`, this or `target_projection` | `icrs\|fk5\|fk4\|galactic` — builds a fresh output grid from scratch, fit to the infiles' own combined footprint, instead of aligning onto a reference file. Omitted (with `target_projection=` given), the infiles' own system is kept; they must all be in the same one, or the run stops and lists each file's system. |
| `target_projection` | `sin` | with no `reffile=`, this or `target_system` | `tan\|sin\|car\|ait\|hpx`, any case — the projection of the fresh grid. `tan`/`sin` suit a pointed field; `car`/`ait`/`hpx` a wide or all-sky one. `sin` samples the sky more coarsely away from the centre (by 1/cos of the angle from it) and cannot reach beyond 90° from it, so a wide combined footprint is better served by `car`/`ait`/`hpx`. Not allowed together with `reffile=` (a hard error, like every `target_*` key): the reffile alone defines the output grid, projection included. `hpx` is different from the other four: it is HEALPix-aligned, not a free-centre/free-scale field — see `target_hpx_nside=` below, which also lists what becomes an error with it. |
| `target_center` | `infiles_centroid` | no | `infiles_centroid`: new field centre is the vector-averaged centroid of every infile's own pixel-grid centre (not naive lon/lat averaging). `user_defined`: centre given via `target_lon=`/`target_lat=` (required together). Not allowed with `target_projection=hpx` (see `target_hpx_nside=`). |
| `target_lon` / `target_lat` | — | with `target_center=user_defined` | Coordinate value in `target_system`. Accepts plain decimal degrees, colon-sexagesimal (`12:34:56.7`), or letter-suffixed sexagesimal (`12h34m56.7s`, `-45d30m00s`, `180.5d`) — same notation as `scripts/casda_fetch.py`'s own `--lon`/`--lat`. A bare colon-sexagesimal `target_lon` reads as hour angle for `icrs`/`fk5`/`fk4`, as degrees for `galactic`. Not allowed with `target_projection=hpx` (see `target_hpx_nside=`). |
| `target_cdelt_lon` / `target_cdelt_lat` | smallest \|CDELT\| across the infiles, per axis | no (give both or neither) | Pixel scale (degrees/pixel) for the new grid's two axes. The default is the smallest \|CDELT\| found across every infile, taken separately for the longitude and latitude axes, so at the grid centre the output samples every input at least as finely as that input's own grid. Away from the centre the projection stretches or compresses the sky per pixel (`sin` coarsens it by 1/cos of the angle from the centre; see `target_projection`). With one infile that is its own scale. Where the infiles' scales differ, the run prints both the value chosen and the spread. Cells coarser than an infile's own are refused when its beam spans fewer than 2 of them; see [The output grid must sample the beam](#the-output-grid-must-sample-the-beam). Not allowed with `target_projection=hpx` (see `target_hpx_nside=`). |
| `target_hpx_nside` | closest to the finest infile pixel | only with `target_projection=hpx` | Power of 2. An HPX grid's centre, pixel scale and reference pixel all follow from Nside alone — HEALPix tiles the whole sky at a fixed orientation once Nside is chosen, so there is no field to centre the way the other four projections have one; `target_center=`/`target_lon=`/`target_lat=`/`target_cdelt_lon=`/`target_cdelt_lat=` are errors with `target_projection=hpx`, not silently ignored. Omitted, the default Nside is the one whose cell is closest to the finest infile pixel; it may only come out coarser than that if every infile's beam spans ≥ 2.5 of the chosen cell (never for a no-beam infile, which is refused instead; see [The output grid must sample the beam](#the-output-grid-must-sample-the-beam)). Like every `target_*` key, this applies only with no `reffile=`: a supplied `reffile=` is used exactly as given, HPX or not, with no check that it is a genuine, Nside-aligned grid. |
| `interp` | `lanczos3` | no | `nearest\|linear\|sinc\|sincsinc\|sinccos\|sincgauss\|blockave\|gauss\|somb\|sombcos\|lanczos2\|lanczos3\|lanczos\|exact` — interpolation kernel for resampling, or `exact` for flux-conserving overlap reprojection (not a kernel — see "Choosing `interp=`" below, which also covers the measured accuracy of every kernel). Avoid `sincsinc` (a confirmed AST bug truncates it before its own window reaches zero); `sinccos` is the clean alternative. |
| `interp_params` | scheme-specific | no | Up to 2 comma-separated numbers: neighbour/support-radius count (`<=0`: let AST choose), then the scheme's own shape value where it has one (sinc width for `sincsinc`/`sinccos`/`sombcos`, Gaussian FWHM for `sincgauss`/`gauss`). `sinc`/`somb`/`sincsinc`/`sinccos`/`sombcos` use the first number as that neighbour/support-radius count (an omitted value lets AST choose it automatically); `blockave` uses it as its own box half-width instead (no automatic choice — see below). `gauss`/`sincgauss`/`blockave` take a width AST cannot pick itself; omitting it uses a measured default and the value is printed (see "Choosing `interp=`"). Required only for `lanczos`, the general form, where no default `a` exists; rejected for `lanczos2`/`lanczos3` (support is fixed). Unused only for `nearest`/`linear`. |
| `mem_frac_ram` | `0.25` | no | Fraction (0,0.95] of system RAM budgeted for one read/resample/write block of planes (see [Shared Parameters](#shared-parameters) above). Parallelism is scoped to a single block — blocks run strictly sequentially, never overlapped — so on a machine where available memory is tight enough to shrink a block below your thread count, you silently lose most of your OpenMP speedup too, not just I/O efficiency. A WARNING prints if this looks likely. |
| `io_overlap` | `n` | no | `y`: write each block on a background thread, overlapped with the *next* block's read+resample. Only one background write is ever in flight. |
| `nwriters` | `1` | no | `>1`: split each block's own write across this many disjoint writer threads (clamped to `[1, OMP_NUM_THREADS]`) instead of one. Only useful alongside `io_overlap=y` on a disk that benefits from concurrent writes (e.g. NVMe) — see `dry_run` below for a per-machine suggestion. |
| `log_level` | `info` | no | `error\|warn\|info\|debug`. At `info` and above, a per-file stage-timing summary (seconds and % of that file's own total, per stage — e.g. `block_read`, `block_convolve`/`reproject_compute`, `block_write`) is always printed after each file finishes, regardless of `timing_enabled`. `debug` additionally logs one `tile_thread` line per block per worker thread (read/process/write start and end times) — fine-grained enough to drive `scripts/plot_tile_async_swimlane.py`'s per-thread timing plot, useful for seeing exactly where threads stall or serialize. |
| `timing_enabled` | `n` | no | Print a whole-run (all files combined) stage timing summary at exit, in addition to the always-on per-file summary described under `log_level` above. |
| `log_output_file` | empty (stdout) | no | Non-empty: append log/timing output to this file instead. |
| `dry_run` | `n` | no | `y`: check the first `infiles=` entry's own target disk (rotational vs. SSD/NVMe), print a suggested `io_overlap`/`nwriters` starting point, write it to `reproject_cubes_dryrun.cfg`, and report the disk-space check below — then exit. No pixel data read, no output written. |
| `min_free_pct_after_write` | `10` | no | Percent `[0,100]` of the output disk's total size that must still be free once every output file is written. Checked before any file is created: if the run would leave less, it refuses and prints the space needed, free now, what would remain, and the minimum. |

What an HPX grid's output actually looks like (the 45° rotation,
diamond footprint, and why) is shown in
[docs/user/EXAMPLES.md §8a](EXAMPLES.md#8a-target_projectionhpx--square-cells-diamond-footprint).

### Output files

One reprojected cube per input, named `<infile-without-extension>_REPROJ.FITS` (not configurable). A pre-existing file at that path aborts the run before anything is touched.

### Choosing `interp=`

What reprojection preserves is the brightness at a given sky position:
each output pixel should hold the sky's own value at that pixel's sky
coordinate. Jy/beam passes through unscaled, source positions follow
the WCS mapping, the beam is unchanged, and integrated flux follows
from those. What is left is interpolation error, which is inherent in
resampling a sampled image — the choice of kernel sets how large it is.

Measured against analytic ground truth, as mean RMS error, each scheme
as a user gets it with no `interp_params=`, on two source widths —
`compact` (sigma 0.8-0.9 pixels, the narrowest this measurement can
represent without itself losing flux at the input-sampling stage —
see the caveat below) and `extended` (sigma 3-4 pixels, comfortably
resolved) — reported separately, not pooled into one number (pooling
would let one width's own effect distort the other's ranking; see the
caveat below for why). Grids: a same-scale shifted grid, a rotated
grid, a 2x finer grid and a 2x coarser one; on the coarser grid only
`extended` is included (see the beam-sampling rule below) — `compact`'s
own narrowest beam does not clear the 2.0-cell sampling floor once the
output cells double, so it is refused there by design, the same as any
other under-sampled input.

**Interpolation accuracy versus conservative-remapping accuracy —
what the RMS and flux-ratio columns each actually measure, and why
`compact` is the cleaner RMS read.** These are two different,
established operations, not one scheme failing where another
succeeds — the distinction geophysical regridding tools (ESMF/xESMF/
CDO) draw between "interpolation" and "conservative remapping", not a
framing specific to this tool:

- **Interpolation** (every scheme in the table below) treats the input
  as point samples of a continuous, band-limited signal and
  reconstructs that signal's value *at a point*, via the
  Shannon-Whittaker sampling theorem a sinc-family kernel implements.
  **RMS error** measures exactly this: at every output pixel, the
  analytic source function's true value is evaluated at that pixel's
  own exact sky coordinate (never area-averaged), and RMS error is the
  root-mean-square of (reprojected output pixel) minus (that point
  value), over every pixel.
- **Conservative remapping** (`interp=exact`, not in the table below —
  see below) treats each input pixel as a piecewise-constant cell and
  computes the exact geometric overlap between input and output cells,
  which by construction exactly preserves the integral (flux) of that
  piecewise-constant field — not an approximation to get closer to, a
  provable identity. **Flux ratio**, not RMS error, is what measures
  this — see below.

A scheme built to answer the first question is under no obligation to
also win the second, and vice versa: a sinc-family kernel's negative
side lobes are what let it resolve fine structure precisely at a
point, but those same negative weights can leave a small, same-signed
residual once many output pixels are summed into a flux — a
structural property of interpolation, not a bug, and not something a
better-tuned kernel removes.

This is also why `extended`'s own RMS column is not a clean
interpolation-accuracy read on its own: for a genuinely curved,
several-pixel-wide source, the true area average over a pixel's
footprint differs measurably from the true point value at that
pixel's centre, even for a hypothetically perfect reconstruction — so
`extended`'s RMS mixes real interpolation error with this point-vs-
area gap, which exists independent of which scheme is used and is not
a property of the scheme being graded. `compact`'s narrower source has
much less curvature across one pixel, so its RMS column is closer to
isolating interpolation accuracy alone — closer, not exact: a
genuinely zero-width delta function would remove the gap entirely, but
is not representable by this test's own point-sampling input generator
without itself losing flux (a narrower attempt measured flux ratios of
0.67-0.89 for every scheme, an input-discretization artifact, not an
interpolation result) — `compact`'s width is practically, not
theoretically, the narrowest clean option here.

| `interp=` | `compact` RMS | `extended` RMS | `interp_params=` |
|---|---|---|---|
| `sincgauss` | 2.2e-4 | 2.9e-5 | optional (default width below) |
| `lanczos3` (default) | 7.4e-4 | 2.7e-4 | not accepted |
| `lanczos2` | 2.0e-3 | 4.5e-4 | not accepted |
| `somb` | 2.2e-3 | 1.5e-3 | optional |
| `linear` | 4.7e-3 | 1.4e-3 | none |
| `sincsinc` | 4.7e-3 | 1.6e-3 | optional |
| `sombcos` | 4.9e-3 | 2.1e-3 | optional |
| `sinccos` | 5.4e-3 | 2.7e-3 | optional |
| `gauss` | 5.6e-3 | 3.0e-3 | optional (default width below) |
| `sinc` | 8.1e-3 | 3.8e-3 | optional |
| `blockave` | 8.1e-3 | 5.3e-3 | optional (default width below) |
| `nearest` | 1.2e-2 | 1.3e-2 | none |

`lanczos3` is the default. On these sources `sincgauss` at its default
width measured lower still, on both widths.

The rows for `sincgauss`, `gauss` and `blockave` depend on the width.
At widths other than their defaults they measure differently: at a
Gaussian FWHM of 2.0 pixels, `sincgauss` 2.8e-3 (`compact`) / 6.4e-4
(`extended`) and `gauss` 1.1e-2 / 4.4e-3; at a box half-width of 2
pixels, `blockave` 1.8e-2 / 9.4e-3. The other nine schemes take no
width.

Two things the table does not cover:

- At a sharp brightness discontinuity no scheme has an advantage —
  every non-trivial one lands in the same error band.
- Error falls as the beam is better sampled. Across beams from 3×3 to
  20×20 pixels, `lanczos3` runs from 0.04% to 0.001% RMS of the source
  peak and `linear` from 0.36% to 0.06%. The kernel matters at every
  sampling, and most in absolute terms where sampling is marginal.

`interp=nearest` is the only scheme that introduces no value absent
from the input, which makes it the only defined choice for data that
is not continuous — see "what not to reproject" below.

Each scheme's support radius also sets how far blanked regions spread,
which the next section measures.

**Kernel widths for `gauss`, `sincgauss` and `blockave`.** These three
need a width that AST cannot choose for itself. Omitting
`interp_params=` uses a measured default and prints the value:

| `interp=` | default width | meaning |
|---|---|---|
| `gauss` | 0.8 px | Gaussian FWHM |
| `sincgauss` | 8.0 px | Gaussian FWHM |
| `blockave` | 1.0 px | box half-width |

Each is the width with the lowest worst-case error across beam
samplings and output geometries. The widths are in pixels and do not
scale with the beam: `gauss` and `blockave` grow less accurate as the
width grows, so their best width is the narrowest that still resamples,
while `sincgauss`'s Gaussian tapers a sinc rather than doing the
interpolating and so runs the other way. Pass `interp_params=` to
override — for `gauss`/`sincgauss` as `0,<FWHM>`, for `blockave` as
`<half-width>`.

**Flux ratio, same two sources, no `interp_params=`.** A different
question from the RMS table above — total integrated flux recovered,
not pointwise accuracy (see `interp=exact` below for the method built
specifically to get this exactly right regardless of source width).
On `extended`, every scheme lands within a narrow 1.0352-1.0379 band —
a near scheme-independent offset at this source width, not something
choosing a different kernel changes. `compact` differentiates schemes
more clearly:

| `interp=` | `compact` flux ratio |
|---|---|
| `sinc` | 1.032919 |
| `somb` | 1.033811 |
| `sinccos` | 1.035100 |
| `lanczos3` (default) | 1.035622 |
| `sincgauss` | 1.035626 |
| `sombcos` | 1.035680 |
| `lanczos2` | 1.035692 |
| `gauss` | 1.036110 |
| `sincsinc` | 1.036466 |
| `linear` | 1.037275 |
| `blockave` | 1.049172 |
| `nearest` | 1.054616 |

`lanczos3`'s flux ratio is close to the group even though it is not
the flux-conserving method — `exact` is the one built to hold this at
1.0 regardless of source width or scheme choice.

**`interp=exact`: conservative remapping, not interpolation.** Every scheme above
convolves: each output pixel is built from nearby input samples
weighted by a kernel, and what is conserved is the sky's own value at
a point (surface brightness in Jy/beam, with the integral following
from that — see above). `interp=exact` conserves flux a different
way: each output pixel is the overlap-area-weighted mean of every
input pixel it covers on the sky, treating every pixel as a spherical
quadrilateral with great-circle edges — the same method as astropy's
`reproject_exact`. The sum of value × pixel solid angle (a source's
flux, for Jy/beam data) is carried over exactly wherever the output
and input pixel sizes differ, which the kernel schemes above only do
approximately. The cost is resolution: averaging over an output
pixel's whole area adds that pixel's own size to the effective beam,
so `exact` is the right choice when getting the total right matters
more than resolving structure within a pixel — not a general
substitute for the kernels above. `interp_params=` is rejected (there
is no kernel width to set). The per-file overlap weights are built
once and kept for the run; their memory comes out of the
`mem_frac_ram` budget, and both the estimate and the actual size
built are printed.

The same convention applies at a coverage edge in both `exact` and
`reproject_exact`: where an infile's coverage ends partway across an
output pixel (a covered pixel next to a blank one), that pixel is
given the mean of its covered part over its *whole* area, by both
methods alike. Summed over a whole plane this adds a small, shared
excess to both tools' total flux (one measured case: +2.0e-5 of a
plane's total, the same excess in both) — not something `exact` gets
wrong relative to astropy, since both follow the identical convention.

`exact`'s own run time is dominated by building the overlap weights, paid
once per input file and comparable to a kernel pass once built.

### The output grid must sample the beam

Reprojection interpolates. It reproduces structure the output grid can
sample and aliases structure it cannot, so going onto coarser cells is
safe only while those cells still sample the beam. Both
`reproject_cubes` and `match_cubes` check this before writing anything.

The rule applies when the output cells are larger than an input's own
cells. The beam's minor axis (`BMIN`), from the header or, with
`CASAMBM=T`, the narrowest channel in the `BEAMS` table, is compared
with the larger side of an output cell:

| `BMIN` / output cell | Result |
|---|---|
| 2.5 or more | runs |
| 2.0 to 2.5 | runs, with a WARNING that the sampling is barely Nyquist |
| under 2.0 | refused: an ERROR names each such input, the output cell size, the smallest acceptable beam, and the `match_cubes` and `convolve_cubes` commands that convolve to a suitable beam first |

Every cell size here, input or output, is measured on the sky at that
grid's reference pixel (`CRPIX`), where the projection adds no
distortion. For `tan`, `sin`, `car` and `ait` that is the pixel scale
the header sets, \|CDELT\| (including any `CD`/`PC` scaling). For `hpx`,
the latitude side is \|CDELT\| × 8/(3π), about 0.85 \|CDELT\|, from that
projection's equal-area geometry. Away from the reference pixel, the
sky a pixel covers changes with the projection (`sin` pixels grow, `tan`
pixels shrink, `car` pixels narrow in longitude only). The WCS describes
that and reprojection follows it, so it is the projection's geometry,
not the pixel scale you set, and is not counted. An output with the
same pixel scale as an input is therefore never treated as coarser,
wherever its reference point lies.

`BMIN` alone is used, so the beam's position angle does not matter: a
Gaussian beam is nowhere narrower than its minor axis.

On a grid no coarser than an input's own, nothing new can alias, so the
run goes ahead; if that input already undersamples its own beam, a
WARNING says so. An input with no beam information (no `BMAJ`/`BMIN`,
no `BEAMS` table) is accepted onto such a grid and refused onto a
coarser one, since aliasing cannot then be ruled out.

The tools never widen a beam themselves. In `match_cubes stages=both`
with `order=convolve_reproject`, the convolution beam (derived or given
as `target_bmaj`/`target_bmin`) is held to the same rule and refused,
not enlarged, when too narrow: to choose it, give `target_bmaj`/
`target_bmin`/`target_bpa`.

### Blanked pixels, and what not to reproject

Bad pixels travel as NaN inside the data cube itself — there is no
separate mask to keep in step with it. Reprojection is strict about
them: if any NaN falls within the interpolation kernel's support, the
output pixel is NaN too, so a blank region grows by roughly the
kernel's support radius. It is never filled in with invented data.
How much it grows depends on `interp=`; for an 8×8 blank block
reprojected onto a half-pixel-shifted grid:

| `interp=` | blanked pixels out | bounding box |
|---|---|---|
| `nearest` | 64 | 8×8 (no growth) |
| `linear` | 81 | 9×9 |
| `lanczos3` (default) | 169 | 13×13 |

If a wider blank margin matters more to you than interpolation
accuracy, `interp=linear` (or `nearest`) keeps it tighter — see
`interp` above for what you give up.

`convolve_cubes`/`match_cubes`'s convolve stage treats blanks
differently, and deliberately so: convolution mixes neighbours by
design, so an output pixel is kept when at least 50% of the kernel's
*weight* came from valid input, and re-blanked otherwise. Neither
operation invents a value for a flagged pixel, but convolution will
fill one in when most of its support was valid, where reprojection
will not.

**Do not reproject categorical products.** Flag cubes, masks and
`COUNT` maps (integer `BITPIX` — this project's own `MASK` cube is
`BITPIX=8`) hold category labels, not measurements, and interpolating
a label is meaningless. Reprojecting a pure 0/1 mask onto a
half-pixel-shifted grid gives:

| `interp=` | output range | fractional pixels |
|---|---|---|
| `nearest` | exactly 0 or 1 | none |
| `linear` | 0 to 1 | yes |
| `lanczos3` (default) | −0.111 to 1.111 | yes |

`lanczos3` and the other sinc-family kernels have negative side lobes,
so they overshoot both ends: a flag of −0.111 or 1.111 means no more
than a flag of 0.5 does. No tool refuses this, so it is yours to
avoid: use `interp=nearest` if such a cube must go onto another grid.
The toolchain itself never does — `rm_synthesis` writes the `MASK` cube
and `rmclean_cubes` reads it, both downstream of any reprojection.

### Example

```bash
make reproject_cubes
bin/reproject_cubes mode=intersection reffile=ref.fits infiles=a.fits,b.fits

# No reference file -- build a fresh grid from the infiles themselves.
# Pixel scale defaults to the smallest |CDELT| across the infiles:
bin/reproject_cubes target_system=icrs infiles=a.fits,b.fits
```

---

## 3. `convolve_cubes`

Convolves all channels, across all input files, to one common angular
resolution — needed whenever channels don't already share one
resolution, whether that's because of combining bands observed at
different resolutions, or because a single band's own native beam
varies channel-to-channel with frequency (the normal case for any cube
with channel-wise varying beams). Unless an explicit target is given
(`target_bmaj`/`target_bmin`/`target_bpa`), that common beam is derived
automatically: each channel's own elliptical beam is sampled around its
edge, pooled across every channel, and fit to the smallest ellipse
containing all of them (see `npts`/`khachiyan_tol` below). Each channel
is then brought to that beam in one FFT step — multiply by (target
beam's Fourier transform ÷ that channel's own beam's Fourier
transform), inverse-transform — which removes the original PSF and
applies the target one together, rather than as two separate
real-space passes. Reads
per-channel beams from a CASA-style `BEAMS` binary table
(`beamfiles=auto`) or a portable ASCII/CSV beam log (see
`cfg/example_beamLog.txt`/`.csv`). A channel missing
from the beam log, or with `BMAJ`/`BMIN`=0, is written as an all-NaN
plane rather than convolved (and is automatically excluded by
`rm_synthesis`'s own NaN detection later). Each input must have 2 sky
axes and exactly one other non-degenerate axis, FREQ; any further axis
(for example STOKES) must have extent 1, so run separate Stokes slices
as separate `infiles`.

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
| `badchan_file` | none | no | See [Shared Parameters](#shared-parameters) above. `convolve_cubes` additionally detects bad channels automatically from each input's own CASA BEAMS table (a degenerate near-zero beam entry) — `badchan_file` is not the only source of rejected channels here. |
| `target_bmaj` / `target_bmin` / `target_bpa` | none (auto-derive) | no | Explicit target beam (arcsec/arcsec/deg) — give all three together to skip automatic common-beam derivation entirely. |
| `max_common_bmaj` | none | no | If the AUTO-derived common beam's BMAJ exceeds this (arcsec), refuse to proceed. Ignored when an explicit target is given. |
| `mem_frac_ram` | `0.25` | no | Fraction (0,0.95] of system RAM budgeted for one read/convolve/write block of planes (see [Shared Parameters](#shared-parameters) above). |
| `npts` | `2000` | no | How finely each channel's own beam ellipse is sampled (points around its edge) before fitting the smallest common beam that contains all of them — higher is a more faithful approximation of each beam's true shape; the default is already generous and this rarely needs changing. Below 12 is a hard parse-time error. |
| `khachiyan_tol` | `1.0e-5` | no | Stopping tolerance for the iterative fit (`npts` above) that derives the common target beam — it keeps refining until its error drops below this; smaller means a tighter fit at more iterations. The default rarely needs changing. |
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
intermediate FITS file written — avoiding that copy matters once cubes
reach 200GB+.
Default chain order is convolve-then-reproject (low-pass-filtering
before resampling avoids baking interpolation/aliasing error into a
marginally-sampled native PSF, and is usually cheaper too).
`stages=convolve`/`both` inherit `convolve_cubes`'s stricter axis scope
(2 sky axes, one non-degenerate FREQ axis, any other axis of extent 1); `stages=reproject` alone keeps the
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
    [target_system=<icrs|fk5|fk4|galactic>] [target_projection=<tan|sin|car|ait|hpx>]
    [target_center=<infiles_centroid|user_defined>]
    [target_lon=<lon> target_lat=<lat>] [target_cdelt_lon=<deg> target_cdelt_lat=<deg>]
    [beamfiles=<spec1>[,<spec2>...]] [badchan_file=<file1>[,<file2>...]]
    [target_bmaj=<arcsec> target_bmin=<arcsec> target_bpa=<deg>]
    [max_common_bmaj=<arcsec>] [mem_frac_ram=<fraction>] [outsuffix=<suffix>]
    [interp=<scheme>] [interp_params=<p1>[,<p2>]]
    [npts=<n>] [khachiyan_tol=<tol>] [manifest=<path>] [io_overlap=y|n] [nwriters=<n>]
    [log_level=<level>] [timing_enabled=y|n] [log_output_file=<path>] [dry_run=y|n]
    [min_free_pct_after_write=<percent>]
bin/match_cubes --config <cfgfile>
bin/match_cubes --help | -h
```

Same `key=value`-only convention as the other two. When `stages`
includes `reproject`, `footprint_mode=`+`reffile=` must be given
together (align onto an existing reference), or both omitted and
`target_system=` builds a fresh output grid from the infiles' own
combined footprint instead (no natural reference band needed) — same
choice, same semantics, as `reproject_cubes`'s own `mode=`/`reffile=`
vs. `target_system=`.

### Parameters

| Key | Default | Required | Meaning |
|---|---|---|---|
| `stages` | — | yes | `reproject \| convolve \| both`. |
| `order` | `convolve_reproject` | no | `convolve_reproject \| reproject_convolve`; only meaningful when `stages=both`. |
| `infiles` | — | yes | 1–50 comma-separated FITS cube paths. |
| `footprint_mode` | — | with `reffile=`, if `stages` includes `reproject` | `intersection \| union \| reference` — same semantics as `reproject_cubes`'s `mode`. |
| `reffile` | — | with `footprint_mode=`, if `stages` includes `reproject` | Reference cube path. Mutually exclusive with every `target_*` key below. |
| `target_system` | the infiles' own system | if `stages` includes `reproject` and `reffile=` is absent, this or `target_projection` | `icrs\|fk5\|fk4\|galactic` — builds a fresh output grid fit to the infiles' own combined footprint. Same semantics as `reproject_cubes`'s own `target_system=` — see its own entry above for `target_center=`/`target_lon=`/`target_lat=`/`target_cdelt_lon=`/`target_cdelt_lat=`. |
| `target_projection` | `sin` | if `stages` includes `reproject` and `reffile=` is absent, this or `target_system` | `tan\|sin\|car\|ait\|hpx` — same semantics as `reproject_cubes`'s own `target_projection=`. |
| `target_center` / `target_lon` / `target_lat` / `target_cdelt_lon` / `target_cdelt_lat` / `target_hpx_nside` | see `reproject_cubes` above | see `reproject_cubes` above | Only meaningful with `target_system=`/`target_projection=`. Same semantics as `reproject_cubes`'s own keys of the same name — `target_hpx_nside=` only with `target_projection=hpx`, the other five errors with it. |
| `beamfiles` | `auto` per input | no | Used when `stages` includes `convolve` — same semantics as `convolve_cubes`. |
| `badchan_file` | none | no | Used when `stages` includes `convolve` — see [Shared Parameters](#shared-parameters) above; same semantics as `convolve_cubes`. |
| `target_bmaj` / `target_bmin` / `target_bpa` | none | no | Used when `stages` includes `convolve` — same semantics as `convolve_cubes`. |
| `max_common_bmaj` | none | no | Used when `stages` includes `convolve` — same semantics as `convolve_cubes`. |
| `mem_frac_ram` | `0.25` | no | Fraction (0,0.95] of RAM budgeted for one read/process/write block (see [Shared Parameters](#shared-parameters) above). |
| `outsuffix` | `_REPROJ.FITS` / `_CONV.FITS` / `_MATCHED.FITS` | no | Default depends on `stages` if not given explicitly. |
| `interp` / `interp_params` | `lanczos3` / scheme-specific | no | Used when `stages` includes `reproject` — same semantics as `reproject_cubes`'s own keys of the same name (see above). |
| `npts` | `2000` | no | Used when `stages` includes `convolve` — same as `convolve_cubes`. |
| `khachiyan_tol` | `1.0e-5` | no | Used when `stages` includes `convolve` — same as `convolve_cubes`. |
| `manifest` | none | no | Path to write a tab-separated `<infile> SKIPPED\|PROCESSED <effective_path>` record, one line per input. Aborts if the manifest path already exists. |
| `io_overlap` | `n` | no | `y`: write each block on a background thread, overlapped with the *next* block's read+process. |
| `nwriters` | `1` | no | `>1`: split each block's own write across this many disjoint writer threads (clamped to `[1, OMP_NUM_THREADS]`) instead of one. Only useful alongside `io_overlap=y` on a disk that benefits from concurrent writes (e.g. NVMe) — see `dry_run` below for a per-machine suggestion. |
| `log_level` | `info` | no | `error\|warn\|info\|debug`. At `info` and above, a per-file stage-timing summary (seconds and % of that file's own total, per stage — e.g. `block_read`, `convolve_compute`, `reproject_compute`, `block_write`) is always printed after each file finishes, regardless of `timing_enabled`. `debug` additionally logs one `tile_thread` line per block per worker thread (read/process/write start and end times) — fine-grained enough to drive `scripts/plot_tile_async_swimlane.py`'s per-thread timing plot, useful for seeing exactly where threads stall or serialize. |
| `timing_enabled` | `n` | no | Print a whole-run (all files combined) stage timing summary at exit, in addition to the always-on per-file summary described under `log_level` above. |
| `log_output_file` | empty (stdout) | no | Non-empty: append log/timing output to this file instead. |
| `dry_run` | `n` | no | `y`: check the first `infiles=` entry's own target disk (rotational vs. SSD/NVMe), print a suggested `io_overlap`/`nwriters` starting point, write it to `match_cubes_dryrun.cfg`, and (when `stages` includes `reproject`) report the disk-space check below — then exit. No pixel data read, no output written. |
| `min_free_pct_after_write` | `10` | no | Used when `stages` includes `reproject` — same semantics as `reproject_cubes`'s own key of the same name. Not applied to `stages=convolve` alone. |

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

One output cube per input (unless skipped as already matched — see above), `<infile-without-extension><outsuffix>`. `manifest=<path>` additionally writes the skip/process record described above.

### Example

```bash
make match_cubes
bin/match_cubes stages=both order=convolve_reproject \
  footprint_mode=intersection reffile=ref.fits infiles=bandA.fits,bandB.fits \
  mem_frac_ram=0.25

# No reference file -- build a fresh grid from the infiles themselves.
# Pixel scale defaults to the smallest |CDELT| across the infiles:
bin/match_cubes stages=both target_system=icrs infiles=bandA.fits,bandB.fits \
  beamfiles=auto,auto
```

---

## 5. `rmclean_cubes`

Standalone RM-CLEAN (Högbom-style deconvolution), driving `rmclean_mod`
(`src/rmclean.f90`, pure computation, no FITS I/O of its own) against a
*dirty* AMP/PHA cube pair that `rm_synthesis` itself wrote, plus
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
    [lsq_ref_compute_mode=native|zero|mid|centroid|min|max|fixed] [lsq_ref_compute_value=<v>]
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
| `abs_flux_floor` | none (opt-in) | Stop the instant a pixel's own peak amplitude drops to/below this literal fixed flux value. Accepts a bare number (native AMP-cube units) or a value suffixed with `Jy`/`mJy`/`uJy` (e.g. `abs_flux_floor=10mJy`; a space before the unit is also accepted), converted against the AMP cube's own `BUNIT`. No noise/baseline adjustment — a pure "brightest remaining feature below X" comparison. |
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
| `lsq_ref_compute_mode` | `mid` (recommended) | `native\|zero\|mid\|centroid\|min\|max\|fixed` — the reference this program's own RMSF table/CLEAN computation uses internally. `mid` minimizes per-pixel escalation-search cost and does not change the RM grid size. |
| `lsq_ref_compute_value` | — | Used only with `lsq_ref_compute_mode=fixed`. |

> The `LSQREF` keyword on every CLEAN/RESID/RESTORED output records
> `lsq_ref_report` — the reference the phase on disk is
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
| `mask_pattern_cache_max` | `4096` | Pixels sharing the same valid-channel mask pattern share one RMSF table: it's built once, on that pattern's first encounter during the per-tile CLEAN pass (not during the lookahead pre-scan below, which only records which patterns exist and where — it never builds a table itself), and reused by every later pixel with the same pattern. Past this many distinct patterns resident at once, additional new patterns fall back to a one-off table built and discarded per pixel (a performance safety valve, not a correctness issue). |
| `cache_eviction_policy` | `belady` | `belady\|hitcount` — which cached pattern to evict once `mask_pattern_cache_max` is reached, rather than simply refusing new entries. `belady`: offline-optimal, informed by a one-time whole-mask-cube pre-scan (Pass 0) that evicts whichever cached pattern isn't needed again for the longest time (or ever) — no tuning parameter, provably the best possible for a given cache size. `hitcount`: a plain least-frequently-used fallback needing no pre-scan — cheaper to start, but not competitive with `belady`; an opt-out, not the recommended choice. Pass 0 has its own safety valve: if distinct patterns exceed 10% of total image pixels (too little reuse for a lookahead cache to be worth its cost), it aborts with a warning and falls back to `hitcount` for that run. |
| `min_valid_chan_frac` | `0.0` (off) | Pixels whose own valid-channel fraction (summed across ALL bands, out of the full channel count) falls below this are never CLEANed — output NaN across every RM-bin instead. `0.0` means every pixel with at least 1 valid channel is CLEANed (original behaviour); e.g. `0.7` skips any pixel with less than 70% of channels valid. Sparse-coverage pixels have a coarser/noisier RMSF regardless of compute spent on them, so skipping them saves time — the skipped pixels are also never registered in Pass 0's pattern registry or the runtime cache. |

When `cache_eviction_policy=belady`, Pass 0 finishes by printing a
**build-time advisory** (to both stdout and `log_output_file` if set):
a per-block, per-total estimate of table-*generation* time only (never
CLEAN time, which depends on the sky and can't be predicted from the
mask alone), derived from a full Belady-cache simulation plus a small
per-block timing sample — not just a naive "first occurrence" count.
Informational only; it does not affect the run. It also writes
`<outfile>.advisory.csv` (one row per block) — fine-grained enough to
drive `scripts/plot_rmclean_advisory.py`'s predicted-vs-actual plot,
useful for checking a run's progress against what Pass 0 expected
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

| File | Type | Units | Meaning |
|---|---|---|---|
| `<outfile>.NITER.MAP.FITS` | int16 | count | CLEAN iterations used at that pixel. |
| `<outfile>.STOP_REASON.MAP.FITS` | int8 | code | Why CLEAN stopped: 0=not cleaned (no valid channels, or below `min_valid_chan_frac`), 1=hit `niter`, 2=stopped at `abs_flux_floor`, 3=stopped at `auto_nsigma`. The code mapping is also written into the FITS header as HISTORY cards, so the file is self-describing without this table. |
| `<outfile>.RESID_PEAK.MAP.FITS` | float32 | same as `CLEAN.AMP` | Peak \|residual\| amplitude across RM bins after CLEAN. NaN where not cleaned. |
| `<outfile>.RESID_RMS.MAP.FITS` | float32 | same as `CLEAN.AMP` | RMS of \|residual\| amplitude across RM bins after CLEAN. NaN where not cleaned. |
| `<outfile>.N_COMPONENTS.MAP.FITS` | int16 | count | Number of nonzero CLEAN components at that pixel. |
| `<outfile>.COMP_RM_SPREAD.MAP.FITS` | float32 | rad/m^2 | Amplitude-weighted RM dispersion (std dev) of the CLEAN components -- a Faraday-complexity metric, distinct from `N_COMPONENTS` alone since several small hits near one true peak shouldn't read as "complex." NaN for 0 or 1 components (undefined/not meaningful, not zero). |
| `<outfile>.TOTAL_POL_FLUX.MAP.FITS` | float32 | same as `CLEAN.AMP` | Unweighted sum of CLEAN component amplitudes -- total polarized flux recovered by CLEAN. 0.0 for a cleaned pixel with zero components (a well-defined "nothing extracted"); NaN only where not cleaned at all. |

Also `<outfile>.advisory.csv`, written when `cache_eviction_policy=belady`
(the default) and any pixel has valid channels: one row per block with mask-pattern cache statistics and the
predicted RMSF-table build time, for `scripts/plot_rmclean_advisory.py`.

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
