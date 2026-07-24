# Changelog

All notable changes to this project are documented in this file.

## [5.0] - in preparation on `multi-band-tomography`

Multi-band Faraday tomography milestone — by far the largest single body
of work this project has shipped at once. Three parts, each usable on
its own but designed to work as a pipeline: `rm_synthesis` itself can now
merge frequency channels from several input files into one RM synthesis
run; two new standalone tools (`reproject_cubes`, `convolve_cubes`)
prepare real, mismatched-geometry/mismatched-resolution bands to actually
be combined that way; and `rm_synthesis`'s own output metadata was fixed
to faithfully carry beam information through the whole chain rather than
silently dropping it. Full design rationale, decisions recorded with the
user, and ticket-by-ticket verification evidence for all of this lives in
`planning/MULTI_BAND_TOMOGRAPHY_PLAN.md` (tickets T0-T12) — this entry is
a summary, not a replacement for that record.

### Added — multi-band RM synthesis (`rm_synthesis`, tickets T1-T9)
- Comma-separated-list config schema for every per-band key (`infileQ`,
  `infileU`, `resiQ`/`slopeQ`/`resiU`/`slopeU`, `infileI`/`path_I`,
  `badchan_file`, and — new this release — `chan_blc`/`chan_trc`/
  `chan_inc`): band count is derived from list length, no separate
  `nbands` key. A config with no commas anywhere behaves exactly as
  before — this is additive, not a breaking change to any existing cfg.
- Unified N-band geometry validation (RA/Dec WCS, NAXIS, frequency-axis
  index) against a `reference_band`, loudly refusing before any compute
  on mismatch — the same exact-equality philosophy the existing Q-vs-U
  check already used, generalized rather than replaced.
- Multi-band frequency/λ² merge: every band's channels concatenated into
  one merged spectrum (no deduplication in overlaps, no sort required —
  the DFT kernel is order-independent) and run through the existing
  single-band RM-synthesis compute path unchanged.
- `δRM`/`max RM scale`/per-band `ΔRM` diagnostic, logged (not
  auto-applied) for multi-band runs, since `use_auto_rm_range=1`'s
  existing heuristic is unsafe across bands — verified against a
  thesis-published table (Table 6.1) to within ~1%.
- Multi-tile multi-band runs (previously single-tile only) — verified
  bit-identical to the single-tile result of the same data.
- Per-band channel sub-range selection (`chan_blc`/`chan_trc`/
  `chan_inc`), independent per band — e.g. reject bad edge channels or
  hand-pick a good sub-range per band.
- Per-band bad-channel files — each band flags its own bad channels via
  its own required file (same required-key convention as `infileQ`).
- GPU offload for multi-band (no compute-path changes needed — the same
  kernel already used by CPU and GPU; verified on real GPU hardware,
  including the two-level VRAM staging path).
- `io_read_threads>1`/`io_overlap` enabled for multi-band (previously
  blocked entirely — found by code inspection, on direct challenge, that
  the restriction was unnecessarily conservative, the same pattern as
  several tickets before it).
- A real thesis-scenario regression (`tests/check_thesis_scenario.py`):
  point-source recovery, Faraday-thick component reveal (P+L combined
  ~9x the P-alone peak amplitude in its own RM window), and F2/F3
  resolved-vs-blended behaviour, reproduced from published Table 6.1
  bands (P: 300/30 MHz, L: 1200/120 MHz).

### Added — cross-band preprocessing toolchain (tickets T10-T11)
- `reproject_cubes`: new standalone tool (own binary, `make
  reproject_cubes`) reprojecting two or more FITS cubes onto one common
  sky grid via Starlink AST + `astResampleR`, with three footprint modes
  (`intersection`/`union`/`reference`), full WCS/header propagation
  (including `CROTA`/`PCi_j`/`CDi_j` sky rotation), `mem_frac_ram`-budgeted
  block I/O, and OpenMP parallelism across planes.
- `gaussft_mod` (`src/gaussft.f90`): new pure-computation module for
  elliptical-Gaussian FFT-domain beam-matching convolution (deconvolve
  from a source PSF, reconvolve to a target PSF), thread-safe for OpenMP
  via a plan-once/execute-many split (FFTW's planner is not thread-safe;
  a single plan's "new-array execute" form is, verified directly).
- `commonbeam_mod` (`src/commonbeam.f90`): new module finding the
  smallest common beam every one of N per-channel PSFs can be
  deconvolved from (convex hull + minimum-volume-enclosing-ellipse via
  Khachiyan's algorithm + Sault/MIRIAD deconvolvability validation),
  verified against the `radio_beam` Python package on real ASKAP data.
- `convolve_cubes`: new standalone tool (own binary, `make
  convolve_cubes`) driving `gaussft_mod`/`commonbeam_mod` to convolve
  cubes — across one or several input files together — to one common
  angular resolution. Reads per-channel beams from a CASA-style `BEAMS`
  binary table or a portable ASCII/CSV beam log (`cfg/
  example_beamLog.txt`/`.csv`, ready-to-adapt examples included); a
  channel is bad if missing from the beam source or listed with BMAJ or
  BMIN equal to 0. `max_common_bmaj` guards against silently convolving
  to an unexpectedly coarse auto-derived resolution.

### Added — beam-metadata propagation (ticket T12)
- `rm_synthesis` now propagates `BMAJ`/`BMIN`/`BPA` to all 8 output
  products (previously propagated none at all). If the input has
  `CASAMBM=T` (a genuinely per-channel-varying beam not yet run through
  `convolve_cubes`), the flux-derived outputs (AMP/PHA, and the
  PEAK/RMPEAK/ANGPEAK/SNR maps when `cubestat=y`) additionally get
  `CASAMBM=T` plus the input's own real per-channel `BEAMS` table
  attached as an extension, plus a `HISTORY` note — deliberately not
  MASK/NVALID, which are validity bookkeeping, not flux data. In
  multi-band mode, every band's own beam metadata is now cross-checked
  against the reference band's, with a runtime warning on mismatch.

### Fixed
- `rm_synthesis` opened its own Q/U/I/mask input cubes `READWRITE`
  despite never writing to any of them (confirmed by grep: no write-type
  CFITSIO call anywhere in the file targets those units), an unnecessary
  risk to irreplaceable input data. Now opened `READONLY`, matching how
  this file's own parallel tile-reader threads for the same files
  already worked.
- `convolve_cubes`' bad-channel detection (both the CASA `BEAMS`-table
  and ASCII/CSV readers) only checked BMAJ for a degenerate (zero) beam
  entry; a channel with BMAJ present but BMIN equal to 0 was silently
  treated as good. Now checks both.
- A latent bug in the per-band channel-count bookkeeping, surfaced (not
  triggered — `subim` was blocked outright for multi-band until the same
  ticket that found it) while adding per-band channel sub-range
  selection: the reference band's own selected-channel count was being
  computed from its raw NAXIS3 rather than its actual selected range.

### Validation
- All 4 build flavours (`scratch/make_all.sh`) clean; full
  `tests/run_tests.sh` 49/49 pass (up from 28 at the start of this
  branch), re-run clean after every change in this release.
- Multi-band `rm_synthesis`: `nbands=1` bit-identical sweep (140/140
  FITS outputs) held after every single ticket in this release, with no
  exceptions — the explicit bar this whole effort was held to throughout.
  Multi-tile-vs-single-tile, per-band-channel-sub-range, and per-band-
  bad-channel-file all verified bit-identical against known-good
  references, not merely "doesn't crash".
- `reproject_cubes`: byte-identical to independently-computed
  (Python/astropy) ground truth at spot-checked pixels; a real
  `FTGSVE` axis-order bug caught by a non-adjacent-sky-axis fixture and
  fixed; 25 repeated stress runs, no failures.
- `gaussft_mod`: identity round-trip and asymmetric-beam cross-check
  against an independent Python implementation; 16-OpenMP-thread
  shared-plan concurrency test bit-identical to serial.
- `commonbeam_mod`: matches `radio_beam` 0.3.9 on a real 286-channel
  ASKAP `BEAMS` table to within 0.003 arcsec (BMAJ/BMIN) and mod-180
  degrees (PA); independently confirmed deconvolvable from all 286 real
  beams.
- `convolve_cubes`: bit-exact identity check (target beam == a
  channel's own native beam reproduces that channel's input exactly,
  validating the SKY-to-PIXEL BPA convention conversion end-to-end);
  smoke-tested against a real cutout of ASKAP data with no NaN/Inf.
- `rm_synthesis` beam propagation: injected real BMAJ/BMIN/BPA and
  `CASAMBM=T`/`BEAMS` cases and confirmed exact propagation to the
  correct output subset only; injected mismatched per-band beams in a
  multi-band run and confirmed the cross-band warning fires correctly
  (and stays silent when bands genuinely match); fed a real
  `convolve_cubes`-produced NaN bad-channel plane into `rm_synthesis`
  with no `badchan_file` and confirmed automatic exclusion via existing
  NaN detection.

### Not yet done
- A full run of the preprocessing toolchain against the complete 23GB
  real ASKAP cube this work targets (only cutouts and synthetic data
  verified so far).
- This branch has not yet merged to `develop`/`main`; `5.0` is not yet
  an actual git tag.

## [4.1] - 2026-07-20

Diagnostics milestone: closes a real gap in the run-timing picture found
while looking into why a real large-cube run's swim-lane plot showed an
oddly long first read and first write, and drops a second build path that
had quietly stopped working.

### Added
- Two new timed stages, `io_read_init` and `io_write_init`, cover work
  that always ran but was never counted anywhere: reading the input
  cubes' true dimensions before anything else can happen, and (when
  `io_write_threads>1`) staking out the output files' full size on disk
  the first time they're closed. Both now get their own timer, their own
  line in the run log (using the same read/write categories and colours
  the swim-lane plot already draws, so existing plots need no changes),
  and their own `io_read_init_sec`/`io_write_init_sec` columns in the
  timing CSV — with real byte counts, not placeholders: true FITS header
  bytes read (via `FTGHSP`) for the read side, and the true declared
  AMP+PHA size for the write side.
- The printed "Timing summary" and its "Macro timing breakdown" now fold
  these two stages into the read/write totals, so the reported total
  finally reconciles against the sum of its own named stages instead of
  leaving an unexplained gap.

### Removed
- A second, CMake-based way of building the project (`CMakeLists.txt`,
  `cmake_build.sh`) that had drifted out of sync with how the code is
  actually built (it never learned the preprocessor switch the real
  build now depends on) and no longer worked at all.
- Three duplicate fixed-form source files (`src/myfits_info.f`,
  `src/printerror.f`, `src/rm_synthesis.f`) left over from before the
  current free-form build path existed; only their `.f90` counterparts
  were ever actually compiled.

### Validation
- Clean 4-variant rebuild (0 new warnings), full 28/28 test suite,
  including the two runs that specifically re-check bit-identical output
  under `io_overlap` and `io_write_threads>1` — the exact code path the
  new write-init timing sits next to.
- New log markers hand-checked on a real run: they appear first, in
  correct time order, ahead of tile 1's own read/write, with byte counts
  that check out exactly against the cube's declared dimensions.

## [4.0] - 2026-07-19

Maintainability and documentation milestone rather than a new-capability
one: the two fortran files carrying almost all of rmtool's logic were
restructured around derived types with zero change in observable
behaviour, a real swim-lane plotter bug was found and fixed, and the
README gained a new "Motivation" section explaining rmtool's parallelism
model for a research (not HPC-expert) audience. 

### Added
- README "Motivation" section: explains tiling for memory, the
  disk-layout-driven tile shape, the CPU-side "corner turn" (and why the
  GPU path skips it), read/write/compute overlap, parallel I/O channels,
  and GPU offload — in plain language.

### Changed
- All ~56 config values (paths, tiling/memory-planning knobs, RM sampling,
  masking, GPU, I/O-parallelism, logging/timing keys) are now bundled into
  one `rmsynth_config_t` derived type and read directly as `cfg%field` at
  every use site, replacing a flat scope of loose local variables that
  `read_cfg_keyval` used to populate through a ~56-argument signature (now
  just `(cfgfile, cfg, status)`). One value, one place it lives, instead
  of a config-parsed copy and a separately-mutated local that could in
  principle drift apart.
- The RAM/VRAM tile-size planner (auto-tiling policy, safety-shrink loop,
  VRAM sub-block sizing) is now a named `plan_tile` routine operating on a
  `tile_plan_t` bundle, instead of ~150 lines of inline arithmetic in the
  middle of the main tile loop.
- The per-tile write dispatch's field assembly is now a named
  `populate_write_job` call; the synchronisation code around it (the
  join-before-reuse/join-before-dispatch invariants documented in
  `docs/ARCHITECTURE.md`, the two safeguards behind a real historical
  production SIGSEGV) is untouched, in its exact original order.
- Read-side byte-count and thread-split arithmetic in the tile loop moved
  into two small named helpers (`compute_tile_read_bytes`,
  `split_channels_across_threads`).

### Validation
- Every ticket independently gated: clean 4-variant rebuild (0 errors, 0
  new warnings), 28/28 tests, and a bit-identical sweep of all 140
  archived FITS outputs against a frozen baseline.
- Additionally validated on a real production-scale run (`io_overlap=y`,
  `io_read_threads=4`, `io_write_threads=2`, real ASKAP-style cube) —
  confirmed data integrity and passed the structural
  no-overlapping-tile-writes check, the same invariant the historical
  postmortem in `docs/ARCHITECTURE.md` was written to guard.

### Fixed
- Swim-lane plotter (`scripts/plot_tile_async_swimlane.py`): the
  legend/info side panel could silently vanish entirely on plots with
  enough content (e.g. a run with both GPU staging slots and a full
  info block) that no candidate layout in a fixed search grid happened
  to fit — every artist got removed and nothing was drawn, with no error.
  A related case rendered both boxes but let their rounded corners
  visually overlap by a few pixels, because the layout check measured
  the info box's raw text extent and missed the padded border box drawn
  around it. Replaced the search entirely with a deterministic layout:
  the side panel is split into a fixed bottom region (legend, narrow
  column count so it can't dwarf the info box) and a top region (info
  box, continuously font-sized to fill whatever the legend left behind).
  The two regions are stacked, not independently placed, so they cannot
  collide by construction, and the info box now uses the space it's
  given instead of stopping at an arbitrary integer point size.

## [3.0] - 2026-07-18

IO-efficiency milestone: parallel reads, genuine parallel writes, async
tile-write overlap, and the crash/correctness work that came with
building them. All planned tickets (T0-T6) are done and validated
end-to-end on real Setonix production hardware, including T6's actual
write-throughput gain (see Validation below). See
`docs/RELEASE_NOTES_3.0.md` for the full writeup.

### Added
- `io_read_threads` cfg key: N independent read-only CFITSIO handles per
  input cube, reading disjoint channel ranges concurrently.
- `io_write_threads` cfg key: N independent Fortran STREAM I/O units
  write disjoint RM-bin byte ranges of the AMP/PHA output cubes directly,
  bypassing CFITSIO's `ftpsse`/handle machinery for pixel data entirely
  (T6) -- see Fixed/Known limitations below for the two-stage history of
  why this replaced an earlier, unsafe design rather than being the
  first approach shipped.
- `io_overlap` cfg key: tile N's write runs on a background POSIX thread
  concurrent with tile N+1's read/mask/prep/compute/cubestat. Uses a raw
  pthread rather than an OpenMP task specifically so it cannot silently
  nest/collapse the existing OpenMP parallel regions used by
  `io_read_threads` and the compute kernel.
- Auto-tiler RAM planning is aware of `io_overlap`'s doubled output-side
  buffers, so `tile_dec` is planned smaller automatically under the same
  `mem_frac_ram` -- no new user-facing memory configuration.
- Swim-lane plotter: I/O read and I/O write render as separate lanes
  (previously one shared lane distinguished only by colour), and every
  plot now includes a stage-time-totals bar panel (seconds and % of wall
  time per stage, largest first).
- `tile_read`/`tile_write` log lines now carry a `bytes=<N>` field, and
  the swim-lane plotter renders a new I/O throughput (MB/s) panel from
  it -- stacked directly below the swim-lane/thread panel, sharing its
  time axis so a dip or spike lines up with the Gantt bar above it.
  Absent (no empty panel) for logs predating this field.
- New regression tests: bit-identical `io_overlap=n` vs `y` comparison, a
  structural "no two tile writes ever overlap" invariant check, and a
  bit-identical `io_write_threads=1` vs `=4` comparison across all 8
  output products (`tests/run_tests.sh` §13-14).
- Full, sectioned cfg reference in `README.md`: every key the parser
  accepts, marked required/required-if/optional with its real default,
  cross-checked against `read_cfg_keyval`'s case statements rather than
  written from memory.

### Fixed
- int64-safe flattened tile indices throughout the tile/scatter/mask
  loops (`ipix_tile`, `pix_base`, `src_idx`, `dst_idx`, `ipix_full`,
  `ipix_sub`, `idx_wts`, and the equivalents in `prepare_cpu_data`/
  `prepare_gpu_data`/`tile_extract_gpu_rm_blocked`/
  `cubestat_tail_quantile_maps`) -- same INT32_MAX overflow class as the
  2.0-era `allocate()` fix, previously missed in the runtime index
  arithmetic that actually reads/writes those buffers.
- `io_write_threads>1` root-caused as unsafe (original design): CFITSIO
  aliases repeat read-write opens of an already-open file onto one shared
  internal buffer (`fits_already_open()`), so the "independent" handles
  corrupt each other under concurrent `ftpsse` calls. Produced a real
  SIGSEGV in testing. Hard-clamped to 1 as an interim measure, then fixed
  permanently by replacing the mechanism entirely (see T6 below).
- `io_overlap` write-vs-write race: an initial version only guarded
  buffer-slot reuse (two tiles apart), not whether the *immediately
  preceding* write (a different slot) had finished -- both share the same
  single FITS handle, so two pthreads could call `ftpsse` on it
  concurrently. Produced a real SIGSEGV on a production-scale run (a
  small leftover tile at the bottom of an image whose height wasn't an
  exact multiple of the tile size raced ahead of the previous tile's
  write). Fixed by unconditionally joining any outstanding write before
  dispatching a new one.
- **T6: genuine write-throughput parallelism.** `io_write_threads>1` no
  longer hard-clamped -- each RM-chunk now writes via an independent
  Fortran STREAM I/O unit directly to its byte offset (from `FTGHAD`),
  bypassing CFITSIO for AMP/PHA pixel data entirely, relying only on the
  POSIX guarantee that concurrent writes to disjoint byte ranges of one
  file are safe. Found and fixed a second bug during this work, before it
  ever shipped: leaving CFITSIO's handle for these files open until
  program exit (as the serial path always had) caused CFITSIO's own
  data-fill-check, at final close, to treat the raw-written pixel data as
  "past its own stale end-of-file" bookkeeping (never updated, since the
  raw writer bypassed CFITSIO) and zero-fill over it -- silent data loss,
  not a crash, and only visible in the final on-disk state after close.
  Fixed by closing CFITSIO's handle for AMP/PHA immediately after
  fetching the byte offset it provides, before any raw write happens.
- Swim-lane plotter: the "CPU stage" row (CPU thread-detail view) was
  silently missing its compute segment -- filtered out on the assumption
  that the per-thread lanes above already covered it, which left the row
  summing to less than the tile's actual non-I/O time and a legend entry
  ("CPU compute") that never had a corresponding bar. Restored; the two
  views aren't redundant (the stage row shows *when* the stage ran as a
  whole, the thread lanes show *how* it was parallelised).
- Swim-lane plotter: the GPU pipeline view's synchronous-fallback path
  (a tile that fits in one VRAM sub-block, so there's no async
  double-buffering) expected an old `send N/M` log format the Fortran
  code no longer emits -- current single-shot GPU compute logs plain
  `gpu send`/`gpu recv` notes instead, so non-staged GPU runs were
  silently rendering with no GPU lane at all. Fixed to recognize the
  current format.

### Validation
- Full build matrix remained successful (`OMP/GPU` combinations, zero
  compiler warnings).
- Test suite green (`28/28`), including a full bit-for-bit
  `io_write_threads=1` vs `=4` output comparison and manual verification
  that `io_write_threads>1` combined with `io_overlap=y` also produces
  bit-identical output.
- End-to-end production-scale validation on real Setonix hardware and
  ASKAP/EMU data (13308x11870, 288 channels): the exact case that
  originally crashed now completes without error, `io_read_threads`/
  `io_overlap` confirmed on the real workload, and T6's write-throughput
  gain measured directly -- `io_write_threads=8` dropped write from
  2479.9s (96% of wall time) to 108.3s (6%), a ~23x reduction, taking
  total wall time from 2586.7s to 1945.4s (~25% faster end-to-end). Full
  before/after table in `docs/RELEASE_NOTES_3.0.md`.

## [2.0] - 2026-07-17

### Added
- Formalized release-cycle documentation:
  - Added this `CHANGELOG.md`.
  - Added `docs/RELEASE_NOTES_2.0.md`.

### Changed
- Tile planner memory accounting split in `src/rm_synthesis.f90`:
  - Host RAM tile budget now uses `bytes_per_tile_pixel_ram`.
  - GPU VRAM sub-block budget now uses `bytes_per_vram_pixel`.
- Updated docs to reflect planner behaviour and measured outcomes:
  - CPU full-image Jennifer benchmark improved after planner split.
  - GPU path remained correct but showed a slight regression on the tested environment.
- Extended swim-lane interpretation notes for CPU thread-detail view:
  - Clarified that single-RM-chunk runs (`nrm_out <= nrm_block_size`) show only odd/non-hatched `cpu_extract` traces.

### Validation
- Full build matrix remained successful (`OMP/GPU` combinations).
- Test suite remained green (`22/22`) during these updates.

## [1.1] - 2026-07-16

### Added
- Expanded observability and timeline diagnostics:
  - Structured stage/tile logging enhancements.
  - Swim-lane plotting workflow and example artifacts.

### Changed
- Host OpenMP performance improvements:
  - Staged gather/scatter loop parallelization in `src/rm_synthesis.f90`.
  - Pack/copy parallelization in `src/rm_synthesis_mod.f90` with nested-region guard (`omp_in_parallel`).
- Documentation refresh across build, parallelism, and design notes.

### Fixed
- Docker tag-build helper compatibility for shell invocation:
  - `docker/build_push_from_tag.sh` now re-execs under bash when invoked with `sh`.
