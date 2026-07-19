# Release Notes 4.0

Status: **tagged `4.0` on `main`, 2026-07-19** (merged from `develop`).

## Summary

4.0 is a maintainability and documentation milestone, not a new-capability
one. The two files carrying almost all of rmtool's logic —
`src/rm_synthesis.f90` and `src/rm_synthesis_mod.f90` — were restructured
around derived types, config threading, and named helper routines, with
zero change in observable behaviour at every step (verified bit-identical
against a frozen pre-refactor baseline throughout). Alongside that, a real
bug in the swim-lane plotter was found and fixed, and the README gained a
new "Motivation" section explaining rmtool's parallelism model in plain
language for a research audience — itself fact-checked claim-by-claim
against the implementation during review, catching three real overclaims
before they shipped (see "README Motivation section" below).

## Highlights

### Encapsulation refactor (`planning/ENCAPSULATION_REFACTOR_PLAN.md`, tickets T0-T5)

- **Config encapsulation.** All ~56 config values (paths, tiling/memory
  knobs, RM sampling, masking, GPU, I/O-parallelism, logging/timing) are
  bundled into one `rmsynth_config_t` derived type and read directly as
  `cfg%field` at every use site in `rm_synthesis.f90`. `read_cfg_keyval`'s
  signature shrank from ~56 individual arguments to `(cfgfile, cfg,
  status)`. One value now lives in one place, instead of a cfg-parsed
  copy and a separately-mutated local that could in principle drift
  apart — the exact coupling this ticket removes for `tile_ra`/`tile_dec`
  and several other auto-tuned keys.
- **Tile/VRAM planner encapsulation.** The RAM byte-budget accounting,
  auto-tiling policy, safety-shrink loop, and VRAM sub-block sizing —
  previously ~150 lines of inline arithmetic in the middle of the tile
  loop — are now a named `plan_tile` routine operating on a `tile_plan_t`
  bundle in `rm_synthesis_mod.f90`.
- **I/O orchestration cleanup.** Byte-count and read-thread-split
  arithmetic moved into two small named helpers
  (`compute_tile_read_bytes`, `split_channels_across_threads`). The
  per-tile write dispatch's ~47-line field assembly moved into a named
  `populate_write_job` call — the synchronisation code immediately
  around it (the join-before-reuse/join-before-dispatch invariants that
  guard against the two historical production SIGSEGVs documented in
  `docs/ARCHITECTURE.md`) was left completely untouched, in its exact
  original order, confirmed via `git diff` inspection before any test
  was run.

### Swim-lane plotter fix (`scripts/plot_tile_async_swimlane.py`)

The legend/info side panel could silently disappear entirely on plots
with enough content that no candidate layout in a fixed search grid
happened to fit (found on a real run combining GPU two-slot staging with
a full info block) — every artist got removed and nothing was drawn, with
no error raised. A related, subtler case rendered both boxes but let
their rounded corners visually overlap by a few pixels, traced to the
layout check measuring the info box's raw text extent while missing the
padded border box drawn around it (`Text.get_window_extent()` doesn't
include an attached `bbox=` patch). Both are fixed by replacing the
search-based layout with a deterministic one: the side panel is split
into a fixed bottom region for the legend (narrow column count, so it
can't grow wider than the info box next to it) and a top region for the
info box, continuously font-sized to fill whatever the legend left
behind. The two regions are stacked, not independently placed and then
checked for overlap, so they cannot collide by construction.

### README "Motivation" section

A new section explaining rmtool's tiling and parallelism story for a
research audience that isn't necessarily expert in HPC computing —
memory-budget tiling (and that it's a fraction of *total* system memory,
not whatever's momentarily free, and why that's deliberate), why tiles
are shaped as full-width strips (matching the FITS RA-fastest on-disk
layout), the CPU-side "corner turn" needed to get per-pixel spectra
contiguous for the compute kernel (and that the GPU path's preferred
layout already matches the disk order, so it does a masked copy instead
of a transpose), the write/next-tile-read-and-compute overlap, parallel
I/O channels, and an honest account of GPU acceleration being currently
bounded by host-device (PCIe) transfer bandwidth rather than a finished,
fully-optimised path.

Three claims in an early draft were caught and corrected during review,
each against the actual implementation rather than left as plausible-
sounding assertions:
- A "read tile N+1 / compute tile N / write tile N-1, three concurrent
  tiles" framing, when the real mechanism (`io_overlap`) is a 2-way
  overlap — one tile's write on a background thread concurrent with the
  *next* tile's entire read-through-compute sequence on the main thread,
  gated by a hard join-before-dispatch, never three tiles at three
  distinct stages simultaneously.
- A "corner turn" description generalised to all compute, when
  `prepare_gpu_data` (`rm_synthesis_mod.f90`) was confirmed to perform a
  masked copy in the disk-native pixel-fastest order — no transpose at
  all. Only `prepare_cpu_data` transposes, because only the CPU's
  per-pixel inner loop needs channel-fastest, stride-1 access.
- A "no separate code path to maintain per machine" framing that implied
  one GPU-capable binary is safe to run on hardware with no physical GPU
  at all. Checked against `rm_synthesis.f90:470-497`: `use_gpu_actual` is
  set purely from the compile-time `USE_GPU` macro and the runtime
  `cfg%use_gpu` flag, with no `omp_get_num_devices()` or other runtime
  device-count check anywhere in the codebase. The CPU-only-binary
  fallback (warn and proceed on CPU) *is* code-handled and tested; a
  GPU-capable binary requesting GPU on hardware with no device at all is
  genuinely untested and depends entirely on the OpenMP/libgomp runtime's
  own default behaviour, which this project has not verified.

## Validation

- Every encapsulation ticket independently gated: clean 4-variant
  rebuild (0 errors, 0 new warnings), the full 28-test suite, and a
  bit-identical sweep of all 140 archived FITS outputs against a frozen
  pre-refactor (`T0`) baseline, for every ticket from `T0` through the
  final `T5c` batch.
- Additionally validated on a real production-scale run (`io_overlap=y`,
  `io_read_threads=4`, `io_write_threads=2`, real ASKAP-style cube):
  confirmed data integrity, and passed the structural
  no-overlapping-tile-writes check (9/9 tile writes, zero overlapping
  start/done windows) — the same invariant the historical SIGSEGV
  postmortem in `docs/ARCHITECTURE.md` was written to guard, now
  empirically re-confirmed at real scale rather than only by code
  inspection.
- Swim-lane fix verified against two real logs: one that previously
  rendered no legend/info panel at all (now renders both, no overlap),
  and one that already rendered correctly before the fix (unchanged,
  confirmed via the printed layout metrics).

## Compatibility and behaviour notes

- Zero functional/observable behaviour change from the encapsulation
  refactor: every cfg key means exactly what it meant before, with the
  same defaults, same validation rules, same numerical output.
- No cfg keys were added, removed, or renamed.
- The swim-lane plotter change affects only the legend/info panel's
  layout in generated PNGs; it does not change what data is plotted or
  how log files are parsed.

## What shipped in this tag

- `rmsynth_config_t` (config), `tile_plan_t`/`plan_tile` (tile/VRAM
  planner), `populate_write_job`, `compute_tile_read_bytes`,
  `split_channels_across_threads` (I/O orchestration helpers) in
  `rm_synthesis_mod.f90`; `rm_synthesis.f90` reads every config value as
  `cfg%field` directly, with the T1 shallow-encapsulation unpack bridge
  fully retired.
- The swim-lane legend/info panel layout fix.
- The README "Motivation" section.

## What's next (beyond this tag)

- No further encapsulation tickets are queued; `planning/ENCAPSULATION_REFACTOR_PLAN.md`'s
  full ticket list (T0-T5) is complete.
- GPU staging/async pipeline code remains explicitly out of scope for
  refactoring until a structural race-detection test exists for it (see
  the plan's "Explicitly Out of Scope" section) — unchanged from the
  scope decision made before this cycle began.
- Whether a GPU-capable binary can run safely on hardware with no
  physical GPU device is still untested (see "README Motivation section"
  above) — worth a real test if that deployment scenario ever comes up.
