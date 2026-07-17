# Release Notes 2.0

Release date: 2026-07-17

## Summary

rmtool 2.0 focuses on planner correctness and operational clarity for large CPU/GPU workloads.
The major internal change separates host RAM tile budgeting from GPU VRAM sub-block budgeting,
with corresponding documentation and interpretation updates for swim-lane diagnostics.

## Highlights

- Planner refinement in `src/rm_synthesis.f90`:
  - Split memory accounting into:
    - `bytes_per_tile_pixel_ram` (host RAM tile budget)
    - `bytes_per_vram_pixel` (GPU VRAM sub-block budget)
  - Prevents staging-only terms from over-shrinking CPU/non-staging RAM tiles.

- Benchmark observations used during release qualification:
  - CPU full-image Jennifer run improved after planner split.
  - GPU path remained correct with active offload, but showed a small runtime regression on the tested environment.

- Documentation improvements:
  - Updated architecture/timeline design note for planner split behavior.
  - Clarified swim-lane odd/even `cpu_extract` semantics for single-RM-block CPU runs.
  - Added project changelog and explicit release note references from README.

## Compatibility and behavior notes

- Configuration keys and public CLI/build entry points are unchanged.
- RM synthesis numerical behavior remains unchanged by planner split; the change targets tile-size planning.
- Swim-lane visual parity (`odd` solid vs `even` hatched) is expected only when multiple RM blocks are present.

## Validation snapshot

- Build matrix: CPU serial, CPU OMP, GPU offload, GPU offload + host OMP built successfully.
- Tests: full suite passing (`22/22`) during this release cycle.
