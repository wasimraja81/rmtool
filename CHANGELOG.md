# Changelog

All notable changes to this project are documented in this file.

## [2.0] - 2026-07-17

### Added
- Formalized release-cycle documentation:
  - Added this `CHANGELOG.md`.
  - Added `docs/RELEASE_NOTES_2.0.md`.

### Changed
- Tile planner memory accounting split in `src/rm_synthesis.f90`:
  - Host RAM tile budget now uses `bytes_per_tile_pixel_ram`.
  - GPU VRAM sub-block budget now uses `bytes_per_vram_pixel`.
- Updated docs to reflect planner behavior and measured outcomes:
  - CPU full-image Jennifer benchmark improved after planner split.
  - GPU path remained correct but showed a slight regression on the tested environment.
- Extended swim-lane interpretation notes for CPU thread-detail view:
  - Clarified that single-RM-block runs (`nrm_out <= nrm_block_size`) show only odd/non-hatched `cpu_extract` traces.

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
