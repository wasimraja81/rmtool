# Cubestat Implementation TODO

Status legend:
- [ ] not started
- [x] completed
- [-] in progress

## Scope
Implement cubestat outputs from RM cubes with:
- user switch: `cubestat` (on/off)
- hardcoded sigma method: `tail_quantile`
- maps: peak power, RM-at-peak, angle-at-peak, SNR
- host-side optimization for OMP CPU and OMP+GPU(host)

## Tasks

1. [x] Add cubestat config switch parsing
- Parse `cubestat` from cfg.
- Default to off when not provided.
- Thread through runtime flags.
Acceptance:
- Program starts with/without key.
- `cubestat=y` enables feature path; `cubestat=n` disables it.

2. [x] Define four cubestat FITS outputs
- Add output names/handles for:
  - `.PEAK.MAP.FITS`
  - `.RM_PEAK.MAP.FITS`
  - `.ANG_PEAK.MAP.FITS`
  - `.SNR.MAP.FITS`
- Ensure creation/open/close flow mirrors existing outputs.
Acceptance:
- Files are created only when cubestat is enabled.

3. [x] Implement tail-quantile sigma estimator
- Per pixel compute q50 and q16 from RM power profile.
- Sigma formula: `(q50 - q16) / 0.67449`.
Acceptance:
- Deterministic sigma for fixed input.

4. [x] Compute peak, RM-at-peak, angle-at-peak, and SNR maps
- `Pmax = max(P(RM))`
- `RM_peak = RM(argmax(P(RM)))`
- `ANG_peak = ANG(argmax(P(RM)))`
- `SNR = Pmax / sigma`
Acceptance:
- Outputs match analytical checks on synthetic profiles.

5. [x] Wire cubestat into tile pipeline
- Execute immediately after tile RM extraction while tile arrays are in memory.
- Write tile blocks to 2D map outputs.
Acceptance:
- No extra reread of RM cube from disk.

6. [x] Apply HOST_OMP host parallelization
- Parallelize pixel loop when `HOST_OMP=1`.
- Keep serial behavior when `HOST_OMP=0`.
Acceptance:
- OMP build uses host threads for cubestat loops.

7. [x] Add FITS headers and metadata
- Add method metadata to cubestat products.
- Record sigma definition and enable flag.
Acceptance:
- Header keywords present and readable via FITS tools.

8. [x] Handle invalid pixels and sigma guards
- Define behavior for fully invalid/masked pixels.
- Add sigma floor to avoid divide-by-zero.
Acceptance:
- No NaN explosions or divide-by-zero crashes.

9. [x] Extend tests for cubestat outputs
- Add tests to validate all four maps.
- Validate consistency across serial/OMP/GPU-host paths.
Acceptance:
- Test suite passes with cubestat enabled and disabled.

10. [x] Benchmark OMP and GPU-hostomp paths
- Measure added cubestat runtime overhead.
- Compare OMP=1,GPU=0 vs OMP=1,GPU=1 host-side behavior.
Acceptance:
- Benchmark report with timing deltas and conclusions.

## Definition of Done
- All tasks checked [x].
- Existing tests pass.
- New cubestat tests pass.
- Docs updated with usage and outputs.
