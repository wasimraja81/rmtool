#!/usr/bin/env bash
# run_tests.sh  –  RM-synthesis validation test suite
#
# Run from the repository root:
#   bash tests/run_tests.sh
#
# Tests performed (in order):
#   1. Generate synthetic Q/U FITS cubes (two sources at known RMs)
#   2. Build CPU-serial binary  (make GPU=0 OMP=0)
#   3. Build CPU-OpenMP binary  (make OMP=1 GPU=0)
#   4. Build GPU binary         (make GPU=1)  [skipped if no suitable compiler]
#   5. Run serial binary  → check RM peaks at source positions
#   6. Run OMP   binary   → bit-identical to serial reference
#   7. Run GPU   binary   → within rtol=1e-4 of serial reference
#   8-9. Staging tests (GPU only, if available)
#  10. Bad channel masking – per-channel NaN and fully-masked pixel handling
#      (run for serial, OMP, and GPU binaries)
#  11. Cubestat outputs – peak/RM-peak/angle-peak/SNR map validation
#      (serial path, cubestat=y)
#  12. Timing report + CSV validation
#  13. io_overlap (async tile write) – bit-identical to io_overlap=n across
#      a 7-tile run (odd tile count exercises the ping-pong buffer join
#      and end-of-loop cleanup, not just a single-tile no-op)
#  14. nwriters>1 – bit-identical to nwriters=1 across a
#      7-tile run (T6 raw-write path bypasses CFITSIO for AMP/PHA pixel
#      writes; guards against the CFITSIO handle-aliasing bug and the
#      stale-buffer-at-close bug that raw-write mode replaced it with)
#  15. Multi-band tomography (T1) – comma-separated-list infileQ/infileU
#      config schema: matched-geometry two-band config validates and runs
#      to completion (T2 replaced the earlier "not yet implemented" stop);
#      mismatched-geometry config is loudly refused; inconsistent per-band
#      list lengths are rejected at config-parse time
#      (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md)
#  16. Multi-band frequency merge (T2) – Sec 10 thesis-grounded scenario
#      (Raja 2014 Table 6.1/6.2): P-band alone, L-band alone, and P+L
#      combined, for a point source + Faraday-thick top-hat + F2/F3 pair.
#      Also T4: multi-tile multi-band produces bit-identical output to
#      the single-tile P+L run above (tiling must not change the answer)
#  17. Split-band identity test (T5) – a contiguous 2-band split of the
#      primary test cube (no gap) must reproduce the undivided cube's own
#      output bit-for-bit; the most direct mechanical regression check
#      for the frequency-merge architecture itself
#  18. Per-band channel sub-range selection (T6) – a 2-band run where band 2
#      is the full undivided primary test cube restricted via per-band
#      subim_chan_blc/trc down to exactly the T5 split's high-channel range
#      must reproduce T5's own bit-identical result; exercises the T6
#      per-band offset/count/z1-shift arithmetic specifically
#  19. Per-band bad-channel files (T7) – flagging raw channel 150 via a
#      single-band badchan_file on the undivided primary test cube must
#      reproduce, bit-identically, a 2-band split run that flags the same
#      raw channel via band 2's own badchan_file (channel 50 in band 2's
#      own numbering, since band 2 starts at original channel 101)
#  20. GPU offload for multi-band (T8) – the T5 split-band config run
#      through the GPU binary: non-staged within rtol=2e-3 of the CPU
#      reference (matching every existing GPU test's own tolerance) plus
#      RM-peak validation; staged (VRAM sub-block path, gpu_vram_mib=1)
#      bit-identical to the non-staged GPU run
#  21. io_overlap for multi-band (T9) – the T5 split-band config forced
#      into 7 tiles with io_overlap=y, bit-identical to the existing
#      single-tile split-band reference; confirms no overlapping tile
#      writes
#  22. io_read_threads>1 for multi-band (T9) – same config, 7 tiles,
#      io_read_threads=4, bit-identical to the same single-tile reference
#  23. RM-CLEAN thesis scenario (docs/dev/RMCLEAN_INTEGRATION_PLAN.md T1) –
#      single line-of-sight reproduction of Raja (2014) Chapter 6 Figures
#      6.1/6.2/6.3 (Table 6.1/6.2 exact): point source + Faraday-thick
#      top-hat, cleaned from P-band alone, L-band alone, and P+L combined.
#      Standalone Fortran program (rmclean_mod has no FITS I/O/binary yet,
#      T2), compiled fresh here; skipped if FFTW3 isn't available. Also
#      renders PNGs in the thesis figures' own panel style (tests/output/
#      rmclean_plots/) via tests/plot_thesis_scenario_rmclean.py, skipped
#      gracefully if ~/venv/rmtool's python3 isn't found.
#  24. RM-CLEAN lsq_ref flexibility + derotate_to_lsq_ref – a single
#      point source with a NONZERO intrinsic angle, cleaned twice (once
#      at lsq_ref=0, once at a band-mean reference requiring a far
#      coarser RM grid per get_drm), confirming
#      derotate_to_lsq_ref recovers the same intrinsic polarization
#      angle at lambda_sq=0 in both cases, within a tolerance derived
#      from each case's own grid resolution (not an arbitrary constant).
#  25. RM-CLEAN get_lsq_ref_compute – confirms the midpoint of a
#      channel set's own l_sq extent (not a channel-count-weighted mean,
#      not any single band's own centroid) minimizes get_drm's bound, on
#      a deliberately imbalanced P-band(61ch)/L-band(121ch) combination
#      where those candidates diverge.
#  26. RM-CLEAN get_drm oversample floor – demonstrates WHY the enforced
#      floor is oversample>=2, not the bare 2-point Nyquist value of 1:
#      oversample=1 (computed by hand, since get_drm itself now refuses
#      it) recovers a WRONG chi0 on a single point source, not just an
#      imprecise one; oversample=2 already recovers it correctly.
#  27. RM-CLEAN get_drm enforcement – confirms get_drm actually REFUSES
#      (stops with a nonzero exit) an unsafe oversample<2, rather than
#      merely producing a bad answer if a caller passes one. Inverted
#      pass logic: a nonzero exit IS the expected, correct outcome here.
#
# A summary of PASS/FAIL is printed at the end.
# Exit code: 0 = all passed, 1 = at least one failure.

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TESTS_DIR="$REPO_ROOT/tests"
DATA_DIR="$TESTS_DIR/data"
OUT_DIR="$TESTS_DIR/output"
TRUTH="$DATA_DIR/truth.json"
TEMPLATE="$TESTS_DIR/rmsynth-test.cfg.template"
TIMING_TEMPLATE="$TESTS_DIR/rmsynth-timing.cfg.template"

BIN_SERIAL="$REPO_ROOT/bin/rm_synthesis_release_cpu_serial"
BIN_OMP="$REPO_ROOT/bin/rm_synthesis_release_cpu_omp"
BIN_GPU="$REPO_ROOT/bin/rm_synthesis_release_gpu_offload"
BIN_GPU_HOSTOMP="$REPO_ROOT/bin/rm_synthesis_release_gpu_offload_hostomp"

PASS=0
FAIL=0
SKIP=0

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
pass() { echo "[PASS] $*"; PASS=$((PASS + 1)); }
fail() { echo "[FAIL] $*"; FAIL=$((FAIL + 1)); }
skip() { echo "[SKIP] $*"; SKIP=$((SKIP + 1)); }
section() { echo; echo "──────────────────────────────────────────────────────"; echo "$*"; echo "──────────────────────────────────────────────────────"; }

make_cfg() {
    local tag="$1" use_gpu="$2" extra="${3:-}"
    local out_prefix="$OUT_DIR/$tag"
    local cfg="$OUT_DIR/${tag}.cfg"
    sed -e "s|__DATADIR__|${DATA_DIR}|g" \
        -e "s|__OUTPREFIX__|${out_prefix}|g" \
        -e "s|__USE_GPU__|${use_gpu}|g" \
        "$TEMPLATE" > "$cfg"
    if [[ -n "${extra}" ]]; then
        printf '%s\n' "${extra}" >> "$cfg"
    fi
    echo "$cfg"
}

make_timing_cfg() {
    local tag="$1" use_gpu="$2"
    local out_prefix="$OUT_DIR/$tag"
    local cfg="$OUT_DIR/${tag}.cfg"
    local csv_file="$OUT_DIR/${tag}.timing.csv"
    sed -e "s|__DATADIR__|${DATA_DIR}|g" \
        -e "s|__OUTPREFIX__|${out_prefix}|g" \
        -e "s|__USE_GPU__|${use_gpu}|g" \
        -e "s|__TIMINGCSV__|${csv_file}|g" \
        "$TIMING_TEMPLATE" > "$cfg"
    echo "$cfg"
}

run_binary() {
    local binary="$1" cfg="$2" logfile="$3"
    if "$binary" "$cfg" > "$logfile" 2>&1; then
        return 0
    else
        echo "  Binary exited non-zero; last lines of log:"
        tail -20 "$logfile" | sed 's/^/    /'
        return 1
    fi
}

require_timing_markers() {
    local logfile="$1"
    grep -q "Run summary:" "$logfile" && \
    grep -q "Timing summary (seconds):" "$logfile" && \
    grep -q "Macro timing breakdown:" "$logfile"
}

require_timing_csv_row() {
    local csv_file="$1"
    [[ -f "$csv_file" ]] || return 1
    local nlines
    nlines=$(wc -l < "$csv_file")
    [[ "$nlines" -ge 2 ]] || return 1
    local header_cols data_cols
    header_cols=$(head -1 "$csv_file" | awk -F',' '{print NF}')
    data_cols=$(tail -1 "$csv_file" | awk -F',' '{print NF}')
    [[ "$header_cols" -eq "$data_cols" ]] || return 1
}

# io_overlap must serialize all tile writes against each other (this test
# runs with the nwriters=1 default, so all tiles share a single
# FITS handle), even though writes overlap in time with the *next* tile's
# read/mask/prep/compute. This checks that structural invariant directly from the
# tile_write start/done log markers, since it's a timing-dependent race:
# bit-identical output comparisons on small/fast test data can pass even
# when the underlying dispatch logic would crash on production-scale data
# (this is exactly how a real double-dispatch bug reached production
# before being caught -- see the git history around this function).
require_no_overlapping_tile_writes() {
    local logfile="$1"
    python3 - "$logfile" <<'PYEOF'
import re, sys
from datetime import datetime

path = sys.argv[1]
starts, ends = [], []
pat = re.compile(r'^(\S+) \[\w+\] \[tile_write\] \[tid=\d+\] tile write (start|done)')
with open(path) as f:
    for line in f:
        m = pat.match(line)
        if not m:
            continue
        ts = datetime.strptime(m.group(1), "%Y-%m-%dT%H:%M:%S.%f")
        (starts if m.group(2) == "start" else ends).append(ts)

if len(starts) != len(ends) or len(starts) == 0:
    print(f"[FAIL] expected matched start/done pairs, got {len(starts)} starts, {len(ends)} ends")
    sys.exit(1)

# Writes are dispatched/logged in order, so pair them positionally and
# check each write's start is not before the previous write's done.
overlaps = []
for i in range(1, len(starts)):
    if starts[i] < ends[i-1]:
        overlaps.append((i, starts[i], ends[i-1]))

if overlaps:
    for i, s, prev_end in overlaps:
        print(f"[FAIL] write {i} started at {s} before write {i-1} finished at {prev_end}")
    sys.exit(1)

print(f"[OK] {len(starts)} tile writes, no overlapping start/done windows")
sys.exit(0)
PYEOF
}

# ---------------------------------------------------------------------------
# 0. Prepare output directory
# ---------------------------------------------------------------------------
section "0. Preparing directories"
mkdir -p "$OUT_DIR"

# Clean previous test outputs (binary refuses to overwrite)
rm -f "$OUT_DIR"/serial.*.FITS "$OUT_DIR"/omp.*.FITS "$OUT_DIR"/gpu.*.FITS
rm -f "$OUT_DIR"/mb_match.*.FITS "$OUT_DIR"/mb_mismatch.*.FITS "$OUT_DIR"/mb_lenmismatch.*.FITS
rm -f "$OUT_DIR"/rmc_*.FITS
rm -f "$OUT_DIR"/*.timing.csv
rm -f "$OUT_DIR"/*.cfg "$OUT_DIR"/*.log

# ---------------------------------------------------------------------------
# 1. Generate synthetic test data
# ---------------------------------------------------------------------------
section "1. Generating synthetic Q/U cubes"
if python3 "$TESTS_DIR/make_test_cubes.py"; then
    pass "Synthetic cubes generated"
else
    fail "make_test_cubes.py failed"
    echo "Cannot continue without test data."
    exit 1
fi

# ---------------------------------------------------------------------------
# 2. Build CPU-serial binary
# ---------------------------------------------------------------------------
section "2. Building CPU-serial binary  (make GPU=0 OMP=0)"
cd "$REPO_ROOT"
if make GPU=0 OMP=0 2>&1 | tail -5; then
    if [[ -x "$BIN_SERIAL" ]]; then
        pass "Serial binary built: $BIN_SERIAL"
    else
        fail "Serial binary not found after make: $BIN_SERIAL"
    fi
else
    fail "make GPU=0 OMP=0 failed"
fi

# ---------------------------------------------------------------------------
# 3. Build CPU-OpenMP binary
# ---------------------------------------------------------------------------
section "3. Building CPU-OpenMP binary  (make OMP=1 GPU=0)"
if make OMP=1 GPU=0 2>&1 | tail -5; then
    if [[ -x "$BIN_OMP" ]]; then
        pass "OMP binary built: $BIN_OMP"
    else
        fail "OMP binary not found after make: $BIN_OMP"
    fi
else
    fail "make OMP=1 GPU=0 failed"
fi

# ---------------------------------------------------------------------------
# 4. Build GPU binary  (best-effort)
# ---------------------------------------------------------------------------
section "4. Building GPU binary  (make GPU=1 OMP=0)"
BUILD_GPU=0
if make GPU=1 OMP=0 2>&1 | tail -5; then
    if [[ -x "$BIN_GPU" ]]; then
        pass "GPU binary built: $BIN_GPU"
        BUILD_GPU=1
    else
        fail "GPU binary not found after make: $BIN_GPU"
    fi
else
    skip "make GPU=1 OMP=0 failed (no GPU compiler?); GPU test will be skipped"
fi

# ---------------------------------------------------------------------------
# 4b. Build GPU+HostOMP binary  (best-effort)
# ---------------------------------------------------------------------------
section "4b. Building GPU+HostOMP binary  (make GPU=1 OMP=1)"
BUILD_GPU_HOSTOMP=0
if make GPU=1 OMP=1 2>&1 | tail -5; then
    if [[ -x "$BIN_GPU_HOSTOMP" ]]; then
        pass "GPU+HostOMP binary built: $BIN_GPU_HOSTOMP"
        BUILD_GPU_HOSTOMP=1
    else
        fail "GPU+HostOMP binary not found after make: $BIN_GPU_HOSTOMP"
    fi
else
    skip "make GPU=1 OMP=1 failed; GPU+HostOMP test will be skipped"
fi

# ---------------------------------------------------------------------------
# 5. Run serial binary  → RM peak check
# ---------------------------------------------------------------------------
section "5. Serial binary – RM peak validation"
if [[ -x "$BIN_SERIAL" ]]; then
    cfg_serial=$(make_cfg "serial" "n")
    log_serial="$OUT_DIR/serial.log"
    rm -f "$OUT_DIR"/serial.*.FITS
    if run_binary "$BIN_SERIAL" "$cfg_serial" "$log_serial"; then
        amp_serial="$OUT_DIR/serial.AMP.RMCUBE.FITS"
        if [[ -f "$amp_serial" ]]; then
            if python3 "$TESTS_DIR/check_rm_peak.py" "$amp_serial" "$TRUTH"; then
                pass "Serial: RM peaks at correct positions"
            else
                fail "Serial: RM peak check failed"
            fi
        else
            fail "Serial: AMP output cube not found: $amp_serial"
        fi
    else
        fail "Serial binary did not complete successfully (see $log_serial)"
    fi
else
    skip "Serial binary not available; skipping run"
fi

# ---------------------------------------------------------------------------
# 6. OMP binary  → bit-identical to serial reference
# ---------------------------------------------------------------------------
section "6. OMP binary – bit-identical comparison with serial"
if [[ -x "$BIN_OMP" && -f "${OUT_DIR}/serial.AMP.RMCUBE.FITS" ]]; then
    cfg_omp=$(make_cfg "omp" "n")
    log_omp="$OUT_DIR/omp.log"
    rm -f "$OUT_DIR"/omp.*.FITS
    if run_binary "$BIN_OMP" "$cfg_omp" "$log_omp"; then
        amp_omp="$OUT_DIR/omp.AMP.RMCUBE.FITS"
        if [[ -f "$amp_omp" ]]; then
            if python3 "$TESTS_DIR/compare_cubes.py" \
                    "$OUT_DIR/serial.AMP.RMCUBE.FITS" "$amp_omp" --exact; then
                pass "OMP AMP: bit-identical to serial"
            else
                # FP reassociation in parallel reductions can cause tiny diffs;
                # fall back to a tight relative tolerance
            if python3 "$TESTS_DIR/compare_cubes.py" \
                    "$OUT_DIR/serial.AMP.RMCUBE.FITS" "$amp_omp" \
                    --rtol 1e-4; then
                    pass "OMP AMP: matches serial within rtol=1e-4 (FP reassociation)"
                else
                    fail "OMP AMP: differs from serial beyond rtol=1e-4"
                fi
            fi
        else
            fail "OMP: AMP output cube not found: $amp_omp"
        fi
    else
        fail "OMP binary did not complete successfully (see $log_omp)"
    fi
else
    skip "OMP binary or serial reference not available; skipping comparison"
fi

# ---------------------------------------------------------------------------
# 7. GPU binary  → within rtol=1e-4 of serial reference
# ---------------------------------------------------------------------------
section "7. GPU binary – tolerance comparison with serial"
if [[ "$BUILD_GPU" -eq 1 && -x "$BIN_GPU" && -f "${OUT_DIR}/serial.AMP.RMCUBE.FITS" ]]; then
    # Disable mandatory offload so test runs on host if no physical GPU
    export OMP_TARGET_OFFLOAD="${OMP_TARGET_OFFLOAD:-DISABLED}"
    cfg_gpu=$(make_cfg "gpu" "y")
    log_gpu="$OUT_DIR/gpu.log"
    rm -f "$OUT_DIR"/gpu.*.FITS
    if run_binary "$BIN_GPU" "$cfg_gpu" "$log_gpu"; then
        amp_gpu="$OUT_DIR/gpu.AMP.RMCUBE.FITS"
        if [[ -f "$amp_gpu" ]]; then
            # Also check RM peaks in GPU output
            if python3 "$TESTS_DIR/check_rm_peak.py" "$amp_gpu" "$TRUTH"; then
                pass "GPU: RM peaks at correct positions"
            else
                fail "GPU: RM peak check failed"
            fi
            if python3 "$TESTS_DIR/compare_cubes.py" \
                    "$OUT_DIR/serial.AMP.RMCUBE.FITS" "$amp_gpu" \
                    --rtol 2e-3; then
                pass "GPU AMP: matches serial within rtol=2e-3 (ffast-math vs IEEE)"
            else
                fail "GPU AMP: differs from serial beyond rtol=2e-3"
            fi
        else
            fail "GPU: AMP output cube not found: $amp_gpu"
        fi
    else
        fail "GPU binary did not complete successfully (see $log_gpu)"
    fi
else
    skip "GPU binary or serial reference not available; skipping GPU test"
fi

# ---------------------------------------------------------------------------
# 7b. GPU+HostOMP binary  → within rtol=1e-4 of serial reference
# ---------------------------------------------------------------------------
section "7b. GPU+HostOMP binary – tolerance comparison with serial"
if [[ "$BUILD_GPU_HOSTOMP" -eq 1 && -x "$BIN_GPU_HOSTOMP" && -f "${OUT_DIR}/serial.AMP.RMCUBE.FITS" ]]; then
    # Disable mandatory offload so test runs on host if no physical GPU
    export OMP_TARGET_OFFLOAD="${OMP_TARGET_OFFLOAD:-DISABLED}"
    cfg_gpu_hostomp=$(make_cfg "gpu_hostomp" "y")
    log_gpu_hostomp="$OUT_DIR/gpu_hostomp.log"
    rm -f "$OUT_DIR"/gpu_hostomp.*.FITS
    if run_binary "$BIN_GPU_HOSTOMP" "$cfg_gpu_hostomp" "$log_gpu_hostomp"; then
        amp_gpu_hostomp="$OUT_DIR/gpu_hostomp.AMP.RMCUBE.FITS"
        if [[ -f "$amp_gpu_hostomp" ]]; then
            # Also check RM peaks in GPU+HostOMP output
            if python3 "$TESTS_DIR/check_rm_peak.py" "$amp_gpu_hostomp" "$TRUTH"; then
                pass "GPU+HostOMP: RM peaks at correct positions"
            else
                fail "GPU+HostOMP: RM peak check failed"
            fi
            if python3 "$TESTS_DIR/compare_cubes.py" \
                    "$OUT_DIR/serial.AMP.RMCUBE.FITS" "$amp_gpu_hostomp" \
                    --rtol 2e-3; then
                pass "GPU+HostOMP AMP: matches serial within rtol=2e-3 (ffast-math vs IEEE)"
            else
                fail "GPU+HostOMP AMP: differs from serial beyond rtol=2e-3"
            fi
        else
            fail "GPU+HostOMP: AMP output cube not found: $amp_gpu_hostomp"
        fi
    else
        fail "GPU+HostOMP binary did not complete successfully (see $log_gpu_hostomp)"
    fi
else
    skip "GPU+HostOMP binary or serial reference not available; skipping GPU+HostOMP test"
fi

# ---------------------------------------------------------------------------
# 8. Auto tiling shape – full-RA Dec strips
#    With a budget that fits >=1 full RA row but not the whole image, the
#    auto planner must produce tile_ra == nx (full RA) and tile_dec < ny
#    (a Dec strip), NOT a square sub-tile. This keeps each plane read
#    contiguous on disk (RA is FITS NAXIS1, fastest-varying).
# ---------------------------------------------------------------------------
section "8. Auto tiling shape – full-RA Dec strips"
if [[ -x "$BIN_SERIAL" ]]; then
    NX=$(python3 -c "import json;print(json.load(open('$TRUTH'))['nx'])")
    # Small mem_frac_ram -> partial tile; tile_auto=y -> Dec-strip policy.
    cfg_auto=$(make_cfg "autotile" "n" "tile_auto=y
mem_frac_ram=0.00003
dry_run=y")
    log_auto="$OUT_DIR/autotile.log"
    rm -f "$OUT_DIR"/autotile.*.FITS
    if "$BIN_SERIAL" "$cfg_auto" > "$log_auto" 2>&1; then
        # Parse "tile_ra x tile_dec (RAM read px):   <ra>   <dec>"
        read -r T_RA T_DEC < <(awk -F: '/tile_ra x tile_dec/{print $2; exit}' "$log_auto")
        if [[ "${T_RA}" -eq "${NX}" && "${T_DEC}" -lt "${NX}" && "${T_DEC}" -ge 1 ]]; then
            pass "Auto tiling: full-RA Dec strip (tile_ra=${T_RA}=nx, tile_dec=${T_DEC}<ny)"
        else
            fail "Auto tiling: expected full-RA strip (nx=${NX}), got tile_ra=${T_RA} tile_dec=${T_DEC}"
        fi
    else
        fail "Auto tiling dry-run did not complete (see $log_auto)"
    fi
else
    skip "Serial binary not available; skipping auto-tiling test"
fi

# ---------------------------------------------------------------------------
# 9. Two-level VRAM sub-block staging (GPU) – bit-identical to non-staged GPU
#    Staging is GPU-only (use_staging requires use_gpu_actual=true).
#    Forcing a tiny gpu_vram_mib makes the RAM block subdivide into
#    Dec-strip sub-blocks, exercising the gather/extract/scatter path.
#    Output must be bit-identical to the single-level GPU reference (test 7),
#    since both paths use the same tile_extract_gpu_rm_blocked kernel.
# ---------------------------------------------------------------------------
section "9. VRAM sub-block staging (GPU) – bit-identical to non-staged GPU"
if [[ "$BUILD_GPU" -eq 1 && -x "$BIN_GPU" && -f "${OUT_DIR}/gpu.AMP.RMCUBE.FITS" ]]; then
    export OMP_TARGET_OFFLOAD="${OMP_TARGET_OFFLOAD:-DISABLED}"
    # gpu_vram_mib=1 (MiB) forces ny_sub << tile_dec -> staging path on.
    cfg_stg=$(make_cfg "stage" "y" "gpu_vram_mib=1
mem_frac_vram=0.10")
    log_stg="$OUT_DIR/stage.log"
    rm -f "$OUT_DIR"/stage.*.FITS
    if run_binary "$BIN_GPU" "$cfg_stg" "$log_stg"; then
        if grep -q "Staging sub-blocks:  T" "$log_stg"; then
            amp_stg="$OUT_DIR/stage.AMP.RMCUBE.FITS"
            if [[ -f "$amp_stg" ]]; then
                if python3 "$TESTS_DIR/compare_cubes.py" \
                        "$OUT_DIR/gpu.AMP.RMCUBE.FITS" "$amp_stg" \
                        --exact; then
                    pass "Staging AMP: bit-identical to non-staged GPU"
                else
                    fail "Staging AMP: differs from non-staged GPU (gather/scatter bug?)"
                fi
            else
                fail "Staging: AMP output cube not found: $amp_stg"
            fi
        else
            fail "Staging path was NOT activated (check planner logic)"
        fi
    else
        fail "Staging run did not complete (see $log_stg)"
    fi
else
    skip "GPU binary or GPU reference not available; skipping staging test"
fi

# ---------------------------------------------------------------------------
# 10. Bad channel masking – per-channel NaN and fully-masked pixel handling
#     Tests that:
#     - Pixels with one bad channel still produce valid RM values
#     - Fully-masked pixels output NaN in RM cube
#     - Mask cube correctly shows per-channel masking
#
#     Runs for all three binaries: serial, OMP, GPU
# ---------------------------------------------------------------------------
section "10. Bad channel masking – Serial binary"
if [[ -x "$BIN_SERIAL" ]]; then
    cfg_badchan=$(make_cfg "badchan_serial" "n")
    # Update config to use the bad channel test data
    sed -i 's|TEST\.Q\.FITSCUBE|TEST_BADCHAN.Q.FITSCUBE|g' "$cfg_badchan"
    sed -i 's|TEST\.U\.FITSCUBE|TEST_BADCHAN.U.FITSCUBE|g' "$cfg_badchan"
    log_badchan="$OUT_DIR/badchan_serial.log"
    rm -f "$OUT_DIR"/badchan_serial.*.FITS
    if run_binary "$BIN_SERIAL" "$cfg_badchan" "$log_badchan"; then
        # Validate using Python script
        if python3 "$TESTS_DIR/check_bad_channel_masking.py" "badchan_serial"; then
            pass "Bad channel masking (serial): per-channel NaN handling correct"
        else
            fail "Bad channel masking (serial): validation failed (see above)"
        fi
    else
        fail "Bad channel test (serial) did not complete successfully (see $log_badchan)"
    fi
else
    skip "Serial binary not available; skipping bad channel test"
fi

section "10. Bad channel masking – OMP binary"
if [[ -x "$BIN_OMP" ]]; then
    cfg_badchan=$(make_cfg "badchan_omp" "n")
    # Update config to use the bad channel test data
    sed -i 's|TEST\.Q\.FITSCUBE|TEST_BADCHAN.Q.FITSCUBE|g' "$cfg_badchan"
    sed -i 's|TEST\.U\.FITSCUBE|TEST_BADCHAN.U.FITSCUBE|g' "$cfg_badchan"
    log_badchan="$OUT_DIR/badchan_omp.log"
    rm -f "$OUT_DIR"/badchan_omp.*.FITS
    if run_binary "$BIN_OMP" "$cfg_badchan" "$log_badchan"; then
        # Validate using Python script
        if python3 "$TESTS_DIR/check_bad_channel_masking.py" "badchan_omp"; then
            pass "Bad channel masking (OMP): per-channel NaN handling correct"
        else
            fail "Bad channel masking (OMP): validation failed (see above)"
        fi
    else
        fail "Bad channel test (OMP) did not complete successfully (see $log_badchan)"
    fi
else
    skip "OMP binary not available; skipping bad channel test"
fi

section "10. Bad channel masking – GPU binary"
if [[ "$BUILD_GPU" -eq 1 && -x "$BIN_GPU" ]]; then
    cfg_badchan=$(make_cfg "badchan_gpu" "y")
    # Update config to use the bad channel test data
    sed -i 's|TEST\.Q\.FITSCUBE|TEST_BADCHAN.Q.FITSCUBE|g' "$cfg_badchan"
    sed -i 's|TEST\.U\.FITSCUBE|TEST_BADCHAN.U.FITSCUBE|g' "$cfg_badchan"
    log_badchan="$OUT_DIR/badchan_gpu.log"
    rm -f "$OUT_DIR"/badchan_gpu.*.FITS
    # Disable mandatory offload so test runs on host if no physical GPU
    export OMP_TARGET_OFFLOAD="${OMP_TARGET_OFFLOAD:-DISABLED}"
    if run_binary "$BIN_GPU" "$cfg_badchan" "$log_badchan"; then
        # Validate using Python script
        if python3 "$TESTS_DIR/check_bad_channel_masking.py" "badchan_gpu"; then
            pass "Bad channel masking (GPU): per-channel NaN handling correct"
        else
            fail "Bad channel masking (GPU): validation failed (see above)"
        fi
    else
        fail "Bad channel test (GPU) did not complete successfully (see $log_badchan)"
    fi
else
    skip "GPU binary not available or not built; skipping bad channel GPU test"
fi

# ---------------------------------------------------------------------------
# 11. Cubestat outputs – peak/RM-peak/angle-peak/SNR maps
# ---------------------------------------------------------------------------
section "11. Cubestat outputs – serial validation"
if [[ -x "$BIN_SERIAL" ]]; then
    cfg_cubestat=$(make_cfg "cubestat_serial" "n" "cubestat=y")
    log_cubestat="$OUT_DIR/cubestat_serial.log"
    rm -f "$OUT_DIR"/cubestat_serial.*.FITS
    if run_binary "$BIN_SERIAL" "$cfg_cubestat" "$log_cubestat"; then
        peak_map="$OUT_DIR/cubestat_serial.PEAK.MAP.FITS"
        rm_peak_map="$OUT_DIR/cubestat_serial.RM_PEAK.MAP.FITS"
        ang_peak_map="$OUT_DIR/cubestat_serial.ANG_PEAK.MAP.FITS"
        snr_map="$OUT_DIR/cubestat_serial.SNR.MAP.FITS"
        if [[ -f "$peak_map" && -f "$rm_peak_map" && -f "$ang_peak_map" && -f "$snr_map" ]]; then
            if python3 - "$TRUTH" "$peak_map" "$rm_peak_map" "$ang_peak_map" "$snr_map" <<'PY'
import json, sys, numpy as np
from astropy.io import fits
truth = json.load(open(sys.argv[1]))
peak = fits.getdata(sys.argv[2]).squeeze()
rm_peak = fits.getdata(sys.argv[3]).squeeze()
ang_peak = fits.getdata(sys.argv[4]).squeeze()
snr = fits.getdata(sys.argv[5]).squeeze()
ok = True
for src in truth["sources"]:
    x, y, rm_exp = src["x"], src["y"], float(src["rm"])
    p = float(peak[y, x])
    r = float(rm_peak[y, x])
    a = float(ang_peak[y, x])
    s = float(snr[y, x])
    if not np.isfinite(p) or p <= 0.0:
        print(f"[FAIL] {src['name']}: peak invalid ({p})")
        ok = False
    if not np.isfinite(r):
        print(f"[FAIL] {src['name']}: RM_peak invalid ({r})")
        ok = False
    if abs(r - rm_exp) > 2.0:
        print(f"[FAIL] {src['name']}: RM_peak mismatch expected {rm_exp:+.1f}, got {r:+.2f}")
        ok = False
    if not np.isfinite(a):
        print(f"[FAIL] {src['name']}: ANG_peak invalid ({a})")
        ok = False
    if not np.isfinite(s) or s <= 0.0:
        print(f"[FAIL] {src['name']}: SNR invalid ({s})")
        ok = False
    if ok:
        print(f"[OK] {src['name']}: peak={p:.4g}, rm_peak={r:+.2f}, snr={s:.3f}")
sys.exit(0 if ok else 1)
PY
            then
                pass "Cubestat maps (serial): files present and source values valid"
            else
                fail "Cubestat maps (serial): value validation failed"
            fi
        else
            fail "Cubestat maps (serial): one or more output files missing"
        fi
    else
        fail "Cubestat run (serial) did not complete successfully (see $log_cubestat)"
    fi
else
    skip "Serial binary not available; skipping cubestat map test"
fi

# ---------------------------------------------------------------------------
# 12. Timing report and CSV validation (Phase 7)
# ---------------------------------------------------------------------------
section "12. Timing report + CSV validation"

if [[ -x "$BIN_SERIAL" ]]; then
    cfg_timing_serial=$(make_timing_cfg "timing_serial" "n")
    log_timing_serial="$OUT_DIR/timing_serial.log"
    csv_timing_serial="$OUT_DIR/timing_serial.timing.csv"
    rm -f "$OUT_DIR"/timing_serial.*.FITS "$csv_timing_serial"
    if run_binary "$BIN_SERIAL" "$cfg_timing_serial" "$log_timing_serial"; then
        if require_timing_markers "$log_timing_serial"; then
            pass "Timing markers present (serial)"
        else
            fail "Timing markers missing (serial)"
        fi
        if require_timing_csv_row "$csv_timing_serial"; then
            pass "Timing CSV emitted (serial)"
        else
            fail "Timing CSV missing/invalid (serial)"
        fi
    else
        fail "Timing run failed (serial)"
    fi
else
    skip "Serial binary not available; skipping timing validation"
fi

if [[ "$BUILD_GPU" -eq 1 && -x "$BIN_GPU" ]]; then
    export OMP_TARGET_OFFLOAD="${OMP_TARGET_OFFLOAD:-DISABLED}"
    cfg_timing_gpu=$(make_timing_cfg "timing_gpu" "y")
    log_timing_gpu="$OUT_DIR/timing_gpu.log"
    csv_timing_gpu="$OUT_DIR/timing_gpu.timing.csv"
    rm -f "$OUT_DIR"/timing_gpu.*.FITS "$csv_timing_gpu"
    if run_binary "$BIN_GPU" "$cfg_timing_gpu" "$log_timing_gpu"; then
        if require_timing_markers "$log_timing_gpu"; then
            pass "Timing markers present (GPU)"
        else
            fail "Timing markers missing (GPU)"
        fi
        if require_timing_csv_row "$csv_timing_gpu"; then
            pass "Timing CSV emitted (GPU)"
        else
            fail "Timing CSV missing/invalid (GPU)"
        fi
        if [[ -f "$csv_timing_serial" && -f "$csv_timing_gpu" ]]; then
            mode_serial=$(tail -1 "$csv_timing_serial" | awk -F',' '{print $2}')
            mode_gpu=$(tail -1 "$csv_timing_gpu" | awk -F',' '{print $2}')
            if [[ "$mode_serial" == "cpu_serial" && "$mode_gpu" == "gpu_offload" ]]; then
                pass "CPU vs GPU timing CSV mode labels valid"
            else
                fail "CPU vs GPU timing CSV mode labels unexpected: serial=$mode_serial gpu=$mode_gpu"
            fi
        else
            skip "CPU/GPU CSV pair not available; skipping mode label comparison"
        fi
    else
        fail "Timing run failed (GPU)"
    fi
else
    skip "GPU binary not available; skipping timing GPU validation"
fi

# ---------------------------------------------------------------------------
# 13. io_overlap (async tile write) – must match io_overlap=n exactly
# ---------------------------------------------------------------------------
# tile_ra/tile_dec force 7 tiles (32 Dec rows / 5 per tile, uneven remainder)
# so the ping-pong buffer join logic and the odd-tile-count end-of-loop
# cleanup both get exercised, not just a single-tile no-op path.
section "13. io_overlap – bit-identical to io_overlap=n (async tile write)"
if [[ -x "$BIN_OMP" ]]; then
    cfg_ovl_n=$(make_cfg "ovl_n" "n" "tile_auto=n
tile_ra=32
tile_dec=5
cubestat=y
io_overlap=n")
    cfg_ovl_y=$(make_cfg "ovl_y" "n" "tile_auto=n
tile_ra=32
tile_dec=5
cubestat=y
io_overlap=y
log_level=debug")
    log_ovl_n="$OUT_DIR/ovl_n.log"
    log_ovl_y="$OUT_DIR/ovl_y.log"
    rm -f "$OUT_DIR"/ovl_n.*.FITS "$OUT_DIR"/ovl_y.*.FITS

    if run_binary "$BIN_OMP" "$cfg_ovl_n" "$log_ovl_n" && \
       run_binary "$BIN_OMP" "$cfg_ovl_y" "$log_ovl_y"; then
        n_tiles_n=$(grep -c "Doing tile" "$log_ovl_n" || true)
        n_tiles_y=$(grep -c "Doing tile" "$log_ovl_y" || true)
        if [[ "$n_tiles_n" -gt 1 && "$n_tiles_y" -eq "$n_tiles_n" ]]; then
            pass "io_overlap: multi-tile run confirmed (${n_tiles_y} tiles)"
        else
            fail "io_overlap: expected >1 matching tile count, got n=$n_tiles_n y=$n_tiles_y"
        fi

        all_match=1
        for suffix in AMP.RMCUBE PHA.RMCUBE MASK.CUBE NVALID.MAP \
                      PEAK.MAP RM_PEAK.MAP ANG_PEAK.MAP SNR.MAP; do
            f_n="$OUT_DIR/ovl_n.${suffix}.FITS"
            f_y="$OUT_DIR/ovl_y.${suffix}.FITS"
            if [[ -f "$f_n" && -f "$f_y" ]]; then
                if ! python3 "$TESTS_DIR/compare_cubes.py" "$f_n" "$f_y" --exact \
                        > /dev/null 2>&1; then
                    all_match=0
                    fail "io_overlap: ${suffix} differs from io_overlap=n"
                fi
            else
                all_match=0
                fail "io_overlap: ${suffix} output missing (expected $f_n and $f_y)"
            fi
        done
        if [[ "$all_match" -eq 1 ]]; then
            pass "io_overlap: all 8 output products bit-identical to io_overlap=n"
        fi

        if require_no_overlapping_tile_writes "$log_ovl_y" > /dev/null 2>&1; then
            pass "io_overlap: tile writes never overlap (single-handle serialization holds)"
        else
            fail "io_overlap: two tile writes overlapped in time -- concurrent use of the" \
                 "same FITS handle, will SIGSEGV on production-scale data (see $log_ovl_y)"
        fi
    else
        fail "io_overlap: OMP run failed (see $log_ovl_n / $log_ovl_y)"
    fi
else
    skip "OMP binary not available; skipping io_overlap test"
fi

# ---------------------------------------------------------------------------
# 14. nwriters>1 (T6 raw-write) – bit-identical to nwriters=1
# ---------------------------------------------------------------------------
# nwriters>1 used to open N read-write FTOPEN handles onto the SAME
# output file; CFITSIO aliases them onto one shared buffer
# (fits_already_open()), so concurrent ftpsse() calls on them corrupted that
# buffer -- a real SIGSEGV on a Setonix run. That mechanism is gone: N>1 now
# bypasses CFITSIO for AMP/PHA pixel writes entirely via independent raw
# STREAM-I/O writes to disjoint byte ranges (write_rm_chunk_raw), with the
# CFITSIO handle for those two files closed immediately after FTGHAD -- see
# the "Parallel write handle setup" comment in rm_synthesis.f90 for why the
# handle can't be left open (ffclos's own data-fill-check machinery treats
# out-of-band writes it never saw as "past EOF" and zero-fills over them).
# tile_ra/tile_dec force 7 tiles (32 Dec rows / 5 per tile, uneven
# remainder) so the RM-chunk split logic is exercised across a tile whose
# width equals the full output width (fast path) with an uneven trailing
# tile, not just a single-tile no-op.
section "14. nwriters>1 – bit-identical to nwriters=1 (T6)"
if [[ -x "$BIN_OMP" ]]; then
    cfg_wt1=$(make_cfg "wt1" "n" "tile_auto=n
tile_ra=32
tile_dec=5
cubestat=y
nwriters=1")
    cfg_wt4=$(make_cfg "wt4" "n" "tile_auto=n
tile_ra=32
tile_dec=5
cubestat=y
nwriters=4")
    log_wt1="$OUT_DIR/wt1.log"
    log_wt4="$OUT_DIR/wt4.log"
    rm -f "$OUT_DIR"/wt1.*.FITS "$OUT_DIR"/wt4.*.FITS

    if run_binary "$BIN_OMP" "$cfg_wt1" "$log_wt1" && \
       run_binary "$BIN_OMP" "$cfg_wt4" "$log_wt4"; then
        pass "nwriters=1 and =4: both runs completed without crashing"

        if grep -q "using raw stream writes" "$log_wt4"; then
            pass "nwriters=4: raw-write path confirmed taken"
        else
            fail "nwriters=4: expected raw-write startup message not found in $log_wt4"
        fi

        all_match=1
        for suffix in AMP.RMCUBE PHA.RMCUBE MASK.CUBE NVALID.MAP \
                      PEAK.MAP RM_PEAK.MAP ANG_PEAK.MAP SNR.MAP; do
            f1="$OUT_DIR/wt1.${suffix}.FITS"
            f4="$OUT_DIR/wt4.${suffix}.FITS"
            if [[ -f "$f1" && -f "$f4" ]]; then
                if ! python3 "$TESTS_DIR/compare_cubes.py" "$f1" "$f4" --exact \
                        > /dev/null 2>&1; then
                    all_match=0
                    fail "nwriters=4: ${suffix} differs from nwriters=1"
                fi
            else
                all_match=0
                fail "nwriters=4: ${suffix} output missing (expected $f1 and $f4)"
            fi
        done
        if [[ "$all_match" -eq 1 ]]; then
            pass "nwriters=4: all 8 output products bit-identical to nwriters=1"
        fi
    else
        fail "nwriters test: OMP run failed (see $log_wt1 / $log_wt4)"
    fi
else
    skip "OMP binary not available; skipping nwriters test"
fi

# ---------------------------------------------------------------------------
# 15. Multi-band tomography (T1) – comma-list config schema + geometry
#     validation (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md). Fortran's bare
#     `stop` always exits 0 in this codebase (every existing error path
#     uses it, matching the pattern already relied on elsewhere in this
#     script), so these checks are log-content based, not exit-code based.
# ---------------------------------------------------------------------------
section "15. Multi-band config schema – geometry validation + frequency merge (T1/T2)"

if [[ -x "$BIN_SERIAL" ]]; then
    mb_match_cfg="$OUT_DIR/mb_match.cfg"
    mb_match_log="$OUT_DIR/mb_match.log"
    cat > "$mb_match_cfg" <<CFGEOF
path                = ${DATA_DIR}/
infileQ             = TEST.Q.FITSCUBE,TEST_BAND2.Q.FITSCUBE
infileU             = TEST.U.FITSCUBE,TEST_BAND2.U.FITSCUBE
outfile             = ${OUT_DIR}/mb_match
remove_badchan      = n
global_badchan_file = /dev/null,/dev/null
subim               = n
rem_mean            = 0
remove_qu_bias      = n
resiQ               = 0.0,0.0
slopeQ              = 0.0,0.0
resiU               = 0.0,0.0
slopeU              = 0.0,0.0
ofac                = 1
fac                 = 3.14159265358979
use_auto_rm_range   = 0
beg_rm              = -50.0
end_rm              = 50.0
nrm                 = 201
output_mode         = ap
ap_angle_mode       = phase
write_mask_output   = y
write_nvalid_output = y
use_gpu             = n
CFGEOF
    "$BIN_SERIAL" "$mb_match_cfg" > "$mb_match_log" 2>&1
    if grep -q "Multi-band geometry validated successfully across" "$mb_match_log" && \
       grep -q "rm_synthesis run completed" "$mb_match_log"; then
        pass "Multi-band matched-geometry config: validated 2 bands, ran to completion (T2)"
    else
        fail "Multi-band matched-geometry config: expected validation/completion messages not found (see $mb_match_log)"
    fi

    mb_match_amp="$OUT_DIR/mb_match.AMP.RMCUBE.FITS"
    if [[ -f "$mb_match_amp" ]]; then
        if python3 "$TESTS_DIR/check_rm_peak.py" "$mb_match_amp" "$TRUTH" > /dev/null 2>&1; then
            pass "Multi-band matched-geometry config (T2): src_A/src_B recovered at correct RM from merged P+L-analogue bands"
        else
            fail "Multi-band matched-geometry config (T2): RM peak(s) not recovered correctly (see $mb_match_amp)"
        fi
    else
        fail "Multi-band matched-geometry config (T2): expected output $mb_match_amp not found"
    fi

    mb_mismatch_cfg="$OUT_DIR/mb_mismatch.cfg"
    mb_mismatch_log="$OUT_DIR/mb_mismatch.log"
    sed -e "s|TEST_BAND2\.Q\.FITSCUBE|TEST_BAND2_MISMATCH.Q.FITSCUBE|" \
        -e "s|${OUT_DIR}/mb_match|${OUT_DIR}/mb_mismatch|" \
        "$mb_match_cfg" > "$mb_mismatch_cfg"
    "$BIN_SERIAL" "$mb_mismatch_cfg" > "$mb_mismatch_log" 2>&1
    if grep -q "ERROR: RA WCS mismatch for band" "$mb_mismatch_log" && \
       ! grep -q "Multi-band geometry validated successfully" "$mb_mismatch_log"; then
        pass "Multi-band mismatched-geometry config: loudly refused before compute"
    else
        fail "Multi-band mismatched-geometry config: expected loud-refuse message not found (see $mb_mismatch_log)"
    fi

    mb_lenmismatch_cfg="$OUT_DIR/mb_lenmismatch.cfg"
    mb_lenmismatch_log="$OUT_DIR/mb_lenmismatch.log"
    sed -e "s|infileU             = TEST.U.FITSCUBE,TEST_BAND2.U.FITSCUBE|infileU             = TEST.U.FITSCUBE|" \
        -e "s|${OUT_DIR}/mb_match|${OUT_DIR}/mb_lenmismatch|" \
        "$mb_match_cfg" > "$mb_lenmismatch_cfg"
    "$BIN_SERIAL" "$mb_lenmismatch_cfg" > "$mb_lenmismatch_log" 2>&1
    if grep -q "infileQ/infileU band-count mismatch" "$mb_lenmismatch_log"; then
        pass "Multi-band inconsistent list-length config: rejected at parse time"
    else
        fail "Multi-band inconsistent list-length config: expected band-count-mismatch message not found (see $mb_lenmismatch_log)"
    fi
else
    skip "Serial binary not available; skipping multi-band config tests"
fi

# ---------------------------------------------------------------------------
# 16. Multi-band frequency merge (T2) – Sec 10 thesis-grounded scenario
#     (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md; Raja 2014 Table 6.1/6.2):
#     P-band (300/30 MHz) alone, L-band (1200/120 MHz) alone, and the P+L
#     combined synthesis, for a point source + Faraday-thick top-hat +
#     the F2/F3 close-pair addition.
# ---------------------------------------------------------------------------
section "16. Multi-band frequency merge – Sec 10 thesis scenario (T2)"

if [[ -x "$BIN_SERIAL" ]]; then
    if python3 "$TESTS_DIR/make_thesis_scenario_cubes.py"; then
        pass "Sec 10 thesis-scenario cubes generated (P-band, L-band)"
    else
        fail "make_thesis_scenario_cubes.py failed"
    fi

    thesis_truth="$DATA_DIR/thesis_scenario_truth.json"
    rm -f "$OUT_DIR"/thesis_p.*.FITS "$OUT_DIR"/thesis_l.*.FITS "$OUT_DIR"/thesis_pl.*.FITS

    make_thesis_cfg() {
        local tag="$1" infileQ="$2" infileU="$3" resi_list="$4"
        local cfg="$OUT_DIR/${tag}.cfg"
        # T7: badchan_file is a required per-band list, same cardinality as
        # resi_list -- reuse its comma count rather than hard-coding it.
        local badchan_list
        badchan_list=$(echo "$resi_list" | sed 's/[^,]*/\/dev\/null/g')
        cat > "$cfg" <<CFGEOF
path                = ${DATA_DIR}/
infileQ             = ${infileQ}
infileU             = ${infileU}
outfile             = ${OUT_DIR}/${tag}
remove_badchan      = n
global_badchan_file = ${badchan_list}
subim               = n
rem_mean            = 0
remove_qu_bias      = n
resiQ               = ${resi_list}
slopeQ              = ${resi_list}
resiU               = ${resi_list}
slopeU              = ${resi_list}
ofac                = 1
fac                 = 3.14159265358979
use_auto_rm_range   = 0
beg_rm              = -500.0
end_rm              = 500.0
nrm                 = 501
output_mode         = ap
ap_angle_mode       = phase
write_mask_output   = y
write_nvalid_output = y
use_gpu             = n
CFGEOF
        echo "$cfg"
    }

    thesis_p_cfg=$(make_thesis_cfg thesis_p THESIS_P.Q.FITSCUBE THESIS_P.U.FITSCUBE "0.0")
    thesis_l_cfg=$(make_thesis_cfg thesis_l THESIS_L.Q.FITSCUBE THESIS_L.U.FITSCUBE "0.0")
    thesis_pl_cfg=$(make_thesis_cfg thesis_pl "THESIS_P.Q.FITSCUBE,THESIS_L.Q.FITSCUBE" \
        "THESIS_P.U.FITSCUBE,THESIS_L.U.FITSCUBE" "0.0,0.0")

    thesis_p_log="$OUT_DIR/thesis_p.log"
    thesis_l_log="$OUT_DIR/thesis_l.log"
    thesis_pl_log="$OUT_DIR/thesis_pl.log"

    if run_binary "$BIN_SERIAL" "$thesis_p_cfg" "$thesis_p_log" && \
       run_binary "$BIN_SERIAL" "$thesis_l_cfg" "$thesis_l_log" && \
       run_binary "$BIN_SERIAL" "$thesis_pl_cfg" "$thesis_pl_log"; then
        pass "Sec 10 scenario: P-alone/L-alone/P+L runs all completed"

        thesis_p_amp="$OUT_DIR/thesis_p.AMP.RMCUBE.FITS"
        thesis_l_amp="$OUT_DIR/thesis_l.AMP.RMCUBE.FITS"
        thesis_pl_amp="$OUT_DIR/thesis_pl.AMP.RMCUBE.FITS"
        if python3 "$TESTS_DIR/check_thesis_scenario.py" \
                "$thesis_p_amp" "$thesis_l_amp" "$thesis_pl_amp" "$thesis_truth"; then
            pass "Sec 10 scenario: point source, thick-component washout/reveal, and F2/F3 P-alone/L-alone behaviour all match thesis-grounded expectations"
        else
            fail "Sec 10 scenario: one or more expected behaviours not observed (see check_thesis_scenario.py output above)"
        fi

        # T4 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): multi-tile multi-band
        # must produce bit-identical output to the single-tile run above --
        # tiling must not change the scientific answer.
        thesis_plmt_cfg="$OUT_DIR/thesis_pl_multitile.cfg"
        thesis_plmt_log="$OUT_DIR/thesis_pl_multitile.log"
        rm -f "$OUT_DIR"/thesis_pl_multitile.*.FITS
        { sed -e "s|outfile             = ${OUT_DIR}/thesis_pl\$|outfile             = ${OUT_DIR}/thesis_pl_multitile|" \
              "$thesis_pl_cfg" | grep -v '^tile_'; \
          echo "tile_ra = 16"; echo "tile_dec = 16"; echo "tile_auto = n"; \
        } > "$thesis_plmt_cfg"
        if run_binary "$BIN_SERIAL" "$thesis_plmt_cfg" "$thesis_plmt_log"; then
            if grep -q "Multi-band run spanning.*4  tile(s)" "$thesis_plmt_log"; then
                pass "Multi-tile multi-band (T4): confirmed 4-tile run (tile_ra=tile_dec=16 on a 32x32 image)"
            else
                fail "Multi-tile multi-band (T4): expected 4-tile message not found (see $thesis_plmt_log)"
            fi
            all_match=1
            for suffix in AMP.RMCUBE PHA.RMCUBE MASK.CUBE NVALID.MAP; do
                f1="$thesis_pl_amp"
                [[ "$suffix" != "AMP.RMCUBE" ]] && f1="${thesis_pl_amp/AMP.RMCUBE/$suffix}"
                f2="$OUT_DIR/thesis_pl_multitile.${suffix}.FITS"
                if [[ -f "$f1" && -f "$f2" ]]; then
                    if ! python3 "$TESTS_DIR/compare_cubes.py" "$f1" "$f2" --exact > /dev/null 2>&1; then
                        all_match=0
                        fail "Multi-tile multi-band (T4): ${suffix} differs from single-tile output"
                    fi
                else
                    all_match=0
                    fail "Multi-tile multi-band (T4): ${suffix} output missing (expected $f1 and $f2)"
                fi
            done
            if [[ "$all_match" -eq 1 ]]; then
                pass "Multi-tile multi-band (T4): all 4 output products bit-identical to single-tile multi-band"
            fi
        else
            fail "Multi-tile multi-band (T4): run failed (see $thesis_plmt_log)"
        fi
    else
        fail "Sec 10 scenario: one or more runs failed (see $thesis_p_log / $thesis_l_log / $thesis_pl_log)"
    fi
else
    skip "Serial binary not available; skipping Sec 10 thesis scenario"
fi

# ---------------------------------------------------------------------------
# 17. Split-band identity test (T5, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md):
#     the single most direct, mechanical regression check for the
#     frequency-merge architecture -- splitting TEST.Q/U.FITSCUBE into two
#     CONTIGUOUS halves (no gap) and running them as a 2-band multi-band
#     config must reproduce the undivided cube's own output (section 5's
#     "serial" run) bit-for-bit, not just approximately. Unlike the Sec 10
#     scientific scenario, this has no qualitative-interpretation caveat --
#     a mismatch here always means a real regression, never expected
#     physics.
# ---------------------------------------------------------------------------
section "17. Split-band identity test – contiguous split == undivided cube (T5)"

if [[ -x "$BIN_SERIAL" ]]; then
    split_cfg="$OUT_DIR/split_identity.cfg"
    split_log="$OUT_DIR/split_identity.log"
    rm -f "$OUT_DIR"/split_identity.*.FITS
    cat > "$split_cfg" <<CFGEOF
path                = ${DATA_DIR}/
infileQ             = TEST_SPLIT_LO.Q.FITSCUBE,TEST_SPLIT_HI.Q.FITSCUBE
infileU             = TEST_SPLIT_LO.U.FITSCUBE,TEST_SPLIT_HI.U.FITSCUBE
outfile             = ${OUT_DIR}/split_identity
remove_badchan      = n
global_badchan_file = /dev/null,/dev/null
subim               = n
rem_mean            = 0
remove_qu_bias      = n
resiQ               = 0.0,0.0
slopeQ              = 0.0,0.0
resiU               = 0.0,0.0
slopeU              = 0.0,0.0
ofac                = 1
fac                 = 3.14159265358979
use_auto_rm_range   = 0
beg_rm              = -50.0
end_rm              = 50.0
nrm                 = 201
output_mode         = ap
ap_angle_mode       = phase
write_mask_output   = y
write_nvalid_output = y
use_gpu             = n
CFGEOF
    if run_binary "$BIN_SERIAL" "$split_cfg" "$split_log"; then
        all_match=1
        for suffix in AMP.RMCUBE PHA.RMCUBE MASK.CUBE NVALID.MAP; do
            f1="$OUT_DIR/serial.${suffix}.FITS"
            f2="$OUT_DIR/split_identity.${suffix}.FITS"
            if [[ -f "$f1" && -f "$f2" ]]; then
                if ! python3 "$TESTS_DIR/compare_cubes.py" "$f1" "$f2" --exact \
                        > /dev/null 2>&1; then
                    all_match=0
                    fail "Split-band identity (T5): ${suffix} differs from undivided-cube output"
                fi
            else
                all_match=0
                fail "Split-band identity (T5): ${suffix} output missing (expected $f1 and $f2)"
            fi
        done
        if [[ "$all_match" -eq 1 ]]; then
            pass "Split-band identity (T5): contiguous 2-band split bit-identical to undivided single-band cube"
        fi
    else
        fail "Split-band identity (T5): run failed (see $split_log)"
    fi
else
    skip "Serial binary not available; skipping split-band identity test"
fi

# ---------------------------------------------------------------------------
# 18. Per-band channel sub-range selection (T6, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md):
#     band 1 = TEST_SPLIT_LO (unrestricted, its own full 100-channel range),
#     band 2 = the FULL undivided TEST.Q/U.FITSCUBE (200 channels) with
#     subim_chan_blc/trc restricting it, per-band, down to exactly channels
#     101-200 -- the same channels T5's TEST_SPLIT_HI file contains. Output
#     must be bit-identical to T5's own already-passing result (itself
#     bit-identical to the undivided single-band run). Only passes if the
#     per-band offset/count/z1-shift arithmetic T6 adds is exactly right.
# ---------------------------------------------------------------------------
section "18. Per-band channel sub-range selection – restricted full cube == T5 split (T6)"

if [[ -x "$BIN_SERIAL" ]]; then
    chan_cfg="$OUT_DIR/split_identity_chan.cfg"
    chan_log="$OUT_DIR/split_identity_chan.log"
    rm -f "$OUT_DIR"/split_identity_chan.*.FITS
    cat > "$chan_cfg" <<CFGEOF
path                = ${DATA_DIR}/
infileQ             = TEST_SPLIT_LO.Q.FITSCUBE,TEST.Q.FITSCUBE
infileU             = TEST_SPLIT_LO.U.FITSCUBE,TEST.U.FITSCUBE
outfile             = ${OUT_DIR}/split_identity_chan
remove_badchan      = n
global_badchan_file = /dev/null,/dev/null
subim               = y
subim_chan_blc      = 0,101
subim_chan_trc      = 0,200
rem_mean            = 0
remove_qu_bias      = n
resiQ               = 0.0,0.0
slopeQ              = 0.0,0.0
resiU               = 0.0,0.0
slopeU              = 0.0,0.0
ofac                = 1
fac                 = 3.14159265358979
use_auto_rm_range   = 0
beg_rm              = -50.0
end_rm              = 50.0
nrm                 = 201
output_mode         = ap
ap_angle_mode       = phase
write_mask_output   = y
write_nvalid_output = y
use_gpu             = n
CFGEOF
    if run_binary "$BIN_SERIAL" "$chan_cfg" "$chan_log"; then
        all_match=1
        for suffix in AMP.RMCUBE PHA.RMCUBE MASK.CUBE NVALID.MAP; do
            f1="$OUT_DIR/serial.${suffix}.FITS"
            f2="$OUT_DIR/split_identity_chan.${suffix}.FITS"
            if [[ -f "$f1" && -f "$f2" ]]; then
                if ! python3 "$TESTS_DIR/compare_cubes.py" "$f1" "$f2" --exact \
                        > /dev/null 2>&1; then
                    all_match=0
                    fail "Per-band chan sub-range (T6): ${suffix} differs from undivided-cube output"
                fi
            else
                all_match=0
                fail "Per-band chan sub-range (T6): ${suffix} output missing (expected $f1 and $f2)"
            fi
        done
        if [[ "$all_match" -eq 1 ]]; then
            pass "Per-band chan sub-range (T6): restricted full cube bit-identical to undivided single-band cube"
        fi
    else
        fail "Per-band chan sub-range (T6): run failed (see $chan_log)"
    fi
else
    skip "Serial binary not available; skipping per-band channel sub-range test"
fi

# ---------------------------------------------------------------------------
# 19. Per-band bad-channel files (T7, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md):
#     raw channel 150 (1-indexed, into the undivided 200-channel TEST cube)
#     flagged via a plain single-band badchan_file must reproduce,
#     bit-identically, a 2-band split run (TEST_SPLIT_LO + TEST_SPLIT_HI)
#     that flags the *same* raw channel via band 2's own badchan_file --
#     channel 50 in band 2's own numbering, since band 2 (TEST_SPLIT_HI)
#     starts at original channel 101 (150 - 100 = 50). Exercises the T7
#     per-band bad-channel read/apply code specifically, composed with the
#     existing frequency-merge architecture.
# ---------------------------------------------------------------------------
section "19. Per-band bad-channel files – flagged split-band run == flagged undivided cube (T7)"

if [[ -x "$BIN_SERIAL" ]]; then
    badchan7_ref_list="$OUT_DIR/badchan7_ref_list.txt"
    badchan7_band2_list="$OUT_DIR/badchan7_band2_list.txt"
    badchan7_none_list="$OUT_DIR/badchan7_none_list.txt"
    printf '150\n' > "$badchan7_ref_list"
    printf '50\n' > "$badchan7_band2_list"
    : > "$badchan7_none_list"

    badchan7_ref_cfg="$OUT_DIR/badchan7_ref.cfg"
    badchan7_ref_log="$OUT_DIR/badchan7_ref.log"
    rm -f "$OUT_DIR"/badchan7_ref.*.FITS
    cat > "$badchan7_ref_cfg" <<CFGEOF
path                = ${DATA_DIR}/
infileQ             = TEST.Q.FITSCUBE
infileU             = TEST.U.FITSCUBE
outfile             = ${OUT_DIR}/badchan7_ref
remove_badchan      = y
global_badchan_file = ${badchan7_ref_list}
subim               = n
rem_mean            = 0
remove_qu_bias      = n
resiQ               = 0.0
slopeQ              = 0.0
resiU               = 0.0
slopeU              = 0.0
ofac                = 1
fac                 = 3.14159265358979
use_auto_rm_range   = 0
beg_rm              = -50.0
end_rm              = 50.0
nrm                 = 201
output_mode         = ap
ap_angle_mode       = phase
write_mask_output   = y
write_nvalid_output = y
use_gpu             = n
CFGEOF

    badchan7_split_cfg="$OUT_DIR/badchan7_split.cfg"
    badchan7_split_log="$OUT_DIR/badchan7_split.log"
    rm -f "$OUT_DIR"/badchan7_split.*.FITS
    cat > "$badchan7_split_cfg" <<CFGEOF
path                = ${DATA_DIR}/
infileQ             = TEST_SPLIT_LO.Q.FITSCUBE,TEST_SPLIT_HI.Q.FITSCUBE
infileU             = TEST_SPLIT_LO.U.FITSCUBE,TEST_SPLIT_HI.U.FITSCUBE
outfile             = ${OUT_DIR}/badchan7_split
remove_badchan      = y
global_badchan_file = ${badchan7_none_list},${badchan7_band2_list}
subim               = n
rem_mean            = 0
remove_qu_bias      = n
resiQ               = 0.0,0.0
slopeQ              = 0.0,0.0
resiU               = 0.0,0.0
slopeU              = 0.0,0.0
ofac                = 1
fac                 = 3.14159265358979
use_auto_rm_range   = 0
beg_rm              = -50.0
end_rm              = 50.0
nrm                 = 201
output_mode         = ap
ap_angle_mode       = phase
write_mask_output   = y
write_nvalid_output = y
use_gpu             = n
CFGEOF

    if run_binary "$BIN_SERIAL" "$badchan7_ref_cfg" "$badchan7_ref_log" && \
       run_binary "$BIN_SERIAL" "$badchan7_split_cfg" "$badchan7_split_log"; then
        all_match=1
        for suffix in AMP.RMCUBE PHA.RMCUBE MASK.CUBE NVALID.MAP; do
            f1="$OUT_DIR/badchan7_ref.${suffix}.FITS"
            f2="$OUT_DIR/badchan7_split.${suffix}.FITS"
            if [[ -f "$f1" && -f "$f2" ]]; then
                if ! python3 "$TESTS_DIR/compare_cubes.py" "$f1" "$f2" --exact \
                        > /dev/null 2>&1; then
                    all_match=0
                    fail "Per-band bad-channel files (T7): ${suffix} differs between flagged-undivided and flagged-split runs"
                fi
            else
                all_match=0
                fail "Per-band bad-channel files (T7): ${suffix} output missing (expected $f1 and $f2)"
            fi
        done
        if grep -q "Number of Bad Channels (band" "$badchan7_split_log"; then
            pass "Per-band bad-channel files (T7): per-band badchan_file read path confirmed taken"
        else
            all_match=0
            fail "Per-band bad-channel files (T7): expected per-band bad-channel log line not found (see $badchan7_split_log)"
        fi
        if [[ "$all_match" -eq 1 ]]; then
            pass "Per-band bad-channel files (T7): flagged split-band run bit-identical to flagged undivided cube"
        fi
    else
        fail "Per-band bad-channel files (T7): one or more runs failed (see $badchan7_ref_log / $badchan7_split_log)"
    fi
else
    skip "Serial binary not available; skipping per-band bad-channel files test"
fi

# ---------------------------------------------------------------------------
# 20. GPU offload for multi-band (T8, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md):
#     the same split-band 2-band config as T5/section 17, run through the
#     GPU binary instead of serial. Mirrors the existing single-band GPU
#     tests exactly: non-staged run checked via rtol=2e-3 against the CPU
#     reference (test 7's own "ffast-math vs IEEE" tolerance, not a
#     multi-band-specific relaxation) plus RM-peak-position validation;
#     staged run (gpu_vram_mib=1 forces VRAM sub-block staging, the one
#     path never previously exercised with a multi-band tile) checked
#     bit-identical against the non-staged GPU run (test 9's own pattern).
#     PHA.RMCUBE is deliberately not rtol-checked here, matching every
#     existing GPU test -- phase angle near low-amplitude bins is known to
#     amplify small ffast-math differences well past 2e-3, in the
#     single-band case too (confirmed by direct comparison during this
#     ticket's investigation), so it was never part of the GPU tolerance
#     gate to begin with.
# ---------------------------------------------------------------------------
section "20. GPU offload for multi-band – split-band run via GPU (T8)"

if [[ "$BUILD_GPU" -eq 1 && -x "$BIN_GPU" && -f "${OUT_DIR}/split_identity.AMP.RMCUBE.FITS" ]]; then
    mbgpu_cfg="$OUT_DIR/mb_gpu.cfg"
    mbgpu_log="$OUT_DIR/mb_gpu.log"
    rm -f "$OUT_DIR"/mb_gpu.*.FITS
    cat > "$mbgpu_cfg" <<CFGEOF
path                = ${DATA_DIR}/
infileQ             = TEST_SPLIT_LO.Q.FITSCUBE,TEST_SPLIT_HI.Q.FITSCUBE
infileU             = TEST_SPLIT_LO.U.FITSCUBE,TEST_SPLIT_HI.U.FITSCUBE
outfile             = ${OUT_DIR}/mb_gpu
remove_badchan      = n
global_badchan_file = /dev/null,/dev/null
subim               = n
rem_mean            = 0
remove_qu_bias      = n
resiQ               = 0.0,0.0
slopeQ              = 0.0,0.0
resiU               = 0.0,0.0
slopeU              = 0.0,0.0
ofac                = 1
fac                 = 3.14159265358979
use_auto_rm_range   = 0
beg_rm              = -50.0
end_rm              = 50.0
nrm                 = 201
output_mode         = ap
ap_angle_mode       = phase
write_mask_output   = y
write_nvalid_output = y
use_gpu             = y
CFGEOF
    if run_binary "$BIN_GPU" "$mbgpu_cfg" "$mbgpu_log"; then
        mbgpu_amp="$OUT_DIR/mb_gpu.AMP.RMCUBE.FITS"
        if [[ -f "$mbgpu_amp" ]]; then
            if python3 "$TESTS_DIR/check_rm_peak.py" "$mbgpu_amp" "$TRUTH" > /dev/null 2>&1; then
                pass "Multi-band GPU (T8): RM peaks at correct positions"
            else
                fail "Multi-band GPU (T8): RM peak check failed (see check_rm_peak.py output)"
            fi
            if python3 "$TESTS_DIR/compare_cubes.py" \
                    "$OUT_DIR/split_identity.AMP.RMCUBE.FITS" "$mbgpu_amp" \
                    --rtol 2e-3 > /dev/null 2>&1; then
                pass "Multi-band GPU AMP (T8): matches CPU reference within rtol=2e-3 (ffast-math vs IEEE)"
            else
                fail "Multi-band GPU AMP (T8): differs from CPU reference beyond rtol=2e-3"
            fi
        else
            fail "Multi-band GPU (T8): AMP output cube not found: $mbgpu_amp"
        fi
    else
        fail "Multi-band GPU (T8): run failed (see $mbgpu_log)"
    fi

    # Staging: force VRAM sub-block subdivision (gpu_vram_mib=1), same
    # pattern as test 9. OMP_TARGET_OFFLOAD=DISABLED for determinism when
    # comparing bit-for-bit against the non-staged run above -- both are
    # the same -ffast-math-compiled kernel, so host-fallback vs real-device
    # dispatch of that same compiled code is expected to be bit-identical
    # (no cross-thread reduction in this kernel), exactly as test 9 already
    # relies on for the single-band case.
    export OMP_TARGET_OFFLOAD="${OMP_TARGET_OFFLOAD:-DISABLED}"
    mbgpu_stg_cfg="$OUT_DIR/mb_gpu_stage.cfg"
    mbgpu_stg_log="$OUT_DIR/mb_gpu_stage.log"
    rm -f "$OUT_DIR"/mb_gpu_stage.*.FITS
    { cat "$mbgpu_cfg"; echo "gpu_vram_mib=1"; echo "mem_frac_vram=0.10"; } | \
        sed "s|outfile             = ${OUT_DIR}/mb_gpu\$|outfile             = ${OUT_DIR}/mb_gpu_stage|" \
        > "$mbgpu_stg_cfg"
    if run_binary "$BIN_GPU" "$mbgpu_stg_cfg" "$mbgpu_stg_log"; then
        if grep -q "Staging sub-blocks:  T" "$mbgpu_stg_log"; then
            mbgpu_stg_amp="$OUT_DIR/mb_gpu_stage.AMP.RMCUBE.FITS"
            if [[ -f "$mbgpu_stg_amp" ]]; then
                if python3 "$TESTS_DIR/compare_cubes.py" \
                        "$OUT_DIR/mb_gpu.AMP.RMCUBE.FITS" "$mbgpu_stg_amp" \
                        --exact > /dev/null 2>&1; then
                    pass "Multi-band GPU staging AMP (T8): bit-identical to non-staged multi-band GPU"
                else
                    fail "Multi-band GPU staging AMP (T8): differs from non-staged multi-band GPU (gather/scatter bug?)"
                fi
            else
                fail "Multi-band GPU staging (T8): AMP output cube not found: $mbgpu_stg_amp"
            fi
        else
            fail "Multi-band GPU staging (T8): staging path was NOT activated (check planner logic)"
        fi
    else
        fail "Multi-band GPU staging (T8): run failed (see $mbgpu_stg_log)"
    fi
else
    skip "GPU binary or multi-band CPU reference not available; skipping multi-band GPU test"
fi

# ---------------------------------------------------------------------------
# 21. io_overlap for multi-band (T9, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md):
#     the T5 split-band 2-band config, forced into 7 tiles (32 Dec rows / 5
#     per tile, uneven remainder -- same shape as test 13's single-band
#     case) with io_overlap=y. Must be bit-identical to the existing
#     single-tile split_identity reference (test 17/T5) -- the ping-pong
#     output buffers io_overlap double-buffers are sized from the already-
#     merged nz_out and never touch how many bands contributed to it, so
#     this exercises exactly the "genuinely does not care about band
#     count" claim from the T9 investigation, not a new code path.
# ---------------------------------------------------------------------------
section "21. io_overlap for multi-band – multi-tile async write == single-tile reference (T9)"

if [[ -x "$BIN_OMP" && -f "${OUT_DIR}/split_identity.AMP.RMCUBE.FITS" ]]; then
    mbovl_cfg="$OUT_DIR/mb_io_overlap.cfg"
    mbovl_log="$OUT_DIR/mb_io_overlap.log"
    rm -f "$OUT_DIR"/mb_io_overlap.*.FITS
    cat > "$mbovl_cfg" <<CFGEOF
path                = ${DATA_DIR}/
infileQ             = TEST_SPLIT_LO.Q.FITSCUBE,TEST_SPLIT_HI.Q.FITSCUBE
infileU             = TEST_SPLIT_LO.U.FITSCUBE,TEST_SPLIT_HI.U.FITSCUBE
outfile             = ${OUT_DIR}/mb_io_overlap
remove_badchan      = n
global_badchan_file = /dev/null,/dev/null
subim               = n
rem_mean            = 0
remove_qu_bias      = n
resiQ               = 0.0,0.0
slopeQ              = 0.0,0.0
resiU               = 0.0,0.0
slopeU              = 0.0,0.0
ofac                = 1
fac                 = 3.14159265358979
use_auto_rm_range   = 0
beg_rm              = -50.0
end_rm              = 50.0
nrm                 = 201
output_mode         = ap
ap_angle_mode       = phase
write_mask_output   = y
write_nvalid_output = y
use_gpu             = n
tile_auto           = n
tile_ra             = 32
tile_dec            = 5
io_overlap          = y
log_level           = debug
CFGEOF
    if run_binary "$BIN_OMP" "$mbovl_cfg" "$mbovl_log"; then
        n_tiles_mb=$(grep -c "Doing tile" "$mbovl_log" || true)
        if [[ "$n_tiles_mb" -gt 1 ]]; then
            pass "io_overlap multi-band (T9): multi-tile run confirmed (${n_tiles_mb} tiles)"
        else
            fail "io_overlap multi-band (T9): expected >1 tile, got $n_tiles_mb"
        fi
        all_match=1
        for suffix in AMP.RMCUBE PHA.RMCUBE MASK.CUBE NVALID.MAP; do
            f1="$OUT_DIR/split_identity.${suffix}.FITS"
            f2="$OUT_DIR/mb_io_overlap.${suffix}.FITS"
            if [[ -f "$f1" && -f "$f2" ]]; then
                if ! python3 "$TESTS_DIR/compare_cubes.py" "$f1" "$f2" --exact \
                        > /dev/null 2>&1; then
                    all_match=0
                    fail "io_overlap multi-band (T9): ${suffix} differs from single-tile reference"
                fi
            else
                all_match=0
                fail "io_overlap multi-band (T9): ${suffix} output missing (expected $f1 and $f2)"
            fi
        done
        if [[ "$all_match" -eq 1 ]]; then
            pass "io_overlap multi-band (T9): bit-identical to single-tile split-band reference"
        fi
        if require_no_overlapping_tile_writes "$mbovl_log" > /dev/null 2>&1; then
            pass "io_overlap multi-band (T9): tile writes never overlap (single-handle serialization holds)"
        else
            fail "io_overlap multi-band (T9): two tile writes overlapped in time (see $mbovl_log)"
        fi
    else
        fail "io_overlap multi-band (T9): run failed (see $mbovl_log)"
    fi
else
    skip "OMP binary or multi-band CPU reference not available; skipping io_overlap multi-band test"
fi

# ---------------------------------------------------------------------------
# 22. io_read_threads>1 for multi-band (T9): same split-band 2-band config,
#     forced into 7 tiles, with io_read_threads=4. Must be bit-identical to
#     the same single-tile split_identity reference -- the parallel
#     channel-split read only ever touches the reference band's own
#     channel range/buffer offset; the other bands' reads happen
#     afterward, sequentially, outside that parallel region entirely.
# ---------------------------------------------------------------------------
section "22. io_read_threads>1 for multi-band – parallel read == single-tile reference (T9)"

if [[ -x "$BIN_OMP" && -f "${OUT_DIR}/split_identity.AMP.RMCUBE.FITS" ]]; then
    mbrt_cfg="$OUT_DIR/mb_io_read_threads.cfg"
    mbrt_log="$OUT_DIR/mb_io_read_threads.log"
    rm -f "$OUT_DIR"/mb_io_read_threads.*.FITS
    sed -e "s|outfile             = ${OUT_DIR}/mb_io_overlap\$|outfile             = ${OUT_DIR}/mb_io_read_threads|" \
        -e "s|io_overlap          = y|io_overlap          = n\nio_read_threads     = 4|" \
        "$mbovl_cfg" > "$mbrt_cfg"
    if run_binary "$BIN_OMP" "$mbrt_cfg" "$mbrt_log"; then
        if grep -qE "Parallel FITS IO: opened +4 +handles/file" "$mbrt_log"; then
            pass "io_read_threads>1 multi-band (T9): parallel-read path confirmed taken (4 handles)"
        else
            fail "io_read_threads>1 multi-band (T9): expected 4-handle parallel-read log line not found (see $mbrt_log)"
        fi
        all_match=1
        for suffix in AMP.RMCUBE PHA.RMCUBE MASK.CUBE NVALID.MAP; do
            f1="$OUT_DIR/split_identity.${suffix}.FITS"
            f2="$OUT_DIR/mb_io_read_threads.${suffix}.FITS"
            if [[ -f "$f1" && -f "$f2" ]]; then
                if ! python3 "$TESTS_DIR/compare_cubes.py" "$f1" "$f2" --exact \
                        > /dev/null 2>&1; then
                    all_match=0
                    fail "io_read_threads>1 multi-band (T9): ${suffix} differs from single-tile reference"
                fi
            else
                all_match=0
                fail "io_read_threads>1 multi-band (T9): ${suffix} output missing (expected $f1 and $f2)"
            fi
        done
        if [[ "$all_match" -eq 1 ]]; then
            pass "io_read_threads>1 multi-band (T9): bit-identical to single-tile split-band reference"
        fi
    else
        fail "io_read_threads>1 multi-band (T9): run failed (see $mbrt_log)"
    fi
else
    skip "OMP binary or multi-band CPU reference not available; skipping io_read_threads>1 multi-band test"
fi

# ---------------------------------------------------------------------------
# 23. RM-CLEAN thesis scenario (T1) – Chapter 6 Figures 6.1-6.3
# ---------------------------------------------------------------------------
section "23. RM-CLEAN thesis scenario – Chapter 6 Figures 6.1/6.2/6.3 (T1)"

RMCLEAN_BUILD_DIR="$OUT_DIR/rmclean_build"
RMCLEAN_CSV_DIR="$OUT_DIR/rmclean_csv"
RMCLEAN_PLOT_DIR="$OUT_DIR/rmclean_plots"
mkdir -p "$RMCLEAN_BUILD_DIR" "$RMCLEAN_CSV_DIR" "$RMCLEAN_PLOT_DIR"
rmclean_o="$RMCLEAN_BUILD_DIR/rmclean.o"
rmclean_mod_dir="$RMCLEAN_BUILD_DIR"
rmclean_bin="$RMCLEAN_BUILD_DIR/thesis_scenario_rmclean"
rmclean_log="$OUT_DIR/thesis_scenario_rmclean.log"
RMCLEAN_VENV_PY="$HOME/venv/rmtool/bin/python3"

if gfortran -cpp -std=gnu -fallow-argument-mismatch -ffree-line-length-none \
        -O3 -fopenmp -J"$rmclean_mod_dir" -c "$REPO_ROOT/src/rmclean.f90" \
        -o "$rmclean_o" 2>"$OUT_DIR/rmclean_mod_build.log" \
    && gfortran -cpp -std=gnu -fallow-argument-mismatch -ffree-line-length-none \
        -O3 -fopenmp -I"$rmclean_mod_dir" -J"$rmclean_mod_dir" \
        "$TESTS_DIR/thesis_scenario_rmclean.f90" "$rmclean_o" \
        -o "$rmclean_bin" -lfftw3 2>>"$OUT_DIR/rmclean_mod_build.log"; then
    if "$rmclean_bin" "$RMCLEAN_CSV_DIR" > "$rmclean_log" 2>&1; then
        while IFS= read -r line; do
            case "$line" in
                *"[PASS]"*) pass "${line#*\[PASS\] }" ;;
                *"[FAIL]"*) fail "${line#*\[FAIL\] }" ;;
            esac
        done < "$rmclean_log"
        if [[ -x "$RMCLEAN_VENV_PY" ]]; then
            if "$RMCLEAN_VENV_PY" "$TESTS_DIR/plot_thesis_scenario_rmclean.py" \
                    "$RMCLEAN_CSV_DIR" "$RMCLEAN_PLOT_DIR" \
                    > "$OUT_DIR/rmclean_plot.log" 2>&1; then
                pass "RM-CLEAN thesis scenario: plots written to $RMCLEAN_PLOT_DIR"
            else
                fail "RM-CLEAN thesis scenario: plotting failed (see $OUT_DIR/rmclean_plot.log)"
            fi
        else
            skip "~/venv/rmtool python3 not found; skipping RM-CLEAN thesis scenario plots"
        fi
    else
        fail "RM-CLEAN thesis scenario: program exited non-zero (see $rmclean_log)"
    fi
else
    skip "FFTW3 (or gfortran) not available for rmclean_mod; skipping RM-CLEAN thesis scenario test (see $OUT_DIR/rmclean_mod_build.log)"
fi

# ---------------------------------------------------------------------------
# 24. RM-CLEAN lsq_ref flexibility + derotate_to_lsq_ref
# ---------------------------------------------------------------------------
section "24. RM-CLEAN lsq_ref flexibility + derotate_to_lsq_ref"

lsqref_flex_bin="$RMCLEAN_BUILD_DIR/test_rmclean_lsqref_flex"
lsqref_flex_log="$OUT_DIR/test_rmclean_lsqref_flex.log"

if [[ -f "$rmclean_o" ]]; then
    if gfortran -cpp -std=gnu -fallow-argument-mismatch -ffree-line-length-none \
            -O3 -fopenmp -I"$rmclean_mod_dir" -J"$rmclean_mod_dir" \
            "$TESTS_DIR/test_rmclean_lsqref_flex.f90" "$rmclean_o" \
            -o "$lsqref_flex_bin" -lfftw3 2>"$OUT_DIR/rmclean_lsqref_flex_build.log"; then
        if "$lsqref_flex_bin" > "$lsqref_flex_log" 2>&1; then
            while IFS= read -r line; do
                case "$line" in
                    *"[PASS]"*) pass "${line#*\[PASS\] }" ;;
                    *"[FAIL]"*) fail "${line#*\[FAIL\] }" ;;
                esac
            done < "$lsqref_flex_log"
        else
            fail "RM-CLEAN lsq_ref flexibility: program exited non-zero (see $lsqref_flex_log)"
        fi
    else
        fail "RM-CLEAN lsq_ref flexibility: build failed (see $OUT_DIR/rmclean_lsqref_flex_build.log)"
    fi
else
    skip "rmclean_mod.o not built (section 23 skipped); skipping RM-CLEAN lsq_ref flexibility test"
fi

# ---------------------------------------------------------------------------
# 25. RM-CLEAN get_lsq_ref_compute
# ---------------------------------------------------------------------------
section "25. RM-CLEAN get_lsq_ref_compute"

optimal_lsqref_bin="$RMCLEAN_BUILD_DIR/test_optimal_lsq_ref"
optimal_lsqref_log="$OUT_DIR/test_optimal_lsq_ref.log"

if [[ -f "$rmclean_o" ]]; then
    if gfortran -cpp -std=gnu -fallow-argument-mismatch -ffree-line-length-none \
            -O3 -fopenmp -I"$rmclean_mod_dir" -J"$rmclean_mod_dir" \
            "$TESTS_DIR/test_optimal_lsq_ref.f90" "$rmclean_o" \
            -o "$optimal_lsqref_bin" -lfftw3 2>"$OUT_DIR/rmclean_optimal_lsqref_build.log"; then
        if "$optimal_lsqref_bin" > "$optimal_lsqref_log" 2>&1; then
            while IFS= read -r line; do
                case "$line" in
                    *"[PASS]"*) pass "${line#*\[PASS\] }" ;;
                    *"[FAIL]"*) fail "${line#*\[FAIL\] }" ;;
                esac
            done < "$optimal_lsqref_log"
        else
            fail "RM-CLEAN get_lsq_ref_compute: program exited non-zero (see $optimal_lsqref_log)"
        fi
    else
        fail "RM-CLEAN get_lsq_ref_compute: build failed (see $OUT_DIR/rmclean_optimal_lsqref_build.log)"
    fi
else
    skip "rmclean_mod.o not built (section 23 skipped); skipping RM-CLEAN get_lsq_ref_compute test"
fi

# ---------------------------------------------------------------------------
# 26. RM-CLEAN get_drm's oversample floor -- why it's 2, not bare Nyquist
# ---------------------------------------------------------------------------
section "26. RM-CLEAN get_drm oversample floor (why 2, not bare Nyquist)"

drmfloor_bin="$RMCLEAN_BUILD_DIR/test_drm_floor"
drmfloor_log="$OUT_DIR/test_drm_floor.log"

if [[ -f "$rmclean_o" ]]; then
    if gfortran -cpp -std=gnu -fallow-argument-mismatch -ffree-line-length-none \
            -O3 -fopenmp -I"$rmclean_mod_dir" -J"$rmclean_mod_dir" \
            "$TESTS_DIR/test_drm_floor.f90" "$rmclean_o" \
            -o "$drmfloor_bin" -lfftw3 2>"$OUT_DIR/rmclean_drmfloor_build.log"; then
        if "$drmfloor_bin" > "$drmfloor_log" 2>&1; then
            while IFS= read -r line; do
                case "$line" in
                    *"[PASS]"*) pass "${line#*\[PASS\] }" ;;
                    *"[FAIL]"*) fail "${line#*\[FAIL\] }" ;;
                esac
            done < "$drmfloor_log"
        else
            fail "RM-CLEAN get_drm oversample floor: program exited non-zero (see $drmfloor_log)"
        fi
    else
        fail "RM-CLEAN get_drm oversample floor: build failed (see $OUT_DIR/rmclean_drmfloor_build.log)"
    fi
else
    skip "rmclean_mod.o not built (section 23 skipped); skipping RM-CLEAN get_drm oversample floor test"
fi

# ---------------------------------------------------------------------------
# 27. RM-CLEAN get_drm actually REFUSES an unsafe oversample -- inverted
#     pass logic: this program is EXPECTED to be terminated by get_drm's
#     own stop(1), so a NONZERO exit code is the correct outcome here.
# ---------------------------------------------------------------------------
section "27. RM-CLEAN get_drm refuses unsafe oversample (enforcement)"

drmenforce_bin="$RMCLEAN_BUILD_DIR/test_drm_floor_enforcement"
drmenforce_log="$OUT_DIR/test_drm_floor_enforcement.log"

if [[ -f "$rmclean_o" ]]; then
    if gfortran -cpp -std=gnu -fallow-argument-mismatch -ffree-line-length-none \
            -O3 -fopenmp -I"$rmclean_mod_dir" -J"$rmclean_mod_dir" \
            "$TESTS_DIR/test_drm_floor_enforcement.f90" "$rmclean_o" \
            -o "$drmenforce_bin" -lfftw3 2>"$OUT_DIR/rmclean_drmenforce_build.log"; then
        if "$drmenforce_bin" > "$drmenforce_log" 2>&1; then
            fail "get_drm: incorrectly ALLOWED an unsafe oversample to proceed (see $drmenforce_log)"
        else
            pass "get_drm: correctly refused an unsafe oversample (<2) with a nonzero exit"
        fi
    else
        fail "RM-CLEAN get_drm enforcement: build failed (see $OUT_DIR/rmclean_drmenforce_build.log)"
    fi
else
    skip "rmclean_mod.o not built (section 23 skipped); skipping RM-CLEAN get_drm enforcement test"
fi

# ---------------------------------------------------------------------------
# 28. RM-CLEAN matched-filter peak refinement (docs/dev/
#     RMCLEAN_INTEGRATION_PLAN.md T3, T3-Phase 1): validates
#     refine_peak_matched_filter/rmsf_point_direct against
#     peak_interp_parabolic across two independent point-source
#     scenarios and 4 sampling densities, including deliberately coarse
#     (resolution-only, get_drm-refusing) grids where peak_interp_
#     parabolic is known to fail. clean_complex itself is NOT yet
#     modified at this point (T3-Phase 2) -- this section only proves
#     the new subroutines behave correctly in isolation.
# ---------------------------------------------------------------------------
section "28. RM-CLEAN matched-filter peak refinement (T3)"

mfr_bin="$RMCLEAN_BUILD_DIR/test_matched_filter_refine"
mfr_log="$OUT_DIR/test_matched_filter_refine.log"

if [[ -f "$rmclean_o" ]]; then
    if gfortran -cpp -std=gnu -fallow-argument-mismatch -ffree-line-length-none \
            -O3 -fopenmp -I"$rmclean_mod_dir" -J"$rmclean_mod_dir" \
            "$TESTS_DIR/test_matched_filter_refine.f90" "$rmclean_o" \
            -o "$mfr_bin" -lfftw3 2>"$OUT_DIR/rmclean_mfr_build.log"; then
        if "$mfr_bin" > "$mfr_log" 2>&1; then
            while IFS= read -r line; do
                if [[ "$line" == *"[PASS]"* || "$line" == *"[FAIL]"* ]]; then
                    echo "  $line"
                fi
            done < "$mfr_log"
            if grep -q "^\[PASS\] test_matched_filter_refine" "$mfr_log"; then
                pass "RM-CLEAN matched-filter peak refinement: all checks passed"
            else
                fail "RM-CLEAN matched-filter peak refinement: one or more checks failed (see $mfr_log)"
            fi
        else
            fail "RM-CLEAN matched-filter peak refinement: program exited non-zero (see $mfr_log)"
        fi
    else
        fail "RM-CLEAN matched-filter peak refinement: build failed (see $OUT_DIR/rmclean_mfr_build.log)"
    fi
else
    skip "rmclean_mod.o not built (section 23 skipped); skipping RM-CLEAN matched-filter peak refinement test"
fi

# ---------------------------------------------------------------------------
# 29. rmclean_cubes: standalone RM-CLEAN tool driven against a REAL
#     rm_synthesis dirty AMP/PHA/MASK cube (docs/dev/
#     RMCLEAN_INTEGRATION_PLAN.md T2) -- built via its own Makefile
#     target (`make rmclean_cubes`), not the ad hoc gfortran calls
#     sections 22-28 use for rmclean_mod's own unit tests. Also exercises
#     rm_synthesis's own lsq_ref_mode option (this section's own T2
#     follow-up): a cube built with lsq_ref_mode=mid must carry a
#     nonzero LSQREF header rmclean_cubes reads back correctly.
# ---------------------------------------------------------------------------
section "29. rmclean_cubes: end-to-end against a real dirty cube (T2)"

if make rmclean_cubes > "$OUT_DIR/rmclean_cubes_build.log" 2>&1; then
    pass "rmclean_cubes: build succeeded"

    rmc_lsqref_cfg=$(make_cfg "rmc_lsqref" "0" "lsq_ref_mode = mid")
    rmc_lsqref_log="$OUT_DIR/rmc_lsqref_rmsynth.log"
    rmc_lsqref_amp="$OUT_DIR/rmc_lsqref.AMP.RMCUBE.FITS"
    if run_binary "$BIN_SERIAL" "$rmc_lsqref_cfg" "$rmc_lsqref_log"; then
        lsqref_val=$(python3 -c "
from astropy.io import fits
print(fits.getheader('$rmc_lsqref_amp').get('LSQREF'))
")
        if python3 -c "
import sys
v = float('$lsqref_val') if '$lsqref_val' != 'None' else None
sys.exit(0 if v is not None and abs(v) > 1e-6 else 1)
"; then
            pass "rm_synthesis lsq_ref_mode=mid: LSQREF header written and nonzero ($lsqref_val)"
        else
            fail "rm_synthesis lsq_ref_mode=mid: LSQREF header missing or zero (got '$lsqref_val')"
        fi

        rmc_out="$OUT_DIR/rmc_cleaned"
        rmc_log="$OUT_DIR/rmclean_cubes_run.log"
        if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" \
                phafile="$OUT_DIR/rmc_lsqref.PHA.RMCUBE.FITS" \
                maskfile="$OUT_DIR/rmc_lsqref.MASK.CUBE.FITS" \
                outfile="$rmc_out" abs_flux_floor=0.01 niter=200 gain=0.1 \
                > "$rmc_log" 2>&1; then
            pass "rmclean_cubes: ran to completion on a real lsq_ref_mode=mid cube"
            if grep -q "^Gate 0 OK" "$rmc_log"; then
                pass "rmclean_cubes: Gate 0 validated the existing RM grid"
            else
                fail "rmclean_cubes: Gate 0 pass marker missing (see $rmc_log)"
            fi
            if python3 "$TESTS_DIR/check_rm_peak.py" \
                    "${rmc_out}.RESTORED.AMP.RMCUBE.FITS" "$TRUTH"; then
                pass "rmclean_cubes: RESTORED.AMP recovers both known point sources' RM"
            else
                fail "rmclean_cubes: RESTORED.AMP peak RM check failed"
            fi
        else
            fail "rmclean_cubes: run failed (see $rmc_log)"
        fi

        # T4a: memory-budgeted tiling -- forcing many small, non-full-width
        # tiles (tile_ra/tile_dec well below the image size, so both the
        # RA-subdivided and multi-row-per-tile paths get exercised) must
        # agree with the default (auto, single whole-image tile) run.
        # NOT byte-identical -- see check_tile_consistency.py's own module
        # docstring for why (a real, pre-existing -O3 -march=native
        # floating-point reassociation effect in clean_complex's own
        # tiered-refinement threshold decision, confirmed via -O0 rebuild;
        # not a tiling logic bug).
        rmc_tiled="$OUT_DIR/rmc_cleaned_tiled"
        rmc_tiled_log="$OUT_DIR/rmclean_cubes_tiled_run.log"
        if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" \
                phafile="$OUT_DIR/rmc_lsqref.PHA.RMCUBE.FITS" \
                maskfile="$OUT_DIR/rmc_lsqref.MASK.CUBE.FITS" \
                outfile="$rmc_tiled" abs_flux_floor=0.01 niter=200 gain=0.1 \
                tile_auto=n tile_ra=3 tile_dec=5 \
                > "$rmc_tiled_log" 2>&1; then
            pass "rmclean_cubes: ran to completion with forced small (3x5) tiles"
            if grep -q "Tile plan: tile_ra x tile_dec = 3 x 5" "$rmc_tiled_log"; then
                pass "rmclean_cubes: forced tile geometry took effect"
            else
                fail "rmclean_cubes: forced tile_ra/tile_dec not reflected in tile plan (see $rmc_tiled_log)"
            fi
            if python3 "$TESTS_DIR/check_tile_consistency.py" "$rmc_out" "$rmc_tiled"; then
                pass "rmclean_cubes: small-tile output numerically agrees with default (single-tile) output"
            else
                fail "rmclean_cubes: small-tile output disagrees with default output beyond tolerance"
            fi
            if python3 "$TESTS_DIR/check_rm_peak.py" \
                    "${rmc_tiled}.RESTORED.AMP.RMCUBE.FITS" "$TRUTH"; then
                pass "rmclean_cubes: small-tile RESTORED.AMP recovers both known point sources' RM"
            else
                fail "rmclean_cubes: small-tile RESTORED.AMP peak RM check failed"
            fi
        else
            fail "rmclean_cubes: forced small-tile run failed (see $rmc_tiled_log)"
        fi

        # T4b: io_read_threads -- parallel readonly chunked tile reads must
        # be BYTE-IDENTICAL to the default (default tile geometry is
        # unchanged here, only the read parallelism, so none of T4a's own
        # floating-point-reassociation caveat applies).
        rt_ok=1
        for rt in 1 2 4; do
            rmc_rt="$OUT_DIR/rmc_rt${rt}"
            rmc_rt_log="$OUT_DIR/rmclean_cubes_rt${rt}_run.log"
            if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" \
                    phafile="$OUT_DIR/rmc_lsqref.PHA.RMCUBE.FITS" \
                    maskfile="$OUT_DIR/rmc_lsqref.MASK.CUBE.FITS" \
                    outfile="$rmc_rt" abs_flux_floor=0.01 niter=200 gain=0.1 \
                    io_read_threads=$rt > "$rmc_rt_log" 2>&1; then
                for suffix in CLEAN.AMP CLEAN.PHA RESID.AMP RESID.PHA \
                              RESTORED.AMP RESTORED.PHA; do
                    if ! cmp -s "${rmc_out}.${suffix}.RMCUBE.FITS" \
                                "${rmc_rt}.${suffix}.RMCUBE.FITS"; then
                        rt_ok=0
                        fail "rmclean_cubes: io_read_threads=$rt differs from io_read_threads=1 on $suffix"
                    fi
                done
            else
                rt_ok=0
                fail "rmclean_cubes: io_read_threads=$rt run failed (see $rmc_rt_log)"
            fi
        done
        [[ "$rt_ok" -eq 1 ]] && \
            pass "rmclean_cubes: io_read_threads=1,2,4 all bit-identical to the default run"

        # T4c: nwriters -- raw stream writes bypassing CFITSIO must
        # be BYTE-IDENTICAL to the default (FTPSSE) write path. Repeated 5x
        # at nwriters=4: this path previously had a genuine,
        # PROBABILISTIC bug (FTGHAD's 3 output arguments come back with
        # their upper 32 bits UNTOUCHED by this system's installed
        # libcfitsio -- a real Fortran-wrapper/library ABI truncation,
        # confirmed with a minimal standalone reproducer outside this
        # codebase entirely -- so an uninitialized receiving variable's
        # leftover stack garbage in the upper bits produced an
        # intermittent, wildly wrong byte offset roughly half the time).
        # Fixed by zero-initializing those 3 variables immediately before
        # every FTGHAD call (open_output_cube's own comment has the full
        # story); a single run would not reliably have caught the
        # original bug, so this repeats the check rather than trusting
        # one pass.
        wt_ok=1
        for wt in 1 2 4; do
            for rep in 1 2 3 4 5; do
                rmc_wt="$OUT_DIR/rmc_wt${wt}_${rep}"
                rmc_wt_log="$OUT_DIR/rmclean_cubes_wt${wt}_${rep}_run.log"
                if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" \
                        phafile="$OUT_DIR/rmc_lsqref.PHA.RMCUBE.FITS" \
                        maskfile="$OUT_DIR/rmc_lsqref.MASK.CUBE.FITS" \
                        outfile="$rmc_wt" abs_flux_floor=0.01 niter=200 gain=0.1 \
                        nwriters=$wt > "$rmc_wt_log" 2>&1; then
                    for suffix in CLEAN.AMP CLEAN.PHA RESID.AMP RESID.PHA \
                                  RESTORED.AMP RESTORED.PHA; do
                        if ! cmp -s "${rmc_out}.${suffix}.RMCUBE.FITS" \
                                    "${rmc_wt}.${suffix}.RMCUBE.FITS"; then
                            wt_ok=0
                            fail "rmclean_cubes: nwriters=$wt (rep $rep) differs from nwriters=1 on $suffix"
                        fi
                    done
                else
                    wt_ok=0
                    fail "rmclean_cubes: nwriters=$wt (rep $rep) run failed (see $rmc_wt_log)"
                fi
            done
        done
        [[ "$wt_ok" -eq 1 ]] && \
            pass "rmclean_cubes: nwriters=1,2,4 all bit-identical to the default run (5 reps each)"

        # T4d: io_overlap -- background-thread tile writes must be
        # BYTE-IDENTICAL to the default (inline write) path, both alone
        # and combined with forced small tiles + io_read_threads +
        # nwriters together (the combined stress case actually
        # exercised, not just each mechanism in isolation).
        rmc_ov_log="$OUT_DIR/rmclean_cubes_ov_run.log"
        ov_ok=1
        for rep in 1 2 3 4 5; do
            rmc_ov="$OUT_DIR/rmc_ov${rep}"
            if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" \
                    phafile="$OUT_DIR/rmc_lsqref.PHA.RMCUBE.FITS" \
                    maskfile="$OUT_DIR/rmc_lsqref.MASK.CUBE.FITS" \
                    outfile="$rmc_ov" abs_flux_floor=0.01 niter=200 gain=0.1 \
                    io_overlap=y > "$rmc_ov_log" 2>&1; then
                for suffix in CLEAN.AMP CLEAN.PHA RESID.AMP RESID.PHA \
                              RESTORED.AMP RESTORED.PHA; do
                    if ! cmp -s "${rmc_out}.${suffix}.RMCUBE.FITS" \
                                "${rmc_ov}.${suffix}.RMCUBE.FITS"; then
                        ov_ok=0
                        fail "rmclean_cubes: io_overlap=y (rep $rep) differs from io_overlap=n on $suffix"
                    fi
                done
            else
                ov_ok=0
                fail "rmclean_cubes: io_overlap=y (rep $rep) run failed (see $rmc_ov_log)"
            fi
        done
        [[ "$ov_ok" -eq 1 ]] && \
            pass "rmclean_cubes: io_overlap=y bit-identical to io_overlap=n (5 reps)"

        rmc_combo_log="$OUT_DIR/rmclean_cubes_combo_run.log"
        combo_ok=1
        for rep in 1 2 3 4 5; do
            rmc_combo="$OUT_DIR/rmc_combo${rep}"
            if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" \
                    phafile="$OUT_DIR/rmc_lsqref.PHA.RMCUBE.FITS" \
                    maskfile="$OUT_DIR/rmc_lsqref.MASK.CUBE.FITS" \
                    outfile="$rmc_combo" abs_flux_floor=0.01 niter=200 gain=0.1 \
                    tile_auto=n tile_ra=3 tile_dec=5 io_read_threads=3 \
                    nwriters=3 io_overlap=y \
                    > "$rmc_combo_log" 2>&1; then
                for suffix in CLEAN.AMP CLEAN.PHA RESID.AMP RESID.PHA \
                              RESTORED.AMP RESTORED.PHA; do
                    if ! cmp -s "${rmc_tiled}.${suffix}.RMCUBE.FITS" \
                                "${rmc_combo}.${suffix}.RMCUBE.FITS"; then
                        combo_ok=0
                        fail "rmclean_cubes: combined small-tile+io_read_threads+nwriters+io_overlap (rep $rep) differs from small-tile-only on $suffix"
                    fi
                done
            else
                combo_ok=0
                fail "rmclean_cubes: combined stress run (rep $rep) failed (see $rmc_combo_log)"
            fi
        done
        [[ "$combo_ok" -eq 1 ]] && \
            pass "rmclean_cubes: combined small-tile+io_read_threads+nwriters+io_overlap bit-identical to small-tile-only (5 reps)"

        # Mask-pattern cache correctness: a deliberately varied mask (many
        # distinct per-pixel valid-channel patterns, forcing both cache
        # hits and the mask_pattern_cache_max overflow fallback) must give
        # BIT-IDENTICAL output whether the cache is generous or disabled
        # entirely (mask_pattern_cache_max=0 -- every pixel a one-off
        # throwaway table) -- the cache is purely a reuse optimization,
        # never allowed to change the result.
        rmc_varied_mask="$OUT_DIR/rmc_varied.MASK.CUBE.FITS"
        python3 - "$OUT_DIR/rmc_lsqref.MASK.CUBE.FITS" "$rmc_varied_mask" <<'PYEOF'
import sys
import numpy as np
from astropy.io import fits

src, dst = sys.argv[1], sys.argv[2]
with fits.open(src) as hdul:
    rng = np.random.default_rng(7)
    data = hdul[0].data
    nchan, ny, nx = data.shape
    for iy in range(ny):
        for ix in range(nx):
            nflip = rng.integers(0, 5)
            if nflip > 0:
                chans = rng.choice(nchan, size=nflip, replace=False)
                data[chans, iy, ix] = 0
    # Preserve the CHANFREQ binary table extension (rmclean_cubes reads
    # it for per-channel L_sq) -- only the primary array's data is
    # perturbed above, in place.
    hdul.writeto(dst, overwrite=True)
PYEOF
        rmc_cache_big="$OUT_DIR/rmc_cache_big"
        rmc_cache_zero="$OUT_DIR/rmc_cache_zero"
        if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" \
                phafile="$OUT_DIR/rmc_lsqref.PHA.RMCUBE.FITS" \
                maskfile="$rmc_varied_mask" outfile="$rmc_cache_big" \
                abs_flux_floor=0.01 niter=200 gain=0.1 mask_pattern_cache_max=4096 \
                > "$OUT_DIR/rmc_cache_big.log" 2>&1 && \
           bin/rmclean_cubes ampfile="$rmc_lsqref_amp" \
                phafile="$OUT_DIR/rmc_lsqref.PHA.RMCUBE.FITS" \
                maskfile="$rmc_varied_mask" outfile="$rmc_cache_zero" \
                abs_flux_floor=0.01 niter=200 gain=0.1 mask_pattern_cache_max=0 \
                > "$OUT_DIR/rmc_cache_zero.log" 2>&1; then
            cache_ok=1
            for suffix in CLEAN.AMP CLEAN.PHA RESID.AMP RESID.PHA \
                          RESTORED.AMP RESTORED.PHA; do
                if ! cmp -s "${rmc_cache_big}.${suffix}.RMCUBE.FITS" \
                            "${rmc_cache_zero}.${suffix}.RMCUBE.FITS"; then
                    cache_ok=0
                    fail "rmclean_cubes: cache vs no-cache differ on $suffix"
                fi
            done
            [[ "$cache_ok" -eq 1 ]] && \
                pass "rmclean_cubes: mask-pattern cache (generous) bit-identical to cache disabled"
        else
            fail "rmclean_cubes: mask-pattern cache comparison run(s) failed"
        fi
    else
        fail "rm_synthesis (lsq_ref_mode=mid, for rmclean_cubes test): run failed (see $rmc_lsqref_log)"
    fi
else
    fail "rmclean_cubes: build failed (see $OUT_DIR/rmclean_cubes_build.log)"
fi

section "30. match_cubes: skip-if-already-matched (planning-doc ticket)"

if make match_cubes > "$OUT_DIR/match_cubes_build.log" 2>&1; then
    pass "match_cubes: build succeeded"

    mc_ref="$OUT_DIR/mc_skip_ref.Q.FITSCUBE"
    mc_shifted="$OUT_DIR/mc_skip_shifted.Q.FITSCUBE"
    # strip_fits_ext (src/match_cubes.f90) strips ANY trailing extension,
    # so the tool's own output name drops .FITSCUBE here.
    mc_ref_reproj="${mc_ref%.FITSCUBE}_REPROJ.FITS"
    mc_shifted_reproj="${mc_shifted%.FITSCUBE}_REPROJ.FITS"
    cp tests/data/TEST.Q.FITSCUBE "$mc_ref"
    python3 - "$mc_ref" "$mc_shifted" <<'PYEOF'
import sys
import shutil
from astropy.io import fits
shutil.copy(sys.argv[1], sys.argv[2])
with fits.open(sys.argv[2], mode="update") as hdul:
    hdul[0].header["CRVAL1"] = float(hdul[0].header["CRVAL1"]) + 0.05
    hdul.flush()
PYEOF

    # --- Positive: reffile=itself already matches -> skip, no output,
    # manifest says SKIPPED ---
    rm -f "$mc_ref_reproj" "$OUT_DIR/mc_manifest_pos.txt"
    mc_pos_log="$OUT_DIR/mc_skip_positive.log"
    if bin/match_cubes stages=reproject footprint_mode=reference reffile="$mc_ref" \
            infiles="$mc_ref" manifest="$OUT_DIR/mc_manifest_pos.txt" \
            > "$mc_pos_log" 2>&1; then
        if [[ -f "$mc_ref_reproj" ]]; then
            fail "match_cubes: positive skip test wrote an output file (should have skipped)"
        elif ! grep -q "^SKIP: $mc_ref " "$mc_pos_log"; then
            fail "match_cubes: positive skip test missing SKIP: message (see $mc_pos_log)"
        elif ! grep -qP "^${mc_ref}\tSKIPPED\t${mc_ref}$" "$OUT_DIR/mc_manifest_pos.txt"; then
            fail "match_cubes: positive skip test manifest missing/wrong SKIPPED line"
        else
            pass "match_cubes: already-matching input skipped (no output, manifest correct)"
        fi
    else
        fail "match_cubes: positive skip test run failed (see $mc_pos_log)"
    fi

    # --- Negative: genuinely offset CRVAL -> processed as normal,
    # manifest says PROCESSED ---
    rm -f "$mc_shifted_reproj" "$OUT_DIR/mc_manifest_neg.txt"
    mc_neg_log="$OUT_DIR/mc_skip_negative.log"
    if bin/match_cubes stages=reproject footprint_mode=reference reffile="$mc_ref" \
            infiles="$mc_shifted" manifest="$OUT_DIR/mc_manifest_neg.txt" \
            > "$mc_neg_log" 2>&1; then
        if [[ ! -s "$mc_shifted_reproj" ]]; then
            fail "match_cubes: negative skip test did not write the expected output"
        elif ! grep -qP "^${mc_shifted}\tPROCESSED\t${mc_shifted_reproj}$" "$OUT_DIR/mc_manifest_neg.txt"; then
            fail "match_cubes: negative skip test manifest missing/wrong PROCESSED line"
        else
            pass "match_cubes: genuinely mismatched input processed as normal (regression guard)"
        fi
    else
        fail "match_cubes: negative skip test run failed (see $mc_neg_log)"
    fi

    # --- Safety: a pre-existing output path always aborts the whole run,
    # regardless of what this run's own skip decision would have been
    # (never silently reused or overwritten -- also the regression guard
    # for the clobber-mode FTINIT bugfix: this is exactly the "run twice
    # without cleanup" scenario that previously silently overwrote) ---
    mc_stale_out_log="$OUT_DIR/mc_skip_stale_output.log"
    rm -f "$OUT_DIR/mc_manifest_stale1.txt"
    if bin/match_cubes stages=reproject footprint_mode=reference reffile="$mc_ref" \
            infiles="$mc_shifted" manifest="$OUT_DIR/mc_manifest_stale1.txt" \
            > "$mc_stale_out_log" 2>&1; then
        fail "match_cubes: rerun with a pre-existing output should have aborted, but succeeded"
    elif grep -q "already exists, refusing to proceed" "$mc_stale_out_log"; then
        pass "match_cubes: pre-existing output path aborts the run (clobber-bug regression guard)"
    else
        fail "match_cubes: rerun aborted but not with the expected message (see $mc_stale_out_log)"
    fi

    # --- Safety: a pre-existing manifest path always aborts the whole
    # run too ---
    rm -f "$mc_shifted_reproj"
    mc_stale_manifest_log="$OUT_DIR/mc_skip_stale_manifest.log"
    if bin/match_cubes stages=reproject footprint_mode=reference reffile="$mc_ref" \
            infiles="$mc_shifted" manifest="$OUT_DIR/mc_manifest_neg.txt" \
            > "$mc_stale_manifest_log" 2>&1; then
        fail "match_cubes: rerun with a pre-existing manifest should have aborted, but succeeded"
    elif grep -q "manifest already exists" "$mc_stale_manifest_log"; then
        pass "match_cubes: pre-existing manifest path aborts the run"
    else
        fail "match_cubes: manifest rerun aborted but not with the expected message (see $mc_stale_manifest_log)"
    fi
else
    fail "match_cubes: build failed (see $OUT_DIR/match_cubes_build.log)"
fi

# ---------------------------------------------------------------------------
# 31. rmclean_cubes: threshold unit conversion + threshold_snr auto noise
#     floor (planning-doc ticket, requested alongside the pipeline suite).
#     Reuses section 29's own rmc_lsqref_amp/PHA/MASK cube trio -- a real
#     dirty cube with known BUNIT, no need to rebuild one here.
# ---------------------------------------------------------------------------
section "31. rmclean_cubes: abs_flux_floor units (Jy/mJy/uJy) + auto_nsigma (T8)"

if [[ -s "$rmc_lsqref_amp" ]]; then
    rmc_pha="$OUT_DIR/rmc_lsqref.PHA.RMCUBE.FITS"
    rmc_mask="$OUT_DIR/rmc_lsqref.MASK.CUBE.FITS"

    # --- Unit conversion: abs_flux_floor=10mJy must resolve to the same
    # native-unit value as a hand-computed conversion from the cube's own
    # BUNIT, not just "some number" -- verified against the actual header,
    # not a hardcoded assumption about what BUNIT is. ---
    rmc_expected_mjy=$(python3 -c "
from astropy.io import fits
bunit = fits.getheader('$rmc_lsqref_amp').get('BUNIT', 'Jy').split('/')[0].strip().lower()
scale = {'jy': 1.0, 'mjy': 1.0e-3, 'ujy': 1.0e-6}.get(bunit, 1.0)
s = f'{10.0 * 1.0e-3 / scale:.6f}'
# Fortran's F0.6 edit descriptor omits the leading zero before the
# decimal point (e.g. '.010000', not '0.010000') -- match that here.
print(s[1:] if s.startswith('0.') else s)
")
    rmc_unit_out="$OUT_DIR/rmc_thresh_unit"
    rmc_unit_log="$OUT_DIR/rmc_thresh_unit.log"
    rm -f "${rmc_unit_out}".*.RMCUBE.FITS
    if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" phafile="$rmc_pha" \
            maskfile="$rmc_mask" outfile="$rmc_unit_out" abs_flux_floor=10mJy \
            niter=200 gain=0.1 > "$rmc_unit_log" 2>&1; then
        if grep -qF -- "-> $rmc_expected_mjy (native AMP-cube units" "$rmc_unit_log"; then
            pass "rmclean_cubes: abs_flux_floor=10mJy converts to the expected native-unit value ($rmc_expected_mjy)"
        else
            fail "rmclean_cubes: abs_flux_floor=10mJy did not convert to the expected native-unit value $rmc_expected_mjy (see $rmc_unit_log)"
        fi
    else
        fail "rmclean_cubes: abs_flux_floor=10mJy run failed (see $rmc_unit_log)"
    fi

    # --- auto_nsigma: runs to completion; per-pixel RM-tail sigma has no
    # randomness at all (unlike the old whole-cube noise_seed pre-scan),
    # so two identical runs must be BIT-IDENTICAL, not just "same printed
    # threshold". ---
    rmc_snr_out1="$OUT_DIR/rmc_thresh_snr1"
    rmc_snr_out2="$OUT_DIR/rmc_thresh_snr2"
    rmc_snr_log1="$OUT_DIR/rmc_thresh_snr1.log"
    rmc_snr_log2="$OUT_DIR/rmc_thresh_snr2.log"
    rm -f "${rmc_snr_out1}".*.RMCUBE.FITS "${rmc_snr_out2}".*.RMCUBE.FITS
    if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" phafile="$rmc_pha" \
            maskfile="$rmc_mask" outfile="$rmc_snr_out1" auto_nsigma=5.0 \
            niter=200 gain=0.1 > "$rmc_snr_log1" 2>&1 && \
       bin/rmclean_cubes ampfile="$rmc_lsqref_amp" phafile="$rmc_pha" \
            maskfile="$rmc_mask" outfile="$rmc_snr_out2" auto_nsigma=5.0 \
            niter=200 gain=0.1 > "$rmc_snr_log2" 2>&1; then
        if [[ -s "${rmc_snr_out1}.RESTORED.AMP.RMCUBE.FITS" ]] && \
           grep -q "auto_nsigma: enabled, multiplier=5" "$rmc_snr_log1"; then
            pass "rmclean_cubes: auto_nsigma mode ran to completion"
        else
            fail "rmclean_cubes: auto_nsigma run did not produce expected output/log (see $rmc_snr_log1)"
        fi
        if cmp -s "${rmc_snr_out1}.RESTORED.AMP.RMCUBE.FITS" \
                  "${rmc_snr_out2}.RESTORED.AMP.RMCUBE.FITS"; then
            pass "rmclean_cubes: auto_nsigma is deterministic (bit-identical across two runs -- no whole-cube noise_seed pre-scan anymore, per-pixel tail sigma has no randomness)"
        else
            fail "rmclean_cubes: auto_nsigma runs were not bit-identical"
        fi
    else
        fail "rmclean_cubes: auto_nsigma run(s) failed (see $rmc_snr_log1 / $rmc_snr_log2)"
    fi

    # --- abs_flux_floor= and auto_nsigma= MAY now be combined (a real
    # design change from the old mutually-exclusive threshold=/
    # threshold_snr= pair -- whichever fires first wins per pixel). ---
    rmc_both_log="$OUT_DIR/rmc_thresh_both.log"
    rmc_both_out="$OUT_DIR/rmc_thresh_both"
    rm -f "${rmc_both_out}".*.RMCUBE.FITS
    if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" phafile="$rmc_pha" \
            maskfile="$rmc_mask" outfile="$rmc_both_out" \
            abs_flux_floor=0.01 auto_nsigma=5.0 niter=200 gain=0.1 \
            > "$rmc_both_log" 2>&1; then
        pass "rmclean_cubes: abs_flux_floor= and auto_nsigma= together now succeeds (combinable, not mutually exclusive)"
    else
        fail "rmclean_cubes: abs_flux_floor=+auto_nsigma= together unexpectedly failed (see $rmc_both_log)"
    fi

    # --- Neither given is also now explicitly VALID (niter-only), with
    # an informational NOTE rather than a refusal. ---
    rmc_neither_log="$OUT_DIR/rmc_thresh_neither.log"
    rmc_neither_out="$OUT_DIR/rmc_thresh_neither"
    rm -f "${rmc_neither_out}".*.RMCUBE.FITS
    if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" phafile="$rmc_pha" \
            maskfile="$rmc_mask" outfile="$rmc_neither_out" \
            niter=200 gain=0.1 > "$rmc_neither_log" 2>&1 && \
       grep -q "NOTE: neither abs_flux_floor= nor auto_nsigma=" "$rmc_neither_log"; then
        pass "rmclean_cubes: omitting both abs_flux_floor= and auto_nsigma= is valid (niter-only), with the expected NOTE"
    else
        fail "rmclean_cubes: niter-only run failed or was missing the expected NOTE (see $rmc_neither_log)"
    fi

    rmc_badunit_log="$OUT_DIR/rmc_thresh_badunit.log"
    if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" phafile="$rmc_pha" \
            maskfile="$rmc_mask" outfile="$OUT_DIR/rmc_thresh_badunit" \
            abs_flux_floor=10Foo > "$rmc_badunit_log" 2>&1; then
        fail "rmclean_cubes: abs_flux_floor=10Foo (unrecognised unit) should have been refused, but succeeded"
    elif grep -q "abs_flux_floor must be a number" "$rmc_badunit_log"; then
        pass "rmclean_cubes: unrecognised abs_flux_floor unit is refused"
    else
        fail "rmclean_cubes: bad-unit rejection had the wrong message (see $rmc_badunit_log)"
    fi
else
    fail "rmclean_cubes: abs_flux_floor-unit tests skipped, section 29's own rmc_lsqref_amp cube is missing"
fi

# ---------------------------------------------------------------------------
# 32. gaussft_mod: FFT-friendly padding (found via the real Jennifer e2e
#     run -- 4501x4501 = 7 x 643, 643 prime, made convolve_cubes/
#     match_cubes' own convolution step ~2x slower than necessary; see
#     next_fast_fft_size's own comment in src/gaussft.f90).
# ---------------------------------------------------------------------------
section "32. gaussft_mod: FFT-friendly padding (next_fast_fft_size)"

gfp_bin="$OUT_DIR/test_gaussft_padding"
gfp_log="$OUT_DIR/test_gaussft_padding.log"
if gfortran -cpp -std=gnu -fallow-argument-mismatch -ffree-line-length-none \
        -O3 -fopenmp -J"$OUT_DIR" \
        src/gaussft.f90 "$TESTS_DIR/test_gaussft_padding.f90" \
        -o "$gfp_bin" -lfftw3 2>"$OUT_DIR/gaussft_padding_build.log"; then
    if "$gfp_bin" > "$gfp_log" 2>&1; then
        while IFS= read -r line; do
            if [[ "$line" == *"[PASS]"* || "$line" == *"[FAIL]"* ]]; then
                echo "  $line"
            fi
        done < "$gfp_log"
        if grep -q "^\[PASS\] test_gaussft_padding" "$gfp_log"; then
            pass "gaussft_mod FFT padding: all checks passed"
        else
            fail "gaussft_mod FFT padding: one or more checks failed (see $gfp_log)"
        fi
    else
        fail "gaussft_mod FFT padding: program exited non-zero (see $gfp_log)"
    fi
else
    fail "gaussft_mod FFT padding: build failed (see $OUT_DIR/gaussft_padding_build.log)"
fi

# ---------------------------------------------------------------------------
# 33. convolve_cubes: io_overlap bit-identical to io_overlap=n (planning-
#     doc ticket, added alongside the FFT padding fix -- same real
#     Jennifer e2e run found convolve's own block write fully serial
#     with compute, ~44s/block dead time on real storage).
# ---------------------------------------------------------------------------
section "33. convolve_cubes: io_overlap bit-identical to io_overlap=n"

if [[ -x bin/convolve_cubes ]]; then
    iob_beamfile="$OUT_DIR/iob_beamlog.txt"
    awk 'BEGIN { for (i=1;i<=200;i++) print i, 10.0, 10.0, 0.0 }' > "$iob_beamfile"

    iob_src="$OUT_DIR/iob_src.Q.FITSCUBE"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$iob_src"

    # strip_fits_ext (src/convolve_cubes.f90) strips ANY trailing
    # extension before appending outsuffix, so .FITSCUBE drops here.
    iob_off_out="${iob_src%.FITSCUBE}_off.CONV.FITS"
    iob_on_out="${iob_src%.FITSCUBE}_on.CONV.FITS"
    rm -f "$iob_off_out" "$iob_on_out"

    iob_off_log="$OUT_DIR/convolve_io_overlap_off.log"
    iob_on_log="$OUT_DIR/convolve_io_overlap_on.log"

    if bin/convolve_cubes infiles="$iob_src" beamfiles="$iob_beamfile" \
            outsuffix="_off.CONV.FITS" target_bmaj=20.0 target_bmin=20.0 \
            target_bpa=0.0 mem_frac_ram=0.1 io_overlap=n \
            > "$iob_off_log" 2>&1 && \
       bin/convolve_cubes infiles="$iob_src" beamfiles="$iob_beamfile" \
            outsuffix="_on.CONV.FITS" target_bmaj=20.0 target_bmin=20.0 \
            target_bpa=0.0 mem_frac_ram=0.1 io_overlap=y \
            > "$iob_on_log" 2>&1; then
        if [[ -s "$iob_off_out" && -s "$iob_on_out" ]] && \
           cmp -s "$iob_off_out" "$iob_on_out"; then
            pass "convolve_cubes: io_overlap=y bit-identical to io_overlap=n"
        else
            fail "convolve_cubes: io_overlap=y output differs from io_overlap=n (see $iob_off_log / $iob_on_log)"
        fi
    else
        fail "convolve_cubes: io_overlap on/off run(s) failed (see $iob_off_log / $iob_on_log)"
    fi
else
    skip "bin/convolve_cubes not built; skipping io_overlap bit-identical test"
fi

# ---------------------------------------------------------------------------
# 34. match_cubes: io_overlap bit-identical to io_overlap=n (stages=
#     convolve -- process_one_file_restricted's own separate port of the
#     same io_overlap mechanism convolve_cubes.f90 uses).
# ---------------------------------------------------------------------------
section "34. match_cubes: io_overlap bit-identical to io_overlap=n (stages=convolve)"

if [[ -x bin/match_cubes ]]; then
    mciob_beamfile="$OUT_DIR/mciob_beamlog.txt"
    awk 'BEGIN { for (i=1;i<=200;i++) print i, 10.0, 10.0, 0.0 }' > "$mciob_beamfile"

    mciob_src="$OUT_DIR/mciob_src.Q.FITSCUBE"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$mciob_src"

    # strip_fits_ext (src/match_cubes.f90) strips ANY trailing extension
    # before appending outsuffix, so .FITSCUBE drops here.
    mciob_off_out="${mciob_src%.FITSCUBE}_off.MATCHED.FITS"
    mciob_on_out="${mciob_src%.FITSCUBE}_on.MATCHED.FITS"
    rm -f "$mciob_off_out" "$mciob_on_out"

    mciob_off_log="$OUT_DIR/match_io_overlap_off.log"
    mciob_on_log="$OUT_DIR/match_io_overlap_on.log"

    if bin/match_cubes stages=convolve infiles="$mciob_src" \
            beamfiles="$mciob_beamfile" outsuffix="_off.MATCHED.FITS" \
            target_bmaj=20.0 target_bmin=20.0 target_bpa=0.0 \
            mem_frac_ram=0.1 io_overlap=n > "$mciob_off_log" 2>&1 && \
       bin/match_cubes stages=convolve infiles="$mciob_src" \
            beamfiles="$mciob_beamfile" outsuffix="_on.MATCHED.FITS" \
            target_bmaj=20.0 target_bmin=20.0 target_bpa=0.0 \
            mem_frac_ram=0.1 io_overlap=y > "$mciob_on_log" 2>&1; then
        if [[ -s "$mciob_off_out" && -s "$mciob_on_out" ]] && \
           cmp -s "$mciob_off_out" "$mciob_on_out"; then
            pass "match_cubes: io_overlap=y bit-identical to io_overlap=n (stages=convolve)"
        else
            fail "match_cubes: io_overlap=y output differs from io_overlap=n (see $mciob_off_log / $mciob_on_log)"
        fi
    else
        fail "match_cubes: io_overlap on/off run(s) failed (see $mciob_off_log / $mciob_on_log)"
    fi
else
    skip "bin/match_cubes not built; skipping io_overlap bit-identical test"
fi

# ---------------------------------------------------------------------------
# 35. reproject_cubes: io_overlap bit-identical to io_overlap=n (planning-
#     doc ticket, added alongside the same fix in convolve_cubes.f90/
#     match_cubes.f90 -- "do not forget reprojection").
# ---------------------------------------------------------------------------
section "35. reproject_cubes: io_overlap bit-identical to io_overlap=n"

if [[ -x bin/reproject_cubes ]]; then
    rcio_ref="$DATA_DIR/TEST.Q.FITSCUBE"
    rcio_src="$OUT_DIR/rcio_src.Q.FITSCUBE"
    cp "$DATA_DIR/TEST_BAND2_MISMATCH.Q.FITSCUBE" "$rcio_src"

    # strip_fits_ext (src/reproject_cubes.f90) strips ANY trailing
    # extension before appending "_REPROJ.FITS", so .FITSCUBE drops here.
    rcio_src_reproj="${rcio_src%.FITSCUBE}_REPROJ.FITS"
    rcio_off_out="${rcio_src}_off.FITS"
    rcio_on_out="${rcio_src}_on.FITS"
    rm -f "$rcio_src_reproj" "$rcio_off_out" "$rcio_on_out"

    rcio_off_log="$OUT_DIR/reproject_io_overlap_off.log"
    rcio_on_log="$OUT_DIR/reproject_io_overlap_on.log"

    if bin/reproject_cubes mode=reference reffile="$rcio_ref" \
            infiles="$rcio_src" mem_frac_ram=0.1 io_overlap=n \
            > "$rcio_off_log" 2>&1; then
        mv "$rcio_src_reproj" "$rcio_off_out"
    fi
    if bin/reproject_cubes mode=reference reffile="$rcio_ref" \
            infiles="$rcio_src" mem_frac_ram=0.1 io_overlap=y \
            > "$rcio_on_log" 2>&1; then
        mv "$rcio_src_reproj" "$rcio_on_out"
    fi

    if [[ -s "$rcio_off_out" && -s "$rcio_on_out" ]] && \
       cmp -s "$rcio_off_out" "$rcio_on_out"; then
        pass "reproject_cubes: io_overlap=y bit-identical to io_overlap=n"
    else
        fail "reproject_cubes: io_overlap=y output differs from (or is missing vs) io_overlap=n (see $rcio_off_log / $rcio_on_log)"
    fi
else
    skip "bin/reproject_cubes not built; skipping io_overlap bit-identical test"
fi

# ---------------------------------------------------------------------------
# 36. match_cubes: io_overlap bit-identical to io_overlap=n (stages=
#     reproject -- process_one_file_general's own separate port).
# ---------------------------------------------------------------------------
section "36. match_cubes: io_overlap bit-identical to io_overlap=n (stages=reproject)"

if [[ -x bin/match_cubes ]]; then
    mcrio_ref="$DATA_DIR/TEST.Q.FITSCUBE"
    mcrio_src="$OUT_DIR/mcrio_src.Q.FITSCUBE"
    cp "$DATA_DIR/TEST_BAND2_MISMATCH.Q.FITSCUBE" "$mcrio_src"

    # strip_fits_ext (src/match_cubes.f90) strips ANY trailing extension
    # before appending outsuffix, so .FITSCUBE drops here.
    mcrio_src_reproj="${mcrio_src%.FITSCUBE}_REPROJ.FITS"
    mcrio_off_out="${mcrio_src}_off.FITS"
    mcrio_on_out="${mcrio_src}_on.FITS"
    rm -f "$mcrio_src_reproj" "$mcrio_off_out" "$mcrio_on_out"

    mcrio_off_log="$OUT_DIR/match_reproject_io_overlap_off.log"
    mcrio_on_log="$OUT_DIR/match_reproject_io_overlap_on.log"

    if bin/match_cubes stages=reproject footprint_mode=reference \
            reffile="$mcrio_ref" infiles="$mcrio_src" outsuffix="_REPROJ.FITS" \
            mem_frac_ram=0.1 io_overlap=n > "$mcrio_off_log" 2>&1; then
        mv "$mcrio_src_reproj" "$mcrio_off_out"
    fi
    if bin/match_cubes stages=reproject footprint_mode=reference \
            reffile="$mcrio_ref" infiles="$mcrio_src" outsuffix="_REPROJ.FITS" \
            mem_frac_ram=0.1 io_overlap=y > "$mcrio_on_log" 2>&1; then
        mv "$mcrio_src_reproj" "$mcrio_on_out"
    fi

    if [[ -s "$mcrio_off_out" && -s "$mcrio_on_out" ]] && \
       cmp -s "$mcrio_off_out" "$mcrio_on_out"; then
        pass "match_cubes: io_overlap=y bit-identical to io_overlap=n (stages=reproject)"
    else
        fail "match_cubes: io_overlap=y output differs from (or is missing vs) io_overlap=n (see $mcrio_off_log / $mcrio_on_log)"
    fi
else
    skip "bin/match_cubes not built; skipping io_overlap bit-identical test"
fi

# ---------------------------------------------------------------------------
# 37. RM-CLEAN clean_complex stop-reason (stopped_by_threshold) +
#     per-iteration trace outputs (docs/dev/RMCLEAN_INTEGRATION_PLAN.md,
#     T6-adjacent convergence diagnostic)
# ---------------------------------------------------------------------------
section "37. RM-CLEAN clean_complex stop-reason + iteration trace"

stopreason_bin="$RMCLEAN_BUILD_DIR/test_rmclean_stop_reason"
stopreason_log="$OUT_DIR/test_rmclean_stop_reason.log"

if [[ -f "$rmclean_o" ]]; then
    if gfortran -cpp -std=gnu -fallow-argument-mismatch -ffree-line-length-none \
            -O3 -fopenmp -I"$rmclean_mod_dir" -J"$rmclean_mod_dir" \
            "$TESTS_DIR/test_rmclean_stop_reason.f90" "$rmclean_o" \
            -o "$stopreason_bin" -lfftw3 2>"$OUT_DIR/rmclean_stopreason_build.log"; then
        if "$stopreason_bin" > "$stopreason_log" 2>&1; then
            while IFS= read -r line; do
                case "$line" in
                    *"[PASS]"*) pass "${line#*\[PASS\] }" ;;
                    *"[FAIL]"*) fail "${line#*\[FAIL\] }" ;;
                esac
            done < "$stopreason_log"
        else
            fail "RM-CLEAN stop-reason: program exited non-zero (see $stopreason_log)"
        fi
    else
        fail "RM-CLEAN stop-reason: build failed (see $OUT_DIR/rmclean_stopreason_build.log)"
    fi
else
    skip "rmclean_mod.o not built (section 23 skipped); skipping RM-CLEAN stop-reason test"
fi

# ---------------------------------------------------------------------------
# 38-40. Multi-band preprocessing END-TO-END (docs/dev/
#     MULTI_BAND_TOMOGRAPHY_PLAN.md T15): sections 15/30/33-36 test the
#     preprocessing tools' own internal correctness (geometry validation,
#     skip-if-already-matched, io_overlap consistency) in isolation, but
#     none of them prove the actual point of the toolchain -- take
#     genuinely mismatched multi-band data, fix it, and get the right
#     answer out of rm_synthesis. These three sections do exactly that,
#     each isolating one failure mode: grid-only (reproject_cubes alone),
#     resolution-only (convolve_cubes alone), and both together
#     (match_cubes stages=both). TEST_BAND2_UNMATCHED/band1_beamlog.txt/
#     band2_beamlog.txt (tests/make_test_cubes.py) are dedicated fixtures
#     for this -- band2's own beam (20") is genuinely different from
#     band1's (10"), so convolve_cubes has real smoothing work to do, not
#     a no-op copy.
# ---------------------------------------------------------------------------
section "38. Multi-band preprocessing: reproject_cubes fixes a grid-only mismatch, then rm_synthesis recovers the known sources (T15)"

if [[ -x bin/reproject_cubes && -x "$BIN_SERIAL" ]]; then
    mbpp_grid_band1q="$OUT_DIR/mbpp_grid.TEST.Q.FITSCUBE"
    mbpp_grid_band1u="$OUT_DIR/mbpp_grid.TEST.U.FITSCUBE"
    mbpp_grid_band2q="$OUT_DIR/mbpp_grid.TEST_BAND2_MISMATCH.Q.FITSCUBE"
    mbpp_grid_band2u="$OUT_DIR/mbpp_grid.TEST_BAND2.U.FITSCUBE"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$mbpp_grid_band1q"
    cp "$DATA_DIR/TEST.U.FITSCUBE" "$mbpp_grid_band1u"
    cp "$DATA_DIR/TEST_BAND2_MISMATCH.Q.FITSCUBE" "$mbpp_grid_band2q"
    cp "$DATA_DIR/TEST_BAND2.U.FITSCUBE" "$mbpp_grid_band2u"
    # strip_fits_ext (src/reproject_cubes.f90) strips ANY trailing
    # extension, so the tool's own output name drops .FITSCUBE here --
    # mirror that with bash's own extension strip rather than
    # concatenating onto the untouched input name.
    mbpp_grid_band2q_reproj="${mbpp_grid_band2q%.FITSCUBE}_REPROJ.FITS"
    rm -f "$mbpp_grid_band2q_reproj"

    mbpp_grid_reproj_log="$OUT_DIR/mbpp_grid_reproject.log"
    if bin/reproject_cubes mode=reference reffile="$mbpp_grid_band1q" \
            infiles="$mbpp_grid_band2q" > "$mbpp_grid_reproj_log" 2>&1 && \
       [[ -s "$mbpp_grid_band2q_reproj" ]]; then
        pass "reproject_cubes: fixed a genuine grid-only mismatch, wrote output"
    else
        fail "reproject_cubes: grid-only mismatch fix failed or wrote no output (see $mbpp_grid_reproj_log)"
    fi

    if [[ -s "$mbpp_grid_band2q_reproj" ]]; then
        mbpp_grid_cfg="$OUT_DIR/mbpp_grid_rmsynth.cfg"
        cat > "$mbpp_grid_cfg" <<CFGEOF
path                = ${OUT_DIR}/
infileQ             = $(basename "$mbpp_grid_band1q"),$(basename "$mbpp_grid_band2q_reproj")
infileU             = $(basename "$mbpp_grid_band1u"),$(basename "$mbpp_grid_band2u")
outfile             = ${OUT_DIR}/mbpp_grid_out
remove_badchan      = n
global_badchan_file = /dev/null,/dev/null
subim               = n
rem_mean            = 0
remove_qu_bias      = n
resiQ               = 0.0,0.0
slopeQ              = 0.0,0.0
resiU               = 0.0,0.0
slopeU              = 0.0,0.0
ofac                = 1
fac                 = 3.14159265358979
use_auto_rm_range   = 0
beg_rm              = -50.0
end_rm              = 50.0
nrm                 = 201
output_mode         = ap
ap_angle_mode       = phase
write_mask_output   = y
write_nvalid_output = y
use_gpu             = n
CFGEOF
        mbpp_grid_rmsynth_log="$OUT_DIR/mbpp_grid_rmsynth.log"
        "$BIN_SERIAL" "$mbpp_grid_cfg" > "$mbpp_grid_rmsynth_log" 2>&1
        if [[ -f "${OUT_DIR}/mbpp_grid_out.AMP.RMCUBE.FITS" ]] && \
           python3 "$TESTS_DIR/check_rm_peak.py" "${OUT_DIR}/mbpp_grid_out.AMP.RMCUBE.FITS" "$TRUTH" > /dev/null 2>&1; then
            pass "reproject_cubes fix + rm_synthesis: src_A/src_B recovered after fixing a real grid-only mismatch"
        else
            fail "reproject_cubes fix + rm_synthesis: RM peak(s) not recovered (see $mbpp_grid_rmsynth_log)"
        fi
    else
        fail "reproject_cubes: grid-only mismatch fix did not produce a usable output; skipping downstream rm_synthesis check"
    fi
else
    skip "reproject_cubes or serial rm_synthesis binary not available; skipping grid-only preprocessing test"
fi

# ---------------------------------------------------------------------------
section "39. Multi-band preprocessing: convolve_cubes fixes a resolution-only mismatch, then rm_synthesis recovers the known sources (T15)"

if [[ -x bin/convolve_cubes && -x "$BIN_SERIAL" ]]; then
    mbpp_res_band1q="$OUT_DIR/mbpp_res.TEST.Q.FITSCUBE"
    mbpp_res_band1u="$OUT_DIR/mbpp_res.TEST.U.FITSCUBE"
    mbpp_res_band2q="$OUT_DIR/mbpp_res.TEST_BAND2.Q.FITSCUBE"
    mbpp_res_band2u="$OUT_DIR/mbpp_res.TEST_BAND2.U.FITSCUBE"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$mbpp_res_band1q"
    cp "$DATA_DIR/TEST.U.FITSCUBE" "$mbpp_res_band1u"
    cp "$DATA_DIR/TEST_BAND2.Q.FITSCUBE" "$mbpp_res_band2q"
    cp "$DATA_DIR/TEST_BAND2.U.FITSCUBE" "$mbpp_res_band2u"
    # strip_fits_ext (src/convolve_cubes.f90) strips ANY trailing
    # extension, so the tool's own output name drops .FITSCUBE here.
    mbpp_res_band1q_conv="${mbpp_res_band1q%.FITSCUBE}_CONV.FITS"
    mbpp_res_band1u_conv="${mbpp_res_band1u%.FITSCUBE}_CONV.FITS"
    mbpp_res_band2q_conv="${mbpp_res_band2q%.FITSCUBE}_CONV.FITS"
    mbpp_res_band2u_conv="${mbpp_res_band2u%.FITSCUBE}_CONV.FITS"
    rm -f "$mbpp_res_band1q_conv" "$mbpp_res_band1u_conv" \
          "$mbpp_res_band2q_conv" "$mbpp_res_band2u_conv"

    mbpp_res_convolve_log="$OUT_DIR/mbpp_res_convolve.log"
    if bin/convolve_cubes \
            infiles="$mbpp_res_band1q,$mbpp_res_band1u,$mbpp_res_band2q,$mbpp_res_band2u" \
            beamfiles="$DATA_DIR/band1_beamlog.txt,$DATA_DIR/band1_beamlog.txt,$DATA_DIR/band2_beamlog.txt,$DATA_DIR/band2_beamlog.txt" \
            mem_frac_ram=0.25 > "$mbpp_res_convolve_log" 2>&1 && \
       [[ -s "$mbpp_res_band1q_conv" && -s "$mbpp_res_band2q_conv" ]]; then
        pass "convolve_cubes: fixed a genuine resolution-only mismatch (band 1 smoothed 10\"->~20\"), wrote output"
    else
        fail "convolve_cubes: resolution-only mismatch fix failed or wrote no output (see $mbpp_res_convolve_log)"
    fi

    if [[ -s "$mbpp_res_band1q_conv" && -s "$mbpp_res_band2q_conv" ]]; then
        mbpp_res_cfg="$OUT_DIR/mbpp_res_rmsynth.cfg"
        cat > "$mbpp_res_cfg" <<CFGEOF
path                = ${OUT_DIR}/
infileQ             = $(basename "$mbpp_res_band1q_conv"),$(basename "$mbpp_res_band2q_conv")
infileU             = $(basename "$mbpp_res_band1u_conv"),$(basename "$mbpp_res_band2u_conv")
outfile             = ${OUT_DIR}/mbpp_res_out
remove_badchan      = n
global_badchan_file = /dev/null,/dev/null
subim               = n
rem_mean            = 0
remove_qu_bias      = n
resiQ               = 0.0,0.0
slopeQ              = 0.0,0.0
resiU               = 0.0,0.0
slopeU              = 0.0,0.0
ofac                = 1
fac                 = 3.14159265358979
use_auto_rm_range   = 0
beg_rm              = -50.0
end_rm              = 50.0
nrm                 = 201
output_mode         = ap
ap_angle_mode       = phase
write_mask_output   = y
write_nvalid_output = y
use_gpu             = n
CFGEOF
        mbpp_res_rmsynth_log="$OUT_DIR/mbpp_res_rmsynth.log"
        "$BIN_SERIAL" "$mbpp_res_cfg" > "$mbpp_res_rmsynth_log" 2>&1
        if [[ -f "${OUT_DIR}/mbpp_res_out.AMP.RMCUBE.FITS" ]] && \
           python3 "$TESTS_DIR/check_rm_peak.py" "${OUT_DIR}/mbpp_res_out.AMP.RMCUBE.FITS" "$TRUTH" > /dev/null 2>&1; then
            pass "convolve_cubes fix + rm_synthesis: src_A/src_B recovered after fixing a real resolution-only mismatch"
        else
            fail "convolve_cubes fix + rm_synthesis: RM peak(s) not recovered (see $mbpp_res_rmsynth_log)"
        fi
    else
        fail "convolve_cubes: resolution-only mismatch fix did not produce usable output; skipping downstream rm_synthesis check"
    fi
else
    skip "convolve_cubes or serial rm_synthesis binary not available; skipping resolution-only preprocessing test"
fi

# ---------------------------------------------------------------------------
section "40. Multi-band preprocessing: match_cubes stages=both fixes grid+resolution mismatch together, then rm_synthesis recovers the known sources (T15)"

if [[ -x bin/match_cubes && -x "$BIN_SERIAL" ]]; then
    mbpp_both_band1q="$OUT_DIR/mbpp_both.TEST.Q.FITSCUBE"
    mbpp_both_band1u="$OUT_DIR/mbpp_both.TEST.U.FITSCUBE"
    mbpp_both_band2q="$OUT_DIR/mbpp_both.TEST_BAND2_UNMATCHED.Q.FITSCUBE"
    mbpp_both_band2u="$OUT_DIR/mbpp_both.TEST_BAND2_UNMATCHED.U.FITSCUBE"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$mbpp_both_band1q"
    cp "$DATA_DIR/TEST.U.FITSCUBE" "$mbpp_both_band1u"
    cp "$DATA_DIR/TEST_BAND2_UNMATCHED.Q.FITSCUBE" "$mbpp_both_band2q"
    cp "$DATA_DIR/TEST_BAND2_UNMATCHED.U.FITSCUBE" "$mbpp_both_band2u"
    # strip_fits_ext (src/match_cubes.f90) strips ANY trailing
    # extension, so the tool's own output name drops .FITSCUBE here.
    mbpp_both_band1q_matched="${mbpp_both_band1q%.FITSCUBE}_MATCHED.FITS"
    mbpp_both_band1u_matched="${mbpp_both_band1u%.FITSCUBE}_MATCHED.FITS"
    mbpp_both_band2q_matched="${mbpp_both_band2q%.FITSCUBE}_MATCHED.FITS"
    mbpp_both_band2u_matched="${mbpp_both_band2u%.FITSCUBE}_MATCHED.FITS"
    rm -f "$mbpp_both_band1q_matched" "$mbpp_both_band1u_matched" \
          "$mbpp_both_band2q_matched" "$mbpp_both_band2u_matched"

    mbpp_both_match_log="$OUT_DIR/mbpp_both_match.log"
    if bin/match_cubes stages=both order=convolve_reproject \
            footprint_mode=reference reffile="$mbpp_both_band1q" \
            infiles="$mbpp_both_band1q,$mbpp_both_band1u,$mbpp_both_band2q,$mbpp_both_band2u" \
            beamfiles="$DATA_DIR/band1_beamlog.txt,$DATA_DIR/band1_beamlog.txt,$DATA_DIR/band2_beamlog.txt,$DATA_DIR/band2_beamlog.txt" \
            mem_frac_ram=0.25 > "$mbpp_both_match_log" 2>&1 && \
       [[ -s "$mbpp_both_band1q_matched" && -s "$mbpp_both_band2q_matched" ]]; then
        pass "match_cubes stages=both: fixed a genuine grid+resolution mismatch together, wrote output"
    else
        fail "match_cubes stages=both: grid+resolution mismatch fix failed or wrote no output (see $mbpp_both_match_log)"
    fi

    if [[ -s "$mbpp_both_band1q_matched" && -s "$mbpp_both_band2q_matched" ]]; then
        mbpp_both_cfg="$OUT_DIR/mbpp_both_rmsynth.cfg"
        cat > "$mbpp_both_cfg" <<CFGEOF
path                = ${OUT_DIR}/
infileQ             = $(basename "$mbpp_both_band1q_matched"),$(basename "$mbpp_both_band2q_matched")
infileU             = $(basename "$mbpp_both_band1u_matched"),$(basename "$mbpp_both_band2u_matched")
outfile             = ${OUT_DIR}/mbpp_both_out
remove_badchan      = n
global_badchan_file = /dev/null,/dev/null
subim               = n
rem_mean            = 0
remove_qu_bias      = n
resiQ               = 0.0,0.0
slopeQ              = 0.0,0.0
resiU               = 0.0,0.0
slopeU              = 0.0,0.0
ofac                = 1
fac                 = 3.14159265358979
use_auto_rm_range   = 0
beg_rm              = -50.0
end_rm              = 50.0
nrm                 = 201
output_mode         = ap
ap_angle_mode       = phase
write_mask_output   = y
write_nvalid_output = y
use_gpu             = n
CFGEOF
        mbpp_both_rmsynth_log="$OUT_DIR/mbpp_both_rmsynth.log"
        "$BIN_SERIAL" "$mbpp_both_cfg" > "$mbpp_both_rmsynth_log" 2>&1
        if [[ -f "${OUT_DIR}/mbpp_both_out.AMP.RMCUBE.FITS" ]] && \
           python3 "$TESTS_DIR/check_rm_peak.py" "${OUT_DIR}/mbpp_both_out.AMP.RMCUBE.FITS" "$TRUTH" > /dev/null 2>&1; then
            pass "match_cubes stages=both fix + rm_synthesis: src_A/src_B recovered after fixing a real grid+resolution mismatch"
        else
            fail "match_cubes stages=both fix + rm_synthesis: RM peak(s) not recovered (see $mbpp_both_rmsynth_log)"
        fi
    else
        fail "match_cubes stages=both: grid+resolution mismatch fix did not produce usable output; skipping downstream rm_synthesis check"
    fi
else
    skip "match_cubes or serial rm_synthesis binary not available; skipping grid+resolution preprocessing test"
fi

# ---------------------------------------------------------------------------
# 41. convolve_cubes/match_cubes: badchan_file is genuinely per-band (T16,
#     docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md) -- a real bug found while
#     preparing a real multi-band run: badchan_file used to be a single
#     shared list applied identically to every infile regardless of that
#     infile's own channel numbering, so band A's own bad channel would
#     also get (incorrectly) excluded from band B, and vice versa. Fixed
#     to be a comma list, one entry per infile, exactly like beamfiles.
#     Verifies each band's own designated bad channel is excluded from
#     THAT band's own output, and NOT from the other band's.
# ---------------------------------------------------------------------------
section "41. convolve_cubes/match_cubes: badchan_file is genuinely per-band (T16)"

if [[ -x bin/convolve_cubes ]]; then
    pbbc_band1="$OUT_DIR/pbbc_band1.Q.FITSCUBE"
    pbbc_band2="$OUT_DIR/pbbc_band2.Q.FITSCUBE"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$pbbc_band1"
    cp "$DATA_DIR/TEST_BAND2.Q.FITSCUBE" "$pbbc_band2"
    pbbc_band1_badchan="$OUT_DIR/pbbc_band1_badchan.txt"
    pbbc_band2_badchan="$OUT_DIR/pbbc_band2_badchan.txt"
    printf '13\n' > "$pbbc_band1_badchan"
    printf '77\n' > "$pbbc_band2_badchan"
    # strip_fits_ext strips ANY trailing extension, so .FITSCUBE drops here.
    # target beam (30") deliberately exceeds BOTH band1_beamlog.txt (10")
    # and band2_beamlog.txt (20") native beams, so neither is skipped as
    # already-matching -- both must genuinely convolve and write output.
    pbbc_band1_conv="${pbbc_band1%.FITSCUBE}_CONV.FITS"
    pbbc_band2_conv="${pbbc_band2%.FITSCUBE}_CONV.FITS"
    rm -f "$pbbc_band1_conv" "$pbbc_band2_conv"

    pbbc_convolve_log="$OUT_DIR/pbbc_convolve.log"
    if bin/convolve_cubes infiles="$pbbc_band1,$pbbc_band2" \
            beamfiles="$DATA_DIR/band1_beamlog.txt,$DATA_DIR/band2_beamlog.txt" \
            badchan_file="$pbbc_band1_badchan,$pbbc_band2_badchan" \
            target_bmaj=30.0 target_bmin=30.0 target_bpa=0.0 \
            > "$pbbc_convolve_log" 2>&1; then
        pbbc_band1_bad_planes=$(python3 -c "
from astropy.io import fits
import numpy as np
data = fits.getdata('$pbbc_band1_conv')[0]
print(','.join(str(i+1) for i in range(data.shape[0]) if np.isnan(data[i]).all()))
")
        pbbc_band2_bad_planes=$(python3 -c "
from astropy.io import fits
import numpy as np
data = fits.getdata('$pbbc_band2_conv')[0]
print(','.join(str(i+1) for i in range(data.shape[0]) if np.isnan(data[i]).all()))
")
        if [[ "$pbbc_band1_bad_planes" == "13" && "$pbbc_band2_bad_planes" == "77" ]]; then
            pass "convolve_cubes: badchan_file is genuinely per-band (band1 only chan 13 bad, band2 only chan 77 bad, no cross-contamination)"
        else
            fail "convolve_cubes: badchan_file per-band check failed (band1 bad=[$pbbc_band1_bad_planes], expected 13; band2 bad=[$pbbc_band2_bad_planes], expected 77) (see $pbbc_convolve_log)"
        fi
    else
        fail "convolve_cubes: per-band badchan_file run failed (see $pbbc_convolve_log)"
    fi
else
    skip "bin/convolve_cubes not built; skipping per-band badchan_file test"
fi

if [[ -x bin/match_cubes ]]; then
    pbbcm_band1="$OUT_DIR/pbbcm_band1.Q.FITSCUBE"
    pbbcm_band2="$OUT_DIR/pbbcm_band2.Q.FITSCUBE"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$pbbcm_band1"
    cp "$DATA_DIR/TEST_BAND2.Q.FITSCUBE" "$pbbcm_band2"
    pbbcm_band1_badchan="$OUT_DIR/pbbcm_band1_badchan.txt"
    pbbcm_band2_badchan="$OUT_DIR/pbbcm_band2_badchan.txt"
    printf '13\n' > "$pbbcm_band1_badchan"
    printf '77\n' > "$pbbcm_band2_badchan"
    # strip_fits_ext strips ANY trailing extension, so .FITSCUBE drops here.
    # target beam (30") deliberately exceeds both native beams so neither
    # band is skipped as already-matching.
    pbbcm_band1_conv="${pbbcm_band1%.FITSCUBE}_CONV.FITS"
    pbbcm_band2_conv="${pbbcm_band2%.FITSCUBE}_CONV.FITS"
    rm -f "$pbbcm_band1_conv" "$pbbcm_band2_conv"

    pbbcm_convolve_log="$OUT_DIR/pbbcm_convolve.log"
    if bin/match_cubes stages=convolve infiles="$pbbcm_band1,$pbbcm_band2" \
            beamfiles="$DATA_DIR/band1_beamlog.txt,$DATA_DIR/band2_beamlog.txt" \
            badchan_file="$pbbcm_band1_badchan,$pbbcm_band2_badchan" \
            target_bmaj=30.0 target_bmin=30.0 target_bpa=0.0 \
            > "$pbbcm_convolve_log" 2>&1; then
        pbbcm_band1_bad_planes=$(python3 -c "
from astropy.io import fits
import numpy as np
data = fits.getdata('$pbbcm_band1_conv')[0]
print(','.join(str(i+1) for i in range(data.shape[0]) if np.isnan(data[i]).all()))
")
        pbbcm_band2_bad_planes=$(python3 -c "
from astropy.io import fits
import numpy as np
data = fits.getdata('$pbbcm_band2_conv')[0]
print(','.join(str(i+1) for i in range(data.shape[0]) if np.isnan(data[i]).all()))
")
        if [[ "$pbbcm_band1_bad_planes" == "13" && "$pbbcm_band2_bad_planes" == "77" ]]; then
            pass "match_cubes: badchan_file is genuinely per-band (band1 only chan 13 bad, band2 only chan 77 bad, no cross-contamination)"
        else
            fail "match_cubes: badchan_file per-band check failed (band1 bad=[$pbbcm_band1_bad_planes], expected 13; band2 bad=[$pbbcm_band2_bad_planes], expected 77) (see $pbbcm_convolve_log)"
        fi
    else
        fail "match_cubes: per-band badchan_file run failed (see $pbbcm_convolve_log)"
    fi
else
    skip "bin/match_cubes not built; skipping per-band badchan_file test"
fi

# ---------------------------------------------------------------------------
# 42. convolve_cubes/match_cubes/reproject_cubes: a long infiles= CLI
#     argument (>512 chars combined) is not silently truncated (T17,
#     docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md). Real bug found while
#     running scripts/run_pipeline.sh on real WALLABY+EMU data: its own
#     symlink-redirection scheme (needed to land match_cubes' output on
#     NVMe rather than match_input_path's own disk) lengthens every
#     path enough that a 4-file infiles= argument crossed 512 chars --
#     silently truncated by the old `character(len=512) :: this_arg`
#     CLI-token buffer (Fortran does not error on this, it just cuts
#     the string), producing a nonsense path and a confusing "failed to
#     open" error far from the real cause. Fixed by widening this_arg/
#     cli_val and the raw_infiles/raw_beamfiles/raw_badchan_file CSV
#     staging buffers to 16384 chars in all three tools.
# ---------------------------------------------------------------------------
section "42. convolve_cubes/match_cubes/reproject_cubes: long infiles= argument not truncated (T17)"

lia_dir="$OUT_DIR/long_infiles_arg/a/really/quite/deeply/nested/directory/structure/designed/specifically/to/push/the/combined/infiles/argument/comfortably/past/the/old/five_hundred_twelve/character/limit/for/a/real/regression/test"
mkdir -p "$lia_dir"
cp "$DATA_DIR/TEST.Q.FITSCUBE" "$lia_dir/bandA.fits"
cp "$DATA_DIR/TEST_BAND2.Q.FITSCUBE" "$lia_dir/bandB.fits"
awk 'BEGIN{for(i=1;i<=200;i++) print i, 10.0, 10.0, 0.0}' > "$lia_dir/bandA_beamlog.txt"
awk 'BEGIN{for(i=1;i<=150;i++) print i, 10.0, 10.0, 0.0}' > "$lia_dir/bandB_beamlog.txt"
lia_infiles="$lia_dir/bandA.fits,$lia_dir/bandB.fits"

if [[ ${#lia_infiles} -le 512 ]]; then
    skip "long infiles= test fixture path is not actually long enough (${#lia_infiles} chars) to exercise the old 512-char limit -- OUT_DIR is unusually short on this machine"
else
    if [[ -x bin/convolve_cubes ]]; then
        rm -f "$lia_dir/bandA_CONV.FITS" "$lia_dir/bandB_CONV.FITS"
        lia_convolve_log="$OUT_DIR/lia_convolve.log"
        if bin/convolve_cubes infiles="$lia_infiles" \
                beamfiles="$lia_dir/bandA_beamlog.txt,$lia_dir/bandB_beamlog.txt" \
                target_bmaj=20.0 target_bmin=20.0 target_bpa=0.0 \
                > "$lia_convolve_log" 2>&1 && \
           [[ -s "$lia_dir/bandA_CONV.FITS" && -s "$lia_dir/bandB_CONV.FITS" ]]; then
            pass "convolve_cubes: ${#lia_infiles}-char infiles= argument (>512) not truncated, both bands processed"
        else
            fail "convolve_cubes: long infiles= argument (${#lia_infiles} chars) failed (see $lia_convolve_log)"
        fi
    else
        skip "bin/convolve_cubes not built; skipping long infiles= test"
    fi

    if [[ -x bin/match_cubes ]]; then
        rm -f "$lia_dir/bandA_CONV.FITS" "$lia_dir/bandB_CONV.FITS"
        lia_match_log="$OUT_DIR/lia_match.log"
        if bin/match_cubes stages=convolve infiles="$lia_infiles" \
                beamfiles="$lia_dir/bandA_beamlog.txt,$lia_dir/bandB_beamlog.txt" \
                target_bmaj=20.0 target_bmin=20.0 target_bpa=0.0 \
                > "$lia_match_log" 2>&1 && \
           [[ -s "$lia_dir/bandA_CONV.FITS" && -s "$lia_dir/bandB_CONV.FITS" ]]; then
            pass "match_cubes: ${#lia_infiles}-char infiles= argument (>512) not truncated, both bands processed"
        else
            fail "match_cubes: long infiles= argument (${#lia_infiles} chars) failed (see $lia_match_log)"
        fi
    else
        skip "bin/match_cubes not built; skipping long infiles= test"
    fi

    if [[ -x bin/reproject_cubes ]]; then
        cp "$DATA_DIR/TEST_BAND2_MISMATCH.Q.FITSCUBE" "$lia_dir/bandC.fits"
        lia_infiles_reproj="$lia_dir/bandA.fits,$lia_dir/bandC.fits"
        rm -f "$lia_dir/bandC_REPROJ.FITS"
        lia_reproject_log="$OUT_DIR/lia_reproject.log"
        if bin/reproject_cubes mode=reference reffile="$lia_dir/bandA.fits" \
                infiles="$lia_infiles_reproj" \
                > "$lia_reproject_log" 2>&1 && \
           [[ -s "$lia_dir/bandC_REPROJ.FITS" ]]; then
            pass "reproject_cubes: ${#lia_infiles_reproj}-char infiles= argument (>512) not truncated, genuinely-mismatched band reprojected"
        else
            fail "reproject_cubes: long infiles= argument (${#lia_infiles_reproj} chars) failed (see $lia_reproject_log)"
        fi
    else
        skip "bin/reproject_cubes not built; skipping long infiles= test"
    fi
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
section "Test Summary"
TOTAL=$((PASS + FAIL + SKIP))
echo "Total : $TOTAL"
echo "Pass  : $PASS"
echo "Fail  : $FAIL"
echo "Skip  : $SKIP"

if [[ "$FAIL" -gt 0 ]]; then
    echo ""
    echo "RESULT: FAILED ($FAIL test(s) failed)"
    exit 1
else
    echo ""
    echo "RESULT: ALL PASSED"
    exit 0
fi
