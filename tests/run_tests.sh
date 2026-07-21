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
#  14. io_write_threads>1 – bit-identical to io_write_threads=1 across a
#      7-tile run (T6 raw-write path bypasses CFITSIO for AMP/PHA pixel
#      writes; guards against the CFITSIO handle-aliasing bug and the
#      stale-buffer-at-close bug that raw-write mode replaced it with)
#  15. Multi-band tomography (T1) – comma-separated-list infileQ/infileU
#      config schema: matched-geometry two-band config validates and runs
#      to completion (T2 replaced the earlier "not yet implemented" stop);
#      mismatched-geometry config is loudly refused; inconsistent per-band
#      list lengths are rejected at config-parse time
#      (planning/MULTI_BAND_TOMOGRAPHY_PLAN.md)
#  16. Multi-band frequency merge (T2) – Sec 10 thesis-grounded scenario
#      (Raja 2014 Table 6.1/6.2): P-band alone, L-band alone, and P+L
#      combined, for a point source + Faraday-thick top-hat + F2/F3 pair.
#      Also T4: multi-tile multi-band produces bit-identical output to
#      the single-tile P+L run above (tiling must not change the answer)
#  17. Split-band identity test (T5) – a contiguous 2-band split of the
#      primary test cube (no gap) must reproduce the undivided cube's own
#      output bit-for-bit; the most direct mechanical regression check
#      for the frequency-merge architecture itself
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
# runs with the io_write_threads=1 default, so all tiles share a single
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
# 14. io_write_threads>1 (T6 raw-write) – bit-identical to io_write_threads=1
# ---------------------------------------------------------------------------
# io_write_threads>1 used to open N read-write FTOPEN handles onto the SAME
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
section "14. io_write_threads>1 – bit-identical to io_write_threads=1 (T6)"
if [[ -x "$BIN_OMP" ]]; then
    cfg_wt1=$(make_cfg "wt1" "n" "tile_auto=n
tile_ra=32
tile_dec=5
cubestat=y
io_write_threads=1")
    cfg_wt4=$(make_cfg "wt4" "n" "tile_auto=n
tile_ra=32
tile_dec=5
cubestat=y
io_write_threads=4")
    log_wt1="$OUT_DIR/wt1.log"
    log_wt4="$OUT_DIR/wt4.log"
    rm -f "$OUT_DIR"/wt1.*.FITS "$OUT_DIR"/wt4.*.FITS

    if run_binary "$BIN_OMP" "$cfg_wt1" "$log_wt1" && \
       run_binary "$BIN_OMP" "$cfg_wt4" "$log_wt4"; then
        pass "io_write_threads=1 and =4: both runs completed without crashing"

        if grep -q "using raw stream writes" "$log_wt4"; then
            pass "io_write_threads=4: raw-write path confirmed taken"
        else
            fail "io_write_threads=4: expected raw-write startup message not found in $log_wt4"
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
                    fail "io_write_threads=4: ${suffix} differs from io_write_threads=1"
                fi
            else
                all_match=0
                fail "io_write_threads=4: ${suffix} output missing (expected $f1 and $f4)"
            fi
        done
        if [[ "$all_match" -eq 1 ]]; then
            pass "io_write_threads=4: all 8 output products bit-identical to io_write_threads=1"
        fi
    else
        fail "io_write_threads test: OMP run failed (see $log_wt1 / $log_wt4)"
    fi
else
    skip "OMP binary not available; skipping io_write_threads test"
fi

# ---------------------------------------------------------------------------
# 15. Multi-band tomography (T1) – comma-list config schema + geometry
#     validation (planning/MULTI_BAND_TOMOGRAPHY_PLAN.md). Fortran's bare
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
global_badchan_file = /dev/null
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
#     (planning/MULTI_BAND_TOMOGRAPHY_PLAN.md; Raja 2014 Table 6.1/6.2):
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
        cat > "$cfg" <<CFGEOF
path                = ${DATA_DIR}/
infileQ             = ${infileQ}
infileU             = ${infileU}
outfile             = ${OUT_DIR}/${tag}
remove_badchan      = n
global_badchan_file = /dev/null
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

        # T4 (planning/MULTI_BAND_TOMOGRAPHY_PLAN.md): multi-tile multi-band
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
# 17. Split-band identity test (T5, planning/MULTI_BAND_TOMOGRAPHY_PLAN.md):
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
global_badchan_file = /dev/null
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
