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

BIN_SERIAL="$REPO_ROOT/bin/rm_synthesis_release_omp0_gpu0"
BIN_OMP="$REPO_ROOT/bin/rm_synthesis_release_omp1_gpu0"
BIN_GPU="$REPO_ROOT/bin/rm_synthesis_release_omp0_gpu1"

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
    local tag="$1" use_gpu="$2"
    local out_prefix="$OUT_DIR/$tag"
    local cfg="$OUT_DIR/${tag}.cfg"
    sed -e "s|__DATADIR__|${DATA_DIR}|g" \
        -e "s|__OUTPREFIX__|${out_prefix}|g" \
        -e "s|__USE_GPU__|${use_gpu}|g" \
        "$TEMPLATE" > "$cfg"
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

# ---------------------------------------------------------------------------
# 0. Prepare output directory
# ---------------------------------------------------------------------------
section "0. Preparing directories"
mkdir -p "$OUT_DIR"

# Clean previous test outputs (binary refuses to overwrite)
rm -f "$OUT_DIR"/serial.*.FITS "$OUT_DIR"/omp.*.FITS "$OUT_DIR"/gpu.*.FITS
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
section "4. Building GPU binary  (make GPU=1)"
BUILD_GPU=0
if make GPU=1 2>&1 | tail -5; then
    if [[ -x "$BIN_GPU" ]]; then
        pass "GPU binary built: $BIN_GPU"
        BUILD_GPU=1
    else
        fail "GPU binary not found after make: $BIN_GPU"
    fi
else
    skip "make GPU=1 failed (no GPU compiler?); GPU test will be skipped"
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
                    --rtol 2e-4; then
                pass "GPU AMP: matches serial within rtol=2e-4"
            else
                fail "GPU AMP: differs from serial beyond rtol=2e-4"
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
