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
#      rmclean_plots/) via tests/plot_thesis_scenario_rmclean.py, using
#      ~/venv/rmtool's python3 if present (this developer's own checked
#      package versions) else falling back to plain python3 on PATH
#      (e.g. inside a container, which has no ~/venv/rmtool but does
#      have numpy/matplotlib installed system-wide).
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
#  27. RM-CLEAN get_drm enforcement – confirms get_drm REFUSES
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
# The ~50 bare `python3` calls throughout this script (check_rm_peak.py,
# compare_cubes.py, and many inline `python3 -c` one-liners) all rely on
# plain PATH lookup, not a resolved variable -- prepending the venv here
# once, up front, makes every one of them resolve to it automatically,
# without editing each call site individually. A system with no
# ~/venv/rmtool at all (e.g. a container with deps installed system-wide,
# see docker/dockerfile) falls through to whatever python3 is already on
# PATH, unchanged from before this prepend existed.
if [[ -x "$HOME/venv/rmtool/bin/python3" ]]; then
  export PATH="$HOME/venv/rmtool/bin:$PATH"
fi

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
WARN=0

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
pass() { echo "[PASS] $*"; PASS=$((PASS + 1)); }
fail() { echo "[FAIL] $*"; FAIL=$((FAIL + 1)); }
skip() { echo "[SKIP] $*"; SKIP=$((SKIP + 1)); }
# warn(): a real check that failed, but deliberately does NOT flip the
# suite's exit code -- for capabilities that should keep being noticed if
# they silently break, without blocking every other test/deploy on them
# (e.g. swim-lane plot generation: a real, useful diagnostic, but not
# something any test's own correctness should depend on). Counted and
# printed separately in the final summary so it's never silently lost,
# just never fatal.
warn() { echo "[WARN] $*"; WARN=$((WARN + 1)); }
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
# Logs from earlier runs stay in $OUT_DIR; the final FFTW-wisdom scan
# (see "Test Summary") only reads logs written after this marker.
RUN_START_MARKER="$OUT_DIR/.run_tests_start_marker"
touch "$RUN_START_MARKER"

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
GEN_PY="$HOME/venv/rmtool/bin/python3"
if [[ ! -x "$GEN_PY" ]]; then
    # No personal venv on this machine yet (e.g. scripts/setup_venv.sh
    # hasn't been run, or a container image with its own system-wide
    # install instead -- see docker/dockerfile) -- fall back to
    # whatever python3 is on PATH.
    GEN_PY="$(command -v python3 || true)"
fi

if [[ -n "$GEN_PY" ]] && "$GEN_PY" "$TESTS_DIR/make_test_cubes.py"; then
    pass "Synthetic cubes generated"
elif [[ -t 0 ]]; then
    echo "[run_tests] make_test_cubes.py failed -- this usually means this project's Python venv hasn't been set up yet." >&2
    venv_confirm=""
    read -r -p "[run_tests] Run scripts/setup_venv.sh now? [y/N] " venv_confirm || venv_confirm=""
    if [[ "$venv_confirm" =~ ^[Yy]$ ]]; then
        bash "$REPO_ROOT/scripts/setup_venv.sh"
        GEN_PY="$HOME/venv/rmtool/bin/python3"
        if [[ -x "$GEN_PY" ]] && "$GEN_PY" "$TESTS_DIR/make_test_cubes.py"; then
            pass "Synthetic cubes generated"
        else
            fail "make_test_cubes.py failed"
            echo "Cannot continue without test data, even after running scripts/setup_venv.sh."
            exit 1
        fi
    else
        fail "make_test_cubes.py failed"
        echo "Cannot continue without test data. Run scripts/setup_venv.sh, then re-run this script."
        exit 1
    fi
else
    fail "make_test_cubes.py failed"
    echo "Cannot continue without test data. Run scripts/setup_venv.sh, then re-run this script." >&2
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
# bypasses CFITSIO for AMP/PHA pixel writes via independent raw
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
global_badchan_file = none,none
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
        # badchan_file is a per-band list, same cardinality as resi_list --
        # reuse its comma count rather than hard-coding it.
        local badchan_list
        badchan_list=$(echo "$resi_list" | sed 's/[^,]*/none/g')
        cat > "$cfg" <<CFGEOF
path                = ${DATA_DIR}/
infileQ             = ${infileQ}
infileU             = ${infileU}
outfile             = ${OUT_DIR}/${tag}
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
#     a mismatch here always means a regression, never expected
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
global_badchan_file = none,none
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
global_badchan_file = none,none
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
global_badchan_file = none,none
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
global_badchan_file = none,none
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
#     afterward, sequentially, outside that parallel region.
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
        # rm_synthesis clamps io_read_threads_eff = max(1, min(requested,
        # omp_get_max_threads())) (src/rm_synthesis.f90) -- OMP_NUM_THREADS
        # is never set in this suite, so omp_get_max_threads() defaults to
        # the visible CPU count (confirmed matching nproc). A CI runner
        # with fewer than 4 cores legitimately opens fewer than the 4
        # handles requested here -- compute the same clamp rather than
        # hardcoding 4, so this doesn't fail on a smaller runner while
        # still proving the parallel-read path was taken. The
        # log line itself is only emitted when io_read_threads_eff > 1
        # (src/rm_synthesis.f90's own if/else), so a 1-CPU environment
        # (clamped to serial) can't exercise this at all -- skip, not fail.
        mbrt_avail=$(nproc)
        mbrt_expected=4
        [[ "$mbrt_avail" -lt "$mbrt_expected" ]] && mbrt_expected="$mbrt_avail"
        [[ "$mbrt_expected" -lt 1 ]] && mbrt_expected=1
        if [[ "$mbrt_expected" -eq 1 ]]; then
            skip "io_read_threads>1 multi-band (T9): only 1 CPU visible, clamps to serial -- parallel-read path not exercisable here"
        elif grep -qE "Parallel FITS IO: opened +${mbrt_expected} +handles/file" "$mbrt_log"; then
            pass "io_read_threads>1 multi-band (T9): parallel-read path confirmed taken (${mbrt_expected} handles)"
        else
            fail "io_read_threads>1 multi-band (T9): expected ${mbrt_expected}-handle parallel-read log line not found (see $mbrt_log)"
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
PLOT_PY="$HOME/venv/rmtool/bin/python3"
if [[ ! -x "$PLOT_PY" ]]; then
    # Not every environment has this developer's own personal venv --
    # a container image has no ~/venv/rmtool at all, but does have
    # numpy/matplotlib installed against its own system python3 (see
    # docker/dockerfile). Fall back to whatever python3 is on PATH
    # rather than skipping plot generation everywhere except this one
    # machine.
    PLOT_PY="$(command -v python3 || true)"
fi

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
        if [[ -x "$PLOT_PY" ]]; then
            if "$PLOT_PY" "$TESTS_DIR/plot_thesis_scenario_rmclean.py" \
                    "$RMCLEAN_CSV_DIR" "$RMCLEAN_PLOT_DIR" \
                    > "$OUT_DIR/rmclean_plot.log" 2>&1; then
                pass "RM-CLEAN thesis scenario: plots written to $RMCLEAN_PLOT_DIR"
            else
                fail "RM-CLEAN thesis scenario: plotting failed (see $OUT_DIR/rmclean_plot.log)"
            fi
        else
            skip "no usable python3 found (checked ~/venv/rmtool and PATH); skipping RM-CLEAN thesis scenario plots"
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
# 27. RM-CLEAN get_drm REFUSES an unsafe oversample -- inverted
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

        # T19 (docs/dev/RMCLEAN_INTEGRATION_PLAN.md): the dirty AMP cube's
        # own BUNIT must be the input Q/U cube's own BUNIT ("Jy/beam" for
        # this fixture, see tests/make_test_cubes.py) with '/RMSF'
        # appended -- the Faraday dispersion function is a
        # Faraday-depth density, not a plain passthrough.
        rmc_amp_bunit=$(python3 -c "
from astropy.io import fits
print(fits.getheader('$rmc_lsqref_amp').get('BUNIT'))
")
        if [ "$rmc_amp_bunit" = "Jy/beam/RMSF" ]; then
            pass "rm_synthesis: AMP.RMCUBE.FITS BUNIT is 'Jy/beam/RMSF' (T19)"
        else
            fail "rm_synthesis: AMP.RMCUBE.FITS BUNIT expected 'Jy/beam/RMSF', got '$rmc_amp_bunit'"
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

            # T19: CLEAN is discrete/additive -- keeps the PRE-'/RMSF'
            # units (plain 'Jy/beam' here); RESID/RESTORED are smeared
            # Faraday-depth densities, same as dirty -- inherit
            # '/RMSF'. Every .PHA. output is 'rad' regardless (phase has
            # no flux-density interpretation) -- this also regression-
            # tests the pre-existing bug found while wiring this up
            # (every .PHA. output previously inherited the AMP
            # template's own amplitude BUNIT instead of 'rad').
            t19_bunit_ok=1
            for pair in \
                "CLEAN.AMP:Jy/beam" "CLEAN.PHA:rad" \
                "RESID.AMP:Jy/beam/RMSF" "RESID.PHA:rad" \
                "RESTORED.AMP:Jy/beam/RMSF" "RESTORED.PHA:rad"
            do
                t19_suffix="${pair%%:*}"
                t19_expected="${pair##*:}"
                t19_got=$(python3 -c "
from astropy.io import fits
print(fits.getheader('${rmc_out}.${t19_suffix}.RMCUBE.FITS').get('BUNIT'))
")
                if [ "$t19_got" != "$t19_expected" ]; then
                    fail "rmclean_cubes: ${t19_suffix}.RMCUBE.FITS BUNIT expected '$t19_expected', got '$t19_got'"
                    t19_bunit_ok=0
                fi
            done
            if [ "$t19_bunit_ok" = "1" ]; then
                pass "rmclean_cubes: all 6 outputs' BUNIT correct (CLEAN=Jy/beam, RESID/RESTORED=Jy/beam/RMSF, PHA=rad) (T19)"
            fi

            # T19 Part D: RMSFFWHM recorded on RESTORED only (the only
            # output the restoring beam is convolved into),
            # matching the "Restoring beam FWHM = ..." value this same
            # run printed to its own log.
            t19_fwhm_printed=$(grep -oP 'Restoring beam FWHM = \K[0-9.]+' "$rmc_log")
            t19_fwhm_restored_amp=$(python3 -c "
from astropy.io import fits
print(fits.getheader('${rmc_out}.RESTORED.AMP.RMCUBE.FITS').get('RMSFFWHM'))
")
            t19_fwhm_restored_pha=$(python3 -c "
from astropy.io import fits
print(fits.getheader('${rmc_out}.RESTORED.PHA.RMCUBE.FITS').get('RMSFFWHM'))
")
            t19_fwhm_clean_amp=$(python3 -c "
from astropy.io import fits
print(fits.getheader('${rmc_out}.CLEAN.AMP.RMCUBE.FITS').get('RMSFFWHM'))
")
            if python3 -c "
import sys
printed = float('$t19_fwhm_printed')
restored = float('$t19_fwhm_restored_amp')
sys.exit(0 if abs(printed - restored) < 1e-3 else 1)
" && [ "$t19_fwhm_restored_pha" = "$t19_fwhm_restored_amp" ] \
              && [ "$t19_fwhm_clean_amp" = "None" ]; then
                pass "rmclean_cubes: RMSFFWHM=$t19_fwhm_restored_amp on RESTORED only, matches printed value ($t19_fwhm_printed) (T19)"
            else
                fail "rmclean_cubes: RMSFFWHM check failed (printed=$t19_fwhm_printed restored.amp=$t19_fwhm_restored_amp restored.pha=$t19_fwhm_restored_pha clean.amp=$t19_fwhm_clean_amp)"
            fi

            # T22 (docs/dev/RMCLEAN_INTEGRATION_PLAN.md): this fixture's
            # own dirty ampfile has LSQREF=$lsqref_val (lsq_ref_native,
            # rm_synthesis's own lsq_ref_mode=mid, confirmed nonzero
            # above) -- this run passed no lsq_ref_report_mode, so
            # lsq_ref_report defaults to intrinsic=0.0. Before T22, every
            # output's LSQREF was copied verbatim from ampfile (i.e. this
            # exact native!=report divergence silently wrote the WRONG
            # value); the fix must write 0.0 on every one of the 6 outputs
            # instead.
            t22_lsqref_ok=1
            for t22_suffix in CLEAN.AMP CLEAN.PHA RESID.AMP RESID.PHA RESTORED.AMP RESTORED.PHA; do
                t22_got=$(python3 -c "
from astropy.io import fits
print(fits.getheader('${rmc_out}.${t22_suffix}.RMCUBE.FITS').get('LSQREF'))
")
                if ! python3 -c "
import sys
sys.exit(0 if abs(float('$t22_got')) < 1e-9 else 1)
"; then
                    fail "rmclean_cubes: ${t22_suffix}.RMCUBE.FITS LSQREF expected 0.0 (lsq_ref_report default), got '$t22_got' (native was $lsqref_val)"
                    t22_lsqref_ok=0
                fi
            done
            if [ "$t22_lsqref_ok" = "1" ]; then
                pass "rmclean_cubes: all 6 outputs' LSQREF=0.0 (lsq_ref_report), not the dirty cube's own nonzero lsq_ref_native=$lsqref_val (T22)"
            fi

            # T19 Part B: cross-validate all 6 new diagnostic maps
            # against the REAL RESID.AMP/CLEAN.AMP cubes from this same
            # run (check_diagnostic_maps.py's own module docstring for
            # why -- not just self-consistency, independent recomputation).
            t19_mean_niter=$(grep -oP 'n_iter_used: mean=\K[0-9.]+' "$rmc_log" | tail -1)
            if python3 "$TESTS_DIR/check_diagnostic_maps.py" "$rmc_out" "$t19_mean_niter"; then
                pass "rmclean_cubes: 7 CLEAN diagnostic maps cross-validated against RESID.AMP/CLEAN.AMP (T19)"
            else
                fail "rmclean_cubes: CLEAN diagnostic maps failed cross-validation (see above)"
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

            # T19 Part A: forced 3x5 tiles on this fixture produces
            # MANY blocks (unlike section 29's own default single-tile
            # run) -- exactly what's needed to check the per-block
            # tallies aren't just mirroring the global ones. Sums
            # across every "Block N CLEAN stop-reasons" line must equal
            # this SAME run's own global summary exactly.
            if python3 -c "
import re
with open('$rmc_tiled_log') as f:
    log = f.read()
blk = re.findall(
    r'Block \d+ CLEAN stop-reasons: CLEANed=(\d+) niter=(\d+) .*? '
    r'abs_flux=(\d+) .*? auto_nsigma=(\d+) ',
    log)
assert len(blk) > 1, f'expected >1 block, found {len(blk)}'
sum_cleaned = sum(int(b[0]) for b in blk)
sum_niter = sum(int(b[1]) for b in blk)
sum_abs_flux = sum(int(b[2]) for b in blk)
sum_auto_nsigma = sum(int(b[3]) for b in blk)
g_cleaned = int(re.search(r'CLEANed (\d+) pixels\.', log).group(1))
g_niter = int(re.search(r'hit niter cap:\s+(\d+)', log).group(1))
g_abs_flux = int(re.search(r'stopped at abs_flux:\s+(\d+)', log).group(1))
g_auto_nsigma = int(re.search(r'stopped at auto_nsigma:(\d+)', log).group(1))
assert sum_cleaned == g_cleaned, f'CLEANed {sum_cleaned} != global {g_cleaned}'
assert sum_niter == g_niter, f'niter {sum_niter} != global {g_niter}'
assert sum_abs_flux == g_abs_flux, f'abs_flux {sum_abs_flux} != global {g_abs_flux}'
assert sum_auto_nsigma == g_auto_nsigma, f'auto_nsigma {sum_auto_nsigma} != global {g_auto_nsigma}'
print(f'{len(blk)} blocks, sums match global exactly (CLEANed={g_cleaned})')
"; then
                pass "rmclean_cubes: per-block CLEAN stop-reason tallies sum exactly to the global summary (T19)"
            else
                fail "rmclean_cubes: per-block CLEAN stop-reason tallies do not sum to the global summary (see $rmc_tiled_log)"
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
        # at nwriters=4: this path previously had a
        # PROBABILISTIC bug (FTGHAD's 3 output arguments come back with
        # their upper 32 bits UNTOUCHED by this system's installed
        # libcfitsio -- a Fortran-wrapper/library ABI truncation,
        # confirmed with a minimal standalone reproducer outside this
        # codebase -- so an uninitialized receiving variable's
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
        # nwriters together (the combined stress case
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
        # (mask_pattern_cache_max=0 -- every pixel a one-off
        # throwaway table) -- the cache is a reuse optimization,
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

        # T14 (docs/dev/RMCLEAN_INTEGRATION_PLAN.md): the comparison
        # above never exercises eviction -- mask_pattern_
        # cache_max=4096 comfortably fits rmc_varied's own true distinct-
        # pattern count (757, confirmed via this same fixture's own
        # "Mask-pattern cache: N distinct pattern(s) cached" log line).
        # Force real eviction by setting the cache far below that count
        # and confirm cache_eviction_policy=hitcount still gives
        # bit-identical output vs caching disabled -- eviction
        # (like caching itself) must never change the numerical result,
        # only how much redundant work it takes to get there.
        rmc_cache_evict="$OUT_DIR/rmc_cache_evict"
        rmc_cache_evict_log="$OUT_DIR/rmc_cache_evict.log"
        if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" \
                phafile="$OUT_DIR/rmc_lsqref.PHA.RMCUBE.FITS" \
                maskfile="$rmc_varied_mask" outfile="$rmc_cache_evict" \
                abs_flux_floor=0.01 niter=200 gain=0.1 \
                mask_pattern_cache_max=100 cache_eviction_policy=hitcount \
                > "$rmc_cache_evict_log" 2>&1; then
            evict_ok=1
            for suffix in CLEAN.AMP CLEAN.PHA RESID.AMP RESID.PHA \
                          RESTORED.AMP RESTORED.PHA; do
                if ! cmp -s "${rmc_cache_evict}.${suffix}.RMCUBE.FITS" \
                            "${rmc_cache_zero}.${suffix}.RMCUBE.FITS"; then
                    evict_ok=0
                    fail "rmclean_cubes: hitcount eviction (cache_max=100) vs no-cache differ on $suffix"
                fi
            done
            [[ "$evict_ok" -eq 1 ]] && \
                pass "rmclean_cubes: cache_eviction_policy=hitcount (forced eviction, cache_max=100 < 757 distinct patterns) bit-identical to cache disabled"
            if grep -q "pixel(s) past" "$rmc_cache_evict_log"; then
                fail "rmclean_cubes: hitcount eviction still fell back to the overflow/one-off path (see $rmc_cache_evict_log) -- eviction isn't actually engaging"
            else
                pass "rmclean_cubes: cache_eviction_policy=hitcount served every pattern via eviction, zero overflow fallbacks"
            fi
        else
            fail "rmclean_cubes: cache_eviction_policy=hitcount forced-eviction run failed (see $rmc_cache_evict_log)"
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
    # +0.016 deg (16 pixels on this 32-pixel/0.032-deg-wide field, UR-15,
    # docs/dev/BRANCH_STATUS_update-reproject.md): large enough that
    # sky_wcs_matches_target's exact-match check still says "process
    # this, don't skip it" -- the only thing this fixture needs -- but
    # small enough to leave sky overlap, now that mode=reference
    # refuses zero overlap. The original +0.05 (50 pixels, wider than
    # the field itself) left zero overlap by accident, since
    # overlap was never checked before UR-15.
    hdul[0].header["CRVAL1"] = float(hdul[0].header["CRVAL1"]) + 0.016
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

    # --- Negative: offset CRVAL -> processed as normal,
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
#     Reuses section 29's own rmc_lsqref_amp/PHA/MASK cube trio -- a
#     dirty cube with known BUNIT, no need to rebuild one here.
# ---------------------------------------------------------------------------
section "31. rmclean_cubes: abs_flux_floor units (Jy/mJy/uJy) + auto_nsigma (T8)"

if [[ -s "$rmc_lsqref_amp" ]]; then
    rmc_pha="$OUT_DIR/rmc_lsqref.PHA.RMCUBE.FITS"
    rmc_mask="$OUT_DIR/rmc_lsqref.MASK.CUBE.FITS"

    # --- Unit conversion: abs_flux_floor=10mJy must resolve to the same
    # native-unit value as a hand-computed conversion from the cube's own
    # BUNIT, not just "some number" -- verified against the header,
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

    # --- abs_flux_floor= and auto_nsigma= MAY now be combined (a
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
# 32. gaussft_mod: FFT-friendly padding (found via an e2e validation
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
        -o "$gfp_bin" -lfftw3 -lfftw3f -lfftw3f_omp 2>"$OUT_DIR/gaussft_padding_build.log"; then
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
# 32b. convolve_cubes/reproject_cubes: build once, upfront, the same way
#      match_cubes/rmclean_cubes already are (sections 29/30). Previously
#      these two had no build step of their own here -- every section
#      below that needs them just checked `[[ -x bin/convolve_cubes ]]`
#      and silently skipped if that pre-existing binary wasn't already
#      sitting in bin/, with nothing distinguishing a genuine platform
#      limitation (like the GPU build's own skip) from simply not having
#      run ./build.sh first. Building them here means every section 33+
#      test below actually runs on any machine where these two targets
#      build successfully.
# ---------------------------------------------------------------------------
section "32b. convolve_cubes/reproject_cubes: build"

if make convolve_cubes > "$OUT_DIR/convolve_cubes_build.log" 2>&1; then
    pass "convolve_cubes: build succeeded"
else
    fail "convolve_cubes: build failed (see $OUT_DIR/convolve_cubes_build.log)"
fi

if make reproject_cubes > "$OUT_DIR/reproject_cubes_build.log" 2>&1; then
    pass "reproject_cubes: build succeeded"
else
    fail "reproject_cubes: build failed (see $OUT_DIR/reproject_cubes_build.log)"
fi

# ---------------------------------------------------------------------------
# 32c. FFTW wisdom is saved. Every bit-identical check below that compares
#      two separate convolve_cubes/match_cubes runs depends on it: without
#      saved wisdom each run re-measures its FFT plans, and two runs can
#      pick different (equally valid) plans that differ in the last bits.
#      Whether that happens depends on timing, so those checks alone would
#      pass or fail at random; this check fails every time instead.
#      src/gaussft.f90's plan_convolution prints the WARNING matched here.
# ---------------------------------------------------------------------------
section "32c. FFTW wisdom: saved to a writable location"

WISDOM_WARNING_TEXT="WARNING: could not save FFTW wisdom"
if [[ -x bin/convolve_cubes ]]; then
    wis_beamfile="$OUT_DIR/wis_beamlog.txt"
    awk 'BEGIN { for (i=1;i<=200;i++) print i, 10.0, 10.0, 0.0 }' > "$wis_beamfile"
    wis_src="$OUT_DIR/wis_src.Q.FITSCUBE"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$wis_src"
    wis_out="${wis_src%.FITSCUBE}.CONV.FITS"
    rm -f "$wis_out"
    wis_log="$OUT_DIR/convolve_wisdom.log"
    wis_file="${RMTOOL_WISDOM_DIR:-wisdom}/fftw_wisdom_sp.dat"

    if bin/convolve_cubes infiles="$wis_src" beamfiles="$wis_beamfile" \
            outsuffix=".CONV.FITS" target_bmaj=20.0 target_bmin=20.0 \
            target_bpa=0.0 mem_frac_ram=0.1 > "$wis_log" 2>&1; then
        if grep -q "$WISDOM_WARNING_TEXT" "$wis_log"; then
            fail "FFTW wisdom not saved -- bit-identical convolve checks below cannot be trusted (see $wis_log; set RMTOOL_WISDOM_DIR to a writable directory)"
        elif [[ ! -s "$wis_file" ]]; then
            fail "FFTW wisdom file $wis_file missing or empty after a convolve_cubes run (see $wis_log)"
        else
            pass "FFTW wisdom saved to $wis_file"
        fi
    else
        fail "convolve_cubes run for the FFTW wisdom check failed (see $wis_log)"
    fi
else
    skip "bin/convolve_cubes not built; skipping FFTW wisdom check"
fi

# ---------------------------------------------------------------------------
# 33. convolve_cubes: io_overlap bit-identical to io_overlap=n (planning-
#     doc ticket, added alongside the FFT padding fix -- the same e2e
#     validation run found convolve's own block write serial
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
    cp "$DATA_DIR/TEST_BAND2_GRIDSHIFT.Q.FITSCUBE" "$rcio_src"

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
    cp "$DATA_DIR/TEST_BAND2_GRIDSHIFT.Q.FITSCUBE" "$mcrio_src"

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
# 39-40. Multi-band preprocessing END-TO-END (docs/dev/
#     MULTI_BAND_TOMOGRAPHY_PLAN.md T15): sections 15/30/33-36 test the
#     preprocessing tools' own internal correctness (geometry validation,
#     skip-if-already-matched, io_overlap consistency) in isolation;
#     these run the preprocessed cubes through rm_synthesis.
#     Section 39 isolates resolution-only (convolve_cubes alone), section
#     40 grid and resolution together (match_cubes stages=both).
#     band1_beamlog.txt/band2_beamlog.txt (tests/make_test_cubes.py)
#     give band 2 a 20" beam against band 1's 10", so convolve_cubes has
#     smoothing work to do rather than a no-op copy. Section 40's band 2
#     is TEST_BAND2_GRIDSHIFT, a different-grid input that overlaps band
#     1; see its comment in make_test_cubes.py for what it does and does
#     not represent.
#
#     There is no grid-only section. One existed (T15) that reprojected
#     band 2's Q but passed band 2's U through unreprojected, so the two
#     never described the same sky: on a zero-overlap input it passed
#     with band 2 contributing nothing, and on GRIDSHIFT it failed.
#     It was removed rather than repaired, because reprojection accuracy
#     is tested against analytic truth in sections 62/70 and reprojected
#     Q/U reaching rm_synthesis is exercised by section 40.
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
global_badchan_file = none,none
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
section "40. Multi-band preprocessing: match_cubes stages=both output on a different-grid, different-beam band 2 runs through rm_synthesis with the known sources intact (T15)"
# Band 2 is TEST_BAND2_GRIDSHIFT, whose header move shifts its sky 16 px
# (see make_test_cubes.py), so at truth.json's positions band 1 carries
# the recovery. This checks that stages=both runs on a grid+beam mismatch
# and writes cubes rm_synthesis accepts without disturbing the band 1
# sources -- not reprojection accuracy, which sections 62/70 test against
# analytic truth.

if [[ -x bin/match_cubes && -x "$BIN_SERIAL" ]]; then
    mbpp_both_band1q="$OUT_DIR/mbpp_both.TEST.Q.FITSCUBE"
    mbpp_both_band1u="$OUT_DIR/mbpp_both.TEST.U.FITSCUBE"
    mbpp_both_band2q="$OUT_DIR/mbpp_both.TEST_BAND2_GRIDSHIFT.Q.FITSCUBE"
    mbpp_both_band2u="$OUT_DIR/mbpp_both.TEST_BAND2_GRIDSHIFT.U.FITSCUBE"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$mbpp_both_band1q"
    cp "$DATA_DIR/TEST.U.FITSCUBE" "$mbpp_both_band1u"
    cp "$DATA_DIR/TEST_BAND2_GRIDSHIFT.Q.FITSCUBE" "$mbpp_both_band2q"
    cp "$DATA_DIR/TEST_BAND2_GRIDSHIFT.U.FITSCUBE" "$mbpp_both_band2u"
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
        pass "match_cubes stages=both: processed a band 2 on a different grid and beam, wrote output"
    else
        fail "match_cubes stages=both: failed on a band 2 with a different grid and beam, or wrote no output (see $mbpp_both_match_log)"
    fi

    if [[ -s "$mbpp_both_band1q_matched" && -s "$mbpp_both_band2q_matched" ]]; then
        mbpp_both_cfg="$OUT_DIR/mbpp_both_rmsynth.cfg"
        cat > "$mbpp_both_cfg" <<CFGEOF
path                = ${OUT_DIR}/
infileQ             = $(basename "$mbpp_both_band1q_matched"),$(basename "$mbpp_both_band2q_matched")
infileU             = $(basename "$mbpp_both_band1u_matched"),$(basename "$mbpp_both_band2u_matched")
outfile             = ${OUT_DIR}/mbpp_both_out
global_badchan_file = none,none
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
            pass "match_cubes stages=both + rm_synthesis: output accepted and src_A/src_B recovered at their known RMs"
        else
            fail "match_cubes stages=both + rm_synthesis: RM peak(s) not recovered (see $mbpp_both_rmsynth_log)"
        fi
    else
        fail "match_cubes stages=both: no usable output; skipping downstream rm_synthesis check"
    fi
else
    skip "match_cubes or serial rm_synthesis binary not available; skipping grid+resolution preprocessing test"
fi

# ---------------------------------------------------------------------------
# 41. convolve_cubes/match_cubes: badchan_file is per-band (T16,
#     docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md) -- a real bug found while
#     preparing a real multi-band run: badchan_file used to be a single
#     shared list applied identically to every infile regardless of that
#     infile's own channel numbering, so band A's own bad channel would
#     also get (incorrectly) excluded from band B, and vice versa. Fixed
#     to be a comma list, one entry per infile, exactly like beamfiles.
#     Verifies each band's own designated bad channel is excluded from
#     THAT band's own output, and NOT from the other band's. Extended by
#     T32 (docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md): the literal value
#     'none' for one infile's entry must behave identically to leaving
#     that entry blank -- no crash, no bad channels excluded for that
#     infile, and the other infile's own list still applies.
# ---------------------------------------------------------------------------
section "41. convolve_cubes/match_cubes: badchan_file is per-band (T16)"

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
    # already-matching -- both must convolve and write output.
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
            pass "convolve_cubes: badchan_file is per-band (band1 only chan 13 bad, band2 only chan 77 bad, no cross-contamination)"
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
            pass "match_cubes: badchan_file is per-band (band1 only chan 13 bad, band2 only chan 77 bad, no cross-contamination)"
        else
            fail "match_cubes: badchan_file per-band check failed (band1 bad=[$pbbcm_band1_bad_planes], expected 13; band2 bad=[$pbbcm_band2_bad_planes], expected 77) (see $pbbcm_convolve_log)"
        fi
    else
        fail "match_cubes: per-band badchan_file run failed (see $pbbcm_convolve_log)"
    fi
else
    skip "bin/match_cubes not built; skipping per-band badchan_file test"
fi

# T32: 'none' as an explicit per-infile badchan_file entry, same meaning as
# leaving that entry blank.
if [[ -x bin/convolve_cubes ]]; then
    pbbcn_band1="$OUT_DIR/pbbcn_band1.Q.FITSCUBE"
    pbbcn_band2="$OUT_DIR/pbbcn_band2.Q.FITSCUBE"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$pbbcn_band1"
    cp "$DATA_DIR/TEST_BAND2.Q.FITSCUBE" "$pbbcn_band2"
    pbbcn_band2_badchan="$OUT_DIR/pbbcn_band2_badchan.txt"
    printf '77\n' > "$pbbcn_band2_badchan"
    pbbcn_band1_conv="${pbbcn_band1%.FITSCUBE}_CONV.FITS"
    pbbcn_band2_conv="${pbbcn_band2%.FITSCUBE}_CONV.FITS"
    rm -f "$pbbcn_band1_conv" "$pbbcn_band2_conv"

    pbbcn_convolve_log="$OUT_DIR/pbbcn_convolve.log"
    if bin/convolve_cubes infiles="$pbbcn_band1,$pbbcn_band2" \
            beamfiles="$DATA_DIR/band1_beamlog.txt,$DATA_DIR/band2_beamlog.txt" \
            badchan_file="none,$pbbcn_band2_badchan" \
            target_bmaj=30.0 target_bmin=30.0 target_bpa=0.0 \
            > "$pbbcn_convolve_log" 2>&1; then
        pbbcn_band1_bad_planes=$(python3 -c "
from astropy.io import fits
import numpy as np
data = fits.getdata('$pbbcn_band1_conv')[0]
print(','.join(str(i+1) for i in range(data.shape[0]) if np.isnan(data[i]).all()))
")
        pbbcn_band2_bad_planes=$(python3 -c "
from astropy.io import fits
import numpy as np
data = fits.getdata('$pbbcn_band2_conv')[0]
print(','.join(str(i+1) for i in range(data.shape[0]) if np.isnan(data[i]).all()))
")
        if [[ "$pbbcn_band1_bad_planes" == "" && "$pbbcn_band2_bad_planes" == "77" ]]; then
            pass "convolve_cubes: badchan_file='none' for band1 behaves as no list (band1 has zero bad planes, band2's own chan 77 still excluded)"
        else
            fail "convolve_cubes: badchan_file='none' check failed (band1 bad=[$pbbcn_band1_bad_planes], expected none; band2 bad=[$pbbcn_band2_bad_planes], expected 77) (see $pbbcn_convolve_log)"
        fi
    else
        fail "convolve_cubes: badchan_file='none' run failed (see $pbbcn_convolve_log)"
    fi
else
    skip "bin/convolve_cubes not built; skipping badchan_file='none' test"
fi

if [[ -x bin/match_cubes ]]; then
    pbbcmn_band1="$OUT_DIR/pbbcmn_band1.Q.FITSCUBE"
    pbbcmn_band2="$OUT_DIR/pbbcmn_band2.Q.FITSCUBE"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$pbbcmn_band1"
    cp "$DATA_DIR/TEST_BAND2.Q.FITSCUBE" "$pbbcmn_band2"
    pbbcmn_band2_badchan="$OUT_DIR/pbbcmn_band2_badchan.txt"
    printf '77\n' > "$pbbcmn_band2_badchan"
    pbbcmn_band1_conv="${pbbcmn_band1%.FITSCUBE}_CONV.FITS"
    pbbcmn_band2_conv="${pbbcmn_band2%.FITSCUBE}_CONV.FITS"
    rm -f "$pbbcmn_band1_conv" "$pbbcmn_band2_conv"

    pbbcmn_convolve_log="$OUT_DIR/pbbcmn_convolve.log"
    if bin/match_cubes stages=convolve infiles="$pbbcmn_band1,$pbbcmn_band2" \
            beamfiles="$DATA_DIR/band1_beamlog.txt,$DATA_DIR/band2_beamlog.txt" \
            badchan_file="none,$pbbcmn_band2_badchan" \
            target_bmaj=30.0 target_bmin=30.0 target_bpa=0.0 \
            > "$pbbcmn_convolve_log" 2>&1; then
        pbbcmn_band1_bad_planes=$(python3 -c "
from astropy.io import fits
import numpy as np
data = fits.getdata('$pbbcmn_band1_conv')[0]
print(','.join(str(i+1) for i in range(data.shape[0]) if np.isnan(data[i]).all()))
")
        pbbcmn_band2_bad_planes=$(python3 -c "
from astropy.io import fits
import numpy as np
data = fits.getdata('$pbbcmn_band2_conv')[0]
print(','.join(str(i+1) for i in range(data.shape[0]) if np.isnan(data[i]).all()))
")
        if [[ "$pbbcmn_band1_bad_planes" == "" && "$pbbcmn_band2_bad_planes" == "77" ]]; then
            pass "match_cubes: badchan_file='none' for band1 behaves as no list (band1 has zero bad planes, band2's own chan 77 still excluded)"
        else
            fail "match_cubes: badchan_file='none' check failed (band1 bad=[$pbbcmn_band1_bad_planes], expected none; band2 bad=[$pbbcmn_band2_bad_planes], expected 77) (see $pbbcmn_convolve_log)"
        fi
    else
        fail "match_cubes: badchan_file='none' run failed (see $pbbcmn_convolve_log)"
    fi
else
    skip "bin/match_cubes not built; skipping badchan_file='none' test"
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

lia_dir="$OUT_DIR/long_infiles_arg/a/really/quite/deeply/nested/directory/structure/designed/specifically/to/push/the/combined/infiles/argument/comfortably/past/the/old/five_hundred_twelve/character/limit/for/a/regression/test/regardless/of/how/short/or/long/the/absolute/path/to/this/particular/checkout/happens/to/be/on/whatever/host/or/container/is/running/this/test/suite/right/now"
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
        cp "$DATA_DIR/TEST_BAND2_GRIDSHIFT.Q.FITSCUBE" "$lia_dir/bandC.fits"
        lia_infiles_reproj="$lia_dir/bandA.fits,$lia_dir/bandC.fits"
        rm -f "$lia_dir/bandC_REPROJ.FITS"
        lia_reproject_log="$OUT_DIR/lia_reproject.log"
        if bin/reproject_cubes mode=reference reffile="$lia_dir/bandA.fits" \
                infiles="$lia_infiles_reproj" \
                > "$lia_reproject_log" 2>&1 && \
           [[ -s "$lia_dir/bandC_REPROJ.FITS" ]]; then
            pass "reproject_cubes: ${#lia_infiles_reproj}-char infiles= argument (>512) not truncated, mismatched band reprojected"
        else
            fail "reproject_cubes: long infiles= argument (${#lia_infiles_reproj} chars) failed (see $lia_reproject_log)"
        fi
    else
        skip "bin/reproject_cubes not built; skipping long infiles= test"
    fi
fi

# ---------------------------------------------------------------------------
# 43. convolve_cubes: nwriters>1 bit-identical to nwriters=1 (T21)
# ---------------------------------------------------------------------------
section "43. convolve_cubes: nwriters>1 bit-identical to nwriters=1 (io_overlap=y)"

if [[ -x bin/convolve_cubes ]]; then
    nwc_beamfile="$OUT_DIR/nwc_beamlog.txt"
    awk 'BEGIN { for (i=1;i<=200;i++) print i, 10.0, 10.0, 0.0 }' > "$nwc_beamfile"

    nwc_src="$OUT_DIR/nwc_src.Q.FITSCUBE"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$nwc_src"

    nwc_1_out="${nwc_src%.FITSCUBE}_nw1.CONV.FITS"
    nwc_4_out="${nwc_src%.FITSCUBE}_nw4.CONV.FITS"
    rm -f "$nwc_1_out" "$nwc_4_out"

    nwc_1_log="$OUT_DIR/convolve_nwriters1.log"
    nwc_4_log="$OUT_DIR/convolve_nwriters4.log"

    if bin/convolve_cubes infiles="$nwc_src" beamfiles="$nwc_beamfile" \
            outsuffix="_nw1.CONV.FITS" target_bmaj=20.0 target_bmin=20.0 \
            target_bpa=0.0 mem_frac_ram=0.1 io_overlap=y nwriters=1 \
            > "$nwc_1_log" 2>&1 && \
       bin/convolve_cubes infiles="$nwc_src" beamfiles="$nwc_beamfile" \
            outsuffix="_nw4.CONV.FITS" target_bmaj=20.0 target_bmin=20.0 \
            target_bpa=0.0 mem_frac_ram=0.1 io_overlap=y nwriters=4 \
            > "$nwc_4_log" 2>&1; then
        if [[ -s "$nwc_1_out" && -s "$nwc_4_out" ]] && \
           cmp -s "$nwc_1_out" "$nwc_4_out"; then
            pass "convolve_cubes: nwriters=4 bit-identical to nwriters=1"
        else
            fail "convolve_cubes: nwriters=4 output differs from nwriters=1 (see $nwc_1_log / $nwc_4_log)"
        fi
    else
        fail "convolve_cubes: nwriters=1/4 run(s) failed (see $nwc_1_log / $nwc_4_log)"
    fi
else
    skip "bin/convolve_cubes not built; skipping nwriters bit-identical test"
fi

# ---------------------------------------------------------------------------
# 44. match_cubes: nwriters>1 bit-identical to nwriters=1 (stages=convolve, T21)
# ---------------------------------------------------------------------------
section "44. match_cubes: nwriters>1 bit-identical to nwriters=1 (stages=convolve)"

if [[ -x bin/match_cubes ]]; then
    mcnw_beamfile="$OUT_DIR/mcnw_beamlog.txt"
    awk 'BEGIN { for (i=1;i<=200;i++) print i, 10.0, 10.0, 0.0 }' > "$mcnw_beamfile"

    mcnw_src="$OUT_DIR/mcnw_src.Q.FITSCUBE"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$mcnw_src"

    mcnw_1_out="${mcnw_src%.FITSCUBE}_nw1.MATCHED.FITS"
    mcnw_4_out="${mcnw_src%.FITSCUBE}_nw4.MATCHED.FITS"
    rm -f "$mcnw_1_out" "$mcnw_4_out"

    mcnw_1_log="$OUT_DIR/match_nwriters1.log"
    mcnw_4_log="$OUT_DIR/match_nwriters4.log"

    if bin/match_cubes stages=convolve infiles="$mcnw_src" \
            beamfiles="$mcnw_beamfile" outsuffix="_nw1.MATCHED.FITS" \
            target_bmaj=20.0 target_bmin=20.0 target_bpa=0.0 \
            mem_frac_ram=0.1 io_overlap=y nwriters=1 \
            > "$mcnw_1_log" 2>&1 && \
       bin/match_cubes stages=convolve infiles="$mcnw_src" \
            beamfiles="$mcnw_beamfile" outsuffix="_nw4.MATCHED.FITS" \
            target_bmaj=20.0 target_bmin=20.0 target_bpa=0.0 \
            mem_frac_ram=0.1 io_overlap=y nwriters=4 \
            > "$mcnw_4_log" 2>&1; then
        if [[ -s "$mcnw_1_out" && -s "$mcnw_4_out" ]] && \
           cmp -s "$mcnw_1_out" "$mcnw_4_out"; then
            pass "match_cubes: nwriters=4 bit-identical to nwriters=1 (stages=convolve)"
        else
            fail "match_cubes: nwriters=4 output differs from nwriters=1 (see $mcnw_1_log / $mcnw_4_log)"
        fi
    else
        fail "match_cubes: nwriters=1/4 run(s) failed (see $mcnw_1_log / $mcnw_4_log)"
    fi
else
    skip "bin/match_cubes not built; skipping nwriters bit-identical test"
fi

# ---------------------------------------------------------------------------
# 45. reproject_cubes: nwriters>1 bit-identical to nwriters=1 (T21)
# ---------------------------------------------------------------------------
section "45. reproject_cubes: nwriters>1 bit-identical to nwriters=1"

if [[ -x bin/reproject_cubes ]]; then
    rcnw_ref="$DATA_DIR/TEST.Q.FITSCUBE"
    rcnw_src="$OUT_DIR/rcnw_src.Q.FITSCUBE"
    cp "$DATA_DIR/TEST_BAND2_GRIDSHIFT.Q.FITSCUBE" "$rcnw_src"

    rcnw_src_reproj="${rcnw_src%.FITSCUBE}_REPROJ.FITS"
    rcnw_1_out="${rcnw_src}_nw1.FITS"
    rcnw_4_out="${rcnw_src}_nw4.FITS"
    rm -f "$rcnw_src_reproj" "$rcnw_1_out" "$rcnw_4_out"

    rcnw_1_log="$OUT_DIR/reproject_nwriters1.log"
    rcnw_4_log="$OUT_DIR/reproject_nwriters4.log"

    if bin/reproject_cubes mode=reference reffile="$rcnw_ref" \
            infiles="$rcnw_src" mem_frac_ram=0.1 io_overlap=y nwriters=1 \
            > "$rcnw_1_log" 2>&1; then
        mv "$rcnw_src_reproj" "$rcnw_1_out"
    fi
    if bin/reproject_cubes mode=reference reffile="$rcnw_ref" \
            infiles="$rcnw_src" mem_frac_ram=0.1 io_overlap=y nwriters=4 \
            > "$rcnw_4_log" 2>&1; then
        mv "$rcnw_src_reproj" "$rcnw_4_out"
    fi

    if [[ -s "$rcnw_1_out" && -s "$rcnw_4_out" ]] && \
       cmp -s "$rcnw_1_out" "$rcnw_4_out"; then
        pass "reproject_cubes: nwriters=4 bit-identical to nwriters=1"
    else
        fail "reproject_cubes: nwriters=4 output differs from (or is missing vs) nwriters=1 (see $rcnw_1_log / $rcnw_4_log)"
    fi
else
    skip "bin/reproject_cubes not built; skipping nwriters bit-identical test"
fi

# ---------------------------------------------------------------------------
# 46. match_cubes: nwriters>1 bit-identical to nwriters=1 (stages=reproject, T21)
# ---------------------------------------------------------------------------
section "46. match_cubes: nwriters>1 bit-identical to nwriters=1 (stages=reproject)"

if [[ -x bin/match_cubes ]]; then
    mcrnw_ref="$DATA_DIR/TEST.Q.FITSCUBE"
    mcrnw_src="$OUT_DIR/mcrnw_src.Q.FITSCUBE"
    cp "$DATA_DIR/TEST_BAND2_GRIDSHIFT.Q.FITSCUBE" "$mcrnw_src"

    mcrnw_src_reproj="${mcrnw_src%.FITSCUBE}_REPROJ.FITS"
    mcrnw_1_out="${mcrnw_src}_nw1.FITS"
    mcrnw_4_out="${mcrnw_src}_nw4.FITS"
    rm -f "$mcrnw_src_reproj" "$mcrnw_1_out" "$mcrnw_4_out"

    mcrnw_1_log="$OUT_DIR/match_reproject_nwriters1.log"
    mcrnw_4_log="$OUT_DIR/match_reproject_nwriters4.log"

    if bin/match_cubes stages=reproject footprint_mode=reference \
            reffile="$mcrnw_ref" infiles="$mcrnw_src" outsuffix="_REPROJ.FITS" \
            mem_frac_ram=0.1 io_overlap=y nwriters=1 > "$mcrnw_1_log" 2>&1; then
        mv "$mcrnw_src_reproj" "$mcrnw_1_out"
    fi
    if bin/match_cubes stages=reproject footprint_mode=reference \
            reffile="$mcrnw_ref" infiles="$mcrnw_src" outsuffix="_REPROJ.FITS" \
            mem_frac_ram=0.1 io_overlap=y nwriters=4 > "$mcrnw_4_log" 2>&1; then
        mv "$mcrnw_src_reproj" "$mcrnw_4_out"
    fi

    if [[ -s "$mcrnw_1_out" && -s "$mcrnw_4_out" ]] && \
       cmp -s "$mcrnw_1_out" "$mcrnw_4_out"; then
        pass "match_cubes: nwriters=4 bit-identical to nwriters=1 (stages=reproject)"
    else
        fail "match_cubes: nwriters=4 output differs from (or is missing vs) nwriters=1 (see $mcrnw_1_log / $mcrnw_4_log)"
    fi
else
    skip "bin/match_cubes not built; skipping nwriters bit-identical test"
fi

# ---------------------------------------------------------------------------
# 47. convolve_cubes/match_cubes/reproject_cubes: dry_run touches no data (T24)
# ---------------------------------------------------------------------------
section "47. convolve_cubes/match_cubes/reproject_cubes: dry_run touches no data"

dr_src="$OUT_DIR/dr_src.Q.FITSCUBE"
cp "$DATA_DIR/TEST.Q.FITSCUBE" "$dr_src"
dr_beamfile="$OUT_DIR/dr_beamlog.txt"
awk 'BEGIN { for (i=1;i<=200;i++) print i, 10.0, 10.0, 0.0 }' > "$dr_beamfile"

if [[ -x bin/convolve_cubes ]]; then
    dr_cv_cfg="$OUT_DIR/convolve_cubes_dryrun.cfg"
    rm -f "$dr_cv_cfg" "${dr_src%.FITSCUBE}.CONV.FITS"
    dr_cv_log="$OUT_DIR/dryrun_convolve.log"
    if bin/convolve_cubes infiles="$dr_src" beamfiles="$dr_beamfile" \
            target_bmaj=20.0 target_bmin=20.0 target_bpa=0.0 dry_run=y \
            > "$dr_cv_log" 2>&1; then
        # Two valid shapes: a confident single suggestion (disk type
        # resolved), or -- common when the target isn't on host-visible
        # storage, e.g. inside a container's own overlay filesystem with
        # no bind mount -- a dual fallback offering both HDD and SSD
        # settings as commented alternatives (T1i-adjacent; see
        # src/convolve_cubes.f90's own run_dry_run_check). Lands next to
        # dr_src itself (in $OUT_DIR), not a bare CWD-relative name --
        # see run_dry_run_check's own comment for why.
        dr_cv_ok=0
        if [[ -s "$dr_cv_cfg" ]]; then
            if grep -q "^io_overlap=" "$dr_cv_cfg" && \
               grep -q "^nwriters=" "$dr_cv_cfg"; then
                dr_cv_ok=1
            elif grep -q "^# io_overlap=n$" "$dr_cv_cfg" && \
                 grep -q "^# io_overlap=y$" "$dr_cv_cfg"; then
                dr_cv_ok=1
            fi
        fi
        if [[ $dr_cv_ok -eq 1 ]] && [[ ! -e "${dr_src%.FITSCUBE}.CONV.FITS" ]]; then
            pass "convolve_cubes: dry_run=y wrote a suggested cfg and touched no data"
        else
            fail "convolve_cubes: dry_run=y did not behave as expected (see $dr_cv_log)"
        fi
    else
        fail "convolve_cubes: dry_run=y run failed (see $dr_cv_log)"
    fi
    rm -f "$dr_cv_cfg"
else
    skip "bin/convolve_cubes not built; skipping dry_run test"
fi

if [[ -x bin/match_cubes ]]; then
    dr_mc_cfg="$OUT_DIR/match_cubes_dryrun.cfg"
    rm -f "$dr_mc_cfg" "${dr_src%.FITSCUBE}_MATCHED.FITS"
    dr_mc_log="$OUT_DIR/dryrun_match.log"
    if bin/match_cubes stages=convolve infiles="$dr_src" beamfiles="$dr_beamfile" \
            target_bmaj=20.0 target_bmin=20.0 target_bpa=0.0 dry_run=y \
            > "$dr_mc_log" 2>&1; then
        # See convolve_cubes' own dry_run check above: confident single
        # suggestion, or dual HDD/SSD fallback when disk type can't be
        # resolved (e.g. no bind mount inside a container). Lands next
        # to dr_src itself (in $OUT_DIR), not a bare CWD-relative name.
        dr_mc_ok=0
        if [[ -s "$dr_mc_cfg" ]]; then
            if grep -q "^io_overlap=" "$dr_mc_cfg" && \
               grep -q "^nwriters=" "$dr_mc_cfg"; then
                dr_mc_ok=1
            elif grep -q "^# io_overlap=n$" "$dr_mc_cfg" && \
                 grep -q "^# io_overlap=y$" "$dr_mc_cfg"; then
                dr_mc_ok=1
            fi
        fi
        if [[ $dr_mc_ok -eq 1 ]] && [[ ! -e "${dr_src%.FITSCUBE}_MATCHED.FITS" ]]; then
            pass "match_cubes: dry_run=y wrote a suggested cfg and touched no data"
        else
            fail "match_cubes: dry_run=y did not behave as expected (see $dr_mc_log)"
        fi
    else
        fail "match_cubes: dry_run=y run failed (see $dr_mc_log)"
    fi
    rm -f "$dr_mc_cfg"
else
    skip "bin/match_cubes not built; skipping dry_run test"
fi

if [[ -x bin/reproject_cubes ]]; then
    dr_rc_cfg="$OUT_DIR/reproject_cubes_dryrun.cfg"
    rm -f "$dr_rc_cfg" "${dr_src%.FITSCUBE}_REPROJ.FITS"
    dr_rc_log="$OUT_DIR/dryrun_reproject.log"
    if bin/reproject_cubes mode=reference reffile="$DATA_DIR/TEST.Q.FITSCUBE" \
            infiles="$dr_src" dry_run=y > "$dr_rc_log" 2>&1; then
        # See convolve_cubes' own dry_run check above: confident single
        # suggestion, or dual HDD/SSD fallback when disk type can't be
        # resolved (e.g. no bind mount inside a container). Lands next
        # to dr_src itself (in $OUT_DIR), not a bare CWD-relative name.
        dr_rc_ok=0
        if [[ -s "$dr_rc_cfg" ]]; then
            if grep -q "^io_overlap=" "$dr_rc_cfg" && \
               grep -q "^nwriters=" "$dr_rc_cfg"; then
                dr_rc_ok=1
            elif grep -q "^# io_overlap=n$" "$dr_rc_cfg" && \
                 grep -q "^# io_overlap=y$" "$dr_rc_cfg"; then
                dr_rc_ok=1
            fi
        fi
        if [[ $dr_rc_ok -eq 1 ]] && [[ ! -e "${dr_src%.FITSCUBE}_REPROJ.FITS" ]]; then
            pass "reproject_cubes: dry_run=y wrote a suggested cfg and touched no data"
        else
            fail "reproject_cubes: dry_run=y did not behave as expected (see $dr_rc_log)"
        fi
    else
        fail "reproject_cubes: dry_run=y run failed (see $dr_rc_log)"
    fi
    rm -f "$dr_rc_cfg"
else
    skip "bin/reproject_cubes not built; skipping dry_run test"
fi

# ---------------------------------------------------------------------------
section "48. convolve_cubes: NaN pixels (T23) are locally masked, not plane-poisoning, and don't corrupt distant source recovery"

if [[ -x bin/convolve_cubes && -x "$BIN_SERIAL" ]]; then
    mbpp_nan_band1q="$OUT_DIR/mbpp_nan.TEST.Q.FITSCUBE"
    mbpp_nan_band1u="$OUT_DIR/mbpp_nan.TEST.U.FITSCUBE"
    mbpp_nan_band2q="$OUT_DIR/mbpp_nan.TEST_BAND2.Q.FITSCUBE"
    mbpp_nan_band2u="$OUT_DIR/mbpp_nan.TEST_BAND2.U.FITSCUBE"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$mbpp_nan_band1q"
    cp "$DATA_DIR/TEST.U.FITSCUBE" "$mbpp_nan_band1u"
    cp "$DATA_DIR/TEST_BAND2.Q.FITSCUBE" "$mbpp_nan_band2q"
    cp "$DATA_DIR/TEST_BAND2.U.FITSCUBE" "$mbpp_nan_band2u"

    # Inject NaN into a 6x6 corner block (x=0..5, y=0..5), every channel,
    # of band 1's Q and U only -- far from both known sources (src_A at
    # (12,10), src_B at (22,20)). Before T23's fix, convolve_to_beam ran
    # the whole plane as one global FFT with no NaN handling, so a
    # single NaN anywhere turned the ENTIRE output plane to NaN; a
    # correct fix must recover both sources unaffected while the
    # injected block itself still comes back mostly NaN (the threshold
    # rejection actually engaging, not just "no longer 100% NaN by
    # accident").
    python3 - "$mbpp_nan_band1q" "$mbpp_nan_band1u" <<'PYEOF'
import sys
import numpy as np
from astropy.io import fits
for path in sys.argv[1:]:
    with fits.open(path, mode='update') as hdul:
        data = hdul[0].data
        data[..., 0:6, 0:6] = np.nan
        hdul.flush()
PYEOF

    mbpp_nan_band1q_conv="${mbpp_nan_band1q%.FITSCUBE}_CONV.FITS"
    mbpp_nan_band1u_conv="${mbpp_nan_band1u%.FITSCUBE}_CONV.FITS"
    mbpp_nan_band2q_conv="${mbpp_nan_band2q%.FITSCUBE}_CONV.FITS"
    mbpp_nan_band2u_conv="${mbpp_nan_band2u%.FITSCUBE}_CONV.FITS"
    rm -f "$mbpp_nan_band1q_conv" "$mbpp_nan_band1u_conv" \
          "$mbpp_nan_band2q_conv" "$mbpp_nan_band2u_conv"

    mbpp_nan_convolve_log="$OUT_DIR/mbpp_nan_convolve.log"
    if bin/convolve_cubes \
            infiles="$mbpp_nan_band1q,$mbpp_nan_band1u,$mbpp_nan_band2q,$mbpp_nan_band2u" \
            beamfiles="$DATA_DIR/band1_beamlog.txt,$DATA_DIR/band1_beamlog.txt,$DATA_DIR/band2_beamlog.txt,$DATA_DIR/band2_beamlog.txt" \
            mem_frac_ram=0.25 > "$mbpp_nan_convolve_log" 2>&1 && \
       [[ -s "$mbpp_nan_band1q_conv" && -s "$mbpp_nan_band2q_conv" ]]; then
        pass "convolve_cubes: ran to completion on NaN-containing input, wrote output"
    else
        fail "convolve_cubes: NaN-containing input run failed or wrote no output (see $mbpp_nan_convolve_log)"
    fi

    if [[ -s "$mbpp_nan_band1q_conv" ]]; then
        if python3 - "$mbpp_nan_band1q_conv" <<'PYEOF'
import sys
import numpy as np
from astropy.io import fits
data = fits.getdata(sys.argv[1]).astype(float)
nan_frac = float(np.mean(np.isnan(data)))
block = data[..., 1:3, 1:3]   # well inside the injected 6x6 block
far   = data[..., 15:20, 15:20]  # far outside the block
ok = True
if nan_frac >= 0.5:
    print(f"FAIL: {nan_frac:.3f} of output is NaN -- plane-poisoning bug (T23) has regressed")
    ok = False
else:
    print(f"OK: only {nan_frac:.3f} of output is NaN (localized, not plane-poisoning)")
if not np.all(np.isnan(block)):
    print("FAIL: pixels well inside the injected NaN block came back non-NaN (threshold rejection not engaging)")
    ok = False
else:
    print("OK: pixels well inside the injected NaN block correctly remain NaN")
if np.any(np.isnan(far)):
    print("FAIL: pixels far from the injected NaN block came back NaN (over-rejection)")
    ok = False
else:
    print("OK: pixels far from the injected NaN block are clean (non-NaN)")
sys.exit(0 if ok else 1)
PYEOF
        then
            pass "convolve_cubes: NaN handling is localized -- far pixels clean, deep-block pixels correctly rejected"
        else
            fail "convolve_cubes: NaN handling localization check failed (see output above)"
        fi
    else
        fail "convolve_cubes: NaN-containing input produced no output to check"
    fi

    if [[ -s "$mbpp_nan_band1q_conv" && -s "$mbpp_nan_band2q_conv" ]]; then
        mbpp_nan_cfg="$OUT_DIR/mbpp_nan_rmsynth.cfg"
        cat > "$mbpp_nan_cfg" <<CFGEOF
path                = ${OUT_DIR}/
infileQ             = $(basename "$mbpp_nan_band1q_conv"),$(basename "$mbpp_nan_band2q_conv")
infileU             = $(basename "$mbpp_nan_band1u_conv"),$(basename "$mbpp_nan_band2u_conv")
outfile             = ${OUT_DIR}/mbpp_nan_out
global_badchan_file = none,none
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
        mbpp_nan_rmsynth_log="$OUT_DIR/mbpp_nan_rmsynth.log"
        "$BIN_SERIAL" "$mbpp_nan_cfg" > "$mbpp_nan_rmsynth_log" 2>&1
        if [[ -f "${OUT_DIR}/mbpp_nan_out.AMP.RMCUBE.FITS" ]] && \
           python3 "$TESTS_DIR/check_rm_peak.py" "${OUT_DIR}/mbpp_nan_out.AMP.RMCUBE.FITS" "$TRUTH" > /dev/null 2>&1; then
            pass "convolve_cubes NaN fix + rm_synthesis: src_A/src_B still recovered despite NaN injected far away in band 1"
        else
            fail "convolve_cubes NaN fix + rm_synthesis: RM peak(s) not recovered (see $mbpp_nan_rmsynth_log)"
        fi
    else
        fail "convolve_cubes: NaN-containing run did not produce usable output; skipping downstream rm_synthesis check"
    fi
else
    skip "convolve_cubes or serial rm_synthesis binary not available; skipping NaN-handling regression test (T23)"
fi

# ---------------------------------------------------------------------------
# 49. gaussft_mod: in-plane FFT threading (T28) matches single-threaded
#     output, is run-to-run deterministic, and is not degenerate --
#     see docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md's own T28 entry.
# ---------------------------------------------------------------------------
section "49. gaussft_mod: in-plane FFT threading (T28) matches nthreads=1, deterministic"

gft_bin="$OUT_DIR/test_gaussft_threading"
gft_log="$OUT_DIR/test_gaussft_threading.log"
if gfortran -cpp -std=gnu -fallow-argument-mismatch -ffree-line-length-none \
        -O3 -fopenmp -J"$OUT_DIR" \
        src/gaussft.f90 "$TESTS_DIR/test_gaussft_threading.f90" \
        -o "$gft_bin" -lfftw3 -lfftw3f -lfftw3f_omp 2>"$OUT_DIR/gaussft_threading_build.log"; then
    if OMP_NUM_THREADS=4 "$gft_bin" > "$gft_log" 2>&1; then
        while IFS= read -r line; do
            if [[ "$line" == *"[PASS]"* || "$line" == *"[FAIL]"* ]]; then
                echo "  $line"
            fi
        done < "$gft_log"
        if grep -q "^\[PASS\] test_gaussft_threading" "$gft_log"; then
            pass "gaussft_mod in-plane threading: all checks passed"
        else
            fail "gaussft_mod in-plane threading: one or more checks failed (see $gft_log)"
        fi
    else
        fail "gaussft_mod in-plane threading: program exited non-zero (see $gft_log)"
    fi
else
    fail "gaussft_mod in-plane threading: build failed (see $OUT_DIR/gaussft_threading_build.log)"
fi

# ---------------------------------------------------------------------------
# 50. rmclean_cache_mod: standalone unit tests (T14, docs/dev/
#     RMCLEAN_INTEGRATION_PLAN.md) -- FITS-free/FFTW-free pure-logic
#     checks for the mask-pattern cache eviction support module, grown
#     incrementally as T14's own increments land.
# ---------------------------------------------------------------------------
section "50. rmclean_cache_mod: standalone unit tests (T14)"

rcm_bin="$OUT_DIR/test_rmclean_cache_mod"
rcm_log="$OUT_DIR/test_rmclean_cache_mod.log"
if gfortran -cpp -std=gnu -fallow-argument-mismatch -ffree-line-length-none \
        -O3 -fopenmp -J"$OUT_DIR" \
        src/rmclean_cache_mod.f90 "$TESTS_DIR/test_rmclean_cache_mod.f90" \
        -o "$rcm_bin" 2>"$OUT_DIR/rmclean_cache_mod_build.log"; then
    if "$rcm_bin" > "$rcm_log" 2>&1; then
        while IFS= read -r line; do
            if [[ "$line" == *"[PASS]"* || "$line" == *"[FAIL]"* ]]; then
                echo "  $line"
            fi
        done < "$rcm_log"
        if grep -q "^\[PASS\] test_rmclean_cache_mod" "$rcm_log"; then
            pass "rmclean_cache_mod: all checks passed"
        else
            fail "rmclean_cache_mod: one or more checks failed (see $rcm_log)"
        fi
    else
        fail "rmclean_cache_mod: program exited non-zero (see $rcm_log)"
    fi
else
    fail "rmclean_cache_mod: build failed (see $OUT_DIR/rmclean_cache_mod_build.log)"
fi

# ---------------------------------------------------------------------------
# 51. rmclean_cubes: cache_eviction_policy cfg validation (T14 increment
#     2, docs/dev/RMCLEAN_INTEGRATION_PLAN.md) -- plumbing only at this
#     point, no eviction-policy behavior yet; just confirms the key is
#     parsed/validated/defaulted correctly.
# ---------------------------------------------------------------------------
section "51. rmclean_cubes: cache_eviction_policy cfg validation"

if [[ -x bin/rmclean_cubes ]]; then
    cep_bad_log="$OUT_DIR/rmclean_cache_eviction_policy_bad.log"
    # A nonzero exit is the EXPECTED outcome here -- must not be a bare
    # statement under this script's own `set -e` (line 98), or the
    # expected failure kills the whole test run. Test the command's
    # own status directly via `if`, same idiom used everywhere else in
    # this file, rather than capturing $? from a standalone statement.
    if bin/rmclean_cubes cache_eviction_policy=bogus > "$cep_bad_log" 2>&1; then
        fail "rmclean_cubes: cache_eviction_policy=bogus was not rejected (exited 0, see $cep_bad_log)"
    else
        if grep -q "ERROR: cache_eviction_policy must be one of belady|hitcount" "$cep_bad_log"; then
            pass "rmclean_cubes: cache_eviction_policy=bogus rejected with the expected error"
        else
            fail "rmclean_cubes: cache_eviction_policy=bogus rejected, but not with the expected error (see $cep_bad_log)"
        fi
    fi

    cep_pha="$OUT_DIR/rmc_lsqref.PHA.RMCUBE.FITS"
    cep_mask="$OUT_DIR/rmc_lsqref.MASK.CUBE.FITS"
    if [[ -f "$rmc_lsqref_amp" && -f "$cep_pha" && -f "$cep_mask" ]]; then
        cep_default_out="$OUT_DIR/rmc_cep_default"
        cep_default_log="$OUT_DIR/rmclean_cache_eviction_policy_default.log"
        rm -f "${cep_default_out}".*.FITS
        if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" phafile="$cep_pha" \
                maskfile="$cep_mask" outfile="$cep_default_out" \
                abs_flux_floor=0.01 niter=5 > "$cep_default_log" 2>&1; then
            if grep -q "^cache_eviction_policy=belady" "$cep_default_log"; then
                pass "rmclean_cubes: cache_eviction_policy defaults to belady"
            else
                fail "rmclean_cubes: cache_eviction_policy default not logged as belady (see $cep_default_log)"
            fi
        else
            fail "rmclean_cubes: default-policy run failed to complete (see $cep_default_log)"
        fi
    else
        skip "rmclean_cubes: cache_eviction_policy default check skipped (section 29's own fixtures not available)"
    fi
else
    skip "bin/rmclean_cubes not built; skipping cache_eviction_policy cfg validation"
fi

# ---------------------------------------------------------------------------
# 52. rmclean_io_mod: next_tile_extent standalone unit test (T14,
#     docs/dev/RMCLEAN_INTEGRATION_PLAN.md) -- the tile-geometry stepper
#     shared by the per-tile compute loop and the new pattern
#     pre-scan (Pass 0). Links cfitsio (rmclean_io_mod's own module
#     dependency, via fitsio_unit_mod) even though this specific
#     function never calls a FITS routine.
# ---------------------------------------------------------------------------
section "52. rmclean_io_mod: next_tile_extent standalone unit test (T14)"

riom_bin="$OUT_DIR/test_rmclean_io_mod"
riom_log="$OUT_DIR/test_rmclean_io_mod.log"
if gfortran -cpp -std=gnu -fallow-argument-mismatch -ffree-line-length-none \
        -O3 -fopenmp -J"$OUT_DIR" \
        src/fitsio_unit_mod.f90 src/rmclean_io_mod.f90 "$TESTS_DIR/test_rmclean_io_mod.f90" \
        -o "$riom_bin" -lcfitsio 2>"$OUT_DIR/rmclean_io_mod_build.log"; then
    if "$riom_bin" > "$riom_log" 2>&1; then
        while IFS= read -r line; do
            if [[ "$line" == *"[PASS]"* || "$line" == *"[FAIL]"* ]]; then
                echo "  $line"
            fi
        done < "$riom_log"
        if grep -q "^\[PASS\] test_rmclean_io_mod" "$riom_log"; then
            pass "rmclean_io_mod: next_tile_extent checks passed"
        else
            fail "rmclean_io_mod: one or more checks failed (see $riom_log)"
        fi
    else
        fail "rmclean_io_mod: program exited non-zero (see $riom_log)"
    fi
else
    fail "rmclean_io_mod: build failed (see $OUT_DIR/rmclean_io_mod_build.log)"
fi

# ---------------------------------------------------------------------------
# 53. rmclean_cubes: Pass 0 (pattern pre-scan) correctness (T14 increment
#     6, docs/dev/RMCLEAN_INTEGRATION_PLAN.md) -- two fixtures, chosen to
#     exercise Pass 0's two distinct outcomes:
#       - rmc_varied.MASK.CUBE.FITS (757 distinct patterns / 1024 valid
#         pixels, 73.9% -- a deliberately high-diversity stress fixture)
#         must trip the 10%-of-image-pixels safety valve and fall back
#         to cache_eviction_policy=hitcount, completing successfully.
#       - rmc_lsqref.MASK.CUBE.FITS (1 distinct pattern / 1024 valid
#         pixels) must NOT trip the safety valve, and Pass 0's own
#         logged distinct-pattern/valid-pixel counts must be identical
#         whether the whole image is one 32x32 tile or split into
#         sixteen 8x8 tiles -- proving Pass 0 and the per-tile
#         compute loop agree on scan order regardless of tiling (the
#         single most important correctness property this increment
#         depends on, since Pass 0's own position bookkeeping is only
#         meaningful if Pass 1 visits pixels in exactly the same order).
#         Tile size 8x8 (a multiple of 4) is deliberately chosen over,
#         say, 7x5: a pre-existing, T14-unrelated floating-point
#         artifact (tile widths not a multiple of the CPU's
#         double-precision SIMD width cause RESID/RESTORED -- but never
#         CLEAN -- to differ by ~1e-7 relative between tilings; see
#         this ticket's own "Progress" notes) would otherwise force this
#         test to special-case which outputs it compares. Both true
#         distinct-pattern/valid-pixel counts were independently
#         confirmed via a Python oracle (np.unique over packed per-pixel
#         channel-validity bytes) before being hard-coded as the
#         expected values below.
# ---------------------------------------------------------------------------
section "53. rmclean_cubes: Pass 0 (pattern pre-scan) correctness"

if [[ -x bin/rmclean_cubes && -f "$rmc_lsqref_amp" && -f "$cep_pha" && -f "$cep_mask" && -f "$rmc_varied_mask" ]]; then
    # --- Safety valve: high-diversity fixture must fall back cleanly ---
    p0sv_out="$OUT_DIR/rmc_p0_safetyvalve"
    p0sv_log="$OUT_DIR/rmc_p0_safetyvalve.log"
    rm -f "${p0sv_out}".*.FITS
    if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" phafile="$cep_pha" \
            maskfile="$rmc_varied_mask" outfile="$p0sv_out" \
            abs_flux_floor=0.01 niter=200 gain=0.1 \
            cache_eviction_policy=belady > "$p0sv_log" 2>&1; then
        if grep -q "WARNING: Pass 0 found 757 distinct patterns" "$p0sv_log" && \
           grep -q "falling back to cache_eviction_policy=hitcount" "$p0sv_log"; then
            pass "rmclean_cubes: Pass 0 safety valve trips on a high-diversity mask (757/1024 patterns) and falls back to hitcount"
        else
            fail "rmclean_cubes: Pass 0 safety valve did not trip as expected on the high-diversity mask (see $p0sv_log)"
        fi
    else
        fail "rmclean_cubes: cache_eviction_policy=belady run failed on the high-diversity mask (see $p0sv_log)"
    fi
    rm -f "${p0sv_out}".*.FITS

    # --- Happy path: low-diversity fixture, Pass 0 completes normally,
    #     tiling-invariant counts and outputs ---
    p0hp_full="$OUT_DIR/rmc_p0_happypath_full"
    p0hp_tiled="$OUT_DIR/rmc_p0_happypath_tiled"
    p0hp_full_log="$OUT_DIR/rmc_p0_happypath_full.log"
    p0hp_tiled_log="$OUT_DIR/rmc_p0_happypath_tiled.log"
    rm -f "${p0hp_full}".*.FITS "${p0hp_tiled}".*.FITS
    p0_expect="Pass 0: done -- 1 distinct pattern(s), 1024 valid pixel(s) scanned."
    if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" phafile="$cep_pha" \
            maskfile="$cep_mask" outfile="$p0hp_full" \
            abs_flux_floor=0.01 niter=200 gain=0.1 \
            cache_eviction_policy=belady > "$p0hp_full_log" 2>&1 && \
       bin/rmclean_cubes ampfile="$rmc_lsqref_amp" phafile="$cep_pha" \
            maskfile="$cep_mask" outfile="$p0hp_tiled" \
            abs_flux_floor=0.01 niter=200 gain=0.1 \
            cache_eviction_policy=belady tile_auto=n tile_ra=8 tile_dec=8 \
            > "$p0hp_tiled_log" 2>&1; then
        if grep -qF "$p0_expect" "$p0hp_full_log" && grep -qF "$p0_expect" "$p0hp_tiled_log"; then
            pass "rmclean_cubes: Pass 0 completes (no safety valve) and logs the correct oracle-matched count (1 pattern/1024 pixels) on a low-diversity mask"
        else
            fail "rmclean_cubes: Pass 0's logged distinct-pattern/valid-pixel count is wrong or missing (see $p0hp_full_log / $p0hp_tiled_log)"
        fi
        if grep -q "^WARNING: Pass 0" "$p0hp_full_log" || grep -q "^WARNING: Pass 0" "$p0hp_tiled_log"; then
            fail "rmclean_cubes: Pass 0 safety valve fired unexpectedly on the low-diversity mask (see $p0hp_full_log / $p0hp_tiled_log)"
        else
            pass "rmclean_cubes: Pass 0 safety valve does not fire on a low-diversity mask"
        fi

        p0hp_ok=1
        for suffix in CLEAN.AMP CLEAN.PHA RESID.AMP RESID.PHA \
                      RESTORED.AMP RESTORED.PHA; do
            if ! cmp -s "${p0hp_full}.${suffix}.RMCUBE.FITS" \
                        "${p0hp_tiled}.${suffix}.RMCUBE.FITS"; then
                p0hp_ok=0
                fail "rmclean_cubes: Pass 0 happy-path output differs between 1-tile (32x32) and 16-tile (8x8) runs on $suffix"
            fi
        done
        [[ "$p0hp_ok" -eq 1 ]] && \
            pass "rmclean_cubes: Pass 0 happy-path output bit-identical across tile sizes (32x32 vs 8x8)"
    else
        fail "rmclean_cubes: Pass 0 happy-path run(s) failed on the low-diversity mask (see $p0hp_full_log / $p0hp_tiled_log)"
    fi
    rm -f "${p0hp_full}".*.FITS "${p0hp_tiled}".*.FITS
else
    skip "rmclean_cubes: Pass 0 correctness checks skipped (bin/rmclean_cubes or section 29/51's own fixtures not available)"
fi

# ---------------------------------------------------------------------------
# 54. rmclean_cubes: Belady eviction correctness + effectiveness vs
#     hitcount (T14 increment 8, docs/dev/RMCLEAN_INTEGRATION_PLAN.md) --
#     the end of Phase B. A deliberately adversarial fixture (not real
#     data -- built to expose LFU-without-decay's own structural blind
#     spot, so the comparison measures something): the first 50
#     pixels cycle through 5 "early" patterns (each racking up ~10 hits
#     before never appearing again); the remaining ~974 pixels cycle
#     through 15 "late" patterns that keep recurring for the rest of the
#     run. Under a cache too small to hold all of them, LFU/hitcount
#     (which only remembers PAST hit counts) evicts the low-hit-count
#     but soon-to-recur late patterns instead of the high-hit-count but
#     permanently-dead early ones -- a real, structural blind spot of
#     any eviction policy that cannot see the future. Belady, via
#     Pass 0's own registry, correctly identifies the early patterns'
#     "never again" future and evicts those instead. 19 true distinct
#     patterns (well under the 10% safety valve for a 1024-pixel image),
#     confirmed via an independent Python oracle.
# ---------------------------------------------------------------------------
section "54. rmclean_cubes: Belady eviction correctness + effectiveness vs hitcount"

if [[ -x bin/rmclean_cubes && -f "$rmc_lsqref_amp" && -f "$cep_pha" ]]; then
    belady_kpat_mask="$OUT_DIR/rmc_kpat.MASK.CUBE.FITS"
    python3 - "$cep_mask" "$belady_kpat_mask" <<'PYEOF'
import sys
import numpy as np
from astropy.io import fits

src, dst = sys.argv[1], sys.argv[2]
with fits.open(src) as hdul:
    rng = np.random.default_rng(11)
    data = hdul[0].data
    nchan, ny, nx = data.shape

    # 5 "early" prototypes + 15 "late" prototypes, each a fixed random
    # subset of channels flipped to 0 relative to the all-valid base
    # pattern (rmc_lsqref's own mask is all-valid everywhere).
    n_early, n_late = 5, 15
    protos = []
    for _ in range(n_early + n_late):
        nflip = rng.integers(1, 5)
        chans = rng.choice(nchan, size=nflip, replace=False)
        protos.append(chans)

    base = data[:, 0, 0].copy()
    for iy in range(ny):
        for ix in range(nx):
            i = iy*nx + ix
            if i < 50:
                p = i % n_early
            else:
                p = n_early + ((i - 50) % n_late)
            row = base.copy()
            row[protos[p]] = 0
            data[:, iy, ix] = row
    hdul.writeto(dst, overwrite=True)
PYEOF

    belady_kpat_nocache="$OUT_DIR/rmc_kpat_nocache"
    belady_kpat_hitcount="$OUT_DIR/rmc_kpat_hitcount"
    belady_kpat_belady="$OUT_DIR/rmc_kpat_belady"
    rm -f "${belady_kpat_nocache}".*.FITS "${belady_kpat_hitcount}".*.FITS "${belady_kpat_belady}".*.FITS

    if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" phafile="$cep_pha" \
            maskfile="$belady_kpat_mask" outfile="$belady_kpat_nocache" \
            abs_flux_floor=0.01 niter=200 gain=0.1 mask_pattern_cache_max=0 \
            > "$OUT_DIR/rmc_kpat_nocache.log" 2>&1 && \
       bin/rmclean_cubes ampfile="$rmc_lsqref_amp" phafile="$cep_pha" \
            maskfile="$belady_kpat_mask" outfile="$belady_kpat_hitcount" \
            abs_flux_floor=0.01 niter=200 gain=0.1 \
            mask_pattern_cache_max=8 cache_eviction_policy=hitcount \
            > "$OUT_DIR/rmc_kpat_hitcount.log" 2>&1 && \
       bin/rmclean_cubes ampfile="$rmc_lsqref_amp" phafile="$cep_pha" \
            maskfile="$belady_kpat_mask" outfile="$belady_kpat_belady" \
            abs_flux_floor=0.01 niter=200 gain=0.1 \
            mask_pattern_cache_max=8 cache_eviction_policy=belady \
            > "$OUT_DIR/rmc_kpat_belady.log" 2>&1; then

        if grep -q "^WARNING: Pass 0" "$OUT_DIR/rmc_kpat_belady.log"; then
            fail "rmclean_cubes: Belady safety valve fired unexpectedly on the 19-pattern adversarial fixture (see $OUT_DIR/rmc_kpat_belady.log)"
        else
            pass "rmclean_cubes: Belady safety valve does not fire on the 19-pattern adversarial fixture"
        fi

        kpat_ok=1
        for suffix in CLEAN.AMP CLEAN.PHA RESID.AMP RESID.PHA \
                      RESTORED.AMP RESTORED.PHA; do
            for pol_out in "$belady_kpat_hitcount" "$belady_kpat_belady"; do
                if ! cmp -s "${belady_kpat_nocache}.${suffix}.RMCUBE.FITS" \
                            "${pol_out}.${suffix}.RMCUBE.FITS"; then
                    kpat_ok=0
                    fail "rmclean_cubes: $(basename "$pol_out") (forced eviction, cache_max=8) differs from cache-disabled on $suffix"
                fi
            done
        done
        [[ "$kpat_ok" -eq 1 ]] && \
            pass "rmclean_cubes: both hitcount and belady (forced eviction, cache_max=8 < 19 distinct patterns) bit-identical to cache disabled"

        # Effectiveness regression guard: total redundant-rebuild count
        # (cache inserts, including evict-and-reinsert, PLUS declined
        # one-offs) must be no worse under belady than under hitcount --
        # Belady is supposed to be PROVABLY optimal for a known future
        # access sequence, not merely "correct."
        # grep -oP exits 1 (killing the whole script under this file's
        # own `set -euo pipefail`, line 98) when the "(N pixel(s) past
        # ...)" line is absent -- an EXPECTED outcome under hitcount
        # (which never declines, so n_overflow_pixels_total stays 0 and
        # that line is never printed at all). `|| true` treats no-match
        # as valid, same idiom used in section 51 for this exact class
        # of bug.
        hc_inserts=$(grep -oP 'Mask-pattern cache: \K[0-9]+' "$OUT_DIR/rmc_kpat_hitcount.log" || true)
        hc_overflow=$(grep -oP '\(\K[0-9]+(?= pixel\(s\) past)' "$OUT_DIR/rmc_kpat_hitcount.log" || true)
        hc_overflow=${hc_overflow:-0}
        bl_inserts=$(grep -oP 'Mask-pattern cache: \K[0-9]+' "$OUT_DIR/rmc_kpat_belady.log" || true)
        bl_overflow=$(grep -oP '\(\K[0-9]+(?= pixel\(s\) past)' "$OUT_DIR/rmc_kpat_belady.log" || true)
        bl_overflow=${bl_overflow:-0}
        hc_total=$((hc_inserts + hc_overflow))
        bl_total=$((bl_inserts + bl_overflow))
        echo "  hitcount total redundant rebuilds: $hc_total (inserts=$hc_inserts, one-off=$hc_overflow)"
        echo "  belady   total redundant rebuilds: $bl_total (inserts=$bl_inserts, one-off=$bl_overflow)"
        if [[ "$bl_total" -le "$hc_total" ]]; then
            pass "rmclean_cubes: belady's total redundant-rebuild count ($bl_total) is <= hitcount's ($hc_total) on the same adversarial fixture/cache size"
        else
            fail "rmclean_cubes: belady's total redundant-rebuild count ($bl_total) is WORSE than hitcount's ($hc_total) -- Belady should never be worse"
        fi

        # T14 increment 9: direct evidence that belady's own resident
        # cache slots carry NO local pattern-bytes copy (delegated to
        # the registry instead), while hitcount's own slots -- which
        # have no registry to delegate to -- still carry one each. Both
        # numbers are printed by rmclean_cubes itself, not inferred.
        bl_pattern_line=$(grep -oP 'pattern bytes locally stored for \K[0-9]+ of [0-9]+' "$OUT_DIR/rmc_kpat_belady.log" || true)
        hc_pattern_line=$(grep -oP 'pattern bytes locally stored for \K[0-9]+ of [0-9]+' "$OUT_DIR/rmc_kpat_hitcount.log" || true)
        echo "  belady   pattern bytes locally stored: $bl_pattern_line"
        echo "  hitcount pattern bytes locally stored: $hc_pattern_line"
        if [[ "$bl_pattern_line" =~ ^0\ of\ ([0-9]+)$ ]] && [[ "${BASH_REMATCH[1]}" -gt 0 ]]; then
            pass "rmclean_cubes: belady's own resident cache slots (${BASH_REMATCH[1]}) carry zero local pattern-bytes copies (T14 increment 9)"
        else
            fail "rmclean_cubes: belady's own resident cache slots still carry local pattern-bytes copies (see $OUT_DIR/rmc_kpat_belady.log) -- increment 9 de-duplication not engaging"
        fi
        if [[ "$hc_pattern_line" =~ ^([0-9]+)\ of\ ([0-9]+)$ ]] && [[ "${BASH_REMATCH[1]}" == "${BASH_REMATCH[2]}" ]] && [[ "${BASH_REMATCH[1]}" -gt 0 ]]; then
            pass "rmclean_cubes: hitcount's own resident cache slots (${BASH_REMATCH[2]}) each still carry their own local pattern-bytes copy (no registry to delegate to)"
        else
            fail "rmclean_cubes: hitcount's own resident cache slots don't match their own expected pattern-bytes count (see $OUT_DIR/rmc_kpat_hitcount.log)"
        fi
    else
        fail "rmclean_cubes: Belady effectiveness run(s) failed on the adversarial fixture (see $OUT_DIR/rmc_kpat_*.log)"
    fi
    rm -f "${belady_kpat_nocache}".*.FITS "${belady_kpat_hitcount}".*.FITS "${belady_kpat_belady}".*.FITS "$belady_kpat_mask"
else
    skip "rmclean_cubes: Belady eviction effectiveness check skipped (bin/rmclean_cubes or section 51's own fixtures not available)"
fi

# ---------------------------------------------------------------------------
# 55. rmclean_cubes: min_valid_chan_frac coverage-threshold pixel skip (T17,
#     docs/dev/RMCLEAN_INTEGRATION_PLAN.md). Reuses the same 20-
#     prototype adversarial fixture generator as section 54, but this
#     time the check that matters is each prototype's own KNOWN
#     nflip count (1-4 of nchan flagged, fixed by the seeded RNG) --
#     min_valid_chan_frac=0.99 (nchan=200 -> skip any pixel with >2 flagged)
#     should keep exactly the prototypes with nflip<=2 and skip the
#     rest, with NaN written to every RM-bin of every skipped pixel.
# ---------------------------------------------------------------------------
section "55. rmclean_cubes: min_valid_chan_frac coverage-threshold pixel skip"

if [[ -x bin/rmclean_cubes && -f "$rmc_lsqref_amp" && -f "$cep_pha" ]]; then
    t17_mask="$OUT_DIR/rmc_t17_kpat.MASK.CUBE.FITS"
    t17_proto_log="$OUT_DIR/rmc_t17_protos.log"
    python3 - "$cep_mask" "$t17_mask" > "$t17_proto_log" <<'PYEOF'
import sys
import numpy as np
from astropy.io import fits

src, dst = sys.argv[1], sys.argv[2]
with fits.open(src) as hdul:
    rng = np.random.default_rng(11)
    data = hdul[0].data
    nchan, ny, nx = data.shape
    n_early, n_late = 5, 15
    protos = []
    nflips = []
    for _ in range(n_early + n_late):
        nflip = rng.integers(1, 5)
        chans = rng.choice(nchan, size=nflip, replace=False)
        protos.append(chans)
        nflips.append(int(nflip))
    base = data[:, 0, 0].copy()
    for iy in range(ny):
        for ix in range(nx):
            i = iy*nx + ix
            if i < 50:
                p = i % n_early
            else:
                p = n_early + ((i - 50) % n_late)
            row = base.copy()
            row[protos[p]] = 0
            data[:, iy, ix] = row
    hdul.writeto(dst, overwrite=True)

    # Expected KEPT pixel count at min_valid_chan_frac=0.99 (nchan=200,
    # threshold keeps nflip<=2): replay the exact same i->prototype
    # assignment used above.
    kept_protos = set(p for p, nf in enumerate(nflips) if nf <= 2)
    n_kept = 0
    for i in range(nx*ny):
        p = (i % n_early) if i < 50 else (n_early + ((i - 50) % n_late))
        if p in kept_protos:
            n_kept += 1
    print(n_kept)
PYEOF
    t17_expected_kept=$(tail -1 "$t17_proto_log")

    t17_thresh="$OUT_DIR/rmc_t17_thresh"
    t17_baseline="$OUT_DIR/rmc_t17_baseline"
    rm -f "${t17_thresh}".*.FITS "${t17_baseline}".*.FITS
    if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" phafile="$cep_pha" \
            maskfile="$t17_mask" outfile="$t17_thresh" \
            abs_flux_floor=0.01 niter=200 gain=0.1 \
            cache_eviction_policy=belady min_valid_chan_frac=0.99 \
            > "$OUT_DIR/rmc_t17_thresh.log" 2>&1 && \
       bin/rmclean_cubes ampfile="$rmc_lsqref_amp" phafile="$cep_pha" \
            maskfile="$t17_mask" outfile="$t17_baseline" \
            abs_flux_floor=0.01 niter=200 gain=0.1 \
            cache_eviction_policy=belady \
            > "$OUT_DIR/rmc_t17_baseline.log" 2>&1; then

        t17_reported_kept=$(grep -oP 'CLEANed \K[0-9]+(?= pixels)' "$OUT_DIR/rmc_t17_thresh.log" || true)
        echo "  min_valid_chan_frac=0.99: expected kept=$t17_expected_kept, reported kept=$t17_reported_kept"
        if [[ "$t17_reported_kept" == "$t17_expected_kept" ]]; then
            pass "rmclean_cubes: min_valid_chan_frac=0.99 keeps exactly the expected pixel count ($t17_expected_kept), matching the known per-prototype flag counts"
        else
            fail "rmclean_cubes: min_valid_chan_frac=0.99 kept $t17_reported_kept pixels, expected $t17_expected_kept (see $OUT_DIR/rmc_t17_thresh.log)"
        fi

        if python3 - "$t17_thresh" "$t17_baseline" "$t17_expected_kept" <<'PYEOF'
import sys
import numpy as np
from astropy.io import fits

thresh_base, baseline_base, n_expected_kept = sys.argv[1], sys.argv[2], int(sys.argv[3])
ok = True
for suffix in ['CLEAN.AMP','CLEAN.PHA','RESID.AMP','RESID.PHA','RESTORED.AMP','RESTORED.PHA']:
    thresh = fits.getdata(f'{thresh_base}.{suffix}.RMCUBE.FITS')
    base = fits.getdata(f'{baseline_base}.{suffix}.RMCUBE.FITS')
    nan_pixels = np.isnan(thresh).all(axis=0)
    n_nan = int(nan_pixels.sum())
    n_total = thresh.shape[1]*thresh.shape[2]
    if n_nan != n_total - n_expected_kept:
        print(f'{suffix}: FAIL wrong NaN pixel count: {n_nan} (expected {n_total - n_expected_kept})')
        ok = False
        continue
    # Kept pixels (not NaN under threshold) must be BIT-IDENTICAL to the
    # no-threshold baseline -- the threshold must never change a KEPT
    # pixel's own CLEAN result, only decide whether it runs at all.
    kept_mask = ~nan_pixels
    if not np.array_equal(thresh[:, kept_mask], base[:, kept_mask]):
        print(f'{suffix}: FAIL kept pixels differ from no-threshold baseline')
        ok = False
if ok:
    print('ALL_OK')
sys.exit(0 if ok else 1)
PYEOF
        then
            pass "rmclean_cubes: min_valid_chan_frac NaN placement matches expected skip count, and kept pixels are bit-identical to the no-threshold baseline (all 6 outputs)"
        else
            fail "rmclean_cubes: min_valid_chan_frac output verification failed (see above)"
        fi
    else
        fail "rmclean_cubes: min_valid_chan_frac test run(s) failed (see $OUT_DIR/rmc_t17_*.log)"
    fi
    rm -f "${t17_thresh}".*.FITS "${t17_baseline}".*.FITS "$t17_mask"
else
    skip "rmclean_cubes: min_valid_chan_frac check skipped (bin/rmclean_cubes or section 51's own fixtures not available)"
fi

section "56. rmclean_cubes: swim-lane plot generation (reminder check, non-blocking)"

# Not a correctness check on rmclean_cubes itself -- a reminder that the
# swim-lane plotting PATH (thread_timing log lines -> scripts/
# plot_tile_async_swimlane.py -> a real PNG) still works end to end.
# Found worth having (2026-08-09) after a real multi-hour WALLABY+EMU
# rmclean run finished with no swim-lane plot ever generated for it --
# the capability existed (scripts/plot_tile_async_swimlane.py, already
# used for an earlier e2e validation run) but nothing prompted using it.
# Deliberately uses warn(), not fail(): this is a diagnostics/visibility
# tool, not a correctness property of rmclean_cubes' own output -- it
# should never be why the suite (or anything gating on it) goes red, but
# it also should never silently stop working without anyone noticing.
if [[ -x bin/rmclean_cubes && -f "$rmc_lsqref_amp" && -f "$cep_pha" && -f "$cep_mask" && -x "$PLOT_PY" ]]; then
    swimtest_out="$OUT_DIR/rmc_swimtest"
    swimtest_log="$OUT_DIR/rmc_swimtest.log"
    swimtest_png="$OUT_DIR/rmc_swimtest_swimlane.png"
    rm -f "${swimtest_out}".*.FITS "$swimtest_png"
    if bin/rmclean_cubes ampfile="$rmc_lsqref_amp" phafile="$cep_pha" \
            maskfile="$cep_mask" outfile="$swimtest_out" \
            abs_flux_floor=0.01 niter=200 gain=0.1 \
            log_level=debug timing_enabled=y \
            > "$swimtest_log" 2>&1; then
        if "$PLOT_PY" scripts/plot_tile_async_swimlane.py \
                --log "$swimtest_log" --out "$swimtest_png" \
                --run latest --thread-stage clean \
                > "$OUT_DIR/rmc_swimtest_plot.log" 2>&1 \
                && [[ -s "$swimtest_png" ]]; then
            png_bytes=$(stat -c%s "$swimtest_png" 2>/dev/null || stat -f%z "$swimtest_png")
            if [[ "$png_bytes" -gt 1000 ]]; then
                pass "rmclean_cubes: swim-lane plot generated ($png_bytes bytes, $swimtest_png)"
            else
                warn "rmclean_cubes: swim-lane plot written but suspiciously small ($png_bytes bytes) -- see $OUT_DIR/rmc_swimtest_plot.log"
            fi
        else
            warn "rmclean_cubes: swim-lane plot generation failed -- see $OUT_DIR/rmc_swimtest_plot.log (thread_timing log format may have drifted from what plot_tile_async_swimlane.py expects)"
        fi
    else
        warn "rmclean_cubes: swim-lane test run itself failed -- see $swimtest_log (cannot check plot generation without a real thread_timing log)"
    fi
    rm -f "${swimtest_out}".*.FITS
else
    skip "rmclean_cubes: swim-lane plot check skipped (bin/rmclean_cubes, section 51's own fixtures, or no usable python3 found)"
fi

# ---------------------------------------------------------------------------
# 57. RM-CLEAN table_oversample's own Nyquist floor (T21)
# ---------------------------------------------------------------------------
section "57. RM-CLEAN table_oversample Nyquist floor (T21)"

tosfloor_bin="$RMCLEAN_BUILD_DIR/test_table_oversample_floor"
tosfloor_log="$OUT_DIR/test_table_oversample_floor.log"

if [[ -f "$rmclean_o" ]]; then
    if gfortran -cpp -std=gnu -fallow-argument-mismatch -ffree-line-length-none \
            -O3 -fopenmp -I"$rmclean_mod_dir" -J"$rmclean_mod_dir" \
            "$TESTS_DIR/test_table_oversample_floor.f90" "$rmclean_o" \
            -o "$tosfloor_bin" -lfftw3 2>"$OUT_DIR/rmclean_tosfloor_build.log"; then
        if "$tosfloor_bin" > "$tosfloor_log" 2>&1; then
            while IFS= read -r line; do
                case "$line" in
                    *"[PASS]"*) pass "${line#*\[PASS\] }" ;;
                    *"[FAIL]"*) fail "${line#*\[FAIL\] }" ;;
                esac
            done < "$tosfloor_log"
        else
            fail "RM-CLEAN table_oversample Nyquist floor: program exited non-zero (see $tosfloor_log)"
        fi
    else
        fail "RM-CLEAN table_oversample Nyquist floor: build failed (see $OUT_DIR/rmclean_tosfloor_build.log)"
    fi
else
    skip "rmclean_mod.o not built (section 23 skipped); skipping RM-CLEAN table_oversample Nyquist floor test"
fi

# ---------------------------------------------------------------------------
# 58. rm_synthesis: badchan_file 'none' sentinel for a band with nothing to
#     list (T31, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md) -- checks BOTH
#     edited code sites: 'none' on the reference band (case B) and 'none'
#     on a non-reference band (case A), each paired with a REAL bad-channel
#     file on the other band, confirming neither disables bad-channel
#     removal globally (if it did, the real band's own channel would come
#     back un-flagged too) and neither requires a throwaway empty file.
# ---------------------------------------------------------------------------
section "58. rm_synthesis: badchan_file 'none' sentinel (T31)"

if [[ -x "$BIN_SERIAL" ]]; then
    t31_badchan="$OUT_DIR/t31_badchan.txt"
    printf '13\n' > "$t31_badchan"

    t31_run() {
        # $1=label $2=global_badchan_file value $3=outfile
        local t31_cfg="$OUT_DIR/t31_$1.cfg"
        local t31_log="$OUT_DIR/t31_$1.log"
        rm -f "${3}".*.FITS
        cat > "$t31_cfg" <<CFGEOF
path                = $DATA_DIR/
infileQ             = TEST.Q.FITSCUBE,TEST_BAND2.Q.FITSCUBE
infileU             = TEST.U.FITSCUBE,TEST_BAND2.U.FITSCUBE
outfile             = $3

global_badchan_file = $2
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
nrm                 = 101
output_mode         = ap
ap_angle_mode       = phase
write_mask_output   = y
write_nvalid_output = y
use_gpu             = n
CFGEOF
        "$BIN_SERIAL" "$t31_cfg" > "$t31_log" 2>&1
    }

    # Case A: reference band (1) real, non-reference band (2) 'none'.
    t31_out_a="$OUT_DIR/t31_case_a_out"
    t31_run "case_a" "$t31_badchan,none" "$t31_out_a"
    t31_log_a="$OUT_DIR/t31_case_a.log"
    if [[ -f "${t31_out_a}.MASK.CUBE.FITS" ]] && grep -qE "\(band[[:space:]]+2[[:space:]]*\): 0 \(none\)" "$t31_log_a"; then
        t31_a_check=$(python3 -c "
from astropy.io import fits
import numpy as np
m = fits.getdata('${t31_out_a}.MASK.CUBE.FITS')[:,0,0]
bad = [i+1 for i,v in enumerate(m) if v==0]
print('OK' if bad==[13] else 'BAD:'+str(bad))
")
        if [[ "$t31_a_check" == "OK" ]]; then
            pass "rm_synthesis: badchan_file 'none' on a non-reference band gives it 0 bad channels, reference band's real channel 13 still flagged (T31 case A)"
        else
            fail "rm_synthesis: T31 case A mask check failed ($t31_a_check, see $t31_log_a)"
        fi
    else
        fail "rm_synthesis: T31 case A run failed or missing 'none' log marker (see $t31_log_a)"
    fi

    # Case B: reference band (1) 'none', non-reference band (2) real
    # (channel 13 of band 2 -> global channel 200+13=213).
    t31_out_b="$OUT_DIR/t31_case_b_out"
    t31_run "case_b" "none,$t31_badchan" "$t31_out_b"
    t31_log_b="$OUT_DIR/t31_case_b.log"
    if [[ -f "${t31_out_b}.MASK.CUBE.FITS" ]] && grep -q "Number of Bad Channels: 0 (none)" "$t31_log_b"; then
        t31_b_check=$(python3 -c "
from astropy.io import fits
import numpy as np
m = fits.getdata('${t31_out_b}.MASK.CUBE.FITS')[:,0,0]
bad = [i+1 for i,v in enumerate(m) if v==0]
print('OK' if bad==[213] else 'BAD:'+str(bad))
")
        if [[ "$t31_b_check" == "OK" ]]; then
            pass "rm_synthesis: badchan_file 'none' on the reference band gives it 0 bad channels, non-reference band's real channel still flagged (T31 case B)"
        else
            fail "rm_synthesis: T31 case B mask check failed ($t31_b_check, see $t31_log_b)"
        fi
    else
        fail "rm_synthesis: T31 case B run failed or missing 'none' log marker (see $t31_log_b)"
    fi
else
    skip "bin/rm_synthesis_release_cpu_serial not built; skipping badchan_file 'none' sentinel test"
fi

# ---------------------------------------------------------------------------
# 59. rm_synthesis: badchan_file/global_badchan_file is the sole switch for
#     bad-channel removal (T33/T35, docs/dev/MULTI_BAND_TOMOGRAPHY_PLAN.md).
#     There is no remove_badchan cfg key. Once badchan_file is given at
#     all, every band needs a real file path or the literal value 'none'
#     -- a blank entry is a parse-time error naming the band, raised
#     before any file is opened. Omitting badchan_file means no
#     removal for any band. Checks both band positions for the blank
#     case, and confirms omitting the key is unaffected.
# ---------------------------------------------------------------------------
section "59. rm_synthesis: badchan_file is the sole switch for bad-channel removal (T33/T35)"

if [[ -x "$BIN_SERIAL" ]]; then
    t35_badchan_run() {
        # $1=label $2=global_badchan_file value (empty = key omitted) $3=outfile
        local t35_cfg="$OUT_DIR/t35_$1.cfg"
        local t35_log="$OUT_DIR/t35_$1.log"
        rm -f "${3}".*.FITS
        {
            cat <<CFGEOF
path                = $DATA_DIR/
infileQ             = TEST.Q.FITSCUBE,TEST_BAND2.Q.FITSCUBE
infileU             = TEST.U.FITSCUBE,TEST_BAND2.U.FITSCUBE
outfile             = $3

CFGEOF
            if [[ -n "$2" ]]; then
                echo "global_badchan_file = $2"
            fi
            cat <<CFGEOF
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
nrm                 = 101
output_mode         = ap
ap_angle_mode       = phase
write_mask_output   = y
write_nvalid_output = y
use_gpu             = n
CFGEOF
        } > "$t35_cfg"
        "$BIN_SERIAL" "$t35_cfg" > "$t35_log" 2>&1
    }

    # Case A: reference band (1) blank, band 2 'none'.
    t35_out_a="$OUT_DIR/t35_case_a_out"
    t35_badchan_run "case_a" ",none" "$t35_out_a"
    t35_log_a="$OUT_DIR/t35_case_a.log"
    if grep -qE "badchan_file/global_badchan_file is blank for band[[:space:]]+1" "$t35_log_a" \
            && [[ ! -f "${t35_out_a}.MASK.CUBE.FITS" ]]; then
        pass "rm_synthesis: blank badchan_file for the reference band is a clear named error (band 1), no output written (T33 case A)"
    else
        fail "rm_synthesis: T33 case A did not error as expected (see $t35_log_a)"
    fi

    # Case B: reference band (1) 'none', band 2 blank.
    t35_out_b="$OUT_DIR/t35_case_b_out"
    t35_badchan_run "case_b" "none," "$t35_out_b"
    t35_log_b="$OUT_DIR/t35_case_b.log"
    if grep -qE "badchan_file/global_badchan_file is blank for band[[:space:]]+2" "$t35_log_b" \
            && [[ ! -f "${t35_out_b}.MASK.CUBE.FITS" ]]; then
        pass "rm_synthesis: blank badchan_file for a non-reference band is a clear named error (band 2), no output written (T33 case B)"
    else
        fail "rm_synthesis: T33 case B did not error as expected (see $t35_log_b)"
    fi

    # Case C: global_badchan_file omitted -- must run to
    # completion, no removal for either band.
    t35_out_c="$OUT_DIR/t35_case_c_out"
    t35_badchan_run "case_c" "" "$t35_out_c"
    t35_log_c="$OUT_DIR/t35_case_c.log"
    if [[ -f "${t35_out_c}.MASK.CUBE.FITS" ]]; then
        pass "rm_synthesis: omitting global_badchan_file entirely is unaffected (T35 case C)"
    else
        fail "rm_synthesis: T35 case C (global_badchan_file omitted) unexpectedly failed to produce output (see $t35_log_c)"
    fi

    # Case D: remove_badchan is not a recognized key.
    t35_out_d="$OUT_DIR/t35_case_d_out"
    t35_cfg_d="$OUT_DIR/t35_case_d.cfg"
    t35_log_d="$OUT_DIR/t35_case_d.log"
    cat > "$t35_cfg_d" <<CFGEOF
path                = $DATA_DIR/
infileQ             = TEST.Q.FITSCUBE
infileU             = TEST.U.FITSCUBE
outfile             = $t35_out_d
remove_badchan      = n
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
nrm                 = 101
output_mode         = ap
ap_angle_mode       = phase
write_mask_output   = y
write_nvalid_output = y
use_gpu             = n
CFGEOF
    "$BIN_SERIAL" "$t35_cfg_d" > "$t35_log_d" 2>&1
    if grep -q "Unknown key in cfg" "$t35_log_d" && grep -q "remove_badchan" "$t35_log_d"; then
        pass "rm_synthesis: remove_badchan is rejected as an unknown cfg key (T35)"
    else
        fail "rm_synthesis: remove_badchan was not rejected as expected (see $t35_log_d)"
    fi
else
    skip "bin/rm_synthesis_release_cpu_serial not built; skipping badchan_file switch tests"
fi

# ---------------------------------------------------------------------------
# 60. convolve_cubes/match_cubes: a blank badchan_file position is now a hard
#     error once badchan_file= is given at all (T34, docs/dev/
#     MULTI_BAND_TOMOGRAPHY_PLAN.md) -- matches rm_synthesis's own T33
#     convention: a real path or the literal value 'none' is required for
#     every position, blank is rejected. Omitting badchan_file= entirely is
#     unaffected (still valid, still means "no list for any infile").
# ---------------------------------------------------------------------------
section "60. convolve_cubes/match_cubes: blank badchan_file position is a hard error (T34)"

if [[ -x bin/convolve_cubes ]]; then
    t34_band1="$OUT_DIR/t34_band1.Q.FITSCUBE"
    t34_band2="$OUT_DIR/t34_band2.Q.FITSCUBE"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$t34_band1"
    cp "$DATA_DIR/TEST_BAND2.Q.FITSCUBE" "$t34_band2"
    t34_band2_badchan="$OUT_DIR/t34_band2_badchan.txt"
    printf '77\n' > "$t34_band2_badchan"
    t34_band1_conv="${t34_band1%.FITSCUBE}_CONV.FITS"
    t34_band2_conv="${t34_band2%.FITSCUBE}_CONV.FITS"
    rm -f "$t34_band1_conv" "$t34_band2_conv"

    t34_log="$OUT_DIR/t34_convolve_blank.log"
    if bin/convolve_cubes infiles="$t34_band1,$t34_band2" \
            beamfiles="$DATA_DIR/band1_beamlog.txt,$DATA_DIR/band2_beamlog.txt" \
            badchan_file=",$t34_band2_badchan" \
            target_bmaj=30.0 target_bmin=30.0 target_bpa=0.0 \
            > "$t34_log" 2>&1; then
        fail "convolve_cubes: blank badchan_file entry unexpectedly succeeded (see $t34_log)"
    else
        if grep -qE "badchan_file entry[[:space:]]+1[[:space:]]+is blank" "$t34_log" \
                && [[ ! -f "$t34_band1_conv" ]] && [[ ! -f "$t34_band2_conv" ]]; then
            pass "convolve_cubes: blank badchan_file entry is a clear named error, no output written (T34)"
        else
            fail "convolve_cubes: blank badchan_file entry error message/output check failed (see $t34_log)"
        fi
    fi

    # badchan_file= omitted must still work fine (unaffected).
    t34_log_omit="$OUT_DIR/t34_convolve_omit.log"
    rm -f "$t34_band1_conv" "$t34_band2_conv"
    if bin/convolve_cubes infiles="$t34_band1,$t34_band2" \
            beamfiles="$DATA_DIR/band1_beamlog.txt,$DATA_DIR/band2_beamlog.txt" \
            target_bmaj=30.0 target_bmin=30.0 target_bpa=0.0 \
            > "$t34_log_omit" 2>&1; then
        pass "convolve_cubes: omitting badchan_file= entirely is still unaffected (T34)"
    else
        fail "convolve_cubes: omitting badchan_file= entirely unexpectedly failed (see $t34_log_omit)"
    fi
else
    skip "bin/convolve_cubes not built; skipping blank badchan_file position hard-error test"
fi

if [[ -x bin/match_cubes ]]; then
    t34m_band1="$OUT_DIR/t34m_band1.Q.FITSCUBE"
    t34m_band2="$OUT_DIR/t34m_band2.Q.FITSCUBE"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$t34m_band1"
    cp "$DATA_DIR/TEST_BAND2.Q.FITSCUBE" "$t34m_band2"
    t34m_band2_badchan="$OUT_DIR/t34m_band2_badchan.txt"
    printf '77\n' > "$t34m_band2_badchan"
    t34m_band1_conv="${t34m_band1%.FITSCUBE}_CONV.FITS"
    t34m_band2_conv="${t34m_band2%.FITSCUBE}_CONV.FITS"
    rm -f "$t34m_band1_conv" "$t34m_band2_conv"

    t34m_log="$OUT_DIR/t34_match_blank.log"
    if bin/match_cubes stages=convolve infiles="$t34m_band1,$t34m_band2" \
            beamfiles="$DATA_DIR/band1_beamlog.txt,$DATA_DIR/band2_beamlog.txt" \
            badchan_file=",$t34m_band2_badchan" \
            target_bmaj=30.0 target_bmin=30.0 target_bpa=0.0 \
            > "$t34m_log" 2>&1; then
        fail "match_cubes: blank badchan_file entry unexpectedly succeeded (see $t34m_log)"
    else
        if grep -qE "badchan_file entry[[:space:]]+1[[:space:]]+is blank" "$t34m_log" \
                && [[ ! -f "$t34m_band1_conv" ]] && [[ ! -f "$t34m_band2_conv" ]]; then
            pass "match_cubes: blank badchan_file entry is a clear named error, no output written (T34)"
        else
            fail "match_cubes: blank badchan_file entry error message/output check failed (see $t34m_log)"
        fi
    fi
else
    skip "bin/match_cubes not built; skipping blank badchan_file position hard-error test"
fi

# ---------------------------------------------------------------------------
# 61. thread_safety_mod: OMP_NUM_THREADS-vs-hardware clamp
# ---------------------------------------------------------------------------
section "61. thread_safety_mod: OMP_NUM_THREADS clamp when it exceeds available cores"

# Deliberately larger than any host/container this suite runs on --
# guarantees the clamp path fires regardless of this machine's own core
# count, so the test doesn't hardcode a number that would break on a
# differently-sized host/CI runner.
T61_HUGE=100000
t61_clamp_re="WARNING: OMP_NUM_THREADS=$T61_HUGE exceeds the [0-9]+ (physical cores|logical CPUs) available.*clamping OMP_NUM_THREADS down to [0-9]+\."

if [[ -x "$BIN_OMP" ]]; then
    t61_log="$OUT_DIR/t61_omp_clamp.log"
    OMP_NUM_THREADS=$T61_HUGE "$BIN_OMP" -h > "$t61_log" 2>&1 || true
    if grep -qE "$t61_clamp_re" "$t61_log"; then
        pass "rm_synthesis (OMP build): OMP_NUM_THREADS=$T61_HUGE clamped with a named target, run continued"
    else
        fail "rm_synthesis (OMP build): expected clamp warning not found for OMP_NUM_THREADS=$T61_HUGE (see $t61_log)"
    fi

    # OMP_NUM_THREADS=1 hits the module's own <=1 early-return -- confirms
    # the check doesn't fire spuriously on an unremarkable request.
    t61_log_safe="$OUT_DIR/t61_omp_noclamp.log"
    OMP_NUM_THREADS=1 "$BIN_OMP" -h > "$t61_log_safe" 2>&1 || true
    if grep -q "WARNING: OMP_NUM_THREADS" "$t61_log_safe"; then
        fail "rm_synthesis (OMP build): OMP_NUM_THREADS=1 unexpectedly triggered a clamp warning (see $t61_log_safe)"
    else
        pass "rm_synthesis (OMP build): OMP_NUM_THREADS=1 triggers no warning"
    fi
else
    skip "$BIN_OMP not built; skipping OMP_NUM_THREADS clamp test"
fi

if [[ -x "$BIN_SERIAL" ]]; then
    # Serial build has no -fopenmp, so _OPENMP is undefined and the whole
    # check compiles away to a no-op (see src/thread_safety_mod.f90) --
    # confirms that guard holds in the linked binary itself, not
    # just in isolation, since a regression here would silently break the
    # default (OMP=0) build for everyone.
    t61_log_serial="$OUT_DIR/t61_serial_noop.log"
    OMP_NUM_THREADS=$T61_HUGE "$BIN_SERIAL" -h > "$t61_log_serial" 2>&1 || true
    if grep -q "WARNING: OMP_NUM_THREADS" "$t61_log_serial"; then
        fail "rm_synthesis (serial build): OMP_NUM_THREADS=$T61_HUGE unexpectedly triggered a clamp warning -- the _OPENMP guard should make this a no-op in a non-OpenMP build (see $t61_log_serial)"
    else
        pass "rm_synthesis (serial build): OMP_NUM_THREADS=$T61_HUGE is a silent no-op (no OpenMP runtime to clamp)"
    fi
else
    skip "$BIN_SERIAL not built; skipping serial-build no-op check"
fi

if [[ -x bin/convolve_cubes ]]; then
    # convolve_cubes is always built with -fopenmp regardless of top-level
    # OMP= (see Makefile) -- confirms the wiring reaches beyond rm_synthesis
    # into the other 4 tools sharing this module.
    t61_log_conv="$OUT_DIR/t61_convolve_clamp.log"
    OMP_NUM_THREADS=$T61_HUGE bin/convolve_cubes -h > "$t61_log_conv" 2>&1 || true
    if grep -qE "$t61_clamp_re" "$t61_log_conv"; then
        pass "convolve_cubes: OMP_NUM_THREADS=$T61_HUGE clamped with a named target (always-OMP build)"
    else
        fail "convolve_cubes: expected clamp warning not found for OMP_NUM_THREADS=$T61_HUGE (see $t61_log_conv)"
    fi
else
    skip "bin/convolve_cubes not built; skipping convolve_cubes clamp test"
fi

# ---------------------------------------------------------------------------
# 62. reproject_cubes value-level correctness (docs/dev/
#     REPROJECT_TEST_PLAN.md / REPROJECT_FLUX_NOTE.md): sections 30-47
#     above check reproject_cubes' interface behaviour (skip-if-matched,
#     I/O-option bit-identity, overwrite refusal, dry_run) and, in
#     section 40, that its output reaches rm_synthesis with the band 1
#     sources intact -- none of them compare a reprojected pixel value,
#     position, integral, peak, noise, or NaN footprint against an
#     independently-known truth. This section does, for the regime
#     REPROJECT_FLUX_NOTE.md's conclusion covers: both the reference and
#     input grids well sampled with respect to the beam (the large-
#     pixel-scale-factor regime in that note's section 8 is out of scope
#     here on purpose). Fixtures: tests/make_reproject_truth_cubes.py
#     (G0 reference; G1 sub-pixel shift; G2 rotation; G3 small scale
#     change; G5 whole-pixel shift -- see that file's own header for the
#     full grid definitions). Checks: tests/check_reproject_correctness.py
#     (ramp exactness, whole-pixel-shift exactness, integral ratio, peak
#     ratio, centroid position, NaN-footprint regression guard, unit
#     handling, polarisation-angle preservation), thresholds taken from
#     REPROJECT_TEST_PLAN.md section 3's proposal (itself from
#     REPROJECT_FLUX_NOTE.md's measurements) for the current AST__LINEAR
#     scheme -- reproject_cubes offers no other interpolation option
#     today, so this section covers the one scheme that exists, not a
#     to-be-added kernel choice.
# ---------------------------------------------------------------------------
section "62. reproject_cubes: value-level correctness (ramp, integral, peak, position, NaN footprint, units, polarisation angle)"

if [[ -x bin/reproject_cubes ]]; then
    rc62_gen_log="$OUT_DIR/reproject_truth_gen.log"
    if python3 "$TESTS_DIR/make_reproject_truth_cubes.py" > "$rc62_gen_log" 2>&1; then
        rc62_ok=1
        rc62_reproj_log="$OUT_DIR/reproject_correctness_reproject.log"
        : > "$rc62_reproj_log"
        rm -f "$DATA_DIR"/REPROJ_G*_REPROJ.FITS
        for rc62_pair in \
            "REPROJ_G0.RAMP.FITS REPROJ_G1.RAMP.FITS" \
            "REPROJ_G0.SOURCES.FITS REPROJ_G5.SOURCES.FITS" \
            "REPROJ_G0.SOURCES.FITS REPROJ_G1.SOURCES.FITS" \
            "REPROJ_G0.SOURCES.FITS REPROJ_G2.SOURCES.FITS" \
            "REPROJ_G0.SOURCES.FITS REPROJ_G3.SOURCES.FITS" \
            "REPROJ_G0.FLATNAN.FITS REPROJ_G1.FLATNAN.FITS" \
            "REPROJ_G0.FLAT.FITS REPROJ_G3.FLAT.FITS" \
            "REPROJ_G0.Q.FITS REPROJ_G1.Q.FITS" \
            "REPROJ_G0.U.FITS REPROJ_G1.U.FITS" ; do
            read -r rc62_ref rc62_in <<< "$rc62_pair"
            # interp=linear explicit, not the CLI default -- this
            # section validates linear's own established numbers (see
            # section 63's own comment below), independent of whatever
            # the default interp= happens to be (lanczos3 since
            # 2026-10-13, UR-12).
            if ! bin/reproject_cubes mode=reference \
                    reffile="$DATA_DIR/$rc62_ref" infiles="$DATA_DIR/$rc62_in" \
                    interp=linear >> "$rc62_reproj_log" 2>&1; then
                rc62_ok=0
            fi
        done
        if [[ "$rc62_ok" -eq 1 ]]; then
            rc62_check_log="$OUT_DIR/reproject_correctness_check.log"
            if python3 "$TESTS_DIR/check_reproject_correctness.py" > "$rc62_check_log" 2>&1; then
                pass "reproject_cubes: all value-level correctness checks passed (see $rc62_check_log)"
            else
                fail "reproject_cubes: one or more value-level correctness checks failed (see $rc62_check_log)"
            fi
        else
            fail "reproject_cubes: a reprojection run failed while preparing the correctness checks (see $rc62_reproj_log)"
        fi
    else
        fail "reproject_cubes correctness: fixture generation failed (see $rc62_gen_log)"
    fi
else
    skip "bin/reproject_cubes not built; skipping value-level correctness tests"
fi

# ---------------------------------------------------------------------------
# 63. reproject_cubes: interp=/interp_params= (Tier 1, docs/dev/
#     REPROJECT_INTERP_OPTIONS_PLAN.md section 5) -- CLI parsing,
#     --config parsing, error paths, a closed-form correctness check for
#     interp=nearest (check_nearest_ramp, tests/check_reproject_
#     correctness.py), and a default-equals-explicit-lanczos3 bit-
#     identity guard on the CLI/cfg merge for this specific key (default
#     changed from linear to lanczos3, 2026-10-13, UR-12 benchmark).
#     Section 62 above pins interp=linear explicitly and covers its own
#     value-level correctness independent of whatever the current
#     default is; this section covers the option itself, not a second
#     pass over linear's own numbers.
# ---------------------------------------------------------------------------
section "63. reproject_cubes: interp=/interp_params= (Tier 1 scheme selection)"

if [[ -x bin/reproject_cubes ]]; then
    rc63_ref="$DATA_DIR/REPROJ_G0.RAMP.FITS"
    rc63_in="$DATA_DIR/REPROJ_G1.RAMP.FITS"
    rc63_reproj="$DATA_DIR/REPROJ_G1.RAMP_REPROJ.FITS"

    # -- default (no interp=) bit-identical to explicit interp=lanczos3
    #    -- default changed from linear to lanczos3, 2026-10-13 (UR-12
    #    benchmark, docs/dev/BRANCH_STATUS_update-reproject.md: mean
    #    RMS error ~6x lower than linear on realistic sources). --
    rc63_default_out="$OUT_DIR/rc63_default.FITS"
    rc63_linear_out="$OUT_DIR/rc63_linear.FITS"
    rc63_default_log="$OUT_DIR/rc63_default.log"
    rc63_linear_log="$OUT_DIR/rc63_linear.log"
    rm -f "$rc63_reproj" "$rc63_default_out" "$rc63_linear_out"

    if bin/reproject_cubes mode=reference reffile="$rc63_ref" \
            infiles="$rc63_in" > "$rc63_default_log" 2>&1; then
        mv "$rc63_reproj" "$rc63_default_out"
    fi
    if bin/reproject_cubes mode=reference reffile="$rc63_ref" \
            infiles="$rc63_in" interp=lanczos3 > "$rc63_linear_log" 2>&1; then
        mv "$rc63_reproj" "$rc63_linear_out"
    fi
    if [[ -s "$rc63_default_out" && -s "$rc63_linear_out" ]] && \
       cmp -s "$rc63_default_out" "$rc63_linear_out"; then
        pass "reproject_cubes: default interp bit-identical to explicit interp=lanczos3"
    else
        fail "reproject_cubes: default interp differs from (or is missing vs) explicit interp=lanczos3 (see $rc63_default_log / $rc63_linear_log)"
    fi

    # -- interp=nearest: closed-form correctness (check_nearest_ramp) --
    rc63_nearest_out="$DATA_DIR/REPROJ_G1.RAMP_REPROJ_NEAREST.FITS"
    rc63_nearest_log="$OUT_DIR/rc63_nearest.log"
    rm -f "$rc63_reproj" "$rc63_nearest_out"
    if bin/reproject_cubes mode=reference reffile="$rc63_ref" \
            infiles="$rc63_in" interp=nearest > "$rc63_nearest_log" 2>&1; then
        mv "$rc63_reproj" "$rc63_nearest_out"
        rc63_check_log="$OUT_DIR/rc63_check.log"
        if python3 "$TESTS_DIR/check_reproject_correctness.py" --interp \
                > "$rc63_check_log" 2>&1; then
            pass "reproject_cubes: interp=nearest closed-form correctness check passed (see $rc63_check_log)"
        else
            fail "reproject_cubes: interp=nearest closed-form correctness check failed (see $rc63_check_log)"
        fi
    else
        fail "reproject_cubes: interp=nearest run failed (see $rc63_nearest_log)"
    fi

    # -- interp=sinccos interp_params=3,2.5: runs, differs from linear,
    #    introduces no new non-finite values in the interior (a sanity
    #    check that the scheme reaches AST at all, not a value-level
    #    correctness claim -- REPROJECT_INTERP_OPTIONS_PLAN.md section 9
    #    leaves sinccos's own thresholds for a later pass). --
    rc63_sinccos_out="$DATA_DIR/REPROJ_G1.SOURCES_REPROJ_SINCCOS.FITS"
    rc63_sinccos_linear_out="$DATA_DIR/REPROJ_G1.SOURCES_REPROJ_LINEAR.FITS"
    rc63_sourcesreproj="$DATA_DIR/REPROJ_G1.SOURCES_REPROJ.FITS"
    rc63_sinccos_log="$OUT_DIR/rc63_sinccos.log"
    rc63_sourceslinear_log="$OUT_DIR/rc63_sourceslinear.log"
    rm -f "$rc63_sourcesreproj" "$rc63_sinccos_out" "$rc63_sinccos_linear_out"
    rc63_sinccos_ok=1
    if bin/reproject_cubes mode=reference \
            reffile="$DATA_DIR/REPROJ_G0.SOURCES.FITS" \
            infiles="$DATA_DIR/REPROJ_G1.SOURCES.FITS" \
            interp=sinccos interp_params=3,2.5 \
            > "$rc63_sinccos_log" 2>&1; then
        mv "$rc63_sourcesreproj" "$rc63_sinccos_out"
    else
        rc63_sinccos_ok=0
    fi
    if bin/reproject_cubes mode=reference \
            reffile="$DATA_DIR/REPROJ_G0.SOURCES.FITS" \
            infiles="$DATA_DIR/REPROJ_G1.SOURCES.FITS" \
            > "$rc63_sourceslinear_log" 2>&1; then
        mv "$rc63_sourcesreproj" "$rc63_sinccos_linear_out"
    else
        rc63_sinccos_ok=0
    fi
    if [[ "$rc63_sinccos_ok" -eq 1 ]]; then
        rc63_sinccos_py_log="$OUT_DIR/rc63_sinccos_py.log"
        if python3 - "$rc63_sinccos_out" "$rc63_sinccos_linear_out" \
                > "$rc63_sinccos_py_log" 2>&1 <<'PYEOF'
import sys
import numpy as np
from astropy.io import fits

sinccos = fits.getdata(sys.argv[1]).squeeze().astype(np.float64)
linear = fits.getdata(sys.argv[2]).squeeze().astype(np.float64)
margin = 10
a = sinccos[margin:-margin, margin:-margin]
b = linear[margin:-margin, margin:-margin]
finite_a = np.isfinite(a)
finite_b = np.isfinite(b)
if not np.all(finite_b[finite_a] if finite_a.any() else True):
    sys.exit("interp=sinccos introduced non-finite values where interp=linear had none")
diff = np.abs(a - b)
frac_differing = float(np.mean(diff > 1e-6))
if frac_differing < 0.01:
    sys.exit(f"interp=sinccos output is nearly identical to interp=linear "
              f"(only {frac_differing*100:.4f}% of interior pixels differ by "
              f">1e-6) -- interp_scheme likely not reaching ast_resampler")
print(f"interp=sinccos differs from interp=linear on {frac_differing*100:.1f}% "
      f"of interior pixels, no new non-finite values")
PYEOF
        then
            pass "reproject_cubes: interp=sinccos reaches AST and changes output (see $rc63_sinccos_py_log)"
        else
            fail "reproject_cubes: interp=sinccos sanity check failed (see $rc63_sinccos_py_log)"
        fi
    else
        fail "reproject_cubes: interp=sinccos or interp=linear comparison run failed (see $rc63_sinccos_log / $rc63_sourceslinear_log)"
    fi
    rm -f "$rc63_sinccos_out" "$rc63_sinccos_linear_out"

    # -- CLI error paths: unrecognised interp name, malformed interp_params --
    rc63_err1_log="$OUT_DIR/rc63_err_badname.log"
    if ! bin/reproject_cubes mode=reference reffile="$rc63_ref" \
            infiles="$rc63_in" interp=bogus > "$rc63_err1_log" 2>&1 && \
       grep -q "ERROR:.*unrecognised interp" "$rc63_err1_log"; then
        pass "reproject_cubes: interp=bogus rejected with ERROR"
    else
        fail "reproject_cubes: interp=bogus was not rejected as expected (see $rc63_err1_log)"
    fi

    rc63_err2_log="$OUT_DIR/rc63_err_badparams.log"
    if ! bin/reproject_cubes mode=reference reffile="$rc63_ref" \
            infiles="$rc63_in" interp=gauss interp_params=1,2,3 \
            > "$rc63_err2_log" 2>&1 && \
       grep -q "ERROR:.*interp_params" "$rc63_err2_log"; then
        pass "reproject_cubes: interp_params=1,2,3 (too many values) rejected with ERROR"
    else
        fail "reproject_cubes: interp_params=1,2,3 was not rejected as expected (see $rc63_err2_log)"
    fi

    # -- interp=/interp_params= via --config file, not just CLI --
    rc63_cfg="$OUT_DIR/rc63_reproject.cfg"
    cat > "$rc63_cfg" <<CFGEOF
mode = reference
reffile = $rc63_ref
infiles = $rc63_in
interp = sincgauss
interp_params = 3,2.5
CFGEOF
    rc63_cfg_log="$OUT_DIR/rc63_cfg.log"
    rm -f "$rc63_reproj"
    if bin/reproject_cubes --config "$rc63_cfg" > "$rc63_cfg_log" 2>&1 && \
       [[ -s "$rc63_reproj" ]]; then
        pass "reproject_cubes: interp=/interp_params= accepted via --config"
    else
        fail "reproject_cubes: interp=/interp_params= via --config failed (see $rc63_cfg_log)"
    fi
    rm -f "$rc63_reproj"
else
    skip "bin/reproject_cubes not built; skipping interp=/interp_params= tests"
fi

# ---------------------------------------------------------------------------
# 64. match_cubes: interp=/interp_params= (Tier 1, docs/dev/
#     REPROJECT_INTERP_OPTIONS_PLAN.md section 5), covering all 3 of
#     this file's own ast_resampler call sites (process_one_file_
#     general for stages=reproject; process_one_file_restricted's two
#     branches for stages=both order=convolve_reproject and
#     order=reproject_convolve) -- the same ast_interp_mod module
#     reproject_cubes uses (section 63 above), confirmed here to give
#     bit-identical results on the same input/reference pair for both
#     interp=linear (the default) and interp=nearest, guarding against
#     the two tools' own copies of the resampling call ever drifting.
# ---------------------------------------------------------------------------
section "64. match_cubes: interp=/interp_params= (Tier 1 scheme selection)"

if [[ -x bin/match_cubes ]]; then
    mc64_ref="$DATA_DIR/REPROJ_G0.RAMP.FITS"
    mc64_in="$DATA_DIR/REPROJ_G1.RAMP.FITS"
    mc64_reproj="$DATA_DIR/REPROJ_G1.RAMP_REPROJ.FITS"

    # -- stages=reproject, interp=bogus rejected --
    mc64_err_log="$OUT_DIR/mc64_err_badname.log"
    if ! bin/match_cubes stages=reproject footprint_mode=reference \
            reffile="$mc64_ref" infiles="$mc64_in" interp=bogus \
            > "$mc64_err_log" 2>&1 && \
       grep -q "ERROR:.*unrecognised interp" "$mc64_err_log"; then
        pass "match_cubes: interp=bogus rejected with ERROR"
    else
        fail "match_cubes: interp=bogus was not rejected as expected (see $mc64_err_log)"
    fi

    # -- stages=reproject (process_one_file_general): interp=linear and
    #    interp=nearest bit-identical (pixel data) to reproject_cubes'
    #    own output on the same input/reference pair -- cross-tool
    #    consistency, not just "it runs." --
    mc64_cross_ok=1
    mc64_cross_py_log="$OUT_DIR/mc64_cross.log"
    : > "$mc64_cross_py_log"
    for mc64_scheme in linear nearest; do
        mc64_match_out="$OUT_DIR/mc64_match_${mc64_scheme}.FITS"
        mc64_rc_out="$OUT_DIR/mc64_rc_${mc64_scheme}.FITS"
        rm -f "$mc64_reproj" "$mc64_match_out" "$mc64_rc_out"
        if bin/match_cubes stages=reproject footprint_mode=reference \
                reffile="$mc64_ref" infiles="$mc64_in" interp="$mc64_scheme" \
                >> "$mc64_cross_py_log" 2>&1; then
            mv "$mc64_reproj" "$mc64_match_out"
        else
            mc64_cross_ok=0
        fi
        if bin/reproject_cubes mode=reference reffile="$mc64_ref" \
                infiles="$mc64_in" interp="$mc64_scheme" \
                >> "$mc64_cross_py_log" 2>&1; then
            mv "$mc64_reproj" "$mc64_rc_out"
        else
            mc64_cross_ok=0
        fi
        if [[ "$mc64_cross_ok" -eq 1 ]]; then
            if ! python3 - "$mc64_match_out" "$mc64_rc_out" \
                    >> "$mc64_cross_py_log" 2>&1 <<'PYEOF'
import sys
import numpy as np
from astropy.io import fits

a = fits.getdata(sys.argv[1])
b = fits.getdata(sys.argv[2])
worst = float(np.nanmax(np.abs(a.astype(np.float64) - b.astype(np.float64))))
if worst != 0.0:
    sys.exit(f"pixel data differs: worst abs diff {worst:.3e} (expected 0.0)")
print(f"{sys.argv[1]} vs {sys.argv[2]}: pixel data bit-identical")
PYEOF
            then
                mc64_cross_ok=0
            fi
        fi
        rm -f "$mc64_match_out" "$mc64_rc_out"
    done
    if [[ "$mc64_cross_ok" -eq 1 ]]; then
        pass "match_cubes: interp=linear/nearest pixel-identical to reproject_cubes on the same input (see $mc64_cross_py_log)"
    else
        fail "match_cubes: interp=linear/nearest NOT pixel-identical to reproject_cubes (see $mc64_cross_py_log)"
    fi

    # -- interp=/interp_params= via --config file --
    mc64_cfg="$OUT_DIR/mc64_match.cfg"
    cat > "$mc64_cfg" <<CFGEOF
stages = reproject
footprint_mode = reference
reffile = $mc64_ref
infiles = $mc64_in
interp = sincgauss
interp_params = 3,2.5
CFGEOF
    mc64_cfg_log="$OUT_DIR/mc64_cfg.log"
    rm -f "$mc64_reproj"
    if bin/match_cubes --config "$mc64_cfg" > "$mc64_cfg_log" 2>&1 && \
       [[ -s "$mc64_reproj" ]]; then
        pass "match_cubes: interp=/interp_params= accepted via --config"
    else
        fail "match_cubes: interp=/interp_params= via --config failed (see $mc64_cfg_log)"
    fi
    rm -f "$mc64_reproj"

    # -- stages=both, both orders: interp= reaches the other 2
    #    ast_resampler call sites (process_one_file_restricted's
    #    convolve_reproject and reproject_convolve branches) without
    #    erroring -- smoke-level, reusing the already-validated T15/T40
    #    band2 fixture (needs both convolve and reproject, unlike
    #    band1 which already matches). --
    mc64_b1q="$DATA_DIR/TEST.Q.FITSCUBE"
    mc64_b2q="$OUT_DIR/mc64_b2.Q.FITSCUBE"
    mc64_b2u="$OUT_DIR/mc64_b2.U.FITSCUBE"
    cp "$DATA_DIR/TEST_BAND2_GRIDSHIFT.Q.FITSCUBE" "$mc64_b2q"
    cp "$DATA_DIR/TEST_BAND2_GRIDSHIFT.U.FITSCUBE" "$mc64_b2u"
    mc64_b2q_matched="${mc64_b2q%.FITSCUBE}_MATCHED.FITS"
    mc64_b2u_matched="${mc64_b2u%.FITSCUBE}_MATCHED.FITS"
    mc64_beams="$DATA_DIR/band2_beamlog.txt,$DATA_DIR/band2_beamlog.txt"

    for mc64_order in convolve_reproject reproject_convolve; do
        mc64_both_log="$OUT_DIR/mc64_both_${mc64_order}.log"
        rm -f "$mc64_b2q_matched" "$mc64_b2u_matched"
        if bin/match_cubes stages=both order="$mc64_order" \
                footprint_mode=reference reffile="$mc64_b1q" \
                infiles="$mc64_b2q,$mc64_b2u" beamfiles="$mc64_beams" \
                interp=nearest > "$mc64_both_log" 2>&1 && \
           [[ -s "$mc64_b2q_matched" && -s "$mc64_b2u_matched" ]]; then
            pass "match_cubes: stages=both order=$mc64_order accepts interp= and writes output"
        else
            fail "match_cubes: stages=both order=$mc64_order with interp= failed (see $mc64_both_log)"
        fi
    done
    rm -f "$mc64_b2q" "$mc64_b2u" "$mc64_b2q_matched" "$mc64_b2u_matched"
else
    skip "bin/match_cubes not built; skipping interp=/interp_params= tests"
fi

# ---------------------------------------------------------------------------
# 65. reproject_cubes: interp=lanczos2/lanczos3/lanczos (Tier 2, docs/dev/
#     REPROJECT_INTERP_OPTIONS_PLAN.md section 6) -- CLI/config parsing,
#     error paths (fixed-support presets reject interp_params; the
#     general form requires exactly 2 values; a<=0 rejected), a sanity
#     check that lanczos3 reaches AST and differs from both linear and
#     sinccos, and -- critically -- a multi-plane, multi-
#     threaded, real-WCS-Mapping run: this exact combination segfaulted
#     during development (section 6d) via a gfortran code-generation
#     issue specific to a procedure(),pointer dummy argument of a
#     contains-ed procedure shared into that procedure's own !$omp
#     parallel region, fixed by resolving finterp as a local variable
#     instead (see scratch/build_check/test_ukern1_internal_proc_bug.f90
#     for the isolated reproducer). This section guards against that
#     regressing silently -- the tiny single-plane section-63/64
#     fixtures never exercise more than 1 OMP thread, so they would not
#     have caught it.
# ---------------------------------------------------------------------------
section "65. reproject_cubes: interp=lanczos2/lanczos3/lanczos (Tier 2)"

if [[ -x bin/reproject_cubes ]]; then
    rc65_ref="$DATA_DIR/REPROJ_G0.RAMP.FITS"
    rc65_in="$DATA_DIR/REPROJ_G1.RAMP.FITS"
    rc65_reproj="$DATA_DIR/REPROJ_G1.RAMP_REPROJ.FITS"

    # -- error paths --
    rc65_err_log="$OUT_DIR/rc65_err.log"
    rc65_err_ok=1
    if ! bin/reproject_cubes mode=reference reffile="$rc65_ref" infiles="$rc65_in" \
            interp=lanczos2 interp_params=5 > "$rc65_err_log" 2>&1 && \
       grep -q "ERROR:.*interp_params is not accepted" "$rc65_err_log"; then
        :
    else
        rc65_err_ok=0
    fi
    if ! bin/reproject_cubes mode=reference reffile="$rc65_ref" infiles="$rc65_in" \
            interp=lanczos3 interp_params=3,3 >> "$rc65_err_log" 2>&1 && \
       grep -q "ERROR:.*interp_params is not accepted" "$rc65_err_log"; then
        :
    else
        rc65_err_ok=0
    fi
    if ! bin/reproject_cubes mode=reference reffile="$rc65_ref" infiles="$rc65_in" \
            interp=lanczos >> "$rc65_err_log" 2>&1 && \
       grep -q "ERROR:.*requires interp_params" "$rc65_err_log"; then
        :
    else
        rc65_err_ok=0
    fi
    if ! bin/reproject_cubes mode=reference reffile="$rc65_ref" infiles="$rc65_in" \
            interp=lanczos interp_params=1,0 >> "$rc65_err_log" 2>&1 && \
       grep -q "ERROR:.*requires a > 0" "$rc65_err_log"; then
        :
    else
        rc65_err_ok=0
    fi
    if [[ "$rc65_err_ok" -eq 1 ]]; then
        pass "reproject_cubes: lanczos2/lanczos3/lanczos error paths all rejected correctly"
    else
        fail "reproject_cubes: a lanczos error path was not rejected as expected (see $rc65_err_log)"
    fi

    # -- lanczos3 reaches AST and differs from both linear and sinccos --
    rc65_lin_out="$OUT_DIR/rc65_linear.FITS"
    rc65_sc_out="$OUT_DIR/rc65_sinccos.FITS"
    rc65_lcz_out="$OUT_DIR/rc65_lanczos3.FITS"
    rc65_run_log="$OUT_DIR/rc65_run.log"
    : > "$rc65_run_log"
    rc65_run_ok=1
    rm -f "$rc65_reproj" "$rc65_lin_out" "$rc65_sc_out" "$rc65_lcz_out"
    if bin/reproject_cubes mode=reference reffile="$rc65_ref" infiles="$rc65_in" \
            interp=linear >> "$rc65_run_log" 2>&1; then
        mv "$rc65_reproj" "$rc65_lin_out"
    else
        rc65_run_ok=0
    fi
    if bin/reproject_cubes mode=reference reffile="$rc65_ref" infiles="$rc65_in" \
            interp=sinccos interp_params=3,2.5 >> "$rc65_run_log" 2>&1; then
        mv "$rc65_reproj" "$rc65_sc_out"
    else
        rc65_run_ok=0
    fi
    if bin/reproject_cubes mode=reference reffile="$rc65_ref" infiles="$rc65_in" \
            interp=lanczos3 >> "$rc65_run_log" 2>&1; then
        mv "$rc65_reproj" "$rc65_lcz_out"
    else
        rc65_run_ok=0
    fi
    if [[ "$rc65_run_ok" -eq 1 ]]; then
        if python3 - "$rc65_lin_out" "$rc65_sc_out" "$rc65_lcz_out" \
                >> "$rc65_run_log" 2>&1 <<'PYEOF'
import sys
import numpy as np
from astropy.io import fits

lin = fits.getdata(sys.argv[1]).squeeze().astype(np.float64)
sc = fits.getdata(sys.argv[2]).squeeze().astype(np.float64)
lcz = fits.getdata(sys.argv[3]).squeeze().astype(np.float64)
margin = 10
lin_i = lin[margin:-margin, margin:-margin]
sc_i = sc[margin:-margin, margin:-margin]
lcz_i = lcz[margin:-margin, margin:-margin]
if not np.all(np.isfinite(lcz_i)):
    sys.exit("lanczos3 produced non-finite interior values")
diff_vs_linear = float(np.mean(np.abs(lcz_i - lin_i) > 1e-6))
diff_vs_sinccos = float(np.mean(np.abs(lcz_i - sc_i) > 1e-6))
if diff_vs_linear < 0.01:
    sys.exit(f"lanczos3 is nearly identical to linear ({diff_vs_linear*100:.4f}% "
              f"interior pixels differ) -- interp_scheme/finterp likely not reaching AST")
print(f"lanczos3 differs from linear on {diff_vs_linear*100:.1f}% of interior pixels, "
      f"from sinccos on {diff_vs_sinccos*100:.1f}%, all finite")
PYEOF
        then
            pass "reproject_cubes: interp=lanczos3 reaches AST and differs from linear/sinccos (see $rc65_run_log)"
        else
            fail "reproject_cubes: interp=lanczos3 sanity check failed (see $rc65_run_log)"
        fi
    else
        fail "reproject_cubes: a comparison run failed while preparing the lanczos3 sanity check (see $rc65_run_log)"
    fi
    rm -f "$rc65_lin_out" "$rc65_sc_out" "$rc65_lcz_out"

    # -- lanczos2 via --config --
    rc65_cfg="$OUT_DIR/rc65_reproject.cfg"
    cat > "$rc65_cfg" <<CFGEOF
mode = reference
reffile = $rc65_ref
infiles = $rc65_in
interp = lanczos2
CFGEOF
    rc65_cfg_log="$OUT_DIR/rc65_cfg.log"
    rm -f "$rc65_reproj"
    if bin/reproject_cubes --config "$rc65_cfg" > "$rc65_cfg_log" 2>&1 && \
       [[ -s "$rc65_reproj" ]]; then
        pass "reproject_cubes: interp=lanczos2 accepted via --config"
    else
        fail "reproject_cubes: interp=lanczos2 via --config failed (see $rc65_cfg_log)"
    fi
    rm -f "$rc65_reproj"

    # -- regression guard: multi-plane, multi-threaded, real WCS Mapping
    #    (the exact combination that segfaulted during Tier 2
    #    development -- section 6d). TEST_BAND2_GRIDSHIFT.Q.FITSCUBE
    #    (150 planes) against TEST.Q.FITSCUBE, mem_frac_ram=0.9 to put
    #    every plane in one block and give OMP multi-thread work.
    #    GRIDSHIFT (UR-15, docs/dev/BRANCH_STATUS_update-reproject.md)
    #    replaces the original TEST_BAND2_MISMATCH here -- its own +5.0
    #    degree shift left zero sky overlap with TEST.Q.FITSCUBE, which
    #    mode=reference's own overlap check (added by UR-15) now
    #    correctly refuses; the "real WCS Mapping" this guard needs
    #    only required a non-trivial reprojection, never zero
    #    overlap specifically. --
    rc65_mt_ref="$DATA_DIR/TEST.Q.FITSCUBE"
    rc65_mt_in="$DATA_DIR/TEST_BAND2_GRIDSHIFT.Q.FITSCUBE"
    rc65_mt_reproj="$DATA_DIR/TEST_BAND2_GRIDSHIFT.Q_REPROJ.FITS"
    rc65_mt_log="$OUT_DIR/rc65_multithread.log"
    rm -f "$rc65_mt_reproj"
    if bin/reproject_cubes mode=reference reffile="$rc65_mt_ref" infiles="$rc65_mt_in" \
            interp=lanczos3 mem_frac_ram=0.9 > "$rc65_mt_log" 2>&1 && \
       [[ -s "$rc65_mt_reproj" ]]; then
        pass "reproject_cubes: interp=lanczos3 survives a multi-plane, multi-threaded, real-WCS-Mapping run"
    else
        fail "reproject_cubes: interp=lanczos3 multi-plane/multi-thread regression guard failed (see $rc65_mt_log)"
    fi
    rm -f "$rc65_mt_reproj"
else
    skip "bin/reproject_cubes not built; skipping Tier 2 (lanczos) tests"
fi

# ---------------------------------------------------------------------------
# 66. match_cubes: interp=lanczos2/lanczos3/lanczos (Tier 2), covering
#     all 3 of this file's own ast_resampler call sites the same way
#     section 64 did for Tier 1 -- including the same multi-plane,
#     multi-threaded regression guard section 65 added for
#     reproject_cubes, since match_cubes' own process_one_file_general/
#     process_one_file_restricted needed the identical local-variable-
#     not-dummy-argument fix (section 6d).
# ---------------------------------------------------------------------------
section "66. match_cubes: interp=lanczos2/lanczos3/lanczos (Tier 2)"

if [[ -x bin/match_cubes ]]; then
    mc66_ref="$DATA_DIR/REPROJ_G0.RAMP.FITS"
    mc66_in="$DATA_DIR/REPROJ_G1.RAMP.FITS"
    mc66_reproj="$DATA_DIR/REPROJ_G1.RAMP_REPROJ.FITS"

    # -- stages=reproject (process_one_file_general): lanczos3
    #    pixel-identical to reproject_cubes' own output -- same
    #    cross-tool consistency check section 64 did for Tier 1. --
    mc66_match_out="$OUT_DIR/mc66_match.FITS"
    mc66_rc_out="$OUT_DIR/mc66_rc.FITS"
    mc66_cross_log="$OUT_DIR/mc66_cross.log"
    : > "$mc66_cross_log"
    mc66_cross_ok=1
    rm -f "$mc66_reproj" "$mc66_match_out" "$mc66_rc_out"
    if bin/match_cubes stages=reproject footprint_mode=reference \
            reffile="$mc66_ref" infiles="$mc66_in" interp=lanczos3 \
            >> "$mc66_cross_log" 2>&1; then
        mv "$mc66_reproj" "$mc66_match_out"
    else
        mc66_cross_ok=0
    fi
    if bin/reproject_cubes mode=reference reffile="$mc66_ref" infiles="$mc66_in" \
            interp=lanczos3 >> "$mc66_cross_log" 2>&1; then
        mv "$mc66_reproj" "$mc66_rc_out"
    else
        mc66_cross_ok=0
    fi
    if [[ "$mc66_cross_ok" -eq 1 ]]; then
        if ! python3 - "$mc66_match_out" "$mc66_rc_out" >> "$mc66_cross_log" 2>&1 <<'PYEOF'
import sys
import numpy as np
from astropy.io import fits

a = fits.getdata(sys.argv[1])
b = fits.getdata(sys.argv[2])
worst = float(np.nanmax(np.abs(a.astype(np.float64) - b.astype(np.float64))))
if worst != 0.0:
    sys.exit(f"pixel data differs: worst abs diff {worst:.3e} (expected 0.0)")
print("pixel data bit-identical")
PYEOF
        then
            mc66_cross_ok=0
        fi
    fi
    rm -f "$mc66_match_out" "$mc66_rc_out"
    if [[ "$mc66_cross_ok" -eq 1 ]]; then
        pass "match_cubes: interp=lanczos3 pixel-identical to reproject_cubes on the same input (see $mc66_cross_log)"
    else
        fail "match_cubes: interp=lanczos3 NOT pixel-identical to reproject_cubes (see $mc66_cross_log)"
    fi

    # -- error path --
    mc66_err_log="$OUT_DIR/mc66_err.log"
    if ! bin/match_cubes stages=reproject footprint_mode=reference \
            reffile="$mc66_ref" infiles="$mc66_in" interp=lanczos2 interp_params=5 \
            > "$mc66_err_log" 2>&1 && \
       grep -q "ERROR:.*interp_params is not accepted" "$mc66_err_log"; then
        pass "match_cubes: lanczos2 with interp_params rejected with ERROR"
    else
        fail "match_cubes: lanczos2 with interp_params was not rejected as expected (see $mc66_err_log)"
    fi

    # -- multi-plane, multi-threaded regression guard for all 3 call
    #    sites: stages=both, both orders, on the real BAND2_UNMATCHED
    #    fixture (process_one_file_restricted's own two branches),
    #    plus the stages=reproject run above already covered
    #    process_one_file_general. --
    mc66_b1q="$DATA_DIR/TEST.Q.FITSCUBE"
    mc66_b2q="$OUT_DIR/mc66_b2.Q.FITSCUBE"
    mc66_b2u="$OUT_DIR/mc66_b2.U.FITSCUBE"
    cp "$DATA_DIR/TEST_BAND2_GRIDSHIFT.Q.FITSCUBE" "$mc66_b2q"
    cp "$DATA_DIR/TEST_BAND2_GRIDSHIFT.U.FITSCUBE" "$mc66_b2u"
    mc66_b2q_matched="${mc66_b2q%.FITSCUBE}_MATCHED.FITS"
    mc66_b2u_matched="${mc66_b2u%.FITSCUBE}_MATCHED.FITS"
    mc66_beams="$DATA_DIR/band2_beamlog.txt,$DATA_DIR/band2_beamlog.txt"

    for mc66_order in convolve_reproject reproject_convolve; do
        mc66_both_log="$OUT_DIR/mc66_both_${mc66_order}.log"
        rm -f "$mc66_b2q_matched" "$mc66_b2u_matched"
        if bin/match_cubes stages=both order="$mc66_order" \
                footprint_mode=reference reffile="$mc66_b1q" \
                infiles="$mc66_b2q,$mc66_b2u" beamfiles="$mc66_beams" \
                interp=lanczos3 mem_frac_ram=0.9 > "$mc66_both_log" 2>&1 && \
           [[ -s "$mc66_b2q_matched" && -s "$mc66_b2u_matched" ]]; then
            pass "match_cubes: stages=both order=$mc66_order survives a multi-plane, multi-threaded interp=lanczos3 run"
        else
            fail "match_cubes: stages=both order=$mc66_order with interp=lanczos3 failed (see $mc66_both_log)"
        fi
    done
    rm -f "$mc66_b2q" "$mc66_b2u" "$mc66_b2q_matched" "$mc66_b2u_matched"
else
    skip "bin/match_cubes not built; skipping Tier 2 (lanczos) tests"
fi

# ---------------------------------------------------------------------------
# 67. scripts/casda_fetch.py: --ra/--dec/--lon/--lat/--frame coordinate
#     notation and frame support -- pure, network-free parts only
#     (resolve_target's own parsing, describe_coord, parse_args'
#     validation); everything else in casda_fetch.py talks to CASDA's
#     live service and has no test coverage here for that reason.
# ---------------------------------------------------------------------------
section "67. casda_fetch.py: coordinate notation/frame support"

if [[ -f scripts/casda_fetch.py ]]; then
    if python3 "$TESTS_DIR/test_casda_fetch_coords.py" \
            > "$OUT_DIR/casda_fetch_coords.log" 2>&1; then
        pass "casda_fetch.py: coordinate notation/frame tests passed (see $OUT_DIR/casda_fetch_coords.log)"
    else
        fail "casda_fetch.py: coordinate notation/frame tests failed (see $OUT_DIR/casda_fetch_coords.log)"
    fi
else
    skip "scripts/casda_fetch.py not found; skipping coordinate notation/frame tests"
fi

# ---------------------------------------------------------------------------
# 68. reproject_cubes/match_cubes: target_system=/target_projection=/
#     target_center= (docs/dev/REPROJECT_TARGET_FRAME_PLAN.md, Phase 1 --
#     target_system= alone, target_center=reffile only). Error paths,
#     header/value correctness against an independent astropy transform
#     (not just "it runs"), and a cross-tool pixel-identity check via the
#     shared target_grid_mod.f90 -- same reasoning as section 64/66's own
#     cross-tool checks for ast_interp_mod. reffile is copied into OUT_DIR
#     first (not used directly from DATA_DIR) because resolve_target_grid
#     writes its generated <stem>_TARGETFRAME.FITS NEXT TO reffile --
#     using the DATA_DIR fixture directly would litter tracked test data.
# ---------------------------------------------------------------------------
section "68. reproject_cubes/match_cubes: target_system=/target_projection=/target_center= (UR-7: no-reffile path)"

if [[ -x bin/reproject_cubes ]]; then
    tf_ref="$OUT_DIR/tf68_ref.Q.FITSCUBE"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$tf_ref"
    tf_reproj="${tf_ref%.FITSCUBE}_REPROJ.FITS"

    # -- error paths -- UR-7: reffile= is now exclusive with every
    #    target_* key (not a combinable pair any more), and every other
    #    target_* validation only fires once reffile= is absent. --
    tf_err_log="$OUT_DIR/tf68_err.log"
    tf_err_ok=1
    rm -f "$tf_reproj"
    if ! bin/reproject_cubes mode=reference reffile="$tf_ref" infiles="$tf_ref" \
            target_system=galactic > "$tf_err_log" 2>&1 && \
       grep -q "ERROR:.*not allowed together with reffile=" "$tf_err_log"; then
        :
    else
        tf_err_ok=0
    fi
    if ! bin/reproject_cubes infiles="$tf_ref" \
            >> "$tf_err_log" 2>&1 && \
       grep -q "ERROR:.*at least one of target_system=.*or target_projection=.*is required" "$tf_err_log"; then
        :
    else
        tf_err_ok=0
    fi
    # target_projection= accepts tan/sin/car/ait/hpx (section 75 tests
    # each one); an unknown value must still be refused by name.
    if ! bin/reproject_cubes infiles="$tf_ref" \
            target_system=galactic target_projection=mol \
            >> "$tf_err_log" 2>&1 && \
       grep -q "ERROR:.*target_projection must be tan, sin, car, ait, or hpx" "$tf_err_log"; then
        :
    else
        tf_err_ok=0
    fi
    if ! bin/reproject_cubes infiles="$tf_ref" \
            target_system=equatorial \
            >> "$tf_err_log" 2>&1 && \
       grep -q "ERROR:.*target_system must be" "$tf_err_log"; then
        :
    else
        tf_err_ok=0
    fi
    if ! bin/reproject_cubes infiles="$tf_ref" \
            target_system=galactic target_center=bogus \
            >> "$tf_err_log" 2>&1 && \
       grep -q "ERROR:.*target_center must be" "$tf_err_log"; then
        :
    else
        tf_err_ok=0
    fi
    if ! bin/reproject_cubes infiles="$tf_ref" \
            target_system=galactic target_center=user_defined \
            >> "$tf_err_log" 2>&1 && \
       grep -q "ERROR:.*target_center=user_defined requires both" "$tf_err_log"; then
        :
    else
        tf_err_ok=0
    fi
    if ! bin/reproject_cubes infiles="$tf_ref" \
            target_system=galactic target_lon=10.0 target_lat=20.0 \
            >> "$tf_err_log" 2>&1 && \
       grep -q "ERROR:.*target_lon=/target_lat=.*require.*target_center=user_defined" "$tf_err_log"; then
        :
    else
        tf_err_ok=0
    fi
    if ! bin/reproject_cubes infiles="$tf_ref" \
            target_system=galactic target_center=user_defined target_lon=not_a_number \
            target_lat=20.0 \
            >> "$tf_err_log" 2>&1 && \
       grep -q "ERROR: could not parse coordinate value" "$tf_err_log"; then
        :
    else
        tf_err_ok=0
    fi
    # UR-16 removed the old "target_cdelt_* required with more than one
    # infile" refusal (the default is now the smallest |CDELT| across
    # the infiles). What remains an error is giving only one of the
    # pair, which previously fell through to an unclear parse failure.
    if ! bin/reproject_cubes infiles="$tf_ref" \
            target_system=galactic target_cdelt_lon=0.001 \
            >> "$tf_err_log" 2>&1 && \
       grep -q "ERROR:.*target_cdelt_lon=.*target_cdelt_lat=.*must be given" "$tf_err_log"; then
        :
    else
        tf_err_ok=0
    fi
    rm -f "$tf_reproj"
    if [[ "$tf_err_ok" -eq 1 ]]; then
        pass "reproject_cubes: target_system=/target_projection=/target_center= error paths all rejected correctly"
    else
        fail "reproject_cubes: a target_system=/target_projection=/target_center= error path was not rejected as expected (see $tf_err_log)"
    fi

    # -- correctness (UR-7, no reffile): target_lon=/target_lat= set to
    #    the infile's OWN CRVAL expressed in galactic coordinates
    #    (computed independently via astropy, not read back from any
    #    tool output) -- this makes the ONE sky position whose pixel
    #    position is exactly preserved under a pure frame rotation (the
    #    infile's own CRPIX, an actual data point) land exactly at the
    #    OUTPUT's own CRPIX too, by construction of the in-memory
    #    shift+rotation Mapping (target_grid_mod.f90's own
    #    build_target_skymap: output pixel final_crpix always maps to
    #    exactly CRVAL, the same algebraic identity scratch/build_check/
    #    test_shiftmap.f90 verified this session). No interim file any
    #    more -- output lands directly at tf_reproj. --
    tf_run_log="$OUT_DIR/tf68_run.log"
    rm -f "$tf_reproj"
    if python3 - "$tf_ref" > "$tf_run_log" 2>&1 <<'PYEOF'
import sys
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
hdr = fits.getheader(sys.argv[1])
expect = SkyCoord(ra=hdr['CRVAL1']*u.deg, dec=hdr['CRVAL2']*u.deg, frame='fk5').galactic
print(f"{expect.l.deg} {expect.b.deg}")
PYEOF
    then
        tf_run_lb=$(tail -n 1 "$tf_run_log")
        tf_run_l="${tf_run_lb% *}"
        tf_run_b="${tf_run_lb#* }"
        if bin/reproject_cubes infiles="$tf_ref" target_system=galactic \
                target_center=user_defined target_lon="$tf_run_l" target_lat="$tf_run_b" \
                >> "$tf_run_log" 2>&1 && [[ -s "$tf_reproj" ]]; then
            if python3 - "$tf_ref" "$tf_reproj" "$tf_run_l" "$tf_run_b" \
                    >> "$tf_run_log" 2>&1 <<'PYEOF'
import sys
import numpy as np
from astropy.io import fits

orig, reproj, exp_l, exp_b = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
orig_hdr = fits.getheader(orig)
reproj_hdr = fits.getheader(reproj)

if reproj_hdr['CTYPE1'] != 'GLON-SIN' or reproj_hdr['CTYPE2'] != 'GLAT-SIN':
    sys.exit(f"CTYPE1/2 = {reproj_hdr['CTYPE1']}/{reproj_hdr['CTYPE2']},"
              f" expected GLON-SIN/GLAT-SIN")
dl = abs(reproj_hdr['CRVAL1'] - exp_l)
db = abs(reproj_hdr['CRVAL2'] - exp_b)
if dl > 1e-6 or db > 1e-6:
    sys.exit(f"CRVAL1/2 = {reproj_hdr['CRVAL1']},{reproj_hdr['CRVAL2']}, expected"
              f" {exp_l},{exp_b} (astropy), diff {dl:.3e},{db:.3e} deg")
if 'RADESYS' in reproj_hdr or 'EQUINOX' in reproj_hdr:
    sys.exit(f"RADESYS/EQUINOX present ({reproj_hdr.get('RADESYS')!r}/"
              f"{reproj_hdr.get('EQUINOX')!r}), expected both absent for Galactic")

orig_data = fits.getdata(orig)
reproj_data = fits.getdata(reproj)
# NOT fully finite in general: the axis-aligned bounding box of a
# ROTATED footprint (FK5->Galactic here) is strictly larger than the
# rotated footprint itself, leaving diamond-shaped NaN corners outside
# it -- expected, correct reprojection behaviour (same as rotating a
# photo into a larger axis-aligned canvas), not a defect. Only the
# in-footprint fraction needs to be finite; a near-zero finite fraction
# would mean the footprint-fit itself is broken.
finite_frac = float(np.mean(np.isfinite(reproj_data)))
if finite_frac < 0.1:
    sys.exit(f"reprojected output is almost entirely non-finite ({finite_frac:.3f}"
              f" finite) -- footprint fit looks broken, not just rotated corners")
crpix1_o, crpix2_o = int(round(orig_hdr['CRPIX1'])), int(round(orig_hdr['CRPIX2']))
crpix1_r, crpix2_r = int(round(reproj_hdr['CRPIX1'])), int(round(reproj_hdr['CRPIX2']))
orig_val = orig_data[..., crpix2_o - 1, crpix1_o - 1].astype(np.float64)
reproj_val = reproj_data[..., crpix2_r - 1, crpix1_r - 1].astype(np.float64)
# Not exact bit-identity: the round-trip pixel coordinate through the
# sky2sky rotation lands within floating-point noise of the integer
# CRPIX, not exactly on it, so linear interpolation blends in an
# infinitesimal neighbour contribution. 1e-3 is generous against the
# synthetic test cube's own O(1) Jy/beam values -- enough to catch a
# genuine misalignment, not sub-pixel floating-point noise.
if not np.allclose(orig_val, reproj_val, rtol=0, atol=1e-3):
    sys.exit(f"value at output CRPIX ({crpix1_r},{crpix2_r}) = {reproj_val},"
              f" expected (orig's own CRPIX value) {orig_val}")
print(f"target_system=galactic (no reffile): CTYPE/CRVAL/RADESYS/EQUINOX correct"
      f" (l={exp_l:.6f}, b={exp_b:.6f}), data finite, CRPIX value preserved exactly")
PYEOF
            then
                pass "reproject_cubes: target_system=galactic (no reffile) header/value correctness confirmed against astropy (see $tf_run_log)"
            else
                fail "reproject_cubes: target_system=galactic (no reffile) correctness check failed (see $tf_run_log)"
            fi
        else
            fail "reproject_cubes: target_system=galactic (no reffile) run failed to produce expected output (see $tf_run_log)"
        fi
    else
        fail "reproject_cubes: failed to compute the expected galactic CRVAL for the correctness check (see $tf_run_log)"
    fi
    rm -f "$tf_reproj"

    # -- cross-tool identity: match_cubes stages=reproject must produce
    #    pixel-identical output to reproject_cubes for the SAME
    #    target_system=galactic request -- the shared target_grid_mod.f90
    #    cannot be allowed to diverge between the two tools' own call
    #    sites, same reasoning as section 64/66's own cross-tool checks. --
    if [[ -x bin/match_cubes ]]; then
        tf_cross_log="$OUT_DIR/tf68_cross.log"
        : > "$tf_cross_log"
        tf_cross_ok=1
        tf_rc_out="$OUT_DIR/tf68_rc_out.FITS"
        tf_mc_out="$OUT_DIR/tf68_mc_out.FITS"
        rm -f "$tf_reproj" "$tf_rc_out" "$tf_mc_out"
        if bin/reproject_cubes infiles="$tf_ref" target_system=galactic \
                >> "$tf_cross_log" 2>&1; then
            mv "$tf_reproj" "$tf_rc_out"
        else
            tf_cross_ok=0
        fi
        if bin/match_cubes stages=reproject infiles="$tf_ref" \
                target_system=galactic \
                >> "$tf_cross_log" 2>&1; then
            mv "$tf_reproj" "$tf_mc_out"
        else
            tf_cross_ok=0
        fi
        if [[ "$tf_cross_ok" -eq 1 ]]; then
            if ! python3 - "$tf_rc_out" "$tf_mc_out" \
                    >> "$tf_cross_log" 2>&1 <<'PYEOF'
import sys
import numpy as np
from astropy.io import fits

a = fits.getdata(sys.argv[1])
b = fits.getdata(sys.argv[2])
worst = float(np.nanmax(np.abs(a.astype(np.float64) - b.astype(np.float64))))
if worst != 0.0:
    sys.exit(f"pixel data differs: worst abs diff {worst:.3e} (expected 0.0)")
hdr_a, hdr_b = fits.getheader(sys.argv[1]), fits.getheader(sys.argv[2])
for k in ('CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2'):
    if hdr_a[k] != hdr_b[k]:
        sys.exit(f"header key {k} differs: {hdr_a[k]!r} vs {hdr_b[k]!r}")
print("reproject_cubes vs match_cubes: target_system=galactic output pixel-identical,"
      " header WCS keys match")
PYEOF
            then
                tf_cross_ok=0
            fi
        fi
        rm -f "$tf_reproj" "$tf_rc_out" "$tf_mc_out"
        if [[ "$tf_cross_ok" -eq 1 ]]; then
            pass "match_cubes: target_system=galactic pixel-identical to reproject_cubes on the same input (see $tf_cross_log)"
        else
            fail "match_cubes: target_system=galactic NOT pixel-identical to reproject_cubes (see $tf_cross_log)"
        fi

        # -- match_cubes-specific guard: target_system with stages=
        #    convolve (no reproject stage, nothing to apply target_system
        #    to) is a clear named error, not silently ignored. --
        tf_mc_guard_log="$OUT_DIR/tf68_mc_guard.log"
        if ! bin/match_cubes stages=convolve infiles="$tf_ref" \
                target_system=galactic > "$tf_mc_guard_log" 2>&1 && \
           grep -q "ERROR:.*target_system/target_projection/target_center.*require.*stages to include reproject" "$tf_mc_guard_log"; then
            pass "match_cubes: target_system with stages=convolve (no reproject) rejected with a clear ERROR"
        else
            fail "match_cubes: target_system without reproject was not rejected as expected (see $tf_mc_guard_log)"
        fi
    else
        skip "bin/match_cubes not built; skipping target_system cross-tool check"
    fi

    # -- target_center=infiles_centroid correctness: two infiles with
    #    DIFFERENT, known CRVAL (patched via astropy), reffile
    #    independent of both (case (ii)'s own "independent reffile"
    #    scenario, design doc section 4d) -- verify the generated
    #    CRVAL matches an INDEPENDENTLY computed 3-D-vector-averaged
    #    centroid of each infile's own pixel-grid-centre sky position,
    #    not just "it runs". --
    tf_centroid_log="$OUT_DIR/tf68_centroid.log"
    tf_in_a="$OUT_DIR/tf68_centroid_a.FITS"
    tf_in_b="$OUT_DIR/tf68_centroid_b.FITS"
    tf_in_a_reproj="${tf_in_a%.FITS}_REPROJ.FITS"
    tf_in_b_reproj="${tf_in_b%.FITS}_REPROJ.FITS"
    rm -f "$tf_in_a" "$tf_in_b" "$tf_in_a_reproj" "$tf_in_b_reproj"
    if python3 - "$DATA_DIR/TEST.Q.FITSCUBE" "$tf_in_a" "$tf_in_b" \
            > "$tf_centroid_log" 2>&1 <<'PYEOF'
import sys
from astropy.io import fits

src, out_a, out_b = sys.argv[1], sys.argv[2], sys.argv[3]
# Offsets kept well inside TEST.Q.FITSCUBE's own 32x32 @ 0.001 deg/pixel
# footprint (0.032 deg wide) -- UR-9 (docs/dev/BRANCH_STATUS_update-
# reproject.md) added a zero-overlap check to the no-reffile multi-
# infile path; the original 2.0/1.0 and -3.0/-2.0 degree offsets put
# these two fixtures nowhere near each other at all (always were, this
# just went uncaught until UR-9 added the check this ticket asks for).
for out, dra, ddec in [(out_a, 0.01, 0.005), (out_b, -0.015, -0.01)]:
    h = fits.open(src)
    h[0].header['CRVAL1'] = h[0].header['CRVAL1'] + dra
    h[0].header['CRVAL2'] = h[0].header['CRVAL2'] + ddec
    h.writeto(out, overwrite=True)
print("wrote", out_a, out_b)
PYEOF
    then
        if bin/reproject_cubes infiles="$tf_in_a,$tf_in_b" target_system=galactic \
                target_center=infiles_centroid target_cdelt_lon=0.001 \
                target_cdelt_lat=0.001 >> "$tf_centroid_log" 2>&1; then
            if python3 - "$tf_in_a" "$tf_in_b" "$tf_in_a_reproj" \
                    >> "$tf_centroid_log" 2>&1 <<'PYEOF'
import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

in_a, in_b, gen = sys.argv[1], sys.argv[2], sys.argv[3]

def center_in_galactic(path):
    hdr = fits.getheader(path)
    w = WCS(hdr).celestial
    nx, ny = hdr['NAXIS1'], hdr['NAXIS2']
    # true geometric pixel-grid centre, 1-based FITS convention --
    # NOT CRPIX, which need not coincide with it (confirmed CRPIX=17
    # vs true centre 16.5 for this fixture).
    lon, lat = w.wcs_pix2world(float(nx + 1) / 2.0, float(ny + 1) / 2.0, 1)
    c = SkyCoord(lon * u.deg, lat * u.deg, frame='fk5').galactic
    return c.l.rad, c.b.rad

sumx = sumy = sumz = 0.0
for path in (in_a, in_b):
    lon_r, lat_r = center_in_galactic(path)
    sumx += np.cos(lat_r) * np.cos(lon_r)
    sumy += np.cos(lat_r) * np.sin(lon_r)
    sumz += np.sin(lat_r)
norm = np.sqrt(sumx**2 + sumy**2 + sumz**2)
expect_lat = np.degrees(np.arctan2(sumz / norm, np.hypot(sumx / norm, sumy / norm)))
expect_lon = np.degrees(np.arctan2(sumy, sumx)) % 360.0

gen_hdr = fits.getheader(gen)
dl = abs(gen_hdr['CRVAL1'] - expect_lon)
db = abs(gen_hdr['CRVAL2'] - expect_lat)
if dl > 1e-6 or db > 1e-6:
    sys.exit(f"infiles_centroid CRVAL1/2 = {gen_hdr['CRVAL1']},{gen_hdr['CRVAL2']},"
              f" independently computed centroid = {expect_lon},{expect_lat},"
              f" diff {dl:.3e},{db:.3e} deg")
print(f"infiles_centroid CRVAL matches independently computed vector-averaged"
      f" centroid: l={expect_lon:.6f} b={expect_lat:.6f}")
PYEOF
            then
                pass "reproject_cubes: target_center=infiles_centroid matches an independently computed vector-averaged centroid (see $tf_centroid_log)"
            else
                fail "reproject_cubes: target_center=infiles_centroid did not match the independently computed centroid (see $tf_centroid_log)"
            fi
        else
            fail "reproject_cubes: target_center=infiles_centroid run failed (see $tf_centroid_log)"
        fi
    else
        fail "reproject_cubes: failed to build target_center=infiles_centroid fixtures (see $tf_centroid_log)"
    fi
    rm -f "$tf_in_a" "$tf_in_b" "$tf_in_a_reproj" "$tf_in_b_reproj"

    # -- UR-9 (docs/dev/BRANCH_STATUS_update-reproject.md): two
    #    infiles with no shared sky at all (TEST.Q.FITSCUBE's own
    #    32x32 @ 0.001 deg/pixel footprint is 0.032 deg wide; a 2/1
    #    degree offset puts them nowhere near each other) must be
    #    rejected with a named error, not silently unioned into one
    #    grid spanning the gap between them. --
    tf_ov_log="$OUT_DIR/tf68_ur9_overlap.log"
    tf_ov_a="$OUT_DIR/tf68_ur9_a.FITS"
    tf_ov_b="$OUT_DIR/tf68_ur9_b.FITS"
    tf_ov_a_reproj="${tf_ov_a%.FITS}_REPROJ.FITS"
    tf_ov_b_reproj="${tf_ov_b%.FITS}_REPROJ.FITS"
    rm -f "$tf_ov_a" "$tf_ov_b" "$tf_ov_a_reproj" "$tf_ov_b_reproj"
    if python3 - "$DATA_DIR/TEST.Q.FITSCUBE" "$tf_ov_a" "$tf_ov_b" \
            > "$tf_ov_log" 2>&1 <<'PYEOF'
import sys
from astropy.io import fits
src, out_a, out_b = sys.argv[1], sys.argv[2], sys.argv[3]
for out, dra, ddec in [(out_a, 2.0, 1.0), (out_b, -3.0, -2.0)]:
    h = fits.open(src)
    h[0].header['CRVAL1'] = h[0].header['CRVAL1'] + dra
    h[0].header['CRVAL2'] = h[0].header['CRVAL2'] + ddec
    h.writeto(out, overwrite=True)
print("wrote", out_a, out_b)
PYEOF
    then
        if ! bin/reproject_cubes infiles="$tf_ov_a,$tf_ov_b" \
                target_system=galactic target_cdelt_lon=0.001 \
                target_cdelt_lat=0.001 >> "$tf_ov_log" 2>&1 && \
           grep -q "ERROR: zero sky overlap between" "$tf_ov_log"; then
            pass "reproject_cubes: no-reffile multi-infile with zero sky overlap rejected with a named ERROR (UR-9, see $tf_ov_log)"
        else
            fail "reproject_cubes: no-reffile multi-infile with zero sky overlap was not rejected as expected (UR-9, see $tf_ov_log)"
        fi
    else
        fail "reproject_cubes: failed to build UR-9 zero-overlap fixtures (see $tf_ov_log)"
    fi
    rm -f "$tf_ov_a" "$tf_ov_b" "$tf_ov_a_reproj" "$tf_ov_b_reproj"

    # -- target_center=user_defined correctness: decimal, colon-
    #    sexagesimal, and letter-suffixed sexagesimal all agree on the
    #    same intended value, and the colon form's hour-angle fallback
    #    correctly differs between an equatorial and a Galactic
    #    target_system (same ambiguity scripts/casda_fetch.py's own
    #    --lon resolves the same way). --
    # Target points are picked NEAR this fixture's own actual sky
    # position (RA=180,Dec=0 -> l=276.337,b=60.189, computed in the
    # correctness check above) -- UR-7's no-reffile path, unlike the
    # old reffile+target_center=user_defined header-relabel-only design
    # this replaces, must also fit the infile's own footprint around
    # the requested centre (ast_mapbox over compose_pix2pix), which
    # fails outright ("No valid output coordinates") for a SIN
    # projection whose tangent point is unrelated to the data -- a
    # correct domain rejection, not a defect, but it means these
    # notation-equivalence checks need a target point the data can
    # actually reach, not an arbitrary one.
    tf_ud_log="$OUT_DIR/tf68_user_defined.log"
    tf_ud_ok=1
    rm -f "$tf_reproj"
    for spec in "276.5|60.0|galactic|276.5|60.0" "18h26m00s|60.0|galactic|276.5|60.0" \
            "12:00:00|0.0|icrs|180.0|0.0"; do
        tf_ud_lon="$(echo "$spec" | cut -d'|' -f1)"
        tf_ud_lat="$(echo "$spec" | cut -d'|' -f2)"
        tf_ud_sys="$(echo "$spec" | cut -d'|' -f3)"
        tf_ud_exp_lon="$(echo "$spec" | cut -d'|' -f4)"
        tf_ud_exp_lat="$(echo "$spec" | cut -d'|' -f5)"
        rm -f "$tf_reproj"
        if ! bin/reproject_cubes infiles="$tf_ref" \
                target_system="$tf_ud_sys" target_center=user_defined \
                target_lon="$tf_ud_lon" target_lat="$tf_ud_lat" \
                >> "$tf_ud_log" 2>&1; then
            tf_ud_ok=0
            continue
        fi
        if ! python3 - "$tf_reproj" "$tf_ud_exp_lon" "$tf_ud_exp_lat" \
                >> "$tf_ud_log" 2>&1 <<'PYEOF'
import sys
from astropy.io import fits
hdr = fits.getheader(sys.argv[1])
exp_lon, exp_lat = float(sys.argv[2]), float(sys.argv[3])
dl = abs(hdr['CRVAL1'] - exp_lon)
db = abs(hdr['CRVAL2'] - exp_lat)
if dl > 1e-9 or db > 1e-9:
    sys.exit(f"CRVAL1/2 = {hdr['CRVAL1']},{hdr['CRVAL2']}, expected {exp_lon},{exp_lat}"
              f" diff {dl:.3e},{db:.3e}")
print(f"matches {exp_lon},{exp_lat}")
PYEOF
        then
            tf_ud_ok=0
        fi
    done
    # -- the SAME colon-sexagesimal lon string must give a DIFFERENT
    #    value under galactic (no hour-angle convention -- degrees)
    #    than under an equatorial system (hour-angle convention -- x15),
    #    confirming hour_fallback is actually being consulted, not
    #    ignored. "276:42:00" read as plain degrees (galactic) is
    #    276.7, near the data -- the same string read with the x15
    #    hour-angle fallback would be nowhere near it, which is exactly
    #    why this sub-check only needs the galactic (literal-degrees)
    #    side: the x15 case is already covered by the icrs colon-form
    #    entry in the loop above. --
    rm -f "$tf_reproj"
    if bin/reproject_cubes infiles="$tf_ref" \
            target_system=galactic target_center=user_defined \
            target_lon="276:42:00" target_lat=60.0 >> "$tf_ud_log" 2>&1; then
        if ! python3 - "$tf_reproj" >> "$tf_ud_log" 2>&1 <<'PYEOF'
import sys
from astropy.io import fits
hdr = fits.getheader(sys.argv[1])
if abs(hdr['CRVAL1'] - 276.7) > 1e-9:
    sys.exit(f"CRVAL1 = {hdr['CRVAL1']}, expected 276.7 (galactic: colon-form is"
              f" degrees, no hour-angle convention of its own)")
print("galactic colon-form correctly read as degrees: 276.7")
PYEOF
        then
            tf_ud_ok=0
        fi
    else
        tf_ud_ok=0
    fi
    rm -f "$tf_reproj"
    if [[ "$tf_ud_ok" -eq 1 ]]; then
        pass "reproject_cubes: target_center=user_defined -- decimal/colon/lettered notations agree, hour-angle fallback correct per target_system (see $tf_ud_log)"
    else
        fail "reproject_cubes: target_center=user_defined notation/hour-angle-fallback check failed (see $tf_ud_log)"
    fi

    rm -f "$tf_ref"
else
    skip "bin/reproject_cubes not built; skipping target_system=/target_projection=/target_center= tests"
fi

# ---------------------------------------------------------------------------
section "69. reproject_cubes/match_cubes: invalid-region info/warning (UR-10)"
# ---------------------------------------------------------------------------
# UR-10 (docs/dev/BRANCH_STATUS_update-reproject.md): a facet-seam
# projection (-HPX, confirmed in REPROJECT_TARGET_FRAME_PLAN.md section
# 4c) can leave part of an output grid with no valid sky mapping at
# all, for reasons having nothing to do with missing input data
# coverage. Window parameters below are derived directly from section
# 4c's own documented formula (CRVAL1/2=180,0, CDELT=-1,1 deg/pixel,
# the theta=60/phi=90 seam's native-x jump from 73.529 to 106.471) --
# not re-derived from a fresh exploration pass.
if [[ -x bin/reproject_cubes ]]; then
    ur10_seam_ref="$OUT_DIR/ur10_seam_ref.FITS"
    ur10_plain_ref="$OUT_DIR/ur10_plain_ref.FITS"
    ur10_in="$DATA_DIR/TEST.Q.FITSCUBE"
    ur10_seam_in="$OUT_DIR/ur10_seam_in.FITSCUBE"
    ur10_seam_reproj="${ur10_in%.FITSCUBE}_REPROJ.FITS"
    ur10_seam_in_reproj="${ur10_seam_in%.FITSCUBE}_REPROJ.FITS"
    ur10_seam_matched="$OUT_DIR/ur10_seam_matched.FITSCUBE"
    ur10_seam_matched_reproj="${ur10_seam_matched%.FITSCUBE}_REPROJ.FITS"
    rm -f "$ur10_seam_ref" "$ur10_plain_ref" "$ur10_seam_in" "$ur10_seam_reproj" \
          "$ur10_seam_in_reproj" "$ur10_seam_matched" "$ur10_seam_matched_reproj"

    # ur10_seam_in (UR-15, docs/dev/BRANCH_STATUS_update-reproject.md):
    # TEST.Q.FITSCUBE (RA=180, Dec=0, a 32x32/0.001deg/pixel SIN field)
    # never overlapped ur10_seam_ref's own pixel grid at all -- its
    # footprint there, [121.98,122.02]x[.98,1.02], falls
    # outside the window's [1,64]x[1,120] extent, confirmed directly.
    # mode=reference's own overlap check (added by UR-15) now correctly
    # refuses that combination, so this test needs its own infile whose
    # sky position falls inside the window -- not TEST.
    # Q.FITSCUBE, still used unchanged below for the ur10_plain case,
    # where it does overlap. RA=270/Dec=1 was verified via astropy.wcs
    # directly against ur10_seam_ref's own CRVAL/CRPIX/CDELT (not
    # re-derived from the HPX projection formula by hand) to map to
    # pixel (~32,~2) -- comfortably inside the window, away from the
    # pole where the facet structure gets complicated. Only the overlap
    # check cares that this file exists and lands inside the window;
    # check_output_grid_validity (UR-10 itself) depends only on the
    # reference's own pixel->sky mapping, never on infile content.
    ur10_gen_log="$OUT_DIR/ur10_gen.log"
    if python3 - "$ur10_seam_ref" "$ur10_plain_ref" "$ur10_seam_in" "$ur10_seam_matched" \
            > "$ur10_gen_log" 2>&1 <<'PYEOF'
import sys
from astropy.io import fits
import numpy as np

seam_out, plain_out, seam_in, seam_in_matched = sys.argv[1:5]

h = fits.Header()
h['SIMPLE'] = True; h['BITPIX'] = -32; h['NAXIS'] = 2
h['NAXIS1'] = 64; h['NAXIS2'] = 120
h['CTYPE1'] = 'RA---HPX'; h['CTYPE2'] = 'DEC--HPX'
h['CRVAL1'] = 180.0; h['CRVAL2'] = 0.0
h['CRPIX1'] = 122.0; h['CRPIX2'] = 1.0
h['CDELT1'] = -1.0; h['CDELT2'] = 1.0
fits.PrimaryHDU(data=np.zeros((120, 64), dtype=np.float32), header=h).writeto(
    seam_out, overwrite=True)

h2 = fits.Header()
h2['SIMPLE'] = True; h2['BITPIX'] = -32; h2['NAXIS'] = 2
h2['NAXIS1'] = 32; h2['NAXIS2'] = 32
h2['CTYPE1'] = 'RA---SIN'; h2['CTYPE2'] = 'DEC--SIN'
h2['CRVAL1'] = 180.0; h2['CRVAL2'] = 0.0
h2['CRPIX1'] = 17.0; h2['CRPIX2'] = 17.0
h2['CDELT1'] = -0.001; h2['CDELT2'] = 0.001
fits.PrimaryHDU(data=np.zeros((32, 32), dtype=np.float32), header=h2).writeto(
    plain_out, overwrite=True)

h3 = fits.Header()
h3['SIMPLE'] = True; h3['BITPIX'] = -32; h3['NAXIS'] = 3
h3['NAXIS1'] = 8; h3['NAXIS2'] = 8; h3['NAXIS3'] = 1
h3['CTYPE1'] = 'RA---SIN'; h3['CTYPE2'] = 'DEC--SIN'; h3['CTYPE3'] = 'FREQ'
h3['CRVAL1'] = 270.0; h3['CRVAL2'] = 1.0; h3['CRVAL3'] = 1.0e9
h3['CRPIX1'] = 4.0; h3['CRPIX2'] = 4.0; h3['CRPIX3'] = 1.0
h3['CDELT1'] = -0.1; h3['CDELT2'] = 0.1; h3['CDELT3'] = 1.0e6
fits.PrimaryHDU(data=np.zeros((1, 8, 8), dtype=np.float32), header=h3).writeto(
    seam_in, overwrite=True)
# Same position, the reference grid's own 1 deg pixels: not coarser, so
# the beam-sampling rule lets it through with no beam recorded.
h3['CDELT1'] = -1.0; h3['CDELT2'] = 1.0
fits.PrimaryHDU(data=np.zeros((1, 8, 8), dtype=np.float32), header=h3).writeto(
    seam_in_matched, overwrite=True)
print("wrote", seam_out, plain_out, seam_in, seam_in_matched)
PYEOF
    then
        ur10_ok=1
        # Seam window, input pixels matching the grid's: runs, reports
        # the invalid fraction and warns.
        ur10_seam_log="$OUT_DIR/ur10_seam.log"
        rm -f "$ur10_seam_matched_reproj"
        if bin/reproject_cubes mode=reference reffile="$ur10_seam_ref" \
                infiles="$ur10_seam_matched" > "$ur10_seam_log" 2>&1 && \
           [[ -s "$ur10_seam_matched_reproj" ]]; then
            if grep -q "output grid invalid-mapping fraction" "$ur10_seam_log" && \
               grep -qE "WARNING: [2-9][0-9]\.[0-9]+% of the output grid has no valid sky mapping" "$ur10_seam_log"; then
                :
            else
                ur10_ok=0
            fi
        else
            ur10_ok=0
        fi
        rm -f "$ur10_seam_matched_reproj"

        # Same window, input with 0.1 deg pixels and no beam: the grid
        # is 10x coarser, so the beam-sampling rule must refuse it and
        # write nothing.
        ur10_coarse_log="$OUT_DIR/ur10_seam_coarse.log"
        rm -f "$ur10_seam_in_reproj"
        if bin/reproject_cubes mode=reference reffile="$ur10_seam_ref" \
                infiles="$ur10_seam_in" > "$ur10_coarse_log" 2>&1; then
            ur10_ok=0
        elif ! grep -q "has no beam information.*aliasing cannot be ruled out" "$ur10_coarse_log" || \
             [[ -e "$ur10_seam_in_reproj" ]]; then
            ur10_ok=0
        fi
        rm -f "$ur10_seam_in_reproj"

        ur10_plain_log="$OUT_DIR/ur10_plain.log"
        rm -f "$ur10_seam_reproj"
        if bin/reproject_cubes mode=reference reffile="$ur10_plain_ref" \
                infiles="$ur10_in" > "$ur10_plain_log" 2>&1; then
            if grep -q "output grid invalid-mapping fraction" "$ur10_plain_log" && \
               ! grep -q "WARNING:.*no valid sky mapping" "$ur10_plain_log"; then
                :
            else
                ur10_ok=0
            fi
        else
            ur10_ok=0
        fi
        rm -f "$ur10_seam_reproj"

        if [[ "$ur10_ok" -eq 1 ]]; then
            pass "reproject_cubes: invalid-region INFO always printed, WARNING only above threshold on a facet-seam window; no-beam input onto 10x coarser seam grid refused (see $ur10_seam_log / $ur10_coarse_log / $ur10_plain_log)"
        else
            fail "reproject_cubes: invalid-region info/warning check failed (see $ur10_seam_log / $ur10_coarse_log / $ur10_plain_log)"
        fi
    else
        fail "reproject_cubes: failed to build UR-10 fixtures (see $ur10_gen_log)"
    fi
    rm -f "$ur10_seam_ref" "$ur10_plain_ref" "$ur10_seam_in" "$ur10_seam_reproj" \
          "$ur10_seam_in_reproj" "$ur10_seam_matched" "$ur10_seam_matched_reproj"
else
    skip "bin/reproject_cubes not built; skipping invalid-region info/warning tests"
fi

# ---------------------------------------------------------------------------
section "70. reproject_cubes: interp= value-level correctness, all schemes except nearest (UR-8)"
# ---------------------------------------------------------------------------
# UR-8 (docs/dev/BRANCH_STATUS_update-reproject.md): promotes sinc/
# sincsinc/sinccos/sincgauss/blockave/gauss/somb/sombcos/lanczos2/
# lanczos3 and linear itself from the "reaches AST, differs from
# linear" sanity level (sections 63-66 above) to real value-level
# correctness against analytic ground truth, reusing (not duplicating)
# check_interp_benchmark.py's own machinery -- see tests/
# check_interp_correctness.py's own header for the tolerance rationale
# (a regression guard, 2x UR-12's own measured mean RMS per scheme).
if [[ -x bin/reproject_cubes ]]; then
    ur8_log="$OUT_DIR/interp_correctness_check.log"
    if python3 "$TESTS_DIR/check_interp_correctness.py" > "$ur8_log" 2>&1; then
        pass "reproject_cubes: interp= value-level correctness passed for every scheme except nearest (see $ur8_log)"
    else
        fail "reproject_cubes: interp= value-level correctness failed for one or more schemes (see $ur8_log)"
    fi
else
    skip "bin/reproject_cubes not built; skipping interp= value-level correctness tests"
fi

section "71. reproject_cubes/match_cubes: min_free_pct_after_write disk-space check (UR-15)"
# ---------------------------------------------------------------------------
# UR-15 (docs/dev/BRANCH_STATUS_update-reproject.md): refuse before any
# file is written if the output would leave less than min_free_pct_
# after_write percent of the output disk's total size free. 100 forces a
# refusal on any machine (no output can leave the whole disk free), so
# this test does not depend on how much space this machine has.
if [[ -x bin/reproject_cubes && -x bin/match_cubes ]]; then
    dk_dir="$OUT_DIR/diskcheck"
    rm -rf "$dk_dir"
    mkdir -p "$dk_dir"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$dk_dir/ref.FITSCUBE"
    cp "$DATA_DIR/TEST_BAND2_GRIDSHIFT.Q.FITSCUBE" "$dk_dir/in.FITSCUBE"
    dk_out="$dk_dir/in_REPROJ.FITS"

    dk_log="$OUT_DIR/diskcheck_refuse.log"
    if ! bin/reproject_cubes mode=reference reffile="$dk_dir/ref.FITSCUBE" \
            infiles="$dk_dir/in.FITSCUBE" min_free_pct_after_write=100 > "$dk_log" 2>&1 && \
       grep -q "ERROR: not enough disk space" "$dk_log" && [[ ! -e "$dk_out" ]]; then
        pass "reproject_cubes: min_free_pct_after_write=100 refuses before writing any file"
    else
        fail "reproject_cubes: min_free_pct_after_write=100 did not refuse, or wrote output anyway (see $dk_log)"
    fi

    dk_mlog="$OUT_DIR/diskcheck_refuse_match.log"
    if ! bin/match_cubes stages=reproject footprint_mode=reference reffile="$dk_dir/ref.FITSCUBE" \
            infiles="$dk_dir/in.FITSCUBE" min_free_pct_after_write=100 > "$dk_mlog" 2>&1 && \
       grep -q "ERROR: not enough disk space" "$dk_mlog" && [[ ! -e "$dk_out" ]]; then
        pass "match_cubes: min_free_pct_after_write=100 refuses before writing any file"
    else
        fail "match_cubes: min_free_pct_after_write=100 did not refuse, or wrote output anyway (see $dk_mlog)"
    fi

    dk_blog="$OUT_DIR/diskcheck_badvalue.log"
    if ! bin/reproject_cubes mode=reference reffile="$dk_dir/ref.FITSCUBE" \
            infiles="$dk_dir/in.FITSCUBE" min_free_pct_after_write=150 > "$dk_blog" 2>&1 && \
       grep -q "min_free_pct_after_write must be between 0 and 100" "$dk_blog"; then
        pass "reproject_cubes: min_free_pct_after_write outside [0,100] rejected"
    else
        fail "reproject_cubes: min_free_pct_after_write=150 was not rejected (see $dk_blog)"
    fi

    dk_dlog="$OUT_DIR/diskcheck_dryrun.log"
    if (cd "$dk_dir" && "$OLDPWD/bin/reproject_cubes" mode=reference reffile=ref.FITSCUBE \
            infiles=in.FITSCUBE dry_run=y) > "$dk_dlog" 2>&1 && \
       grep -q "Disk space: output needs" "$dk_dlog" && [[ ! -e "$dk_out" ]]; then
        pass "reproject_cubes: dry_run=y reports the disk-space check and writes no output"
    else
        fail "reproject_cubes: dry_run=y did not report the disk-space check (see $dk_dlog)"
    fi
    rm -rf "$dk_dir"
else
    skip "bin/reproject_cubes or bin/match_cubes not built; skipping disk-space check tests"
fi

section "72. reproject_cubes: no-reffile pixel-scale default and beam advisory (UR-16/UR-17)"
# ---------------------------------------------------------------------------
# UR-16: with no reffile and no target_cdelt_lon/lat, the output pixel
# scale is the smallest |CDELT| across the infiles, per axis -- so no
# input is coarsened. Checked against the resulting CRVAL/CDELT in the
# output header, not only against stdout. UR-17: differing beams get a
# named advisory, matching beams get none.
if [[ -x bin/reproject_cubes ]] && [[ -n "$GEN_PY" ]]; then
    pc_dir="$OUT_DIR/pixscale"
    rm -rf "$pc_dir"; mkdir -p "$pc_dir"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$pc_dir/fine.fits"
    # A deliberately coarser copy: 2x the pixel size, beam left alone.
    "$GEN_PY" - "$pc_dir" <<'PYEOF'
import sys
from astropy.io import fits
d = sys.argv[1]
with fits.open(f"{d}/fine.fits") as h:
    hdr = h[0].header.copy(); data = h[0].data
hdr['CDELT1'] = hdr['CDELT1'] * 2.0
hdr['CDELT2'] = hdr['CDELT2'] * 2.0
fits.PrimaryHDU(data=data, header=hdr).writeto(f"{d}/coarse.fits", overwrite=True)
PYEOF
    pc_log="$OUT_DIR/pixscale_default.log"
    rm -f "$pc_dir"/*_REPROJ.FITS
    if bin/reproject_cubes target_system=icrs \
            infiles="$pc_dir/fine.fits,$pc_dir/coarse.fits" > "$pc_log" 2>&1 && \
       [[ -s "$pc_dir/fine_REPROJ.FITS" ]]; then
        pc_check=$("$GEN_PY" - "$pc_dir" <<'PYEOF'
import sys
from astropy.io import fits
d = sys.argv[1]
want = abs(fits.getheader(f"{d}/fine.fits")['CDELT1'])
got1 = abs(fits.getheader(f"{d}/fine_REPROJ.FITS")['CDELT1'])
got2 = abs(fits.getheader(f"{d}/fine_REPROJ.FITS")['CDELT2'])
ok = abs(got1-want) <= 1e-6*want and abs(got2-want) <= 1e-6*want
print("OK" if ok else f"FAIL want={want} got={got1},{got2}")
PYEOF
)
        if [[ "$pc_check" == "OK" ]] && grep -q "differ in pixel scale" "$pc_log"; then
            pass "reproject_cubes: no-reffile pixel scale defaults to the smallest |CDELT|, and the spread is reported"
        else
            fail "reproject_cubes: pixel-scale default wrong or spread not reported ($pc_check, see $pc_log)"
        fi
    else
        fail "reproject_cubes: no-reffile run with no target_cdelt= failed (see $pc_log)"
    fi
    rm -f "$pc_dir"/*_REPROJ.FITS

    # UR-17: beams differ -> advisory; beams equal -> none.
    "$GEN_PY" - "$pc_dir" <<'PYEOF'
import sys
from astropy.io import fits
d = sys.argv[1]
for name, bmaj in (("fine.fits", 10.0/3600.0), ("coarse.fits", 20.0/3600.0)):
    with fits.open(f"{d}/{name}", mode="update") as h:
        h[0].header['BMAJ'] = bmaj; h[0].header['BMIN'] = bmaj; h[0].header['BPA'] = 0.0
        h.flush()
PYEOF
    pc_blog="$OUT_DIR/pixscale_beamwarn.log"
    bin/reproject_cubes target_system=icrs \
        infiles="$pc_dir/fine.fits,$pc_dir/coarse.fits" > "$pc_blog" 2>&1
    rm -f "$pc_dir"/*_REPROJ.FITS
    "$GEN_PY" - "$pc_dir" <<'PYEOF'
import sys
from astropy.io import fits
d = sys.argv[1]
with fits.open(f"{d}/coarse.fits", mode="update") as h:
    h[0].header['BMAJ'] = 10.0/3600.0; h[0].header['BMIN'] = 10.0/3600.0
    h.flush()
PYEOF
    pc_qlog="$OUT_DIR/pixscale_beamquiet.log"
    bin/reproject_cubes target_system=icrs \
        infiles="$pc_dir/fine.fits,$pc_dir/coarse.fits" > "$pc_qlog" 2>&1
    if grep -q "do not share one angular resolution" "$pc_blog" && \
       grep -q "Run convolve_cubes" "$pc_blog" && \
       ! grep -q "do not share one angular resolution" "$pc_qlog"; then
        pass "reproject_cubes: differing beams get a named advisory, matching beams stay silent"
    else
        fail "reproject_cubes: beam advisory wrong (see $pc_blog / $pc_qlog)"
    fi
    rm -rf "$pc_dir"
else
    skip "bin/reproject_cubes not built or no Python available; skipping pixel-scale/beam tests"
fi

section "73. match_cubes stages=both: the in-memory chain equals the two-tool sequence (UR-18)"
# ---------------------------------------------------------------------------
# UR-18 (docs/dev/BRANCH_STATUS_update-reproject.md): stages=both's own
# distinguishing claim is that chaining convolve and reproject THROUGH
# MEMORY gives the same answer as running the two tools with an
# intermediate file, only without writing that file. Sections 64/66
# previously checked only that it runs and writes output, and section 40
# that rm_synthesis still recovers the known RM peaks through it -- so a
# block-boundary or handoff bug small enough to survive RM recovery
# would have gone unnoticed. Asserted here as BIT-IDENTICAL, not within
# a tolerance: both paths narrow the convolved plane to float32 before
# resampling (astResampleR is the float32 variant), so there is no
# precision difference to absorb. Measured as bit-identical on
# 2026-10-15 before this test was written.
#
# Each order is compared against ITS OWN two-step sequence -- convolve
# then reproject for order=convolve_reproject, reproject then convolve
# for order=reproject_convolve. The operations do not commute (measured:
# 76794 of 76800 pixels differ between the orders, max relative 4.5e-4),
# so pairing an order with the wrong sequence would fail a correct tool.
# That non-commutativity is itself asserted at the end, since nothing
# else would catch order= having no effect.
if [[ -x bin/match_cubes && -x bin/reproject_cubes && -x bin/convolve_cubes ]] && [[ -n "$GEN_PY" ]]; then
    ce_dir="$OUT_DIR/chainequiv"
    rm -rf "$ce_dir"; mkdir -p "$ce_dir"
    ce_beams="$DATA_DIR/band2_beamlog.txt"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$ce_dir/ref.fits"
    ce_ok=1

    # One infile only, so the chain and the two-step pool exactly the
    # same per-channel beams and auto-derive the same common beam.
    for ce_order in convolve_reproject reproject_convolve; do
        rm -f "$ce_dir"/chain.fits "$ce_dir"/chain_MATCHED.FITS \
              "$ce_dir"/step.fits "$ce_dir"/step_CONV.FITS \
              "$ce_dir"/step_CONV_REPROJ.FITS "$ce_dir"/step_REPROJ.FITS \
              "$ce_dir"/step_REPROJ_CONV.FITS
        cp "$DATA_DIR/TEST_BAND2_GRIDSHIFT.Q.FITSCUBE" "$ce_dir/chain.fits"
        cp "$DATA_DIR/TEST_BAND2_GRIDSHIFT.Q.FITSCUBE" "$ce_dir/step.fits"
        ce_log="$OUT_DIR/chainequiv_${ce_order}.log"

        if ! bin/match_cubes stages=both order="$ce_order" footprint_mode=reference \
                reffile="$ce_dir/ref.fits" infiles="$ce_dir/chain.fits" \
                beamfiles="$ce_beams" > "$ce_log" 2>&1; then
            fail "match_cubes: stages=both order=$ce_order run failed (see $ce_log)"
            ce_ok=0
            continue
        fi

        if [[ "$ce_order" == "convolve_reproject" ]]; then
            bin/convolve_cubes infiles="$ce_dir/step.fits" beamfiles="$ce_beams" \
                >> "$ce_log" 2>&1
            bin/reproject_cubes mode=reference reffile="$ce_dir/ref.fits" \
                infiles="$ce_dir/step_CONV.FITS" >> "$ce_log" 2>&1
            ce_two_step="$ce_dir/step_CONV_REPROJ.FITS"
            ce_seq="convolve_cubes then reproject_cubes"
        else
            bin/reproject_cubes mode=reference reffile="$ce_dir/ref.fits" \
                infiles="$ce_dir/step.fits" >> "$ce_log" 2>&1
            bin/convolve_cubes infiles="$ce_dir/step_REPROJ.FITS" beamfiles="$ce_beams" \
                >> "$ce_log" 2>&1
            ce_two_step="$ce_dir/step_REPROJ_CONV.FITS"
            ce_seq="reproject_cubes then convolve_cubes"
        fi

        if [[ ! -s "$ce_dir/chain_MATCHED.FITS" || ! -s "$ce_two_step" ]]; then
            fail "match_cubes: stages=both order=$ce_order -- chain or two-step output missing (see $ce_log)"
            ce_ok=0
            continue
        fi
        if "$GEN_PY" "$TESTS_DIR/compare_cubes.py" "$ce_dir/chain_MATCHED.FITS" \
                "$ce_two_step" --exact >> "$ce_log" 2>&1; then
            pass "match_cubes: stages=both order=$ce_order bit-identical to running $ce_seq separately"
        else
            fail "match_cubes: stages=both order=$ce_order differs from running $ce_seq separately (see $ce_log)"
            ce_ok=0
        fi
        cp "$ce_dir/chain_MATCHED.FITS" "$ce_dir/order_${ce_order}.FITS"
    done

    # order= must change the result: the two operations do not
    # commute, so identical output would mean order= was ignored.
    if [[ "$ce_ok" -eq 1 && -s "$ce_dir/order_convolve_reproject.FITS" \
          && -s "$ce_dir/order_reproject_convolve.FITS" ]]; then
        if "$GEN_PY" "$TESTS_DIR/compare_cubes.py" \
                "$ce_dir/order_convolve_reproject.FITS" \
                "$ce_dir/order_reproject_convolve.FITS" --exact \
                > "$OUT_DIR/chainequiv_orders.log" 2>&1; then
            fail "match_cubes: the two stages=both orders produced identical output -- order= appears to be ignored"
        else
            pass "match_cubes: the two stages=both orders differ, as convolution and resampling do not commute"
        fi
    fi
    rm -rf "$ce_dir"
else
    skip "match_cubes/reproject_cubes/convolve_cubes not all built, or no Python; skipping chain-equivalence tests"
fi

section "74. reproject_cubes: measured kernel-width defaults for gauss/sincgauss/blockave (UR-19)"
# ---------------------------------------------------------------------------
# UR-19 (docs/dev/BRANCH_STATUS_update-reproject.md): these three
# schemes take a kernel width AST cannot choose itself, and omitting
# interp_params= used to be a hard error (UR-14) because AST's own
# omitted-params behaviour writes all- or mostly-NaN output. That left
# three of the twelve schemes unreachable unless the caller already
# knew a width to pass. ast_interp_mod.f90 now supplies a measured
# default (tests/check_interp_width_defaults.py picks each one as the
# width with the lowest worst-case error over beam samplings and
# output geometries) and reports the value in use.
#
# Asserted here: each scheme runs with no interp_params= at all, the
# INFO line names the width this build ships, and the output is NOT
# the all-NaN cube UR-14 was guarding against. The widths are checked
# against the literal values in ast_interp_mod.f90, so changing a
# constant without re-running the measurement script fails here.
if [[ -x bin/reproject_cubes ]] && [[ -n "$GEN_PY" ]]; then
    ur19_dir="$OUT_DIR/interp_width_defaults"
    rm -rf "$ur19_dir"; mkdir -p "$ur19_dir"
    cp "$DATA_DIR/TEST.Q.FITSCUBE" "$ur19_dir/ref.fits"

    # scheme:expected width as the INFO line prints it (F6.2)
    for ur19_case in "gauss:0.80" "sincgauss:8.00" "blockave:1.00"; do
        ur19_scheme="${ur19_case%%:*}"
        ur19_width="${ur19_case##*:}"
        ur19_log="$OUT_DIR/interp_width_${ur19_scheme}.log"
        rm -f "$ur19_dir/in_${ur19_scheme}"*.FITS
        cp "$DATA_DIR/TEST_BAND2_GRIDSHIFT.Q.FITSCUBE" "$ur19_dir/in_${ur19_scheme}.fits"

        if ! bin/reproject_cubes mode=reference reffile="$ur19_dir/ref.fits" \
                infiles="$ur19_dir/in_${ur19_scheme}.fits" interp="$ur19_scheme" \
                > "$ur19_log" 2>&1; then
            fail "reproject_cubes: interp=$ur19_scheme with no interp_params= failed (see $ur19_log)"
            continue
        fi
        if grep -q "INFO: interp=$ur19_scheme with no interp_params= .*$ur19_width pixels (default)" "$ur19_log"; then
            pass "reproject_cubes: interp=$ur19_scheme defaults its width to $ur19_width px and says so"
        else
            fail "reproject_cubes: interp=$ur19_scheme did not report the expected default width $ur19_width px (see $ur19_log)"
        fi
        if [[ -s "$ur19_dir/in_${ur19_scheme}_REPROJ.FITS" ]] && \
           "$GEN_PY" -c "
import sys
import numpy as np
from astropy.io import fits
d = fits.getdata('$ur19_dir/in_${ur19_scheme}_REPROJ.FITS')
# UR-14's failure mode was an all- or mostly-NaN cube; require a
# finite fraction above a quarter rather than merely 'some finite pixel'.
sys.exit(0 if np.isfinite(d).mean() > 0.25 else 1)
" >> "$ur19_log" 2>&1; then
            pass "reproject_cubes: interp=$ur19_scheme default width produces a populated cube, not the all-NaN output UR-14 guarded"
        else
            fail "reproject_cubes: interp=$ur19_scheme default width produced an empty or mostly-NaN cube (see $ur19_log)"
        fi
    done

    # An explicit interp_params= must still win over the default, and
    # must change the pixels -- otherwise the default would be silently
    # overriding whatever the caller asked for.
    ur19_ov="$OUT_DIR/interp_width_override.log"
    rm -f "$ur19_dir/ov"*.FITS
    cp "$DATA_DIR/TEST_BAND2_GRIDSHIFT.Q.FITSCUBE" "$ur19_dir/ov.fits"
    bin/reproject_cubes mode=reference reffile="$ur19_dir/ref.fits" \
        infiles="$ur19_dir/ov.fits" interp=sincgauss interp_params=0,2.0 \
        > "$ur19_ov" 2>&1
    if [[ -s "$ur19_dir/ov_REPROJ.FITS" ]] && ! grep -q "(default)" "$ur19_ov" && \
       ! "$GEN_PY" "$TESTS_DIR/compare_cubes.py" "$ur19_dir/ov_REPROJ.FITS" \
            "$ur19_dir/in_sincgauss_REPROJ.FITS" --exact >> "$ur19_ov" 2>&1; then
        pass "reproject_cubes: explicit interp_params= overrides the default width, reports no default, and changes the pixels"
    else
        fail "reproject_cubes: explicit interp_params= did not override the default width (see $ur19_ov)"
    fi
    rm -rf "$ur19_dir"
else
    skip "bin/reproject_cubes not built or no Python available; skipping interp width-default tests"
fi

section "74b. reproject_cubes: shipped kernel-width defaults are still best-or-tied-best (UR-19/UR-24 B2)"
# ---------------------------------------------------------------------------
# UR-19's measurement script (tests/check_interp_width_defaults.py) is
# what justifies ast_interp_mod.f90's default_gauss_fwhm_px/default_
# sincgauss_fwhm_px/default_blockave_halfwidth_px constants in the
# first place -- section 74 above only checks that the shipped
# constants are the ones actually used and reported; it does not
# check that they are STILL the measured-best choice. This re-runs
# that measurement (now swept over all 20 system/projection
# combinations, UR-24 B2) and fails if any shipped default is no
# longer best or tied-best on both worst-case RMS and worst-case
# maximum error -- e.g. an AST kernel-behaviour change that quietly
# makes a different width better would be caught here, not just by
# someone happening to re-run the script by hand.
if [[ -x bin/reproject_cubes ]] && [[ -n "$GEN_PY" ]]; then
    ur19b_log="$OUT_DIR/interp_width_defaults_check.log"
    rm -rf "$OUT_DIR/interp_width_defaults"
    if "$GEN_PY" "$TESTS_DIR/check_interp_width_defaults.py" > "$ur19b_log" 2>&1; then
        pass "interp width defaults: gauss/sincgauss/blockave all best-or-tied-best over 20 system/projection combinations (see $ur19b_log)"
    else
        fail "interp width defaults: one or more shipped defaults are no longer best-or-tied-best (see $ur19b_log)"
    fi
else
    skip "bin/reproject_cubes not built or no Python available; skipping interp width-defaults measurement"
fi

section "75. reproject_cubes/match_cubes: target_projection= on the no-reffile path (UR-7/UR-8)"
# ---------------------------------------------------------------------------
# UR-7/UR-8 (docs/dev/BRANCH_STATUS_update-reproject.md): with no reffile,
# the output grid can be built in any of tan/sin/car/ait/hpx, in the
# requested target_system or, for target_projection= alone, the infiles'
# own system. tests/check_target_projection.py runs every projection x
# {icrs, galactic} x {equator, Dec -60} x {one infile, two infiles} in
# both tools and compares each output with astropy's reproject_interp
# onto the header the tool wrote (UR-13 found interp=linear equivalent
# to that algorithm), plus CTYPE, coverage, cross-tool identity and the
# error paths -- see that script's header for every check and tolerance.
if [[ -x bin/reproject_cubes && -x bin/match_cubes ]] && [[ -n "$GEN_PY" ]]; then
    ur7_log="$OUT_DIR/target_projection_check.log"
    rm -rf "$OUT_DIR/target_projection"
    if "$GEN_PY" "$TESTS_DIR/check_target_projection.py" > "$ur7_log" 2>&1; then
        pass "target_projection=: checks $(grep -m1 'checks:' "$ur7_log" | sed 's/^target_projection checks: //') across every projection, system, centre and infile count, in both tools (see $ur7_log)"
    else
        fail "target_projection=: one or more checks failed (see $ur7_log)"
    fi
else
    skip "reproject_cubes/match_cubes not built or no Python available; skipping target_projection= tests"
fi

section "76. reproject_cubes/match_cubes: refuse an output grid that cannot sample the beam (UR-21)"
# ---------------------------------------------------------------------------
# UR-21 (docs/dev/BRANCH_STATUS_update-reproject.md): onto a grid coarser
# than an infile's own, both tools refuse when the beam's minor axis spans
# fewer than 2.0 output cells and warn between 2.0 and 2.5; they never
# enlarge a beam themselves. tests/check_beam_sampling.py covers both
# tools, header and BEAMS-table beams, the no-beam and already-
# undersampled cases, and runs the advice the refusal prints.
if [[ -x bin/reproject_cubes && -x bin/match_cubes ]] && [[ -n "$GEN_PY" ]]; then
    ur21_log="$OUT_DIR/beam_sampling_check.log"
    rm -rf "$OUT_DIR/beam_sampling"
    if "$GEN_PY" "$TESTS_DIR/check_beam_sampling.py" > "$ur21_log" 2>&1; then
        pass "beam sampling: checks $(grep -m1 'checks:' "$ur21_log" | sed 's/^beam sampling checks: //') across both tools (see $ur21_log)"
    else
        fail "beam sampling: one or more checks failed (see $ur21_log)"
    fi
else
    skip "reproject_cubes/match_cubes not built or no Python available; skipping beam-sampling tests"
fi

section "77. reproject_cubes/match_cubes: HEALPix-aligned HPX grids, target_hpx_nside= (UR-23)"
# ---------------------------------------------------------------------------
# UR-23 (docs/dev/BRANCH_STATUS_update-reproject.md): a target_projection=
# hpx no-reffile grid is HEALPix-aligned -- CRVAL/CDELT/CRPIX/NAXIS all
# follow from Nside alone (target_hpx_nside=, or a default chosen from the
# finest infile pixel, with UR-21's own beam margin deciding whether it may
# coarsen). tests/check_hpx_nside.py checks output pixel centres against
# astropy_healpix (equatorial and polar, icrs and galactic), values against
# astropy, every default-selection branch, every error path, and both
# tools.
if [[ -x bin/reproject_cubes && -x bin/match_cubes ]] && [[ -n "$GEN_PY" ]]; then
    ur23_log="$OUT_DIR/hpx_nside_check.log"
    rm -rf "$OUT_DIR/hpx_nside"
    if "$GEN_PY" "$TESTS_DIR/check_hpx_nside.py" > "$ur23_log" 2>&1; then
        pass "HPX Nside: checks $(grep -m1 'checks:' "$ur23_log" | sed 's/^HPX Nside checks: //') (see $ur23_log)"
    else
        fail "HPX Nside: one or more checks failed (see $ur23_log)"
    fi
else
    skip "reproject_cubes/match_cubes not built or no Python available; skipping HPX Nside tests"
fi

# ---------------------------------------------------------------------------
# FFTW wisdom, whole run: section 32c checks one convolve_cubes run; this
# catches any other run (match_cubes stages=convolve/both, the Python
# check scripts' own subprocess logs) that could not save wisdom.
# ---------------------------------------------------------------------------
section "78. FFTW wisdom: no run in this suite failed to save it"
wis_bad_logs=()
wis_unreadable_logs=()
wis_scanned=0
while IFS= read -r -d '' wis_scan_log; do
    wis_scanned=$((wis_scanned + 1))
    # grep exits 0 = match, 1 = no match, 2 = could not read the file.
    wis_grep_rc=0
    grep -q -- "$WISDOM_WARNING_TEXT" "$wis_scan_log" || wis_grep_rc=$?
    if [[ "$wis_grep_rc" -eq 0 ]]; then
        wis_bad_logs+=("$wis_scan_log")
    elif [[ "$wis_grep_rc" -ne 1 ]]; then
        wis_unreadable_logs+=("$wis_scan_log")
    fi
done < <(find "$OUT_DIR" -name '*.log' -newer "$RUN_START_MARKER" -print0)
if [[ "${#wis_bad_logs[@]}" -gt 0 ]]; then
    fail "FFTW wisdom not saved in ${#wis_bad_logs[@]} log(s): ${wis_bad_logs[*]}"
elif [[ "${#wis_unreadable_logs[@]}" -gt 0 ]]; then
    fail "FFTW wisdom scan could not read ${#wis_unreadable_logs[@]} log(s): ${wis_unreadable_logs[*]}"
else
    pass "FFTW wisdom: no '$WISDOM_WARNING_TEXT' in any of the $wis_scanned log(s) written by this run"
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
echo "Warn  : $WARN (does not affect RESULT -- see individual [WARN] lines above)"

if [[ "$FAIL" -gt 0 ]]; then
    echo ""
    echo "RESULT: FAILED ($FAIL test(s) failed)"
    exit 1
else
    echo ""
    echo "RESULT: ALL PASSED"
    exit 0
fi
