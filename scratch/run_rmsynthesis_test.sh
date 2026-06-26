#!/usr/bin/env bash
set -euo pipefail

# RM-synthesis RM-cube generation runner
# Requires: config file (mandatory), optional thread count

usage() {
  cat <<EOF
Usage: $(basename "$0") <config_file> [num_threads]

Positional Arguments:
  <config_file>    Config file name (relative to cfg/); required
  [num_threads]    Number of OMP threads (default: 6)

Examples:
  $(basename "$0") rmsynth-casa.fullim.cfg
  $(basename "$0") rmsynth-casa.fullim.cfg 8
  $(basename "$0") cfg/benchmark.cfg 4

Environment:
  OMP_PROC_BIND=close  (hardcoded to avoid hyperthreading)
  OMP_PLACES=cores     (hardcoded to avoid hyperthreading)

Binary:
  Uses: bin/rm_synthesis_release_omp1 (release + OpenMP)

EOF
}

if [[ $# -eq 0 ]]; then
  echo "ERROR: missing required argument <config_file>" >&2
  usage >&2
  exit 1
fi

if [[ "${1}" == "-h" || "${1}" == "--help" ]]; then
  usage
  exit 0
fi

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXE="${ROOT_DIR}/bin/rm_synthesis_release_omp1"
CFG_NAME="${1}"
OMP_NUM_THREADS="${2:-6}"
CFG_PATH="${ROOT_DIR}/cfg/${CFG_NAME}"
CFG_ARG="../cfg/${CFG_NAME}"

cd "${ROOT_DIR}/scratch"

echo "[runFile] Working directory: $PWD"
echo "[runFile] Executable: ${EXE}"
echo "[runFile] Config: ${CFG_PATH}"

if [[ ! -x "${EXE}" ]]; then
  echo "[runFile] ERROR: executable not found or not executable: ${EXE}" >&2
  echo "[runFile] Build first with: make all" >&2
  exit 1
fi

if [[ ! -f "${CFG_PATH}" ]]; then
  echo "[runFile] ERROR: config not found: ${CFG_PATH}" >&2
  exit 1
fi

DATA_PATH="$(awk -F= '/^path=/{print $2; exit}' "${CFG_PATH}" | xargs)"
Q_FILE="$(awk -F= '/^infileQ=/{print $2; exit}' "${CFG_PATH}" | xargs)"
U_FILE="$(awk -F= '/^infileU=/{print $2; exit}' "${CFG_PATH}" | xargs)"
OUT_BASE="$(awk -F= '/^outfile=/{print $2; exit}' "${CFG_PATH}" | xargs)"
AP_MODE="$(awk -F= '/^ap_angle_mode=/{print $2; exit}' "${CFG_PATH}" | xargs)"

if [[ -z "${DATA_PATH}" || -z "${Q_FILE}" || -z "${U_FILE}" || -z "${OUT_BASE}" ]]; then
  echo "[runFile] ERROR: cfg missing one of path/infileQ/infileU/outfile" >&2
  exit 1
fi

if [[ -z "${AP_MODE}" ]]; then
  AP_MODE="phase"
fi

if [[ ! -f "${DATA_PATH}${Q_FILE}" ]]; then
  echo "[runFile] ERROR: Q cube not found: ${DATA_PATH}${Q_FILE}" >&2
  exit 1
fi

if [[ ! -f "${DATA_PATH}${U_FILE}" ]]; then
  echo "[runFile] ERROR: U cube not found: ${DATA_PATH}${U_FILE}" >&2
  exit 1
fi

ANGLE_FILE="${OUT_BASE}.PHA.RMCUBE.FITS"
if [[ "${AP_MODE}" == "pol" ]]; then
  ANGLE_FILE="${OUT_BASE}.POLA.RMCUBE.FITS"
fi

rm -f "${OUT_BASE}.AMP.RMCUBE.FITS" "${ANGLE_FILE}"

# Run with OMP settings and timing
export OMP_NUM_THREADS="${OMP_NUM_THREADS}"
export OMP_PROC_BIND=close
export OMP_PLACES=cores

echo "[runFile] Running with OMP_NUM_THREADS=${OMP_NUM_THREADS}"
/usr/bin/time -v "${EXE}" "${CFG_ARG}"

if [[ ! -s "${OUT_BASE}.AMP.RMCUBE.FITS" || ! -s "${ANGLE_FILE}" ]]; then
  echo "[runFile] ERROR: run did not produce expected output cubes." >&2
  echo "[runFile] Expected: ${PWD}/${OUT_BASE}.AMP.RMCUBE.FITS and ${PWD}/${ANGLE_FILE}" >&2
  exit 2
fi

echo "[runFile] Test run finished."
