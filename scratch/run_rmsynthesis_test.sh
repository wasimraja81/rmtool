#!/usr/bin/env bash
set -euo pipefail

# Run RM-synthesis CASA Q/U test from scratch/
# This script uses cfg/test_my_casa_q_u.cfg

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXE="${ROOT_DIR}/bin/rm_synthesis"
CFG_NAME="test_my_casa_q_u.cfg"
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

# Run standard test
"${EXE}" "${CFG_ARG}"

if [[ ! -s "${OUT_BASE}.AMP.RMCUBE.FITS" || ! -s "${ANGLE_FILE}" ]]; then
  echo "[runFile] ERROR: run did not produce expected output cubes." >&2
  echo "[runFile] Expected: ${PWD}/${OUT_BASE}.AMP.RMCUBE.FITS and ${PWD}/${ANGLE_FILE}" >&2
  exit 2
fi

echo "[runFile] Test run finished."
