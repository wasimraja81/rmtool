#!/usr/bin/env bash
set -euo pipefail

# RM-synthesis RM-cube generation runner
# Requires: config file (mandatory), optional thread count and backend

usage() {
  cat <<EOF
Usage: $(basename "$0") <config_file> [num_threads] [backend]

Positional Arguments:
  <config_file>    Config file name (relative to cfg/) or absolute path; required
  [num_threads]    Number of OMP threads (default: 6)
  [backend]        auto|cpu|gpu (default: auto)

Examples:
  $(basename "$0") rmsynth-casa.fullim.cfg
  $(basename "$0") rmsynth-casa.fullim.cfg 8
  $(basename "$0") cfg/benchmark.cfg 4 auto
  $(basename "$0") rmsynth-subim.cfg 1 gpu

Environment:
  CPU mode sets:
    OMP_PROC_BIND=close
    OMP_PLACES=cores
  GPU mode defaults (if unset):
    OMP_TARGET_OFFLOAD=MANDATORY
    OMP_DEFAULT_DEVICE=0

Binary:
  auto mode selects executable from cfg key use_gpu/use_gpus:
    use_gpu=y  -> bin/rm_synthesis_release_omp0_gpu1
    use_gpu=n  -> bin/rm_synthesis_release_omp1_gpu0

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
CFG_NAME="${1}"
OMP_NUM_THREADS="${2:-6}"
BACKEND="${3:-auto}"

if [[ "${BACKEND}" != "auto" && "${BACKEND}" != "cpu" && "${BACKEND}" != "gpu" ]]; then
  echo "[runFile] ERROR: backend must be one of auto|cpu|gpu" >&2
  usage >&2
  exit 1
fi

if [[ -f "${CFG_NAME}" ]]; then
  CFG_PATH="$(realpath "${CFG_NAME}")"
else
  CFG_PATH="${ROOT_DIR}/cfg/${CFG_NAME}"
fi

cd "${ROOT_DIR}/scratch"

if [[ "${CFG_PATH}" == "${ROOT_DIR}"/* ]]; then
  CFG_ARG="..${CFG_PATH#${ROOT_DIR}}"
else
  CFG_ARG="${CFG_PATH}"
fi

echo "[runFile] Config: ${CFG_PATH}"

if [[ ! -f "${CFG_PATH}" ]]; then
  echo "[runFile] ERROR: config not found: ${CFG_PATH}" >&2
  exit 1
fi

USE_GPU_CFG="$(awk -F= '
  /^use_gpu[[:space:]]*=/  {gsub(/[[:space:]]/,"",$2); print tolower($2); found=1; exit}
  /^use_gpus[[:space:]]*=/ {gsub(/[[:space:]]/,"",$2); print tolower($2); found=1; exit}
  END {if (!found) print "n"}
' "${CFG_PATH}")"

if [[ "${BACKEND}" == "cpu" ]]; then
  EXE="${ROOT_DIR}/bin/rm_synthesis_release_omp1_gpu0"
  USE_GPU_MODE="n"
elif [[ "${BACKEND}" == "gpu" ]]; then
  EXE="${ROOT_DIR}/bin/rm_synthesis_release_omp0_gpu1"
  USE_GPU_MODE="y"
else
  if [[ "${USE_GPU_CFG}" == "y" || "${USE_GPU_CFG}" == "yes" || "${USE_GPU_CFG}" == "1" ]]; then
    EXE="${ROOT_DIR}/bin/rm_synthesis_release_omp0_gpu1"
    USE_GPU_MODE="y"
  else
    EXE="${ROOT_DIR}/bin/rm_synthesis_release_omp1_gpu0"
    USE_GPU_MODE="n"
  fi
fi

echo "[runFile] Working directory: $PWD"
echo "[runFile] Backend: ${BACKEND} (cfg use_gpu=${USE_GPU_CFG}, effective=${USE_GPU_MODE})"
echo "[runFile] Executable: ${EXE}"

if [[ ! -x "${EXE}" ]]; then
  echo "[runFile] ERROR: executable not found or not executable: ${EXE}" >&2
  if [[ "${USE_GPU_MODE}" == "y" ]]; then
    echo "[runFile] Build first with: make GPU=1" >&2
  else
    echo "[runFile] Build first with: make GPU=0 OMP=1" >&2
  fi
  exit 1
fi

DATA_PATH="$(awk -F= '/^path=/{print $2; exit}' "${CFG_PATH}" | xargs)"
Q_FILE="$(awk -F= '/^infileQ=/{print $2; exit}' "${CFG_PATH}" | xargs)"
U_FILE="$(awk -F= '/^infileU=/{print $2; exit}' "${CFG_PATH}" | xargs)"
OUT_BASE="$(awk -F= '/^outfile=/{print $2; exit}' "${CFG_PATH}" | xargs)"
AP_MODE="$(awk -F= '/^ap_angle_mode=/{print $2; exit}' "${CFG_PATH}" | xargs)"
OUTPUT_MODE="$(awk -F= '/^output_mode=/{print $2; exit}' "${CFG_PATH}" | xargs)"

if [[ -z "${DATA_PATH}" || -z "${Q_FILE}" || -z "${U_FILE}" || -z "${OUT_BASE}" ]]; then
  echo "[runFile] ERROR: cfg missing one of path/infileQ/infileU/outfile" >&2
  exit 1
fi

if [[ -z "${AP_MODE}" ]]; then AP_MODE="phase"; fi
if [[ -z "${OUTPUT_MODE}" ]]; then OUTPUT_MODE="ap"; fi

if [[ ! -f "${DATA_PATH}${Q_FILE}" ]]; then
  echo "[runFile] ERROR: Q cube not found: ${DATA_PATH}${Q_FILE}" >&2
  exit 1
fi

if [[ ! -f "${DATA_PATH}${U_FILE}" ]]; then
  echo "[runFile] ERROR: U cube not found: ${DATA_PATH}${U_FILE}" >&2
  exit 1
fi

if [[ "${OUTPUT_MODE}" == "ri" ]]; then
  OUT_FILE_1="${OUT_BASE}.REAL.RMCUBE.FITS"
  OUT_FILE_2="${OUT_BASE}.IMAG.RMCUBE.FITS"
else
  OUT_FILE_1="${OUT_BASE}.AMP.RMCUBE.FITS"
  OUT_FILE_2="${OUT_BASE}.PHA.RMCUBE.FITS"
  if [[ "${AP_MODE}" == "pol" ]]; then
    OUT_FILE_2="${OUT_BASE}.POLA.RMCUBE.FITS"
  fi
fi

rm -f "${OUT_FILE_1}" "${OUT_FILE_2}"

# Run with OMP settings and timing
export OMP_NUM_THREADS="${OMP_NUM_THREADS}"
if [[ "${USE_GPU_MODE}" == "y" ]]; then
  export OMP_TARGET_OFFLOAD="${OMP_TARGET_OFFLOAD:-MANDATORY}"
  export OMP_DEFAULT_DEVICE="${OMP_DEFAULT_DEVICE:-0}"
else
  export OMP_PROC_BIND=close
  export OMP_PLACES=cores
fi

echo "[runFile] Running with OMP_NUM_THREADS=${OMP_NUM_THREADS}"
if [[ "${USE_GPU_MODE}" == "y" ]]; then
  echo "[runFile] OMP_TARGET_OFFLOAD=${OMP_TARGET_OFFLOAD} OMP_DEFAULT_DEVICE=${OMP_DEFAULT_DEVICE}"
fi
/usr/bin/time -v "${EXE}" "${CFG_ARG}"

DRY_RUN_CFG="$(awk -F= '/^dry_run[[:space:]]*=/{gsub(/[[:space:]]/,"",$2); print tolower($2); exit}' "${CFG_PATH}")"

if [[ "${DRY_RUN_CFG}" == "y" || "${DRY_RUN_CFG}" == "yes" || "${DRY_RUN_CFG}" == "1" ]]; then
  echo "[runFile] Dry-run mode: no output cubes expected."
  echo "[runFile] Dry-run finished. Check scratch/tile_autotune.cfg and scratch/runtime_estimate.txt"
else
  if [[ ! -s "${OUT_FILE_1}" || ! -s "${OUT_FILE_2}" ]]; then
    echo "[runFile] ERROR: run did not produce expected output cubes." >&2
    echo "[runFile] Expected: ${PWD}/${OUT_FILE_1} and ${PWD}/${OUT_FILE_2}" >&2
    exit 2
  fi
  echo "[runFile] Test run finished."
fi
