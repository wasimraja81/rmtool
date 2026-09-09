#!/usr/bin/env bash
set -euo pipefail

# scripts/setup_venv.sh -- one-command Python environment setup for
# rmtool's own helper scripts (scripts/casda_fetch.py and friends) and
# test helpers (tests/*.py). Creates ~/venv/rmtool only if it doesn't
# already exist -- never recreates or touches an existing one.
# Dependencies always install via that venv's own pip binary directly
# (never system pip), so system Python is never touched.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VENV_DIR="${HOME}/venv/rmtool"

if [[ -d "${VENV_DIR}" ]]; then
  echo "[setup_venv] ${VENV_DIR} already exists -- leaving it as-is."
else
  echo "[setup_venv] Creating ${VENV_DIR}..."
  python3 -m venv "${VENV_DIR}"
fi

echo "[setup_venv] Installing dependencies into ${VENV_DIR} (system Python untouched)..."
"${VENV_DIR}/bin/pip" install --upgrade pip
"${VENV_DIR}/bin/pip" install -r "${ROOT_DIR}/requirements.txt"

echo "[setup_venv] Done. Use ${VENV_DIR}/bin/python3 to run rmtool's Python scripts, e.g.:"
echo "  ${VENV_DIR}/bin/python3 scripts/casda_fetch.py --help"
