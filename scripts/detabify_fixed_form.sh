#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <file1.f> [file2.f ...]" >&2
  exit 2
fi

for f in "$@"; do
  if [[ ! -f "$f" ]]; then
    echo "Missing file: $f" >&2
    exit 2
  fi

  tmp="${f}.tmp.$$"
  # Fixed-form Fortran should be treated with tab stops at 8 columns.
  expand -t 8 "$f" > "$tmp"
  mv "$tmp" "$f"
  echo "Detabified: $f"
done
