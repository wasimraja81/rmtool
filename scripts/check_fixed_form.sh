#!/usr/bin/env bash
set -euo pipefail

strict_block=0
if [[ ${1:-} == "--strict-block" ]]; then
  strict_block=1
  shift
fi

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 [--strict-block] <file1.f> [file2.f ...]" >&2
  exit 2
fi

status=0

for f in "$@"; do
  if [[ ! -f "$f" ]]; then
    echo "Missing file: $f" >&2
    status=2
    continue
  fi

  awk -v file="$f" -v strict_block="$strict_block" '
function ltrim(s) { sub(/^[ ]+/, "", s); return s }
BEGIN { depth = 0; err = 0 }
{
  line = $0

  if (index(line, "\t") > 0) {
    printf("%s:%d: tab character found\n", file, NR)
    err = 1
  }

  if (length(line) > 72) {
    printf("%s:%d: WARNING: line exceeds column 72 (%d chars)\n", file, NR, length(line))
    # not setting err=1: treated as warning only
  }

  # Blank line
  if (line ~ /^[ ]*$/) next

  c1 = substr(line, 1, 1)
  first = match(line, /[^ ]/)
  if (first == 0) next
  firstch = substr(line, first, 1)
  label = substr(line, 1, 5)

  # Fixed-form comment styles
  if (c1 == "c" || c1 == "C" || c1 == "*" || c1 == "!") next
  if (firstch == "!") next

  # Continuation style used in this project: '-' as first non-space at column 6.
  if (firstch == "-") {
    if (first != 6) {
      printf("%s:%d: continuation marker '-' must be in column 6 (found %d)\n", file, NR, first)
      err = 1
    }
    next
  }

  # Numeric statement labels are allowed in columns 1-5.
  if (label ~ /[0-9]/) {
    next
  }

  if (first < 7) {
    printf("%s:%d: statement starts before column 7 (found %d)\n", file, NR, first)
    err = 1
    next
  }

  if (strict_block == 1) {
    code = substr(line, 7)
    if (code ~ /^[ ]*$/) next

    indent = match(code, /[^ ]/) - 1
    stmt = ltrim(code)
    up = toupper(stmt)

    checkDepth = depth
    is_end = (up ~ /^(END[ ]*IF|ENDIF|END[ ]*DO|ENDDO|END[ ]*SELECT|ENDSELECT)\b/)
    is_else = (up ~ /^(ELSE\b|ELSE[ ]*IF\b|ELSEIF\b)/)

    if (is_end || is_else) {
      checkDepth = depth - 1
      if (checkDepth < 0) checkDepth = 0
    }

    expected = checkDepth * 2
    if (indent != expected) {
      printf("%s:%d: indent %d spaces after col7, expected %d\n", file, NR, indent, expected)
      err = 1
    }

    if (is_end) {
      depth = depth - 1
      if (depth < 0) depth = 0
    }

    is_if_open = (up ~ /^IF[ ]*\(.*\)[ ]*THEN\b/)
    is_do_open = (up ~ /^DO([ ]+[0-9]+|[ ]+|$)/)
    is_select_open = (up ~ /^SELECT[ ]+CASE\b/)

    if (is_if_open || is_do_open || is_select_open) {
      depth = depth + 1
    }
  }
}
END {
  if (err) exit 1
}
  ' "$f" || status=1
done

exit $status
