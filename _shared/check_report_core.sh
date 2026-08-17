#!/bin/sh
# Verify that every app's embedded copy of the report core is byte-identical
# to the canonical Apps/_shared/report_core.js.
#
#   sh _shared/check_report_core.sh          (run from the Apps folder)
#
# The apps are single-file by design, so the core is pasted into each one
# between the sentinel comments:
#     /* == report-core v1 ==   ...   /* == /report-core == */
# This script extracts the block from each app and diffs it against the
# canonical file.  Exit status 0 = all copies match.

cd "$(dirname "$0")/.." || exit 2
CANON=_shared/report_core.js
[ -f "$CANON" ] || { echo "missing $CANON"; exit 2; }

# temp files live beside the canonical copy: $TMPDIR is not always writable
REF=_shared/.core_ref.$$
TMP=_shared/.core_cmp.$$
trap 'rm -f "$REF" "$TMP"' EXIT INT TERM
awk '/^\/\* == report-core/,/^\/\* == \/report-core == \*\//' "$CANON" > "$REF"
[ -s "$REF" ] || { echo "no sentinel block in $CANON"; exit 2; }

fail=0
found=0
for app in */*.html; do
  grep -q '== report-core' "$app" || continue
  found=$((found + 1))
  awk '/^\/\* == report-core/,/^\/\* == \/report-core == \*\//' "$app" > "$TMP"
  if cmp -s "$REF" "$TMP"; then
    echo "  ok    $app"
  else
    echo "  DIFF  $app"
    diff "$REF" "$TMP" | head -20
    fail=1
  fi
done

echo "---"
if [ "$found" -eq 0 ]; then
  echo "no app carries the report core yet"
  exit 0
fi
if [ "$fail" -eq 0 ]; then
  echo "$found copies, all identical to $CANON"
else
  echo "copies out of sync — re-paste from $CANON"
fi
exit $fail
