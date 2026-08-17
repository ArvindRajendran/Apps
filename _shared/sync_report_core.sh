#!/bin/sh
# Replace the embedded report-core block of one or more files with the
# canonical Apps/_shared/report_core.js.
#
#   sh _shared/sync_report_core.sh DAC/DAC.html BAAM/BAAM.html ...
#
# Each target must already contain the sentinel pair
#     /* == report-core v1 ==   ...   /* == /report-core == */
# (paste the canonical file once by hand when adding the core to a new app).
# For apps assembled from parts, sync the SOURCE part, not the built file.

cd "$(dirname "$0")/.." || exit 2
CANON=_shared/report_core.js
[ -f "$CANON" ] || { echo "missing $CANON"; exit 2; }
[ $# -gt 0 ] || { echo "usage: sync_report_core.sh <file> [file ...]"; exit 2; }

for f in "$@"; do
  [ -f "$f" ] || { echo "  skip  $f (not found)"; continue; }
  grep -q '== report-core' "$f" || { echo "  skip  $f (no sentinel block)"; continue; }
  grep -q '== /report-core == \*/' "$f" || { echo "  SKIP  $f (no closing sentinel)"; continue; }
  awk -v canon="$CANON" '
    /^\/\* == report-core/ { print "\x01MARK\x01"; skip = 1 }
    skip && /^\/\* == \/report-core == \*\// { skip = 0; next }
    !skip { print }
  ' "$f" | awk -v canon="$CANON" '
    /\x01MARK\x01/ { while ((getline line < canon) > 0) print line; close(canon); next }
    { print }
  ' > "$f.tmp" && mv "$f.tmp" "$f" && echo "  synced  $f"
done

echo "---"
sh _shared/check_report_core.sh
