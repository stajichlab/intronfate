#!/bin/sh
# Example 3: MALIN intron tables in, scored here, and back out unchanged.
# Usage: run.sh [output-dir]      (set INTRONFATE to the binary if it is not on PATH)
set -eu
cd "$(dirname "$0")"
BIN="${INTRONFATE:-intronfate}"
OUT="${1:-out}"
rm -rf "$OUT"
# One matrix.tsv per alignment named in the table's #MAP lines.
"$BIN" convert --from-malin tiny.table -o "$OUT/matrices"
for k in K1 K2; do
    "$BIN" score --matrix "$OUT/matrices/$k/matrix.tsv" -t tiny.nwk -m dollo,mp,counts -o "$OUT/score_$k"
done
# And back: the data rows round-trip exactly (`?` cells included).
"$BIN" convert --to-malin "$OUT/matrices/K1/matrix.tsv" "$OUT/matrices/K2/matrix.tsv" -o "$OUT/roundtrip.table"
grep -v '^#' tiny.table > "$OUT/original.rows"
grep -v '^#' "$OUT/roundtrip.table" > "$OUT/roundtrip.rows"
cmp "$OUT/original.rows" "$OUT/roundtrip.rows" && echo "round trip: identical"
