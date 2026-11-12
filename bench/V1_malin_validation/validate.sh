#!/usr/bin/env bash
# V1: validate intron-projector against MALIN's own published test data.
#
# Runs on the HPCC (see run.sbatch). Reads the MALIN distribution's test files, never writes outside
# its own results/ folder, and is purgeable: `rm -rf results` removes everything it makes.
#
#   SMOKE=1 ./validate.sh      tiny synthetic fixtures (tests/data, examples/03_malin); seconds, any machine
#   ./validate.sh              Malin/test/{rogozin_introns_table,Rogozin-ecdysozoa.tre,Chromo23.*}
#
# Environment: INTRON_PROJECTOR (binary), MALIN_TEST (directory holding the MALIN test files).
set -euo pipefail
cd "$(dirname "$0")"
HERE=$PWD
REPO=$(cd ../.. && pwd)
BIN=${INTRON_PROJECTOR:-$REPO/target/release/intron-projector}
OUT=${RESULTS:-$HERE/results}
rm -rf "$OUT"; mkdir -p "$OUT"
TIMES="$OUT/timings.tsv"; printf 'step\tseconds\n' > "$TIMES"

timed() {  # timed <label> <command...>
    local label=$1; shift
    local t0; t0=$(date +%s.%N)
    "$@"
    awk -v a="$t0" -v b="$(date +%s.%N)" -v l="$label" 'BEGIN{printf "%s\t%.2f\n", l, b-a}' >> "$TIMES"
}

if [ "${SMOKE:-0}" = 1 ]; then
    BINARY_TABLE=$REPO/tests/data/malin_binary.table; BINARY_TREE=$REPO/examples/03_malin/tiny.nwk
    MAP_TABLE=$REPO/tests/data/malin_tiny.table;      MAP_TREE=$REPO/examples/03_malin/tiny.nwk
else
    M=${MALIN_TEST:?set MALIN_TEST to the directory containing rogozin_introns_table, Chromo23.table, ...}
    BINARY_TABLE=$M/rogozin_introns_table; BINARY_TREE=$M/Rogozin-ecdysozoa.tre
    MAP_TABLE=$M/Chromo23.table;           MAP_TREE=$M/Chromo23.newick
fi
for f in "$BIN" "$BINARY_TABLE" "$BINARY_TREE" "$MAP_TABLE" "$MAP_TREE"; do
    [ -e "$f" ] || { echo "missing input: $f" >&2; exit 1; }
done
echo "binary: $("$BIN" --version)   smoke=${SMOKE:-0}"

# ---- 1. Rogozin et al. 2003 table (0/1 only, no #MAP) --------------------------------------------
R=$OUT/rogozin
timed rogozin_convert "$BIN" convert --from-malin "$BINARY_TABLE" -o "$R/matrices"
RM=$(ls "$R"/matrices/*/matrix.tsv)
timed rogozin_score "$BIN" score --matrix "$RM" -t "$BINARY_TREE" -m dollo,mp,counts -o "$R/score"
# Independent second importer (Python) must produce the same matrix body.
python3 -I "$REPO/scripts/malin_to_matrix.py" "$BINARY_TABLE" "$R/python_matrix.tsv" > "$R/python_convert.log"
diff <(grep -v '^#' "$RM") <(grep -v '^#' "$R/python_matrix.tsv") > "$R/importer_parity.diff" \
    && echo "importer parity (rogozin): identical" || { echo "importer parity (rogozin): DIFFERS, see $R/importer_parity.diff" >&2; exit 1; }

# ---- 2. Chromo23 (phase codes + #MAP over many alignments) ---------------------------------------
C=$OUT/chromo23
timed chromo23_convert "$BIN" convert --from-malin "$MAP_TABLE" -o "$C/matrices"
mkdir -p "$C/score"
t0=$(date +%s.%N)
for d in "$C"/matrices/*/; do
    k=$(basename "$d")
    "$BIN" score --matrix "$d/matrix.tsv" -t "$MAP_TREE" -m dollo,mp,counts -o "$C/score/$k" > /dev/null
done
awk -v a="$t0" -v b="$(date +%s.%N)" 'BEGIN{printf "chromo23_score_all\t%.2f\n", b-a}' >> "$TIMES"
timed chromo23_export "$BIN" convert --to-malin "$C"/matrices/*/matrix.tsv -o "$C/roundtrip.table"

# ---- 3. Compare, summarise -----------------------------------------------------------------------
python3 -I "$HERE/summarize.py" "$OUT" "$BINARY_TABLE" "$MAP_TABLE" | tee "$OUT/REPORT.txt"
