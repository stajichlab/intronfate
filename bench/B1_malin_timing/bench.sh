#!/usr/bin/env bash
# B1: wall time and peak memory of intronfate vs MALIN's command-line programs on the same matrix.
#
# Compared (the same data, the same tree, REPS repeats each):
#   dollo       intronfate `score -m dollo`            vs  MALIN DolloParsimony
#   likelihood  intronfate `score -m likelihood`       vs  MALIN calculateRASRates (constant rates)
# MALIN has no command-line projection step (its tables are built in the GUI), so projection is timed for
# intronfate only (bench/V3) and the table export is timed here as `export_malin_table`.
# Convergence rules differ, so every run keeps its output and the final log-likelihoods are compared in README.
#
# Purgeable: everything is written under $OUT (default ./results). Run it on the HPCC (see run.sbatch).
#   SMOKE=1 ./bench.sh        tiny fixtures, seconds, any machine; MALIN is skipped if MALIN_JAR is unset
#
# Environment: INTRONFATE (binary), MALIN_JAR, MATRIX (one matrix.tsv), TREE (Newick), REPS (default 3), OUT.
set -euo pipefail
cd "$(dirname "$0")"
HERE=$PWD; REPO=$(cd ../.. && pwd)
BIN=${INTRONFATE:-$REPO/target/release/intronfate}
OUT=${OUT:-$HERE/results}; REPS=${REPS:-3}
rm -rf "$OUT"; mkdir -p "$OUT"
TSV="$OUT/timings.tsv"; printf 'tool\ttask\trep\twall_s\tmax_rss_kb\texit\n' > "$TSV"

if [ "${SMOKE:-0}" = 1 ]; then
    "$BIN" convert --from-malin "$REPO/tests/data/malin_binary.table" -o "$OUT/smoke_matrix" > /dev/null
    MATRIX=$(ls "$OUT"/smoke_matrix/*/matrix.tsv | head -1); TREE=$REPO/examples/03_malin/tiny.nwk; REPS=1
fi
: "${MATRIX:?set MATRIX to one matrix.tsv (for example OUT/combined/matrix.tsv from analyze --combine)}"
: "${TREE:?set TREE}"
for f in "$BIN" "$MATRIX" "$TREE"; do [ -e "$f" ] || { echo "missing: $f" >&2; exit 1; }; done
echo "intronfate: $("$BIN" --version)   matrix: $MATRIX   reps: $REPS"

timed() {  # timed <tool> <task> <rep> <stdout-file> <command...>
    local tool=$1 task=$2 rep=$3 out=$4; shift 4
    local tf; tf=$(mktemp); local rc=0
    /usr/bin/time -f '%e\t%M' -o "$tf" "$@" > "$out" 2> "$out.err" || rc=$?
    printf '%s\t%s\t%s\t%s\t%s\n' "$tool" "$task" "$rep" "$(tail -1 "$tf")" "$rc" >> "$TSV"; rm -f "$tf"
}
# the awk below moves `exit` to the last column: timed rows are tool task rep wall rss exit
TABLE=$OUT/matrix.malin.table
timed intronfate export_malin_table 1 "$OUT/export.out" "$BIN" convert --to-malin "$MATRIX" -o "$TABLE"

for r in $(seq 1 "$REPS"); do
    timed intronfate dollo "$r" "$OUT/if_dollo.$r.out" "$BIN" score --matrix "$MATRIX" -t "$TREE" -m dollo -o "$OUT/if_dollo.$r"
    timed intronfate likelihood "$r" "$OUT/if_lik.$r.out" "$BIN" score --matrix "$MATRIX" -t "$TREE" -m likelihood -o "$OUT/if_lik.$r"
done

if [ -n "${MALIN_JAR:-}" ] && [ -e "$MALIN_JAR" ]; then
    P=ca.umontreal.iro.evolution.introns
    J="java -Xmx${MALIN_HEAP:-8g} -cp $MALIN_JAR"
    for r in $(seq 1 "$REPS"); do
        timed malin dollo "$r" "$OUT/malin_dollo.$r.out" $J $P.DolloParsimony "$TABLE" "$TREE"
        timed malin likelihood "$r" "$OUT/malin_rates.$r.out" $J $P.calculateRASRates -opt-eps "${MALIN_EPS:-1e-6}" "$TABLE" "$TREE"
    done
else
    echo "MALIN_JAR not set or missing: MALIN rows skipped" >&2
fi
python3 -I "$HERE/summarize.py" "$TSV" | tee "$OUT/summary.md"
