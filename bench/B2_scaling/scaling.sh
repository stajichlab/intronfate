#!/usr/bin/env bash
# B2: how intronfate's import, Dollo and likelihood scoring scale with taxa and sites, on synthetic data (see make_scaling_data.py).
# For each (TAXA, SITES) cell: generate a MALIN table + tree (not timed), then time, REPS repeats each:
#   import   `convert --from-malin`            (ease-of-import cost: MALIN table -> matrix.tsv)
#   dollo    `score -m dollo`
#   mp       `score -m mp`
#   likelihood `score -m likelihood`           (skipped above LIK_MAX_CELLS cells; killed after CELL_TIMEOUT seconds)
# Purgeable: everything is under $OUT (default ./results); data are deleted after each cell unless KEEP=1.
#   SMOKE=1 ./scaling.sh runs a tiny grid: TAXA_LIST="6 12" SITES_LIST="30 60" REPS=1, seconds, any machine.
# Environment: INTRONFATE (binary), TAXA_LIST, SITES_LIST, REPS (3), OUT, CELL_TIMEOUT (1800), LIK_MAX_CELLS (2e7), KEEP.
set -euo pipefail
cd "$(dirname "$0")"
HERE=$PWD; REPO=$(cd ../.. && pwd)
BIN=${INTRONFATE:-$REPO/target/release/intronfate}
OUT=${OUT:-$HERE/results}; REPS=${REPS:-3}
if [ "${SMOKE:-0}" = 1 ]; then TAXA_LIST="6 12"; SITES_LIST="30 60"; REPS=1; fi
TAXA_LIST=${TAXA_LIST:-"10 30 100 300 1000"}; SITES_LIST=${SITES_LIST:-"1000 10000 100000"}
CELL_TIMEOUT=${CELL_TIMEOUT:-1800}; LIK_MAX_CELLS=${LIK_MAX_CELLS:-20000000}
[ -x "$BIN" ] || { echo "missing binary: $BIN" >&2; exit 1; }
rm -rf "$OUT"; mkdir -p "$OUT"
TSV=$OUT/timings.tsv; printf 'task\ttaxa\tsites\trep\twall_s\tmax_rss_kb\texit\n' > "$TSV"
echo "intronfate: $("$BIN" --version)  node: $(hostname)  reps: $REPS"

timed() {  # timed <task> <taxa> <sites> <rep> <command...>
    local task=$1 t=$2 s=$3 rep=$4; shift 4
    local tf rc=0; tf=$(mktemp)
    /usr/bin/time -f '%e\t%M' -o "$tf" timeout "$CELL_TIMEOUT" "$@" > /dev/null 2> "$OUT/last.err" || rc=$?
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$task" "$t" "$s" "$rep" "$(tail -1 "$tf")" "$rc" >> "$TSV"; rm -f "$tf"
}

for T in $TAXA_LIST; do
  for S in $SITES_LIST; do
    D=$OUT/data_${T}_${S}; mkdir -p "$D"
    python3 -I "$HERE/make_scaling_data.py" "$D/d" "$T" "$S"
    echo "cell taxa=$T sites=$S"
    for r in $(seq 1 "$REPS"); do
        rm -rf "$D/m"
        timed import "$T" "$S" "$r" "$BIN" convert --from-malin "$D/d.table" -o "$D/m"
    done
    M=$(ls "$D"/m/*/matrix.tsv | head -1)
    for r in $(seq 1 "$REPS"); do
        timed dollo "$T" "$S" "$r" "$BIN" score --matrix "$M" -t "$D/d.nwk" -m dollo -o "$D/o"
        timed mp "$T" "$S" "$r" "$BIN" score --matrix "$M" -t "$D/d.nwk" -m mp -o "$D/o"
        if [ $((T * S)) -le "$LIK_MAX_CELLS" ]; then
            timed likelihood "$T" "$S" "$r" "$BIN" score --matrix "$M" -t "$D/d.nwk" -m likelihood -o "$D/o"
        fi
    done
    [ "${KEEP:-0}" = 1 ] || rm -rf "$D"
  done
done
python3 -I "$HERE/summarize.py" "$TSV" | tee "$OUT/summary.md"
