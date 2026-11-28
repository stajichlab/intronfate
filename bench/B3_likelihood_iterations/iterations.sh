#!/usr/bin/env bash
# B3: why does `score -m likelihood` take so long on the B2 simulated data, and how much do looser tolerances save?
# For each (TAXA, SITES) cell: generate the B2 data, then run the likelihood with each tolerance pair, one run each (the run is long),
# recording wall time, peak RSS, EM iterations, final lnL and whether it converged. `-vv` keeps the per-cycle trace (E-step and M-step
# time) of the default run in $OUT/trace_<T>_<S>.txt. Purgeable: everything is under $OUT.
#   SMOKE=1 ./iterations.sh      tiny cells, seconds
# Environment: INTRONFATE, CELLS ("taxa:sites ..."), TOLS ("tol:param_tol ..."), OUT, CELL_TIMEOUT (1800).
set -euo pipefail
cd "$(dirname "$0")"
HERE=$PWD; REPO=$(cd ../.. && pwd)
BIN=${INTRONFATE:-$REPO/target/release/intronfate}
OUT=${OUT:-$HERE/results}
CELLS=${CELLS:-"10:10000 30:1000 30:10000 100:1000 300:1000"}
TOLS=${TOLS:-"1e-6:1e-7 1e-4:1e-5 1e-2:1e-3"}
if [ "${SMOKE:-0}" = 1 ]; then CELLS="6:30 12:60"; fi
CELL_TIMEOUT=${CELL_TIMEOUT:-1800}
[ -x "$BIN" ] || { echo "missing binary: $BIN" >&2; exit 1; }
rm -rf "$OUT"; mkdir -p "$OUT"
TSV=$OUT/iterations.tsv
printf 'taxa\tsites\ttol\tparam_tol\twall_s\tmax_rss_kb\texit\titerations\tlnl\tconverged\n' > "$TSV"
echo "intronfate: $("$BIN" --version)  node: $(hostname)"
for c in $CELLS; do
  T=${c%%:*}; S=${c##*:}; D=$OUT/data_${T}_${S}; mkdir -p "$D"
  python3 -I "$HERE/../B2_scaling/make_scaling_data.py" "$D/d" "$T" "$S"
  "$BIN" convert --from-malin "$D/d.table" -o "$D/m" > /dev/null
  M=$(ls "$D"/m/*/matrix.tsv | head -1)
  for p in $TOLS; do
    tol=${p%%:*}; ptol=${p##*:}; o=$D/o_$tol; rc=0; tf=$(mktemp)
    /usr/bin/time -f '%e\t%M' -o "$tf" timeout "$CELL_TIMEOUT" "$BIN" score --matrix "$M" -t "$D/d.nwk" -m likelihood \
        --ml-tol "$tol" --ml-param-tol "$ptol" -vv -o "$o" > "$o.out" 2> "$o.err" || rc=$?
    [ "$p" = "${TOLS%% *}" ] && cp "$o.err" "$OUT/trace_${T}_${S}.txt"
    it=$(python3 -I "$HERE/read_summary.py" "$o" ml_iterations 2>/dev/null || echo NA)
    ll=$(python3 -I "$HERE/read_summary.py" "$o" ml_log_likelihood 2>/dev/null || echo NA)
    conv=$(python3 -I "$HERE/read_summary.py" "$o" ml_converged 2>/dev/null || echo NA)
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$T" "$S" "$tol" "$ptol" "$(tail -1 "$tf")" "$rc" "$it" "$ll" "$conv" >> "$TSV"
    rm -f "$tf"
  done
  rm -rf "$D"/m "$D"/d.table
done
cat "$TSV"
