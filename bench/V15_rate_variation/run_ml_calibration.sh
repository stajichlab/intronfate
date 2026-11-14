#!/usr/bin/env bash
# Fit the likelihood model to one V13 calibration matrix on every placement tree of its control; prints one line per tree:
#   control  variant  rep  categories  tree  lnL  converged  shape
#   run_ml_calibration.sh CONTROL VARIANT REP MATRIX TREES_DIR CATEGORIES [ERROR_RATES]      (INTRONFATE = the binary)
set -euo pipefail
x=$1 var=$2 rep=$3 m=$4 trees=$5 k=$6 err=${7:-}
tmp=$(mktemp -d "${SCRATCH:-/tmp}/v15cal.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
extra=(); [ -n "$err" ] && extra=(--error-rates "$err")
for t in "$trees"/place_*.nwk; do
  id=$(basename "$t" .nwk)
  "$INTRONFATE" score --matrix "$m" -t "$t" -m likelihood --ml-rate-categories "$k" "${extra[@]}" -o "$tmp/$id" \
      > /dev/null 2> "$tmp/err" || { cat "$tmp/err" >&2; exit 1; }
  s=$tmp/$id/summary_score.tsv
  get() { awk -F'\t' -v k="$1" '$1==k{print $2}' "$s"; }
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$x" "$var" "$rep" "$k" "$id" "$(get ml_log_likelihood)" "$(get ml_converged)" \
      "$(get ml_gamma_shape)"
  rm -rf "${tmp:?}/$id"
done
