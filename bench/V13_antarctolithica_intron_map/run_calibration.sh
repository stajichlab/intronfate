#!/usr/bin/env bash
# Score one calibration matrix on every placement tree of its control; prints control, variant, rep, mode, tree, dollo_losses.
#   run_calibration.sh CONTROL VARIANT REP MODE MATRIX TREES_DIR      (one line of calib/jobs.txt; INTRONFATE = the binary)
set -euo pipefail
x=$1 var=$2 rep=$3 mode=$4 m=$5 trees=$6
tmp=$(mktemp -d "${SCRATCH:-/tmp}/v13cal.XXXXXX"); trap 'rm -rf "$tmp"' EXIT
for t in "$trees"/place_*.nwk; do
  id=$(basename "$t" .nwk)
  "$INTRONFATE" score --matrix "$m" -t "$t" -m dollo -o "$tmp/$id" > /dev/null 2> "$tmp/err" || { cat "$tmp/err" >&2; exit 1; }
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$x" "$var" "$rep" "$mode" "$id" "$(awk -F'\t' '$1=="dollo_losses"{print $2}' "$tmp/$id/summary_score.tsv")"
  rm -rf "$tmp/$id"
done
