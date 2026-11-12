#!/bin/sh
# Example 1: project, build the matrix and score it, in one pass.
# Usage: run.sh [output-dir]      (set INTRONFATE to the binary if it is not on PATH)
set -eu
cd "$(dirname "$0")"
BIN="${INTRONFATE:-intronfate}"
OUT="${1:-out}"
rm -rf "$OUT"
"$BIN" analyze --gff family.gff3 -p proteins.faa --msa family.aln.fasta -t family.nwk \
    -m dollo,mp,counts --reconstructions -o "$OUT"
