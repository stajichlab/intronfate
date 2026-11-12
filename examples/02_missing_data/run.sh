#!/bin/sh
# Example 2: an alignment gap is missing data, not an intron loss.
#
# A, B and C share a phase-1 intron inside residue 6 (column 5). Taxon D either
#   gapped: has no residues in columns 4-5, so it cannot show that intron  -> `?`, no loss
#   full:   is aligned straight through column 5 and has no intron there   -> `0`, one loss
# Usage: run.sh [output-dir]      (set INTRON_PROJECTOR to the binary if it is not on PATH)
set -eu
cd "$(dirname "$0")"
BIN="${INTRON_PROJECTOR:-intron-projector}"
OUT="${1:-out}"
rm -rf "$OUT"
for variant in gapped full; do
    cat abc.gff3 d.$variant.gff3 > "family.$variant.gff3"
    "$BIN" analyze --gff "family.$variant.gff3" -p "proteins.$variant.faa" \
        --msa "family.$variant.aln.fasta" -t family.nwk -m dollo,mp,counts \
        -o "$OUT/$variant"
    rm -f "family.$variant.gff3"
done
