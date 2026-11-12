#!/usr/bin/env bash
# One protein FASTA per species for OrthoFinder, named <label>.faa (OrthoFinder's species name is the file stem, and
# `species.nwk` and the BED files use the same labels). Reads the BFD funannotate proteomes read-only via the data/ symlinks (link_data.sh).
set -euo pipefail
HERE=${HERE:-$PWD}; cd "$HERE"
D=${D:-$HERE/data}
mkdir -p proteomes
# Written through temp files and renamed, so two jobs preparing at once cannot hand each other a half-written file.
printf 'label\tproteins\tmean_len\tinternal_stop_fixed\n' > "proteomes/manifest.tsv.tmp$$"
grep -v '^#' genomes.tsv | while IFS=$'\t' read -r label dir group; do
    src=$D/$label/$label.proteins.fa.gz
    [ -s "$src" ] || { echo "missing $src (run link_data.sh first)" >&2; exit 1; }
    # funannotate writes one protein per gene (-T1); strip the terminal '*' and turn any internal '*' into X
    # (a premature stop); drop empty records.
    zcat "$src" | awk -v lab="$label" '
        function flush() { if (id != "" && seq != "") { n++; tot += length(seq); sub(/\*+$/, "", seq); if (seq ~ /\*/) { fixed++; gsub(/\*/, "X", seq) } print ">" id "\n" seq } }
        /^>/ { flush(); id = substr($1, 2); seq = ""; next }
        { seq = seq $0 }
        END { flush(); printf "%s\t%d\t%.0f\t%d\n", lab, n, tot / n, fixed + 0 > "/dev/stderr" }' \
        > "proteomes/$label.faa.tmp$$" 2>> "proteomes/manifest.tsv.tmp$$"
    mv "proteomes/$label.faa.tmp$$" "proteomes/$label.faa"
done
mv "proteomes/manifest.tsv.tmp$$" proteomes/manifest.tsv
cat proteomes/manifest.tsv
