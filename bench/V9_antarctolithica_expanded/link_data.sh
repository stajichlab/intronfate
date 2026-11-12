#!/usr/bin/env bash
# Link each genome's funannotate GFF3 and proteins (read-only, from BFD) into data/<label>/ for this folder.
# Own folder (not V2's) because V6 uses S. cerevisiae S288C, V2 uses AMP_1a.
set -euo pipefail
HERE=${HERE:-$PWD}; cd "$HERE"
BFD=${BFD:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation}
grep -v '^#' genomes.tsv | while IFS=$'\t' read -r label dir group; do
    src=$BFD/$dir/predict_results
    for k in gff3 proteins.fa; do
        [ -s "$src/$dir.$k.gz" ] || { echo "missing $src/$dir.$k.gz" >&2; exit 1; }
        mkdir -p data/$label; ln -sf "$src/$dir.$k.gz" "data/$label/$label.$k.gz"
    done
done
ls data
