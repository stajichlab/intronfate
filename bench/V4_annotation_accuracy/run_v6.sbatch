#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 1 --mem 8gb --time 00:30:00
#SBATCH --job-name=V4_v6_accuracy
#SBATCH --output=logs/V4_v6_accuracy.%j.log

# Intron accuracy over the genes in each V6 tier (A_strict, B_partial, C_hog) and the loss ceiling with the measured Sn.
# Submit from ~/projects/intron_projector/bench/V4_annotation_accuracy/ :  mkdir -p logs && sbatch -p highclock run_v6.sbatch
# Needs the V6 tier job to have run (families/, results_projection/ in bench/V6_orthology_nine_species). Absolute path, not dirname "$0".
set -euo pipefail
cd ~/projects/intron_projector/bench/V4_annotation_accuracy
SCORING=${SCORING:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/rc6_vs_pre_rc6}
V6=../V6_orthology_nine_species
for tier in A_strict B_partial C_hog; do
    python3 busco_gene_accuracy.py "$SCORING" "$V6" --genomes "$V6/genomes.tsv" --taxon-map "$V6/families/$tier/taxon_map.tsv" \
        --projection "$V6/results_projection/$tier" --out "results_v6/$tier" > "results_v6_$tier.txt"
done
ls results_v6/*
