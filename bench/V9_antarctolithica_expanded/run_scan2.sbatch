#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 16gb --time 04:00:00
#SBATCH --job-name=V9_bootstrap_scan
#SBATCH --output=logs/V9_bootstrap_scan.%j.log

# Placement scan with a family bootstrap, on all sites and on sites where an Antarctolithica strain carries the intron, for the gene-tree backbone
# (trees/) and a backbone with Taphrinomycotina forced monophyletic (trees_taphmono/). Needs run_all.sbatch. Submit: sbatch -p highclock run_scan2.sbatch
set -euo pipefail
: "${SCRATCH:?run inside a Slurm job}"
HERE=~/projects/intron_projector/bench/V9_antarctolithica_expanded
cd "$HERE"
BIN=~/projects/intron_projector/target/release/intronfate
python3 make_trees.py results/plain/Species_Tree/SpeciesTree_rooted.txt --group-prefix Alith --genomes genomes.tsv --root-on Mmucedo --out trees_taphmono --taph-monophyletic
mkdir -p scan
for tier in A B; do for bb in trees trees_taphmono; do for mode in all carrier; do
    extra=(); [ "$mode" = carrier ] && extra=(--require-carrier Alith)
    python3 bootstrap_scan.py --matrix results_projection/$tier/combined/matrix.tsv --trees $bb --bin "$BIN" --work "$SCRATCH/s_${tier}_${bb}_${mode}" \
        --out scan/${tier}_${bb}_${mode}.tsv --threads "${SLURM_CPUS_PER_TASK:-8}" "${extra[@]}" > /dev/null
    echo "== tier $tier, backbone $bb, sites: $mode"; head -9 scan/${tier}_${bb}_${mode}.tsv | cut -f1-8
done; done; done
