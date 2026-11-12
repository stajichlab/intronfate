#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 16gb --time 04:00:00
#SBATCH --job-name=V8_calibration
#SBATCH --output=logs/V8_calibration.%j.log

# Calibration of the placement scan: treat a lineage whose position is known as the unknown (Antarctolithica removed), scan all its candidate placements
# on sites (a) all (b) where it carries the intron, and ask whether the true position is recovered. Targets are intron-poor lineages like Antarctolithica
# and, as controls, intron-rich ones. Needs run_all.sbatch. Submit: sbatch -p highclock run_calibration.sbatch. Absolute path, not dirname "$0".
set -euo pipefail
: "${SCRATCH:?run inside a Slurm job}"
HERE=~/projects/intron_projector/bench/V8_antarctolithica_placement
cd "$HERE"
BIN=~/projects/intron_projector/target/release/intronfate
TARGETS=${TARGETS:-"Protomyces Taphrina Savitreella Spombe Yarrowia Lipomyces Smaydis Cneoformans Neolecta Saitoella Pyronema Orbilia Botrytis Cryomyces"}
for x in $TARGETS; do
    python3 make_trees.py results/plain/Species_Tree/SpeciesTree_rooted.txt --group-prefix $x --drop-prefix Alith --genomes genomes.tsv --out loo/$x > /dev/null
    for tier in A B; do for mode in all carrier; do
        extra=(); [ "$mode" = carrier ] && extra=(--require-carrier $x)
        python3 bootstrap_scan.py --matrix results_projection/$tier/combined/matrix.tsv --trees loo/$x --bin "$BIN" --work "$SCRATCH/loo_${x}_${tier}_${mode}" \
            --out loo/$x/${tier}_${mode}.tsv --threads "${SLURM_CPUS_PER_TASK:-8}" --reps 200 "${extra[@]}" > /dev/null
    done; done
done
for tier in A B; do echo "== tier $tier"; python3 loo_calibration.py loo $tier | tee loo/summary_$tier.tsv; done
