#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 16gb --time 04:00:00
#SBATCH --job-name=V7_error_sensitivity
#SBATCH --output=logs/V7_error_sensitivity.%j.log

# V7: how much of the Dollo gains/losses could annotation error produce? Perturb each V6 tier's matrices at the measured Sn/Pr and re-score.
# Needs the V6 tier results (bench/V6_orthology_nine_species/results_projection/<tier>) and the V4 V6-gene accuracy (bench/V4_annotation_accuracy/results_v6/<tier>).
# Submit from ~/projects/intron_projector/bench/V7_annotation_error_sensitivity/ :  mkdir -p logs && sbatch -p highclock run.sbatch
# Heavy temporary files (perturbed matrices, score outputs) live in $SCRATCH; only the tables in results/ are kept. Absolute path, not dirname "$0".
set -euo pipefail
: "${SCRATCH:?run inside a Slurm job}"
HERE=~/projects/intron_projector/bench/V7_annotation_error_sensitivity
B=~/projects/intron_projector/bench
cd "$HERE"
python3 tests/test_perturb.py
BIN=~/projects/intron_projector/target/release/intronfate
for tier in A_strict B_partial C_hog; do
    python3 run_sensitivity.py --proj "$B/V6_orthology_nine_species/results_projection/$tier" --tree "$B/V6_orthology_nine_species/species.nwk" \
        --accuracy "$B/V4_annotation_accuracy/results_v6/$tier/busco_accuracy.tsv" --bin "$BIN" --work "$SCRATCH/v7_$tier" \
        --out "results/$tier" --reps "${REPS:-20}" --threads "${SLURM_CPUS_PER_TASK:-16}" > "results/$tier.log" 2>&1 || { mkdir -p results; tail -20 "results/$tier.log"; exit 1; }
    tail -25 "results/$tier.log"
done
