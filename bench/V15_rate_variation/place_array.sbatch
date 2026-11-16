#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 16gb --time 12:00:00 -p epyc,intel
#SBATCH --job-name=V15_place
#SBATCH --output=logs/V15place.%A_%a.log

# Part B of run_v15.sbatch as an array, one placement tree per task (for an error set whose fits are too slow to run 4 at a time):
#   ls trees > todo.txt; sbatch --export=ALL,SET=measured,CATS=4,TODO=place_todo.txt --array=1-N place_array.sbatch
set -euo pipefail
HERE=~/projects/intron_projector/bench/V15_rate_variation
V9=~/projects/intron_projector/bench/V9_antarctolithica_expanded
V14=~/projects/intron_projector/bench/V14_likelihood_alith
REPO=${REPO:-$HOME/projects/intron_projector/bench/D1_howto_public_genomes/repo}
M=$V9/results_projection/A/combined/matrix.tsv
cd "$HERE"
t=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TODO")
id=$(basename "$t" .nwk)
d=place/${SET}_k$CATS
[ -s "$d/$id/likelihood_sites.tsv" ] && { echo "skip $id"; exit 0; }
start=$(date +%s)
RAYON_NUM_THREADS=$SLURM_CPUS_PER_TASK "$REPO/target/release/intronfate" score --matrix "$M" -t "$t" -m likelihood \
    --error-rates "$V14/errors/$SET.tsv" --ml-rate-categories "$CATS" -o "$d/$id" > "$d/$id.log" 2>&1
echo "$id: $(( $(date +%s) - start )) s; $(grep -E '^ml_(log_likelihood|converged)' "$d/$id/summary_score.tsv" | tr '\n' ' ')"
