#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 8gb --time 24:00:00 -p epyc,intel
#SBATCH --job-name=V15_cal
#SBATCH --output=logs/V15cal.%A_%a.log

# V15 likelihood calibration: one array task per V13 calibration matrix ("all" mode: every position; the carrier-only matrices
# break the model's conditioning on what is observed and are not used). Each matrix is fitted on all placement trees of its
# control, without annotation error, with CATS rate categories (default 1). Resumable: finished outputs are skipped.
# JOBS: the job list (default calib_jobs.txt, V13's "all" lines); REPO: the built copy; ML_OPTS: extra `score` options.
#   sbatch --array=1-165%40 calib_ml.sbatch
#   sbatch --export=ALL,CATS=4,JOBS=calib_jobs_k4.txt,REPO=$PWD/repo,ML_OPTS="--ml-param-tol 1e-5" --array=1-60%30 calib_ml.sbatch
set -euo pipefail
HERE=~/projects/intron_projector/bench/V15_rate_variation
REPO=${REPO:-$HOME/projects/intron_projector/bench/D1_howto_public_genomes/repo}
export INTRONFATE=$REPO/target/release/intronfate RAYON_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8} ML_OPTS=${ML_OPTS:-}
JOBS=${JOBS:-calib_jobs.txt}
CATS=${CATS:-1}
cd "$HERE"
line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBS")
IFS=$'\t' read -r x var rep mode m trees <<< "$line"
V13=~/projects/intron_projector/bench/V13_antarctolithica_intron_map    # V13's job list has paths relative to its folder
[[ $m = /* ]] || m=$V13/$m
[[ $trees = /* ]] || trees=$V13/$trees
mkdir -p calib
for k in $CATS; do
  out=calib/${x}_${var}_${rep}_k$k.tsv
  [ -s "$out" ] && { echo "skip $out"; continue; }
  start=$(date +%s)
  ./run_ml_calibration.sh "$x" "$var" "$rep" "$m" "$trees" "$k" > "$out.part"
  mv "$out.part" "$out"
  echo "$out: $(wc -l < "$out") trees, $(( $(date +%s) - start )) s"
done
