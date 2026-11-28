#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 32gb --time 12:00:00 -p batch
#SBATCH --job-name=B2_scaling
#SBATCH --output=logs/B2_scaling.%j.log

# B2: scaling of import and scoring in taxa and sites on synthetic data. Submit from this folder in the fresh clone:
#   mkdir -p logs && sbatch --export=ALL,SMOKE=1 run.sbatch      # smoke test (seconds)
#   sbatch run.sbatch                                            # the full grid
# Absolute paths: Slurm copies the script.
set -euo pipefail
REPO=~/projects/intronfate_main
cd "$REPO/bench/B2_scaling"
module load rust >/dev/null 2>&1 || true
( cd "$REPO" && cargo build --release --locked )
python3 -I test_make_scaling_data.py && python3 -I test_summarize.py
export INTRONFATE=$REPO/target/release/intronfate
echo "node: $(hostname)  cpu: $(lscpu | sed -n 's/^Model name: *//p')"
./scaling.sh
