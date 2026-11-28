#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 16gb --time 12:00:00 -p batch
#SBATCH --job-name=B3_lik_iter
#SBATCH --output=logs/B3_lik_iter.%j.log

# B3: likelihood iteration counts and tolerance. Submit from this folder in the fresh clone:
#   mkdir -p logs && sbatch --export=ALL,SMOKE=1 run.sbatch ; sbatch run.sbatch
set -euo pipefail
REPO=~/projects/intronfate_main
cd "$REPO/bench/B3_likelihood_iterations"
module load rust >/dev/null 2>&1 || true
( cd "$REPO" && git pull -q --ff-only 2>/dev/null || true; cargo build --release --locked )
python3 -I test_read_summary.py && python3 -I ../B2_scaling/test_make_scaling_data.py
export INTRONFATE=$REPO/target/release/intronfate
echo "node: $(hostname)  cpu: $(lscpu | sed -n 's/^Model name: *//p')"
./iterations.sh
