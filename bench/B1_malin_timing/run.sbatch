#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 32gb --time 04:00:00 -p batch
#SBATCH --job-name=B1_malin_timing
#SBATCH --output=logs/B1_malin_timing.%j.log

# B1: intronfate vs MALIN timing. Submit from the folder holding this script, in the fresh clone:
#   mkdir -p logs && sbatch --export=ALL,SMOKE=1 run.sbatch                      # smoke test (seconds)
#   sbatch --export=ALL,MATRIX=/path/combined/matrix.tsv,TREE=/path/species.nwk,MALIN_JAR=/path/Malin.jar run.sbatch
# Both tools run on the same node, one after the other, so they see the same hardware. Absolute path: Slurm copies the script.
set -euo pipefail
HERE=~/projects/intronfate_main/bench/B1_malin_timing
REPO=~/projects/intronfate_main
cd "$HERE"
module load rust >/dev/null 2>&1 || true
module load java >/dev/null 2>&1 || true
( cd "$REPO" && cargo build --release --locked )
python3 -I test_summarize.py
export INTRONFATE=$REPO/target/release/intronfate
echo "node: $(hostname)  cpu: $(lscpu | sed -n 's/^Model name: *//p')  java: $(java -version 2>&1 | head -1)"
./bench.sh
