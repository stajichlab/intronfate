#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 16gb --time 02:00:00 -p batch
#SBATCH --job-name=V1_malin_validation
#SBATCH --output=logs/V1_malin_validation.%j.log

# V1: validate intronfate against MALIN's published test data (read-only inputs, purgeable output).
# Submit from ~/projects/intron_projector/bench/V1_malin_validation/ :
#     mkdir -p logs && MALIN_TEST=/path/to/Malin/test sbatch --export=ALL run.sbatch
set -euo pipefail
module load rust >/dev/null 2>&1 || true
REPO=$(cd ../.. && pwd)
find "$REPO/crates" -name '*.rs' -exec touch {} +   # see bench/V2_gff_projection/run.sbatch
( cd "$REPO" && cargo build --release --locked )
: "${MALIN_TEST:?set MALIN_TEST}"
export INTRONFATE=$REPO/target/release/intronfate
./validate.sh
