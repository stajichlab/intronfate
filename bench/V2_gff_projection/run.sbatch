#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 16gb --time 02:00:00
#SBATCH --job-name=V2_gff_projection
#SBATCH --output=logs/V2_gff_projection.%j.log

# V2: real-annotation test of the GFF -> protein -> alignment projection on nine fungal genomes.
# Submit from ~/projects/intron_projector/bench/V2_gff_projection/ :
#     mkdir -p logs && sbatch -p epyc run.sbatch
# Purgeable: everything is under data/ (symlinks), results/, families/, alignments/, taxon_map.tsv.
set -euo pipefail
module load rust >/dev/null 2>&1 || true
REPO=$(cd ../.. && pwd)
# rsync keeps the laptop's file times, which can be older than a binary built a minute earlier on the
# cluster; cargo would then skip the rebuild and run a stale binary. Touch the sources first.
find "$REPO/crates" -name '*.rs' -exec touch {} +
( cd "$REPO" && cargo build --release --locked )
export INTRON_PROJECTOR=$REPO/target/release/intron-projector
. ./config.sh
./prepare.sh
# BUSCO per genome, up to 2 at a time (4 cpus each); skips genomes already done.
grep -v '^#' genomes.tsv | cut -f1 | xargs -P 2 -I{} ./busco_one.sh {}
python3 make_families.py "$MIN_GENOMES"
./align.sh
./run_projection.sh
