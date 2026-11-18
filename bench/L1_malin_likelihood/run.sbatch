#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 16gb --time 02:00:00 -p highclock
#SBATCH --job-name=L1_malin_likelihood
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/L1_malin_likelihood/logs/L1.%j.log

# L1: compare the likelihood model (-m likelihood, one class, no annotation error) with MALIN's own constant-rates estimates.
# Layout on the HPCC (all purgeable): bench/L1_malin_likelihood/{build_src (a copy of Cargo.toml, Cargo.lock, crates/), data (MALIN's .rates files)}.
# The build goes to $SCRATCH (large temporary files). Absolute paths (no dirname $0).
set -euo pipefail
module load rust >/dev/null 2>&1 || true
D=/rhome/jstajich/projects/intron_projector/bench/L1_malin_likelihood
V1=/rhome/jstajich/projects/intron_projector/bench/V1_malin_validation/data
export CARGO_TARGET_DIR=$SCRATCH/target
find "$D/build_src/crates" -name '*.rs' -exec touch {} +   # rsync keeps laptop times; a stale binary was run once (see V2)
cd "$D/build_src"
INTRONFATE_MALIN_TABLE=$V1/rogozin_introns_table INTRONFATE_MALIN_RATES=$D/data/Rog-ecdysozoa.rates \
  cargo test --release --locked -p intronfate-core --test malin_rates -- --ignored --nocapture 2>&1 | tee "$D/results_rog.txt"
