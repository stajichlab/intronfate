#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 16gb --time 08:00:00 -p highclock
#SBATCH --job-name=L1_rates
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/L1_malin_likelihood/logs/L1_rates.%j.log

# L1c: the free-class likelihood model against MALIN's rate-variation estimates (calculateRASRates): Rog-ecdysozoa1000.rates (Rogozin
# table, two loss levels) and Chromo23.rates (Chromo23 table, 23 taxa; MALIN reports 7,030 sites and log-likelihood -28430.9249).
# Layout as run.sbatch (build_src, data). Build in $SCRATCH. Absolute paths.
set -euo pipefail
module load rust >/dev/null 2>&1 || true
D=/rhome/jstajich/projects/intron_projector/bench/L1_malin_likelihood
V1=/rhome/jstajich/projects/intron_projector/bench/V1_malin_validation/data
export CARGO_TARGET_DIR=$SCRATCH/target
# NS=global models the never-seen positions as MALIN does (SPEC §5.5); default is this tool's per-pattern conditioning
export INTRONFATE_MALIN_NEVER_SEEN=${NS:-per-pattern}
find "$D/build_src/crates" -name '*.rs' -exec touch {} +
cd "$D/build_src"
T() { set +e; cargo test --release --locked -p intronfate-core --test malin_rates -- --ignored --nocapture 2>&1 | grep -v -E "Compiling|Finished|Running|^$|^warning|running 1 test"; set -e; return 0; }
{
echo "== Rogozin table, Rog-ecdysozoa1000.rates"
INTRONFATE_MALIN_TABLE=$V1/rogozin_introns_table INTRONFATE_MALIN_RATES=$D/data/Rog-ecdysozoa1000.rates T
echo "== Chromo23: positions kept by each minimum-observed filter (MALIN used 7,030)"
INTRONFATE_MALIN_COUNT_ONLY=1 INTRONFATE_MALIN_TABLE=$V1/Chromo23.table INTRONFATE_MALIN_RATES=$D/data/Chromo23.rates T
# Chromo23 fits at the minimum-observed filters named in MINOBS (space separated; set after reading the counts above)
for n in ${MINOBS:-}; do
  echo "== Chromo23 table, Chromo23.rates, min observed taxa $n (MALIN: 7,030 positions, log-likelihood -28430.9249)"
  INTRONFATE_MALIN_MIN_OBSERVED=$n INTRONFATE_MALIN_LL=-28430.9249 INTRONFATE_MALIN_TABLE=$V1/Chromo23.table INTRONFATE_MALIN_RATES=$D/data/Chromo23.rates T
done
} 2>&1 | tee "$D/results_rates_${TAG:-1}.txt"
