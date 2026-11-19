#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 16gb --time 04:00:00 -p highclock
#SBATCH --job-name=L1_classes
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/L1_malin_likelihood/logs/L1_classes.%j.log

# L1b: rate classes on the Rogozin table (7 taxa, 7,221 sites): none, shared factor, separate loss / gain classes. Compare ml_log_likelihood.
# MALIN's published rate-variation run on this table (Rog-ecdysozoa1000.rates) has two loss classes and one gain class.
# Layout as run.sbatch (build_src, data); build in $SCRATCH; absolute paths.
set -euo pipefail
module load rust >/dev/null 2>&1 || true
D=/rhome/jstajich/projects/intron_projector/bench/L1_malin_likelihood
V1=/rhome/jstajich/projects/intron_projector/bench/V1_malin_validation/data
export CARGO_TARGET_DIR=$SCRATCH/target
find "$D/build_src/crates" -name '*.rs' -exec touch {} +
( cd "$D/build_src" && cargo build --release --locked )
BIN=$CARGO_TARGET_DIR/release/intronfate
mkdir -p "$D/classes" && cd "$D/classes"
"$BIN" convert --from-malin "$V1/rogozin_introns_table" -o matrices
M=$(ls matrices/*/matrix.tsv | head -1)
run() {  # name, options...
  n=$1; shift
  /usr/bin/time -f "$n %e s" "$BIN" score --matrix "$M" -t "$V1/Rogozin-ecdysozoa.tre" -m likelihood --min-taxa 1 "$@" -o "out_$n" > /dev/null
  s=out_$n/summary_score.tsv
  printf '%s\t' "$n"; grep -E "^ml_(log_likelihood|iterations|converged|rate_categories|gamma_shape|gain_categories|gain_shape|loss_categories|loss_shape)\b" "$s" | tr '\n' ' '; echo
}
{
run none
run shared2 --ml-rate-categories 2
run shared4 --ml-rate-categories 4
run loss2 --ml-loss-categories 2
run loss3 --ml-loss-categories 3
run loss4 --ml-loss-categories 4
run gain2 --ml-gain-categories 2
run gain2_loss2 --ml-gain-categories 2 --ml-loss-categories 2
run gain3_loss3 --ml-gain-categories 3 --ml-loss-categories 3
} 2>&1 | tee "$D/results_classes.txt"
