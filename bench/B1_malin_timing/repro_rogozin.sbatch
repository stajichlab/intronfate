#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 16gb --time 02:00:00 -p batch
#SBATCH --job-name=B1_rog_repro
#SBATCH --output=logs/B1_rog_repro.%j.log

# B1: make MALIN's command line reproduce its published Rogozin fit. For each variant of calculateRASRates switches, fit MALIN on
# the Rogozin table, then evaluate intronfate's log-likelihood at MALIN's parameters (bench/L1's test). The published
# Rog-ecdysozoa.rates is the control: L1 gets -17817.3652 at those parameters (intronfate's own maximum is -17817.3633).
# Needs results_rog/{matrix.malin.table,tree.labelled.nwk} from run_lnl.sbatch. Output: repro_rog.tsv. Absolute paths.
set -euo pipefail
H=~/projects/intronfate_main; HERE=$H/bench/B1_malin_timing; R=$HERE/results_rog
PUB=~/projects/intron_projector/bench/L1_malin_likelihood/data/Rog-ecdysozoa.rates
PUB1000=~/projects/intron_projector/bench/L1_malin_likelihood/data/Rog-ecdysozoa1000.rates   # calculateRASRates, 2 loss levels (L1: lnL -17796.1286 at its parameters)
JAR=$H/Malin/Malin.jar; P=ca.umontreal.iro.evolution.introns
module load rust >/dev/null 2>&1 || true; module load java >/dev/null 2>&1 || true
export CARGO_TARGET_DIR=$SCRATCH/target
cd "$H"
evalrates() {  # evalrates <rates file> -> "max<TAB>at_malin"
  INTRONFATE_MALIN_TABLE=$R/matrix.malin.table INTRONFATE_MALIN_RATES=$1 INTRONFATE_MALIN_NEVER_SEEN=global \
    cargo test --release --locked -p intronfate-core --test malin_rates -- --ignored --nocapture 2>&1 \
    | sed -n "s/.*maximum \(-[0-9.]*\).*MALIN's parameters \(-[0-9.]*\).*/\1\t\2/p" | head -1
}
printf 'variant\tswitches\twall_s\tintronfate_max\tlnL_at_malin\n' > "$HERE/repro_rog.tsv"
printf 'published_control\t(none, evaluate Rog-ecdysozoa.rates)\t-\t%s\n' "$(evalrates "$PUB")" >> "$HERE/repro_rog.tsv"
printf 'published1000_control\t(none, evaluate Rog-ecdysozoa1000.rates)\t-\t%s\n' "$(evalrates "$PUB1000")" >> "$HERE/repro_rog.tsv"
run() {  # run <name> <switches...>
  local name=$1; shift
  local out=$R/rates.$name.out t0; t0=$(date +%s.%N)
  java -Xmx8g -cp "$JAR" $P.calculateRASRates "$@" "$R/matrix.malin.table" "${TREEFILE:-$R/tree.labelled.nwk}" > "$out" 2> "$out.err" || echo "variant $name exited non-zero" >&2
  local w; w=$(awk -v a="$t0" -v b="$(date +%s.%N)" 'BEGIN{printf "%.1f", b-a}')
  printf '%s\t%s\t%s\t%s\n' "$name" "$*" "$w" "$(evalrates "$out")" >> "$HERE/repro_rog.tsv"
}
run default
run eps1e-10 -opt-eps 1e-10
run round1000 -opt-round 1000
run eps_round -opt-eps 1e-10 -opt-round 1000
run start_pub -start-with "$PUB" -opt-eps 1e-10 -opt-round 1000
run lowloss1 -low-loss-levels 1
# MALIN takes branch lengths in the Newick file as STARTING values; a tree without them leaves internal branches stuck near 0.
# (a) the published fit's lengths, (b) a uniform guess on every branch.
sed -n 's/^#tree //p' "$PUB" | head -1 | tr -d ' ' > "$R/tree.publen.nwk"
python3 - "$R/tree.labelled.nwk" "$R/tree.uniform.nwk" <<'PY'
import re, sys
t = open(sys.argv[1]).read().strip().rstrip(';')
open(sys.argv[2], 'w').write(re.sub(r'(?<=[A-Za-z0-9])(?=[,)])', ':0.5', t) + ':0.5;\n')
PY
TREEFILE=$R/tree.publen.nwk run publen -opt-eps 1e-8 -opt-round 1000
TREEFILE=$R/tree.uniform.nwk run uniform -opt-eps 1e-8 -opt-round 1000
TREEFILE=$R/tree.uniform.nwk run uniform_lowloss1 -low-loss-levels 1 -opt-eps 1e-8 -opt-round 1000
cat "$HERE/repro_rog.tsv"
