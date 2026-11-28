#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 16gb --time 02:00:00 -p batch
#SBATCH --job-name=B1_lnl
#SBATCH --output=logs/B1_lnl.%j.log

# B1 follow-up: is the speed comparison like for like? Fit MALIN (calculateRASRates) on a tree with named internal nodes,
# then evaluate intronfate's log-likelihood at MALIN's fitted parameters and compare with intronfate's own maximum, as bench/L1 did.
# Submit:  sbatch --export=ALL,RUN=<bench.sh OUT dir with matrix.malin.table>,TREE=<tree>,MALIN_JAR=<Malin.jar> run_lnl.sbatch
# Uses the ignored test crates/intronfate-core/tests/malin_rates.rs. Build goes to $SCRATCH. Absolute paths (Slurm copies the script).
set -euo pipefail
REPO=~/projects/intronfate_main
HERE=$REPO/bench/B1_malin_timing
: "${RUN:?}" "${TREE:?}" "${MALIN_JAR:?}"
module load rust >/dev/null 2>&1 || true
module load java >/dev/null 2>&1 || true
export CARGO_TARGET_DIR=$SCRATCH/target
python3 -I "$HERE/test_label_internal_nodes.py"
python3 -I "$HERE/label_internal_nodes.py" "$TREE" "$RUN/tree.labelled.nwk"
P=ca.umontreal.iro.evolution.introns
# MALIN `-ambiguous d` ignores columns with more than d*n question marks (n = 9 taxa); the matching intronfate filter is
# INTRONFATE_MALIN_MIN_OBSERVED = taxa observed >= n - floor(d*n). Pairs "d:min_observed" (9 taxa): d=0 -> 9, d=0.5 -> 5.
cd "$REPO"
for pair in ${PAIRS:-0:9 0.5:5}; do
  d=${pair%%:*}; k=${pair##*:}
  out=$RUN/malin_rates.amb$d${TAG:-}.out
  ( cd "$HERE" && /usr/bin/time -f "malin calculateRASRates -ambiguous $d: wall %e s, max RSS %M kB" \
      java -Xmx8g -cp "$MALIN_JAR" $P.calculateRASRates -ambiguous "$d" ${MALIN_EXTRA:-} -opt-eps "${MALIN_EPS:-1e-6}" "$RUN/matrix.malin.table" "$RUN/tree.labelled.nwk" > "$out" )
  for mode in per_pattern global; do
    echo "=== -ambiguous $d (min observed $k), never-seen: $mode"
    INTRONFATE_MALIN_TABLE=$RUN/matrix.malin.table INTRONFATE_MALIN_RATES=$out \
    INTRONFATE_MALIN_NEVER_SEEN=$mode INTRONFATE_MALIN_MIN_OBSERVED=$k \
      cargo test --release --locked -p intronfate-core --test malin_rates -- --ignored --nocapture 2>&1 \
      | grep -v -E "^\s*(Compiling|Finished|Running)" | tee "$HERE/lnl_amb${d}${TAG:-}_$mode.txt"
  done
done
