#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 48gb --time 04:00:00
#SBATCH -p highclock
#SBATCH --job-name=S1_sim
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/S1_soft_evidence_simulation/logs/S1.%j.log
# S1: simulate annotation errors and junction reads on the real strict-family matrix's sites (known truth), then fit the hard, hard + true error
# rates, and soft (k = 3, k = 0.2, and p0 0.5 / k 0.2) likelihoods and score them against the truth. Needs bench/J1_junction_support/work/run and J2b results.
# Submit: mkdir -p logs results && sbatch run.sbatch      (SEEDS="1 2 3" SCENARIOS="none low mid high" READS="empirical clean")
set -euo pipefail
: "${SCRATCH:?run inside a Slurm job}"
B=~/projects/intron_projector/bench
S1=$B/S1_soft_evidence_simulation; REPO=~/projects/intron_projector
J1W=$B/J1_junction_support/work
module load rust >/dev/null 2>&1 || true
module load gcc/12.2.0 >/dev/null 2>&1 || true
find "$REPO/crates" -name '*.rs' -exec touch {} +
export CARGO_TARGET_DIR=$SCRATCH/target
( cd "$REPO" && cargo build --release --locked )
export BIN=$SCRATCH/target/release/intronfate
cd "$S1"
python3 -I simulate.py --selftest
python3 -I evaluate.py --selftest
export NWK=$B/V6_orthology_nine_species/species.nwk
# positions a false intron could land on, seen or not: matrix sites + the constant-rate fit's expected never-seen positions
export U=$(awk -F"\t" '$1=="sites"{s=$2} $1=="ml_unobserved_positions"{u=$2} END{printf "%d", s+u}' $J1W/score_rc1_hard/summary_score.tsv)
echo "universe $U"
W=$SCRATCH/s1; rm -rf "$W"; mkdir -p "$W" results
export W S1 J1W U
SEEDS=${SEEDS:-"1 2 3"}; SCENARIOS=${SCENARIOS:-"none low mid high"}; READS=${READS:-"empirical clean"}
jobs=()
for sc in $SCENARIOS; do for rd in $READS; do
  { [ "$sc" = none ] && [ "$rd" = clean ]; } && continue          # no errors: the read model only changes true cells, so one is enough
  for sd in $SEEDS; do jobs+=("$sc $rd $sd"); done
done; done
one() {   # scenario reads seed: simulate, fit five models, evaluate
  set -euo pipefail
  set -- $1; sc=$1; rd=$2; sd=$3; d=$W/${sc}_${rd}_$sd; mkdir -p "$d"
  python3 -I "$S1/simulate.py" "$J1W/run/combined/matrix.tsv" "$J1W/score_rc1_hard/likelihood_branches.tsv" "$S1/../J2b_matrix_cells/results/J2b_bins.tsv" "$d" "$sc" "$rd" "$sd" "$U"
  SC=(--matrix "$d/matrix.tsv" -t "$NWK" -m likelihood --ml-origins --ml-rate-categories 1 --ml-max-iter 20000)
  export RAYON_NUM_THREADS=2
  "$BIN" score "${SC[@]}" -o "$d/hard" > "$d/hard.out" 2>&1
  "$BIN" score "${SC[@]}" --error-rates "$d/errors.tsv" -o "$d/hard_err" > "$d/hard_err.out" 2>&1
  "$BIN" score "${SC[@]}" --ml-soft-evidence "$d/soft_0.8_3.tsv" -o "$d/soft_k3" > "$d/soft_k3.out" 2>&1
  "$BIN" score "${SC[@]}" --ml-soft-evidence "$d/soft_0.8_0.2.tsv" -o "$d/soft_k0.2" > "$d/soft_k0.2.out" 2>&1
  "$BIN" score "${SC[@]}" --ml-soft-evidence "$d/soft_0.5_0.2.tsv" -o "$d/soft_p0.5_k0.2" > "$d/soft_p0.5_k0.2.out" 2>&1
  for m in hard hard_err soft_k3 soft_k0.2 soft_p0.5_k0.2; do
    python3 -I "$S1/evaluate.py" "$d/truth.tsv" "$d/$m" "$m" "$sc" "$rd" "$sd" > "$d/$m.row"
  done
}
export -f one
printf '%s\n' "${jobs[@]}" | xargs -P 8 -I{} bash -c 'one "{}"'
{ python3 -I evaluate.py --header; cat "$W"/*/*.row; } > results/S1_runs.tsv
cat results/S1_runs.tsv
