#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 48gb --time 06:00:00
#SBATCH -p highclock
#SBATCH --job-name=J3b_creep
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/J3b_creeping_parameters/logs/J3b.%j.log
# J3 follow-up: which parameters are still moving when a 4-class fit stops at an iteration cap, and what tolerance is safe?
# The same 4-class fit (hard, and soft p0 0.8 / k 3) is stopped at 500 ... 50,000 iterations and with looser tolerances; the saved branch tables, summaries
# and per-site origin tables are compared with the 50,000-iteration fit by analyse.py. Needs bench/J1_junction_support/work (run + diag_soft.tsv).
# Submit: mkdir -p logs results && sbatch run.sbatch
set -euo pipefail
: "${SCRATCH:?run inside a Slurm job}"
H=~/projects/intron_projector/bench; W=$H/J1_junction_support/work; HERE=$H/J3b_creeping_parameters
module load rust >/dev/null 2>&1 || true; module load gcc/12.2.0 >/dev/null 2>&1 || true
export CARGO_TARGET_DIR=$SCRATCH/target
( cd ~/projects/intron_projector && cargo build --release --locked )
export BIN=$SCRATCH/target/release/intronfate RAYON_NUM_THREADS=1 W H
python3 -I "$HERE/analyse.py" --selftest
mkdir -p "$SCRATCH/j3b" "$HERE/results"; export OUT=$SCRATCH/j3b
cd "$H/V6_orthology_nine_species"
fit() {  # EVIDENCE NAME MAXIT TOL PTOL
  ev=$1; nm=$2; d=$OUT/${ev}__$nm; mkdir -p "$d"
  args=(--matrix "$W/run/combined/matrix.tsv" -t species.nwk -m likelihood --ml-origins --ml-rate-categories 4 --ml-max-iter "$3" --ml-tol "$4" --ml-param-tol "$5" -o "$d")
  [ "$ev" = soft ] && args+=(--ml-soft-evidence "$W/diag_soft.tsv")
  "$BIN" score "${args[@]}" > "$d/stdout" 2> "$d/stderr" || echo "FAILED $ev $nm" >> "$OUT/failed"
}
export -f fit
{ for ev in hard soft; do
    for it in 500 1000 2000 5000 10000 20000 50000; do echo "$ev cap$it $it 1e-6 1e-7"; done
    echo "$ev tol1e-5 50000 1e-5 1e-6"; echo "$ev tol1e-4 50000 1e-4 1e-5"; echo "$ev tol1e-3 50000 1e-3 1e-4"
  done; } | xargs -P 16 -L 1 bash -c 'fit $0 $1 $2 $3 $4'
python3 -I "$HERE/analyse.py" "$OUT" > "$HERE/results/creep.tsv"
cat "$HERE/results/creep.tsv"; [ -f "$OUT/failed" ] && cat "$OUT/failed" || true
