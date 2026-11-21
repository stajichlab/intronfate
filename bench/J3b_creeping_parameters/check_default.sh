#!/usr/bin/bash -l
# Checks that the default tolerance for a 4-class fit equals the explicit `--ml-tol 1e-5 --ml-param-tol 1e-6` of J3b (same iterations and ln L),
# and that a constant-rate fit keeps the old default. Run inside a job: srun -p highclock -c 4 bash check_default.sh
set -euo pipefail
: "${SCRATCH:?run inside a Slurm job}"
H=~/projects/intron_projector/bench; W=$H/J1_junction_support/work
module load rust gcc/12.2.0 >/dev/null 2>&1 || true
export CARGO_TARGET_DIR=$SCRATCH/target; ( cd ~/projects/intron_projector && cargo build --release --locked 2>&1 | tail -1 )
BIN=$SCRATCH/target/release/intronfate; cd "$H/V6_orthology_nine_species"
run() { "$BIN" score --matrix "$W/run/combined/matrix.tsv" -t species.nwk -m likelihood "${@:2}" -o "$SCRATCH/chk_$1" > /dev/null 2>&1 || true
        echo "$1: $(awk -F'\t' '$1=="ml_log_likelihood"||$1=="ml_iterations"||$1=="ml_converged"{printf "%s=%s ", $1, $2}' "$SCRATCH/chk_$1/summary_score.tsv")"; }
run default4 --ml-rate-categories 4
run explicit4 --ml-rate-categories 4 --ml-tol 1e-5 --ml-param-tol 1e-6
run old4 --ml-rate-categories 4 --ml-tol 1e-6 --ml-param-tol 1e-7
run default1
run old1 --ml-tol 1e-6 --ml-param-tol 1e-7
