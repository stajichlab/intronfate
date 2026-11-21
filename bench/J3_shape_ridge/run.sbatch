#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 24gb --time 04:00:00
#SBATCH -p highclock
#SBATCH --job-name=J3_ridge
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/J3_shape_ridge/logs/J3.%j.log
# Profile of the log-likelihood over the gamma shape (4 classes) on the J1 combined matrix, hard and soft (p0 0.8, k 3).
# Each shape is FIXED (--ml-gamma-shape), so the flat ridge cannot slow the fit; plus the free-shape fit and constant rates.
# Needs bench/J1_junction_support/work (run J1 with ONLY_PREP=1 first); work/diag_soft.tsv comes from diag_soft_convergence.sbatch.
# Output: results/profile.tsv (one row per fit), work under $SCRATCH.
set -euo pipefail
: "${SCRATCH:?run inside a Slurm job}"
H=~/projects/intron_projector/bench
W=$H/J1_junction_support/work
BIN=$SCRATCH/target/release/intronfate
if [ ! -x "$BIN" ]; then module load rust gcc/12.2.0 >/dev/null 2>&1; export CARGO_TARGET_DIR=$SCRATCH/target; ( cd ~/projects/intron_projector && cargo build --release --locked ); fi
mkdir -p "$H/J3_shape_ridge/results" "$SCRATCH/j3"
export RAYON_NUM_THREADS=1
MAXIT=${MAXIT:-20000}
cd "$H/V6_orthology_nine_species"
fit() {  # fit EVIDENCE SHAPE   (SHAPE: a number, "free", or "const")
  ev=$1; sh=$2; d=$SCRATCH/j3/${ev}_${sh}; mkdir -p "$d"
  args=(--matrix "$W/run/combined/matrix.tsv" -t species.nwk -m likelihood --ml-max-iter "$MAXIT" -o "$d")
  case $sh in const) ;; free) args+=(--ml-rate-categories 4) ;; *) args+=(--ml-rate-categories 4 --ml-gamma-shape "$sh") ;; esac
  [ "$ev" = soft ] && args+=(--ml-soft-evidence "$W/diag_soft.tsv")
  "$BIN" score "${args[@]}" > "$d/stdout" 2> "$d/stderr" || echo "FAILED $ev $sh" >> "$SCRATCH/j3/failed"
}
export -f fit; export W BIN SCRATCH MAXIT
{ for ev in hard soft; do for sh in 0.5 1 2 3 4 5 6 8 12 20 50 free const; do echo "$ev $sh"; done; done; } | xargs -P 8 -L 1 bash -c 'fit $0 $1'
printf 'evidence\tshape\tlnL\tconverged\titerations\tfitted_shape\n' > "$H/J3_shape_ridge/results/profile.tsv"
for d in "$SCRATCH"/j3/*_*/; do
  n=$(basename "$d"); ev=${n%%_*}; sh=${n#*_}; s=$d/summary_score.tsv
  [ -s "$s" ] || { echo "no summary for $n" >&2; continue; }
  g() { awk -F'\t' -v k=$1 '$1==k{print $2}' "$s"; }
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$ev" "$sh" "$(g ml_log_likelihood)" "$(g ml_converged)" "$(g ml_iterations)" "$(g ml_gamma_shape)" >> "$H/J3_shape_ridge/results/profile.tsv"
done
sort -k1,1 -k2,2V "$H/J3_shape_ridge/results/profile.tsv" -o "$H/J3_shape_ridge/results/profile.tsv.sorted" && head -1 "$H/J3_shape_ridge/results/profile.tsv" > "$H/J3_shape_ridge/results/p.tmp" && grep -v '^evidence' "$H/J3_shape_ridge/results/profile.tsv.sorted" >> "$H/J3_shape_ridge/results/p.tmp" && mv "$H/J3_shape_ridge/results/p.tmp" "$H/J3_shape_ridge/results/profile.tsv" && rm "$H/J3_shape_ridge/results/profile.tsv.sorted"
cat "$H/J3_shape_ridge/results/profile.tsv"; [ -f "$SCRATCH/j3/failed" ] && cat "$SCRATCH/j3/failed" || true
