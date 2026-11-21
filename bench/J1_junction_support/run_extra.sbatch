#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 32gb --time 02:00:00
#SBATCH -p highclock
#SBATCH --job-name=J1_extra
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/J1_junction_support/logs/J1x.%j.log
# J1 follow-up: more (p0, k) settings on the cache that run.sbatch left in work/ (BEDs, tier A projection, hard fit).
# Review of J2 (Fable) asked for k = 0.2, the value the one-read agreement in J2 points to. Constant rates (converged, fast).
# Submit: sbatch run_extra.sbatch        EXTRA="0.8 0.2;0.8 0.5" TAG=rc1_ RATECAT=1 (defaults)
set -euo pipefail
: "${SCRATCH:?run inside a Slurm job}"
J1=~/projects/intron_projector/bench/J1_junction_support
V6=~/projects/intron_projector/bench/V6_orthology_nine_species
REPO=~/projects/intron_projector
W=$J1/work
module load rust >/dev/null 2>&1 || true
module load gcc/12.2.0 >/dev/null 2>&1 || true
find "$REPO/crates" -name '*.rs' -exec touch {} +
export CARGO_TARGET_DIR=$SCRATCH/target
( cd "$REPO" && cargo build --release --locked )
BIN=$SCRATCH/target/release/intronfate
TAG=${TAG:-rc1_}; RATECAT=${RATECAT:-1}
[ -s "$W/run/combined/matrix.tsv" ] && [ -d "$W/score_${TAG}hard" ] || { echo "run run.sbatch first (TAG=$TAG RATECAT=$RATECAT)"; exit 1; }
cd "$V6"
jargs=(); for l in $(grep -v '^#' genomes.tsv | cut -f1); do [ -s "$W/$l.bed" ] && jargs+=(--junctions "$l=$W/$l.bed"); done
SC=(--matrix "$W/run/combined/matrix.tsv" -t species.nwk -m likelihood --ml-origins --ml-rate-categories "$RATECAT" --ml-max-iter "${MAXITER:-50000}")
IFS=';' read -ra PKS <<< "${EXTRA:-0.8 0.2;0.8 0.5}"
for pk in "${PKS[@]}"; do
  set -- $pk; p=$1; k=$2
  "$BIN" junctions --introns "$W/run/introns.tsv" --sites-dir "$W/run" --matrix "$W/run/combined/matrix.tsv" --combined-labels --junction-format bed \
     "${jargs[@]}" --unsupported-p "$p" --half-reads "$k" -o "$W/soft_${p}_${k}.tsv" 2> "$W/junctions_${p}_${k}.stderr"
  "$BIN" score "${SC[@]}" --ml-soft-evidence "$W/soft_${p}_${k}.tsv" -o "$W/score_${TAG}${p}_${k}" > "$W/score_${TAG}${p}_${k}.stdout"
done
python3 -I "$J1/summarize.py" "$W" "$TAG" | tee "$J1/results/summary${TAG:+.$TAG}.txt"
cp "$W"/junctions_*.stderr "$J1/results/"
