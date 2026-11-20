#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 16gb --time 01:00:00
#SBATCH -p highclock
#SBATCH --job-name=J1_diag
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/J1_junction_support/logs/J1diag.%j.log
# Why do the soft fits stop at the iteration cap while the hard fit converges? Traces (-vv) of the hard and one soft fit, 3000 iterations.
# Needs work/run from `ONLY_PREP=1 sbatch run.sbatch`. Output: results/trace_hard.txt, results/trace_soft.txt (every 100th line + the last 20).
set -euo pipefail
: "${SCRATCH:?run inside a Slurm job}"
J1=~/projects/intron_projector/bench/J1_junction_support
V6=~/projects/intron_projector/bench/V6_orthology_nine_species
BIN=$SCRATCH/target/release/intronfate
[ -x "$BIN" ] || { module load rust gcc/12.2.0 >/dev/null 2>&1; export CARGO_TARGET_DIR=$SCRATCH/target; ( cd ~/projects/intron_projector && cargo build --release --locked ); }
W=$J1/work; cd "$V6"
SC=(--matrix "$W/run/combined/matrix.tsv" -t species.nwk -m likelihood --ml-rate-categories 4 --ml-max-iter 3000)
jargs=(); for l in Ccinerea Cneoformans Spombe Scerevisiae Ncrassa Anidulans; do jargs+=(--junctions "$l=$W/$l.bed"); done
"$BIN" junctions --introns "$W/run/introns.tsv" --sites-dir "$W/run" --matrix "$W/run/combined/matrix.tsv" --combined-labels --junction-format bed "${jargs[@]}" -o "$W/diag_soft.tsv" 2> /dev/null
trim() { awk 'NR%100==1' "$1"; echo ...; tail -20 "$1"; }
"$BIN" -vv score "${SC[@]}" -o "$SCRATCH/diag_hard" > /dev/null 2> "$SCRATCH/trace_hard.raw" || true
"$BIN" -vv score "${SC[@]}" --ml-soft-evidence "$W/diag_soft.tsv" -o "$SCRATCH/diag_soft" > /dev/null 2> "$SCRATCH/trace_soft.raw" || true
wc -l "$SCRATCH"/trace_*.raw
trim "$SCRATCH/trace_hard.raw" > "$J1/results/trace_hard.txt"
trim "$SCRATCH/trace_soft.raw" > "$J1/results/trace_soft.txt"
echo "== hard"; tail -6 "$SCRATCH/trace_hard.raw"; echo "== soft"; tail -6 "$SCRATCH/trace_soft.raw"
