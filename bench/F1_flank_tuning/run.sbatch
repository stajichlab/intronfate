#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 16gb --time 01:00:00
#SBATCH -p highclock
#SBATCH --job-name=F1_flank
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/F1_flank_tuning/logs/F1.%j.log
# F1: --flank sweep on the V6 nine-species tier A alignments. Outputs go to $SCRATCH; only the small tables are kept in results/.
# Submit: mkdir -p logs results && sbatch run.sbatch      (absolute paths, not dirname "$0")
set -euo pipefail
: "${SCRATCH:?run inside a Slurm job}"
F1=~/projects/intron_projector/bench/F1_flank_tuning
V6=~/projects/intron_projector/bench/V6_orthology_nine_species
REPO=~/projects/intron_projector
module load rust >/dev/null 2>&1 || true
module load gcc/12.2.0 >/dev/null 2>&1 || true
find "$REPO/crates" -name '*.rs' -exec touch {} +
export CARGO_TARGET_DIR=$SCRATCH/target
( cd "$REPO" && cargo build --release --locked )
BIN=$SCRATCH/target/release/intronfate
cd "$V6"
tier=A_strict
gffs=(); prots=(); for l in $(grep -v '^#' genomes.tsv | cut -f1); do gffs+=(--gff "data/$l/$l.gff3.gz"); prots+=(-p "data/$l/$l.proteins.fa.gz"); done
mkdir -p "$F1/results"
for fl in ${FLANKS:-0 1 2 3 5 8}; do
  out=$SCRATCH/flank$fl; rm -rf "$out"
  "$BIN" analyze "${gffs[@]}" "${prots[@]}" --msa-dir alignments/$tier --msa-ext aln.fasta --combine --combined-only \
     --taxon-map families/$tier/taxon_map.tsv -t species.nwk -m dollo,counts --flank $fl -o "$out" > /dev/null
  cp "$out/combined/matrix.tsv" "$SCRATCH/matrix.flank$fl.tsv"
  cp "$out/combined/dollo_branches.tsv" "$F1/results/dollo_branches.flank$fl.tsv"
  head -1 "$out/combined/matrix.tsv" > "$F1/results/matrix_header.flank$fl.txt"
done
python3 -I "$F1/gap_distance.py" --matrix "$SCRATCH/matrix.flank0.tsv" --aln-dir alignments/$tier --aln-ext aln.fasta \
   --taxon-map families/$tier/taxon_map.tsv --out "$F1/results/absence_by_gap_distance.tsv"
python3 -I "$F1/summarize.py" "$F1/results" > "$F1/results/summary.txt"
cat "$F1/results/summary.txt"
