#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 48gb --time 06:00:00
#SBATCH -p highclock
#SBATCH --job-name=J1_junctions
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/J1_junction_support/logs/J1.%j.log
# J1: RNA-seq junction support as soft evidence on the V6 nine-species tier A families.
# Junctions come from the funannotate training alignments (hisat2 BAMs of the normalised RNA-seq) that BFD already holds.
# Large files go to $SCRATCH; only the small tables are kept in results/.
# Submit: mkdir -p logs results && sbatch run.sbatch      (ONLY_PREP=1 stops after the projection; MAXITER=n caps the EM)      (absolute paths, not dirname "$0")
set -euo pipefail
: "${SCRATCH:?run inside a Slurm job}"
J1=~/projects/intron_projector/bench/J1_junction_support
V6=~/projects/intron_projector/bench/V6_orthology_nine_species
REPO=~/projects/intron_projector
GA=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation
module load rust >/dev/null 2>&1 || true
module load gcc/12.2.0 >/dev/null 2>&1 || true
module load samtools >/dev/null 2>&1 || true
find "$REPO/crates" -name '*.rs' -exec touch {} +
export CARGO_TARGET_DIR=$SCRATCH/target
( cd "$REPO" && cargo build --release --locked )
BIN=$SCRATCH/target/release/intronfate
python3 -I "$REPO/bench/test_junction_scripts.py" "$BIN"          # the fixture test first (coordinates end to end)

# $J1/work is the purgeable cache (the tier A projection and the tables); delete it to start over
W=$J1/work; mkdir -p "$W" "$J1/results"
cd "$V6"
# 1. junctions from the BAMs that exist
jargs=(); : > "$J1/results/junction_inputs.tsv"; echo -e "label\tbam\tjunctions\ttotal_reads" > "$J1/results/junction_inputs.tsv"
for l in $(grep -v '^#' genomes.tsv | cut -f1); do
  d=$(grep -v '^#' genomes.tsv | awk -v l=$l '$1==l{print $2}')
  bam=$GA/$d/training/funannotate_train.coordSorted.bam
  [ -s "$bam" ] || { echo "no BAM for $l"; continue; }
  samtools view "$bam" | python3 -I "$J1/bam_junctions.py" > "$W/$l.bed"
  echo -e "$l\t$bam\t$(wc -l < $W/$l.bed)\t$(awk '{s+=$5} END{print s+0}' $W/$l.bed)" >> "$J1/results/junction_inputs.tsv"
  jargs+=(--junctions "$l=$W/$l.bed")
done
cat "$J1/results/junction_inputs.tsv"
# 2. the tier A projection with one folder per family (sites.tsv) and the combined matrix
tier=A_strict
gffs=(); prots=(); for l in $(grep -v '^#' genomes.tsv | cut -f1); do gffs+=(--gff "data/$l/$l.gff3.gz"); prots+=(-p "data/$l/$l.proteins.fa.gz"); done
if [ ! -s "$W/run/combined/matrix.tsv" ]; then
  rm -rf "$W/run"
  "$BIN" analyze "${gffs[@]}" "${prots[@]}" --msa-dir alignments/$tier --msa-ext aln.fasta --combine \
     --taxon-map families/$tier/taxon_map.tsv -t species.nwk -m counts -o "$W/run" > /dev/null
fi
[ -z "${ONLY_PREP:-}" ] || { echo "prep done"; exit 0; }
rm -rf "$W"/score_* "$W"/soft_* "$W"/junctions_*
# 3. hard fit, then soft fits over a grid of (p0, k)
SC=(--matrix "$W/run/combined/matrix.tsv" -t species.nwk -m likelihood --ml-origins --ml-rate-categories 4 --ml-max-iter "${MAXITER:-50000}")
"$BIN" score "${SC[@]}" -o "$W/score_hard" > "$W/score_hard.stdout"
for pk in "0.8 3" "0.5 3" "0.95 3" "0.8 1" "0.8 10"; do
  set -- $pk; p=$1; k=$2
  "$BIN" junctions --introns "$W/run/introns.tsv" --sites-dir "$W/run" --matrix "$W/run/combined/matrix.tsv" --combined-labels --junction-format bed \
     "${jargs[@]}" --unsupported-p "$p" --half-reads "$k" -o "$W/soft_${p}_${k}.tsv" 2> "$W/junctions_${p}_${k}.stderr"
  "$BIN" score "${SC[@]}" --ml-soft-evidence "$W/soft_${p}_${k}.tsv" -o "$W/score_${p}_${k}" > "$W/score_${p}_${k}.stdout"
done
python3 -I "$J1/summarize.py" "$W" | tee "$J1/results/summary.txt"
cp "$W"/junctions_*.stderr "$J1/results/" 
