#!/usr/bin/env bash
# What are the unmatched / no_cds records of the how-to run? Re-downloads the fungi5 genomes into $SCRATCH and runs `intronfate project`.
#   srun -p short -c 2 --mem 8gb --time 0:20:00 bash check_unmatched.sh
set -euo pipefail
: "${SCRATCH:?run inside a Slurm job}"
HERE=~/projects/intron_projector/bench/D1_howto_public_genomes; REPO=$HERE/repo
cd "$SCRATCH"; mkdir -p genomes
grep -v '^#' "$REPO/examples/04_public_genomes/fungi5.tsv" | while IFS=$'\t' read -r label asm organism; do
  url=https://ftp.ncbi.nlm.nih.gov/genomes/all/${asm:0:3}/${asm:4:3}/${asm:7:3}/${asm:10:3}/$asm/$asm
  curl -fsSL -o genomes/$label.gff.gz ${url}_genomic.gff.gz; curl -fsSL -o genomes/$label.protein.faa.gz ${url}_protein.faa.gz
done
"$REPO/target/release/intronfate" project --gff-dir genomes --gff-ext gff --proteins-dir genomes --proteins-ext protein.faa -o p > /dev/null 2>&1
echo "== unmatched.tsv header and first rows"; head -4 p/unmatched.tsv
echo "== unmatched by source file and reason"; python3 -I - p/unmatched.tsv <<'PY'
import csv, sys, collections
r = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
print(collections.Counter(tuple(v for k, v in row.items() if k in ("source", "file", "reason", "kind")) for row in r).most_common(10))
PY
echo "== no_cds by feature type"; grep no_cds p/diagnostics_detail.tsv | wc -l
: # no_cds rows are tRNA, rRNA and ncRNA transcripts (checked 2026-11-13)
echo "== the CDS rows of unmatched transcripts (first per genome)"
for g in Anidulans Spombe Scerevisiae; do
  t=$(awk -F'\t' -v s="genomes/$g.gff.gz" '$5==s {print $1; exit}' p/unmatched.tsv)
  echo "-- $g $t"; gzip -dc genomes/$g.gff.gz | awk -F'\t' -v t="Parent=$t;" 'index($9,t) {print $3"\t"$9}' | head -3 | cut -c1-300 || true
done
