#!/usr/bin/env bash
# Link each genome's GFF3 and protein FASTA (read-only, from the BFD run) into data/<label>/ and verify them.
set -euo pipefail
cd "$(dirname "$0")"; . ./config.sh
mkdir -p data
printf 'label\tgff3\tproteins\tn_proteins\tn_cds_features\n' > data/manifest.tsv
grep -v '^#' genomes.tsv | while IFS=$'\t' read -r label dir group; do
    src=$BFD/$dir/predict_results
    gff=$src/$dir.gff3.gz; prot=$src/$dir.proteins.fa.gz
    [ -s "$gff" ] && [ -s "$prot" ] || { echo "missing files for $label ($src)" >&2; exit 1; }
    mkdir -p "data/$label"; ln -sf "$gff" "data/$label/$label.gff3.gz"; ln -sf "$prot" "data/$label/$label.proteins.fa.gz"
    np=$(zcat "$prot" | grep -c '^>'); nc=$(zcat "$gff" | awk -F'\t' '$3=="CDS"' | wc -l)
    printf '%s\t%s\t%s\t%s\t%s\n' "$label" "$gff" "$prot" "$np" "$nc" >> data/manifest.tsv
done
cat data/manifest.tsv
