#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 16gb --time 02:00:00
#SBATCH --job-name=V8_screen_density
#SBATCH --output=logs/V8_screen_density.%j.log

# Intron density (introns per protein, fraction of proteins with introns) for candidate genomes, from `intronfate project` on their funannotate GFF3.
# Candidates: screen_patterns.txt (directory-name prefixes under BFD genome_annotation). Output: screen_density.tsv. Work in $SCRATCH.
# Submit: sbatch -p highclock screen_density.sbatch       Absolute path, not dirname "$0".
set -euo pipefail
: "${SCRATCH:?run inside a Slurm job}"
HERE=~/projects/intron_projector/bench/V8_antarctolithica_placement
BFD=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation
BIN=~/projects/intron_projector/target/release/intronfate
cd "$HERE"
printf 'prefix\tgenome\tproteins\tintrons\tintrons_per_protein\tfraction_proteins_with_introns\n' > screen_density.tsv
grep -v '^#' screen_patterns.txt | while read -r pre; do
    d=""
    for c in "$BFD"/${pre}*; do
        n=$(basename "$c")
        [ -s "$c/predict_results/$n.gff3.gz" ] && [ -s "$c/predict_results/$n.proteins.fa.gz" ] && { d=$n; break; }
    done
    [ -n "$d" ] || { echo "no complete genome for $pre" >&2; continue; }
    rm -rf "$SCRATCH/scr"; mkdir -p "$SCRATCH/scr"
    "$BIN" project --gff "$BFD/$d/predict_results/$d.gff3.gz" -p "$BFD/$d/predict_results/$d.proteins.fa.gz" -o "$SCRATCH/scr" > /dev/null 2>&1 || { echo "project failed for $d" >&2; continue; }
    np=$(zcat "$BFD/$d/predict_results/$d.proteins.fa.gz" | grep -c '^>')
    awk -F'\t' -v pre="$pre" -v d="$d" -v np="$np" 'NR==1{for(i=1;i<=NF;i++) if($i=="protein_id") c=i; next} {n++; s[$c]=1} END{k=0; for(x in s) k++; printf "%s\t%s\t%d\t%d\t%.2f\t%.2f\n", pre, d, np, n, n/np, k/np}' "$SCRATCH/scr/introns.tsv" >> screen_density.tsv
done
sort -t$'\t' -k5,5 -nr screen_density.tsv | column -t -s$'\t' | cut -c1-150
