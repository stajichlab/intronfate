#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 8gb --time 00:20:00
#SBATCH --job-name=V6_check_tutorial
#SBATCH --output=logs/V6_check_tutorial.%j.log

# Runs the section 6 command of docs/TUTORIAL.md exactly as written (list for the GFFs, proteomes folder, alignments folder, --combine --combined-only)
# on the V6 tier A and compares the result with the earlier run. Output in $SCRATCH. Absolute path, not dirname "$0".
set -euo pipefail
: "${SCRATCH:?run inside a Slurm job}"
HERE=~/projects/intron_projector/bench/V6_orthology_nine_species
cd "$HERE"
export PATH=$HOME/projects/intron_projector/target/release:$PATH
ls data/*/*.gff3.gz > "$SCRATCH/gffs.txt"; cp "$SCRATCH/gffs.txt" gffs.txt       # the list lives in the working folder, as in the tutorial
rm -rf "$SCRATCH/tut"
intronfate analyze \
    --gff-list gffs.txt --proteins-dir proteomes --proteins-ext faa \
    --msa-dir alignments/A_strict --msa-ext aln.fasta \
    --taxon-map families/A_strict/taxon_map.tsv \
    -t species.nwk -m dollo,mp,counts \
    --combine --combined-only -o "$SCRATCH/tut" > /dev/null
rm -f gffs.txt
for f in matrix.tsv dollo_branches.tsv dollo_sites.tsv; do
    cmp "results_projection/A_strict/combined/$f" "$SCRATCH/tut/combined/$f" && echo "identical to the earlier run: combined/$f"
done
cat "$SCRATCH/tut/combined/dollo_branches.tsv" | cut -f1-5 | head -20
