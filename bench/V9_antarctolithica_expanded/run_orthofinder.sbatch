#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 48gb --time 06:00:00
#SBATCH --job-name=V9_orthofinder
#SBATCH --output=logs/V9_orthofinder.%j.log

# OrthoFinder 3.1.0 on the proteomes in genomes.tsv with MSA-based gene trees; no species tree is given, so OrthoFinder infers one from the genes (STAG) (default: famsa + FastTree).
# VARIANT=plain (default) or VARIANT=y (adds -y: split paralogous clades below the root of a HOG into separate HOGs).
# Submit from the bench folder:  mkdir -p logs && sbatch -p highclock --export=ALL,VARIANT=plain run_orthofinder.sbatch
# Heavy temporary files live in $SCRATCH; only the result tables and trees are copied back to results/<variant>/.
# NB: absolute path, not dirname "$0": under Slurm $0 is a copy in the spool directory.
set -euo pipefail
HERE=~/projects/intron_projector/bench/V9_antarctolithica_expanded
cd "$HERE"
: "${SCRATCH:?run inside a Slurm job}"
VARIANT=${VARIANT:-plain}
module load orthofinder/3.1.0
OF=/opt/linux/rocky/8.x/x86_64/pkgs/orthofinder/3.1.0/env/bin/orthofinder

./prepare_proteomes.sh
rm -rf "$SCRATCH/of_in" "$SCRATCH/of_out"; mkdir -p "$SCRATCH/of_in"
cp proteomes/*.faa "$SCRATCH/of_in/"
extra=(); [ "$VARIANT" = y ] && extra=(-y)
date; time "$OF" -f "$SCRATCH/of_in" -t "${SLURM_CPUS_PER_TASK:-16}" -a 4 -o "$SCRATCH/of_out" -n "$VARIANT" "${extra[@]}"
date
res=$(ls -d "$SCRATCH"/of_out/Results_"$VARIANT")
mkdir -p "results/$VARIANT"
# Keep what later steps read; skip WorkingDirectory (all-vs-all search results, MSAs): large and reproducible.
rsync -a --exclude WorkingDirectory "$res"/ "results/$VARIANT/"
du -sh "results/$VARIANT"; ls "results/$VARIANT"
