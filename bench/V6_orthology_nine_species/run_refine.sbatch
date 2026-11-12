#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 16gb --time 02:00:00
#SBATCH --job-name=V6_orthosynassign
#SBATCH --output=logs/V6_orthosynassign.%j.log

# Synteny refinement of OrthoFinder orthogroups with orthoSynAssign (pixi env in this folder), at several window sizes.
# Submit:  mkdir -p logs && sbatch -p highclock --export=ALL,VARIANT=plain run_refine.sbatch
# NB: absolute path, not dirname "$0" (under Slurm $0 is a copy in the spool directory).
set -euo pipefail
export PATH=$HOME/.pixi/bin:$PATH
HERE=~/projects/intron_projector/bench/V6_orthology_nine_species
cd "$HERE"
VARIANT=${VARIANT:-plain}
D=$HERE/data
R=results/$VARIANT
mkdir -p bed "refined/$VARIANT"
# BED per species (labels equal the OrthoFinder species names), from the funannotate GFF3s
grep -v '^#' genomes.tsv | cut -f1 | while read -r l; do
    python3 make_bed.py "$D/$l/$l.gff3.gz" > "bed/$l.bed"
done
wc -l bed/*.bed | tail -1
for w in ${WINDOWS:-8 30}; do
    pixi run orthosynassign --og_file "$R/Orthogroups/Orthogroups.tsv" --bed bed/*.bed -w "$w" \
        -t "${SLURM_CPUS_PER_TASK:-8}" -o "refined/$VARIANT/Refined_SOGs.w$w.tsv" 2>&1 | grep -v -E 'Progress:|WARN' | tail -3
done
python3 orthogroup_summary.py "$R" refined/$VARIANT/Refined_SOGs.w*.tsv | tee "refined/$VARIANT/summary.txt"
