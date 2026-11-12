#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 1 --mem 8gb --time 00:30:00
#SBATCH --job-name=V4_busco_accuracy
#SBATCH --output=logs/V4_busco_accuracy.%j.log

# V4: intron accuracy over the genes actually used in the V2 alignments (read-only on the shared scoring folder).
# Submit from ~/projects/intron_projector/bench/V4_annotation_accuracy/ :  mkdir -p logs && sbatch -p highclock run.sbatch
set -euo pipefail
SCORING=${SCORING:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/rc6_vs_pre_rc6}
python3 busco_gene_accuracy.py "$SCORING" ../V2_gff_projection | tee results_report.txt
