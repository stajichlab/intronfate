#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 1 --mem 4gb --time 00:20:00
#SBATCH -p highclock
#SBATCH --job-name=F2_missing
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/F2_column_ambiguity/logs/F2.%j.log
# F2: how much of the missing data is a short alignment gap. Reads the V6 tier A matrix (flank 1) and alignments; no heavy output.
# Submit: mkdir -p logs results && sbatch run.sbatch
set -euo pipefail
module load gcc/12.2.0 >/dev/null 2>&1 || true
F2=~/projects/intron_projector/bench/F2_column_ambiguity
V6=~/projects/intron_projector/bench/V6_orthology_nine_species
cd "$V6"
python3 -I "$F2/missing_cells.py" --matrix results_projection/A_strict/combined/matrix.tsv --aln-dir alignments/A_strict \
   --taxon-map families/A_strict/taxon_map.tsv --flank 1 | tee "$F2/results/missing_cells.txt"
