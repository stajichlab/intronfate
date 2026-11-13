#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 48gb --time 02:00:00
#SBATCH --job-name=V12_composition
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/V12_antarctolithica_composition/logs/V12_composition.%j.log
# Needs V10's selected_genomes.tsv and pfam_counts.tsv on the HPCC. Absolute paths.
set -euo pipefail
module load gcc/12.2.0 duckdb/1.1.3
cd /rhome/jstajich/projects/intron_projector/bench/V12_antarctolithica_composition
duckdb < gene_counts.sql
python3 composition_analysis.py | tee out/composition.log
