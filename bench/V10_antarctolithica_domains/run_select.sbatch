#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 48gb --time 03:00:00
#SBATCH --job-name=V10_pfam_counts
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/V10_antarctolithica_domains/logs/V10_pfam_counts.%j.log
# Select genomes, then count Pfam domains per genome (domain i-Evalue < 1e-5) from BFD tables/pfam.parquet. Absolute paths (no dirname $0).
set -euo pipefail
module load gcc/12.2.0 duckdb/1.1.3
cd /rhome/jstajich/projects/intron_projector/bench/V10_antarctolithica_domains
python3 select_genomes.py ../V9_antarctolithica_expanded selected_genomes.tsv
duckdb < pfam_counts.sql
wc -l pfam_counts.tsv
