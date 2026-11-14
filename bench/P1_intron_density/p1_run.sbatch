#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 48gb --time 06:00:00 -p stajichlab
#SBATCH --job-name=P1_intron_density
#SBATCH --output=logs/P1_intron_density.%j.log

# P1 pilot (a): per-genome intron density from the BFD tables (read-only).
# Submit from ~/projects/intron_projector/bench/P1_intron_density/ . Purgeable.
set -euo pipefail
mkdir -p logs
python3 p1_intron_density.py intron_density_per_genome.tsv
