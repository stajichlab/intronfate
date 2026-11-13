#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 96gb --time 06:00:00 -p stajichlab
#SBATCH --job-name=P1_marker
#SBATCH --output=logs/P1_marker.%j.log
# P1 pilot (b2): marker-gene intron density. Submit from bench/P1_intron_density/. Purgeable.
set -euo pipefail
mkdir -p logs
python3 p1_marker_density.py --selftest
python3 p1_marker_density.py intron_density_per_genome.tsv .
