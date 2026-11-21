#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 32gb --time 01:00:00
#SBATCH -p stajichlab
#SBATCH --job-name=J2_support
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/J2_support_vs_reference/logs/J2.%j.log
# Read support vs RefSeq agreement for the six genomes with a J1 junction BED. Needs bench/J1_junction_support/work/*.bed.
set -euo pipefail
cd ~/projects/intron_projector/bench/J2_support_vs_reference
mkdir -p logs results
python3 -I support_vs_reference.py --selftest
python3 -I support_vs_reference.py results/J2
