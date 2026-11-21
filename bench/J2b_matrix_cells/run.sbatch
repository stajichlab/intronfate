#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 32gb --time 01:00:00
#SBATCH -p stajichlab
#SBATCH --job-name=J2b_cells
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/J2b_matrix_cells/logs/J2b.%j.log
# J2 on the matrix cells with an exact / near3 / near12 RefSeq match. Needs bench/J1_junction_support/work (run.sbatch there).
set -euo pipefail
cd ~/projects/intron_projector/bench/J2b_matrix_cells
mkdir -p logs results
J1=../J1_junction_support/work
python3 -I matrix_cells.py --selftest
python3 -I matrix_cells.py $J1 $J1/run results/J2b
