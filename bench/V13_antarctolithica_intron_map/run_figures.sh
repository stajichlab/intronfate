#!/usr/bin/env bash
# V13 figures with the pixi environment (numpy, matplotlib). Run on the HPCC after run_v13.sbatch:
#   srun -p short -c 2 --mem 8gb --time 0:30:00 bash run_figures.sh
set -euo pipefail
HERE=~/projects/intron_projector/bench/V13_antarctolithica_intron_map
V9=~/projects/intron_projector/bench/V9_antarctolithica_expanded
cd "$HERE"
pixi run python -I v13_figures.py --trees "$V9/trees" --genomes "$V9/genomes.tsv" --score score --results results --out figures
