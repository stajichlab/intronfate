#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 64gb --time 03:00:00 -p stajichlab
#SBATCH --job-name=P1c_mark
#SBATCH --output=P1c_mark.%j.log
set -euo pipefail
python3 p1c_markers.py --selftest
python3 p1c_markers.py ../P1_intron_density/marker_pfams.tsv .
