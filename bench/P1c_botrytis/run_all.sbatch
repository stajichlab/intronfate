#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 64gb --time 04:00:00 -p stajichlab
#SBATCH --job-name=P1c_bot
#SBATCH --output=P1c_bot.%j.log
set -euo pipefail
python3 p1c_ref_introns.py .
python3 p1c_miss_profile.py .
python3 p1c_markers.py ../P1_intron_density/marker_pfams.tsv .
