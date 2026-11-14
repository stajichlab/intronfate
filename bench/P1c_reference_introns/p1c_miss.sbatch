#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 48gb --time 02:00:00 -p stajichlab
#SBATCH --job-name=P1c_miss
#SBATCH --output=P1c_miss.%j.log
set -euo pipefail
python3 p1c_miss_profile.py --selftest
python3 p1c_miss_profile.py .
