#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 48gb --time 04:00:00 -p stajichlab
#SBATCH --job-name=P1d_audit
#SBATCH --output=P1d_audit.%j.log
set -euo pipefail
python3 p1d_version_audit.py --selftest
python3 p1d_version_audit.py version_audit.tsv 8
