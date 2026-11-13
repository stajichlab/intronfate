#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 16gb --time 01:00:00
#SBATCH --job-name=V10_pfam_analysis
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/V10_antarctolithica_domains/logs/V10_pfam_analysis.%j.log
set -euo pipefail
module load gcc/12.2.0
cd /rhome/jstajich/projects/intron_projector/bench/V10_antarctolithica_domains
python3 pfam_analysis.py
