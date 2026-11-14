#!/bin/bash
#SBATCH -p epyc
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:20:00
#SBATCH -J P7_lineage
#SBATCH -o P7_lineage.%j.log
cd "$SLURM_SUBMIT_DIR"
python3 lineage_summary.py --selftest && python3 lineage_summary.py
