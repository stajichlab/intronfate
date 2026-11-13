#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 32gb --time 02:00:00 -p stajichlab
#SBATCH --job-name=P1e_classify
#SBATCH --output=P1e_classify.%j.log
set -euo pipefail
python3 p1e_classify.py --selftest
python3 p1e_classify.py . GCF_000011425.1_ASM1142v1 GCF_000219625.1_MYCGR_v2.0 GCF_000146045.2_R64 GCF_001661235.1_Picme2 GCF_000151645.1_ASM15164v1 GCF_000143535.2_ASM14353v4 GCF_000182925.2_NC12 GCF_000149245.1_CNA3 GCF_000002945.2_ASM294v3
