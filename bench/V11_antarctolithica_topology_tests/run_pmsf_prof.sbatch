#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 32 --mem 320gb --time 48:00:00
#SBATCH --job-name=V11_pmsf_prof
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests/logs/V11_pmsf_prof.%j.log
# PMSF step 1: fit LG+C60+F+G4 on the guide tree (the unconstrained LG+F+G4 ML tree) once and write the site-specific frequency profile
# (IQ-TREE asks for ~245 GB for 200,000 sites; the first attempt with 32 GB per hypothesis task ran out of memory in seconds).
set -euo pipefail
module load iqtree/2.2.2.6
cd /rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests
mkdir -p pmsf
iqtree2 -s supermatrix.fasta -m LG+C60+F+G4 -ft hyp/free.treefile -n 0 -T $SLURM_CPUS_PER_TASK --prefix pmsf/prof -seed 1 -redo
ls -l pmsf
