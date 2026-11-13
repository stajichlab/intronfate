#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 32gb --time 06:00:00
#SBATCH --job-name=V11_prep
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests/logs/V11_prep.%j.log
# Per-site rates on the unconstrained ML tree (IQ-TREE --rate), then the filtered / recoded alignments in variants/.
set -euo pipefail
module load iqtree/2.2.2.6
cd /rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests
mkdir -p variants
iqtree2 -s supermatrix.fasta -m LG+F+G4 -te hyp/free.treefile --rate -T $SLURM_CPUS_PER_TASK --prefix variants/rates -seed 1 -redo
python3 site_filter.py drop-fast supermatrix.fasta variants/rates.rate 0.25 variants/fast25.fasta
python3 site_filter.py drop-fast supermatrix.fasta variants/rates.rate 0.50 variants/fast50.fasta
python3 site_filter.py recode supermatrix.fasta variants/dayhoff6.fasta
ls -l variants
