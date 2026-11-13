#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 32gb --time 12:00:00
#SBATCH --job-name=V11_au
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests/logs/V11_au.%j.log
# Tree tests over the per-hypothesis ML trees: IQ-TREE -z (RELL: bp-RELL, KH, SH, wKH, wSH, c-ELW) and -au; then a RAxML-NG likelihood cross-check.
# RAxML-NG has no AU/KH/SH test, so the tests are IQ-TREE's; RAxML-NG re-evaluates the same trees as an independent check of the likelihood ranking.
set -euo pipefail
module load iqtree/2.2.2.6
cd /rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests
: > all_trees.nwk; : > all_trees.names
for f in hyp/*.treefile; do cat $f >> all_trees.nwk; basename $f .treefile >> all_trees.names; done
iqtree2 -s supermatrix.fasta -m LG+F+G4 -z all_trees.nwk -zb 10000 -zw -au -n 0 -T $SLURM_CPUS_PER_TASK --prefix au -seed 1 -redo
module purge; module load raxml-ng/2.0.2
mkdir -p rx
raxml-ng --evaluate --msa supermatrix.fasta --model LG+G4+F --tree all_trees.nwk --prefix rx/eval --threads $SLURM_CPUS_PER_TASK --seed 1 --redo
