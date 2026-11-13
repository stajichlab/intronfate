#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 32gb --time 24:00:00
#SBATCH --array=0-8
#SBATCH --job-name=V11_hyp
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests/logs/V11_hyp.%A_%a.log
# One ML search per hypothesis (task 0 = unconstrained), IQ-TREE LG+F+G4 on the supermatrix, constrained by hypotheses/<H>.constraint.nwk.
set -euo pipefail
module load iqtree/2.2.2.6
cd /rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests
mkdir -p hyp
HYPS=(free $(ls hypotheses/*.constraint.nwk | xargs -n1 basename | sed 's/.constraint.nwk//'))
H=${HYPS[$SLURM_ARRAY_TASK_ID]}
if [ "$H" = free ]; then CON=""; else CON="-g hypotheses/$H.constraint.nwk"; fi
iqtree2 -s supermatrix.fasta -m LG+F+G4 $CON -T $SLURM_CPUS_PER_TASK --prefix hyp/$H -seed 1 -redo
