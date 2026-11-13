#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 32gb --time 24:00:00
#SBATCH --job-name=V11_varau
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests/logs/V11_varau_%x.%j.log
# Tree tests (RELL tests + AU) over the per-hypothesis trees of one variant.
set -euo pipefail
module load iqtree/2.2.2.6
cd /rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests
source variant_settings.sh
: > var_$VARIANT/all_trees.nwk; : > var_$VARIANT/all_trees.names
for f in var_$VARIANT/*.treefile; do cat $f >> var_$VARIANT/all_trees.nwk; basename $f .treefile >> var_$VARIANT/all_trees.names; done
iqtree2 -s $ALN -m $MODEL $EXTRA -z var_$VARIANT/all_trees.nwk -zb 10000 -zw -au -n 0 -T $SLURM_CPUS_PER_TASK --prefix var_$VARIANT/au -seed 1 -redo
