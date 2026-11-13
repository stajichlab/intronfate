#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 32gb --time 12:00:00
#SBATCH --array=0-8
#SBATCH --job-name=V11_var
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests/logs/V11_var_%x.%A_%a.log
# One ML search per hypothesis for one variant (VARIANT exported by submit_variants.sh); same hypotheses as run_hyp.sbatch.
set -euo pipefail
module load iqtree/2.2.2.6
cd /rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests
source variant_settings.sh
mkdir -p var_$VARIANT
HYPS=(free $(ls hypotheses/*.constraint.nwk | xargs -n1 basename | sed 's/.constraint.nwk//'))
H=${HYPS[$SLURM_ARRAY_TASK_ID]}
if [ "$H" = free ]; then CON=""; else CON="-g hypotheses/$H.constraint.nwk"; fi
iqtree2 -s $ALN -m $MODEL $EXTRA $CON -T $SLURM_CPUS_PER_TASK --prefix var_$VARIANT/$H -seed 1 -redo
