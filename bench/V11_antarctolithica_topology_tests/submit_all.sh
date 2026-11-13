#!/usr/bin/bash
# Run from /rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests on the HPCC:  bash submit_all.sh
set -euo pipefail
cd /rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests
J1=$(sbatch --parsable -p highclock run_supermatrix.sbatch)
J2=$(sbatch --parsable -p highclock --dependency=afterok:$J1 run_hyp.sbatch)
J3=$(sbatch --parsable -p highclock --dependency=afterok:$J2 run_au.sbatch)
echo "supermatrix $J1, hypotheses $J2, tests $J3"
