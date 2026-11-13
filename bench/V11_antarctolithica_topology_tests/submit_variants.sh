#!/usr/bin/bash
# Run from /rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests on the HPCC (after submit_all.sh finished):  bash submit_variants.sh [variant ...]   (default: pmsf fast25 fast50 dayhoff6)
set -euo pipefail
cd /rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests
V=("$@"); [ ${#V[@]} -eq 0 ] && V=(pmsf fast25 fast50 dayhoff6)
P=$(sbatch --parsable -p highclock run_variants_prep.sbatch)
for v in "${V[@]}"; do
  T=12:00:00; [ "$v" = pmsf ] && T=48:00:00
  A=$(sbatch --parsable -p highclock --dependency=afterok:$P --time=$T --job-name=V11_var_$v --export=ALL,VARIANT=$v run_variant.sbatch)
  U=$(sbatch --parsable -p highclock --dependency=afterok:$A --job-name=V11_au_$v --export=ALL,VARIANT=$v run_variant_au.sbatch)
  echo "$v: prep $P, searches $A, tests $U"
done
