#!/usr/bin/bash
# Restricted-matrix runs for the site-heterogeneous (PMSF C60) test: keep the slowest 50,000 of the 200,000 sites (variants/rates.rate from run_variants_prep.sbatch),
# then the nine hypothesis searches and the tree tests under LG+C60+F+G4 (PMSF) and, for comparison, under LG+F+G4 (fast75).
# Run from /rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests on the HPCC:  bash submit_pmsf50k.sh
set -euo pipefail
cd /rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests
F=$(sbatch --parsable -p highclock --mem 4gb --time 00:20:00 --job-name=V11_fast75 --output=/rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests/logs/V11_fast75.%j.log \
     --wrap "cd /rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests && python3 site_filter.py drop-fast supermatrix.fasta variants/rates.rate 0.75 variants/fast75.fasta")
A=$(sbatch --parsable -p highclock --dependency=afterok:$F --time=48:00:00 --mem=110gb --job-name=V11_var_pmsf50k --export=ALL,VARIANT=pmsf50k run_variant.sbatch)
U=$(sbatch --parsable -p highclock --dependency=afterok:$A --mem=110gb --job-name=V11_au_pmsf50k --export=ALL,VARIANT=pmsf50k run_variant_au.sbatch)
B=$(sbatch --parsable -p highclock --dependency=afterok:$F --time=12:00:00 --job-name=V11_var_fast75 --export=ALL,VARIANT=fast75 run_variant.sbatch)
V=$(sbatch --parsable -p highclock --dependency=afterok:$B --job-name=V11_au_fast75 --export=ALL,VARIANT=fast75 run_variant_au.sbatch)
echo "fast75 filter $F; pmsf50k searches $A tests $U; fast75 searches $B tests $V"
