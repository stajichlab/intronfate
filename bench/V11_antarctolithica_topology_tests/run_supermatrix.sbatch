#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 1 --mem 8gb --time 01:00:00
#SBATCH --job-name=V11_supermatrix
#SBATCH --output=/rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests/logs/V11_supermatrix.%j.log
# Supermatrix from the V9 tier-A alignments + the constraint trees for the candidate Alith placements. Absolute paths (no dirname $0).
set -euo pipefail
cd /rhome/jstajich/projects/intron_projector/bench/V11_antarctolithica_topology_tests
V9=../V9_antarctolithica_expanded
python3 build_supermatrix.py --aln-dir $V9/alignments/A --taxon-map $V9/families/A/taxon_map.tsv --genomes $V9/genomes.tsv --out supermatrix \
    --min-lineages 35 --max-gap 0.5 --min-cols 100 --max-sites 200000
python3 make_hypotheses.py $V9/genomes.tsv hypotheses
ls hypotheses
