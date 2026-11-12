#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 8gb --time 01:00:00
#SBATCH --job-name=V9_analysis
#SBATCH --output=logs/V9_analysis.%j.log

# Intron density, family coverage per genome, exclusive sharing and clade-specific intron retention for the V9 tiers (needs run_all.sbatch).
# Submit: sbatch -p highclock run_analysis.sbatch      Absolute path, not dirname "$0".
set -euo pipefail
HERE=~/projects/intron_projector/bench/V9_antarctolithica_expanded
cd "$HERE"
mkdir -p summary
python3 intron_density.py results_projection/A/introns.tsv proteomes genomes.tsv > summary/intron_density.tsv
python3 family_coverage.py families genomes.tsv > summary/family_coverage.tsv
for t in A B; do
    python3 exclusive_sharing.py results_projection/$t/combined/matrix.tsv genomes.tsv > summary/exclusive_sharing_$t.tsv
    python3 clade_retention.py results_projection/$t/combined/matrix.tsv genomes.tsv > summary/clade_retention_$t.tsv
done
python3 tests/test_clade_retention.py; python3 tests/test_make_trees.py; python3 tests/test_scan_stats.py; python3 tests/test_make_families_pruned.py | tail -1
echo analysis done
