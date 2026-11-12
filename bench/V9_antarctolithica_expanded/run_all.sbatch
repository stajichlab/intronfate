#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 32gb --time 06:00:00
#SBATCH --job-name=V9_families_scan
#SBATCH --output=logs/V9_families_scan.%j.log

# Families (A: >=22 lineages, B: 18-21; multi-copy genomes treated as missing; the Antarctolithica strains count as one lineage) -> MAFFT ->
# intronfate analyze --combine -> every candidate placement of the Antarctolithica clade scored.
# Needs results/plain (run_orthofinder.sbatch). Submit: sbatch -p highclock --export=ALL,VARIANT=plain run_all.sbatch
# Purgeable: families/ alignments/ results_projection/ trees/. Absolute path, not dirname "$0".
set -euo pipefail
HERE=~/projects/intron_projector/bench/V9_antarctolithica_expanded
REPO=~/projects/intron_projector
cd "$HERE"
VARIANT=${VARIANT:-plain}
module load rust >/dev/null 2>&1 || true
module load mafft >/dev/null 2>&1 || true
find "$REPO/crates" -name '*.rs' -exec touch {} +
( cd "$REPO" && cargo build --release --locked )
BIN=$REPO/target/release/intronfate
python3 tests/test_make_families_pruned.py

rm -rf families; python3 make_families_pruned.py --of "results/$VARIANT" --data data --out families --group-prefix Alith --tiers A=${TIER_A:-30},B=${TIER_B:-24}
python3 make_trees.py "results/$VARIANT/Species_Tree/SpeciesTree_rooted.txt" --group-prefix Alith --genomes genomes.tsv --root-on Mmucedo --out trees
ls data/*/*.gff3.gz > gffs.txt
mkdir -p results_projection
for tier in A B; do
    mkdir -p alignments/$tier
    ls families/$tier/families/*.faa | xargs -P "${SLURM_CPUS_PER_TASK:-8}" -I{} sh -c \
        'b=$(basename {} .faa); [ -s alignments/'$tier'/$b.aln.fasta ] || mafft --auto --thread 1 --quiet {} > alignments/'$tier'/$b.aln.fasta'
    rm -rf results_projection/$tier
    /usr/bin/time -f "$tier\t%e s\t%M kB" "$BIN" -v analyze --gff-list gffs.txt --proteins-dir proteomes --proteins-ext faa \
        --msa-dir alignments/$tier --msa-ext aln.fasta --taxon-map families/$tier/taxon_map.tsv \
        -t trees/gene_tree.nwk -m dollo,mp,counts --combine --combined-only -o results_projection/$tier \
        > results_projection/$tier.stdout 2> results_projection/$tier.stderr || { tail -20 results_projection/$tier.stderr; exit 1; }
    tail -3 results_projection/$tier.stderr
    python3 run_scan.py --matrix results_projection/$tier/combined/matrix.tsv --trees trees --bin "$BIN" --work "$SCRATCH/scan_$tier" \
        --out results_projection/$tier/placements.tsv --threads "${SLURM_CPUS_PER_TASK:-8}" | head -30
done
