#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 4 --mem 8gb --time 00:20:00
#SBATCH --job-name=V6_combined_only
#SBATCH --output=logs/V6_combined_only.%j.log

# `analyze --combine --combined-only` on the V6 tier A: same combined/ files as the full run, far fewer output files.
# Needs the V6 tier job (families, alignments, results_projection/A_strict/combined). Output goes to $SCRATCH; only a diff and sizes are kept.
# Submit: sbatch -p highclock check_combined_only.sbatch      Absolute path, not dirname "$0".
set -euo pipefail
: "${SCRATCH:?run inside a Slurm job}"
HERE=~/projects/intron_projector/bench/V6_orthology_nine_species
cd "$HERE"
REPO=~/projects/intron_projector
module load rust >/dev/null 2>&1 || true
find "$REPO/crates" -name '*.rs' -exec touch {} +      # rsync keeps laptop mtimes; make cargo rebuild
( cd "$REPO" && cargo build --release --locked )
BIN=$REPO/target/release/intronfate
tier=A_strict
gffs=(); prots=(); for l in $(grep -v '^#' genomes.tsv | cut -f1); do gffs+=(--gff "data/$l/$l.gff3.gz"); prots+=(-p "data/$l/$l.proteins.fa.gz"); done
rm -rf "$SCRATCH/lean"
/usr/bin/time -f "combined-only\t%e s\t%M kB" "$BIN" analyze "${gffs[@]}" "${prots[@]}" --msa-dir alignments/$tier --msa-ext aln.fasta --combine --combined-only \
    --taxon-map families/$tier/taxon_map.tsv -t species.nwk -m dollo,mp,counts -o "$SCRATCH/lean" > /dev/null
echo "files in the full run: $(find results_projection/$tier -type f | wc -l)   in the combined-only run: $(find "$SCRATCH/lean" -type f | wc -l)"
echo "size (MB) full: $(du -sm results_projection/$tier | cut -f1)   combined-only: $(du -sm "$SCRATCH/lean" | cut -f1)"
for f in matrix.tsv dollo_branches.tsv dollo_sites.tsv mp_sites.tsv counts_sites.tsv pairwise.tsv; do
    cmp "results_projection/$tier/combined/$f" "$SCRATCH/lean/combined/$f" && echo "identical: combined/$f"
done
