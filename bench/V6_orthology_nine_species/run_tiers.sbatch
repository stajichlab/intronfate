#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 32gb --time 04:00:00
#SBATCH --job-name=V6_tiers
#SBATCH --output=logs/V6_tiers.%j.log

# Tiered families from the OrthoFinder results -> MAFFT -> intron-projector analyze, per tier (A_strict, B_partial, C_hog, D_synteny).
# Needs results/$VARIANT (run_orthofinder.sbatch) and refined/$VARIANT/Refined_SOGs.w$WINDOW.tsv (run_refine.sbatch).
# Submit:  sbatch -p highclock --export=ALL,VARIANT=plain run_tiers.sbatch
# Purgeable: families/, alignments/, results_projection/ in this folder. Absolute path, not dirname "$0" (Slurm copies the script).
set -euo pipefail
HERE=~/projects/intron_projector/bench/V6_orthology_nine_species
REPO=~/projects/intron_projector
cd "$HERE"
VARIANT=${VARIANT:-plain}; WINDOW=${WINDOW:-8}; MIN=${MIN:-7}
module load rust >/dev/null 2>&1 || true
module load mafft >/dev/null 2>&1 || true
find "$REPO/crates" -name '*.rs' -exec touch {} +      # rsync keeps laptop mtimes; make cargo rebuild
( cd "$REPO" && cargo build --release --locked )
BIN=$REPO/target/release/intron-projector
python3 -m pytest -q tests >/dev/null 2>&1 || python3 tests/test_make_families.py     # the fixture test first

rm -rf families; python3 make_families.py --of "results/$VARIANT" --data data --min "$MIN" --out families \
    --refined "refined/$VARIANT/Refined_SOGs.w$WINDOW.tsv"
labels=$(grep -v '^#' genomes.tsv | cut -f1)
gffs=(); prots=(); for l in $labels; do gffs+=(--gff "data/$l/$l.gff3.gz"); prots+=(-p "data/$l/$l.proteins.fa.gz"); done
for tier in A_strict B_partial C_hog D_synteny; do
    n=$(ls families/$tier/families/*.faa 2>/dev/null | wc -l || true); [ "$n" -gt 0 ] || { echo "$tier: no families"; continue; }
    mkdir -p alignments/$tier results_projection
    ls families/$tier/families/*.faa | xargs -P "${SLURM_CPUS_PER_TASK:-8}" -I{} sh -c \
        'b=$(basename {} .faa); [ -s alignments/'$tier'/$b.aln.fasta ] || mafft --auto --thread 1 --quiet {} > alignments/'$tier'/$b.aln.fasta'
    alns=(); for f in alignments/$tier/*.aln.fasta; do alns+=(--msa "$f"); done
    rm -rf results_projection/$tier
    /usr/bin/time -f "$tier\t%e s\t%M kB" "$BIN" -v analyze "${gffs[@]}" "${prots[@]}" "${alns[@]}" \
        --taxon-map families/$tier/taxon_map.tsv -t species.nwk -m dollo,mp,counts -o results_projection/$tier \
        > results_projection/$tier.stdout 2> results_projection/$tier.stderr || { tail -20 results_projection/$tier.stderr; exit 1; }
    echo "== $tier: $n families"; tail -3 results_projection/$tier.stderr
done
