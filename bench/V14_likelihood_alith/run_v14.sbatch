#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 32gb --time 12:00:00 -p epyc,intel
#SBATCH --job-name=V14_ml
#SBATCH --output=logs/V14.%j.log

# V14: the likelihood model with annotation error (-m likelihood) on the V9 tier-A matrix: the gene tree under several error sets
# (with a family bootstrap for the measured set), and all 70 Antarctolithica placements compared by fitted lnL and RELL.
# REPO: a copy of the branch with -m likelihood, built (default: the D1 copy). Absolute paths, not $0.
set -euo pipefail
HERE=~/projects/intron_projector/bench/V14_likelihood_alith
V9=~/projects/intron_projector/bench/V9_antarctolithica_expanded
REPO=${REPO:-$HOME/projects/intron_projector/bench/D1_howto_public_genomes/repo}
IF=$REPO/target/release/intronfate
M=$V9/results_projection/A/combined/matrix.tsv
BOOT=${BOOT:-100}
PLACE_SETS=${PLACE_SETS:-none measured alith_sn0.70}
cd "$HERE"
python3 -I tests/test_v14.py
git -C "$REPO" log --oneline -1 2>/dev/null || echo "REPO $REPO (no git metadata)"
python3 -I make_error_sets.py --v4 "$REPO/bench/V4_annotation_accuracy/results_v6/A_strict/busco_accuracy.tsv" \
    --genomes "$V9/genomes.tsv" --out errors

echo "== gene tree, every error set ($(date +%T))"
mkdir -p score
for e in errors/*.tsv; do
  set_=$(basename "$e" .tsv)
  extra=(); [ "$set_" = measured ] && extra=(--ml-bootstrap "$BOOT")
  start=$(date +%s)
  "$IF" score --matrix "$M" -t "$V9/trees/gene_tree.nwk" -m dollo,likelihood --error-rates "$e" "${extra[@]}" \
      -o "score/gene_$set_" > "score/gene_$set_.log" 2>&1
  echo "$set_: $(( $(date +%s) - start )) s; $(grep -E '^ml_(log_likelihood|iterations|converged)' score/gene_$set_/summary_score.tsv | tr '\n' ' ')"
done

echo "== placements ($(date +%T))"
for set_ in $PLACE_SETS; do
  mkdir -p "place/$set_"
  ls "$V9"/trees/place_*.nwk "$V9/trees/gene_tree.nwk" | RAYON_NUM_THREADS=4 xargs -P 4 -I{} sh -c '
    id=$(basename {} .nwk)
    [ -s "place/'"$set_"'/$id/likelihood_sites.tsv" ] || "'"$IF"'" score --matrix "'"$M"'" -t {} -m likelihood --error-rates "errors/'"$set_"'.tsv" \
        -o "place/'"$set_"'/$id" > "place/'"$set_"'/$id.log" 2>&1'
  ids=$(ls -d place/$set_/*/ | xargs -n1 basename | paste -sd, -)
  mkdir -p results
  python3 -I rell.py --score "place/$set_" --ids "$ids" --reps 1000 > "results/rell_$set_.tsv"
  grep -h ml_converged place/$set_/*/summary_score.tsv | sort | uniq -c
  echo "$set_: $(sed -n 3,6p results/rell_$set_.tsv | cut -f1-4 | tr '\n' ' ')"
done
echo "== done ($(date +%T))"
