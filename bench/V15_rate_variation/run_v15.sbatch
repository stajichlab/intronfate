#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 32gb --time 48:00:00 -p epyc,intel
#SBATCH --job-name=V15_ml
#SBATCH --output=logs/V15.%j.log

# V15: the likelihood model with rate variation across positions (--ml-rate-categories) on the V9 tier-A matrix, as V14 without it:
#   A. the gene tree under the V14 error sets none / measured / alith_sn0.70, with 1 and CATS rate categories (bootstrap BOOT
#      replicates for measured with CATS categories; per-cycle trace of the CATS fits in score/*.log);
#   B. all 70 Antarctolithica placements with CATS categories under PLACE_SETS, compared by fitted lnL and RELL (V14's rell.py).
# REPO: a built copy of the branch (default: the D1 copy). Absolute paths, not $0.
set -euo pipefail
HERE=~/projects/intron_projector/bench/V15_rate_variation
V9=~/projects/intron_projector/bench/V9_antarctolithica_expanded
V14=~/projects/intron_projector/bench/V14_likelihood_alith
REPO=${REPO:-$HOME/projects/intron_projector/bench/D1_howto_public_genomes/repo}
IF=$REPO/target/release/intronfate
M=$V9/results_projection/A/combined/matrix.tsv
CATS=${CATS:-4}
BOOT=${BOOT:-100}
PLACE_SETS=${PLACE_SETS:-none measured}
GENE_SETS=${GENE_SETS:-none measured alith_sn0.70}
cd "$HERE"
python3 -I tests/test_v15.py
"$IF" --version || true

echo "== A. gene tree ($(date +%T))"
mkdir -p score results
for set_ in $GENE_SETS; do
  for k in 1 $CATS; do
    out=score/gene_${set_}_k$k
    [ -s "$out/summary_score.tsv" ] && { echo "skip $out"; continue; }
    extra=(); [ "$set_" = measured ] && [ "$k" != 1 ] && extra=(--ml-bootstrap "$BOOT")
    start=$(date +%s)
    "$IF" -vv score --matrix "$M" -t "$V9/trees/gene_tree.nwk" -m dollo,likelihood --error-rates "$V14/errors/$set_.tsv" \
        --ml-rate-categories "$k" "${extra[@]}" -o "$out" > "$out.log" 2>&1
    echo "$set_ k$k: $(( $(date +%s) - start )) s; $(grep -E '^ml_(log_likelihood|iterations|converged|gamma_shape|expected_losses|root_expected_present)' "$out/summary_score.tsv" | tr '\n' ' ')"
  done
done

echo "== B. placements, $CATS categories ($(date +%T))"
for set_ in $PLACE_SETS; do
  d=place/${set_}_k$CATS
  mkdir -p "$d"
  ls "$V9"/trees/place_*.nwk "$V9/trees/gene_tree.nwk" | RAYON_NUM_THREADS=4 xargs -P 4 -I{} sh -c '
    id=$(basename {} .nwk)
    [ -s "'"$d"'/$id/likelihood_sites.tsv" ] || "'"$IF"'" score --matrix "'"$M"'" -t {} -m likelihood --error-rates "'"$V14"'/errors/'"$set_"'.tsv" \
        --ml-rate-categories '"$CATS"' -o "'"$d"'/$id" > "'"$d"'/$id.log" 2>&1'
  ids=$(ls -d "$d"/*/ | xargs -n1 basename | paste -sd, -)
  python3 -I "$V14/rell.py" --score "$d" --ids "$ids" --reps 1000 > "results/rell_${set_}_k$CATS.tsv"
  grep -h ml_converged "$d"/*/summary_score.tsv | sort | uniq -c
  echo "$set_: $(sed -n 3,6p "results/rell_${set_}_k$CATS.tsv" | cut -f1-4 | tr '\n' ' ')"
done
echo "== done ($(date +%T))"
