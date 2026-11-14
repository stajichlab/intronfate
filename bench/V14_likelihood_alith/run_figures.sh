#!/usr/bin/env bash
# V14 tables and figures after run_v14.sbatch: RELL over the 70 placements only (the gene tree is the same topology as one of them,
# so it is left out to avoid splitting that placement's wins), then v14_figures.py with V13's pixi environment.
#   srun -p short -c 4 --mem 8gb --time 0:30:00 bash run_figures.sh
set -euo pipefail
HERE=~/projects/intron_projector/bench/V14_likelihood_alith
V9=~/projects/intron_projector/bench/V9_antarctolithica_expanded
V13=~/projects/intron_projector/bench/V13_antarctolithica_intron_map
cd "$HERE"
for d in place/*/; do
  set_=$(basename "$d")
  ids=$(ls -d place/$set_/place_*/ | xargs -n1 basename | paste -sd, -)
  python3 -I rell.py --score "place/$set_" --ids "$ids" --reps 1000 > "results/rell_$set_.tsv"
done
pixi run --manifest-path "$V13/pixi.toml" python -I v14_figures.py --trees "$V9/trees" --genomes "$V9/genomes.tsv" \
    --score score --results results --out figures --v13results "$V13/results"
