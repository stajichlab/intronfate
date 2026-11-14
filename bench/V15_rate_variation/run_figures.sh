#!/usr/bin/env bash
# V15 tables and figures after run_v15.sbatch and the calibration arrays (calib_ml.sbatch): RELL over the 70 placements only (as V14),
# the calibration summaries per number of classes, then v15_figures.py with V13's pixi environment.
#   srun -p short -c 4 --mem 8gb --time 0:30:00 bash run_figures.sh
set -euo pipefail
HERE=~/projects/intron_projector/bench/V15_rate_variation
V9=~/projects/intron_projector/bench/V9_antarctolithica_expanded
V13=~/projects/intron_projector/bench/V13_antarctolithica_intron_map
V14=~/projects/intron_projector/bench/V14_likelihood_alith
K=${K:-4}
cd "$HERE"
mkdir -p results figures
for d in place/*/; do
  set_=$(basename "$d")
  ids=$(ls -d "place/$set_"/place_*/ | xargs -n1 basename | paste -sd, -)
  python3 -I "$V14/rell.py" --score "place/$set_" --ids "$ids" --reps 1000 > "results/rell_$set_.tsv"
done
for k in 1 "$K"; do
  ls calib/*_k$k.tsv > /dev/null 2>&1 && python3 -I summarize_ml_calibration.py --scores "calib/*_k$k.tsv" --calib calib_v13 \
      --out "results/calib_k$k" > "results/calib_k$k.txt"
done
pixi run --manifest-path "$V13/pixi.toml" python -I v15_figures.py --trees "$V9/trees" --genomes "$V9/genomes.tsv" \
    --score score --results results --out figures --v13results "$V13/results" --k "$K"
