#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 32gb --time 08:00:00 -p epyc,intel
#SBATCH --job-name=V13_alith
#SBATCH --output=logs/V13.%j.log

# V13: Antarctolithica intron positions on the V9 tier-A matrix: rescoring on the gene tree and all 39 placements, tables (analyze.py),
# calibration with intron-poor stand-ins (calibrate.py, run_calibration.sh). Figures: run_figures.sh (pixi). Absolute paths, not $0.
# INTRONFATE: a release build of the branch that has `-m dollo_weighted` (default: the D1 copy of it).
set -euo pipefail
HERE=~/projects/intron_projector/bench/V13_antarctolithica_intron_map
V9=~/projects/intron_projector/bench/V9_antarctolithica_expanded
export INTRONFATE=${INTRONFATE:-$HOME/projects/intron_projector/bench/D1_howto_public_genomes/repo/target/release/intronfate}
M=$V9/results_projection/A/combined/matrix.tsv
CONTROLS=${CONTROLS:-Pyronema,Orbilia,Ncrassa,Botrytis,Geoglossum,Anidulans,Lipomyces,Neoaidaea,Babjevia,Yarrowia,Geotrichum,Saitoella,Neolecta,Taphrina,Spombe}
REPS=${REPS:-5}
cd "$HERE"
python3 -I tests/test_v13.py
"$INTRONFATE" --version 2>/dev/null || true

echo "== A. score the gene tree and every placement ($(date +%T))"
mkdir -p score
"$INTRONFATE" score --matrix "$M" -t "$V9/trees/gene_tree.nwk" -m dollo,mp,dollo_weighted --gain-loss-weights 1:1,2:1,3:1,5:1,10:1,20:1 \
    --reconstructions -o score/gene_tree > score/gene_tree.log 2>&1
ls "$V9"/trees/place_*.nwk | xargs -P "$SLURM_CPUS_PER_TASK" -I{} sh -c \
    'id=$(basename {} .nwk); "$INTRONFATE" score --matrix '"$M"' -t {} -m dollo,mp -o score/$id > score/$id.log 2>&1'
ls -d score/place_* | wc -l

echo "== B. tables ($(date +%T))"
python3 -I analyze.py --matrix "$M" --trees "$V9/trees" --genomes "$V9/genomes.tsv" --score score --out results

echo "== C. calibration ($(date +%T))"
rm -rf calib; mkdir -p calib
python3 -I calibrate.py prepare --matrix "$M" --gene-tree "$V9/trees/gene_tree.nwk" --genomes "$V9/genomes.tsv" \
    --gene-sites score/gene_tree/dollo_sites.tsv --gene-branches score/gene_tree/dollo_branches.tsv \
    --controls "$CONTROLS" --reps "$REPS" --make-trees "$V9/make_trees.py" --out calib
printf 'control\tvariant\trep\tmode\ttree\tdollo_losses\n' > calib/scores.tsv
xargs -P "$SLURM_CPUS_PER_TASK" -L 1 bash "$HERE/run_calibration.sh" < calib/jobs.txt >> calib/scores.tsv
wc -l calib/scores.tsv
python3 -I calibrate.py summarize --calib calib --scores calib/scores.tsv --out results
cp calib/plan.tsv results/calibration_plan.tsv
find calib -name '*.matrix.tsv' -type f -delete                     # the thinned matrices are large and reproducible
echo "== done ($(date +%T))"
