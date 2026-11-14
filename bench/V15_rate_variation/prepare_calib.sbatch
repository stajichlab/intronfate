#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 2:00:00 -p epyc,intel
#SBATCH --job-name=V15_prep
#SBATCH --output=logs/V15prep.%j.log

# Regenerate V13's calibration matrices and trees (V13's run_v13.sbatch step C, same arguments; the thinning is seeded) into
# calib_v13/, check that the plan matches V13's, and write the job lists for calib_ml.sbatch:
#   calib_jobs.txt     every "all"-mode matrix (165)
#   calib_jobs_k4.txt  the full matrices and age-matched replicates 1-3 (60), for the fits with rate classes
set -euo pipefail
HERE=~/projects/intron_projector/bench/V15_rate_variation
V9=~/projects/intron_projector/bench/V9_antarctolithica_expanded
V13=~/projects/intron_projector/bench/V13_antarctolithica_intron_map
M=$V9/results_projection/A/combined/matrix.tsv
CONTROLS=Pyronema,Orbilia,Ncrassa,Botrytis,Geoglossum,Anidulans,Lipomyces,Neoaidaea,Babjevia,Yarrowia,Geotrichum,Saitoella,Neolecta,Taphrina,Spombe
cd "$HERE"
rm -rf calib_v13
python3 -I "$V13/calibrate.py" prepare --matrix "$M" --gene-tree "$V9/trees/gene_tree.nwk" --genomes "$V9/genomes.tsv" \
    --gene-sites "$V13/score/gene_tree/dollo_sites.tsv" --gene-branches "$V13/score/gene_tree/dollo_branches.tsv" \
    --controls "$CONTROLS" --reps 5 --make-trees "$V9/make_trees.py" --out calib_v13 > /dev/null
diff calib_v13/plan.tsv "$V13/results/calibration_plan.tsv" && echo "plan identical to V13's"
awk -F'\t' '$4=="all"' calib_v13/jobs.txt > calib_jobs.txt
awk -F'\t' '$2=="full" || ($2=="age" && $3<=3)' calib_jobs.txt > calib_jobs_k4.txt
wc -l calib_jobs.txt calib_jobs_k4.txt
