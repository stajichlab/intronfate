#!/bin/bash
#SBATCH -p stajichlab
#SBATCH -c 32
#SBATCH --mem=96G
#SBATCH -t 10:00:00
#SBATCH -J P6_topo
#SBATCH -o P6_topo.%j.log
# P6: likelihoods of fixed topologies and of constrained searches (FastTree only), on the filtered concatenation. Writes only here ($PWD).
set -euo pipefail
module load singularity 2>/dev/null || module load apptainer
IMG=/bigdata/stajichlab/shared/singularity_cache/ghcr.io-stajichlab-phyling-v2.4.1.2.img
BASE=/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/phylo/results/protein
ALN=$BASE/buildtree/protein-atlas_genus_reps-taxa_1654.fungi_odb10.fa
FT=$BASE/buildtree/fungi_odb10/fasttree/protein-atlas_genus_reps-taxa_1654.fungi_odb10.fasttree.nosupport.treefile
AS=$BASE/tree/fungi_odb10/astral/consensus.nw
run() { singularity exec --no-home --pid -B /bigdata/stajichlab -B "$PWD" "$IMG" env PATH=/opt/conda/bin:/usr/local/bin:/usr/bin:/bin OMP_NUM_THREADS=32 "$@"; }
python3 topology_tests.py --selftest
python3 topology_tests.py constraints cons
ls $BASE/tree/fungi_odb10/ ; GT=$(ls $BASE/tree/fungi_odb10/*.nw 2>/dev/null | head -1 || true)
# gene trees: concatenate the per-marker trees saved by phyling (ASTRAL input)
cp $BASE/tree/fungi_odb10/astral/trees.nw gene_trees.nw
python3 topology_tests.py concordance gene_trees.nw . | tee concordance.txt
# strip support labels and branch lengths are kept; FastTree -intree reads Newick
evalt() { # name tree
  echo "== eval $1"; run FastTreeMP -lg -nome -mllen -intree "$2" "$ALN" > "eval_$1.nw" 2> "eval_$1.log"; grep -i "loglk\|Total time" "eval_$1.log" | tail -3; }
evalt fasttree "$FT"
evalt astral "$AS"
for g in Microsporidia Kickxellomycota; do
  echo "== constrained search $g"
  run FastTreeMP -lg -nosupport -constraints cons/constraint_$g.fa -constraintWeight 100 "$ALN" > "cons_$g.nw" 2> "cons_$g.log"
  tail -2 "cons_$g.log"
  evalt "cons_$g" "cons_$g.nw"
done
echo ALL_DONE
