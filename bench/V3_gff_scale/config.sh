# Shared settings for the V3 scripts (sourced). Override from the environment.
BFD=${BFD:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation}
TOP_K=${TOP_K:-8}                     # largest GFFs, run one at a time
G_STEPS=${G_STEPS:-"10 100 1000"}     # genome counts for the many-genomes run (one process each)
THREADS=${THREADS:-"1 2 4 8 16"}      # thread counts for the scaling run on the largest GFF
SEED=${SEED:-20261111}
PARTITION=${PARTITION:-highclock}
