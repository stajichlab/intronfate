# Shared settings for the V2 scripts (sourced). Override from the environment.
BFD=${BFD:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation}
BUSCO_LINEAGE=${BUSCO_LINEAGE:-/srv/projects/db/BUSCO/v10/lineages/fungi_odb10}
# The `busco` module only loads miniconda and sets BUSCO_LINEAGES; the program is a standalone prefix.
BUSCO_HOME=${BUSCO_HOME:-/opt/linux/rocky/8.x/x86_64/pkgs/busco/5.8.0}
PARTITION=${PARTITION:-highclock}
# A BUSCO family is used when it is single-copy ("Complete") in at least this many of the genomes.
MIN_GENOMES=${MIN_GENOMES:-7}
