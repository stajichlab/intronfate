#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 00:40:00
#SBATCH --job-name=V6_pixi_install
#SBATCH --output=logs/V6_pixi_install.%j.log

# Does pixi build and install orthoSynAssign? Then run the package's own example (read-only copy of its example dir).
set -euo pipefail
export PATH=$HOME/.pixi/bin:$PATH
# An absolute path: under Slurm $0 is a copy in the scheduler spool, so dirname "$0" points at the wrong place.
HERE=~/projects/intron_projector/bench/V6_orthology_nine_species
cd "$HERE"
pixi --version
time pixi install
pixi run orthosynassign --version
pixi run orthosynassign --help | head -8
EX=~/projects/orthoSynAssign/example
if [ -d "$EX" ]; then
    mkdir -p "$SCRATCH/ex" && cd "$SCRATCH/ex"
    pixi run --manifest-path "$HERE/pixi.toml" \
        orthosynassign --og_file "$EX/orthogroups.tsv" --bed "$EX"/*.bed -o Refined_SOGs.tsv -t 2
    wc -l Refined_SOGs.tsv "$EX/orthogroups.tsv"
    # compare with the package's shipped result for the same input
    [ -f "$EX/Refined_SOGs.tsv" ] && { diff <(sort Refined_SOGs.tsv) <(sort "$EX/Refined_SOGs.tsv") > /dev/null && echo "matches the shipped example output" || echo "DIFFERS from the shipped example output"; }
fi
