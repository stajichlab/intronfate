#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 16gb --time 04:00:00
#SBATCH --job-name=D1_howto
#SBATCH --output=logs/D1_howto.%j.log

# Runs docs/HOWTO.md as written: every ```sh block from step 1 on, extracted by scripts/doc_blocks.py, in $SCRATCH.
#   SET=fungi5 sbatch run_howto.sbatch      (Example 1, the default)
#   SET=fungi9 sbatch run_howto.sbatch      (Example 2)
# Only two lines of the page are replaced: the working folder (-> $SCRATCH) and, for fungi9, the example line the page tells the reader to
# change. `apptainer` is mapped to the site's `singularity` (the page says either works). The small results are copied to results/<SET>/.
# REPO is a copy of the repository under test in ./repo (rsync it there, then `cargo build --release` in it); absolute paths, not $0.
set -euo pipefail
: "${SCRATCH:?run inside a Slurm job}"
SET=${SET:-fungi5}
HERE=~/projects/intron_projector/bench/D1_howto_public_genomes
REPO=$HERE/repo
export PATH=$REPO/target/release:$PATH
module load singularity
apptainer() { singularity "$@"; }; export -f apptainer

python3 -I "$REPO/scripts/test_howto_scripts.py"           # the helper scripts' fixture tests first

WORK=$SCRATCH/howto_$SET; mkdir -p "$WORK"
[ -s "$HERE/containers/orthofinder_3.2.sif" ] && cp "$HERE/containers/orthofinder_3.2.sif" "$WORK/"   # skips the 650 MB download (the page's `[ -s ... ] ||`)
subs=(--sub 'mkdir -p ~/intronfate-howto && cd ~/intronfate-howto' "cd $WORK")
[ "$SET" = fungi5 ] || subs+=(--sub 'SET=fungi5; MIN=4; OUTGROUP=Cneoformans' 'SET=fungi9; MIN=7; OUTGROUP=Spunctatus')
python3 -I "$REPO/scripts/doc_blocks.py" "$REPO/docs/HOWTO.md" --from "## 1." --trace "${subs[@]}" > "$SCRATCH/howto.sh"

start=$(date +%s)
REPO=$REPO bash "$SCRATCH/howto.sh" 2>&1 | tee "$SCRATCH/howto.out"
echo "total wall time: $(( $(date +%s) - start )) s"

out=$HERE/results/$SET; rm -rf "$out"; mkdir -p "$out"
cp "$SCRATCH/howto.out" "$out/howto.out"
cd "$WORK"
cp families/summary.tsv "$out/families_summary.tsv"
cp results/A_strict/diagnostics.tsv results/A_strict/summary_project.tsv "$out/" 2>/dev/null || cp results/A_strict/diagnostics.tsv "$out/"
cp results/A_strict/combined/dollo_branches.tsv results/A_strict/combined/summary_score.tsv "$out/"
[ -d results/B_partial ] && cp results/B_partial/combined/dollo_branches.tsv "$out/dollo_branches.B_partial.tsv"
wc -l results/A_strict/unmatched.tsv results/A_strict/combined/dollo_sites.tsv > "$out/counts.txt"
du -sh --exclude=orthofinder_3.2.sif . genomes proteomes orthofinder families alignments results >> "$out/counts.txt"
echo "results in $out"
