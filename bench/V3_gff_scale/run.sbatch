#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 16 --mem 64gb --time 04:00:00
#SBATCH --job-name=V3_gff_scale
#SBATCH --output=logs/V3_gff_scale.%j.log

# V3: how does the projection path perform on genome-scale annotation? (reads BFD GFF3/proteins read-only)
# Submit from ~/projects/intron_projector/bench/V3_gff_scale/ :
#     mkdir -p logs && sbatch -p highclock run.sbatch
# Purgeable: data/ (lists only) and results/ (timings only; the output tables are deleted after each run).
set -euo pipefail
module load rust >/dev/null 2>&1 || true
REPO=$(cd ../.. && pwd)
find "$REPO/crates" -name '*.rs' -exec touch {} +      # rsync keeps laptop mtimes: force the rebuild
( cd "$REPO" && cargo build --release --locked )
export INTRON_PROJECTOR=$REPO/target/release/intron-projector
. ./config.sh
mkdir -p data results
if [ ! -s data/gff_sizes.tsv ]; then       # one slow pass over the shared filesystem, then cached
    # find, not a glob: 23 000 paths would exceed the argument-length limit
    find "$BFD" -mindepth 3 -maxdepth 3 -path '*/predict_results/*.gff3.gz' -printf '%s\t%P\n' \
        | awk -F'\t' '{split($2,p,"/"); print $1 "\t" p[1]}' > data/gff_sizes.tsv
fi
python3 select.py "$BFD" "$TOP_K" "$G_STEPS" "$SEED"

# A. the largest GFF3s, one genome per process
i=0; while IFS=$'\t' read -r label gff prot bytes; do
    printf '%s\t%s\t%s\t%s\n' "$label" "$gff" "$prot" "$bytes" > "data/one_$i.tsv"
    ./run_one.sh "data/one_$i.tsv" "top_$label"; i=$((i+1))
done < data/largest.tsv

# B. thread scaling on the very largest
head -1 data/largest.tsv > data/largest1.tsv
for t in $THREADS; do ./run_one.sh data/largest1.tsv "threads_$t" "$t"; done

# C. many genomes in one process
for n in $G_STEPS; do ./run_one.sh "data/random_$n.tsv" "genomes_$n"; done

# C2. thread scaling with many genomes (files are now read in parallel)
for t in ${G_THREADS:-1 4 16}; do ./run_one.sh data/random_100.tsv "gthreads_$t" "$t"; done

# D. what gzip decoding costs: same genome from a plain copy in $SCRATCH
: "${SCRATCH:?run inside a Slurm job}"
IFS=$'\t' read -r label gff prot bytes < data/largest1.tsv
zcat "$gff" > "$SCRATCH/$label.gff3"
printf '%s\t%s\t%s\t%s\n' "$label" "$SCRATCH/$label.gff3" "$prot" "$bytes" > data/plain1.tsv
./run_one.sh data/largest1.tsv gzip_largest
./run_one.sh data/plain1.tsv plain_largest
echo "uncompressed GFF3 of the largest: $(stat -c %s "$SCRATCH/$label.gff3") bytes" | tee results/uncompressed_largest.txt

python3 summarize_scale.py | tee results/REPORT.txt
