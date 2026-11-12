#!/usr/bin/env bash
# Run `intronfate project` on the genomes listed in a TSV (label, gff, proteins, bytes) and record
# wall time, peak RSS, the tool's own [time]/[size] stage lines, and the sizes of what it wrote.
# Usage: run_one.sh <list.tsv> <tag> [threads]      (writes results/<tag>/)
set -euo pipefail
cd "$(dirname "$0")"
list=$1; tag=$2; threads=${3:-}
REPO=$(cd ../.. && pwd)
BIN=${INTRONFATE:-$REPO/target/release/intronfate}
out=results/$tag; rm -rf "$out"; mkdir -p "$out"
args=()
while IFS=$'\t' read -r label gff prot bytes; do args+=(--gff "$gff" -p "$prot"); done < "$list"
n=$(wc -l < "$list")
envs=(); [ -n "$threads" ] && envs=(env RAYON_NUM_THREADS="$threads")
"${envs[@]}" /usr/bin/time -f 'wall_s\t%e\nuser_s\t%U\nsys_s\t%S\nmax_rss_kb\t%M' -o "$out/time.tsv" \
    "$BIN" -v project "${args[@]}" -o "$out/o" > "$out/stdout.txt" 2> "$out/stderr.txt" || { tail -5 "$out/stderr.txt" >&2; exit 1; }
{
  echo -e "tag\t$tag"; echo -e "genomes\t$n"; echo -e "threads\t${threads:-default}"
  echo -e "gff_bytes_gz\t$(awk -F'\t' '{s+=$4} END{print s}' "$list")"
  echo -e "output_bytes\t$(du -sb "$out/o" | cut -f1)"
  echo -e "introns_tsv_bytes\t$(stat -c %s "$out/o/introns.tsv")"
} > "$out/meta.tsv"
cp "$out/o/summary_project.tsv" "$out/summary_project.tsv" 2>/dev/null || true
rm -rf "$out/o"      # the tables are large; keep only timings and summaries
echo "$tag: $(grep wall_s "$out/time.tsv" | cut -f2) s, $(awk -F'\t' '/max_rss_kb/{printf "%.0f MB", $2/1024}' "$out/time.tsv")"
