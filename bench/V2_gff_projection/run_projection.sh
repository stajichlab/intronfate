#!/usr/bin/env bash
# The step under test: GFF3 -> protein space -> alignment columns -> matrices -> Dollo/MP, on real annotation.
set -euo pipefail
cd "$(dirname "$0")"; . ./config.sh
REPO=$(cd ../.. && pwd)
BIN=${INTRON_PROJECTOR:-$REPO/target/release/intron-projector}
OUT=results/projection; rm -rf "$OUT"; mkdir -p "$OUT"
labels=$(grep -v '^#' genomes.tsv | cut -f1)
gffs=(); prots=()
for l in $labels; do gffs+=(--gff "data/$l/$l.gff3.gz"); prots+=(-p "data/$l/$l.proteins.fa.gz"); done

# (1) GFF parsing + protein-space projection, per genome, timed with peak memory.
mkdir -p "$OUT/per_genome"
printf 'label\twall_s\tmax_rss_kb\n' > "$OUT/per_genome_timing.tsv"
for l in $labels; do
    /usr/bin/time -f "$l\t%e\t%M" -a -o "$OUT/per_genome_timing.tsv" \
        "$BIN" project --gff "data/$l/$l.gff3.gz" -p "data/$l/$l.proteins.fa.gz" -o "$OUT/per_genome/$l" > /dev/null 2> "$OUT/per_genome/$l.log"
done

# (2) The whole family set in one process: every genome, every alignment, species tree via --taxon-map.
alns=(); for f in alignments/*.aln.fasta; do alns+=(--msa "$f"); done
/usr/bin/time -f "analyze_all\t%e\t%M" -a -o "$OUT/per_genome_timing.tsv" \
    "$BIN" -v analyze "${gffs[@]}" "${prots[@]}" "${alns[@]}" --taxon-map taxon_map.tsv \
        -t species.nwk -m dollo,mp,counts -o "$OUT/all" > "$OUT/analyze.stdout" 2> "$OUT/analyze.stderr" || { tail -20 "$OUT/analyze.stderr" >&2; exit 1; }
python3 summarize_projection.py "$OUT" | tee "$OUT/REPORT.txt"
