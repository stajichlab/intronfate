#!/usr/bin/env bash
# BUSCO (protein mode) for one genome. Usage: busco_one.sh <label>. Temporary files go to $SCRATCH.
set -euo pipefail
cd "$(dirname "$0")"; . ./config.sh
label=$1; : "${SCRATCH:?run inside a Slurm job}"
export PATH=$BUSCO_HOME/bin:$PATH
command -v busco >/dev/null || { echo "busco not found under $BUSCO_HOME" >&2; exit 1; }
out=results/busco/$label; mkdir -p "$out"
[ -s "$out/full_table.tsv" ] && { echo "$label: already done"; exit 0; }
zcat "data/$label/$label.proteins.fa.gz" > "$SCRATCH/$label.faa"
cd "$SCRATCH"
busco -i "$label.faa" -m proteins -l "$BUSCO_LINEAGE" --offline -c "${SLURM_CPUS_PER_TASK:-4}" -o "busco_$label" -f > "busco_$label.log" 2>&1
cd - > /dev/null
cp "$SCRATCH"/busco_"$label"/run_*/full_table.tsv "$out/full_table.tsv"
cp "$SCRATCH"/busco_"$label"/short_summary*.txt "$out/" 2>/dev/null || true
tail -3 "$SCRATCH/busco_$label.log"
