#!/usr/bin/env bash
# Align every <family>.faa in a folder with MAFFT (--auto), untrimmed, N families at a time.
#   align_families.sh IN_DIR OUT_DIR [JOBS]      -> OUT_DIR/<family>.aln.fasta
# Do not trim the result (trimAl, Gblocks): intronfate needs each aligned row to equal its protein.
# An alignment that already exists is kept, so an interrupted run can be restarted. Each output is written to a temporary file and
# renamed, so a killed run never leaves a half-written alignment behind. Needs mafft and xargs on PATH (both are in the OrthoFinder
# container).
set -euo pipefail
in=${1:?usage: align_families.sh IN_DIR OUT_DIR [JOBS]}; out=${2:?usage: align_families.sh IN_DIR OUT_DIR [JOBS]}; jobs=${3:-4}
command -v mafft >/dev/null || { echo "mafft is not on PATH" >&2; exit 1; }
mkdir -p "$out"
n=$(find "$in" -maxdepth 1 -name '*.faa' | wc -l)
[ "$n" -gt 0 ] || { echo "no .faa files in $in" >&2; exit 1; }
# xargs appends the file after the fixed arguments, so the output folder travels in the environment.
find "$in" -maxdepth 1 -name '*.faa' -print0 | sort -z | ALIGN_OUT="$out" xargs -0 -P "$jobs" -n 1 sh -c '
    f=$1; o=$ALIGN_OUT/$(basename "$f" .faa).aln.fasta
    [ -s "$o" ] && exit 0
    if mafft --auto --thread 1 --quiet "$f" > "$o.tmp$$"; then mv "$o.tmp$$" "$o"
    else rm -f "$o.tmp$$"; echo "mafft failed on $f" >&2; exit 1; fi' _
echo "aligned $(find "$out" -maxdepth 1 -name '*.aln.fasta' | wc -l) of $n families into $out" >&2
