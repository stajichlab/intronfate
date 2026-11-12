#!/usr/bin/env bash
# MAFFT every family (no trimming: gappy regions are handled by the tool's observability rule).
set -euo pipefail
cd "$(dirname "$0")"; . ./config.sh
module load mafft 2>/dev/null || true
mkdir -p alignments
ls families/*.faa | xargs -P "${SLURM_CPUS_PER_TASK:-4}" -I{} sh -c 'b=$(basename {} .faa); [ -s alignments/$b.aln.fasta ] || mafft --auto --thread 1 --quiet {} > alignments/$b.aln.fasta'
echo "alignments: $(ls alignments/*.aln.fasta | wc -l)"
