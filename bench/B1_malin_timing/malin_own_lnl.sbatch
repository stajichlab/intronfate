#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 00:20:00 -p batch
#SBATCH --job-name=B1_malin_lnl
#SBATCH --output=logs/B1_malin_lnl.%j.log

# B1: ask MALIN itself (computeRASPosteriors) for the likelihood of the Rogozin table under (a) its published rate files and
# (b) the rate files the CLI fitted in repro_rogozin.sbatch. If MALIN scores its CLI fit near -17817, the mismatch is in
# our reading of the file; if it scores it near -30984, the CLI optimiser stops at a bad optimum. Output: malin_own_lnl.txt.
set -euo pipefail
H=~/projects/intronfate_main; HERE=$H/bench/B1_malin_timing; R=$HERE/results_rog
L1=~/projects/intron_projector/bench/L1_malin_likelihood/data
module load java >/dev/null 2>&1 || true
P=ca.umontreal.iro.evolution.introns
: > "$HERE/malin_own_lnl.txt"
for f in "$L1/Rog-ecdysozoa.rates" "$L1/Rog-ecdysozoa1000.rates" "$R/rates.default.out" "$R/rates.eps_round.out" "$R/rates.publen.out"; do
  echo "=== $f" | tee -a "$HERE/malin_own_lnl.txt"
  java -Xmx4g -cp "$H/Malin/Malin.jar" $P.computeRASPosteriors "$R/matrix.malin.table" "$R/tree.labelled.nwk" "$f" 2>&1 \
    | grep -i -E "likelihood|lnL|log|never|no intron|predicted|Exception|rror" | head -8 | tee -a "$HERE/malin_own_lnl.txt" || true
done
