# Sourced by run_variant*.sbatch: VARIANT selects the alignment, model and extra IQ-TREE options.
case "$VARIANT" in
  pmsf)      ALN=supermatrix.fasta;          MODEL="LG+C60+F+G4"; EXTRA="-ft hyp/free.treefile" ;;   # site-heterogeneous mixture, PMSF with the LG+F+G4 ML tree as guide
  fast25)    ALN=variants/fast25.fasta;      MODEL="LG+F+G4";     EXTRA="" ;;                        # the fastest 25 % of sites removed
  fast50)    ALN=variants/fast50.fasta;      MODEL="LG+F+G4";     EXTRA="" ;;                        # the fastest 50 % removed
  fast75)    ALN=variants/fast75.fasta;      MODEL="LG+F+G4";     EXTRA="" ;;                        # the fastest 75 % removed: the slowest 50,000 sites kept
  pmsf50k)   ALN=variants/fast75.fasta;      MODEL="LG+C60+F+G4"; EXTRA="-ft hyp/free.treefile" ;;   # PMSF C60 on the slowest 50,000 sites (~1/4 of the memory of the full matrix: ~61 GB instead of ~245 GB)
  dayhoff6)  ALN=variants/dayhoff6.fasta;    MODEL="MK+G4";       EXTRA="-st MORPH" ;;              # Dayhoff-6 recoding
  *) echo "unknown VARIANT '$VARIANT'" >&2; exit 2 ;;
esac
