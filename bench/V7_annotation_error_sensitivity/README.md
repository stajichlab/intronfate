# V7: how much of the Dollo gains/losses could annotation error explain?

V4 gave a worst-case ceiling (`sites_carried x (1 - Sn)/Sn`). V7 simulates the error instead: it perturbs every V6 matrix at each genome's measured
sensitivity (Sn) and precision (Pr) over the tier's own genes (`../V4_annotation_accuracy/results_v6/<tier>/busco_accuracy.tsv`), re-scores with
`intronfate score -m dollo`, and sums the per-branch Dollo gains/losses over families, 20 replicates per mode.

| File | Role |
|---|---|
| `perturb.py` | the error model (below); `tests/test_perturb.py` is its fixture test (`python3 tests/test_perturb.py`) |
| `run_sensitivity.py` | one tier: perturb, re-score in `$SCRATCH`, write `results/<tier>/{replicates,summary}.tsv` |
| `run.sbatch` | the three tiers (A_strict, B_partial, C_hog); `mkdir -p logs && sbatch -p highclock run.sbatch` |

**Model.** Missed introns: a taxon with `present` introns and sensitivity Sn is missing `M = present (1 - Sn)/Sn`. A miss matters to the Dollo
score only where it hides a loss, so it is placed on a random cell that is a terminal-branch loss (the taxon is `0`, lies inside the clade of the site's
carriers, and a relative in its parent clade carries the site) and set to `1`. `ceiling_fn` fills min(M, candidates) cells (every miss hides a loss);
`typical_fn` scales M by the taxon's share of introns at sites with >= 2 carriers. False introns (`fp`): every `1` becomes `0` with probability 1 - Pr.
`?` cells are never changed. (A first version placed misses on arbitrary `0` cells; that creates sharing with distant lineages and pushes gains
deeper, which raised losses and was wrong for this purpose.)

## Results (tier A, 1,106 single-copy families; B and C in `results/*/summary.tsv`)
Terminal-branch Dollo losses, observed vs. after undoing the missed introns (mean and 2.5-97.5 % of 20 replicates):

| Genome | Observed | Misses undone (ceiling) | Misses undone (typical) |
|---|---|---|---|
| *S. maydis* | 4,055 | 4,020 [4,013-4,030] | 4,033 |
| *S. cerevisiae* | 1,447 | 1,446 | 1,446 |
| *C. neoformans* | 1,091 | 1,047 [1,031-1,062] | 1,063 |
| *S. pombe* | 1,011 | 986 | 994 |
| *N. crassa* | 628 | 612 | 617 |
| *M. mucedo* | 1,516 | 1,347 [1,328-1,359] | 1,386 |
| *C. cinerea* | 514 | 301 [286-314] | 349 |
| *A. nidulans* | 151 | 103 [96-108] | 114 |

- Losses in *S. maydis*, *S. cerevisiae*, *C. neoformans*, *S. pombe* and *N. crassa* move by 0-4 % even at the ceiling: they are not artefacts of missed introns.
- *M. mucedo* moves by 9-11 %. *C. cinerea* (-32 to -41 %) and *A. nidulans* (-25 to -32 %) are the ones annotation error can move a lot, consistent with V4's
  "cannot be told apart" for those two. Tiers B and C give the same ordering.
- Missed introns never fully explain a loss count here, because only cells that look like a terminal loss can be filled, which is fewer than V4's ceiling assumed.
- False introns (`fp`): terminal gains move by -10 % to +17 % with no consistent direction. Removing an intron at a shared site can create a loss or move a gain, and the
  uniform removal is a crude model (a false intron is unlikely to coincide with an orthologous site), so read gains with V4's precision numbers, not this table.
  Do not read the `ceiling` column (both corrections) for losses; its false-intron half raises losses by deleting carriers.

Limits: Sn and Pr are RefSeq-relative and measured at V2/V6 genes; missed introns are assumed to hide terminal losses or to be distributed like annotated ones;
independent draws per family and taxon; Sn/Pr of the outgroup *S. punctatus* are `nan` (not scored by V4 at these genes), so it is left uncorrected.
