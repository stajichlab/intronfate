# V7: how much of the Dollo gains/losses could annotation error explain?

V4 gave a worst-case ceiling (`sites_carried x (1 - Sn)/Sn`). V7 simulates the error instead: it perturbs each V6 tier's **combined** matrix (all families in one
matrix, `analyze --combine`) at each genome's measured sensitivity (Sn) and precision (Pr) over the tier's own genes (`../V4_annotation_accuracy/results_v6/<tier>/busco_accuracy.tsv`),
re-scores with `intronfate score -m dollo`, and reports the per-branch Dollo gains/losses over 100 replicates per mode (mean and the 2.5th/97.5th percentiles).

| File | Role |
|---|---|
| `perturb.py` | the error model (below); `tests/test_perturb.py` is its fixture test (`python3 tests/test_perturb.py`) |
| `run_sensitivity.py` | one tier: score the observed matrix (asserted equal to the tier's own Dollo table), perturb, re-score in `$SCRATCH`, write `results/<tier>/{replicates,summary}.tsv` |
| `run.sbatch` | the three tiers (A_strict, B_partial, C_hog); `mkdir -p logs && sbatch -p highclock run.sbatch` |

**Model.** Missed introns: a taxon with `P` introns in the matrix, precision Pr and sensitivity Sn is missing `M = P*Pr*(1 - Sn)/Sn` (true positives `P*Pr`; FN = TP(1-Sn)/Sn).
The budget is the genome's, so it is pooled over all families. A miss matters to the Dollo score only where it hides a loss, so it is placed on a random *candidate* cell and
set to `1`; the candidates are the cells where the Dollo scorer put a loss on the taxon's own terminal branch (read from the observed run's `dollo_sites.tsv`, so the
rule, including `?` cells and ties, is the scorer's own). `ceiling_fn` fills min(M, candidates) cells (every miss hides a loss); `typical_fn` scales M by the taxon's share
of introns at sites with >= 2 carriers. False introns (`fp`): every `1` becomes `0` with probability 1 - Pr. `?` cells are never changed.
**Correction.** The first version of this analysis capped the budget per family, which discarded misses in families where a taxon had no candidate cell and understated the
ceiling (at the ceiling it removed only 213 of *C. cinerea*'s 514 losses, where all 514 can be explained); an independent code review found it. The first version also drew random cells that were
not terminal losses, which created sharing with distant lineages and raised losses; both are fixed.

## Results (tier A, 1,106 single-copy families; B and C in `results/*/summary.tsv`)
Terminal-branch Dollo losses, observed vs. after undoing the missed introns:

| Genome | Observed | Misses undone (ceiling) | Misses undone (typical) |
|---|---|---|---|
| *S. maydis* | 4,055 | 4,015 | 4,030 |
| *S. cerevisiae* | 1,447 | 1,446 | 1,446 |
| *C. neoformans* | 1,091 | 1,024 | 1,046 |
| *S. pombe* | 1,011 | 965 | 977 |
| *N. crassa* | 628 | 588 | 597 |
| *M. mucedo* | 1,516 | 1,289 | 1,351 |
| *C. cinerea* | 514 | 0 | 140 |
| *A. nidulans* | 151 | 0 | 0 |

- Losses in *S. maydis*, *S. cerevisiae*, *C. neoformans*, *S. pombe* and *N. crassa* move by at most 6 % even at the ceiling: they are not artefacts of missed introns.
- *M. mucedo* moves by 15 %. For *C. cinerea* and *A. nidulans* the ceiling removes **all** their terminal losses, so those counts cannot be told from annotation misses,
  as V4 found (the simulation now agrees with V4's analytic ceilings). Tiers B and C give the same ordering.
- The replicates differ little (each fill removes exactly one loss, so the spread comes only from which cells are drawn), so the percentile intervals are narrow.
- False introns (`fp`): terminal gains move by -10 % to +17 % with no consistent direction. Removing an intron at a shared site can create a loss or move a gain, and the
  uniform removal is a crude model (a false intron is unlikely to coincide with an orthologous site), so read gains with V4's precision numbers, not this table.
  Do not read the `ceiling` column (both corrections) for losses; its false-intron half raises losses by deleting carriers.

Limits: Sn and Pr are RefSeq-relative and measured at V2/V6 genes; missed introns are assumed to hide terminal losses or to be distributed like annotated ones; independent draws per
taxon; *S. cerevisiae* Sn/Pr rest on 39-46 introns in tier A (noisy; its losses barely move either way); in tier C *S. punctatus* has no genes, so it is uncorrected there.
