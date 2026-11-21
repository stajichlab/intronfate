# S1: does soft junction evidence help the convergent-insertion calls? A simulation with known truth

[J2b](../J2b_matrix_cells/README.md) found too little RefSeq disagreement in the reliable genomes to test the soft-evidence likelihood against RefSeq. Here the truth is simulated, at any error rate.

## Design (`simulate.py`, `evaluate.py`; selftests in `bench/test_support_vs_reference.py`; job `run.sbatch`)

- **Sites and truth.** The real strict-family matrix (J1, 15,697 sites, its missing-data mask) is kept. Presence at the nine tips is simulated on the species tree with the
  two-state model the likelihood fits, with the J1 constant-rate hard fit's parameters (root presence, branch gain and loss probabilities). Every site's number of origins
  (root presence + 0 → 1 changes) is known. A site is kept only if an observable tip is *annotated* present (the model's ascertainment).
- **Annotation error.** A true intron is missed with probability `miss`; a true absence is annotated present at a per-taxon rate chosen so false introns are a set share of that
  taxon's annotated-present cells, over all 109,669 positions (matrix sites + the fit's never-seen positions): none (0, 0%), low (5%, 1%), mid (10%, 3%), high (15%, 8%).
- **Reads** (six genomes; the other three stay hard). `empirical`: a cell's read count is drawn from the distribution among J2b cells that agree (true) or disagree (false) with
  RefSeq. This is pessimistic for soft evidence, since some disagreeing cells are real introns, so false introns look better supported than they are. `clean`: a false intron has
  one read with probability 0.05. True cells use the empirical distribution in both.
- **Fits** (constant rates, 20,000 iterations, 3 seeds each): `hard` (no error model); `hard_err` (hard, given the **true** miss and false rates: an oracle competitor);
  `soft_k3` (`junctions` defaults, p0 0.8, k 3); `soft_k0.2` (p0 0.8, k 0.2); `soft_p0.5_k0.2` (p0 0.5, k 0.2: a map that separates unsupported from supported cells more).
- **Scoring.** AUC of ln LR for the sites with ≥ 2 true origins; the calls at p(≥ 2 origins) > 0.5 (count, precision, recall); precision of the top 100 by ln LR; AUC of p(no origin)
  for sites with no true presence anywhere.

## Results (mean of 3 seeds; `results/S1_runs.tsv`, `results/S1_summary.tsv`)

AUC of ln LR for convergent insertion (about 650–720 of 15,697 sites truly have ≥ 2 origins), `empirical` reads:

| error scenario | false-only sites | hard | hard + true error rates | soft k 3 (default) | soft k 0.2 | soft p0 0.5, k 0.2 |
|---|---|---|---|---|---|---|
| none | 0 | 0.828 | 0.828 | 0.762 | 0.821 | 0.799 |
| low | 200 | 0.804 | 0.805 | 0.731 | 0.795 | 0.771 |
| mid | 612 | 0.769 | 0.770 | 0.706 | 0.761 | 0.741 |
| high | 1,581 | 0.743 | 0.745 | 0.697 | 0.738 | 0.722 |

- **Soft evidence did not improve the convergent-insertion ranking in any scenario.** Against the hard fit the AUC was lower in every seed, by 0.05–0.08 for the default (k 3), 0.00–0.01
  for k 0.2 and 0.02–0.04 for p0 0.5. It is worst when there is no error at all, and the gap narrows as the error rate rises but never closes. `clean` reads give the same ordering
  (the soft k 0.2 fit ties the hard fit only in the high-error scenario).
- **Soft fits make no calls.** At p(≥ 2 origins) > 0.5 the soft fits flagged nothing in any scenario, while the hard fit flags 96–174 sites, of which 37–70% are truly convergent (precision
  falls as errors rise: 0.70 with none, 0.37–0.56 at mid and high). With `empirical` reads the top-100 precision is 0.13–0.47 for soft fits against 0.41–0.75 for hard.
- **The error model that knows the true rates is the better remedy, but is no better at ranking.** `hard_err` has the same AUC as `hard`, makes fewer and more precise calls
  (0.45–0.67 at mid and high against 0.37–0.56; it makes only 8 calls in the high `empirical` scenario), and is the only fit that can say a site's presences are all errors (AUC 0.74–0.76 for no-origin sites; hard and soft fits have no
  mechanism to do so: with `--ml-soft-evidence` and no false-intron rate the seen-conditioning forces some carrier to be real, so p(no origin) is 0 and its AUC 0.5 by construction).
- **Why soft evidence loses here.** A soft cell is only partly present, so any site with two or more carriers has a lower chance that both are real; that lowers p(≥ 2 origins) for the
  truly convergent sites as much as for the rest. The reads would help only if they separated false from true cells strongly (the reads-given-error model), but the default map gives an
  unsupported cell 0.8 and a one-read cell 0.85, and k 0.2 gives 0.8 and 0.97: a gap too small to outweigh the discount on every true intron. Tuning to p0 0.5 widens the gap and
  still loses 0.02–0.03 AUC against the hard fit.

## Limits (read before quoting any number)

- Constant rates and independent random errors: false introns are placed uniformly, not as shifted splice sites beside a true intron, and real errors correlate with gene structure and
  expression. The fitted model is the generating model.
- The read models are assumptions anchored on J2b, not measurements; the no-read share of true cells (6%) and the false-cell read distribution drive the result.
- Truth is "≥ 2 origins on this tree", which is hard to detect (many convergent sites look like one origin), so even the oracle ranks at AUC ≈ 0.8.
- The soft k 3 fits at high error with `clean` reads did not converge in 20,000 iterations (0 of 3; every other cell converged except one of three soft k 0.2 fits at low error), so those rows are lower bounds on its fit.
- Three seeds per cell and means only: the differences quoted are consistent across all seeds, but there are no intervals. Ranking is evaluated; the effect on gain and loss totals is not.
- This says nothing against soft evidence for other uses (flagging suspect cells, a fixed per-cell confidence from a better source than read counts); it says the read-count map as built, with
  the default or any tuned setting tried here, does not improve the convergent-insertion calls and can remove them.

## Consequence for the defaults

The `junctions` defaults stay as they are; there is no case for changing them from this. The case for presenting the soft likelihood as an improvement is not supported: it should be
labelled as a sensitivity analysis (it shows which calls depend on trusting the annotation), as in J1.
