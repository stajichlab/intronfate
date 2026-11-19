# L1: the likelihood model against MALIN's own estimates

**Question.** Does `-m likelihood` (SPEC §5.5, one rate class, no annotation error) maximise the same likelihood as MALIN's
`calculateConstantRates`? MALIN's published output for the Rogozin table and the Ecdysozoa tree is `Rog-ecdysozoa.rates`
(MALIN's licensed test set; not committed, `Malin/` is git-ignored).

**Method** (`crates/intronfate-core/tests/malin_rates.rs`, an ignored test; `run.sbatch` runs it on the HPCC, 39 s): fit this tool's model to
`rogozin_introns_table` on the tree in the `.rates` file (strict tolerance 1e-10), evaluate this tool's log-likelihood **at MALIN's printed
parameters** (per-branch `p(gain)`, `p(loss)` and the root presence probability, parsed by `malin::parse_rates`), and print both parameter sets.
Evaluating MALIN's parameters under our likelihood does not depend on how the non-identifiable root split is resolved.

**Result** (`results_rog.txt`, 2026-11-18). 7,221 sites, 12 branches.

| | log-likelihood |
|---|---|
| this tool at its maximum (400 E-steps, converged) | -17817.3633 |
| this tool at MALIN's parameters | -17817.3652 |
| difference | +0.0019 |

MALIN's parameters are within 0.002 log-likelihood units of this tool's maximum, so the two tools maximise the same function, and MALIN's
estimates sit at our optimum (the likelihood is flat in the direction where they differ). Per-branch probabilities agree to three or four
digits on every branch that is separately identifiable (for example Dm gain 0.00291 / 0.00291, loss 0.32189 / 0.32191; Diptera 0.00408 /
0.00407, 0.58176 / 0.58182; Ecdysozoa 0.00368 / 0.00369, 0.31913 / 0.31901). They differ on the two root branches (Crown and Pf) and the root
presence probability (ours 0.0253, MALIN 0.0972): with two root children the root probability and the two branches are not separately
identifiable (SPEC §5.5), this tool fixes the larger child's branch to "no change", and MALIN resolves the split differently; the likelihood
is the same.

## Rate variation (free levels)

`Rog-ecdysozoa1000.rates` (Rogozin table) and `Chromo23.rates` (Chromo23 table, 23 taxa) come from `calculateRASRates`: classes indexed by loss
level and gain level, each level with its own probability and multiplier. Both have two loss levels and one gain level. This tool's free-level model
(SPEC §5.5d, `--ml-loss-levels 2`) is compared with them by `run_rates.sbatch` (`results_rates_3.txt`), using the same test as above
(`INTRONFATE_MALIN_MIN_OBSERVED` selects MALIN's table filter, `INTRONFATE_MALIN_NEVER_SEEN=global` MALIN's never-seen correction; `results_rates_4.txt` is the global run).

| | Rogozin | Chromo23 |
|---|---|---|
| positions | 7,221 | 7,030 (the 9,256-site table filtered to positions observed in at least 18 taxa; MALIN reports 7,030) |
| loss levels, MALIN (probability @ multiplier) | 0.798 @ 1, 0.202 @ 0.402 | 0.845 @ 1, 0.155 @ 0.387 |
| loss levels, this tool, MALIN's never-seen correction | 0.815 @ 1, 0.185 @ 0.401 (root branch free, started at MALIN's parameters) | 0.8446 @ 1, 0.1554 @ 0.3876 |
| lnL at MALIN's parameters | -17796.1286 | **-28430.9249** (MALIN prints -28430.9249; uncorrected -40617.5994, MALIN -40617.5994; p0 0.823339, MALIN 0.82334; never-seen 32,763.79, MALIN 32,763.79) |
| lnL at this tool's maximum | -17796.0627 (+0.066; root branch free) | -28430.9159 (+0.009; the same with the larger root branch fixed or free) |
| lnL at this tool's maximum, larger root branch fixed (default) | -17796.7474 (-0.62) | -28430.9159 (no cost) |

On both tables MALIN's parameters sit at this tool's optimum (within 0.07 and 0.01 log-likelihood units, with the same two-level structure), and on `Chromo23`
the likelihood this tool computes at MALIN's parameters equals the value MALIN prints to four decimals, including its uncorrected value and its
estimate of the never-seen positions. Two details matter. (1) The never-seen correction: `Chromo23` has missing cells, and MALIN corrects the whole table with one
probability `p0` that a position has no intron at any genome (log-likelihood = uncorrected - N ln(1 - p0); its "predicted missing" is N p0 / (1 - p0)); this tool's
default conditions each position on its own pattern of missing cells, which gave an optimum 26 log-likelihood units above its value at MALIN's parameters and a
different quantity from MALIN's printed -28430.9. With `--ml-never-seen global` the two coincide; without missing cells (Rogozin) there is no difference. (2) The larger
root branch: fixing it (the default) cost 0.7 units on the Rogozin table once the classes scale the two root branches differently (SPEC §5.5d) and nothing on `Chromo23`.

## Rate classes on the Rogozin table (`run_classes.sbatch`, `results_classes.txt`)

| model | log-likelihood |
|---|---|
| no rate classes | -17817.36 |
| shared factor, 2 / 4 gamma classes | -17797.64 / -17795.83 |
| loss factor only, 2 / 3 / 4 gamma classes | -17800.89 / -17798.97 / -17797.91 |
| gain factor only, 2 classes | -17817.20 (not converged; the shape drifts to 0.16) |
| gain 2 x loss 2 / gain 3 x loss 3 | -17800.76 / -17798.92 (not converged in 50,000 E-steps) |
| free levels, 2 loss levels (MALIN's structure) | about -17796.7 (root branch fixed), -17796.1 (free) |

Variation in the gain rate adds nothing here (a gain x loss grid improves on loss-only by under 0.2), which matches MALIN's choice of two loss levels and one gain
level; with the same number of gamma classes a shared factor fits slightly better than a loss-only factor.
