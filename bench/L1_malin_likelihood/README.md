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
(`INTRONFATE_MALIN_MIN_OBSERVED` selects MALIN's table filter).

| | Rogozin | Chromo23 |
|---|---|---|
| positions | 7,221 | 7,030 (the 9,256-site table filtered to positions observed in at least 18 taxa; MALIN reports 7,030) |
| loss levels, MALIN (probability @ multiplier) | 0.798 @ 1, 0.202 @ 0.402 | 0.845 @ 1, 0.155 @ 0.387 |
| loss levels, this tool | 0.815 @ 1, 0.185 @ 0.401 (root branch free, started at MALIN's parameters) | 0.849 @ 1, 0.151 @ 0.386 (root branch free, started at MALIN's parameters) |
| lnL at MALIN's parameters | -17796.1286 | -27985.6768 |
| lnL at this tool's maximum, root branch free | -17796.0627 (+0.066) | -27959.6099 (+26.07) |
| lnL at this tool's maximum, larger root branch fixed (default) | -17796.7474 (-0.62) | not run |

On the Rogozin table MALIN's parameters sit at this tool's optimum (within 0.07 log-likelihood units, with the same two-level structure), and the
default fit (larger root branch fixed) is 0.7 lower, the cost of that restriction once the root branches are scaled differently by the classes
(SPEC §5.5d). On `Chromo23` the level structure and most branch probabilities agree (for example Osat loss 0.02682 / 0.02700, Hsap 0.20106 / 0.19901) but
this tool's optimum is 26 units higher than its value at MALIN's parameters, with a few branches clearly different (Chlorophyta gain 0.036 / 0.053,
Alveolata loss 0.038 / 0.015), and MALIN's own reported log-likelihood (-28430.9) is a different quantity from this tool's (-27959.6). A likely cause is the treatment
of missing data: MALIN's output reports one global estimate of the never-seen positions ("predicted missing", p0 0.823), while this tool conditions each position
on its own pattern of missing cells (SPEC §5.5). That was not tested here (it would need a MALIN-style correction as an option). The Rogozin table has no missing
cells, and there the two agree.

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
