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

**Not covered.** `Rog-ecdysozoa1000.rates` and `Chromo23.rates` come from `calculateRASRates` (rate variation): classes indexed by *loss level*
and *gain level* with their own probabilities and rate multipliers (for example two loss classes with probabilities 0.20 / 0.80 and loss
multipliers 0.40 / 1.00, one gain class). That is the separate gain and loss rate factors of SPEC §5.6 (open item 1), and a free-class model
rather than the equal-weight gamma classes of §5.5b, so those two files cannot be compared until that is implemented.
