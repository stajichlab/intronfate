# B1: speed of intronfate vs MALIN (JOSS paper evidence)

Purpose: the paper's speed claim. V3 timed `intronfate` against its own earlier versions only; this folder times it against
MALIN's command-line programs on the same matrix and tree. **Purgeable:** `results/` and `logs/` can be deleted freely.

| Task | intronfate | MALIN |
|---|---|---|
| Dollo parsimony | `score -m dollo` | `DolloParsimony` |
| Likelihood, constant rates | `score -m likelihood` | `calculateRASRates` (default: no rate variation) |
| Table export | `convert --to-malin` | n/a |

Not compared: **projection** (annotation to alignment columns). MALIN has no command-line step for it (its tables come
from the GUI), so it is timed for `intronfate` alone ([V3](../V3_gff_scale/README.md)). Say so in the paper rather than
quote a projection speed-up.

Fairness notes (state these in the paper):
- Convergence rules differ (`-opt-eps` vs `--ml-tol`/`--ml-param-tol`), so a wall-time ratio alone is not a claim. Compare the
  final log-likelihoods too (`results/if_lik.1/`, `results/malin_rates.1.out`) and report time to a similar likelihood.
- The JVM start-up is included in MALIN's time (it is part of using it); `MALIN_HEAP` sets `-Xmx` (default 8g).
- Same node, runs one after the other, `REPS` repeats; medians and ranges in `results/summary.md`.
- `intronfate` thread count is the tool's default; MALIN's CLI is single-threaded.

## Run (HPCC only)

```sh
cd ~/projects/intronfate_main/bench/B1_malin_timing
mkdir -p logs
sbatch --export=ALL,SMOKE=1 run.sbatch        # tiny fixtures; MALIN skipped unless MALIN_JAR is set
sbatch --export=ALL,MATRIX=<OUT>/combined/matrix.tsv,TREE=<species.nwk>,MALIN_JAR=<Malin.jar> run.sbatch
```

Real data: the V6 tier A matrix (1,106 families, nine species), `bench/V6_orthology_nine_species/results_projection/A_strict/combined/matrix.tsv`
with `species.nwk`. Read those paths in place; this folder writes only under its own `results/`.

`Malin.jar` is third-party and git-ignored; it is not in this repository.

## Tests

`python3 -I test_summarize.py` (summary statistics, failed runs excluded) and `python3 -I test_make_smoke_table.py`. `SMOKE=1` runs the whole pipeline on a generated 6-taxon table (`make_smoke_table.py`);


## Results

Job 32512585 (2026-11-27), HPCC node c05 (AMD Opteron 6376), `intronfate` 0.1.0 at `main` `2c8aee4`, MALIN `Malin.jar` as
supplied by the owner, Java 21.0.7 (`-Xmx8g`). Input: V6 tier A combined matrix (1,106 families, 9 species, 15,716 table lines),
`species.nwk`; 3 repeats; raw rows in `results_V6A/timings.tsv` on the HPCC.

| tool | task | median wall (s) | min-max (s) | peak RSS (MB) |
|---|---|---|---|---|
| intronfate | `score -m dollo` | 0.17 | 0.17-0.19 | 37 |
| MALIN | `DolloParsimony` | 0.58 | 0.57-0.59 | 81 |
| intronfate | `score -m likelihood` | 1.50 | 1.49-1.54 | 30 |
| MALIN | `calculateRASRates` | 90.75 | 89.72-90.95 | 175 |
| intronfate | `convert --to-malin` | 0.14 | one run | 28 |

## Likelihood comparability check (2026-11-27): NOT like for like, cause unknown

`run_lnl.sbatch` fits MALIN (`calculateRASRates`) on the same table (the tree needs named internal nodes: `label_internal_nodes.py`),
then evaluates intronfate's log-likelihood at MALIN's printed parameters with the L1 test (`malin_rates.rs`), next to intronfate's own
maximum on the same positions. Jobs 32512872, 32512893, 32512945; outputs `lnl_*.txt` on the HPCC.

| MALIN run | positions scored | intronfate lnL at its maximum | intronfate lnL at MALIN's parameters | gap |
|---|---|---|---|---|
| default (all columns) | 15,697 | -46,709.99 | -71,581.81 | 24,872 |
| `-ambiguous 0` (complete columns only) | 10,926 | -37,224.28 | -58,020.57 | 20,796 |
| `-ambiguous 0.5` | 14,031 | -45,668.16 | -69,872.33 | 24,204 |
| `-ambiguous 0`, `-opt-round 5000` | 10,926 | -37,224.28 | -58,018.80 | 20,795 |

**Rogozin control (job 32513341, `results_rog/`).** The same harness on the Rogozin table and tree, where L1 agrees with MALIN's
*published* `.rates` to 0.002, gives lnL -17,817.36 at intronfate's maximum but **-30,984.14 at the parameters this harness got
from `calculateRASRates`** (default switches, `-ambiguous 0`, 7.7 s). So the fault is in how the harness obtains MALIN's fit, not in the
V6 data and not in the model: the fitted internal branches sit at about 1e-5 (for example Diptera, Ecdysozoa), unlike the published file
(0.004, 0.32). The jar has no `calculateConstantRates` class (L1's published file came from that program or the GUI), only
`calculateRASRates`, whose CLI defaults reach a different, much worse optimum under this likelihood. Until MALIN's CLI is made to reproduce
its published fit on Rogozin (try `-start-with` the published rates, a tighter `-opt-eps`, and check that `-opt-round` reaches the program;
the switch syntax is in the guide, Appendix B.3), none of the MALIN timings above is a like-for-like number.

Reading: the gap is the same size with complete columns only (so missing-data handling, per-pattern vs global never-seen, is not
the cause; both modes agree there) and does not move when MALIN is allowed 5,000 rounds (so it is not early stopping). MALIN
fits one class. On the Rogozin table the same test agrees to 0.002 (L1), so the two tools do not maximise the same function on *this*
table, or the exported table or its tree is read differently by MALIN. Not yet diagnosed. Candidates, none tested: how MALIN reads the
exported phase codes and `#MAP` block for this table; a root with a tip child (`Spunctatus`) and MALIN's resolution of the root split;
MALIN's treatment of columns with an intron in one genome only. Cheapest next step: run this harness on the Rogozin table, where L1
already agrees, to tell a harness fault from a data-specific difference.

MALIN wall time depends strongly on its options (default 90 s; `-ambiguous 0` 23 s; `-ambiguous 0.5` 146 s; `-opt-round 5000` 108 s), so
no speed ratio is quotable until the comparison above is explained.

Open before any of this is quoted:
- **The likelihood ratio (about 60x) is withdrawn pending the check above.** intronfate's `ml_log_likelihood` is -46709.99 (default
  settings, tolerances 1e-6/1e-7). MALIN's `calculateRASRates` output carries no log-likelihood line, so the two fits have not
  been compared. Next step: score MALIN's `.rates` with `intronfate` (as [L1](../L1_malin_likelihood/README.md) did) or run
  `computeRASPosteriors`, and compare the likelihood each stopping rule reaches; also try `-opt-eps` looser and tighter.
- The models differ in detail (MALIN fits branch lengths and a root density; intronfate fixes the root and treats the
  error rates as inputs, SPEC §5.5), so say what is fitted by each.
- Dollo is under a second for both: report it as comparable, not as a win.
- One data set, one node, nine species. A scaling series in taxa and sites is still to do (HANDOFF item 8).
