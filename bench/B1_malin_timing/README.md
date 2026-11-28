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

`python3 -I test_summarize.py` (summary statistics, failed runs excluded). `SMOKE=1` runs the whole pipeline on the fixtures in
`tests/data` and `examples/03_malin`.

## Results

Not run yet. TODO: table from `results/summary.md`, the matching log-likelihoods, node and CPU.
