# B2: scaling of import and scoring in taxa and sites

Question: how do `convert --from-malin` (import), `score -m dollo`, `score -m mp` and `score -m likelihood` scale as the number of taxa and
sites grow? V3 measures projection (GFF3 to protein space) in genomes; this covers everything after projection (HANDOFF item 8). It is
intronfate only: the MALIN likelihood comparison is parked ([B1](../B1_malin_timing/README.md)).

**Data.** Synthetic, deterministic (`make_scaling_data.py`, seed 1): a random binary tree on T taxa; each site has one gain at a random
node and independent loss (p = 0.1) on every branch below it. Written as a MALIN table with `#MAP` lines so the import path is the one
users take. Grid: taxa 10, 30, 100, 300, 1,000 by sites 1,000, 10,000, 100,000. Timings use `/usr/bin/time` (wall, peak RSS), 3 repeats,
median reported; failures (124 = timed out after `CELL_TIMEOUT`) are shown, never averaged in. The likelihood is skipped above
`LIK_MAX_CELLS` (2e7 matrix cells) so one slow cell cannot eat the allocation; raise it to extend the grid.

```sh
mkdir -p logs && sbatch --export=ALL,SMOKE=1 run.sbatch    # seconds: 6-12 taxa, 30-60 sites
sbatch run.sbatch                                          # the full grid (results/summary.md)
python3 -I test_make_scaling_data.py && python3 -I test_summarize.py
```

Purgeable: `results/` and `logs/`. Run on the HPCC, one node type, nothing else on the node.

**What it can and cannot say.** Simulated sites are cleaner than real ones (no alignment gaps, no `?` cells, balanced loss rate), so the
numbers show how cost grows with taxa and sites, not what a real data set costs; the real-data anchors are V3 (projection), V6 tier A
(15,697 sites, 9 taxa: Dollo 0.17 s, likelihood 1.5 s) and L1/V1. The run is single node, one thread requested; any speed-up from the
parallel parts is not measured.

## Results

Not yet run.
