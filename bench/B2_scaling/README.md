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

## Results (job 32520847, node c09, AMD Opteron 6376, 3 repeats, 7 h 16 min in all)

Median wall time (peak RSS), taxa down, sites across. `fail124` = no finish in 1,800 s; `-` = skipped by the likelihood cap (above 2e7 cells).

#### dollo: median wall time (peak RSS), taxa down, sites across

| taxa \ sites | 1000 | 10000 | 100000 |
|---|---|---|---|
| 10 | 0.04 s (5 MB) | 0.11 s (14 MB) | 1.19 s (116 MB) |
| 30 | 0.03 s (6 MB) | 0.22 s (32 MB) | 2.12 s (295 MB) |
| 100 | 0.13 s (11 MB) | 0.60 s (72 MB) | 5.21 s (686 MB) |
| 300 | 0.28 s (38 MB) | 1.31 s (220 MB) | 12.74 s (2169 MB) |
| 1000 | 2.22 s (382 MB) | 4.95 s (587 MB) | 34.03 s (5844 MB) |

#### import: median wall time (peak RSS), taxa down, sites across

| taxa \ sites | 1000 | 10000 | 100000 |
|---|---|---|---|
| 10 | 0.03 s (5 MB) | 0.08 s (15 MB) | 0.63 s (119 MB) |
| 30 | 0.02 s (6 MB) | 0.14 s (26 MB) | 1.43 s (242 MB) |
| 100 | 0.05 s (11 MB) | 0.36 s (67 MB) | 3.71 s (647 MB) |
| 300 | 0.11 s (21 MB) | 0.85 s (184 MB) | 8.34 s (1812 MB) |
| 1000 | 0.27 s (62 MB) | 2.69 s (595 MB) | 25.74 s (5907 MB) |

#### likelihood: median wall time (peak RSS), taxa down, sites across

| taxa \ sites | 1000 | 10000 | 100000 |
|---|---|---|---|
| 10 | 1.88 s (5 MB) | 1.80 s (14 MB) | 29.19 s (116 MB) |
| 30 | 7.35 s (6 MB) | 101.77 s (32 MB) | 1215.52 s (296 MB) |
| 100 | 34.79 s (11 MB) | 335.59 s (72 MB) | fail124 |
| 300 | 39.20 s (39 MB) | 1302.82 s (220 MB) | - |
| 1000 | fail124 | fail124 | - |

#### mp: median wall time (peak RSS), taxa down, sites across

| taxa \ sites | 1000 | 10000 | 100000 |
|---|---|---|---|
| 10 | 0.04 s (7 MB) | 0.10 s (14 MB) | 0.88 s (116 MB) |
| 30 | 0.03 s (6 MB) | 0.22 s (32 MB) | 2.09 s (295 MB) |
| 100 | 0.11 s (11 MB) | 0.54 s (71 MB) | 5.39 s (686 MB) |
| 300 | 0.29 s (38 MB) | 1.42 s (220 MB) | 12.19 s (2170 MB) |
| 1000 | 2.21 s (380 MB) | 5.15 s (587 MB) | 37.16 s (5844 MB) |

Reading it:
- **Import, Dollo and parsimony are linear in cells.** Going from 10 to 1,000 taxa (100x) at 100,000 sites costs about 30-40x in time and 50x in memory; sites 1,000 to 100,000 costs about 15-30x. The largest cell (1,000 taxa x 100,000 sites = 10^8 cells) imports in 26 s and scores in 34-37 s, at 5.8 GB. Memory is the limit, about 58 bytes per cell.
- **Likelihood is the expensive method and its time is not set by size alone.** It takes 29 s for 10 taxa x 100,000 sites, 1,216 s for 30 x 100,000, and does not finish in 1,800 s for 1,000 taxa x 1,000 sites. The cost follows the number of EM iterations to converge as well as the cell count: 300 taxa x 1,000 sites (39 s) is faster than 100 x 1,000 (35 s) only by a little, yet 30 x 10,000 (102 s) is slower than 100 x 1,000. On these simulated sites (one gain per site, 10 % loss per branch, no real alignment structure) the fit is hard; real data is much quicker (V6 tier A, 9 taxa and 15,697 sites: 1.5 s). Treat these as an upper bound for the likelihood, not a typical cost, and state that the likelihood does not scale to hundreds of taxa and tens of thousands of sites at the default tolerance.
- Peak memory for the likelihood equals the Dollo figure, so memory does not grow with iterations.
- Not measured: threads, and the likelihood's per-iteration cost. If the paper claims the likelihood scales, this table says it does not past ~30 taxa x 100,000 sites; Dollo and parsimony are the scalable methods.
