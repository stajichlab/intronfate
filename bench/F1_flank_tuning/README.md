# F1: what does `--flank` do, and should the default (1) change?

MALIN's test data (`Chromo23.table`, Rogozin) holds finished tables, not alignments, so the flank rule cannot be compared with MALIN's own
masking. This test uses the V6 nine-species fungal families (tier A, 1,106 single-copy families, MAFFT alignments).

Purgeable on the HPCC: `results/` is small and kept. Run: `mkdir -p logs results && sbatch run.sbatch` (`FLANKS="0 1 2"` limits the sweep).

| File | What |
|---|---|
| `run.sbatch` | `analyze --combine --combined-only` at flank 0, 1, 2, 3, 5, 8, then the two scripts below |
| `gap_distance.py` | for conserved sites (present in at least two other taxa) the absence rate of a taxon as a function of the number of residue columns between the site and the nearest gap in that taxon's row; flank-0 matrix |
| `summarize.py` | masked presences and Dollo totals per flank |
| `results/` | `summary.txt`, `absence_by_gap_distance.tsv` (+ `.by_taxon`), per-flank `dollo_branches` and matrix headers |

## Result (jobs 32196757, 32196811)

| flank | sites kept | presences masked | Dollo losses |
|---|---|---|---|
| 0 | 16,504 | 584 | 15,101 |
| 1 (default) | 15,697 | 2,059 | 14,316 |
| 2 | 14,442 | 4,288 | 13,142 |
| 3 | 13,447 | 6,148 | 12,194 |
| 5 | 11,756 | 9,363 | 10,555 |
| 8 | 9,793 | 13,218 | 8,590 |

Each step of flank discards sites and turns annotated introns into `?`. The question is whether what it removes is unreliable. If gaps near an
intron produced false losses, the absence rate of conserved introns would be raised close to a gap. **It is not**: within every taxon the rate
is flat from distance 0 to 8 or more (for example *N. crassa* 0.86 / 0.85 / 0.83 / 0.84 / 0.87, *C. cinerea* 0.39 / 0.40 / 0.36 / 0.39 / 0.35;
`absence_by_gap_distance.tsv.by_taxon`). Pooled over taxa it is 0.70–0.72 at every distance, which is dominated by the taxa that have lost most
introns and is why the table is also given per taxon.

**Conclusion: keep 1.** Nothing in these data supports a larger flank (it only discards data), and flank 1 is what MALIN's rule needs for a
site inside a residue. Flank 0 would keep 1,475 more presences, but the rule would no longer need a residue on both sides of a boundary intron
beyond the one-column minimum the code already applies, so it is the least conservative choice rather than a better one.

## Caveats

- 353 of 8,377 present phase-1/2 cells (4%) fall on a gap at the column read from `family_column` (clustered: the same site in 2–3 taxa, e.g.
  OG0001922 column 2219). At flank 0 a phase-1/2 carrier needs a residue there, so the cause is not understood; the distances for those sites
  are therefore uncertain. It does not change the flat trend (the 96% that index correctly show it) but should be traced.
- One data set, nine genomes, strict single-copy families; the result may differ for sparser or more divergent alignments.
