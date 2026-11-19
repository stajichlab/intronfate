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

## Result (jobs 32196757 for the sweep, 32217541 for the gap-distance table)

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
is flat from distance 0 to 8 or more (for example *N. crassa* 0.85 / 0.82 / 0.86 / 0.83 / 0.87, *C. cinerea* 0.35 / 0.40 / 0.34 / 0.40 / 0.35;
`absence_by_gap_distance.tsv.by_taxon`). Pooled over taxa it is 0.70–0.72 at every distance, which is dominated by the taxa that have lost most
introns and is why the table is also given per taxon.

**Conclusion: keep 1.** Nothing in these data supports a larger flank (it only discards data), and flank 1 is what MALIN's rule needs for a
site inside a residue. Flank 0 would keep 1,475 more presences (584 against 2,059 masked), but a boundary intron would still need residues on both sides, and nothing here shows the extra cells are better. The default stays at MALIN's rule with one flank column.

## Caveats

- One data set, nine genomes, strict single-copy families; the result may differ for sparser or more divergent alignments.
- Per taxon the rate varies between bins by a few points (e.g. *S. pombe* 0.82 at distance 0, 0.86 to 0.87 beyond; *N. crassa* 0.82 to 0.87) with 200–400 cells per near bin. If anything the rate is slightly *lower* right at a gap, the opposite of false losses.
- The first run of `gap_distance.py` read the matrix `column` as 1-based (it is 0-based), which showed up as 353 present cells "on a gap".
  Corrected in job 32217541: 0 of 8,377, and the trend above is from the corrected run.
