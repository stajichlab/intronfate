# F2: how much missing data could column ambiguity resolve?

Question: SPEC §5.6b lists "column ambiguity" (an intron beside an indel could sit in any column of the gap run; give it weights instead
of masking the cell) as the next step of probabilistic projection. Is there enough such missing data to justify it?

`missing_cells.py` classifies every `?` cell of the V6 tier A matrix (nine genomes, 1,106 strict single-copy families, `--flank 1`; every
taxon is annotated, so `?` means "not aligned around the site") by the gap run in that taxon's row at the site. `run.sbatch` runs it on the
HPCC (job 32216850, seconds). `results/missing_cells.txt` has the output.

## Result (job 32217540)

141,273 cells, 17,251 missing (12.2%):

| gap at the site | cells | share of missing | share of all cells |
|---|---|---|---|
| short (run of at most 6 columns, residues on both sides) | 3,336 | 19% | 2.4% |
| long run | 11,255 | 65% | 8.0% |
| row begins or ends inside the window | 2,651 | 15% | 1.9% |
| no gap in the window | 0 | 0% | 0% |
| taxon has no row | 9 | 0.1% | 0.0% |

Only the first row is what column ambiguity could resolve, and only 1,195 of 15,697 sites (7.6%) have a carrier and such a cell. Long gaps
(most of the missing data) are regions with no homologous sequence in that taxon: they are correctly `?`, and no choice of column helps.

**Conclusion: not worth building now.** At most about 2.4% of cells would change, and the rule that decides them (how likely an intron at one
end of an insertion is the same as one at the boundary across the gap) has no data behind it. Revisit if a data set with many short indels
(divergent, poorly aligned families) shows a much larger share; rerun this script on it first.

Every missing cell has a gap in its window, so the missing data is fully explained by alignment gaps and row ends. (The first run of this
script read the matrix `column` as 1-based; it is 0-based, which gave 9% "no gap" cells and a few short-gap counts too low. Fixed.)
