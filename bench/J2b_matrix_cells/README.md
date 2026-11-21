# J2b: J2 restricted to the matrix cells, with a relaxed RefSeq match

[J2](../J2_support_vs_reference/README.md) scored every BFD intron genome-wide and required an identical RefSeq intron. The likelihood only sees introns that
are cells of the strict-family matrix, and a splice site shifted by a few bases is the same intron for it. The independent review of J2 asked for this
redo. `matrix_cells.py` (selftest in `bench/test_support_vs_reference.py`) takes the introns from the per-family `sites.tsv` of the J1 work folder,
joins them to `introns.tsv` for coordinates (all joined), and scores agreement with RefSeq at three levels: `exact`, `near3` (both ends within 3 nt) and
`near12` (within 12 nt). Reads are at the BFD coordinates, as in J2. Job `run.sbatch` (17 s); outputs `results/J2b_bins.tsv`, `results/J2b_summary.tsv`.

Limit: there is no RefSeq projection into the alignment columns, so `near3`/`near12` are an approximation of "same column", not the column test itself.

## Results

Cell introns, share with no reads, and agreement with RefSeq (exact / near12) by read bin:

| genome | cell introns | no reads | agree at 0 reads | agree at 1 read | agree at 2 reads | AUC exact / near12 |
|---|---|---|---|---|---|---|
| *C. neoformans* | 6,470 | 5.8% | 87.3% / 88.9% | 99.9% / 99.9% | 99.8% / 99.8% | 0.85 / 0.86 |
| *S. pombe* | 1,387 | 3.7% | 78.8% / 80.8% | 100% / 100% | 99.5% / 99.5% | 0.81 / 0.90 |
| *N. crassa* | 1,841 | 4.3% | 79.7% / 81.0% | 98.7% / 98.7% | 99.3% / 99.7% | 0.77 / 0.81 |
| *S. cerevisiae* | 46 | 15.2% | 1 of 7 | 94.6% (35/37) | 2 of 2 | 0.87 / 0.87 |
| *C. cinerea* | 6,941 | 9.3% | 82.6% / 83.6% | 91.2% / 92.0% | 91.2% / 92.3% | 0.50 / 0.50 |
| *A. nidulans* | 2,594 | 6.9% | 61.7% / 63.9% | 87.8% / 88.9% | 85.5% / 86.8% | 0.57 / 0.57 |

## What this changes

- **Few model cells lack reads.** In the matrix only 4–9% of cell introns have no reads (15% in *S. cerevisiae*, which has only 46 cells), against 20–35% genome-wide in
  J2. The zero-read prior `p0` therefore matters for few cells. J2's alarm about `p0` (0.10 to 0.81 across genomes) came from the genome-wide set, which is dominated by
  ab initio introns outside the strict families; it mostly does not apply to the cells the likelihood uses.
- **Zero-read agreement in the cells is close to the default `p0` of 0.8** for four genomes (79–87%). The outliers are *A. nidulans* (62%) and *S. cerevisiae* (1 of 7: too few cells to say).
- **The one-read step stands.** With one read, agreement is 99.9% (*C. neoformans*), 100% (*S. pombe*), 98.7% (*N. crassa*), against the default map's 0.85 at k = 3.
  In the genomes with a reliable RefSeq the default is cautious, as J2 found; in *C. cinerea* and *A. nidulans* agreement saturates at 88–92%, the weak-label case.
- **Relaxing the match barely matters**, except that it raises the *S. pombe* AUC (0.81 to 0.91) and *N. crassa* (0.77 to 0.81): shifted splice sites explain a part of those disagreements.
  Exact matching was not the main cause of disagreement.
- **Not changed:** the `junctions` defaults. The cell-level agreement is still a RefSeq agreement with reads that trained BFD (circularity), and it says nothing yet about
  whether soft evidence moves gain/loss/origin calls toward the truth (steps 3 and 4).
