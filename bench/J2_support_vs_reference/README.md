# J2: does RNA-seq junction support predict agreement with the RefSeq annotation?

Question: [J1](../J1_junction_support/README.md) turned junction reads into soft evidence with a prior map `p(n)` (`--unsupported-p`, `--half-reads`).
Does read support actually track whether an annotated intron is right? Ground truth here is the RefSeq annotation of the same genome
(the intron sets of [P1f](../P1f_adjudication) and [V4](../V4_annotation_accuracy/README.md)), which is not truth but is an independent
annotation: it was not trained on these reads.

## Method

`support_vs_reference.py` (fixture-tested: `bench/test_support_vs_reference.py`; coordinates, strand rule and AUC are asserted).
For each of the six genomes with a J1 junction BED, every intron of the BFD (funannotate rc.6) annotation is labelled `in_ref` when RefSeq has an
intron with the same sequence, coordinates and strand. Its reads are the reads at exactly those coordinates and strand (a `.` strand
counts for either). Job `run.sbatch`; outputs `results/J2_bins.tsv`, `results/J2_summary.tsv`.

## Results

Fraction of BFD introns that RefSeq also has, by junction reads at the intron (n = BFD introns in the bin):

| genome | 0 reads | 1 | 2 | 3-4 | 5-9 | 10+ |
|---|---|---|---|---|---|---|
| *C. neoformans* | 81.4% (3,623) | 97.8% | 98.7% | 99.2% | 98.8% | 98.3% |
| *S. pombe* | 67.4% (377) | 98.8% | 99.1% | 98.8% | 98.7% | 100% |
| *N. crassa* | 70.6% (3,032) | 96.6% | 98.3% | 98.5% | 98.6% | 100% |
| *S. cerevisiae* | 9.6% (115) | 95.7% | 100% | 100% | – | – |
| *C. cinerea* | 80.2% (14,068) | 87.6% | 87.7% | 86.3% | 86.0% | 87.1% |
| *A. nidulans* | 66.2% (6,844) | 83.1% | 81.7% | 82.2% | 83.3% | 77.8% |

AUC of the read count for agreement with RefSeq: *S. cerevisiae* 0.94, *S. pombe* 0.82, *N. crassa* 0.82, *C. neoformans* 0.79, *A. nidulans* 0.60,
*C. cinerea* 0.54. Fraction with at least one read, by the P1f classes (shared with RefSeq / BFD only / RefSeq only):

| genome | shared | BFD only | RefSeq only |
|---|---|---|---|
| *C. neoformans* | 91.5% | 38.2% | 39.6% |
| *S. pombe* | 94.7% | 30.1% | 43.7% |
| *N. crassa* | 83.4% | 23.2% | 15.7% |
| *S. cerevisiae* | 94.9% | 7.1% | 19.0% |
| *C. cinerea* | 78.5% | 68.5% | 7.4% |
| *A. nidulans* | 73.6% | 52.9% | 2.4% |

## What this shows, and what it does not

- **A single read is the signal; more reads add nothing.** In the three genomes where RefSeq and the annotation agree most (*C. neoformans*, *S. pombe*,
  *N. crassa*), agreement jumps from 67–81% with no reads to 97–99% with one, and stays there. These reads are digitally normalised and shallow
  (J1), so a count is not an expression level; the information is presence, not number.
- **The prior map does not match that shape.** The default `p(n) = 1 − 0.2·3/(n+3)` gives 0.80 at no reads (close to the 0.67–0.81 seen here), but only
  0.85 at one read and 0.95 at ten, where the data say about 0.97–0.99 from one read on. In these data `--half-reads` near 0.2 would fit better. This
  is a fit to six genomes, not a recommendation to change the default; the shallow normalised reads set the shape.
- **RefSeq is a weak label in two genomes.** In *A. nidulans* and *C. cinerea* the BFD-only introns have more read support (53%, 68%) than in the
  other genomes, and the RefSeq-only introns have almost none (2.4%, 7.4%): the reads side with the BFD annotation, so agreement tops out at
  83–88% and the AUC is low. That is a statement about the reference there (P1f reached the same view), not about the reads.
- **Circularity is partial.** The same reads trained the BFD annotation (J1), so the BFD introns are biased towards read support; RefSeq did not use them,
  which is why agreement with RefSeq is informative at all. A fully independent test needs reads that were not used for training.
- **Not shown:** whether soft evidence moves the likelihood's *inferences* (gains, losses, origins) towards the truth. That needs per-cell
  annotation-error calls compared with the RefSeq disagreements.
