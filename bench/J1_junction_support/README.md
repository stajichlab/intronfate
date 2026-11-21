# J1: RNA-seq junction support as soft evidence (V6 tier A families)

Question: how many annotated introns do RNA-seq junction reads support, and do the likelihood and the convergent-insertion calls
(`--ml-origins`) move when each supported cell is given a probability instead of a hard 1 (`intronfate junctions`, SPEC §5.6b)?

## Inputs and method

- **Families:** the 1,106 strict families of [V6](../V6_orthology_nine_species/README.md) (`alignments/A_strict`), nine species, projected with
  `analyze --combine` into one combined matrix (15,697 sites).
- **Junctions:** from the hisat2 alignments that funannotate training left in BFD
  (`genome_annotation/<genome>/training/funannotate_train.coordSorted.bam`), turned into BED by `bam_junctions.py` (fixture-tested, with the
  coordinate convention pinned end to end by `bench/test_junction_scripts.py`). Six species have a BAM: *A. nidulans*, *C. cinerea*,
  *C. neoformans*, *N. crassa*, *S. cerevisiae*, *S. pombe*. *S. punctatus* has empty read files, *M. mucedo* and *M. maydis* have reads but
  no BAM. Cells of those three, and any cell the matrix masks as unobservable, stay hard.
- **Fits:** `score` of the combined matrix, `-m likelihood --ml-origins`, hard (no soft table) and soft for `(p0, k)` in
  (0.8, 3), (0.5, 3), (0.95, 3), (0.8, 1), (0.8, 10); once with 4 gamma rate classes (`results/sensitivity.tsv`) and once with constant rates
  (`TAG=rc1_ RATECAT=1`, `results/sensitivity.rc1_.tsv`). `run.sbatch` runs everything; `summarize.py` writes the tables. `work/` is a purgeable cache.

## What the data are, and are not

- The reads are **digitally normalised** (`*_norm_*` in BFD) and shallow (a BAM of 4–14 MB per species): read counts are not expression levels,
  and the depth differs by species (2% of *A. nidulans* cells have three or more reads, 47% of *C. neoformans*), so `q` here varies with
  sequencing depth as much as with annotation quality.
- The same reads trained the funannotate annotation, so **support is not independent of the annotation**: high support partly confirms the
  annotation against its own training data.
- `p0` and `k` are priors, not measurements.

## Results

Support of annotated introns (cells of the combined matrix, junction reads at exactly the intron's coordinates and strand):

| taxon | cells | ≥1 read | ≥3 reads | ≥10 reads |
|---|---|---|---|---|
| A. nidulans | 2,442 | 93.2% | 2.1% | 0.04% |
| C. cinerea | 6,502 | 90.8% | 31.5% | 1.1% |
| C. neoformans | 6,055 | 94.2% | 47.2% | 0.5% |
| N. crassa | 1,728 | 95.8% | 19.0% | 0% |
| S. cerevisiae | 40 | 82.5% | 0% | 0% |
| S. pombe | 1,279 | 96.5% | 35.3% | 0.6% |

Likelihood and convergent-insertion calls, constant rates (converged; `results/sensitivity.rc1_.tsv`):

| setting | ln L | sites with p(≥2 origins) > 0.5 | Spearman of ln LR with hard | hard top-100 kept |
|---|---|---|---|---|
| hard | −46,710 | 125 | 1 | 100 |
| p0 0.95, k 3 | −45,707 | 0 | 0.98 | 10 |
| p0 0.8, k 1 | −44,673 | 0 | 0.90 | 5 |
| p0 0.8, k 3 | −43,657 | 0 | 0.85 | 3 |
| p0 0.8, k 10 | −42,772 | 0 | 0.83 | 2 |
| p0 0.5, k 3 | −39,680 | 0 | −0.12 | 0 |

What can and cannot be said:

- The count of sites above an ln LR threshold is **not a usable comparison** across fits: ln LR is measured against each fit's own prior odds,
  and those differ (it ranged from 4 to 12,841 here). Use rank agreement and the posterior instead.
- Soft ln L is not comparable with hard ln L (different data); only differences among soft settings on the same data mean anything.
- With cells 85–95% reliable, **no site keeps a posterior above 0.5 for two or more origins**, and the hard-model top candidates are mostly not
  preserved. The model can now explain a scattered presence as an annotation error instead of an independent insertion, so the convergent-insertion
  signal depends on how far the annotation is trusted. In this data set it is not robust to that. This does not show which model is right; it
  shows that the hard-model calls rest entirely on the annotation being correct.
- With four rate classes the fits are slow to converge and reach the iteration cap (the log-likelihood is bounded and climbs 0.006 per iteration; the hard fit behaves the same way); their rankings are less stable (Spearman 0.03–0.80). [J3](../J3_shape_ridge/README.md) shows the gamma shape itself is well determined (the profile is peaked, 95% interval about 3–6 for the soft fit) and that what is slow is the EM, with about 2 log-likelihood units unfinished. Prefer the constant-rate table for rankings, or the converged free-shape fit.
- `p0 = 0.5` is pathological (negative rank correlation, 12,841 flagged sites); disregard it.

## Bugs this run found

1. `junctions` softened cells the matrix had masked as unobservable (and `score` then failed with "no site"); it now requires `--matrix`.
2. The soft-evidence likelihood did not condition on the position having been seen: ln L went positive (+31,884) and the gain rate collapsed.
   Fixed (SPEC §5.6b) with a brute-force check; results from before the fix are invalid.

## Not done

- Whether read support predicts annotation errors: [J2](../J2_support_vs_reference/README.md) compares reads with the RefSeq annotation (one read is the signal; RefSeq is a weak label in two genomes). Whether the likelihood's calls move towards the truth, and non-normalised RNA-seq, are open.
- Junctions at coordinates the annotation lacks are not projected into the alignment.
- *S. punctatus*, *M. mucedo* and *M. maydis* are not covered.
