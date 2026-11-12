# V4 — how far can annotation error move the intron gain/loss results?

**Scope.** This folder is about one question for *this tool*: the shared/lost-intron results rest on each genome's annotated introns, so how
much of a reported loss or gain could be annotation error? It uses another project's accuracy scores as input. The broader
comparison of funannotate gene-model accuracy (for the funannotate manuscript) is a separate project with its own plan, in
`funannotate-manuscript/assessment_reference_comparison/`.

The V2/V3 analyses use the BFD funannotate rc.6 predictions (`genome_annotation/<genome>/predict_results/`). There has been
critique of the released gene models. This folder collects the accuracy evidence and turns it into bounds on what the tool
reports. **No numbers here were measured by this project:** the scores come from a colleague's session
(`monitor-gene-stats-pipeline`, gffcompare, funannotate rc.6 vs NCBI RefSeq GFF3, files in
`BFD/Fungi_BFD_runs/analysis/rc6_vs_pre_rc6/reference_compare.tsv`), received 2026-11-11.

## What is known

CDS mode (CDS segments scored as exons), rc.6 predictions, intron level (gffcompare "intron level"; Pr = precision):

| genome | in V2? | intron Sn % | intron Pr % | exon Sn/Pr % | intron-chain Sn/Pr % | locus Sn/Pr % |
|---|---|---|---|---|---|---|
| *Cryptococcus neoformans* H99 | yes | **94.6** | **97.0** | 92.2 / 94.8 | 78.3 / 85.3 | 83.1 / 83.8 |
| *Neurospora crassa* OR74A | yes | **81.1** | **91.7** | 77.5 / 86.2 | 70.4 / 82.1 | 76.0 / 79.2 |
| *Botrytis cinerea* B05.10 | no | 92.4 | 93.7 | 88.4 / 90.1 | 78.8 / 86.8 | 87.8 / 83.7 |
| *Zymoseptoria tritici* IPO323 | no | 64.8 | 68.7 | 58.2 / 55.4 | 48.6 / 46.2 | 53.5 / 43.8 |

The other V2 genomes were scored afterwards (same colleague, 2026-11-11); see the tables below. All scores are the colleague's gffcompare results, not measurements by this project.

Cross-check: gffcompare's predicted-intron counts for H99 (35 967) and OR74A (14 078) equal the intron counts this tool read from
the same GFF3s in V2, so both see the same intron set.

The earlier (pre-rc.6, funannotate 1.9.0) models were far worse (OR74A intron Sn 2.8 %), so which model set an analysis used
matters; V2 and V3 used rc.6.

Caveats from the scorer, which apply here too: RefSeq is not truth; RefSeq has alternative isoforms while predictions have one
mRNA per gene (caps transcript-level Sn; intron and locus levels are the better figures); the funannotate runs used
`swissprot_fungi.faa` as protein evidence, which may contain proteins of these very species, so scores may be inflated; four
strains only.

## What it implies for the V2 results — measured over the genes the analysis uses

`busco_gene_accuracy.py` (job 32031140; outputs in `results/`) restricts the colleague's per-intron table to the V2 transcripts:
their predicted introns for Pr, and the introns of the RefSeq transcripts they map to in the gffcompare `.tmap` (class codes
`= c j k m n e o`, i.e. the same locus) for Sn. Exact coordinates must match. The restriction matters a lot, because conserved
single-copy (BUSCO) genes are modelled far better than the genome average:

| genome | intron Sn: genome-wide → V2 genes (n ref introns) | intron Pr: genome-wide → V2 genes (n pred introns) |
|---|---|---|
| *S. punctatus* | 88.7 → **94.0** (4 846) | 92.1 → 95.1 (4 934) |
| *M. mucedo* | 79.9 → **95.3** (2 699) | 72.5 → 89.3 (2 975) |
| *M. maydis* | 72.0 → **88.8** (223) | 76.0 → 82.2 (253) |
| *C. cinerea* | 84.9 → **92.8** (4 794) | 85.5 → 90.2 (4 922) |
| *C. neoformans* | 94.6 → **99.0** (4 476) | 97.0 → 98.9 (4 506) |
| *S. pombe* | 94.2 → **97.3** (903) | 96.5 → 99.0 (890) |
| *N. crassa* | 81.1 → **98.7** (1 323) | 91.7 → 99.0 (1 320) |
| *A. nidulans* | 70.2 → **89.3** (1 837) | 77.7 → 88.0 (1 858) |
| *S. cerevisiae* | not scored: the V2 assembly (AMP_1a) has no RefSeq annotation; the S288C proxy has only 312 reference introns | |

Errors are close to uniform across intron phase and position in the gene (`busco_accuracy_by_class.tsv`; for example *N. crassa* Sn
98.2–99.3 % across all phases and positions), so annotation error does not look like it systematically favours one phase or the ends
of genes — which matters because sites are identified by column *and phase*.

**Terminal-branch losses vs. what annotation misses could explain** (`loss_ceiling.tsv`; ceiling = `sites_carried × (1 − Sn)/Sn`,
assuming every missed intron sits at a site other species share, so a true upper bound for this mechanism):

| genome | terminal Dollo losses | missed-intron ceiling (V2-gene Sn) | share | with the earlier genome-wide Sn |
|---|---|---|---|---|
| *N. crassa* | 470 | 16 | **3.5 %** | (61 %) |
| *C. neoformans* | 803 | 42 | **5.3 %** | (30 %) |
| *S. pombe* | 720 | 23 | **3.2 %** | (7 %) |
| *M. maydis* | 2 883 | 29 | **1.0 %** | (3 %) |
| *M. mucedo* | 963 | 137 | **14 %** | (72 %) |
| *C. cinerea* | 351 | 358 | ≥ 100 % | (234 %) |
| *A. nidulans* | 98 | 210 | ≥ 100 % | (759 %) |
| *S. punctatus* | 0 (outgroup) | – | – | – |

So, for the V2 analysis:
- Terminal-branch loss counts for *N. crassa*, *C. neoformans*, *S. pombe*, *M. maydis* (and, less firmly, *M. mucedo*) are **not**
  explained by annotation misses: the ceiling is a few percent. The "up to ~60 % for *N. crassa*" figure given earlier (from genome-wide
  Sn) was wrong for this analysis and is superseded.
- For *C. cinerea* and *A. nidulans* the ceiling is as large as the observed losses, so their terminal losses **cannot be told apart from
  annotation misses**; treat them as unreliable. The same lower precision (Pr ~ 88–90 %) means a tenth of their introns (and of *M. mucedo*'s)
  are not in RefSeq, so their terminal **gains** carry a similar caveat.
- Internal-branch results rest on several species agreeing and are less exposed still.

Limits of this measurement: it is a comparison with RefSeq, which is not truth (a "novel" intron may be real, a "missed" one may be an
annotation error in RefSeq); the funannotate runs used protein evidence that may include these species' own proteins, which can
inflate scores; exact-coordinate matching; the V2 genes are conserved BUSCO genes, which is exactly the population the analysis uses but
not representative of the genome as a whole; *S. maydis* has few reference introns at these genes (223), so its rates are noisy.

## Still to do

1. **Sensitivity run:** perturb the V2 `matrix.tsv` files at each genome's measured `1 − Sn` / `1 − Pr` and report the spread of per-branch
   losses and gains (the ceiling above is the worst case for one mechanism; this gives the typical effect).
2. **Direct comparison on genomes with both annotations:** run the pipeline on the RefSeq GFF3s (paths in the scoring folder's
   `reference_pairs.tsv`) and compare projected `(column, phase)` sites with the funannotate ones.
3. *S. cerevisiae*: the V2 genome (BFD folder `Saccharomyces_cerevisiae`, assembly GCA_029989085.1, strain AMP_1a) has no RefSeq annotation, so
   nothing can score it. The colleague has already scored the reference strain, BFD folder `Saccharomyces_cerevisiae_S288C`
   (GCF_000146045.2_R64): rc.6, CDS mode, intron Sn/Pr 73.3 / 66.0, intron-chain 73.4 / 66.3, locus 89.1 / 96.5, but on only 312 reference
   (329 predicted) introns, so any per-gene, per-phase or BUSCO-restricted breakdown would have very few events (not counted yet). The
   two are different strains and assemblies, so S288C's scores do not describe AMP_1a. Decision for the project owner: swap S288C into the
   V2 set so every genome has a measured accuracy, or keep AMP_1a and state that it is unscored. Files (read-only):
   predictions `BFD/Fungi_BFD_runs/genome_annotation/Saccharomyces_cerevisiae_S288C/predict_results/*.gff3.gz`; RefSeq
   `1KFG/2026/NCBI_fungi/source/NCBI_ASM/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gff.gz`; per-intron table and gffcompare maps under
   `analysis/rc6_vs_pre_rc6/reference_introns/` and `reference_gffcompare/` (`Saccharomyces_cerevisiae_S288C.new.*`).


## Rescoring over the V6 tiers
`run_v6.sbatch` runs `busco_gene_accuracy.py` (now with `--genomes/--taxon-map/--projection/--out`) over each V6 tier's genes into `results_v6/<tier>/`. *S. cerevisiae* S288C is scored there (tier A genes: Sn 97.4 %, Pr 82.6 %, but over only 39 reference and 46 predicted introns, so noisy; 343 reference introns genome-wide), so the earlier "unscored" note no longer applies to the V6 set. Tier C has no *S. punctatus* genes (the N1 clade excludes it), so that row has no reference introns and the script now skips such rows instead of crashing.
