---
title: 'intronfate: projecting annotated introns onto protein alignments and reconstructing their gain and loss'
tags:
  - Rust
  - comparative genomics
  - introns
  - molecular evolution
  - Dollo parsimony
  - fungi
authors:
  - name: Jason E. Stajich
    orcid: 0000-0002-7591-0020
    affiliation: 1
affiliations:
  - name: University of California-Riverside, Riverside, CA, USA   # TODO: add department if JOSS asks
    index: 1
date: TODO
bibliography: paper.bib
---

<!--
OUTLINE STATUS (2026-11-27): evidence map below is current; B2 scaling series submitted; prose not yet written.
DRAFT SKELETON. JOSS limit is roughly 750-1000 words in the body; target word counts are in brackets.
Every number below must be traced to a bench/ README or a log before submission; "TODO(bench)" marks the ones
that still need a run. Do not cite from memory: check each paper.bib entry against the publisher record.
-->

# Summary  [~120 words]

<!-- Introns sit at specific positions in homologous genes; presence/absence across species records history.
intronfate = Rust CLI + library. GFF3/GTF + protein alignments -> intron positions as (alignment column, phase)
-> presence/absence/missing matrix -> gains and losses on a species tree (Dollo, parsimony, weighted Dollo,
likelihood). Imports and exports MALIN tables. -->

# Statement of need  [~180 words]

<!-- Who: comparative genomicists studying intron density, loss and gain across large clades (fungi first).
Problem: the projection step (annotation -> alignment columns) is done with ad hoc scripts; MALIN [@Csuros2008]
starts from a hand-built matrix. Two design commitments: position identity is alignment column + phase
(docs/SPEC.md section 2); missing is not absent (a column in a gap or unannotated region stays '?').
Scale: genome-scale annotation is routine; a folder of GFFs + alignments should give a scored matrix in one step. -->

# State of the field  [~120 words]

<!-- MALIN [@Csuros2008; @Csuros2011] is the reference likelihood implementation; intronfate is complementary,
covers the upstream projection step, and reads/writes MALIN tables (validated on MALIN's published data, bench/V1, L1).
Earlier Dollo analyses of intron positions [@Rogozin2003]. State plainly what other tools do not do. -->

# Software design  [~250 words]

<!-- Figure 1: pipeline (annotation -> projection -> matrix -> scoring -> tables).
Core + CLI crates; AnnotationSource seam. Inputs: file / list-file / folder, stdin and --stdout TABLE,
--taxon-map for one alignment per family, --combine (equals the per-family sum; tested).
Scoring: dollo, mp, dollo_weighted, counts, likelihood (per-genome annotation error, family-block bootstrap,
gamma classes, free class levels). Optional junction-support soft evidence, presented as a sensitivity analysis
(bench/S1: it did not improve ranking). Testing: unit + CLI tests; MALIN validation V1; MALIN likelihood
reproduced to four decimals (L1); doc-sample tests with mutation checks so the pages cannot drift from the binary
(scripts/test_doc_samples.py, scripts/mutate_doc_samples.py). -->

# Research impact statement  [~200 words]

<!-- Three evidence blocks, each tied to a bench result. NO Antarctolithica material (private project).
Ease of use: docs/HOWTO.md run as written on the HPCC with Apptainer (bench/D1) and via conda on macOS/arm64
(within 0.2 percent of D1); give install + run commands and counts of steps.
Speed and scale: bench/V3 (1,000 genomes in 52 s and 7.0 GB; largest GFF3 in the BFD 1.6 s; single-file parsing
is single-threaded, state as a limitation). MALIN likelihood timing is PARKED (bench/B1_malin_timing README): cite L1 for likelihood equivalence, time only Dollo against MALIN, and do not claim a likelihood speedup. TODO(bench): intronfate-only scaling series (taxa x sites).
Use cases: bench/P1 intron density across the BFD; F1/F2 flank and column-ambiguity; V4/V7 annotation-error
sensitivity; J1-J3 junction support vs RefSeq. Future use: intron density and position questions at clade
scale (loss rate by lineage, convergent insertion, annotation-quality screening). -->

# AI usage disclosure

Generative AI tools were used in developing the software, its documentation and this paper.

- **Tools and models.** Claude Code (Anthropic), a coding assistant that runs in the terminal, driven by Claude models (recent Sonnet and
  Opus releases; commits written with its help carry a `Co-Authored-By: Claude` trailer in the repository history).
- **Software.** The assistant wrote and edited Rust code, unit and command-line tests, benchmark harnesses (`bench/`) and the user
  documentation (`docs/`, README). The author specified the design (position identity as alignment column plus phase, missing is not
  absent, the scoring models), directed each change, and reviewed the diffs.
- **Validation.** Correctness was not taken from the assistant. The code is tested against MALIN's published test data and outputs
  (`bench/V1`, `bench/L1`), against hand-computed fixtures, and the documentation examples are executed by tests that fail when a page
  and the binary disagree (`scripts/test_doc_samples.py`, with mutation checks). The author ran the benchmarks on the UCR HPCC and
  checked the reported numbers against the logs.
- **Paper.** The assistant drafted an outline and helped with wording and the benchmark tables; the author wrote and approved the
  final text, and checked every reference against its publisher record. No AI tool is an author, and the author takes responsibility for
  the software and the paper.

<!-- TODO(author): confirm the model names, the extent statement, and the "checked every reference" claim once it is true; compare with
JOSS's current AI policy (https://joss.readthedocs.io/en/latest/policies.html) before submission. -->

# Acknowledgements

<!-- TODO: funding, HPCC (UCR), MALIN authors. -->

# References
