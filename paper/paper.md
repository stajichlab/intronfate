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

<!-- PROSE DRAFT 2026-11-28 (~820 words, body only). Numbers come from: V3 (52 s, 7.0 GB, 1.6 s), B2 (10^8 cells; likelihood 1,216 s),
V6 tier A + B1 (Dollo 0.17 s vs 0.58 s, likelihood 1.5 s), L1 and V1. TODO(author) to verify before submission: the D1 conda/container
'within 0.2 percent' statement; 'rarely documented or tested' (opinion about the field); the use-case paragraph has no numbers yet
(P1, V7 results); the likelihood-scaling sentence will change if B3 finds a cheaper tolerance. -->

# Summary

Spliceosomal introns sit at specific positions in homologous genes, and their presence or absence across species records
gain and loss over evolutionary time. `intronfate` is a Rust command-line tool and library that turns genome annotations
(GFF3 or GTF) and protein alignments into intron positions, expressed as an alignment column and a phase, and builds a
presence, absence and missing matrix with one row per position and one column per genome. It then reconstructs gains and
losses on a species tree with Dollo parsimony, maximum parsimony, weighted Dollo, simple counts or a maximum-likelihood model
that allows annotation error and rate variation. It reads and writes MALIN tables [@Csuros2008], so results move between the two programs.

# Statement of need

Studies of intron density and position across a clade need the same step first: place every annotated intron into a shared
coordinate system. This is usually done with one-off scripts, and the output is a hand-built matrix that is then given to a
likelihood program such as MALIN [@Csuros2008]. The step is where errors enter, and it is rarely documented or tested.
`intronfate` makes it a single command, from a folder of annotations and alignments to a scored matrix.

Two design commitments follow from how annotation-based data fail. First, a position is identified by its alignment column
and phase, so introns from different genomes are the same position when they fall in the same column and phase. Second,
missing is not absent: a column that is in an alignment gap, or in a region a genome does not annotate, stays `?` and is not
counted as a loss. Treating gaps as losses inflates loss estimates in exactly the places where annotations are weakest.
The target users are comparative genomicists working across tens to thousands of genomes, with fungi as the first clade.

# State of the field

MALIN [@Csuros2008; @Csuros2011] is the reference program for likelihood analysis of intron presence and absence, with
rate variation across positions and branches. It starts from a table of intron presence and absence, and its table is
prepared outside the program. Earlier analyses of intron positions used Dollo-style counts [@Rogozin2003]. `intronfate` is
complementary: it produces the table from annotation, then offers its own scoring and exports the table to MALIN. It is
a reimplementation of the scoring models, not a replacement for MALIN's graphical interface or its published methods.

# Software design

The core crate holds the annotation reader, the projection onto alignments, the matrix and the scoring models, and a
second crate provides the command-line interface. Annotation enters through an `AnnotationSource` seam, so other formats can
be added. Input can be one file, a list file, a folder, or standard input, and `--taxon-map` handles one alignment per
gene family. `--combine` joins families into one matrix and equals the sum of the per-family results, which a test checks.

The likelihood model has per-genome annotation error, a family-block bootstrap, gamma or free-level rate classes, and
MALIN's correction for positions never seen with an intron. Documentation examples are executed as tests, with mutation
checks that confirm a changed binary breaks them, so the pages cannot drift from the program.

Validation against MALIN used its published test data. Branch-wise losses from Dollo parsimony on the Rogozin table
match the MALIN guide exactly. For the likelihood, the published Rogozin parameters score within 0.002 log-likelihood units
of this tool's own maximum, and the Chromo23 log-likelihood, with its never-seen correction, reproduces MALIN's printed
value (-28430.9249) to four decimals.

# Research impact statement

**Ease of use.** The documented workflow, from public genomes to a scored matrix, was run as written on a cluster with a
container and again with conda on macOS. The two runs agreed within 0.2 percent. A MALIN table converts to and from the
matrix with one command each (`convert --from-malin`, `--to-malin`).

**Speed and scale.** Projection of 1,000 fungal genomes in one process takes 52 s and 7.0 GB; the largest
annotation in the genome collection tested (165 MB uncompressed) takes 1.6 s. On simulated data, import, Dollo and parsimony
scale linearly in matrix size: a matrix of 1,000 taxa and 100,000 sites (10^8 cells) imports in 26 s and is scored in about 35 s,
using 5.8 GB. The likelihood does not scale like this: on the same simulated data it took 1,216 s for 30 taxa and 100,000
sites, so it is suited to tens of taxa, not hundreds. On a real nine-species set with 15,697 sites it takes 1.5 s. For Dollo
parsimony on that set, `intronfate` took 0.17 s and MALIN's command-line program 0.58 s; both are fast, and we make no
speed claim for the likelihood against MALIN.

**Use cases.** The tool supports questions about intron density and position at clade scale: loss rates by lineage,
convergent insertion at the same position, and screening of annotation quality, with error-sensitivity analyses to show
how much a result depends on annotation misses.

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
