# intronfate

Project annotated intron positions into protein sequences and alignment columns, then score whether those positions are
**shared** across orthologs: Dollo parsimony (`dollo`), Fitch parsimony (`mp`), weighted Dollo (`dollo_weighted`), plain counts
(`counts`), or a likelihood model with per-genome annotation error (`likelihood`).

Rust workspace, no runtime dependencies beyond the binary.

## Start here

1. **See it work (5 minutes):** run the demo ([Try it](#try-it)) and the hand-checked [`examples/`](examples/README.md).
2. **A full run on public genomes:** [`docs/HOWTO.md`](docs/HOWTO.md), copy-paste steps from NCBI download to gains and losses per
   branch, for five or nine fungi. Needs OrthoFinder and MAFFT; run it on a cluster node.
3. **Your own genomes:** [`docs/TUTORIAL.md`](docs/TUTORIAL.md) explains each step of the pipeline and what to check.
4. **Reference:** [`docs/USAGE.md`](docs/USAGE.md) (commands, options, methods, diagnostics) and [`docs/OUTPUTS.md`](docs/OUTPUTS.md)
   (every output file and column).

Also: [`docs/SPEC.md`](docs/SPEC.md) design, assumptions and references · [`docs/TESTING.md`](docs/TESTING.md) the test suite ·
[`bench/README.md`](bench/README.md) validation runs on real data · [`docs/HANDOFF.md`](docs/HANDOFF.md) current state of the work.

## What it does

```text
GFF3/GTF ──► intron ──► cds_offset ──► (codon, phase) ──► MSA column ──► site × taxon matrix ──► scores
             genome      coding bases    0/1/2 nt         homology        1 / 0 / ? per taxon     dollo, mp, dollo_weighted,
             space       into the CDS    into the codon   unit            per position            counts, likelihood
```

The key idea: genome coordinates are useless across species, and residue numbers drift with indels —
**alignment columns** are the homology unit. An intron's position is its CDS nucleotide offset, split
into codon + phase, looked up in the protein alignment.

## Install

Requires Rust 1.85 or newer ([rustup.rs](https://rustup.rs), or `conda install -c conda-forge rust`).

```sh
git clone https://github.com/stajichlab/intronfate.git
cd intronfate
cargo build --release          # binary at target/release/intronfate
cargo test --workspace         # full suite, see docs/TESTING.md
```

The binary has no runtime dependencies; copy it or put `target/release` on your `PATH`.

## Try it

```sh
cargo run --release -p intronfate-cli -- demo -o intron_demo
```

Writes a worked five-species example (inputs + every table) and prints the interesting ones:

```text
site    column  phase  n_taxa  n_observable  pattern  species1  species2  species3  species4  outgroup
4:0     4       0      4       5            11101    1         1         1         0         1
6:0     6       0      2       5            01100    0         1         1         0         0
9:0     9       0      3       5            10011    1         0         0         1         1

site    column  phase  n_present  n_missing  gains  losses  monophyletic  dollo_compatible  gained_at            lost_in
4:0     4       0      4          0          1      1       no            no                root                 species4
6:0     6       0      2          0          1      0       yes           yes               species2,species3
9:0     9       0      3          0          1      1       no            no                root                 species2,species3

site    column  phase  n_present  n_missing  steps  informative  ambiguous
4:0     4       0      4          0          1      no           yes
6:0     6       0      2          0          1      yes          yes
9:0     9       0      3          0          1      yes          yes
```

Three intron positions: two ancestral (each lost once) and one derived insertion on the
species2/species3 stem. The annotated tree shows where the origins were placed:

```text
(outgroup,((species1,(species2,species3)[gains=1])A,species4))[gains=2];
```

## Real data

One pass over many gene-family alignments, with every family also stacked into one matrix (`out/combined/`):

```sh
intronfate analyze --gff-dir genomes/ --proteins-dir genomes/ --proteins-ext faa \
                   --msa-dir alignments/ --msa-ext aln.fasta --taxon-map taxon_map.tsv \
                   -t species.nwk --outgroup Spunctatus --combine -o out
```

`matrix.tsv` is the interchange format: `intronfate score --matrix out/combined/matrix.tsv ...` re-scores it under another tree,
rooting or method in seconds, without re-reading any annotation. The separate `project`, `matrix` and `score` commands, every option
and the input forms (files, list files, folders) are in [`docs/USAGE.md`](docs/USAGE.md); the whole pipeline from download to results
is [`docs/HOWTO.md`](docs/HOWTO.md).

## Output

Read these first: `diagnostics.tsv` and `unmatched.tsv` (what failed to join, and why), then `combined/dollo_branches.tsv` (gains and
losses per branch) or, per family, `<family>/summary_score.tsv` and `dollo_sites.tsv`. Every file and column is described in
[`docs/OUTPUTS.md`](docs/OUTPUTS.md).

## Things that will bite you if you don't read them

- **Phase is part of site identity.** `codon 5 phase 0` and `codon 5 phase 1` are different sites.
  Use `--ignore-phase` deliberately, not by accident.
- **Missing is `?`, never `0`.** An unannotated taxon, a gap around the site, or an alignment row that is not the annotated protein
  is *no observation*; scoring it as "intron lost" would invent losses. `--flank` sets how much alignment around a site a taxon needs.
  Details: [`docs/USAGE.md`](docs/USAGE.md) §6.
- **Classic Dollo minimises losses, not gains.** One gain per site is assumed, so the informative statistic is the loss count;
  `-m dollo_weighted` sweeps the gain:loss cost when you need a real trade-off ([`docs/SPEC.md`](docs/SPEC.md) §5.1–5.2).
- **Dollo gain placement depends on rooting.** Pass `--outgroup`; if the Newick looks unrooted the tool warns instead of guessing.
- **`introns.tsv` numbers are per-protein.** Cross-species comparisons use `matrix.tsv` columns.
- **Annotation and sequence must correspond.** `unmatched.tsv` and the `length_mismatch` diagnostic exist because a join onto the
  wrong paralog otherwise yields plausible-looking tables.
- **The result is only as good as the annotated introns.** A missed intron looks like a loss and a false one like a gain;
  [`bench/V7`](bench/V7_annotation_error_sensitivity/README.md) measures how far that moves the counts, and `-m likelihood`
  with `--error-rates` models it ([`docs/USAGE.md`](docs/USAGE.md) §11).

## Validation

Every step was run on a university HPC cluster; scripts, small result tables and a README per step are in [`bench/`](bench/), indexed
with one-line results in [`bench/README.md`](bench/README.md): agreement with MALIN's published Dollo results, projection on nine
fungal genomes, genome-scale timing, annotation accuracy and its effect on the counts, ortholog tiers, and the how-to run end to end.
These are validation and first-use results, not a finished biological study. The tools have also been applied to a separate
research project.

## Layout

```text
crates/intronfate-core/  library: gff, annotation, seqio, msa, project, matrix, tree, score, likelihood, malin, tables, inputs, ...
crates/intronfate-cli/   the `intronfate` binary
docs/                    HOWTO, TUTORIAL, USAGE, OUTPUTS, SPEC (design + references), TESTING, HANDOFF (state)
examples/                three runnable examples with committed expected output, plus 04_public_genomes (inputs for docs/HOWTO.md)
bench/                   validation runs on the HPCC (index: bench/README.md)
tests/data/              tiny synthetic fixtures
scripts/                 helpers for docs/HOWTO.md (one_protein_per_gene.py, make_families.py, align_families.sh, doc_blocks.py);
                         malin_to_matrix.py, an independent MALIN importer used as a cross-check
```

## Status

Implemented: GFF3/GTF ingestion with frame/exon/length cross-checks, protein-space projection, FASTA/Stockholm alignments,
site × taxon matrices with per-cell missing data, pairwise sharing, the five scoring methods above (`likelihood` with
`--error-rates`, bootstrap intervals and gamma rate classes), TSV/NEXUS output, MALIN import/export (`convert`), and combined
matrices over many families (`--combine`).

Not yet: SQL/Chado/Ensembl annotation sources (the `AnnotationSource` seam is there), separate gain and loss rate factors in the
likelihood model, validation against MALIN's `*.rates`, stdin/stdout support. The full list: [`docs/HANDOFF.md`](docs/HANDOFF.md) and
[`docs/SPEC.md`](docs/SPEC.md) §12.
