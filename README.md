# intronfate

*Formerly `intron-projector`.* Working folders on the laptop and the HPCC are still named `intron_projector`; only the repository, crates, binary and environment variable (`INTRONFATE`) changed name.

Project annotated intron positions into protein sequences and alignment columns, then score whether
those positions are **shared** across orthologs — under Dollo parsimony, maximum parsimony, or plain
counts.

Rust workspace, no runtime dependencies beyond the binary.

**Documentation:** [`docs/USAGE.md`](docs/USAGE.md) user guide ·
[`docs/OUTPUTS.md`](docs/OUTPUTS.md) every output column · [`examples/`](examples/) runnable, hand-checked
examples · [`docs/SPEC.md`](docs/SPEC.md) design, assumptions and references ·
[`docs/TESTING.md`](docs/TESTING.md) the test suite · [`docs/HANDOFF.md`](docs/HANDOFF.md) current state.

## What it does

```
GFF3/GTF ──► intron ──► cds_offset ──► (codon, phase) ──► MSA column ──► site × taxon matrix ──► scores
             genome      coding bases    0/1/2 nt         homology        0/1 per taxon        Dollo / MP /
             space       into the CDS    into the codon   unit            per position         counts
```

The key idea: genome coordinates are useless across species, and residue numbers drift with indels —
**alignment columns** are the homology unit. An intron's position is its CDS nucleotide offset, split
into codon + phase, looked up in the protein alignment.

## Install / build

```sh
cargo build --release          # binary at target/release/intronfate
cargo test --workspace         # full suite, see docs/TESTING.md
```

## Try it

```sh
cargo run -p intronfate-cli -- demo -o intron_demo
```

Writes a worked five-species example (inputs + every table) and prints the interesting ones:

```
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

```
(outgroup,((species1,(species2,species3)[gains=1])A,species4))[gains=2];
```

## Real data

```sh
# 1. intron positions in protein space (no alignment needed yet)
intronfate project --gff species.gff.gz -p species.faa -o out

# 2. project onto an ortholog alignment and build the matrices
intronfate matrix --gff a.gff.gz --gff b.gff.gz \
                        -p a.faa -p b.faa \
                        --msa family.fasta -o out

# 3. score sharing on a tree
intronfate score --matrix out/family/matrix.tsv \
                       --tree family.nwk --outgroup Species_x \
                       --methods dollo,mp -o out/family
```

Or one pass with `analyze`. `matrix.tsv` is the interchange format, so step 3 is cheap and can be
re-run over rootings, tolerances and methods without re-parsing annotation.

## Output

Per alignment, in `out/<alignment>/`: `matrix.tsv` (binary site × taxon), `sites.tsv` (long form with
the full coordinate stack), `taxa.tsv` (what was projected and what was dropped), `pairwise.tsv` +
`sharing_matrix.tsv` (shared counts, Jaccard), `spectrum.tsv` (sites shared by *k* taxa),
`counts_sites.tsv`, `characters.nex` (NEXUS for external phylogenetics). With `score`/`analyze`:
`dollo_sites.tsv`, `mp_sites.tsv`, `reconstructions.tsv`, `tree_gains.nwk`, `summary_score.tsv`.
Globally: `introns.tsv`, `unmatched.tsv`, `diagnostics.tsv` (+ `_detail`).

`demo` writes a focused subset of these so the interesting tables print next to the inputs it
generates (`demo.gff3`, `demo.faa`, `demo.aln.fasta`, `demo.nwk`).

## Things that will bite you if you don't read them

- **Phase is part of site identity.** `codon 5 phase 0` and `codon 5 phase 1` are different sites.
  Use `--ignore-phase` deliberately, not by accident.
- **Missing is `?`, never `0`.** An unannotated taxon, a gap around the site, or an alignment row that
  is not the annotated protein is *no observation*; treating it as "intron lost" invents losses.
  `--flank` controls how much alignment around a site a taxon needs. See `docs/USAGE.md` §6.
- **Classic Dollo minimises losses, not gains.** With gain cost 1 and loss cost 0 and free losses,
  one gain at the root explains any pattern — the min-gain score is degenerate at exactly 1 for every
  site. The informative statistic is the loss count given one origin; `dollo_weighted` sweeps the
  gain/loss ratio when you need a real trade-off.
- **Dollo gain placement depends on rooting.** Pass `--outgroup`; if the Newick looks unrooted the
  tool warns instead of guessing.
- **`introns.tsv` numbers are per-protein.** Cross-species comparisons use `matrix.tsv` columns.
- **Annotation/sequence must correspond.** `unmatched.tsv` and the `length_mismatch` diagnostic exist
  because a join onto the wrong paralog otherwise yields plausible-looking tables.

## Layout

```
crates/intronfate-core/  library: gff, annotation, seqio, msa, project, matrix, tree, score, tables
crates/intronfate-cli/   the `intronfate` binary
docs/                 USAGE, OUTPUTS, TESTING, SPEC (design + references), HANDOFF (state)
examples/             three runnable examples with committed expected output
tests/data/           tiny synthetic fixtures
scripts/              malin_to_matrix.py — independent MALIN importer used as a cross-check
```

## Status

Implemented: GFF3/GTF ingestion with frame/exon/length cross-checks, protein-space projection,
FASTA/Stockholm alignments, site × taxon matrices, pairwise sharing, Dollo (classic + weighted),
Fitch MP with ancestral reconstruction, counts, TSV/NEXUS output.

Also: per-cell missing data (alignment gaps are `?`), MALIN import/export (`convert`), `--min-observed`.

Not yet: SQL/Chado/Ensembl sources (the `AnnotationSource` seam is there), likelihood-based
gain/loss models, benchmarks on real genomes (run those on the HPCC), stdin support.
