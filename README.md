# intronfate

*Formerly `intron-projector`.* Working folders on the laptop and the HPCC are still named `intron_projector`; only the repository, crates, binary and environment variable (`INTRONFATE`) changed name.

Project annotated intron positions into protein sequences and alignment columns, then score whether
those positions are **shared** across orthologs — under Dollo parsimony, maximum parsimony, or plain
counts.

Rust workspace, no runtime dependencies beyond the binary.

**Documentation:** [`docs/HOWTO.md`](docs/HOWTO.md) start here: copy-paste steps on complete public genomes (five or nine fungi, downloaded from NCBI) ·
[`docs/TUTORIAL.md`](docs/TUTORIAL.md) the same pipeline explained, for your own genomes ·
[`docs/USAGE.md`](docs/USAGE.md) user guide ·
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

For a complete run on public genomes, from download to results, follow [`docs/HOWTO.md`](docs/HOWTO.md). The commands in short:

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

Or one pass with `analyze`. For many families or genomes, give lists or folders instead of repeating flags
(`--msa-dir alignments/ --msa-ext aln.fasta`, `--msa-list families.txt`, likewise `--gff-dir/--gff-list` and `--proteins-dir/--proteins-list`), and add
`--combine` to also stack every family into one matrix, scored once (`out/combined/`; see [`docs/USAGE.md`](docs/USAGE.md) §9):

```sh
intronfate analyze --gff-dir genomes/ --proteins-dir genomes/ --proteins-ext faa \
                   --msa-dir alignments/ --msa-ext aln.fasta --taxon-map taxon_map.tsv \
                   -t species.nwk --combine -o out        # out/<family>/... and out/combined/dollo_branches.tsv
```

`matrix.tsv` is the interchange format, so step 3 is cheap and can be
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

## Results so far

Everything below was run on a university HPC cluster; the scripts, the small result tables and a README for each step are in [`bench/`](bench/).
These are validation and first-use results, not a finished biological study.

| Step | What was tested | Result |
|---|---|---|
| [V1](bench/V1_malin_validation/README.md) | Scoring and MALIN interop on MALIN's published test data | Every figure in the spec reproduced. Per-branch Dollo losses match the MALIN guide exactly (Rogozin table: 7,221 sites, 2,770 losses); MALIN tables round-trip losslessly. |
| [V2](bench/V2_gff_projection/README.md) | GFF3 to protein space to alignment columns on nine fungal genomes | 739 BUSCO families (single-copy in at least 7 of 9 genomes), 6,363 protein rows, every row joined to its transcript; all families through `analyze` in 8.4 s and 507 MB. |
| [V3](bench/V3_gff_scale/README.md) | Genome-scale `project` | 1,000 genomes in one process: 294 s and 55.8 GB at the start, **52 s and 7.0 GB** after the optimisations (the largest GFF3 in the BFD, 165 MB uncompressed, takes 1.6 s). A single file does not scale with threads. |
| [V4](bench/V4_annotation_accuracy/README.md) | How accurate are the annotated introns? (funannotate against RefSeq, at the genes actually used) | Intron sensitivity at those genes is 89-99 % for eight of the nine species (*S. cerevisiae* is scored separately, over very few introns); per-species figures and the resulting ceilings on how many apparent losses could be annotation misses are in the README. |
| [V6](bench/V6_orthology_nine_species/README.md) | Ortholog families that are as single-copy as possible (OrthoFinder, synteny refinement, tiers) | 8,766 orthogroups; 1,106 single-copy in all 9 species, 689 more in 7-8, 401 more from splitting multi-copy orthogroups. Synteny refinement added none (166 groups of about two genes). |
| [V7](bench/V7_annotation_error_sensitivity/README.md) | How much could annotation error change the Dollo losses and gains? | Terminal losses in *S. maydis*, *S. cerevisiae*, *C. neoformans*, *S. pombe* and *N. crassa* move by at most 6 % even when every missed intron is assumed to hide a loss; *M. mucedo* by about 15 %; for *C. cinerea* and *A. nidulans* the ceiling removes all of their terminal losses, so those counts cannot be told from annotation misses. |
| [V8](bench/V8_antarctolithica_placement/README.md), [V9](bench/V9_antarctolithica_expanded/README.md) | Do intron positions support the placement of an unplaced Antarctic ascomycete (31, then 42 genomes)? | Whole-matrix parsimony and Dollo are miscalibrated for an intron-poor lineage; a calibrated clade-specific retention test places it outside Pezizomycotina, Saccharomycotina and Taphrinomycetes, without resolving its deeper position. |
| [V10](bench/V10_antarctolithica_domains/README.md) | Pfam domain composition (472 genomes): PCoA, random forest, signature test | A normal fungal repertoire; between the yeasts and the Pezizomycotina, nearest to early-diverging yeast lineages. |
| [V11](bench/V11_antarctolithica_topology_tests/README.md) | Likelihood tests of nine candidate placements (IQ-TREE, RAxML-NG cross-check; C60 PMSF, site removal, recoding) | Same best tree in every variant; the gap to the alternative published placement narrows with a better model but does not close. |
| [V12](bench/V12_antarctolithica_composition/README.md) | Genome traits, amino-acid / codon / Pfam composition, isolation test, random forest | Compact yeast-like genome; composition does not mark a separate group, only the Pfam repertoire is isolated (as reduced-genome lineages are). |

Summary of V4-V12 and what they do and do not show: [`docs/ANTARCTOLITHICA_PLACEMENT.md`](docs/ANTARCTOLITHICA_PLACEMENT.md).

**A worked result.** On the 1,106 single-copy families (nine fungal genomes; 15,697 inferred intron gains and 14,316 losses in total) the largest
terminal-branch loss counts are in *Mycosarcoma maydis* (4,055), *M. mucedo* (1,516), *S. cerevisiae* (1,447) and *C. neoformans* (1,091). For
*M. maydis* and *S. cerevisiae*, whose low intron density is well known, that is the expected direction, and V7 shows these counts are not
artefacts of missed introns. The per-lineage totals come from one command (`analyze --combine`, see below) with nothing to sum afterwards.
The other counts need the limits below before they are interpreted.

**What the results do and do not show.**

- The annotations are funannotate predictions, compared with RefSeq where a RefSeq annotation exists. RefSeq is not truth, and the scores come from
  exact-coordinate matching at conserved single-copy genes, which are the genes the analysis uses but not representative of whole genomes.
- Annotation error is measured, not removed: V7 shows which loss counts survive it and which do not. Terminal gains are less certain than
  losses, because false introns move them in no consistent direction.
- *S. cerevisiae* has so few introns that its scores rest on a few dozen introns and are noisy.
- Nothing here validates the evolutionary model itself (Dollo parsimony assumes one gain per site). Likelihood scoring is not implemented.

## Layout

```
crates/intronfate-core/  library: gff, annotation, seqio, msa, project, matrix, tree, score, tables
crates/intronfate-cli/   the `intronfate` binary
docs/                 USAGE, OUTPUTS, TESTING, SPEC (design + references), HANDOFF (state)
examples/             three runnable examples with committed expected output, plus 04_public_genomes (inputs for docs/HOWTO.md)
tests/data/           tiny synthetic fixtures
scripts/              helpers for docs/HOWTO.md (one_protein_per_gene.py, make_families.py, align_families.sh, doc_blocks.py);
                      malin_to_matrix.py — independent MALIN importer used as a cross-check
```

## Status

Implemented: GFF3/GTF ingestion with frame/exon/length cross-checks, protein-space projection,
FASTA/Stockholm alignments, site × taxon matrices, pairwise sharing, Dollo (classic + weighted),
Fitch MP with ancestral reconstruction, counts, TSV/NEXUS output.

Also: per-cell missing data (alignment gaps are `?`), MALIN import/export (`convert`), `--min-observed`.

Not yet: SQL/Chado/Ensembl sources (the `AnnotationSource` seam is there), likelihood-based
gain/loss models, benchmarks on real genomes (run those on the HPCC), stdin support.
