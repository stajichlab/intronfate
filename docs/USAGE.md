# intronfate — user guide

How to run the tool, what it needs, how to read what it writes, and what to do when something looks
wrong. For *why* it is designed this way see [`SPEC.md`](SPEC.md); for every output column see
[`OUTPUTS.md`](OUTPUTS.md); for worked, runnable examples see [`../examples/`](../examples/).

## 1. The idea in one paragraph

Genome annotation tells you where each intron is on a chromosome. To ask "do these orthologs have an
intron at the *same place*?" you need the position **in the protein**, and you need to compare it
across species **through a multiple alignment**. For every intron the tool computes `cds_offset` (how
many coding bases lie 5' of it), splits that into a codon and a phase (0, 1 or 2 bases into the
codon), looks up the alignment column of that codon in each ortholog, and builds a table of
`site × taxon`: `1` intron present, `0` absent, `?` cannot tell. It then scores the table on a tree.

## 2. Install

```sh
cargo build --release            # binary: target/release/intronfate
cargo test --workspace           # the full suite (see TESTING.md)
```

Requires Rust 1.85 or newer. There are no runtime dependencies.

## 3. What you must provide

| Input | Flag | Requirements |
|---|---|---|
| Annotation | `--gff` (repeatable; also `--gff-list`, `--gff-dir`) | GFF3 or GTF, optionally gzip. Needs `CDS` features grouped by transcript. `exon` features are used to find UTR introns. |
| Proteins | `-p` / `--proteins` (repeatable; also `--proteins-list`, `--proteins-dir`) | FASTA of the proteins those transcripts encode. |
| Alignment | `--msa` (repeatable; also `--msa-list`, `--msa-dir`) | Aligned FASTA or Stockholm of the orthologs. Row ids must equal the protein ids. |
| Tree (to score) | `-t` / `--tree` | Newick. Tip names must equal the alignment row ids exactly. |

### Many input files: lists and folders

A run over thousands of alignments or genomes does not need thousands of flags. Every multi-file input (`--gff`, `-p/--proteins`,
`--msa`, and `--matrix` for `score`) takes the same three forms, which can be mixed:

| Form | Example | Meaning |
|---|---|---|
| files | `--msa a.aln --msa b.aln` | exactly these files |
| list | `--msa-list families.txt` | one path per line; blank lines and `#` comments are skipped; a relative path is relative to the folder of the list file, not the working directory, so a list can travel with its data; a path that does not exist is an error naming the list and line |
| folder | `--msa-dir alignments/ --msa-ext aln.fasta` | the files directly inside the folder whose names match, sorted by name |

An extension (`--gff-ext`, `--proteins-ext`, `--msa-ext`, `--matrix-ext`; comma-separated or repeated) matches a file name that **is** it or
**ends in `.<ext>`**, ignoring case and a trailing `.gz`/`.bgz`: `aln` matches `OG1.aln` and `OG1.ALN.gz`, but not `OGaln` or `OG1.aln.bak`;
`aln.fasta` matches `OG1.aln.fasta`. Defaults when no extension is given: GFF `gff3, gff, gtf`; proteins `faa, fa, fasta, pep`;
alignments `aln, aln.fasta, aln.fa, aln.faa, afa, msa, sto, stockholm` (plain `.fasta` is deliberately not an alignment default,
because proteins use it too; a folder that mixes proteins and alignments needs an explicit `--msa-ext`). Folders are not searched recursively,
except `--matrix-dir` (below).

The result keeps the order *files, then each list's lines, then each folder's sorted matches*, and a path that occurs twice is used once.
Naming a list or folder that yields nothing is an error that says what was looked for (nothing named at all falls through to the usual
"`--msa` is required").

`score` also accepts many matrices: `--matrix-dir out/` searches the folder and its immediate sub-folders for `matrix.tsv` (or
`X.matrix.tsv`), so the output folder of `analyze` or `matrix` can be passed as it is. A single `--matrix FILE` writes its results into `-o`
as before; several matrices (or any list or folder form) are each scored into `<-o>/<family>/`, where `<family>` is the folder holding
`matrix.tsv` or the file name without `.matrix.tsv`; two matrices that would share a folder are an error.

Three requirements that cause most problems:

1. **The alignment row must *be* the protein.** The row's residues (gaps removed) must equal the
   sequence in the protein FASTA. A trimmed row (trimAl, Gblocks), a different isoform, or an edited
   sequence is reported as `row_sequence_mismatch` and the taxon is treated as **missing data**,
   because projecting through it would put every intron in the wrong column.
2. **Annotation and protein must correspond.** The transcript's protein id (`protein_id=` attribute,
   a `polypeptide` child feature, or `--id-map`) must name a sequence in the FASTA. Failures are listed
   in `unmatched.tsv`. The CDS length must be about `3 × protein length` (±3 residues by default,
   `--length-slack`); a larger gap is the signature of joining to the wrong paralog.
3. **Names must match exactly** between alignment rows, protein ids and tree tips. Nothing is fuzzily
   matched, because a wrong match would place introns on the wrong species.

### GFF conventions that are handled

`Parent=` is tried first, then `transcript_id` (so NCBI RefSeq, Ensembl, MAKER, Augustus/BRAKER and
protein2genome layouts all group correctly; change the order with `--transcript-attrs`). Multi-valued
`Parent=t1,t2` feeds every parent. 5'-partial gene models (first CDS with phase 1 or 2) are projected from
the first complete codon. Both strands, gzip, CRLF and stray Latin-1 bytes are fine. Use
`--derive-cds-from-exons` only for GFFs that carry exons but no CDS.

## 4. Commands

```
intronfate project  --gff a.gff3 -p a.faa -o out                    # introns -> protein space
intronfate matrix   --gff a.gff3 -p a.faa --msa fam.aln -o out      # -> alignment columns, site x taxon
intronfate score    --matrix out/fam/matrix.tsv -t fam.nwk          # -> Dollo / MP / counts
intronfate analyze  --gff ... -p ... --msa ... -t ...               # all of the above in one pass
intronfate convert  --from-malin t.table -o dir | --to-malin m.tsv... -o t.table
intronfate demo     -o intron_demo                                  # built-in worked example
```

`project` does not need an alignment. `matrix` writes everything up to the matrix; `score` re-reads a
`matrix.tsv`, so you can re-score under other rootings, filters and methods without re-parsing the
annotation. Give several `--gff`, `-p` and `--msa` for multi-species runs; ids must be unique across
files (collisions are reported).

**Each `--gff` file is treated as an independent annotation** (normally one species) and files are read in
parallel. A transcript id that appears in two files is therefore two records — never one chimera built from both —
and is reported as `transcript_id_collision`; ids should be unique across species (funannotate's locus-tag prefix
makes them so). Likewise each `-p` file is read in parallel and merged in the order given.

### Options that change the result

| Option | Default | Effect |
|---|---|---|
| `--ignore-phase` | off | Treat `codon 5 phase 0` and `codon 5 phase 1` as the same site. Off by default: they are one base apart but different sites, and merging them manufactures "conserved" positions. |
| `--flank N` | 1 | Residues a taxon must have on each side of a site before its *absence* counts. Fewer ⇒ more `0`s (and more losses); more ⇒ more `?`s. See §6. |
| `--clamp-tail` | off | Place an intron anchored past the end of a row on the row's last column instead of dropping it. Piles unrelated C-terminal sites together and is approximate by construction (the codon is not in the row); a clamped intron is exempt from the `--flank` window and counted in `clamped_sites`. Leave off unless your alignment was trimmed. |
| `--tolerance N` | 0 | Pairwise statistics only: match sites within ±N columns. The matrix itself always keeps exact columns. |
| `--length-slack N\|none` | 3 | Allowed difference between protein length and CDS/3 before warning. |
| `--strip-version` | off | Match `ENSP1.2` to `ENSP1`. Refuses a key shared by several sequences. |
| `--id-map FILE` | – | Two-column `annotation_id → protein_id` for stubborn joins. |
| `--taxon-map FILE` | – | Two-column `alignment_row_id → taxon`. **Needed for one-alignment-per-gene-family runs**, where each species' row carries its own protein id but the tree's tips are species. Rows not listed keep their id; two rows of one alignment may not map to one taxon. See §9. |
| `--strict` | off | Abort on a malformed annotation line instead of counting it. |
| `--outgroup NAME` | – | Root the tree on a tip (or a set that forms a split). **Dollo results depend on rooting**; a tree that looks unrooted triggers a warning. |
| `-m dollo,mp,counts` | all | Which scoring methods to run. |
| `--min-taxa N` | 1 | Skip sites carried by fewer than N scored taxa (`2` drops lineage-specific singletons). |
| `--min-observed N` | 0 | Skip sites observed (`1` or `0`, not `?`) in fewer than N taxa. |
| `--reconstructions` | off | Also write per-node ancestral states. |

## 5. Reading the results

Start with `out/<alignment>/summary_score.tsv`, then `dollo_sites.tsv`.

- **`matrix.tsv`** — the heart of it. `pattern` is one character per taxon in file order.
- **Dollo** assumes an intron at an exact position and phase arises **once** (at the common ancestor of
  the taxa that carry it) and is then lost independently any number of times. So `gains` is always 1 per
  site and the informative column is **`losses`**: the minimum number of independent losses. A site with
  `dollo_compatible = yes` needs none. Many losses means either that the intron really was ancestral and
  lost repeatedly, or that it arose more than once (parallel gain), or that the tree is wrong. Dollo is
  known to under-estimate ancestral intron density; treat it as one model among several.
- **Maximum parsimony** (`mp_sites.tsv`) allows gain and loss at equal cost. Compare it with Dollo on
  the same matrix: where they disagree, the Dollo assumption is doing the work.
- **Counts** (`counts_sites.tsv`) need no tree: per-site frequency among *observed* taxa, and the carriers.
- **`pairwise.tsv`** — shared sites, Jaccard and shared/min per taxon pair. Compare pairs only when
  `n_both_observable` is similar: a pair with little overlap looks dissimilar for the wrong reason.

## 6. Missing data: `?` is not `0`

A cell is **observed** only when the taxon is annotated *and* its alignment row has residues around the
site. Otherwise it is `?`, and scoring gives it no cost: **a gap never counts as an intron loss.**

- Phase 1/2 introns sit inside a residue: the taxon needs residues in columns `c-flank … c+flank`.
- Phase 0 introns sit between residues: the taxon needs residues on *both* sides, in adjacent
  columns (`c-m … c+m-1`, `m = max(flank,1)`). An intron beside an indel is therefore `?`, not placed on
  a different column in each species.
- An annotated intron whose taxon is not aligned around it becomes `?` too (counted in
  `masked_present` and reported as `present_not_observable`) — keeping only the presences in poorly
  aligned regions would bias the table towards gains.
- A taxon with no annotation for the gene is `?` at every site.

[`examples/02_missing_data`](../examples/02_missing_data/README.md) shows the difference: the same taxon
scores 0 losses when gapped and 1 when aligned straight through.

`--flank 1` is a sensible default, not a validated one. If your conclusions change between `--flank 0`,
`1` and `3`, they depend on how gappy the alignment is — say so when reporting them.

## 7. Using MALIN

[MALIN](https://www.iro.umontreal.ca/~csuros/introns/malin/) (Csűrös 2008) reads and writes tables of
one character per intron site. `convert` moves between the two:

```sh
intronfate convert --from-malin Chromo23.table -o matrices      # one matrices/<aln>/matrix.tsv each
intronfate score   --matrix matrices/KOG0004/matrix.tsv -t tree.nwk
intronfate convert --to-malin out/*/matrix.tsv -o ours.table    # taxa must be identical
```

A MALIN site is `(column, phase) = (alipos / 3, alipos % 3)`; phase digits are `3`=phase 0, `1`, `2`.
Export refuses matrices built with `--ignore-phase` (no phase to write), and import of a 0/1 table
without `#MAP` lines produces a phase-less matrix (`phase_mode=ignore`). See
[`../examples/03_malin`](../examples/03_malin/README.md) for a lossless round trip.

## 8. Troubleshooting

`diagnostics.tsv` counts every kind of problem; `diagnostics_detail.tsv` lists examples (`-v` prints some).

| Diagnostic | Meaning | What to do |
|---|---|---|
| `protein_not_in_fasta` | A transcript names a protein that is not in the FASTA. | Check ids; try `--strip-version` or `--id-map`. |
| `protein_id_fallback_transcript` | Transcripts with no protein id were joined on their transcript id (one informational record per annotation file, with the count). | Normal for funannotate and most GTF, whose protein FASTA ids are transcript ids; otherwise add `protein_id` or use `--id-map`. |
| `length_mismatch` | Protein length ≠ CDS/3 beyond the slack. | Wrong paralog or isoform? Check the join; raise `--length-slack` only if the annotation really omits stops/UTR. |
| `frame_mismatch` | CDS `phase` values disagree with the derived reading frame. | The CDS set is incomplete or inconsistent; offsets may be unreliable for that transcript. |
| `no_cds` | A transcript has no CDS and was dropped. | `--derive-cds-from-exons` if it is protein2genome-style. |
| `intron_before_first_codon` | An intron lies inside a 5' partial codon. | Informational; it is not projectable. |
| `isoform_intron_conflict` | Transcripts of one protein id differ in intron set. | Give isoforms distinct protein ids. |
| `site_beyond_protein_end` | An intron anchors past the protein's last residue. | Often a boundary intron before the stop codon. |
| `row_sequence_mismatch` | Alignment row ≠ protein sequence. | Use the untrimmed protein for the alignment, or rebuild it. |
| `proteins_not_in_alignment` | Some annotated proteins have no row in this alignment (one informational record per alignment: the count and the first few ids). | Expected when several species and families are run together; otherwise check row ids. |
| `unannotated_taxon` / `protein_not_annotated` | An alignment row has no annotation. | Its cells are `?`. |
| `tail_sites_dropped` | Introns anchor past the aligned part of a row. | Expected for trimmed alignments; `--clamp-tail` is the (blunt) alternative. |
| `present_not_observable` | Annotated introns where the taxon isn't aligned around the site. | Reported as `?`; inspect the alignment there. |
| `invalid_utf8` | Non-UTF-8 bytes in an input. | Read lossily; fix the file if the field matters. |

Exit status: `0` success, `1` the run failed (message on stderr), `2` command-line usage error.

## 9. Several species, many gene families

The usual real-data layout is one alignment per ortholog family, with each species' protein in it:

```
fam0001.aln.fasta   >Smaydis_g123-T1  >Scerevisiae_g88-T1  ...     (row ids = protein ids)
species.nwk         ((Smaydis,...),Scerevisiae)...                  (tips = species)
taxon_map.tsv       Smaydis_g123-T1<TAB>Smaydis   Scerevisiae_g88-T1<TAB>Scerevisiae ...
```

```sh
intronfate analyze --gff sp1.gff3.gz --gff sp2.gff3.gz ... -p sp1.faa.gz -p sp2.faa.gz ... \
    --msa-dir alignments/ --msa-ext aln.fasta \
    --taxon-map taxon_map.tsv -t species.nwk -o out
```

One process handles every alignment and writes `out/<alignment>/…`; alignments are processed one at a time and
released, so memory does not grow with the number of families. Per-family Dollo/MP results are then summed:
`dollo_branches.tsv` of each family gives losses per branch, and adding them gives the total per lineage.
A family that lacks a species simply has no row for it, and that tip is scored as missing.
[`bench/V2_gff_projection`](../bench/V2_gff_projection/README.md) builds exactly this from BUSCO orthologs.

## 10. Limits

**The result is only as good as the annotation's introns.** A real intron the gene model missed looks exactly like a
loss, and a spurious predicted intron looks like a gain, so per-genome intron sensitivity and precision (for example from
`gffcompare` against a reference annotation) bound how far to trust per-branch loss and gain counts — most of all on the terminal
branches of the least accurate genomes, since internal branches need several species to agree. Report the accuracy of the
annotations you analyse alongside the results. A worked example, with measured numbers for eight fungal genomes (accuracy over just the genes analysed, which can differ a lot
from the genome-wide rate), is in `bench/V4_annotation_accuracy/README.md`.

CDS features define protein space; UTR introns are reported but never projected. One protein per id
(isoforms are unioned and flagged). Likelihood-based gain/loss models and rate variation are not implemented (SPEC §5.5).
Analyses on real genomes belong on a cluster: nothing here has been benchmarked at genome scale yet.
