# Tutorial: from genomes to intron gains and losses

This page explains the pipeline for your own genomes: what each step is for, what to choose, and what to check before you believe a
number. The commands themselves are in [`HOWTO.md`](HOWTO.md), which runs the same steps on public genomes and is tested as written;
follow it with your own files in place of the downloads. The worked numbers below come from the nine-fungus run in
[`bench/V6_orthology_nine_species`](../bench/V6_orthology_nine_species/README.md) (funannotate annotation).

What you will have at the end: for every lineage of your species tree, how many intron positions were gained and lost, with the table
of positions behind it. The heavy steps (OrthoFinder, MAFFT) belong on a cluster; `intronfate` itself takes seconds for a thousand
families.

## 0. Try the tool first

Install it ([README](../README.md#install)), then run `intronfate demo -o intron_demo`. The demo builds a five-protein toy family,
projects its introns onto an alignment and prints the Dollo and parsimony tables, so you can see every output before touching real
data. The hand-checked [`examples/`](../examples/README.md) work the same way.

## 1. What you need

For each species, a short **label** (no spaces, e.g. `Ncrassa`) and:

| File | Requirement |
|---|---|
| Annotation | GFF3 (or GTF) with `CDS` features grouped by transcript; gzip is fine |
| Proteins | FASTA of the proteins those transcripts encode, with ids that match the annotation (`protein_id=`, or the transcript id) |

and one **species tree** in Newick whose tips are exactly those labels, for example the taxonomy of the nine fungi:

```text
(Spunctatus,(Mmucedo,((Smaydis,(Ccinerea,Cneoformans)),(Spombe,(Scerevisiae,(Ncrassa,Anidulans))))));
```

`Smaydis` is *Mycosarcoma maydis* (*My. maydis* below); the label is historical and is kept so that results match earlier runs.
The tree is not inferred from the data, and Dollo places gains according to it, so use one you trust and root it (`--outgroup`).

Use annotation you trust: the result is only as good as the introns in it (§7). Keep protein ids unique across species
(funannotate's locus-tag prefix does this); `intronfate` reports collisions rather than guessing.

## 2. One protein per gene

OrthoFinder compares proteins, and two isoforms of one gene look to it like paralogs, which makes fewer families single-copy. So the
proteome given to OrthoFinder must have one protein per gene, named `<label>.faa`:

- NCBI (RefSeq or GenBank) and Ensembl protein files list every isoform. [`scripts/one_protein_per_gene.py`](../scripts/one_protein_per_gene.py)
  keeps the longest per gene, using the `protein_id` on the GFF3's CDS rows ([HOWTO](HOWTO.md) step 3).
- Gene predictors such as funannotate or BRAKER usually write one protein per gene already; copy that FASTA to `<label>.faa`, without
  terminal `*`. An internal `*` (a premature stop) makes the aligned row differ from the protein, so such proteins are left out of the
  families.

`intronfate` itself is later given the **original** annotation and protein files, not this reduced set: every transcript then joins its
protein, and the alignment picks out the one isoform per gene that was used.

## 3. Orthology

OrthoFinder groups the proteomes into orthogroups ([HOWTO](HOWTO.md) step 4). Give it the species tree (`-s`) so it does not have to
infer one. Keep `Orthogroups/Orthogroups.tsv` and `Phylogenetic_Hierarchical_Orthogroups/N*.tsv`.

**Why single-copy matters.** A shared intron is evidence about the history of a *gene* only if the proteins compared are orthologs,
not paralogs. So families are built with at most one protein per species, and multi-copy orthogroups are not guessed at: they are
either split with the hierarchical orthogroups or left out. For the nine fungi, 1,106 of 8,766 orthogroups were single-copy in all
species. Synteny-based refinement ([orthoSynAssign](https://github.com/stajichlab/orthoSynAssign)) is optional; across these distant
fungi it resolved almost nothing.

## 4. Families in tiers

[`scripts/make_families.py`](../scripts/make_families.py) turns the OrthoFinder tables into tiers of families ([HOWTO](HOWTO.md) step 5):

| Tier | Meaning | Nine fungi |
|---|---|---|
| `A_strict` | single-copy in every species | 1,106 |
| `B_partial` | single-copy in at least `--min` species, but not all | 689 (`--min 7`) |
| `C_hog` | from multi-copy orthogroups, split by the hierarchical orthogroups into single-copy groups | 401 |

Each tier has `families/<tier>/families/<id>.faa` and `taxon_map.tsv` (protein id to label). A species missing from a family is
simply absent: `intronfate` treats that as missing data, never as a loss. Start with `A_strict`; the other tiers add families at the
cost of more missing cells.

## 5. Align, and do not trim

Each family is aligned with MAFFT (`scripts/align_families.sh`, [HOWTO](HOWTO.md) step 6). **Do not trim** (trimAl, Gblocks): the
aligned row must equal the protein, otherwise it is rejected as `row_sequence_mismatch`, because projecting through a trimmed row would
put introns in the wrong column. Gappy, poorly aligned regions are handled by the tool: an absence there counts as *missing*, not as a
loss ([USAGE](USAGE.md) §6).

## 6. Run `intronfate`

One `intronfate analyze` call projects every intron, builds the matrices and scores them ([HOWTO](HOWTO.md) step 7). The options that
matter for your own data:

- `--taxon-map families/<tier>/taxon_map.tsv` maps each alignment row (a protein id) to its species, so the tree tips are species.
- The annotation and protein inputs take files, list files or folders (`--gff-list`, `--gff-dir`/`--gff-ext`, likewise for `-p` and
  `--msa`); a folder is not searched recursively, so with one folder per species use a list file ([USAGE](USAGE.md) §3).
- `--combine` stacks every family into one matrix and scores it once, so the totals per lineage need no summing; `--combined-only`
  skips the per-family folders.
- `--outgroup` roots the tree; Dollo placements depend on it.

For the nine fungi, tier A (1,106 families) runs in about 2 seconds with `--combined-only`.

## 7. Read the results

Start with `diagnostics.tsv` and `unmatched.tsv` (what failed to join, and why), then `combined/dollo_branches.tsv`, the gains and
losses per branch summed over all families. Every column of every file is described in [`OUTPUTS.md`](OUTPUTS.md).

Excerpt of `combined/dollo_branches.tsv` for the nine fungi (tier A; branches named by the tips below them):

| branch (tips below) | gains | losses |
|---|---|---|
| all but *S. punctatus* | 831 | 0 |
| *S. pombe*, *S. cerevisiae*, *N. crassa*, *A. nidulans* | 72 | 2,569 |
| *My. maydis* | 147 | 4,055 |
| *S. cerevisiae* | 17 | 1,447 |
| *C. cinerea* | 1,799 | 514 |

How to read it: a **gain** is where an intron position first appears (one origin per position, under Dollo); a **loss** is a
descendant that lacks it. Large losses on *My. maydis* and *S. cerevisiae* are the expected direction, since both are known for few
introns.

**Before you believe a number**, check three things:

1. **Diagnostics.** `row_sequence_mismatch`, `length_mismatch` or many `unmatched.tsv` rows mean the annotation and proteins do not
   correspond ([USAGE](USAGE.md) §8).
2. **Annotation error.** A missed intron looks exactly like a loss. Measure each annotation against a reference where one exists, over
   the genes you analysed. For the nine fungi, [V7](../bench/V7_annotation_error_sensitivity/README.md) shows which terminal-loss counts
   survive the measured error and which cannot be told from missed introns ([V4](../bench/V4_annotation_accuracy/README.md) has the
   accuracy figures).
3. **Terminal gains** are less certain than losses: false introns move them, in no consistent direction.

## 8. Where to go next

- **Re-score** without redoing the projection: `intronfate score --matrix results/A_strict/combined/matrix.tsv -t other.nwk --outgroup X
  -o rescored` tries another tree, rooting or method. `-m dollo_weighted` shows how the counts change when a second origin of a position
  is allowed at a price.
- **Model annotation error** with `-m likelihood --error-rates errors.tsv`: expected gains and losses per branch, with each genome's
  sensitivity and precision taken into account, and bootstrap intervals ([USAGE](USAGE.md) §11).
- **Compare tiers**: repeat with `B_partial` and `C_hog`; the loss pattern was the same for all three in the nine-species run.
- All options and input forms: [`USAGE.md`](USAGE.md). The coordinate model and assumptions: [`SPEC.md`](SPEC.md). MALIN users:
  `convert` imports and exports MALIN tables (USAGE §7).
