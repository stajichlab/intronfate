# Tutorial: from genomes to intron gains and losses

This walks through the whole pipeline on your own species: prepare the data, find orthologs, align them, run `intronfate`, and read the
results. It follows exactly what was done for the nine-fungus set in [`bench/V6_orthology_nine_species`](../bench/V6_orthology_nine_species/README.md)
(scripts there are the reference implementation), so every number quoted below is real. Run the heavy steps (OrthoFinder, MAFFT) on a cluster or
a large machine; `intronfate` itself is light (the nine-genome, 1,106-family run takes seconds).

What you will have at the end: for every lineage of your species tree, how many intron positions were gained and lost, with the full table of
positions behind it.

## 0. Try the tool first (2 minutes)

```sh
cargo build --release                       # binary: target/release/intronfate
target/release/intronfate demo -o intron_demo
```

`demo` builds a five-protein toy family, projects its introns onto an alignment and prints the Dollo and parsimony tables, so you can see every
output before touching real data. The hand-checked examples in [`examples/`](../examples/) are runnable the same way.

## 1. What you need

For each species (a short **label**, no spaces, e.g. `Ncrassa`):

| File | Requirement |
|---|---|
| Annotation | GFF3 (or GTF) with `CDS` features grouped by transcript; gzip is fine |
| Proteins | FASTA of the proteins those transcripts encode, **one per gene**, ids that match the annotation (`protein_id=`, or the transcript id) |

and one **species tree** in Newick whose tips are exactly those labels, e.g. the taxonomy:

```
(Spunctatus,(Mmucedo,((Smaydis,(Ccinerea,Cneoformans)),(Spombe,(Scerevisiae,(Ncrassa,Anidulans))))));
```

Lay the data out as `data/<label>/<label>.gff3.gz` and `data/<label>/<label>.proteins.fa.gz` (the V6 helper scripts expect this).
Use annotation you trust: the result is only as good as the introns in it (see §6).

## 2. Prepare the proteomes

OrthoFinder wants one FASTA per species, named `<label>.faa`, with no terminal `*`. [`prepare_proteomes.sh`](../bench/V6_orthology_nine_species/prepare_proteomes.sh)
does this (and turns a rare internal `*`, a premature stop, into `X`):

```sh
mkdir -p proteomes
for l in Ncrassa Anidulans ...; do
  zcat data/$l/$l.proteins.fa.gz | sed 's/\*$//' > proteomes/$l.faa
done
```

Keep a protein id unique across species (funannotate's locus-tag prefix does this); `intronfate` reports collisions rather than guessing.

## 3. Orthology (OrthoFinder)

```sh
orthofinder -f proteomes -s species.nwk -t 16 -a 4 -o of_out -n run1
```

`-s` gives the species tree (OrthoFinder then does not have to infer it), and the default MSA-based gene trees are used. The nine-species run took
about 7 minutes. Keep `Orthogroups/Orthogroups.tsv` and `Phylogenetic_Hierarchical_Orthogroups/N*.tsv`.

Look at what you got with [`orthogroup_summary.py`](../bench/V6_orthology_nine_species/orthogroup_summary.py): of 8,766 orthogroups, 1,106 were
single-copy in all nine species, 1,795 in at least seven, and 4,785 were multi-copy (67 % of grouped genes).

**Why single-copy matters.** A shared intron is only evidence for homology of the *gene* if you compare orthologs, not paralogs. So families are
built with at most one protein per species. Multi-copy orthogroups are not guessed at: either split them using the hierarchical orthogroups, or leave them out.
(Synteny-based refinement, [orthoSynAssign](https://github.com/stajichlab/orthoSynAssign), is optional; across these distant fungi it resolved almost nothing.)

## 4. Build the families

[`make_families.py`](../bench/V6_orthology_nine_species/make_families.py) turns the OrthoFinder tables into tiers of families:

```sh
python3 make_families.py --of of_out/Results_run1 --data data --min 7 --out families
```

| Tier | Families | Meaning |
|---|---|---|
| `A_strict` | 1,106 | single-copy in all 9 species |
| `B_partial` | 689 | single-copy in 7 or 8 species |
| `C_hog` | 401 | from multi-copy orthogroups, the hierarchical groups that split into single-copy in at least 7 species |

Each tier is `families/<tier>/families/<id>.faa` plus `families/<tier>/taxon_map.tsv` (`protein_id<TAB>label`). A species missing from a family
is simply absent: `intronfate` treats that as missing data, never as a loss. Pick a tier (start with `A_strict`).

## 5. Align the families (untrimmed)

```sh
mkdir -p alignments/A_strict
for f in families/A_strict/families/*.faa; do
  b=$(basename $f .faa); mafft --auto --thread 1 --quiet $f > alignments/A_strict/$b.aln.fasta
done        # or xargs -P 16 / a job array
```

**Do not trim** (trimAl, Gblocks): the alignment row must equal the protein, otherwise the row is rejected as `row_sequence_mismatch`, because projecting
through a trimmed row would put introns in the wrong column. Gappy, poorly aligned regions are handled by the tool: an absence there counts as
*missing*, not as a loss.

## 6. Run `intronfate`

```sh
ls data/*/*.gff3.gz > gffs.txt            # one path per line; relative paths resolve against the list's folder
intronfate analyze \
    --gff-list gffs.txt --proteins-dir proteomes --proteins-ext faa \
    --msa-dir alignments/A_strict --msa-ext aln.fasta \
    --taxon-map families/A_strict/taxon_map.tsv \
    -t species.nwk -m dollo,mp,counts \
    --combine --combined-only -o results/A_strict
```

- `--taxon-map` maps each protein id (alignment row) to its species so the tree tips are species, not proteins.
- `--*-list` and `--*-dir/--*-ext` avoid one flag per file. A folder is searched non-recursively, so with the `data/<label>/` layout use a list
  for the GFFs (as above), or put all GFFs in one folder and use `--gff-dir`.
- `--combine` stacks every family into one matrix and scores it once; `--combined-only` skips the per-family folders (drop it if you want them).
- The nine-genome, 1,106-family run takes about 2 seconds and 165 MB with `--combined-only`.

## 7. Read the results

Everything is in `results/A_strict/` (details of every column: [`OUTPUTS.md`](OUTPUTS.md)):

| File | What it tells you |
|---|---|
| `combined/dollo_branches.tsv` | **The answer**: Dollo gains and losses summed over all families, per branch of the tree |
| `combined/dollo_sites.tsv` | One row per intron position (`<family>:<column>:<phase>`): carriers, gains, losses |
| `combined/matrix.tsv` | The site × species table (`1` present, `0` absent, `?` unknown) behind everything |
| `combined/pairwise.tsv` | How many positions each pair of species shares |
| `introns.tsv` | Every annotated coding intron, with its position in the protein |
| `unmatched.tsv`, `diagnostics.tsv` | What failed to join or was skipped, and why: read these first |

Excerpt of `combined/dollo_branches.tsv` for the nine fungi (families `A_strict`; `tips` shortened):

| branch (tips below) | gains | losses |
|---|---|---|
| all but *S. punctatus* | 831 | 0 |
| *S. pombe*, *S. cerevisiae*, *N. crassa*, *A. nidulans* | 72 | 2,569 |
| *M. maydis* | 147 | 4,055 |
| *S. cerevisiae* | 17 | 1,447 |
| *C. cinerea* | 1,799 | 514 |

How to read it: a **gain** is where an intron position first appears (one origin per position, Dollo); a **loss** is a descendant that lacks it. Large losses on
*M. maydis* and *S. cerevisiae* are the expected direction, since both are known for few introns.

**Before you believe a number**, check three things:

1. **`diagnostics.tsv`**: `row_sequence_mismatch`, `length_mismatch`, many `unmatched` rows mean the annotation and proteins do not correspond.
2. **Annotation error.** A missed intron looks exactly like a loss. For the nine fungi the missed-intron ceiling is 0.1-6.5 % of terminal losses for five species, ~17 % for
   *M. mucedo*, and as large as the losses themselves for *C. cinerea* and *A. nidulans*, which makes those two unreliable
   ([V4](../bench/V4_annotation_accuracy/README.md), [V7](../bench/V7_annotation_error_sensitivity/README.md)). Compare your annotation with a reference if one exists.
3. **Terminal gains** are less certain than losses (false introns inflate them).

To score other rootings, trees or methods without redoing the projection: `intronfate score --matrix results/A_strict/combined/matrix.tsv -t other.nwk --outgroup X -o rescored`.
To compare tiers, repeat §6 with `B_partial` and `C_hog`; the loss pattern was the same for all three in the nine-species run.

## 8. Where to go next

- All options and input forms: [`USAGE.md`](USAGE.md). The coordinate model and assumptions: [`SPEC.md`](SPEC.md). MALIN users: `convert` imports and exports MALIN tables.
- A worked, scripted run of everything above on a cluster, with logs and checks: [`bench/V6_orthology_nine_species`](../bench/V6_orthology_nine_species/README.md).
