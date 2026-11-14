# How-to: intron gains and losses from public genomes, step by step

Every command below can be copied and pasted in order. The data are complete public genomes (NCBI RefSeq annotation and proteins),
downloaded in step 2. Two examples use the same steps:

| Example | Genomes | Time, 8 cores | Disk |
|---|---|---|---|
| **1. five fungi** (`fungi5`) | *Cryptococcus neoformans* (outgroup), *S. pombe*, *S. cerevisiae*, *N. crassa*, *A. nidulans* | 14 min (OrthoFinder 10) | 0.4 GB |
| **2. nine fungi** (`fungi9`) | the five above, plus *Spizellomyces punctatus* (outgroup), *Mucor mucedo*, *Mycosarcoma maydis*, *Coprinopsis cinerea* | 37 min (OrthoFinder 30) | 0.7 GB |

Example 2 is the species set of the worked result in the [README](../README.md) (which used funannotate annotation; here RefSeq is used, so
the numbers differ a little, see step 8). Times and sizes are from the test run on a cluster node; the container adds 0.65 GB. The genome lists and trees are in [`examples/04_public_genomes/`](../examples/04_public_genomes/).

What you get at the end: for every branch of the species tree, how many intron positions were gained and lost, with the table of positions
behind it. For the reasoning behind each step, see [`TUTORIAL.md`](TUTORIAL.md). For every option, see [`USAGE.md`](USAGE.md).

**Use bash.** Keep one terminal open from step 1 to the end, because later steps use the variables set in step 1. If you come back
later, run the step 1 block again first. On a cluster, run the steps on a compute node, either interactively (`srun --pty -c 8 bash` or
your site's equivalent) or as one batch job ([at the end](#running-it-as-one-batch-job)).

## 0. Install (once)

**intronfate** is one Rust binary. If `cargo` is missing, install Rust from <https://rustup.rs>, or use conda: `conda install -c conda-forge rust`.

```sh
git clone https://github.com/stajichlab/intronfate.git
cd intronfate && cargo build --release && cd ..
REPO=$PWD/intronfate
export PATH=$REPO/target/release:$PATH
intronfate --help | head -3
```

The other tools are **OrthoFinder** and **MAFFT**, plus `python3` (3.8 or newer, standard library only) and `curl`. There are two ways
to get OrthoFinder and MAFFT. Choose one in step 1:

- **A container (recommended on a cluster)**: if `apptainer` (or `singularity`) is available, nothing else needs installing. The
  [BioContainers](https://biocontainers.pro) OrthoFinder image already includes MAFFT, DIAMOND and the other tools OrthoFinder needs.
- **conda/mamba**: install both into one environment:

```bash
mamba create -y -n orthofinder -c conda-forge -c bioconda orthofinder=3.2 mafft
```

## 1. Pick an example, a folder and the tools

```sh
mkdir -p ~/intronfate-howto && cd ~/intronfate-howto
SET=fungi5; MIN=4; OUTGROUP=Cneoformans
CPUS=8
EX=$REPO/examples/04_public_genomes
```

For Example 2, change the second line to `SET=fungi9; MIN=7; OUTGROUP=Spunctatus`. `MIN` is the fewest species a partial family may
have (tiers B and C in step 5).

Then set `RUN`, the prefix that runs OrthoFinder and MAFFT. **With the container** (downloads a 650 MB image once):

```sh
[ -s orthofinder_3.2.sif ] || curl -fL -o orthofinder_3.2.sif https://depot.galaxyproject.org/singularity/orthofinder:3.2--hdfd78af_0
RUN="apptainer exec --cleanenv -B $PWD,$REPO $PWD/orthofinder_3.2.sif"
$RUN orthofinder -h > /dev/null && echo "OrthoFinder runs"
```

Where the command is called `singularity`, write `singularity exec` instead of `apptainer exec`. Keep `--cleanenv`: without it, the
cluster's `LD_LIBRARY_PATH` and conda settings leak into the container, and OrthoFinder then fails its MCL check. You can also build the image yourself
with `apptainer pull orthofinder_3.2.sif docker://quay.io/biocontainers/orthofinder:3.2--hdfd78af_0`.

**With conda** instead:

```bash
conda activate orthofinder
RUN=""
```

**With Docker** (a laptop without Apptainer; not tested here, and slow on Apple-silicon Macs, which emulate the x86 image):

```bash
RUN="docker run --rm -u $(id -u):$(id -g) -v $PWD:$PWD -v $REPO:$REPO -w $PWD quay.io/biocontainers/orthofinder:3.2--hdfd78af_0"
```

## 2. Download the genomes

The GFF3 annotation and the protein FASTA of each RefSeq assembly in the list:

```sh
mkdir -p genomes
grep -v '^#' $EX/$SET.tsv | while IFS=$'\t' read -r label asm organism; do
  url=https://ftp.ncbi.nlm.nih.gov/genomes/all/${asm:0:3}/${asm:4:3}/${asm:7:3}/${asm:10:3}/$asm/$asm
  echo "$label  $organism"
  curl -fsSL -o genomes/$label.gff.gz         ${url}_genomic.gff.gz
  curl -fsSL -o genomes/$label.protein.faa.gz ${url}_protein.faa.gz
done
ls -l genomes
```

The label (first column) is the species name used from here on: proteome file names, tree tips and output columns.

## 3. One protein per gene

RefSeq protein files list every isoform. OrthoFinder would treat two isoforms of one gene as paralogs, so keep the longest per gene:

```sh
mkdir -p proteomes
for g in genomes/*.gff.gz; do
  label=$(basename $g .gff.gz)
  python3 $REPO/scripts/one_protein_per_gene.py --gff $g --proteins genomes/$label.protein.faa.gz --out proteomes/$label.faa
done
```

Each line of the output reports how many genes were kept, for example:

```text
Ncrassa: 9757 genes kept from 10812 proteins (1055 other isoforms, 0 internal stop, 0 named in the GFF but not in the FASTA, 0 in the FASTA but not in the GFF)
```

The number of genes kept equals NCBI's count of protein-coding genes for every genome in both examples.

## 4. Find orthologs (OrthoFinder)

```sh
$RUN orthofinder -f proteomes -s $EX/$SET.nwk -t $CPUS -a 4 -o orthofinder -n $SET 2>&1 | tee orthofinder.log
ls orthofinder/Results_$SET/Orthogroups/Orthogroups.tsv
```

`-s` gives OrthoFinder the species tree, so it does not have to infer one. OrthoFinder refuses to write into an existing `-o` folder.
To run this step again, first `rm -r orthofinder`.

If OrthoFinder stops with `Cannot run MSA method 'famsa'`, the CPU is too old for FAMSA, its default aligner (seen on 2012-era AMD
cluster nodes). Add `-A mafft` to the command, or run on a newer node.

## 5. Build single-copy families

```sh
python3 $REPO/scripts/make_families.py --of orthofinder/Results_$SET --proteomes proteomes --min $MIN --out families
```

This prints the families in each tier. Every family has at most one protein per species:

| Tier | Meaning |
|---|---|
| `A_strict` | single-copy in every species: **start here** |
| `B_partial` | single-copy in `MIN` or more species (but not all); a species without a protein is *missing data*, never a loss |
| `C_hog` | from multi-copy orthogroups, split by OrthoFinder's hierarchical orthogroups into single-copy groups |

Each tier has `families/<tier>/families/*.faa` and `families/<tier>/taxon_map.tsv`, which maps each protein to its species.

## 6. Align each family

```sh
$RUN bash $REPO/scripts/align_families.sh families/A_strict/families alignments/A_strict $CPUS
```

This runs MAFFT `--auto` on each family. **Do not trim the alignments** (trimAl, Gblocks). `intronfate` needs each aligned row to be exactly
its protein, and it refuses a row that is not (`row_sequence_mismatch`). Poorly aligned regions are dealt with in step 7: there, an
intron absence counts as missing, not as a loss.

## 7. Project the introns and score gains and losses

```sh
intronfate analyze \
    --gff-dir genomes --gff-ext gff \
    --proteins-dir genomes --proteins-ext protein.faa \
    --msa-dir alignments/A_strict --msa-ext aln.fasta \
    --taxon-map families/A_strict/taxon_map.tsv \
    -t $EX/$SET.nwk --outgroup $OUTGROUP -m dollo,mp,counts \
    --combine --combined-only -o results/A_strict
```

- The **original downloads** are given here (every isoform, the full GFF3). That lets every transcript join its protein, and the
  alignment picks out the one isoform per gene that was used.
- `--taxon-map` makes the tree tips species rather than proteins.
- `--combine --combined-only` stacks all families into one matrix and scores it once. Drop `--combined-only` to also get one folder per family.
- This step takes 4 s (Example 1) or 7 s (Example 2).

## 8. Read the results

```sh
column -t results/A_strict/combined/dollo_branches.tsv | cut -c1-150
cut -f1,2 results/A_strict/diagnostics.tsv
```

Example 1 gives (columns `branch`, `gains`, `losses`; `node`, `n_tips` and `tips` left out):

```text
branch                                gains  losses
root                                  1935   0
Cneoformans                           8342   0
Spombe,Scerevisiae,Ncrassa,Anidulans  344    0
Spombe                                923    965
Scerevisiae,Ncrassa,Anidulans         17     444
Scerevisiae                           23     1747
Ncrassa,Anidulans                     1105   21
Ncrassa                               894    679
Anidulans                             1557   164
```

How to read it: a **gain** is where an intron position first appears (Dollo parsimony allows one origin per position). A **loss** is a
descendant that lacks it.

The outgroup's own branch (here *C. neoformans*, 8,342 gains) holds every position found only in the outgroup. With a single outgroup these
cannot be dated, so do not read them as lineage-specific gains. Add a second outgroup (Example 2) to see more of the deep history.

In Example 2, the RefSeq annotation gives almost the same terminal losses as the funannotate run behind the README's worked result:

| Branch | Losses, this how-to (RefSeq) | Losses, README (funannotate) |
|---|---|---|
| *M. maydis* | 3,864 | 4,055 |
| *S. cerevisiae* | 1,452 | 1,447 |
| *M. mucedo* | 1,479 | 1,516 |
| *C. neoformans* | 1,010 | 1,091 |
| Ascomycota (*S. pombe*, *S. cerevisiae*, *N. crassa*, *A. nidulans*) | 2,455 | 2,569 |

with 1,103 strict single-copy families (README: 1,106). The full tables of both examples are in
[`bench/D1_howto_public_genomes/results/`](../bench/D1_howto_public_genomes/results/).

The other tables:

| File | What it tells you |
|---|---|
| `combined/dollo_branches.tsv` | Dollo gains and losses per branch, summed over all families (above) |
| `combined/dollo_sites.tsv` | one row per intron position (`<family>:<column>:<phase>`): which species carry it, gains, losses |
| `combined/matrix.tsv` | the position × species table (`1` present, `0` absent, `?` unknown) behind everything |
| `combined/pairwise.tsv` | intron positions shared by each pair of species |
| `introns.tsv` | every annotated coding intron of every genome, with its position in the protein |
| `unmatched.tsv`, `diagnostics.tsv` | what failed to join or was skipped, and why. **Read these first.** |

All columns are explained in [`OUTPUTS.md`](OUTPUTS.md). With RefSeq, expect these diagnostics, which are harmless:

| Diagnostic | Example 1 | Why |
|---|---|---|
| `no_cds` | 10,132 | tRNA, rRNA and ncRNA transcripts: they have no CDS |
| `protein_not_in_fasta` (in `unmatched.tsv`) | 101 | pseudogenes (`pseudo=true`): RefSeq gives their CDS rows but no protein |
| `proteins_not_in_alignment` | 1,817 | one record per alignment: the proteins (other genes, other isoforms) that are not in that family |
| `present_not_observable` | 624 | introns in a stretch where that species is not aligned (gaps around the site): scored `?`, not counted |

There are also a few dozen minor records (`site_beyond_protein_end`, `cds_gap_not_exon_gap`, `overlapping_exons`, `frame_mismatch`, ...; see
[`USAGE.md`](USAGE.md) for each). Neither example had a `row_sequence_mismatch` or `length_mismatch`; many of those would mean that the
annotation and the protein file do not match.

**Before you interpret a number:** RefSeq annotation is not truth. A missed intron looks exactly like a loss, and a false intron looks like a
gain. For how large that effect is in these species, see [V4](../bench/V4_annotation_accuracy/README.md) and
[V7](../bench/V7_annotation_error_sensitivity/README.md).

## 9. Next steps

**Score again** without redoing the projection (this takes seconds). Change `-t`, `--outgroup` or `-m` to try another tree, rooting
or method. `dollo_weighted` asks how the counts change when a second origin of a position is allowed at a price (gain:loss costs
1:1, 2:1, 5:1 and 10:1; set others with `--gain-loss-weights`):

```sh
intronfate score --matrix results/A_strict/combined/matrix.tsv -t $EX/$SET.nwk --outgroup $OUTGROUP -m dollo,mp,dollo_weighted -o results/A_strict/rescored
column -t results/A_strict/rescored/dollo_weighted.tsv
```

**Use another tier.** Align and run it the same way, for example tier B:

```sh
$RUN bash $REPO/scripts/align_families.sh families/B_partial/families alignments/B_partial $CPUS
intronfate analyze --gff-dir genomes --gff-ext gff --proteins-dir genomes --proteins-ext protein.faa \
    --msa-dir alignments/B_partial --msa-ext aln.fasta --taxon-map families/B_partial/taxon_map.tsv \
    -t $EX/$SET.nwk --outgroup $OUTGROUP -m dollo,mp,counts --combine --combined-only -o results/B_partial
```

**Use your own species.** Write a list in the same format as [`fungi5.tsv`](../examples/04_public_genomes/fungi5.tsv) (label, then the
RefSeq `GCF_...` accession with its assembly name, as in the folder name under
`https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/`), and write a Newick tree with the labels as tips. Then follow the steps unchanged.

- A **GenBank** (`GCA_`) assembly with annotation works the same way.
- For **Ensembl** genomes, download the GFF3 and the `pep.all` FASTA. `one_protein_per_gene.py` reads Ensembl's `protein_id`.
- For **your own annotation** (funannotate, BRAKER, MAKER), see [`TUTORIAL.md`](TUTORIAL.md) §1-2. If the protein FASTA already has one
  protein per gene, copy it into `proteomes/` instead of running step 3.

## Running it as one batch job

On a cluster, the steps can run unattended. `doc_blocks.py` writes the `sh` blocks of this page, from step 1 on, into one script.
The alternatives shown in `bash` blocks (conda, Docker) are left out, so the script uses the container:

```bash
python3 $REPO/scripts/doc_blocks.py $REPO/docs/HOWTO.md --from "## 1." > howto.sh
sbatch -c 8 --mem 16G --time 2:00:00 --wrap "REPO=$REPO; PATH=$REPO/target/release:\$PATH; source howto.sh"
```

This is how the page itself is tested ([`bench/D1_howto_public_genomes`](../bench/D1_howto_public_genomes/README.md)).
