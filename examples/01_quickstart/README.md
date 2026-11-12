# Example 1 — quickstart

Four species, a 10-residue protein each (`MKVLGAWEHT`), aligned without gaps, on the tree
`((A,B),(C,D));`. Files: `family.gff3`, `proteins.faa`, `family.aln.fasta`, `family.nwk`. Run `./run.sh`.

## The annotation (so you can check it)

| Species | Strand | CDS segments (nt) | Introns at `cds_offset` |
|---|---|---|---|
| A | + | 12, 9, 9 | 12, 21 |
| B | + | 12, 9, 9 | 12, 21 |
| C | − | 21, 9 (transcript order: `601-621` then `401-409`) | 21 |
| D | + | 12, 18 | 12 |

`cds_offset` is the number of coding bases before the intron, in transcript order — so C's intron is
at 21 although it lies at *lower* genomic coordinates than its first exon. Every offset is a multiple
of 3, so all introns have phase 0: offset 12 is after codon 4 (column 4), offset 21 after codon 7
(column 7). Output: `introns.tsv` lists exactly these six rows.

## The matrix

```
site  column phase pattern  A B C D
4:0   4      0     1101     1 1 0 1      carried by A, B, D
7:0   7      0     1110     1 1 1 0      carried by A, B, C
```

## The scores, by hand

Each site is carried by three species and absent from the fourth, and no clade of the tree is
"three species", so neither is monophyletic:

- **Dollo**: site `4:0` — the carriers' common ancestor is the root; C lost it → 1 loss. Site `7:0` — D lost it
  → 1 loss. Total **2 gains, 2 losses**; neither site is Dollo-compatible.
- **Maximum parsimony**: one change per site (a loss on one branch) → **2 steps**.

`summary_score.tsv` reports `dollo_gains 2`, `dollo_losses 2`, `mp_steps 2`. `tree_gains.nwk` shows
both origins placed at the root.

## Things to try

- `--outgroup D`: roots the tree on D (`rooting outgroup:D` in `summary_score.tsv`). Dollo changes to
  **2 gains, 1 loss** — with D as the outgroup, `{A,B,C}` is a clade, so site `7:0` needs no loss and only
  `4:0` does — while MP stays at 2 steps, because Fitch does not depend on the root. Dollo placements
  depend on rooting; MP counts do not.
- `--min-taxa 4`: nothing survives (no site is carried by all four taxa) — `characters_skipped 2`.
- Without `-t`, only `-m counts` can run (`counts_sites.tsv`, `pairwise.tsv`); asking for Dollo or MP
  without a tree is an error that says so.
