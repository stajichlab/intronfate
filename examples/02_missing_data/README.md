# Example 2 — a gap is missing data, not a loss

`A`, `B` and `C` share an intron **inside residue 6** (`cds_offset 16`, so codon 5, **phase 1**, alignment column 5).
Their exons are 16 + 20 nt — not multiples of 3, so the second exon's GFF `phase` is `2`
(`(3 − 16 mod 3) mod 3`), which exercises the phase rule. Taxon `D` is annotated too (its own intron is at
column 9) and the tree is `((A,B),(C,D));`. `run.sh` analyses two versions of `D`:

| Variant | D's alignment row | At column 5 |
|---|---|---|
| `gapped` | `MKVL--WEHTRS` (no residues in columns 4–5) | D cannot show the intron: cell `?` |
| `full` | `MKVLGAWEHTRS` (aligned straight through) | D is observed and lacks the intron: cell `0` |

## Expected, by hand

Site `5:1` is carried by A, B, C. Its common ancestor is the root.

- **gapped**: D is `?`, which costs nothing → **0 losses**. Site `9:0` is D's own (one gain) → total **2 gains, 0 losses**.
- **full**: D is observed absent → a loss on D's branch → **1 loss** (site `9:0` still adds one gain) → **2 gains, 1 loss**.

`matrix.tsv` for `gapped` shows the pattern `111?` for site `5:1`; `full` shows `1110`.
Without the missing-data rule both variants would report one loss.

Try `--flank 0`, which needs only the intron's own residue: the result for this example does not change
because D has no residue in column 5 at all. A taxon gapped only in a *neighbouring* column would
change from `?` to `0` — see `crates/intron-core/tests/per_cell_missing.rs`.
