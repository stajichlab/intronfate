# Output reference

Every file is tab-separated text with one header row, stable row order (diffable between runs) and
absolute 1-based inclusive genomic coordinates. A leading `# key=value ...` line carries metadata where
a file has nowhere else to put it. Files are written atomically (temp file, then rename).

Conventions used below: **column** is the 0-based alignment column; **phase** is `0` (intron between
codons), `1` or `2` (after that many bases of the codon); a **site** is a `(column, phase)` pair, labelled
`<column>:<phase>`; `cds_offset` is the number of coding bases 5' of the intron, counted from the first
complete codon.

## Per run

| File | Content |
|---|---|
| `introns.tsv` | Every coding intron: `protein_id gene_id transcript_id seqid strand intron_rank coding_rank genomic_start genomic_end intron_len cds_offset phase codon residue_before residue_after site_key protein_len in_range`. `intron_rank` counts all introns of the transcript (UTR included), `coding_rank` only coding ones. `codon`, `residue_before/after` are 1-based; `site_key` is `cds_offset`. These numbers are **per protein** — compare species through `matrix.tsv`, not here. |
| `unmatched.tsv` | `transcript_id gene_id attempted_protein_id reason source` for every transcript that failed to join a sequence. |
| `diagnostics.tsv` | `kind count` per kind of problem, covering the annotation reader (`frame_mismatch`, `no_cds`, `strand_unknown`, `protein_id_fallback_transcript`, `intron_before_first_codon`, `transcript_id_collision`, `invalid_utf8`, …), the join (`protein_not_in_fasta`, `length_mismatch`, `site_beyond_protein_end`, …) and the alignment step. Written only when at least one warning or error was raised (informational records alone, such as the transcripts-retained summary, do not create it); a stale file from an earlier run is removed. |
| `diagnostics_detail.tsv` | Retained examples: level, kind, subject, message. |

## Per alignment: `<out>/<alignment>/`

| File | Content |
|---|---|
| `matrix.tsv` | Line 1 `# format=2 alignment= alen= phase_mode= flank= masked_present= taxa= sites= unannotated= residues=`. Header `site column phase n_taxa n_observable pattern <taxon>...`. Cells are `1` present, `0` observed absent, `?` not observed. `n_taxa` = carriers, `n_observable` = taxa with `1` or `0`. `pattern` is the taxon cells concatenated. This is the interchange format read by `score`. Files without `?` and `format=` (older versions) still load; an unannotated taxon's `0`s are read as `?`. |
| `sites.tsv` | Long form audit trail: one row per taxon × projected site with the full coordinate stack (`cds_offset codon residue_before residue_after site_label protein_len residues_aligned`). |
| `taxa.tsv` | Per taxon: `protein_id protein_len residues_aligned n_sites n_intron_records tail_sites clamped_sites annotated reason_if_not`. `reason_if_not` ∈ `protein_not_in_sequence_set`, `protein_not_annotated`, `row_sequence_mismatch`. |
| `presence.tsv` | `site column phase taxon`, presences only. |
| `pairwise.tsv` | Per taxon pair: `n_a n_b shared jaccard fraction_of_smaller both_annotated n_both_observable`. `shared` uses greedy one-to-one matching within `--tolerance` columns and requires equal phase unless `--ignore-phase`. `jaccard = shared/(n_a+n_b-shared)`; `fraction_of_smaller = shared/min(n_a,n_b)`. |
| `sharing_matrix.tsv` | Square matrix of `shared` counts; the diagonal is each taxon's site count. |
| `spectrum.tsv` | `n_taxa_with_site n_sites share_of_sites`: how many sites are carried by exactly *k* taxa. |
| `counts_sites.tsv` | Tree-free: `n_present n_absent n_missing frequency taxa`. `frequency = n_present/(n_present+n_absent)`; missing cells are excluded. |
| `characters.nex` | NEXUS discrete-character block (`0 1`, missing `?`) for external phylogenetics programs. |

## Scoring (`score`, `analyze`)

| File | Content |
|---|---|
| `dollo_sites.tsv` | `n_present n_missing gains losses monophyletic dollo_compatible gained_at lost_in`. Classic Dollo: one gain at the common ancestor of the carriers (`gained_at`), `losses` = minimum independent losses, `lost_in` the branches where they fall, separated by `;` (an internal branch is named by its tips, which are comma-separated). `monophyletic` = the observed carriers form exactly one clade (unobserved tips ignored); `dollo_compatible` = `losses == 0`. |
| `dollo_branches.tsv` | `branch node n_tips gains losses tips`: Dollo origins placed at, and losses on the branch leading to, every node, summed over sites — the table to compare with another tool's per-lineage losses (e.g. MALIN's). `branch` is `root`, a tip name, a non-numeric node label, or the node's tip list (`A,B`; long lists are truncated); where two nodes would print alike — or an internal node is itself labelled `root` — `#<node id>` is appended, so a name identifies exactly one node. `tips` is always the full comma-separated tip list. |
| `mp_sites.tsv` | `n_present n_missing steps informative ambiguous`: Fitch/Hartigan step count with gains and losses equal cost; `informative` = at least two taxa on each side; `ambiguous` = several equally short reconstructions. |
| `dollo_weighted.tsv` | (`-m dollo_weighted`) One row per weighting: `gain_weight loss_weight sites gains losses sites_one_gain sites_several_gains sites_with_losses`. As the gain weight rises, positions move from several independent origins to one origin with losses; where the totals stop changing, the loss counts no longer depend on the weighting. |
| `dollo_weighted_sites.tsv` | (`-m dollo_weighted`) Long form: `site gain_weight loss_weight gains losses cost` per site and weighting. Counts only: where on the tree the gains and losses fall is given by classic Dollo (`dollo_sites.tsv`). |
| `reconstructions.tsv` | (`--reconstructions`) `method site node is_tip state` per node. Internal nodes are named by their tip set; numeric internal labels in the input tree are treated as support values, not names. |
| `tree_gains.nwk` | The tree with `[gains=N]` on the branches where Dollo origins were placed. |
| `summary_score.tsv` | `key value`: method totals (`dollo_gains`, `dollo_losses`, `mp_steps`, `mp_consistency_index`, …), how many tree tips and matrix taxa matched, `characters_scored`, `characters_skipped` (by `--min-taxa` or no carrier), `characters_too_sparse` (by `--min-observed`), and whether the tree looked rooted. |

## Combined matrix: `<out>/combined/` (`--combine`)

Written by `matrix`/`analyze --combine`, or by `score --combine -o DIR`. It holds the matrix-level files of a single alignment
(`matrix.tsv`, `presence.tsv`, `pairwise.tsv`, `sharing_matrix.tsv`, `spectrum.tsv`, `counts_sites.tsv`; plus the scoring files for `analyze` and
`score`), computed over all families at once. Differences from a per-alignment folder:

- `matrix.tsv` line 1 has `alignment=combined families=<n> family_pad=<p>`, and the header is `site column phase family family_column n_taxa
  n_observable pattern <taxon>...`. `site` is `<family>:<family_column>:<phase>`; `column` is `family_column` plus the family's offset.
- In every other file the `site`/`label` is that same `<family>:<column>:<phase>` and `column` is the family-local column.
- `pairwise.tsv`, `sharing_matrix.tsv`: `n_a`, `n_b`, `shared`, `n_both_observable` are sums over families (sites of different families never pair).
- No `taxa.tsv`, `sites.tsv` or `characters.nex`: those are per-alignment projection records.

## MALIN export

`convert --to-malin` writes `#ORGANISM` lines, one row per taxon (`0`, `?`, and `3`/`1`/`2` for a present
intron of phase 0/1/2), then `#MAP site alignment alipos` with `alipos = 3·column + phase`.
