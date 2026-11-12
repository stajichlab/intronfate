# V6: orthology for the nine-species set (OrthoFinder, synteny refinement, tiered families)

Goal: ortholog families that are as single-copy as possible, for `intronfate analyze --taxon-map`. Same nine species as V2 except
*S. cerevisiae* is the reference strain S288C (BFD `Saccharomyces_cerevisiae_S288C`; V2 uses AMP_1a, which has no RefSeq annotation to score).
Everything runs on the HPCC; the laptop only holds the scripts and the small summaries in `reports/`. Purgeable on the HPCC:
`data/ proteomes/ bed/ results/ refined/ families/ alignments/ results_projection/ logs/`.

| Step | Script | Notes |
|---|---|---|
| link GFF3 + proteins | `link_data.sh` | read-only symlinks into BFD |
| proteomes | `prepare_proteomes.sh` | one `<label>.faa` per species, terminal `*` stripped |
| orthologs | `run_orthofinder.sbatch` | OrthoFinder 3.1.0, `-s species.nwk`, MSA gene trees; `VARIANT=y` adds `-y` (changes little) |
| synteny | `run_refine.sbatch`, `make_bed.py`, `pixi.toml` | [orthoSynAssign](https://github.com/stajichlab/orthoSynAssign) via pixi, windows 8 and 30 |
| families, align, project | `run_tiers.sbatch`, `make_families.py` | MAFFT `--auto`, then `analyze` per tier |
| summaries | `orthogroup_summary.py`, `summarize_tiers.py` | output copied to `reports/` |
| test | `tests/test_make_families.py` | 3-species fixture (`python3 tests/test_make_families.py`) |

Tiers (at most one protein per species; a missing species is missing data, never a loss; MIN = 7 of 9; ambiguous paralogs are in no tier):
A_strict = root orthogroups single-copy in all 9; B_partial = single-copy in 7-8; C_hog = from multi-copy orthogroups, hierarchical
orthogroups (`N1`, the clade without *Spunctatus*) single-copy in >= 7; D_synteny = synteny-refined groups single-copy in >= 7.

## Results (job chain 32032336, 32032337, 32032364)
- OrthoFinder: 8,766 orthogroups, 4,785 multi-copy (67 % of grouped genes), 12,996 genes unassigned; 1,106 strict single-copy, 1,795
  single-copy in >= 7/9. `N0` has no HOG file (the orthogroups are the root level); each `N<k>` covers only the species below its node.
- Synteny (window 8): 166 groups, 333 genes (pairs); window 30: none. Too little conserved gene order across these clades, so D_synteny is empty.
- Families: A 1,106 (9,953 proteins), B 689, C 401, D 0. `analyze` takes 10-26 s per tier. See `reports/`.
- Dollo gains/losses are summed per branch in `reports/tier_projection_REPORT.txt`. These are annotation-based counts with no correction
  for annotation misses (V4 gives per-species ceilings), so read the terminal-branch losses with V4 beside them.
- Diagnostics: `no_cds` warnings are mRNAs without CDS in the annotation (none is in a family); `proteins_not_in_alignment` is one record per alignment.
- Pitfall met: `rsync --exclude data` also skipped `tests/fixture/data`; anchor excludes (`/data`).

## Dollo losses vs. annotation error (`loss_ceiling_tiers.py`, `reports/loss_ceiling_tiers.txt`)
*Superseded by V4's rescoring over the V6 genes (`../V4_annotation_accuracy/results_v6/`, which also gives S288C a gene-restricted Sn of 97.4 %, not the 73.3 % used below) and by the simulation in `../V7_annotation_error_sensitivity/`.*
Terminal-branch losses compared with the V4 ceiling (`carried x (1 - Sn)/Sn`) and terminal gains with `carried x (1 - Pr)`, using V4's accuracy at the V2 BUSCO
genes (an approximation: the V6 gene sets overlap V2's but differ; S288C has only genome-wide Sn/Pr from 312 reference introns).
Same split in all three tiers: *N. crassa*, *C. neoformans*, *S. pombe*, *S. maydis* (and *S. cerevisiae*) losses are at most 4-8 % explainable by missed introns
(tier A: 3.9, 6.0, 3.8, 1.3, 1.2 %); *M. mucedo* 16-29 %; *C. cinerea* and *A. nidulans* are at or above 100 %, so their losses cannot be told from misses.
Gains are less safe: ceilings of 24-52 % for *A. nidulans*, *C. cinerea*, *M. mucedo*, *S. maydis* and *S. cerevisiae* (possible novel/false introns).
*S. cerevisiae* carries only 46 introns in tier A (1,447 terminal losses), matching its known intron poverty, but its Sn/Pr is the weakest measured.

## `--combine` and `--combined-only` on the tiers
`run_tiers.sbatch` runs `analyze --combine`: `results_projection/<tier>/combined/` holds all families in one matrix, and `check_combined.py` confirms its Dollo
table equals the per-family tables summed (0 differences in all three tiers; tier A: 15,697 sites, 17 branches). `check_combined_only.sbatch` re-runs tier A with
`--combined-only`: the six compared `combined/` files are byte-identical, the run writes 15 files instead of 15,499 (32 MB instead of 107 MB), and `analyze` takes 2.0 s
instead of about 8 s.
