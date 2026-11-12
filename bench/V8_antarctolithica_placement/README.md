# V8: do intron positions support the placement of *Antarctolithica*?

*Ascomycota* sp. *Antarctolithica obscura* (BFD folders `Ascomycota_sp._Antarctolithica_obscura_Mars-*`, funannotate gene models) has only 4,232-5,668
predicted proteins per strain, far fewer than the black fungi it might be assumed to belong to (*Cryomyces*: 17,347). The BFD genus-representative tree does not
contain it (its genus is unassigned). This asks what its intron positions say about where it belongs. **Purgeable on the HPCC:** `data/ proteomes/ results/
families/ alignments/ results_projection/ trees*/ scan/ loo/ logs/`; the small result tables are in `results_summary/`.

## Data and taxon set (`genomes.tsv`, 31 genomes)
Chosen with the project owner: Taphrinomycotina (6: *S. pombe*, *Taphrina*, *Protomyces*, *Savitreella*, *Saitoella*, *Neolecta*; *Pneumocystis* dropped),
Saccharomycotina (*Yarrowia*, *Lipomyces*; *S. cerevisiae* dropped), 14 Pezizomycotina from ten classes, outgroups *Mucor mucedo*, *M. maydis* and *C. neoformans*
(*Spizellomyces* dropped), and 6 *Antarctolithica* strains (Mars-19.38 and Mars-9.1 have no finished annotation). The backbone for the placement scans is OrthoFinder's gene tree for these genomes (with the strains removed). **Correction:** an earlier version of this README said the
BFD genus-representative tree, pruned to these genera and rooted on *Mucor*, nests the basidiomycetes inside the Ascomycota. That was wrong: it came from a bug in
my re-rooting code (found later by a fixture test in V9; the V8 scripts that were committed do not contain that code). Re-done correctly, the BFD tree gives
(*Mucor*, ((basidiomycetes), Ascomycota)) with the Taphrinomycotina as one clade sister to Saccharomycotina + Pezizomycotina (`../V9_antarctolithica_expanded/bfd_backbone_check.py`).
The OrthoFinder tree splits the Taphrinomycotina (*Saitoella* and *Neolecta* apart from the rest), so the scans were also run with the Taphrinomycotina forced to be one clade
(`trees_taphmono`); the conclusions did not change.

## Steps (all on the HPCC; scripts here)
| Step | Script | Result |
|---|---|---|
| link data, proteomes | `link_data.sh`, `prepare_proteomes.sh` | 31 proteomes |
| orthology | `run_orthofinder.sbatch` | OrthoFinder 3.1.0, no species tree given (30 min): 12,489 orthogroups; the inferred species tree roots on *Mucor*, puts the basidiomycetes outside the Ascomycota, and places the six strains as one tight clade on a long branch (0.71 substitutions per site) between (*Saitoella*, *Neolecta*) and (Saccharomycotina + Pezizomycotina) (`results_summary/gene_tree.nwk`) |
| families | `make_families_pruned.py`, `lineage_counts.py` | Only 16 orthogroups are single-copy in all 26 lineages (the strains count as one). Treating a genome with two or more genes in an orthogroup as **missing** (never a loss): tier A = at least 22 lineages, 2,317 families; tier B = 18-21 lineages, 769 families |
| align + project | `run_all.sbatch` | MAFFT untrimmed, `analyze --combine --combined-only`: one matrix per tier |
| placement scans | `make_trees.py`, `bootstrap_scan.py`, `run_scan2.sbatch` | the strains attached at every branch of the backbone (48 trees), scored by Fitch parsimony and Dollo, bootstrap over families |
| calibration | `run_calibration.sbatch`, `loo_calibration.py` | 14 known lineages treated as the unknown |
| clade-specific retention | `clade_retention.py`, `exclusive_sharing.py`, `intron_density.py` | see below |
| tests | `tests/` | fixture tests for the family builder and the retention test |

## Results
**Intron density** (`results_summary/intron_density.tsv`): *Antarctolithica* has 0.55-0.63 introns per protein (about a third of proteins have any); Pezizomycotina
1.2-2.5; the intron-poor controls *Savitreella* and *Yarrowia* 0.24, *Protomyces* 0.64, *Taphrina* 0.82.

**Whole-matrix parsimony is miscalibrated for an intron-poor lineage, so it cannot be used.** Attaching the strains at every branch of the backbone (scores in
`results_summary/scan/` and `A|B/placements.tsv`):
- on all sites, parsimony chooses the Taphrinomycetes (bootstrap over families 0.98 in tier A) and Dollo chooses *Yarrowia* (1.00): two intron-poor lineages, the two
  criteria disagree, and the pattern is shared loss;
- on sites where a strain carries the intron, parsimony chooses *Neolecta* (1.00 in tier A, 0.55-0.58 in tier B), the most intron-rich genome of all.

Calibration (`results_summary/loo/`): each of 14 known lineages was treated as the unknown. Both modes recover the truth for only part of them and misplace others with
bootstrap 1.00: *Yarrowia* and *M. maydis* are pulled into the Taphrinomycetes on all sites; *C. neoformans*, *Lipomyces* and *S. pombe* are pulled to *Neolecta*, and
*Botrytis*, *Orbilia* and *Pyronema* to *Geoglossum* (both intron-rich), on carrier-only sites. The two *Antarctolithica* results are the same artefacts.

**Clade-specific intron retention** (`clade_retention.py`, `results_summary/*/clade_retention.tsv`, fixture-tested) is calibrated on the known lineages. R = a
lineage's retention of the introns that define a clade / its retention of ancient introns (normalised for its own intron richness; tier A):

| Clade-specific introns (sites) | members | non-members | *Antarctolithica* |
|---|---|---|---|
| Pezizomycotina (~1,700) | 0.98-1.9 | 0.00-0.32 | **0.18** |
| Taphrinomycetes (245) | 7.0-9.0 | 0.00-0.12 | **0.04** |
| Taphrinomycotina (100) | 0.33-6.5 | 0.0-0.7 | **0.00** |
| Saccharomycotina + Pezizomycotina (~1,800) | 0.37-2.1 | 0.00-0.17 (*Neolecta* 0.32) | 0.30 |

## Conclusions, and their limits
- Intron positions **rule out** *Antarctolithica* inside Pezizomycotina (including Dothideomycetes and the black fungi) and inside the Taphrinomycetes. This agrees with its small
  proteome and with the gene tree, and it removes the whole-matrix "Taphrinomycetes" and "*Neolecta*" answers as artefacts.
- They **do not resolve** whether it sits just outside Saccharomycotina + Pezizomycotina (the gene tree) or inside it as a yeast-like branch: its retention of that group's
  introns (0.30) lies between the non-members and the members, and the clade cannot be calibrated for Saccharomycotina with two genomes.
- Power is limited by the lineage's own intron poverty (about 1,765 carried introns in tier A), and the annotations are funannotate predictions. The six strains were analysed as
  one lineage (their gene-tree branches differ by under 0.005 substitutions per site); a strain-by-strain check was not done.
- Not done: more Saccharomycotina and Taphrinomycotina genomes (to calibrate the Saccharomycotina question), the V7-style error sensitivity, and families below 18 lineages.
