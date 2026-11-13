# V9: V8 repeated with 11 more genomes (Saccharomycotina, Taphrinomycotina)

V8 could not calibrate Saccharomycotina-specific introns (two genomes) and left the placement of *Antarctolithica* (Alith) relative to Saccharomycotina + Pezizomycotina
unresolved. V9 adds 11 genomes (`genomes.tsv`, 42 in total): Lipomycetes *Babjevia, Kockiozyma, Myxozyma, Limtongia, Dipodascopsis*, *Neoaidaea* (intron-rich, 1.5-2.2 introns/protein),
*Trigonopsis, Geotrichum, Alloascoidea* (Saccharomycotina, 12 in all), *Taphrina deformans* (`Tdeformans`) and *Schizosaccharomyces japonicus*, chosen from `../V8_antarctolithica_placement/screen_density.tsv`.
**Purgeable on the HPCC:** everything except the scripts and `results_summary/` (`families/ alignments/ results_projection/ trees*/ scan/ loo/ summary/ logs/`; the small tables are copied to `results_summary/`).

Method as V8: OrthoFinder 3.1.0 (46 min), families with a multi-copy genome treated as missing for that genome, the six strains one lineage; tier A = at least 33 of 37 lineages (1,704 families, 29,085 sites),
tier B = 29-32 (1,021 families); MAFFT untrimmed, `analyze --combine --combined-only`; backbone = OrthoFinder gene tree re-rooted on *Mucor* (`make_trees.py --root-on`, fixture-tested), and a variant with the Taphrinomycotina forced monophyletic.
Scans use the tie-group bootstrap and the calibration includes Dollo (both fixed after the code review of V8). Tests: `tests/`.

## Results (tables in `results_summary/`)
**Clade-specific retention** (`summary/clade_retention_{A,B}.tsv`; R as defined in `clade_retention.py`, with its stated biases), tier A:

| Clade-specific introns (sites) | members (n) | non-members (n) | Alith |
|---|---|---|---|
| Pezizomycotina (1,040) | 1.08-2.37 (14) | 0.00-0.19 (22) | 0.14 |
| Taphrinomycetes (458) | 3.96-10.56 (4) | 0.00-0.05 (32) | 0.00 |
| Saccharomycotina (547) | 1.10-1.88 (11) | 0.00-0.33 (25; *Neolecta* highest) | 0.25 |
| Lipomycetes (304) | 0.98-1.84 (7) | 0.00-0.72 (29) | 0.25 |
| Saccharomycotina + Pezizomycotina (1,736) | 0.66-1.83 (25) | 0.00-0.34 (11) | 0.29 |

Alith falls in the non-member range for every clade, with no member below it (tier B agrees: 0.10, 0.02, 0.13, 0.11, 0.32 against members' minimum 0.67-1.80; in tier B the Lipomycetes row is weak,
10 of 29 non-members score above Alith). With 12 Saccharomycotina genomes the Saccharomycotina-specific and the combined Saccharomycotina + Pezizomycotina introns are now calibrated, which V8 could not do,
so the V8 "unresolved" status of Sacch+Pezizo changes: Alith retains them at 0.29, below every one of 25 members (lowest 0.66) and inside the range of the non-members (highest *Neolecta*, 0.34).
Taphrinomycotina-specific introns (222 sites): Alith 0.09 against 0.2-0.6 for *Neolecta, Saitoella, S. pombe, S. japonicus* and 4.9-9.2 for the Taphrinomycetes, which is consistent with Alith not being a Taphrinomycete
but is not decisive for deeper Taphrinomycotina membership, because the intron-poor members themselves are low.

**Whole-matrix scans** (`scan/`) are again not usable: parsimony on all sites puts Alith in the Taphrinomycotina (0.999, tiers A and B), carrier-only parsimony on *Neolecta* (1.00 in A, 0.72 in B); Dollo on all sites prefers the
basidiomycete outgroup (0.70) in A and the Taphrinomycotina (0.99) in B, carrier-only Dollo the gene-tree position (1.00 in A, 0.75 in B; tied with a 25-genome mixed clade).
Calibration over 20 known lineages (`loo/summary_{A,B}.tsv`; true position = sister clade in the gene tree), tier A: parsimony 8/20 (all sites) and 9/20 (carrier sites); Dollo 11/20 and 5/20; tier B: parsimony 11/20 and 8/20, Dollo 14/20 and 8/20.
Because carrier-only Dollo recovers the truth in only 5/20 (tier A), its agreement with the gene tree for Alith is not evidence, as in V8.

## Conclusions and limits
- Together with V8: Alith is not inside Pezizomycotina, Saccharomycotina, Lipomycetes or Taphrinomycetes, and its retention of the introns that define Saccharomycotina + Pezizomycotina is in the non-member range.
  It is not placed more precisely: the intron-poor lineage also scores low for deeper Taphrinomycotina membership. The introns therefore agree with the gene-tree position (outside Saccharomycotina + Pezizomycotina) as far as they can discriminate; they do not by themselves locate it among the earlier-diverging ascomycete lineages or relative to *Neolecta*.
- Limits unchanged from V8: R is biased (members and non-members use different site sets), the lineage has about 1,700 carried introns, annotations are funannotate predictions, strains are merged into one lineage.
- *Cryomyces* is present in only 57 of 1,704 tier-A families (`summary/family_coverage.tsv`, 3 %; probably removed by the multi-copy rule, not checked), so it contributes nothing and the "black fungi" statements rest on *Hortaea* and *Zymoseptoria*.
- Not done: the strain-by-strain check, the V7-style error sensitivity for this set, the BFD species tree placement with an unassigned genus.
