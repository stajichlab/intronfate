# What was tested about the placement of *Antarctolithica* (V4-V12)

*Ascomycota* sp. *Antarctolithica obscura* (BFD folders `Ascomycota_sp._Antarctolithica_obscura_Mars-*`, funannotate gene models, MAG assemblies with no phylum or class in the BFD taxonomy) has 4,232-5,668 predicted proteins per strain
and a long branch (0.71 substitutions/site in the OrthoFinder gene tree). The question: where does it belong, and do intron positions, protein-domain content and likelihood tests agree with a published CAT-model placement at the base of the Pezizomycotina?
Each step has its own README in `bench/`; this page is the summary. All analyses ran on the HPCC; the numbers below are copied from those READMEs.

## 1. Are the inputs trustworthy? (V4, V6, V7)
- **V4** compared funannotate introns with RefSeq (colleague's gffcompare tables): intron sensitivity is 89-99 % at the genes actually used. *S. cerevisiae* S288C genome-wide sensitivity looked low (63.3 %); the cause was **47 of 343 "reference introns" being 1-bp CDS gaps** (start == end, frameshift-type) that no annotation can match. Excluding them gives 73.3 % over 296 introns, the colleague's figure (`busco_gene_accuracy.py --min-ref-intron 10`). Recall by length is flat (about 75 % from 50 to 800 bp), so long introns are not being dropped; 36 novel predicted introns under 50 bp are probably false positives.
- **V6** built the single-copy families (OrthoFinder; synteny refinement added none). **V7** asked how much annotation error could change Dollo loss counts: terminal losses move by at most 6 % for five species, 15 % for *M. mucedo*, and all of the terminal losses of *C. cinerea* and *A. nidulans* can be explained by missed introns.

## 2. Intron positions (V8: 31 genomes; V9: 42 genomes)
- *Antarctolithica* is intron-poor (0.55-0.63 introns/protein). Whole-matrix parsimony and Dollo placement are **miscalibrated** for such a lineage: on all sites parsimony picks the Taphrinomycetes and Dollo picks *Yarrowia* or the basidiomycete outgroup; carrier-only parsimony picks *Neolecta*. In leave-one-lineage-out calibration on 20 known lineages (V9, tier A) parsimony recovers the true position for 8-9 and Dollo for 5-11 of them, with confident wrong answers, so none of these agrees-with-the-gene-tree results is evidence.
- The calibrated test is **clade-specific intron retention** (R = retention of a clade's defining introns / retention of ancient introns, normalised by the lineage's own intron richness). V9, with 12 Saccharomycotina genomes: *Antarctolithica* R is 0.14 for Pezizomycotina, 0.00 for Taphrinomycetes, 0.25 for Saccharomycotina and 0.29 for Saccharomycotina + Pezizomycotina, always in the non-member range (members: 1.08-2.37, 3.96-10.56, 1.10-1.88, 0.66-1.83). It is therefore **not inside** Pezizomycotina, Saccharomycotina, Lipomycetes or Taphrinomycetes, and its retention of the Saccharomycotina + Pezizomycotina introns looks like a non-member's. It does not locate the lineage among the earlier-diverging ascomycetes.
- Corrections made along the way (see the commit history): a re-rooting bug in my own check had produced a false "basidiomycetes nest inside Ascomycota" statement (retracted); the bootstrap penalised placements that had a tied twin tree (fixed with tie groups); V7 undercounted missed introns (pooled budget). Two independent read-only code reviews were run.

## 2b. Where its introns map, and a calibrated placement test (V13)
On the V9 tier-A matrix (42 genomes, 29,085 positions; Dollo and Fitch parsimony, no likelihood; rooted on *Mucor*):
- **Gains and losses.** The Alith stem has 183 gains and 4,185 losses, the largest loss count on the tree. Its ancestor held about 1,367
  positions (Dollo; MP 989) against 5,369 at its parent. Each strain carries 666-798 positions and the six together 1,208; the
  strain-to-strain spread is best read as MAG annotation noise.
- **Retention.** Alith keeps 8-11 % of the positions from every measurable age class (root to Ascomycota). For its intron count that is
  among the lowest, together with *Taphrina* and the black yeasts. Unrelated groups share this, so it is not a placement signal.
- **Placements.** Over all positions, the fewest losses put Alith next to *M. maydis*. Over the positions Alith carries, the fewest put it
  sister to Saccharomycotina + Pezizomycotina (the gene-tree position); sister to Pezizomycotina (+50), *Neolecta* + *Saitoella* (+53)
  and Saccharomycotina (+56) are close behind.
- **Calibration.** 15 known lineages were thinned to Alith's intron content and placed the same way. Over all positions, they were pulled
  to intron-poor genomes (54 of 75 age-matched runs), so the *M. maydis* result is that artefact. Over carried positions, they were never
  placed shallower than their true position, only deeper (50-53 of 75 runs) or correctly. So the carried-positions optimum is an upper
  bound: Alith is sister to Saccharomycotina + Pezizomycotina **or shallower** (for example sister to Pezizomycotina). The introns cannot
  separate these two, and they do not support a position deeper than Saccharomycotina + Pezizomycotina; 22 positions pointing only to
  *Neolecta* are the one contrary signal.

## 2c. The likelihood model with annotation error (V14)
Same matrix, fitted with `-m likelihood` (per-genome sensitivity and false-intron rates: measured for 6 genomes, the median for the
rest, Alith swept from 1.0 to 0.7):
- **Placement:** sister to Saccharomycotina + Pezizomycotina is best under every error set (family-resampled RELL 0.86-0.91), sister
  to Pezizomycotina second (Δ lnL 4.5-4.9, RELL 0.08-0.13), sister to Saccharomycotina third (Δ 5.1-5.2). Deeper positions are 36-47
  worse; next to *M. maydis* 186-197. Unlike Dollo over all positions, the likelihood is not attracted to intron-poor genomes. Not
  calibrated as V13's Dollo tests were.
- **Gains and losses:** the stem has about 6,560 expected losses [6,105-7,024] against Dollo's 4,185, because the likelihood
  reconstructs a richer parent (7,544 positions against 5,369). Alith's ancestor holds about 1,300 positions under both. Assuming
  Alith's annotation misses 30 % of introns moves its stem losses by less than 1 %: its intron poverty is not an annotation artefact.
- **Strains:** with errors modelled, almost no change is placed between the six strains (0-19 per branch), consistent with the
  strain-level Dollo events being annotation noise.
- No rate variation across positions yet, which can overstate ancestral richness and losses.

## 3. Protein-domain content (V10)
Pfam counts for 472 genomes (the V9 genomes plus a BFD sample of up to 3 per order; 10,642 Pfams). *Antarctolithica* carries 89-98 % of the 2,500 core fungal Pfams (other genomes: median 99.4 %), so it looks like a normal fungus. In the PCoA it lies between the yeast cluster (Saccharomycotina, Taphrinomycotina) and the Pezizomycotina; its nearest genomes are early-diverging yeasts (Lipomycetes-type Saccharomycotina, *Taphrina*, *Saitoella*), *Novakomyces* and *Tortispora*. A random forest that classifies known genomes at 99-100 % (leave-one-order-out) does not assign it confidently (Taphrinomycotina 0.41, Saccharomycotina 0.22, Pezizomycotina 0.18, Basidiomycota 0.17). Pfam signatures: it carries 0.11 of the Pezizomycotina signature (members 0.65-1.00) and 0.14 of the Saccharomycotina signature (members 0.53-1.00). Best reading: a yeast-like, small-proteome ascomycete outside both groups; but Pfam content tracks lifestyle and proteome size, so this is not a calibrated placement method.

## 3b. Genomic composition and whether it is a stand-alone group (V12)
Genome-level traits and composition for the same 472 genomes. Scalars: 7.8-10.5 Mb (Saccharomycotina median 14.6, Taphrinomycotina 13.0, Pezizomycotina 33.5), 4,232-5,668 genes, 541-572 genes/Mb (above almost every genome in every group), GC 49-50 % and GC3 52-54 % (typical, not diagnostic), 5-6 % masked. Nearest genomes differ by feature: amino-acid composition *Taphrina*, codon usage lichen-forming Pezizomycotina, Pfam content *Saitoella* and early-diverging yeasts. Isolation test (distance to the nearest genome outside the own order): typical for amino-acid composition (28 % of genomes are more isolated) and codon usage (63 %), but isolated for Pfam content (only 1.5 % of genomes are more isolated; the others are mostly reduced-genome lineages such as *Pneumocystis*). A random forest on amino-acid, codon and scalar features cannot place it (Taphrinomycotina 0.34, Pezizomycotina 0.26, Saccharomycotina 0.17) and cannot even recover known Taphrinomycotina (leave-one-order-out accuracy 0.06). Reading: composition does not mark *Antarctolithica* as a separate group and does not point to one group; the isolated Pfam repertoire together with the small compact genome and low intron density fits a reduced yeast-like lineage, which was not tested directly (for example whether its Pfams are a subset of its neighbours').

## 4. Likelihood tests of candidate placements (V11)
Supermatrix of 270 single-copy families (200,000 columns, 37 lineages, the six strains merged). Nine trees: unconstrained ML plus eight constraint hypotheses (sister to Pezizomycotina; to Saccharomycotina + Pezizomycotina; to Taphrinomycotina; to Saccharomycotina; to all ascomycetes; with *Saitoella* + *Neolecta*; with the Dothideomycetes black fungi; and a control with Taphrinomycotina monophyletic). IQ-TREE 2.2.2.6 for the searches and for the tree tests (RELL/KH/SH/ELW, AU); RAxML-NG 2.0.2 reproduced all nine log-likelihoods to within 0.04 (RAxML-NG has no AU/KH/SH tests and no mixture models, so it is a cross-check only).

| Tree | full, LG+F+G4 | fastest 75 % removed | PMSF LG+C60+F+G4, 50k slowest sites | Dayhoff-6 recoded |
|---|---|---|---|---|
| sister to Saccharomycotina + Pezizomycotina (= unconstrained ML) | 0 | 0 | 0 | 0 |
| sister to Pezizomycotina | 418 | 135 | 102 | 181 |
| sister to Saccharomycotina | 400 | 143 | 106 | 161 |
| sister to Taphrinomycotina | 2,804 | 765 | 650 | 1,758 |
| with the black fungi | 28,269 | 3,619 | 3,062 | 16,770 |
(delta lnL; the 25 % and 50 % removed runs agree). The best tree is the same in every variant. The gap to "sister to Pezizomycotina" shrinks with a more realistic model but never reaches zero; AU/KH reject it everywhere, the conservative SH test does not in the PMSF, site-removed and recoded runs. With this many sites the AU p-values are overconfident, so the delta lnL differences are used as a ranking under each model, not as evidence against a CAT result.

## 5. What can and cannot be said
- Consistent across three kinds of evidence (introns, Pfam content, likelihood under single-matrix and mixture models): *Antarctolithica* is **not a black yeast or a Pezizomycotina**, and not inside Saccharomycotina or the Taphrinomycetes. Its likelihood-best position is sister to Saccharomycotina + Pezizomycotina.
- **Not shown:** a decision between "sister to Pezizomycotina" (the published CAT placement), "sister to Saccharomycotina + Pezizomycotina" and a Taphrinomycotina-side position. The introns lean against being inside Saccharomycotina + Pezizomycotina; the likelihood tests lean toward sister to it; neither is decisive, and the long branch is the situation in which a single-matrix model is least reliable.
- **Limits:** funannotate gene models of MAGs (Mars-9.1 and Mars-19.38 lack finished annotation and are absent); the six strains are one lineage; the retention test has stated biases (members and non-members are scored on different site sets); the C60 test uses only the 50,000 slowest sites (the full matrix needs about 245 GB) with a guide tree from the single-matrix analysis; no full-data C60, partitioned or CAT-GTR (PhyloBayes) analysis on this matrix, and no posterior predictive tests; Pfam analysis is not calibrated for placement.
- **Next steps that would discriminate:** the same hypothesis set under CAT-GTR on this supermatrix with posterior predictive tests on the long branch; a topology test restricted to genes with the fewest gaps; the strain-by-strain check; the V7-style error sensitivity for the V9 set.

## Where things are
`bench/V4_annotation_accuracy`, `V6_orthology_nine_species`, `V7_annotation_error_sensitivity`, `V8_antarctolithica_placement` (31 genomes), `V9_antarctolithica_expanded` (42 genomes), `V10_antarctolithica_domains`, `V11_antarctolithica_topology_tests`, `V12_antarctolithica_composition`, `V13_antarctolithica_intron_map` (figures and the calibrated intron placement test), `V14_likelihood_alith` (the likelihood model with annotation error). Large data stays on the HPCC (purgeable folders are listed in each README); the repository holds the scripts, fixture tests and small result tables.
