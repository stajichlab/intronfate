# V11: likelihood tests of candidate placements of *Antarctolithica* (IQ-TREE; RAxML-NG cross-check)

Why IQ-TREE and not only RAxML-NG: RAxML-NG has no AU/KH/SH tests and no site-heterogeneous mixture models; IQ-TREE has both (C60 mixtures with PMSF approximate the CAT model). RAxML-NG re-evaluates the same trees as a likelihood cross-check.
**Purgeable on the HPCC:** `supermatrix.*`, `hyp/`, `rx/`, `au.*`, `all_trees.*`, `logs/`; small results are in `results_summary/`.

## Method (all on the HPCC; `submit_all.sh` runs the chain)
| Step | Script | Notes |
|---|---|---|
| supermatrix | `build_supermatrix.py` (fixture test `tests/test_topology_setup.py`) | V9 tier-A alignments, the six strains merged into one lineage (`Alith`, per family the strain with most residues), families with >= 35 of 37 lineages, columns gapped in > 50 % of lineages dropped; 270 families, 200,000 columns, 7.6 % gaps (families added best-covered first up to 200,000 columns, so the 1,704 tier-A families are subsampled) |
| hypotheses | `make_hypotheses.py` | constraint trees H1-H8 (sister to Pezizomycotina; to Saccharomycotina + Pezizomycotina; to Taphrinomycotina; to Saccharomycotina; to all sampled ascomycetes; with *Saitoella* + *Neolecta*; with the Dothideomycetes black fungi; and H8, Taphrinomycotina monophyletic with Alith free, a control for the cost of the monophyly H3 needs) |
| ML searches | `run_hyp.sbatch` | IQ-TREE 2.2.2.6, LG+F+G4, one search per hypothesis plus the unconstrained tree |
| tests | `run_au.sbatch` | IQ-TREE `-z` (RELL: bp-RELL, KH, SH, wKH, wSH, c-ELW) and `-au` (10,000 replicates); RAxML-NG 2.0.2 `--evaluate` (LG+G4+F) on the same trees |

## Results (`results_summary/au.iqtree`)
| Tree | lnL | delta lnL | p-AU |
|---|---|---|---|
| unconstrained ML | -9,715,346.2 | 0 | 0.56 |
| H2 sister to Saccharomycotina + Pezizomycotina | -9,715,346.2 | 0.0004 | 0.55 |
| H8 control (Taphrinomycotina monophyletic, Alith free) | -9,715,489.8 | 143.6 | 0.012 |
| H4 sister to Saccharomycotina | -9,715,746.4 | 400.2 | 2e-4 |
| H1 sister to Pezizomycotina | -9,715,764.5 | 418.4 | 2e-52 |
| H5 sister to all sampled ascomycetes | -9,718,006.9 | 2,660.8 | 2e-50 |
| H3 sister to Taphrinomycotina | -9,718,150.5 | 2,804.3 | 1e-46 |
| H6 with *Saitoella* + *Neolecta* | -9,718,461.3 | 3,115.1 | 4e-54 |
| H7 with Dothideomycetes black fungi | -9,743,614.6 | 28,269 | 4e-4 |
The unconstrained ML tree places Alith as sister to the 25 Saccharomycotina + Pezizomycotina genomes (the same topology as H2; *Saitoella* and *Neolecta* are outside). RAxML-NG reproduces all nine log-likelihoods to within 0.04.

## Reading and limits
- Under this model the data prefer sister to Saccharomycotina + Pezizomycotina (as the OrthoFinder gene tree, and as the intron retention test, V8/V9) over sister to Pezizomycotina (delta lnL 418) and sister to Saccharomycotina (400).
- **The AU p-values are not to be taken at face value**: 200,000 sites and a single-matrix model (LG+F+G4, one partition) make every likelihood difference "significant" if the model is misspecified, and the long *Antarctolithica* branch (V8 gene tree: 0.71 substitutions/site) is the case where a single-matrix model is most likely to mislead. Site-heterogeneous (CAT-like) models can move exactly this placement, which is what the submitted CAT result is about. Use the delta lnL differences as a ranking under this model, not as evidence against the CAT placement.
- H8 itself is rejected (delta 144): a monophyletic Taphrinomycotina costs likelihood in this dataset, so H3's poor fit is partly that.
- Not done: LG+C60+F+G4 (PMSF) tests, per-gene or partitioned models, recoding or removing the fastest sites, using all 1,704 families, and the other strains.
