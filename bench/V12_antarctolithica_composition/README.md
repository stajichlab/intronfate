# V12: genomic composition of *Antarctolithica* and whether it is a stand-alone group

Uses BFD tables (`asm_stats`, `aa_freq`, `codon_freq`, `gene_info`, `samples`) for the 472 genomes of V10 (the V9 genomes plus a BFD sample of up to 3 per order) and V10's Pfam counts. **Purgeable on the HPCC:** `gene_counts.tsv`, `logs/`; small result tables and the log are in `out/`.
Scripts: `gene_counts.sql` (duckdb), `composition_analysis.py`, `run_composition.sbatch` (needs `module load gcc/12.2.0`). No fixture test: the script only joins tables and computes distances and a random forest, and it was checked by its own assertions (no NaN distances) and by the V10 analysis it extends.

## Results (`out/`)
**Scalar traits** (`scalars.tsv`, `scalar_percentiles.tsv`): genome 7.8-10.5 Mb (Saccharomycotina median 14.6, Taphrinomycotina 13.0, Pezizomycotina 33.5); 4,232-5,668 genes (Saccharomycotina median 5,923, Taphrinomycotina 5,550, Pezizomycotina 10,294); gene density 541-572 genes/Mb (above almost every genome in every group, medians 310-460); GC 49-50 % and GC3 52-54 % (typical, not diagnostic); 5-6 % masked (low in every group). The assemblies are MAGs.
**Nearest genomes** (`neighbours.tsv`): amino-acid composition (CLR): *Taphrina*, *T. deformans*, *Bifiguratus* (distance 0.24); codon usage (CLR): lichen-forming Pezizomycotina (*Lobaria*, *Arthonia*, *Xylographa*), 0.88; Pfam content (Jaccard): *Saitoella* and early-diverging Saccharomycotina (*Dipodascopsis*, *Kockiozyma*, *Babjevia*), 0.22. The three feature sets disagree.
**Isolation test** (`isolation.tsv`): distance from each genome to its nearest genome outside its own order (class); the *Antarctolithica* strains (nearest non-*Antarctolithica* genome) are placed in that distribution.

| Feature space | Alith nearest-genome distance | share of other genomes more isolated (order / class) | median / 95th percentile of the others |
|---|---|---|---|
| amino-acid composition | 0.237 | 28 % / 35 % | 0.169 / 0.355 |
| codon usage | 0.874 | 63 % / 72 % | 1.064 / 2.530 |
| Pfam content | 0.236 | **1.5 % / 3.9 %** | 0.116 / 0.204 |
The strains are close to each other in every space (mean pairwise distance 0.027, 0.128, 0.085). The genomes most isolated in Pfam space are mostly reduced-genome lineages (*Phlyctochytrium*, *Pneumocystis*, *Hyaloraphidium*).
**Random forest** (amino acids + codons + scalars, leave-one-order-out; `rf.tsv`): accurate for Pezizomycotina (0.99), Basidiomycota (0.92), Saccharomycotina (0.86), but **not for Taphrinomycotina (0.06) or Chytridiomycota (0.00)**, so composition does not identify those groups. *Antarctolithica* mean probabilities: Taphrinomycotina 0.34, Pezizomycotina 0.26, Saccharomycotina 0.17, Basidiomycota 0.08. The Taphrinomycotina value is not informative, because the classifier misassigns known Taphrinomycotina.

## Reading
- Amino-acid and codon composition are ordinary for a fungus and do not mark *Antarctolithica* as a separate group; they also do not point to one group (the nearest genomes differ by feature).
- Its Pfam repertoire is unusually isolated, and the most similar isolated genomes are reduced-genome lineages; together with the small, compact genome (few genes, very high gene density, few repeats, few introns in V8/V9) this fits a reduced-genome yeast-like lineage. Whether its Pfams are a subset of its neighbours' (reduction) or a different repertoire was not tested.
- Limits: compositional distances are sensitive to lifestyle, GC and proteome size, not only ancestry; MAG assemblies and funannotate models; the isolation test is relative to this sampling (3 genomes per order, with Antarctolithica's closest relatives possibly unsampled).
