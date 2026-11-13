# V10: is the Pfam domain composition of *Antarctolithica* consistent with fungi, and what does it suggest about placement?

Question from the project owner: *Antarctolithica* (BFD `Ascomycota_sp._Antarctolithica_obscura_Mars-*`; the assemblies are MAGs with no phylum or class in the BFD taxonomy) is not a black yeast, so what does its protein-domain content say?
No earlier PCoA/random-forest output was found in the repositories or on the HPCC; this is a new analysis. **Purgeable on the HPCC:** `selected_genomes.tsv`, `pfam_counts.tsv` (2.1 M rows), `logs/`; small tables and figures are in `out/`.

## Method (all on the HPCC; scripts here)
| Step | Script | Notes |
|---|---|---|
| genome set | `select_genomes.py` | 472 genomes: the 42 V9 genomes (found by the locus tag in their proteome ids; six annotated *Antarctolithica* strains, each kept separate) + a BFD sample of up to 3 genomes per order (fungal BUSCO complete >= 80 %, duplicated < 10 %, best-first, fixed seed) |
| Pfam counts | `pfam_counts.sql` (duckdb, `run_select.sbatch`) | BFD `tables/pfam.parquet`, domain i-Evalue < 1e-5, per genome and Pfam; 10,642 Pfams |
| analysis | `pfam_analysis.py` (`run_analysis.sbatch`; needs `module load gcc/12.2.0` for pandas) | PCoA (Bray-Curtis on relative counts, Jaccard on presence), nearest neighbours, random forest with leave-one-order-out CV, Pfam signature test, core-Pfam completeness |
Groups (`rf_group`): Pezizomycotina (215 genomes), Saccharomycotina (43), Taphrinomycotina (18), Basidiomycota (132), early-diverging groups; BFD leaves SUBPHYLUM empty, so ascomycete subphyla are assigned from the class (`pfam_analysis.py`); 6 genomes with missing or uncertain class are "unassigned" and not used for training. Pfams seen in < 2 % of the non-*Antarctolithica* genomes are dropped for the distances and the forest.

## Results (`out/`)
- **Consistent with fungi:** the strains carry 89-98 % of the 2,500 core Pfams (present in >= 95 % of the other genomes); other genomes carry a median 99.4 % (1st percentile 90.4 %). The lowest strain (Mars-14.71, 0.887) is also the least complete (BUSCO 78 %). 3,590-4,095 Pfam families per strain.
- **PCoA** (`pcoa_*.png`, `.tsv`; PC1/PC2 27 %/13 % Bray-Curtis, 29 %/11 % Jaccard): the strains lie between the yeast cluster (Saccharomycotina, Taphrinomycotina) and the Pezizomycotina cluster, apart from all clusters.
- **Nearest genomes** (`neighbours.tsv`): early-diverging Saccharomycotina (Lipomycetes-type: *Myxozyma, Kockiozyma, Dipodascopsis, Babjevia, Limtongia*), Taphrinomycotina (*Taphrina, Saitoella*), *Novakomyces*, *Tortispora*; no Pezizomycotina among the closest ten.
- **Random forest** (`rf_cv.tsv`, `rf_alith.tsv`): leave-one-order-out accuracy 99-100 % for the large groups (Kickxellomycota 0.67; 296 genomes with an *Antarctolithica*-sized repertoire: 0.993). The strains are not assigned confidently: mean probabilities Taphrinomycotina 0.41, Saccharomycotina 0.22, Pezizomycotina 0.18, Basidiomycota 0.17. The forest cannot place a profile that matches none of its groups.
- **Pfam signature test** (Pfams carried by >= 80 % of a group and <= 10 % of the rest; `pfam_signature.tsv`): Pezizomycotina (138 Pfams): members 0.65-1.00, *Antarctolithica* 0.11, only 12 of 251 non-members above it; Saccharomycotina (35 Pfams): members 0.53-1.00, *Antarctolithica* 0.14, 11 of 423 non-members above it. The Taphrinomycotina signature has fewer than 20 Pfams and is not tested.

## Reading
The domain content is a normal fungal one and is that of a yeast-like, small-proteome ascomycete outside the Pezizomycotina and outside crown Saccharomycotina, nearest to early-diverging yeast lineages. This agrees with the intron result (V8/V9) and is weak evidence against a position nested in either group.
It does not separate sister-to-Pezizomycotina from sister-to-(Saccharomycotina + Pezizomycotina) or a Taphrinomycotina-side position.
## Limits
- Pfam content tracks lifestyle and proteome size, so yeast-like neighbours can reflect convergent genome reduction instead of ancestry; this is not a calibrated placement method (unlike the intron retention test, no known-lineage calibration was done for the signature test beyond the member / non-member ranges shown).
- The strains are MAGs; Mars-9.1 and Mars-19.38 have no finished annotation and are absent. The BFD sample is not phylogenetically balanced beyond 3 genomes per order.
- Not done: a signature test for Taphrinomycotina, per-strain tests, domain architecture (instead of Pfam presence), and a topology test of the three candidate positions under a site-heterogeneous (CAT) model, which would be the discriminating analysis.
