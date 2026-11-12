# V2 — projection on real annotation (nine fungal genomes)

Purpose: V1 validated scoring and MALIN interop. This tests the tool's main job — GFF3 → protein space →
alignment columns — on real annotation, and gives the first real timings. **Purgeable** (`data/` holds only
symlinks; delete `results/ families/ alignments/ taxon_map.tsv`).

**Data.** Funannotate GFF3 + proteins from the BFD runs (`BFD/Fungi_BFD_runs/genome_annotation/<genome>/predict_results/`),
nine genomes chosen by the project owner (`genomes.tsv`): *Spizellomyces punctatus*, *Mucor mucedo*,
*Mycosarcoma (Ustilago) maydis*, *Coprinopsis cinerea*, *Cryptococcus neoformans*, *Schizosaccharomyces pombe*,
*Saccharomyces cerevisiae*, *Neurospora crassa*, *Aspergillus nidulans*. These are **funannotate predictions**, not
curated gene sets: expect differences from literature intron counts. Mycosarcoma is the current name of
*Ustilago maydis*. (*Rhizopus delemar* was offered as an alternative to *Mucor mucedo*; only one is used.)

**Tree.** `species.nwk` is the accepted taxonomy (Chytridiomycota outgroup; Mucoromycota; Basidiomycota
`(Smaydis,(Ccinerea,Cneoformans))`; Ascomycota `(Spombe,(Scerevisiae,(Ncrassa,Anidulans)))`) with **no branch
lengths** — it is not inferred from these data. Dollo placements depend on it.

**Orthologs.** BUSCO `fungi_odb10` in protein mode per genome; a BUSCO is a family when it is single-copy
("Complete") in at least `MIN_GENOMES` (7) of the 9 genomes, so some families lack a taxon, which exercises
missing data. Aligned with MAFFT `--auto`, **untrimmed** (gappy regions are handled by the observability rule,
`--flank`). Alignment row ids are the protein ids; `taxon_map.tsv` maps each to its genome so the species tree
can be used (`--taxon-map`).

Steps (`run.sbatch`): `prepare.sh` (link + count) → `busco_one.sh` per genome → `make_families.py` →
`align.sh` → `run_projection.sh` (per-genome `project` timing with peak RSS; one `analyze` over every genome
and family) → `summarize_projection.py` (REPORT.txt).

```sh
mkdir -p logs && sbatch -p highclock run.sbatch      # on 2026-11-11 `batch` was down and `epyc` had a month-long queue; highclock ran at once
```

## First result (2026-11-11, job 32023112, `highclock`, 34 s for the projection step; `REPORT.2026-11-11.txt`)

- 739 BUSCO families (single-copy in ≥ 7 of 9 genomes), 6 363 protein rows; every row joined to its transcript
  and only **one** was refused: a *C. cinerea* model with an internal stop codon, which MAFFT dropped from the
  row, so the row no longer equalled the protein (`row_sequence_mismatch`; the taxon was scored missing for that
  family instead of projecting through a shifted sequence). `make_families.py` now skips such proteins.
- No `length_mismatch`; 64 `site_beyond_protein_end` and 5 `tail_sites_dropped` across 739 families.
- Per genome, GFF parse + projection: 0.05–0.32 s and 18–125 MB peak RSS (largest: *C. cinerea*, 74 k CDS rows).
- All 739 families through `analyze` in one process: **8.4 s, 507 MB**. (The same data on a binary built before the
  per-alignment memory release, parallel pairwise, O(k) `mrca` and streaming writer fixes took 34 s and 3.7 GB; two
  changes, not isolated.)
- Biology agrees with the literature on every point checked by eye: phase 0 commonest (0.36–0.46) and phase 2
  rarest or near-rarest (0.22–0.28); *S. cerevisiae* has 407 introns in 335 intron-containing genes against ~5.5
  introns per intron-containing gene in *C. neoformans* and *S. punctatus*; Dollo losses are largest on the
  *S. maydis* (2 883), Ascomycota-stem (1 654) and *S. cerevisiae* (1 030) branches.
- **Annotation accuracy** (colleague's gffcompare scores vs RefSeq, measured over the genes in these alignments; `bench/V4_annotation_accuracy/`):
  intron sensitivity / precision at the V2 genes is 94–99 % for *C. neoformans*, *S. pombe*, *N. crassa* (99 / 99), 93–95 % for *M. mucedo*
  (Pr 89 %), *C. cinerea* and *S. punctatus*, and ~89 % for *A. nidulans* and *M. maydis*; errors do not depend much on phase or position.
  Annotation misses can explain at most a few percent of the terminal-branch losses of *N. crassa*, *C. neoformans*, *S. pombe*,
  *M. maydis* (and ~14 % for *M. mucedo*), but **not rule out** that all of the terminal losses of *C. cinerea* and *A. nidulans* are
  misses — treat those two as unreliable. *S. cerevisiae* (an assembly with no RefSeq annotation) is unscored.
- Caveats: funannotate predictions, not curated gene sets; a taxonomic (not inferred) tree; Dollo only; untrimmed
  alignments; `--flank 1` unvalidated. This shows the projection behaves sensibly on real data; it does not
  validate any particular count.

What a good result looks like: every BUSCO protein joins to its GFF transcript (no `protein_not_in_fasta`
for family members, no `length_mismatch`), no `row_sequence_mismatch`, phase 0 commonest, *S. cerevisiae* the
fewest introns and Dollo losses concentrated on the yeast branches. Anything else is a finding: check the
diagnostics before the biology.
