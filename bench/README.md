# Bench: validation runs on real data

Each folder holds the scripts, small result tables and a README for one question. All runs used the HPCC; large inputs and
intermediate files stay there (each README says which folders are purgeable). The numbers belong to the folder READMEs; the other
docs link here instead of repeating them.

| Run | Question | Result |
|---|---|---|
| [V1](V1_malin_validation/README.md) | Do scoring and MALIN interop reproduce MALIN's published test data? | Every figure checked reproduced. Per-branch Dollo losses match the MALIN guide exactly (Rogozin table: 7,221 sites, 2,770 losses); MALIN tables round-trip losslessly. |
| [V2](V2_gff_projection/README.md) | Does GFF3 → protein space → alignment columns behave on real annotation (nine fungal genomes)? | 739 BUSCO families (single-copy in at least 7 of 9 genomes), 6,363 protein rows, every row joined to its transcript (one refused for an internal stop codon); `analyze` over all families in 8.4 s and 507 MB. |
| [V3](V3_gff_scale/README.md) | Time and memory of `project` at genome scale | 1,000 genomes in one process: 294 s and 55.8 GB at the start, **52 s and 7.0 GB** after the optimisations; the largest GFF3 in the BFD (165 MB uncompressed) takes 1.6 s. A single file does not get faster with more threads. |
| [V4](V4_annotation_accuracy/README.md) | How accurate are the annotated introns at the genes actually used (funannotate against RefSeq)? | Intron sensitivity at those genes is 89–99 % for eight of the nine species, far better than genome-wide; *S. cerevisiae* is scored separately, over very few introns. Gives an analytic ceiling on the losses that missed introns could explain. |
| [V6](V6_orthology_nine_species/README.md) | Ortholog families that are as single-copy as possible (OrthoFinder, synteny refinement, tiers) | 8,766 orthogroups; 1,106 single-copy in all 9 species, 689 more in 7–8, 401 more from splitting multi-copy orthogroups. Synteny refinement added none. Tier A, combined: 15,697 intron positions; the largest terminal losses are in *Mycosarcoma maydis* (4,055), *Mucor mucedo* (1,516), *S. cerevisiae* (1,447) and *C. neoformans* (1,091). |
| [V7](V7_annotation_error_sensitivity/README.md) | How much could annotation error change the Dollo losses and gains? | Simulated at each genome's measured accuracy: terminal losses of *My. maydis*, *S. cerevisiae*, *C. neoformans*, *S. pombe* and *N. crassa* move by 0.1–6.4 % even if every missed intron hides a loss; *M. mucedo* by 15 %; for *C. cinerea* and *A. nidulans* all terminal losses could be misses. Gains move in no consistent direction. |
| [D1](D1_howto_public_genomes/README.md) | Does [`docs/HOWTO.md`](../docs/HOWTO.md) work when pasted in order? | Both examples ran from download to results (14 and 37 min on 8 cores). With RefSeq annotation the nine-fungus terminal losses are within 8 % of V6's funannotate run. |

Notes on the numbering:

- There is no V5. V8–V15 were analyses for a separate research project and moved to that project's own repository.
- The P-series (P1–P7, R1) are exploratory pilot analyses on the BFD annotation tables that scoped this work. They are one-off scripts,
  kept on `main` with their own index (`bench/PILOT_SERIES.md` there), and are not part of this validation set.

Conventions for a new run: a folder `<id>_<what>/` with a README that states the question, how to run it (`sbatch ...`), what is
purgeable on the HPCC, and the result with job ids and dates; only scripts and small tables are committed. Large temporary files go
to `$SCRATCH` inside the job.
