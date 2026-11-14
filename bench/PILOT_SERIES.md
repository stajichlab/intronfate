# BFD pilot series (P1–P7, R1)

Exploratory analyses on the BFD fungal annotation tables that motivate and scope the intron-evolution work.
They are **one-off, purgeable scripts**, not library code: they read the shared BFD tables on the HPCC
(`/bigdata/stajichlab/shared/projects/BFD/...`, read-only) and write only into their own folder here. Only scripts,
sbatch files and small summary tables are tracked; large per-gene tables, logs, trees and alignments are not
(regenerate by rerunning the sbatch on the HPCC, partition `stajichlab` or `epyc`).

Most analysis scripts carry a self-test (grep for `selftest`); the summary, plotting and small helper scripts
(for example `p1_summary.py`, `p1e_summary.py`, `refmap.py`, `r1_figures.py`) do not. Migrating these to tested
code is still open.

| Folder | Question | Entry point |
|---|---|---|
| `P1_intron_density` | Genome-wide and marker-gene intron density per genome; assembly-quality floors | `p1_run.sbatch`, `p1_marker.sbatch` |
| `P1c_reference_introns`, `P1c_botrytis` | BFD intron calls versus RefSeq reference annotations; miss profiles | `p1c_run.sbatch`, `run_all.sbatch` |
| `P1d_version_audit` | Agreement with RefSeq by funannotate version | `p1d_run.sbatch` |
| `P1e_error_classes` | Gene-level error classes and causes | `p1e_run.sbatch` |
| `P1f_adjudication` | RNA-seq adjudication of disagreements | `p1f_run.sbatch` |
| `P2_genus_tree_check` | Taxonomy and support checks on the genus tree; FastTree vs ASTRAL | `tree_check.py`, `compare_trees.py` |
| `P3_reduced_mapping` | Intron-poor lineages mapped on the rooted tree; yeast-clade sensitivity | `reduced_rooted.py`, `yeast_clade.py` |
| `P6_topology_tests` | Gene-tree concordance and FastTree constrained-topology likelihoods | `topology_tests.sbatch` |
| `P7_lineage_choice` | Per-lineage marker-intron counts for choosing site-level loss lineages | `lineage_summary.sbatch` |
| `R1_report_figures` | Figures for the pilot report | `r1_figures.py` |

Decision recorded 2026-11-11 (owner): the BFD annotations are accepted as good; P1c–P1f agreement numbers are
sensitivity covariates, not a reason to rerun annotation.

Known caveats: `gene_introns.codon_position` in the BFD tables is off by one exon and must not be used; the
`Saccharomycotina sp. YJ-2025a/b` genomes (locus tags FF90B97F, F6EB82EF) are excluded from tree work.
