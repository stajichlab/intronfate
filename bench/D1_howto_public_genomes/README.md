# D1: docs/HOWTO.md run as written, on both examples

Purpose: check that [`docs/HOWTO.md`](../../docs/HOWTO.md) works when its commands are pasted in order, on complete public genomes (NCBI
RefSeq), and record the numbers the page quotes. **Purgeable** on the HPCC: `repo/ containers/ logs/` (all work ran in `$SCRATCH`).

| File | What |
|---|---|
| `run_howto.sbatch` | `SET=fungi5` or `SET=fungi9`: runs every `sh` block of the page from step 1 on (`scripts/doc_blocks.py`), in `$SCRATCH`, then copies the small tables to `results/<SET>/` |
| `check_unmatched.sh` | what the `unmatched.tsv` and `no_cds` records of a RefSeq run are |
| `results/fungi5/`, `results/fungi9/` | `dollo_branches.tsv` (tiers A and B), `families_summary.tsv`, `diagnostics.tsv`, `summary_score.tsv`, `counts.txt` (sizes), `howto.out` (full log of the run) |

Only two lines of the page are replaced by the test: the working folder (→ `$SCRATCH`) and, for fungi9, the example line the page tells the
reader to change. `doc_blocks.py` fails if a replaced line is no longer in the page, so the test cannot drift from it. `apptainer` is
mapped to the site's `singularity` (`module load singularity`, SingularityCE 3.9.3). `repo/` is an rsync copy of the branch under test,
built with `cargo build --release`. The fixture tests (`scripts/test_howto_scripts.py`) run first in each job.

```sh
rsync -a --delete --exclude target --exclude .git ./ hpcc:projects/intron_projector/bench/D1_howto_public_genomes/repo/   # from the laptop
cd ~/projects/intron_projector/bench/D1_howto_public_genomes && (cd repo && cargo build --release) && cp repo/bench/D1_howto_public_genomes/run_howto.sbatch .
SET=fungi5 sbatch run_howto.sbatch; SET=fungi9 sbatch run_howto.sbatch
```

## Result (2026-11-13, jobs 32075194 and 32075195, node i33, 8 cores)

Both examples ran from download to the tier B rescore without error.

| | fungi5 | fungi9 |
|---|---|---|
| download (GFF3 + proteins) | 7 s | 13 s |
| one protein per gene | 6 s | 11 s |
| OrthoFinder 3.2 (container) | 593 s | 1,777 s |
| MAFFT, tier A | 166 s | 276 s |
| `intronfate analyze`, tier A | 4 s | 7 s |
| total, incl. tier B | 824 s | 2,223 s |
| disk (without the 650 MB image) | 373 MB | 717 MB |
| families A / B / C | 1,817 / 635 / 94 | 1,103 / 720 / 415 |
| intron positions (tier A, combined) | 15,140 | 15,283 |

- **One protein per gene**: the genes kept equal NCBI's protein-coding gene count for all nine genomes (e.g. *N. crassa* 9,757 of 10,812
  proteins; *C. neoformans* 6,975 of 7,826). No protein was dropped for an internal stop or a missing record.
- **Diagnostics** (`check_unmatched.sh`): the 101 `protein_not_in_fasta` transcripts (63 *A. nidulans*, 32 *S. pombe*, 6 *S. cerevisiae*)
  are pseudogenes (`pseudo=true` CDS rows, no protein); the ~10,000 `no_cds` are tRNA/rRNA/ncRNA. No `row_sequence_mismatch`, no
  `length_mismatch`.
- **fungi9 against the README's funannotate run** (V6 tier A): terminal losses *M. maydis* 3,864 (4,055), *S. cerevisiae* 1,452 (1,447),
  *M. mucedo* 1,479 (1,516), *C. neoformans* 1,010 (1,091), Ascomycota stem 2,455 (2,569); 1,103 strict families (1,106). Two independent
  annotations of the same strains give nearly the same loss pattern. Gains differ more (e.g. *C. cinerea* 1,855 here against 1,799), as
  expected from V4/V7.

## Problems found by running the page (all fixed in the page or the scripts)

1. **MCL check fails inside Singularity**: the cluster's `LD_LIBRARY_PATH` and conda variables leak into the container. The page now uses
   `apptainer exec --cleanenv`.
2. **FAMSA** (OrthoFinder 3.2's default MSA program, in the image) cannot run on the old AMD Abu Dhabi nodes of `batch` (c05). The page
   says to add `-A mafft` or use a newer node; this job asks for `epyc,intel`.
3. The page wrote OrthoFinder's output to a log in the working folder; when it failed, nothing was shown. It now `tee`s the output.
4. `align_families.sh` passed its arguments to `xargs` in the wrong order, and left a temporary file behind when MAFFT failed. Both fixed,
   with a stub-`mafft` test.
5. The rescore example used `-m dollo_weighted`, which the CLI does not offer (the README mentions weighted Dollo; it is in the library
   only). The page now uses `-m dollo,mp`.

The Docker option on the page was not tested.
