# Handoff — state after the 2026-11-11 review and fix pass

**Nothing is committed.** The whole tree is untracked, so there is no baseline to diff against; the first
commit should be made once you have read this and are happy with it (the repo owner decides).

Snapshot: `cargo test --workspace` → **239 tests pass** (180 core unit, 23 core integration, 10 CLI unit,
23 CLI end-to-end, 3 doc); `cargo clippy --workspace --all-targets` → 0 warnings; `cargo fmt` applied;
`python3 -I scripts/test_malin_to_matrix.py` → 3 pass. Real-data validation (MALIN's published data, nine fungal
genomes, genome-scale timing, annotation accuracy, orthology tiers, error sensitivity) is in `bench/V1`–`V7`; per project
policy analysis runs on the HPCC. Inputs can be given as files, list files or folders (`--msa-list`/`--msa-dir` etc., USAGE §3).

Design rationale lives in [`SPEC.md`](SPEC.md). This file is state only.

## Second pass (same day): towards a usable tool

Added: native MALIN import/export (`convert`, `crates/intronfate-core/src/malin.rs`), `--min-observed`, the
`n_both_observable` pairwise column, atomic table writes, `rooting outgroup:` in `summary_score.tsv`, exit
code 1 for run failures (2 stays clap's usage error), the CLI `--transcript-attrs` default aligned with the
library (`Parent` first — the earlier fix had not reached the CLI), three runnable examples with committed
expected output (`examples/`), 13 binary-level tests, and the user documentation
(`docs/USAGE.md`, `OUTPUTS.md`, `TESTING.md`, example READMEs). Every claim in the example READMEs was run.

## What was done (first pass)

Three read-only reviews (a plan review and two code reviews) were run, then their confirmed findings were
fixed with a regression test each. The full list is in SPEC §9 ("Found by the multi-agent review"). The
headline changes:

1. **Per-cell missing data is implemented end to end** (SPEC §2). `Site` has `presence` + `observable`
   bitsets; `Cell {Present, Absent, Missing}`; `cell_observable()` applies MALIN's flank rule
   (`--flank`, default 1); `matrix.tsv` is `format=2` with `?` cells and an `n_observable` column;
   `score` reads `?` back; counts report `n_missing`; NEXUS writes `?`. Acceptance tests:
   `crates/intronfate-core/tests/per_cell_missing.rs`.
2. **Scoring correctness**: Fitch now handles polytomies (Hartigan counting); Dollo and Fitch are
   checked against brute force over every `{0,1,missing}` pattern; weighted Dollo uses 64-bit costs;
   `min_taxa` counts scored carriers; Dollo `monophyletic` ignores unobserved tips.
3. **Projection correctness**: 5'-partial CDS; NCBI and multi-parent GFF3; alignment rows that are not the
   joined protein are treated as missing; `Clamp` only clamps true tail overhangs.
4. **Robustness**: lossy UTF-8 reads, panic fixes, Newick comments/UTF-8, multi-outgroup rooting.
5. **Docs/spec** corrected and cited (SPEC §13); the GFF `phase` direction bug fixed earlier in this session.

## Decisions worth knowing

- `--flank` default **1** is an engineering default, not a validated one (MALIN's guide gives the rule, not a
  value). Tune against `Malin/test` on the HPCC before relying on loss counts.
- A present intron whose taxon is not aligned around it becomes `?` and is counted (`masked_present`), not kept
  as presence — otherwise ambiguous regions would be observed only when they show a gain.
- `strip_version` strips *all* trailing `.<digits>` segments; ambiguous stripped keys refuse to join.
- Default transcript key order is now `Parent`, `transcript_id`, ... (GFF3 first).

## Third pass: reviews and real-data setup

Three Opus reviews (tests via mutation testing, code, performance) and their fixes — see SPEC §9 and
TESTING.md. Highlights: unique branch names and gains by node id; `--ignore-phase` / `--clamp-tail`
observability fixed; exact `ambiguous` flag; `--taxon-map` (needed for one-alignment-per-family runs; the
error message had pointed at a flag that did not exist); randomised reference tests; and the cheap, safe
performance fixes, each checked output-identical:

| Fix | Why |
|---|---|
| `Tree::mrca` O(k·h) with precomputed depths; Dollo finds the MRCA and observed count in one postorder pass | was O(k²·h²) per site (measured 128 ms per call on a 200-tip caterpillar) |
| readable branch names computed once per tree from node ids | was `leaves_under` + formatting per loss per site |
| `build_matrices` processes and releases one alignment at a time | memory grew with the number of alignments (tens of GB at thousands) |
| `pairwise` builds per-taxon column lists and per-taxon observability bitsets once, runs in parallel, and the sharing matrix reuses the pairs | was O(T²·S) bit tests, run twice |
| `write_matrix_tsv` streams rows to a buffered file | `Table` made one `String` per cell (1e8 cells for 1e5 sites × 1e3 taxa) |

**V2** (`bench/V2_gff_projection`) **ran** (job 32023112, `highclock`; `batch` was down, `epyc` had a month-long
queue; BUSCO is the standalone prefix `/opt/linux/rocky/8.x/x86_64/pkgs/busco/5.8.0`, not the module): nine
fungal genomes, 739 BUSCO families, one `analyze` with `--taxon-map` in 8.4 s / 507 MB. Everything joined, one
row correctly refused (internal stop codon), phase and intron-density patterns match the literature. Two lessons:
rsync keeps laptop file times, so a job must `touch` the sources before `cargo build` or it can run a stale
binary (it did, once; fixed in `run.sbatch`); and `protein_not_in_alignment` was raised 58 million times in a
multi-family run (now one informational record per alignment). Results: `bench/V2_gff_projection/README.md`.

## Open / not done (in priority order)

0. **Annotation accuracy** (`bench/V4_annotation_accuracy/`): measured. At the V2 genes intron Sn/Pr (colleague's gffcompare vs RefSeq) are
   88–99 %, far better than the genome-wide rates; annotation misses can explain only a few percent of terminal-branch losses for
   *N. crassa*, *C. neoformans*, *S. pombe*, *M. maydis* (~14 % *M. mucedo*) but not exclude it for *C. cinerea* and *A. nidulans*.
   Still to do: sensitivity run on the V2 matrices, a direct funannotate-vs-RefSeq projection comparison, and an annotated replacement
   for the *S. cerevisiae* AMP_1a assembly (no RefSeq).
1. ~~V2 results~~ — read (see above). Next real-data step: a genome-scale GFF (hundreds of MB), more genomes
   and families, and a curated-annotation comparison (e.g. curated *S. cerevisiae* / *S. pombe* GFF3) to check
   the funannotate-vs-curated difference the V2 README flags.
2. **Genome-scale measurements exist** (`bench/V3_gff_scale`, 2026-11-11). The largest GFF3 in the BFD is 14 MB gz
   (165 MB plain, 108 k transcripts): 1.6 s / 222 MB. 1 000 genomes in one process: **52 s / 7.0 GB** (was
   294 s / 55.8 GB) after streaming `introns.tsv`, reading files in parallel, streaming the annotation into the
   join and slimming the join state. A single file still does not scale with threads. Remaining levers: one
   allocation per protein id (`Arc<str>`), parallelism inside one large file, `zlib-rs`.
3. **Performance backlog** (from the Opus review; none changes results; measure first, on the HPCC):
   word-parallel observability (`SiteMatrix::build` does S·T·(2·flank+1) bit tests and `observability()` one
   allocation per column); per-site `Vec<u8>` states in `prepare` → two bitsets; reusable buffers in
   `dollo_one`/`fitch_one` (`rayon::map_init`); streaming `reconstructions.tsv` (S·N rows); parallel GFF
   parsing and a faster gzip backend (`zlib-rs`); `read_tsv` and MSA reading without whole-file `String`s;
   a thread-count flag (the pool uses every core). Per-alignment
   parallelism (`par_iter` over alignments) is safe but needs a memory cap.
4. **A synthetic generator and a scaling benchmark** (the review's plan): seeded random trees, MSAs with
   indels and gap fraction, Dollo-evolved intron sites, GFF3 per species; axes T ∈ {10…1000}, sites per
   alignment {1e3…1e5}, alignments {1…5000}, species {1…100}, tree shape, flank, threads; metrics per stage
   wall time (needs `-v` stage timers), peak RSS, bytes written, GFF MB/s, and byte-identical output
   checksums as the regression gate. Not built yet; the V2 run gives only whole-process timings.
5. **MALIN FASTA headers** (`/organism=` + `{i … i}`) as an `AnnotationSource`; site-by-site comparison with
   MALIN's own output files.
6. **Dollo reporting**: `gain_ambiguous` when all-missing lineages sit beside the MRCA; stratify `compatible`
   by at least two carriers; per-branch "assessable sites" denominators.
7. Suspected, unconfirmed: a phase-0 intron right before the stop codon is dropped as out of range (needs a
   decision and a test); an intron inside a skipped 5' partial codon is warned about but not listed when the
   file has no exon rows; Stockholm `#=GS SF/SE` may not be a real convention (Pfam uses `name/start-end`);
   `FlyBase` `Derives_from` protein links; `Diagnostics` retains no errors after its cap; GTF `;` inside
   quoted values; a MALIN table using `1` for every presence (the guide says 1/2/3 all mean present) is rejected
   at phase-0 sites.

## Verify

```sh
cargo test --workspace && cargo clippy --workspace --all-targets && cargo fmt --check
python3 -I scripts/test_malin_to_matrix.py
cargo run -p intronfate-cli -- demo -o /tmp/ip     # expect 3 sites, 3 gains, 2 losses, 3 MP steps
cargo build --release && INTRONFATE=$PWD/target/release/intronfate examples/01_quickstart/run.sh
```

## Settled — do not relitigate without a reason

- Position identity is alignment **column + phase**; residue numbers are per-protein bookkeeping only.
- Missing ≠ absent: unannotated taxon, a gap around the site, or a row that is not the joined protein.
- Classic Dollo = one gain at the MRCA of carriers, minimise **losses**; weighted Dollo is the sensitivity tool.
- Re-rooting splits an **edge** and inserts a new root; MP is root-independent, Dollo placements are not.
- Third-party `Malin/` stays git-ignored; reference it by path, add no licensed fixtures to `tests/`.
