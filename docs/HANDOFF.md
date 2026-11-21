# Handoff: state on 2026-11-16

This page is state only: what the branch holds, what is open, and how to check it. Design and rationale are in [`SPEC.md`](SPEC.md),
user documentation in [`USAGE.md`](USAGE.md), validation runs in [`../bench/README.md`](../bench/README.md).

## Branch

`docs/howto-public-genomes` is **not yet merged into `main`**; merging is the owner's decision. It carries everything below.

*Naming:* the project was formerly `intron-projector`. Working folders on the laptop and the HPCC are still named `intron_projector`
(e.g. `~/projects/intron_projector/bench/...` in the bench READMEs); only the repository, crates, binary and environment variable
(`INTRONFATE`) changed name.

## Implemented

- **Projection**: GFF3/GTF → protein space → alignment columns, with frame/exon/length checks, per-cell missing data (`?`) and the
  flank rule; file, list-file and folder inputs; `--taxon-map` for one-alignment-per-family runs.
- **Combining**: `--combine` / `--combined-only` stack all families into one matrix scored once (equal to the per-family sum, tested).
- **Scoring**: `dollo`, `mp` (Fitch/Hartigan), `dollo_weighted` (`--gain-loss-weights` sweep), `counts`.
- **Likelihood** (`-m likelihood`, SPEC §5.5–5.5b): per-genome annotation error (`--error-rates`, sensitivity with precision or false
  rate), bootstrap intervals (`--ml-bootstrap`, family blocks), gamma rate classes (`--ml-rate-categories`, `--ml-gamma-shape`), EM
  stopping rule as options (`--ml-tol`, `--ml-param-tol`), and a per-cycle trace with `-vv`. Two independent reviews; fixes tested.
- **MALIN** import/export (`convert`), validated on MALIN's published data (V1).
- **Docs**: [`HOWTO.md`](HOWTO.md) on public genomes, run as written on the HPCC for both examples ([D1](../bench/D1_howto_public_genomes/README.md));
  the doc set was restructured on 2026-11-16 so each topic has one home (install → README, outputs → OUTPUTS, options → USAGE,
  numbers → the bench READMEs, tests → TESTING).

Analyses for a separate research project (formerly `bench/V8`–`V15`) moved to that project; nothing here depends on them.

## Open items

In rough priority order. Items marked SPEC are described there.

1. ~~Separate gain and loss rate factors~~ and ~~free class levels like MALIN's~~ — done (SPEC §5.5c, §5.5d; `--ml-gain-categories` etc., `--ml-gain-levels`, `--ml-loss-levels`, `--ml-free-root-branch`).
2. ~~Validate the likelihood fit against MALIN's `*.rates`~~ — done ([L1](../bench/L1_malin_likelihood/README.md)): constant rates (within 0.002), Rogozin rate variation (within 0.07), and `Chromo23` rate variation, which reproduces MALIN's printed log-likelihood (-28430.9249) to four decimals at MALIN's parameters once the same never-seen correction is used (`--ml-never-seen global`); this tool's maximum is 0.009 higher.
3. ~~EM `--ml-max-iter` overshoot~~ — fixed: the limit is exact (plain EM steps near the limit; test `max_iter_is_an_exact_limit_on_e_steps`).
4. ~~`--flank` default~~ — examined ([F1](../bench/F1_flank_tuning/README.md)): MALIN's test data has no alignments, so it was done on nine fungal genomes; the absence rate of conserved introns does not depend on the distance to a gap, so a larger flank only discards data. Kept at 1.
5. **SQL/Chado/Ensembl annotation sources** (`AnnotationSource` seam exists; needs a schema decision; SPEC §12 item 4).
6. **HOWTO Docker route is untested** (the page says so; D1 used Apptainer/Singularity and conda is not run by the test either).
7. **MALIN FASTA input** is built (`--malin-fasta`, SPEC §6c; equivalent to the GFF route on the quickstart). Open: a site-by-site
   comparison with a MALIN run on the same files, which needs a MALIN run supplied by the owner (SPEC §12 item 2).
8. **Scaling**: a synthetic generator and `matrix`/`analyze` scaling in taxa and sites; remaining memory levers for `project` (one
   allocation per id, parallel parsing inside one file, `zlib-rs`) are in [V3](../bench/V3_gff_scale/README.md).
9. Python bindings only if needed (stdin inputs and `--stdout TABLE` are done, USAGE §3) (SPEC §12 items 7–8).
10. **Dollo reporting ideas**: `gain_ambiguous` when all-missing lineages sit beside the MRCA; per-branch "assessable sites"
    denominators.
11. **Suspected bugs, checked** (each with a test): resolved: GTF `;` inside quoted values (attributes are split quote-aware); a 0/1 MALIN
    table with `#MAP` lines (MALIN treats `1`, `2`, `3` alike; the phase comes from `alipos` unless the table uses `2`/`3`); an intron inside
    the skipped 5' partial codon, which without exon rows vanished from the intron list and shifted later `intron_rank`s (fixed). Not bugs:
    an intron between the last residue and the stop codon is reported (`site_beyond_protein_end`, `tail_sites_dropped`) and set aside by
    design; `--clamp-tail` places it one codon early (documented as approximate); Stockholm `SF/SE` is this tool's own convention, now
    stated as such (SPEC §3, USAGE). Still unchecked: FlyBase `Derives_from` protein links (no FlyBase file at hand).
12. **Soft-evidence likelihood was wrong before the seen-conditioning fix (SPEC §5.6b)**: results from `--ml-soft-evidence` / `--ml-origins` with soft cells computed earlier are invalid (ln L could exceed 0 and gains collapsed). Fixed and checked against a brute force; the J1 benchmark was rerun with it (`bench/J1_junction_support`).
13. **Junction support as soft evidence** (`intronfate junctions`, USAGE §4, SPEC §5.6b) is built, unit/CLI tested, and run on real RNA-seq
    ([J1](../bench/J1_junction_support/README.md): six genomes, hisat2 BAMs, normalised shallow reads). [J2](../bench/J2_support_vs_reference/README.md)
    compares read support with the RefSeq annotation: one read is the biggest step (agreement 67–81% with none, 97–99% with one in *C. neoformans*, *S. pombe*,
    *N. crassa*; disagreement still halves from one read to two). [J2b](../bench/J2b_matrix_cells/README.md) restricts to the matrix cells: only 4–9% lack reads, and their agreement is near the default `p0`, the default map (`p0` 0.8, `k` 3) gives only 0.85 at one read, and `k` near 0.2 would fit these data; RefSeq is a weak
    label in *A. nidulans* and *C. cinerea*. Defaults not changed (owner's decision). RefSeq has too few disagreeing cells in the reliable genomes to test the calls (J2b step 3). Open: non-normalised RNA-seq, and junctions at coordinates the annotation lacks
    (not projected).
    **Simulation ([S1](../bench/S1_soft_evidence_simulation/README.md)):** on the real matrix's sites with simulated truth and errors, soft junction evidence never improved the
    convergent-insertion ranking (AUC lower than the hard fit in every seed; no calls at p > 0.5), and it cannot flag all-error sites; the hard fit with the true error rates
    is the better remedy. Treat the soft likelihood as a sensitivity analysis, not an improvement.
14. **Four-class gamma fits** are slow to converge, not on a flat shape ridge ([J3](../bench/J3_shape_ridge/README.md): the profile is peaked; about 2 ln L
    units remain at the 20,000-iteration cap). [J3b](../bench/J3b_creeping_parameters/README.md): one loss probability creeps (the branch above *S. cerevisiae*+*N. crassa*+*A. nidulans* in the hard fit); totals and the origin ranking are stable from ~5,000 iterations; `--ml-tol 1e-5 --ml-param-tol 1e-6` halves the hard fit time with no visible cost (default unchanged, owner's decision).

## Settled (see SPEC; do not relitigate without a reason)

- Position identity is alignment column + phase (SPEC §2); missing is not absent (SPEC §2, "Missing data is not absence").
- Classic Dollo = one gain at the MRCA of carriers, minimise losses; weighted Dollo is the sensitivity tool (SPEC §5.1–5.2).
- Re-rooting splits an edge; MP is root-independent, Dollo placements are not (SPEC §9).
- Likelihood: error rates are inputs, not estimated; root branch fixed; no-data and unary branches fixed (SPEC §5.5).
- Third-party `Malin/` stays git-ignored; reference it by path, add no licensed fixtures to `tests/`.

## Verify

```sh
cargo build --release
cargo test --workspace --release            # one test is #[ignore]d: the real-data likelihood check (TESTING.md)
cargo clippy --workspace --all-targets --release && cargo fmt --check
python3 -I scripts/test_malin_to_matrix.py && python3 -I scripts/test_howto_scripts.py
target/release/intronfate demo -o intron_demo   # expect 3 sites, 3 gains, 2 losses, 3 MP steps
```

What each layer covers: [`TESTING.md`](TESTING.md). The HOWTO end-to-end run is on the HPCC ([D1](../bench/D1_howto_public_genomes/README.md)).
