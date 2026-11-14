# Testing

```sh
cargo test --workspace                      # Rust: unit, integration, CLI end-to-end, doc tests
cargo clippy --workspace --all-targets      # must stay warning-free
cargo fmt --check
python3 -I scripts/test_malin_to_matrix.py  # the Python converter's fixture tests
python3 -I scripts/test_howto_scripts.py    # the HOWTO helper scripts' fixture tests
```

## Layers

| Layer | Where | What it proves |
|---|---|---|
| Unit | `crates/intronfate-core/src/*.rs` (`#[cfg(test)]`) | Coordinate model, GFF parsing, CDS merging on both strands, 5'-partial CDS, frame/phase checks, Newick, bitsets, MALIN parsing. |
| Brute-force | `score.rs` tests | Dollo DP and Fitch/Hartigan are compared with exhaustive search over every `{present, absent, missing}` tip pattern on eight tree shapes, polytomies included. Reconstructions must realise exactly the reported step count. |
| Pipeline integration | `crates/intronfate-core/tests/per_cell_missing.rs` | Join → project → matrix → file → score on tiny hand-built families: a gap scores `?` and no loss, the same taxon aligned through scores one loss, `?` survives `matrix.tsv`, adding an all-missing taxon changes no score, trimmed alignment rows are rejected. |
| Combined matrix | `crates/intronfate-core/tests/combined.rs` | A combined matrix scores exactly the sum of its families (totals and every branch), keeps the same column in two families as two sites, never pairs sites across families, keeps missing taxa `?`, survives the file round trip, and refuses duplicate names, mixed flank/phase mode and already-combined inputs. |
| CLI end to end | `crates/intronfate-cli/tests/cli.rs` | The real binary: golden output for every example, hand-calculated values, `matrix`+`score` = `analyze`, gzip input, exit codes, error messages, output-folder collisions, and the list and folder input forms (`--msa-dir`/`--msa-list` give byte-identical results to repeated `--msa`; `--gff-dir`/`--proteins-list`; `score --matrix-dir`, `convert --matrix-dir/--matrix-list`, a `combined` sub-folder skipped). `inputs.rs` unit tests cover the extension rule, list parsing and folder scanning. |
| Examples | `examples/*/` | Each has `run.sh` and a committed `expected/`; the CLI test re-runs them and requires byte-identical output. Their READMEs list the hand-checked numbers. |
| Randomised vs reference | `matrix.rs`, `annotation.rs`, `tree.rs`, `score.rs` tests | Random CDS layouts on both strands (with 5'-partial skips and exon rows) against a naive per-base model; random matrices through `matrix.tsv` and MALIN; random trees for `mrca`, Dollo and Fitch against brute force; weighted Dollo against brute force with missing tips; the parallel pairwise code against the per-pair reference; the `ambiguous` flag against an exhaustive count of optimal labellings. |
| Cross-implementation | `scripts/test_malin_to_matrix.py` vs `malin.rs` tests | Two independent MALIN importers agree on the same fixture (`tests/data/malin_tiny.table`). |
| Likelihood | `likelihood.rs` tests, `cli.rs::likelihood_writes_its_tables_and_reads_error_rates` | Likelihood, ascertainment term, node and pair posteriors equal brute-force enumeration (errors, `?`, polytomy, unscored tip); conditional probabilities of all seen patterns sum to 1; EM never decreases the likelihood; simulation recovery of losses, and of terminal losses only when misses are modelled; the larger root branch stays fixed; family-block bootstrap; CLI tables, error-file forms and refusals. |
| How-to helpers | `scripts/test_howto_scripts.py` | `one_protein_per_gene.py` keeps the longest isoform per gene for NCBI- and Ensembl-style GFF3 (ties, CDS under a gene, internal stops, missing proteins); `make_families.py` builds the tiers from a tiny OrthoFinder layout with plain and gzip proteomes; `doc_blocks.py` extracts only ```` ```sh ```` blocks and refuses a `--sub` line it cannot find. |
| How-to, end to end | `bench/D1_howto_public_genomes/run_howto.sbatch` (HPCC) | `docs/HOWTO.md` run as written, on both examples, from download to `dollo_branches.tsv`. |

Golden files only catch *changes*. The numbers that make the examples trustworthy are the
hand-calculated ones (`quickstart_numbers_match_the_hand_calculation` and the example READMEs).

## Test-suite review (2026-11-11)

An Opus agent applied 44 single-line mutations one at a time to a copy of the repo and ran the suite on each.
Surviving mutants became tests: a hand-written truth table for `cell_observable` (every phase × flank × gap
position, plus the merged `--ignore-phase` window), weighted Dollo against brute force, the clamp boundary,
`read_tsv`/MALIN duplicate-row merging, `write_atomic` leaving no temp file, counts with unannotated taxa, and
thread-count independence of the output. It also found two real bugs: the `ambiguous` flag in `mp_sites.tsv`
was true for every unrooted-ish case, and an empty matrix exported to an unreadable MALIN table. Both are fixed.
Golden files are only regression tests (they were generated by the code); the independent evidence is the
hand-calculated values and the randomised reference tests above.

## Adding a case

- New parsing behaviour: a unit test next to the code with the smallest failing input.
- New science (a scoring rule, observability rule): a hand-built family in
  `per_cell_missing.rs`, and if it is a DP, extend the brute-force check.
- New example: copy an `examples/` folder, keep `run.sh` self-contained (`INTRONFATE` env var,
  relative output dir `out`), generate `expected/` with `./run.sh out && rm -rf expected && mv out expected`,
  **check the numbers by hand**, add the name to the list in `every_example_reproduces_its_expected_output_exactly`.

## What is not tested here

Real-genome behaviour and speed. Per project policy analyses and benchmarks run on the HPCC; a
benchmark plan and the numbers to reproduce are in [`SPEC.md`](SPEC.md) §6 and [`HANDOFF.md`](HANDOFF.md).
The third-party MALIN datasets in `Malin/` are git-ignored and are not used by the automated suite.
