# Examples

Three small, runnable, hand-checkable examples (01–03), plus the inputs of the public-genome how-to (04). Each of 01–03 has:

- the input files,
- `run.sh` — runs the tool (set `INTRONFATE` to the binary, or have it on `PATH`; output goes to `out/`),
- `expected/` — the committed output, which `cargo test` re-generates and compares byte for byte,
- a `README.md` with the numbers worked out by hand.

| Example | Shows |
|---|---|
| [`01_quickstart`](01_quickstart/README.md) | The full pipeline on four species, one on the minus strand, with two shared sites that each need a loss. |
| [`02_missing_data`](02_missing_data/README.md) | An alignment gap is missing data (`?`), not an intron loss; non-multiple-of-3 exons and phases. |
| [`03_malin`](03_malin/README.md) | Importing a MALIN table, scoring it, and exporting it back losslessly. |
| [`04_public_genomes`](04_public_genomes/README.md) | Genome lists and species trees for [`docs/HOWTO.md`](../docs/HOWTO.md) (downloads complete genomes; not part of `cargo test`). |

Build the binary first ([Install](../README.md#install)), then from the repository root:

```sh
INTRONFATE=$PWD/target/release/intronfate examples/01_quickstart/run.sh
```
