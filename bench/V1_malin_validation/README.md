# V1 — validation against MALIN's published test data

Purpose: find out whether the figures quoted in `docs/SPEC.md` §6/§6b (which an earlier session produced
on a laptop and never re-derived) are right, and compare per-branch Dollo losses with the MALIN guide
(p. 27). **Purgeable:** everything this makes is under `results/`; delete it freely.

What it does (`validate.sh`):

1. Rogozin et al. 2003 table (`rogozin_introns_table`, 0/1 only): import with `convert`, score with
   `Rogozin-ecdysozoa.tre` (Dollo, MP, counts), cross-check against the independent Python importer.
2. `Chromo23.table` (phase codes + `#MAP`, 405 alignments in its `#MAP` block): import, score every alignment,
   export, and verify the round trip row by row (order-independent).
3. `summarize.py` prints observed vs. claimed values (it does **not** assert the claims) and writes
   `rogozin_branches.tsv` / `chromo23_branches.tsv` (per-branch gains and losses).

Run:

```sh
SMOKE=1 ./validate.sh                       # seconds, tiny fixtures; use to check the scripts anywhere
# on the HPCC (build + validate):
mkdir -p logs && MALIN_TEST=<dir containing Chromo23.table, rogozin_introns_table, ...> sbatch --export=ALL run.sbatch
```

Result of the first run (2026-11-11, HPCC `epyc`, 44 s): every spec figure reproduced, per-branch Dollo losses
matched the MALIN guide p. 27 exactly, round trip lossless — see `REPORT.2026-11-11.txt` and
`rogozin_branches.2026-11-11.tsv`. (Submit with `-p epyc` if `batch` is unavailable.)

Interpretation guide:

- **The spec's figures differ** → the spec was wrong or the implementation changed; correct the spec, not
  this script.
- **Per-branch losses vs guide p. 27** → MALIN's Dollo is the reference for classic Dollo. Differences on a
  branch point at rooting, the treatment of `?`, or a real bug. (Branch names in `dollo_branches.tsv` are
  tip names, tree labels, or `node<id>(<n>tips)`; map the guide's `Diptera` / `Ecdysozoa` to the listed tips.)
- This validates *scoring and interop* only. Projection from GFF (the tool's main job) has no MALIN
  counterpart and is covered by `examples/` and the unit tests.
