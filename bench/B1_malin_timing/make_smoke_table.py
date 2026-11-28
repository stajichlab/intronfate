#!/usr/bin/env python3
"""Write a small deterministic MALIN intron table and tree for the B1 smoke test (no third-party data).

Six taxa on ((A,B),(C,D),(E,F)); each site is one gain at a random node and independent losses on the branches below it,
so Dollo and the likelihood both have something to find. The table has `#MAP` lines (alignment name, alignment column),
so the phase is kept and `convert --to-malin` can write it back. Usage: make_smoke_table.py OUT_PREFIX [N_SITES] [SEED]
"""
import random
import sys

TAXA = list("ABCDEF")
CLADES = [list("AB"), list("CD"), list("EF"), TAXA]  # nodes a gain can occur at
TREE = "((A,B),(C,D),(E,F));\n"


def build(n_sites=60, seed=1):
    rng = random.Random(seed)
    rows = {t: [] for t in TAXA}
    for _ in range(n_sites):
        carriers = rng.choice(CLADES)
        for t in TAXA:
            present = t in carriers and rng.random() > 0.25  # 25 % loss on the way to each tip
            rows[t].append("1" if present else "0")
    return rows


def write(prefix, n_sites=60, seed=1):
    rows = build(n_sites, seed)
    with open(prefix + ".table", "w") as fh:
        for t in TAXA:
            fh.write(f"#ORGANISM {t} {t}\n")
        for t in TAXA:
            fh.write(f"{t}\t{''.join(rows[t])}\n")
        fh.write("#MAP\tsite\tkog\talipos\n")
        for i in range(n_sites):
            fh.write(f"#MAP\t{i}\tK0\t{3 + 7 * i}\n")
    with open(prefix + ".nwk", "w") as fh:
        fh.write(TREE)


if __name__ == "__main__":
    write(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 60, int(sys.argv[3]) if len(sys.argv) > 3 else 1)
