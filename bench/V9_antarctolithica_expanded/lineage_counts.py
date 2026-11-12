#!/usr/bin/env python3
"""How many single-copy families does the set support when the near-identical Antarctolithica strains count as ONE lineage?
A family is single-copy when no genome has more than one gene in it. Lineages = every non-Alith genome + one for the Alith group (present
when any strain is). Usage: lineage_counts.py Orthogroups.tsv [HOG N1.tsv ...]"""
import csv, sys

def counts(path, first_col):
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        names = r.fieldnames[first_col:]
        alith = [n for n in names if n.startswith("Alith")]
        others = [n for n in names if not n.startswith("Alith")]
        tally = {}
        n_all = len(others) + (1 if alith else 0)
        for row in r:
            c = {n: len([g for g in row[n].replace(", ", ",").split(",") if g]) if row.get(n) else 0 for n in names}
            if any(v > 1 for v in c.values()): continue               # multi-copy somewhere
            L = sum(1 for n in others if c[n]) + (1 if any(c[n] for n in alith) else 0)
            has_alith = any(c[n] for n in alith)
            for k in (n_all, n_all - 2, n_all - 4, n_all - 8, n_all - 12):
                if L >= k: tally.setdefault(k, [0, 0]); tally[k][0] += 1; tally[k][1] += has_alith
        print(f"{path}: {n_all} lineages ({len(others)} genomes + Alith group)")
        for k in sorted(tally, reverse=True): print(f"   single-copy in >= {k} lineages: {tally[k][0]:5d} families, {tally[k][1]:5d} include an Antarctolithica strain")

counts(sys.argv[1], 1)
for p in sys.argv[2:]: counts(p, 3)

# Same count when a genome with >= 2 genes in the family is treated as MISSING there (its paralogs are ambiguous), instead of rejecting the family.
def pruned(path, first_col):
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        names = r.fieldnames[first_col:]
        alith = [n for n in names if n.startswith("Alith")]
        others = [n for n in names if not n.startswith("Alith")]
        tally = {}; n_all = len(others) + 1
        for row in r:
            c = {n: len([g for g in row[n].replace(", ", ",").split(",") if g]) if row.get(n) else 0 for n in names}
            single = {n for n, v in c.items() if v == 1}
            if not (single & set(alith)): continue                    # needs an Antarctolithica strain
            L = len(single & set(others)) + 1
            for k in (n_all, n_all - 4, n_all - 8, n_all - 12, n_all - 16):
                if L >= k: tally[k] = tally.get(k, 0) + 1
        print(f"{path}: multi-copy genomes treated as missing, Alith strain required")
        for k in sorted(tally, reverse=True): print(f"   >= {k} of {n_all} lineages: {tally[k]:5d} families")
pruned(sys.argv[1], 1)
for p in sys.argv[2:]: pruned(p, 3)
