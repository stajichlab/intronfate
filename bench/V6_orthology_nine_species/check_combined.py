#!/usr/bin/env python3
"""Check that `analyze --combine` (results_projection/<tier>/combined/) equals the per-family Dollo gains/losses summed afterwards.
Usage: check_combined.py results_projection     exits 1 on any difference."""
import collections, csv, glob, os, sys
root = sys.argv[1]; bad = 0
def read(path):
    return {r["tips"]: (int(r["gains"]), int(r["losses"])) for r in csv.DictReader(open(path), delimiter="\t")}
for tier in sorted(d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d))):
    summed = collections.defaultdict(lambda: [0, 0]); fams = 0
    for aln in glob.glob(os.path.join(root, tier, "*.aln")):
        fams += 1
        for k, (g, l) in read(os.path.join(aln, "dollo_branches.tsv")).items():
            summed[k][0] += g; summed[k][1] += l
    comb = read(os.path.join(root, tier, "combined", "dollo_branches.tsv"))
    diff = {k: (summed.get(k), comb.get(k)) for k in set(summed) | set(comb) if tuple(summed.get(k, (0, 0))) != comb.get(k, (0, 0))}
    n_sites = sum(1 for l in open(os.path.join(root, tier, "combined", "matrix.tsv")) if not l.startswith("#")) - 1
    print(f"{tier}: {fams} families, combined matrix {n_sites} sites, branches {len(comb)}, differences from the summed per-family tables: {len(diff)}")
    for k, v in diff.items(): print("  DIFF", k, v)
    bad += len(diff)
sys.exit(1 if bad else 0)
