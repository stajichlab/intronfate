#!/usr/bin/env python3
"""Compare candidate trees by their likelihood-model fits: fitted lnL per tree and a RELL bootstrap over families.

RELL (resampling estimated log-likelihoods; Kishino, Miyata & Hasegawa 1990): the per-position log-likelihoods of each tree's own
fit are resampled, by whole family (the `<family>:<column>:<phase>` label prefix), and the tree with the highest resampled sum wins
the replicate. Reported per tree: lnL, delta to the best, and the fraction of replicates it wins. Ties at 1e-9 split the win.
Usage: rell.py --score score_dir --ids gene_tree,place_01,... --reps 1000 --seed 1 > table.tsv
"""
import argparse, collections, csv, os, random, sys


def site_lnl(path):
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        return [(row["site"], float(row["log_lik"])) for row in r]


def family(label):
    parts = label.rsplit(":", 2)
    return parts[0] if len(parts) == 3 else label


def rell(per_tree, reps, seed):
    """per_tree: {tree: [(site, lnl)]} with identical site order. Returns {tree: (lnL, win fraction)}."""
    trees = sorted(per_tree)
    sites = [s for s, _ in per_tree[trees[0]]]
    for t in trees:
        if [s for s, _ in per_tree[t]] != sites:
            sys.exit(f"{t}: sites differ from {trees[0]}")
    fam = collections.OrderedDict()
    for i, s in enumerate(sites):
        fam.setdefault(family(s), []).append(i)
    sums = {t: [sum(per_tree[t][i][1] for i in idx) for idx in fam.values()] for t in trees}
    rng = random.Random(seed)
    n = len(fam)
    wins = collections.Counter()
    for _ in range(reps):
        draw = [rng.randrange(n) for _ in range(n)]
        tot = {t: sum(sums[t][k] for k in draw) for t in trees}
        hi = max(tot.values())
        best = [t for t in trees if tot[t] >= hi - 1e-9]
        for t in best:
            wins[t] += 1 / len(best)
    return {t: (sum(sums[t]), wins[t] / reps) for t in trees}, n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--score", required=True); ap.add_argument("--ids", required=True)
    ap.add_argument("--reps", type=int, default=1000); ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    per = {t: site_lnl(os.path.join(a.score, t, "likelihood_sites.tsv")) for t in a.ids.split(",")}
    res, nfam = rell(per, a.reps, a.seed)
    best = max(v[0] for v in res.values())
    print(f"# RELL over {nfam} families, {a.reps} replicates")
    print("tree\tlnL\tdelta\trell_win")
    for t, (l, w) in sorted(res.items(), key=lambda x: -x[1][0]):
        print(f"{t}\t{l:.3f}\t{best - l:.3f}\t{w:.4f}")


if __name__ == "__main__":
    main()
