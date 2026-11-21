#!/usr/bin/env python3
"""S1: mean over seeds of results/S1_runs.tsv, one row per scenario, read model and fit (standard library only).
Run: summarize.py results/S1_runs.tsv > results/S1_summary.tsv"""
import collections, csv, statistics, sys

ORDER = ["none", "low", "mid", "high"]
MODELS = ["hard", "hard_err", "soft_k3", "soft_k0.2", "soft_p0.5_k0.2"]
COLS = ["true_multi", "true_noorigin", "auc_multi_lnlr", "called", "precision", "recall", "top100_precision", "auc_noorigin"]


def main(path):
    g = collections.defaultdict(list)
    for r in csv.DictReader(open(path), delimiter="\t"):
        g[(r["scenario"], r["reads"], r["model"])].append(r)
    print("\t".join(["scenario", "reads", "model", "seeds", "converged"] + COLS))
    for k in sorted(g, key=lambda k: (ORDER.index(k[0]), k[1], MODELS.index(k[2]))):
        rs = g[k]
        out = []
        for c in COLS:
            v = [float(r[c]) for r in rs if r[c] != "NA"]
            out.append(f"{statistics.mean(v):.3f}" if v and c not in ("true_multi", "true_noorigin", "called") else (f"{statistics.mean(v):.0f}" if v else "NA"))
        print("\t".join([k[0], k[1], k[2], str(len(rs)), f"{sum(r['converged'] == 'yes' for r in rs)}/{len(rs)}"] + out))


if __name__ == "__main__":
    main(sys.argv[1])
