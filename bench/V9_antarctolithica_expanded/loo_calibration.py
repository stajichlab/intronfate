#!/usr/bin/env python3
"""Summarise the calibration runs: for each known lineage treated as unknown, where did the scan put it, and is that its true position?
Reads loo/<X>/true_sibling.txt, loo/<X>/index.tsv and loo/<X>/<tier>_<mode>.tsv (from bootstrap_scan.py). A placement counts as true when it, or a tree tied
with it (identical per-family scores), attaches the lineage next to its true sibling clade; a lineage that is a child of the root has no single true edge (n/a).
Usage: loo_calibration.py loo TIER > table.tsv"""
import csv, os, sys
loo, tier = sys.argv[1:3]
print("target\tmode\ttrue_sibling\tbest_n_tips_below\tbest_tips_below\tbest_is_true\ttied_with_best\tbootstrap_best\tsteps_true_over_best\tdollo_best_n_tips_below\tdollo_best_is_true\tdollo_bootstrap_best")
for x in sorted(os.listdir(loo)):
    d = os.path.join(loo, x)
    if not os.path.isfile(os.path.join(d, "true_sibling.txt")): continue
    true = open(os.path.join(d, "true_sibling.txt")).read().strip()
    idx = {r["id"]: r for r in csv.DictReader(open(os.path.join(d, "index.tsv")), delimiter="\t")}
    for mode in ("all", "carrier"):
        p = os.path.join(d, f"{tier}_{mode}.tsv")
        if not os.path.exists(p): continue
        rows = [r for r in csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t") if r["tree"] != "gene_tree"]
        best = rows[0]
        group = [best["tree"]] + [t for t in best["tied_with"].split(",") if t and t in idx]
        is_true = "n/a" if true == "ROOT" else ("yes" if any(idx[t]["tips_below"] == true for t in group) else "NO")
        tr = [r for r in rows if r["tree"] in idx and idx[r["tree"]]["tips_below"] == true]
        dol = max(rows, key=lambda r: float(r["dollo_best_in_fraction_of_replicates"]))      # most often the best tree by Dollo losses
        dgroup = [dol["tree"]] + [t for t in dol["tied_with"].split(",") if t and t in idx]
        d_true = "n/a" if true == "ROOT" else ("yes" if any(idx[t]["tips_below"] == true for t in dgroup) else "NO")
        short = lambda t: t if len(t) < 40 else t[:37] + "..."
        print(f"{x}\t{mode}\t{short(true)}\t{best['n_tips_below']}\t{short(best['tips_below'])}\t{is_true}\t{len(group) - 1}\t{best['mp_best_in_fraction_of_replicates']}\t"
              f"{tr[0]['steps_over_best'] if tr else '?'}\t{dol['n_tips_below']}\t{d_true}\t{dol['dollo_best_in_fraction_of_replicates']}")
