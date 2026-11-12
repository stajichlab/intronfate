#!/usr/bin/env python3
"""Summarise the calibration runs: for each known lineage treated as unknown, where did the scan put it, and is that its true position?
Reads loo/<X>/true_sibling.txt and loo/<X>/scan_<mode>.tsv (from bootstrap_scan.py). Usage: loo_calibration.py loo TIER > table.tsv"""
import csv, os, sys
loo, tier = sys.argv[1:3]
print("target\tmode\ttrue_sibling\tbest_placement_tips_below\tbest_is_true\tbootstrap_best\tsteps_true_over_best\tdollo_best_tips")
for x in sorted(os.listdir(loo)):
    d = os.path.join(loo, x)
    if not os.path.isfile(os.path.join(d, "true_sibling.txt")): continue
    true = open(os.path.join(d, "true_sibling.txt")).read().strip()
    for mode in ("all", "carrier"):
        p = os.path.join(d, f"{tier}_{mode}.tsv")
        if not os.path.exists(p): continue
        rows = [r for r in csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")]
        rows = [r for r in rows if r["tree"] != "gene_tree"]
        best = rows[0]
        idx = {r["id"]: r["tips_below"] for r in csv.DictReader(open(os.path.join(d, "index.tsv")), delimiter="\t")}
        full = lambda r: idx[r["tree"]]
        tr = [r for r in rows if full(r) == true]
        dol = min(rows, key=lambda r: int(r["dollo_losses"]))
        short = lambda t: t if len(t) < 40 else t[:37] + "..."
        print(f"{x}\t{mode}\t{short(true)}\t{short(full(best))}\t{'yes' if full(best)==true else 'NO'}\t{best['mp_best_in_fraction_of_replicates']}\t{tr[0]['steps_over_best'] if tr else '?'}\t{short(full(dol))}")
