#!/usr/bin/env python3
"""Does the `OG` column of Phylogenetic_Hierarchical_Orthogroups/N*.tsv name the same orthogroup as the row of that id in Orthogroups.tsv?
For every HOG row, its genes must all lie in the Orthogroups.tsv row with that OG id. Usage: check_hog_join.py results/plain"""
import csv, sys, glob, os
res = sys.argv[1]
og = {}
with open(os.path.join(res, "Orthogroups", "Orthogroups.tsv"), newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        og[r["Orthogroup"]] = {g for k, v in r.items() if k != "Orthogroup" for g in v.replace(", ", ",").split(",") if g}
for p in sorted(glob.glob(os.path.join(res, "Phylogenetic_Hierarchical_Orthogroups", "N*.tsv"))):
    bad = n = 0
    with open(p, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            genes = {g for k, v in r.items() if k not in ("HOG", "OG", "Gene Tree Parent Clade") for g in v.replace(", ", ",").split(",") if g}
            n += 1
            if not genes <= og.get(r["OG"], set()): bad += 1
    print(f"{os.path.basename(p)}: {n} HOG rows, {bad} whose genes are NOT all in Orthogroups.tsv row {'OG'}")
