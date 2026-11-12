#!/usr/bin/env python3
"""Score a matrix under every candidate placement tree and rank them.
For each trees/place_NN.nwk: `intronfate score -m dollo,mp`, then total Dollo gains and losses (dollo_branches.tsv) and Fitch steps (mp_sites.tsv).
Fewer steps / losses = better supported. Usage: run_scan.py --matrix M --trees trees --bin intronfate --work DIR --out placements.tsv [--threads 8]"""
import argparse, concurrent.futures as cf, csv, os, shutil, subprocess

ap = argparse.ArgumentParser()
for k in ("matrix", "trees", "bin", "work", "out"): ap.add_argument("--" + k, required=True)
ap.add_argument("--threads", type=int, default=8); a = ap.parse_args()
idx = {r["id"]: r for r in csv.DictReader(open(os.path.join(a.trees, "index.tsv")), delimiter="\t")}
shutil.rmtree(a.work, ignore_errors=True); os.makedirs(a.work)

def run(name, tree):
    d = os.path.join(a.work, name)
    subprocess.run([a.bin, "score", "--matrix", a.matrix, "-t", tree, "-m", "dollo,mp", "-o", d], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    g = l = 0
    for r in csv.DictReader(open(os.path.join(d, "dollo_branches.tsv")), delimiter="\t"): g += int(r["gains"]); l += int(r["losses"])
    steps = sum(int(r["steps"]) for r in csv.DictReader(open(os.path.join(d, "mp_sites.tsv")), delimiter="\t"))
    shutil.rmtree(d)
    return name, g, l, steps

jobs = [(k, os.path.join(a.trees, k + ".nwk")) for k in sorted(idx)] + [("gene_tree", os.path.join(a.trees, "gene_tree.nwk"))]
with cf.ThreadPoolExecutor(a.threads) as ex: res = list(ex.map(lambda j: run(*j), jobs))
rows = []
for name, g, l, s in res:
    i = idx.get(name, {"n_tips_below": "", "group": "(OrthoFinder gene tree)", "tips_below": ""})
    rows.append((name, i["group"], i["n_tips_below"], g, l, s, i["tips_below"]))
rows.sort(key=lambda r: (r[5], r[4]))
best = rows[0][5]
with open(a.out, "w") as f:
    f.write("tree\tgroup\tn_tips_below\tdollo_gains\tdollo_losses\tmp_steps\tsteps_over_best\ttips_below\n")
    for r in rows: f.write("\t".join(map(str, r[:6])) + f"\t{r[5]-best}\t{r[6]}\n")
print(open(a.out).read())
