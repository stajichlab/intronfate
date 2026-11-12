#!/usr/bin/env python3
"""Placement scan with a bootstrap over gene families.

Each candidate tree (trees/place_NN.nwk plus the gene tree) is scored with `intronfate score -m dollo,mp`; per-site Fitch steps and Dollo losses are
summed per family, and families are resampled with replacement to see how often each placement is the best. `--require-carrier PREFIX` keeps only
sites where some taxon whose name starts with PREFIX carries the intron, so that absences in that lineage (shared loss) cannot attract it to other
intron-poor lineages. Usage: bootstrap_scan.py --matrix M --trees DIR --bin intronfate --work DIR --out tsv [--require-carrier Alith] [--reps 1000]"""
import argparse, collections, concurrent.futures as cf, csv, os, shutil, subprocess
import numpy as np

ap = argparse.ArgumentParser()
for k in ("matrix", "trees", "bin", "work", "out"): ap.add_argument("--" + k, required=True)
ap.add_argument("--threads", type=int, default=8); ap.add_argument("--reps", type=int, default=1000); ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--require-carrier", default=""); a = ap.parse_args()
shutil.rmtree(a.work, ignore_errors=True); os.makedirs(a.work)
matrix = a.matrix
if a.require_carrier:
    lines = open(a.matrix).read().split("\n"); out = []; hdr = None
    for l in lines:
        if not l or l.startswith("#"): out.append(l); continue
        if hdr is None: hdr = l.split("\t"); out.append(l); first = hdr.index("pattern") + 1; idx = [i for i in range(first, len(hdr)) if hdr[i].startswith(a.require_carrier)]; continue
        c = l.split("\t")
        if any(c[i] == "1" for i in idx): out.append(l)
    matrix = os.path.join(a.work, "filtered.tsv"); open(matrix, "w").write("\n".join(out) + "\n")
index = {r["id"]: r for r in csv.DictReader(open(os.path.join(a.trees, "index.tsv")), delimiter="\t")}
names = sorted(index) + ["gene_tree"]

def run(name):
    tree = os.path.join(a.trees, name + ".nwk"); d = os.path.join(a.work, name)
    subprocess.run([a.bin, "score", "--matrix", matrix, "-t", tree, "-m", "dollo,mp", "-o", d], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    mp = {r["site"]: int(r["steps"]) for r in csv.DictReader(open(os.path.join(d, "mp_sites.tsv")), delimiter="\t")}
    lo = {r["site"]: int(r["losses"]) for r in csv.DictReader(open(os.path.join(d, "dollo_sites.tsv")), delimiter="\t")}
    shutil.rmtree(d); return name, mp, lo

with cf.ThreadPoolExecutor(a.threads) as ex: res = list(ex.map(run, names))
sites = sorted(res[0][1]); fam = sorted({s.split(":")[0] for s in sites}); fi = {f: i for i, f in enumerate(fam)}
def table(k):
    X = np.zeros((len(fam), len(names)), dtype=np.int64)
    for t, r in enumerate(res):
        for s, v in r[k].items(): X[fi[s.split(":")[0]], t] += v
    return X
MP, LO = table(1), table(2)
rng = np.random.default_rng(a.seed); F = len(fam)
W = rng.multinomial(F, [1.0 / F] * F, size=a.reps)
def wins(X):
    tot = W @ X; best = tot == tot.min(axis=1, keepdims=True)
    return (best / best.sum(axis=1, keepdims=True)).sum(axis=0) / a.reps
wm, wl = wins(MP), wins(LO)
order = np.argsort(MP.sum(axis=0))
with open(a.out, "w") as f:
    f.write(f"# {len(sites)} sites in {F} families; {a.reps} family-bootstrap replicates; carrier filter: {a.require_carrier or 'none'}\n")
    f.write("tree\tgroup\ttips_below\tmp_steps\tsteps_over_best\tdollo_losses\tmp_best_in_fraction_of_replicates\tdollo_best_in_fraction_of_replicates\n")
    for t in order:
        n = names[t]; i = index.get(n, {"group": "(gene tree)", "tips_below": ""})
        f.write(f"{n}\t{i['group'][:60]}\t{i['tips_below'][:60]}\t{MP[:, t].sum()}\t{MP[:, t].sum() - MP.sum(axis=0).min()}\t{LO[:, t].sum()}\t{wm[t]:.3f}\t{wl[t]:.3f}\n")
print(open(a.out).read())
