#!/usr/bin/env python3
"""Re-score a tier's matrices after simulating annotation error (see perturb.py), many replicates, and summarise per-branch Dollo gains/losses.

Usage: run_sensitivity.py --proj <results_projection/TIER> --tree species.nwk --accuracy <busco_accuracy.tsv> --bin <intronfate>
                          --work $SCRATCH/v7 --out results/TIER [--reps 20] [--seed 1] [--threads 8]
Replicate 0 of mode `observed` is the unperturbed matrix (it must equal the tier's own Dollo table). Modes: observed; ceiling_fn / typical_fn (missed introns only, filled at terminal-loss cells); fp (false introns only); ceiling (both).
Writes <out>/replicates.tsv (mode rep branch gains losses) and <out>/summary.tsv (mean and 2.5/97.5 percentiles per mode and branch).
"""
import argparse, collections, concurrent.futures as cf, csv, glob, os, random, shutil, statistics, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import perturb as P

ap = argparse.ArgumentParser()
ap.add_argument("--proj", required=True); ap.add_argument("--tree", required=True); ap.add_argument("--accuracy", required=True)
ap.add_argument("--bin", required=True); ap.add_argument("--work", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--reps", type=int, default=20); ap.add_argument("--seed", type=int, default=1); ap.add_argument("--threads", type=int, default=8)
a = ap.parse_args()

sn, pr = {}, {}
for r in csv.DictReader(open(a.accuracy), delimiter="\t"):
    if r.get("sn_v2_genes"): sn[r["genome"]] = float(r["sn_v2_genes"]) / 100; pr[r["genome"]] = float(r["pr_v2_genes"]) / 100
alns = sorted(glob.glob(os.path.join(a.proj, "*.aln", "matrix.tsv")))
os.makedirs(a.out, exist_ok=True); shutil.rmtree(a.work, ignore_errors=True); os.makedirs(a.work)
print(f"{len(alns)} matrices; Sn/Pr for {len(sn)} genomes: " + " ".join(f"{k}={sn[k]:.3f}/{pr[k]:.3f}" for k in sorted(sn)), flush=True)
loaded = {m: P.read_matrix(m) for m in alns}
tree = P.parse_newick(open(a.tree).read())

def score_one(args):
    mode, rep, i, m = args
    meta, hdr, rows = loaded[m]
    d = os.path.join(a.work, f"{mode}_{rep}", str(i)); os.makedirs(d)
    if mode != "observed":
        rows = P.perturb(hdr, rows, sn, pr, mode, random.Random(f"{a.seed}/{mode}/{rep}/{i}"), tree)
    if not rows: return []
    P.write_matrix(os.path.join(d, "matrix.tsv"), meta, hdr, rows)
    subprocess.run([a.bin, "score", "--matrix", os.path.join(d, "matrix.tsv"), "-t", a.tree, "-m", "dollo", "-o", os.path.join(d, "out")],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    res = [(r["tips"], int(r["gains"]), int(r["losses"])) for r in csv.DictReader(open(os.path.join(d, "out", "dollo_branches.tsv")), delimiter="\t")]
    shutil.rmtree(d)
    return res

jobs = [("observed", 0, i, m) for i, m in enumerate(alns)]
MODES = ("ceiling_fn", "typical_fn", "fp", "ceiling")
for mode in MODES:
    jobs += [(mode, rep, i, m) for rep in range(1, a.reps + 1) for i, m in enumerate(alns)]
tot = collections.defaultdict(lambda: [0, 0]); done = 0
with cf.ThreadPoolExecutor(a.threads) as ex:
    for job, res in zip(jobs, ex.map(score_one, jobs, chunksize=16)):
        for tips, g, l in res:
            tot[(job[0], job[1], tips)][0] += g; tot[(job[0], job[1], tips)][1] += l
        done += 1
        if done % 5000 == 0: print(f"  {done}/{len(jobs)} scored", flush=True)
with open(os.path.join(a.out, "replicates.tsv"), "w") as f:
    f.write("mode\trep\tbranch\tgains\tlosses\n")
    for (mode, rep, tips), (g, l) in sorted(tot.items()): f.write(f"{mode}\t{rep}\t{tips}\t{g}\t{l}\n")

def pct(xs, q):
    xs = sorted(xs); return xs[min(len(xs) - 1, int(q * len(xs)))]
branches = sorted({t for (_, _, t) in tot}, key=lambda t: (-t.count(","), t))
with open(os.path.join(a.out, "summary.tsv"), "w") as f:
    f.write("branch\tn_tips\tobs_gains\tobs_losses\t" + "\t".join(f"{m}_{k}_{w}" for m in MODES for k in ("gains", "losses") for w in ("mean", "p2.5", "p97.5")) + "\n")
    for t in branches:
        row = [t, str(t.count(",") + 1), str(tot[("observed", 0, t)][0]), str(tot[("observed", 0, t)][1])]
        for m in MODES:
            for k in (0, 1):
                xs = [tot[(m, r, t)][k] for r in range(1, a.reps + 1)]
                row += [f"{statistics.mean(xs):.0f}", str(pct(xs, 0.025)), str(pct(xs, 0.975))]
        f.write("\t".join(row) + "\n")
print(open(os.path.join(a.out, "summary.tsv")).read())
