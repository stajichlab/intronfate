#!/usr/bin/env python3
"""Re-score one tier's COMBINED matrix after simulating annotation error (see perturb.py), many replicates, and summarise per-branch Dollo gains/losses.

Usage: run_sensitivity.py --matrix <tier>/combined/matrix.tsv --tree species.nwk --accuracy <busco_accuracy.tsv> --bin intronfate
                          --work $SCRATCH/v7 --out results/TIER [--reps 100] [--seed 1] [--threads 8] [--expect <tier>/combined/dollo_branches.tsv]
Replicate 0 (mode `observed`) is the unperturbed matrix; with --expect it must equal that Dollo table exactly (asserted). Candidate cells for the missed introns
come from the observed run's dollo_sites.tsv. Modes: observed; ceiling_fn / typical_fn (missed introns only); fp (false introns only); ceiling (both).
Writes <out>/replicates.tsv (mode rep branch gains losses) and <out>/summary.tsv (mean and the 2.5th / 97.5th percentiles over replicates per mode and branch)."""
import argparse, collections, concurrent.futures as cf, csv, os, random, shutil, subprocess, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import perturb as P

ap = argparse.ArgumentParser()
for k in ("matrix", "tree", "accuracy", "bin", "work", "out"): ap.add_argument("--" + k, required=True)
ap.add_argument("--reps", type=int, default=100); ap.add_argument("--seed", type=int, default=1); ap.add_argument("--threads", type=int, default=8)
ap.add_argument("--expect", default=""); a = ap.parse_args()

sn, pr = {}, {}
for r in csv.DictReader(open(a.accuracy), delimiter="\t"):
    if r.get("sn_v2_genes") and r["sn_v2_genes"] != "nan":
        sn[r["genome"]] = float(r["sn_v2_genes"]) / 100; pr[r["genome"]] = float(r["pr_v2_genes"]) / 100
os.makedirs(a.out, exist_ok=True); shutil.rmtree(a.work, ignore_errors=True); os.makedirs(a.work)
meta, hdr, rows = P.read_matrix(a.matrix)
taxa = hdr[hdr.index("pattern") + 1:]
print(f"{len(rows)} sites x {len(taxa)} taxa; Sn/Pr for {len(sn)} genomes", flush=True)

def score(name, mrows):
    d = os.path.join(a.work, name); os.makedirs(d)
    P.write_matrix(os.path.join(d, "matrix.tsv"), meta, hdr, mrows)
    subprocess.run([a.bin, "score", "--matrix", os.path.join(d, "matrix.tsv"), "-t", a.tree, "-m", "dollo", "-o", os.path.join(d, "out")],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return os.path.join(d, "out")

# replicate 0: the observed matrix, also the source of the candidate cells
out0 = score("observed", rows)
observed = {r["tips"]: (int(r["gains"]), int(r["losses"])) for r in csv.DictReader(open(os.path.join(out0, "dollo_branches.tsv")), delimiter="\t")}
if a.expect:
    want = {r["tips"]: (int(r["gains"]), int(r["losses"])) for r in csv.DictReader(open(a.expect), delimiter="\t")}
    assert observed == want, "replicate 0 differs from the tier's own Dollo table"
site_index = {r[0]: i for i, r in enumerate(rows)}
cand = P.candidates_from_dollo(os.path.join(out0, "dollo_sites.tsv"), site_index, taxa)
print("terminal-loss candidate cells per taxon: " + " ".join(f"{t}={len(v)}" for t, v in cand.items() if v), flush=True)

def run(job):
    mode, rep = job
    mrows = P.perturb(hdr, rows, sn, pr, mode, random.Random(f"{a.seed}/{mode}/{rep}"), cand)
    out = score(f"{mode}_{rep}", mrows)
    res = [(r["tips"], int(r["gains"]), int(r["losses"])) for r in csv.DictReader(open(os.path.join(out, "dollo_branches.tsv")), delimiter="\t")]
    shutil.rmtree(os.path.join(a.work, f"{mode}_{rep}"))
    return mode, rep, res

MODES = ("ceiling_fn", "typical_fn", "fp", "ceiling")
tot = collections.defaultdict(lambda: [0, 0])
for tips, (g, l) in observed.items(): tot[("observed", 0, tips)] = [g, l]
with cf.ThreadPoolExecutor(a.threads) as ex:
    for mode, rep, res in ex.map(run, [(m, r) for m in MODES for r in range(1, a.reps + 1)]):
        for tips, g, l in res: tot[(mode, rep, tips)] = [g, l]
with open(os.path.join(a.out, "replicates.tsv"), "w") as f:
    f.write("mode\trep\tbranch\tgains\tlosses\n")
    for (mode, rep, tips), (g, l) in sorted(tot.items()): f.write(f"{mode}\t{rep}\t{tips}\t{g}\t{l}\n")
branches = sorted(observed, key=lambda t: (-t.count(","), t))
with open(os.path.join(a.out, "summary.tsv"), "w") as f:
    f.write("branch\tn_tips\tobs_gains\tobs_losses\t" + "\t".join(f"{m}_{k}_{w}" for m in MODES for k in ("gains", "losses") for w in ("mean", "p2.5", "p97.5")) + "\n")
    for t in branches:
        row = [t, str(t.count(",") + 1), str(observed[t][0]), str(observed[t][1])]
        for m in MODES:
            for k in (0, 1):
                xs = np.array([tot[(m, r, t)][k] for r in range(1, a.reps + 1)])
                row += [f"{xs.mean():.0f}", f"{np.percentile(xs, 2.5):.0f}", f"{np.percentile(xs, 97.5):.0f}"]
        f.write("\t".join(row) + "\n")
print(open(os.path.join(a.out, "summary.tsv")).read())
