#!/usr/bin/env python3
"""Calibration of the placement-by-Dollo-losses test with intron-poor stand-ins for Antarctolithica.

Question: if a lineage whose position we know had as few introns as Antarctolithica, would the same test still put it in the right
place, or would it be pulled somewhere else (and where)? For each control lineage X (one genome with a known sister group):
  1. trees: X removed together with Antarctolithica, then attached on every edge of the remaining backbone (V9 make_trees.py logic);
  2. matrices: X's intron presences thinned to Antarctolithica's level, a thinned intron becoming an observed absence:
       full     no thinning (the test on the real genome)
       uniform  every '1' kept with probability (Alith carried positions / X carried positions)
       age      positions not gained on X's own branch kept with probability (Alith's retention of positions gained above its parent
                node / X's own, the parent node being excluded because Dollo makes those carried by construction), and
                positions gained on X's own branch with probability (Alith's own-lineage positions / X's), so both the deep-intron
                retention and the lineage-specific count match Antarctolithica's
     each thinned variant REPS times with different seeds; every matrix also as a carrier-only version (rows where X carries it);
  3. every matrix is scored on every placement tree (`intronfate score -m dollo`, by run_calibration.sh), and the placement with the
     fewest Dollo losses is compared with the true one.
Thinning removes introns at random (within age class); real loss is not random (e.g. 3' or phase biases), so this is a first calibration.

prepare:   calibrate.py prepare --matrix M --gene-tree T --genomes G --gene-sites S --gene-branches B --controls a,b --reps 5 --out calib
summarize: calibrate.py summarize --calib calib --scores calib/scores.tsv --out results
"""
import argparse, collections, csv, os, statistics, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import v13lib as L  # noqa: E402


def prepare(a):
    root = L.parse_newick(open(a.gene_tree).read())
    tips = L.leaves(root)
    alith = sorted(t for t in tips if t.startswith(a.group_prefix))
    branches = L.read_branches(a.gene_branches)
    sites = L.read_dollo_sites(a.gene_sites, L.branch_index(branches))
    taxa, rows = L.read_matrix(a.matrix)
    a_prof = L.origin_profile(sites, rows, alith, root)
    old = lambda prof: [p for p in prof if not p["by_construction"]]      # origins above the parent node (see v13lib.origin_profile)
    a_old_p = sum(p["present"] for p in old(a_prof))
    a_old_a = sum(p["absent"] for p in old(a_prof))
    a_young = a_prof[-2]["present"] + a_prof[-1]["present"]  # gained on the Alith stem or inside
    a_carried = sum(1 for r in rows.values() if L.unit_state(r, alith) == "1")
    a_ret_old = a_old_p / (a_old_p + a_old_a)
    os.makedirs(a.out, exist_ok=True)
    jobs = open(os.path.join(a.out, "jobs.txt"), "w")
    plan = open(os.path.join(a.out, "plan.tsv"), "w")
    plan.write("control\tcarried\told_present\told_absent\tyoung\tp_uniform\tp_old\tp_young\talith_carried\talith_old_retention\talith_young\n")
    for x in a.controls.split(","):
        d = os.path.join(a.out, x)
        os.makedirs(d, exist_ok=True)
        subprocess.run([sys.executable, "-I", a.make_trees, a.gene_tree, "--group-prefix", x, "--drop-prefix", a.group_prefix,
                        "--genomes", a.genomes, "--out", os.path.join(d, "trees")], check=True, stdout=subprocess.DEVNULL)
        prof = L.origin_profile(sites, rows, [x], root)
        old_p = sum(p["present"] for p in old(prof)); old_a = sum(p["absent"] for p in old(prof))
        young = prof[-2]["present"] + prof[-1]["present"]
        carried = sum(1 for r in rows.values() if r[x] == "1")
        p_uni = min(1.0, a_carried / carried)
        p_old = min(1.0, a_ret_old / (old_p / (old_p + old_a)))
        p_young = min(1.0, a_young / young) if young else 1.0
        young_sites = {s for s, v in sites.items() if v["gained_at"] == frozenset([x])}
        plan.write(f"{x}\t{carried}\t{old_p}\t{old_a}\t{young}\t{p_uni:.4f}\t{p_old:.4f}\t{p_young:.4f}\t{a_carried}\t{a_ret_old:.4f}\t{a_young}\n")
        variants = [("full", 0, None)] + [("uniform", r, p_uni) for r in range(1, a.reps + 1)] + \
                   [("age", r, (lambda s, ys=young_sites, po=p_old, py=p_young: py if s in ys else po)) for r in range(1, a.reps + 1)]
        for name, rep, kp in variants:
            m = os.path.join(d, f"{name}_{rep}.matrix.tsv")
            if kp is None:
                if not os.path.exists(m):
                    os.symlink(os.path.abspath(a.matrix), m)
            else:
                seed = sum(map(ord, x)) * 1000 + rep * 10 + len(name)           # deterministic, different per control/variant/rep
                L.thin_matrix(a.matrix, m, x, kp, seed=seed)
            c = os.path.join(d, f"{name}_{rep}.carrier.matrix.tsv")
            L.carrier_only(m, c, x)
            for mode, path in (("all", m), ("carrier", c)):
                jobs.write(f"{x}\t{name}\t{rep}\t{mode}\t{path}\t{os.path.join(d, 'trees')}\n")
    jobs.close(); plan.close()
    print(open(os.path.join(a.out, "plan.tsv")).read(), end="")


def summarize(a):
    scores = collections.defaultdict(dict)                   # (x, variant, rep, mode) -> {place: losses}
    for r in csv.DictReader(open(a.scores), delimiter="\t"):
        scores[(r["control"], r["variant"], r["rep"], r["mode"])][r["tree"]] = int(r["dollo_losses"])
    rows = []
    for (x, var, rep, mode), v in sorted(scores.items()):
        d = os.path.join(a.calib, x, "trees")
        idx = {r["id"]: r for r in L.read_tsv(os.path.join(d, "index.tsv"))}
        true = open(os.path.join(d, "true_sibling.txt")).read().strip()
        backbone = L.parse_newick(open(os.path.join(d, "backbone.nwk")).read())
        lo = min(v.values())
        best = sorted(t for t in v if v[t] == lo)
        true_ids = [t for t in idx if idx[t]["tips_below"] == true]
        hit = "n/a" if not true_ids else ("yes" if true_ids[0] in best else "no")
        dist = min(L.edge_distance(backbone, frozenset(idx[b]["tips_below"].split(",")), frozenset(true.split(","))) for b in best) if true_ids else ""
        gap = (v[true_ids[0]] - lo) if true_ids else ""
        rows.append(dict(control=x, variant=var, rep=rep, mode=mode, true_sibling=true, best=",".join(best), n_tied=len(best), true_is_best=hit,
                         edges_from_true=dist, losses_true_minus_best=gap, best_sister_tips=idx[best[0]]["tips_below"]))
    os.makedirs(a.out, exist_ok=True)
    keys = list(rows[0])
    with open(os.path.join(a.out, "calibration_runs.tsv"), "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in rows:
            f.write("\t".join(str(r[k]) for k in keys) + "\n")
    agg = collections.defaultdict(list)
    for r in rows:
        agg[(r["control"], r["variant"], r["mode"])].append(r)
    with open(os.path.join(a.out, "calibration_summary.tsv"), "w") as f:
        f.write("control\tvariant\tmode\treplicates\ttrue_is_best\tmedian_edges_from_true\tmost_common_best_sister\n")
        for (x, var, mode), rs in sorted(agg.items()):
            hits = sum(r["true_is_best"] == "yes" for r in rs)
            ds = [r["edges_from_true"] for r in rs if r["edges_from_true"] != ""]
            common = collections.Counter(r["best_sister_tips"] for r in rs).most_common(1)[0][0]
            common = common if len(common) < 60 else common[:57] + "..."
            f.write(f"{x}\t{var}\t{mode}\t{len(rs)}\t{hits}\t{statistics.median(ds) if ds else ''}\t{common}\n")
    print(open(os.path.join(a.out, "calibration_summary.tsv")).read(), end="")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    for k in ("matrix", "gene_tree", "genomes", "gene_sites", "gene_branches", "controls", "out", "make_trees"):
        p.add_argument("--" + k.replace("_", "-"), required=True)
    p.add_argument("--reps", type=int, default=5); p.add_argument("--group-prefix", default="Alith")
    s = sub.add_parser("summarize")
    s.add_argument("--calib", required=True); s.add_argument("--scores", required=True); s.add_argument("--out", required=True)
    a = ap.parse_args()
    prepare(a) if a.cmd == "prepare" else summarize(a)
