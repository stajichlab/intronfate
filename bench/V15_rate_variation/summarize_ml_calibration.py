#!/usr/bin/env python3
"""Summarise the V15 likelihood calibration (V13's design, with the likelihood model instead of Dollo losses).

Input: calib/*.tsv from run_ml_calibration.sh (control, variant, rep, categories, tree, lnL, converged, shape) and V13's calibration
folder (each control's placement trees, index.tsv and true_sibling.txt). For every matrix and number of rate categories the
placement with the highest lnL is compared with the control's true position, as V13 did for the fewest Dollo losses:
  true_is_best             the true attachment has the highest lnL (ties within 1e-6 count)
  edges_from_true          edges between the best attachment and the true one
  depth_best_minus_true    negative: the best attachment is closer to the root than the true one (placed too deep)
  lnl_best_minus_true      how much better the best attachment fits than the true one
Writes calibration_runs.tsv and calibration_summary.tsv (per control, variant, categories) to --out and prints the summary.
Usage: summarize_ml_calibration.py --scores 'calib/*.tsv' --calib V13/calib --out results
"""
import argparse, collections, glob, os, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "V13_antarctolithica_intron_map"))
import v13lib as L  # noqa: E402

COLS = ["control", "variant", "rep", "categories", "tree", "lnl", "converged", "shape"]


def read_scores(pattern):
    runs = collections.defaultdict(dict)
    unconverged = 0
    for path in sorted(glob.glob(pattern)):
        if path.endswith(".part"):
            continue
        for line in open(path):
            r = dict(zip(COLS, line.rstrip("\n").split("\t")))
            if r["converged"] != "yes":
                unconverged += 1
            runs[(r["control"], r["variant"], r["rep"], r["categories"])][r["tree"]] = float(r["lnl"])
    return runs, unconverged


def judge(v, trees_dir):
    """v: {tree id: lnL}. Returns the comparison of the best placement with the true one."""
    idx = {r["id"]: r for r in L.read_tsv(os.path.join(trees_dir, "index.tsv"))}
    true = open(os.path.join(trees_dir, "true_sibling.txt")).read().strip()
    backbone = L.parse_newick(open(os.path.join(trees_dir, "backbone.nwk")).read())
    hi = max(v.values())
    best = sorted(t for t in v if v[t] >= hi - 1e-6)
    true_ids = [t for t in idx if idx[t]["tips_below"] == true and t in v]
    tips = lambda t: frozenset(idx[t]["tips_below"].split(","))
    out = dict(true_sibling=true, best=",".join(best), n_tied=len(best), best_sister_tips=idx[best[0]]["tips_below"])
    if not true_ids:
        return dict(out, true_is_best="n/a", edges_from_true="", depth_best_minus_true="", lnl_best_minus_true="")
    t = true_ids[0]
    return dict(out, true_is_best="yes" if t in best else "no",
                edges_from_true=min(L.edge_distance(backbone, tips(b), tips(t)) for b in best),
                depth_best_minus_true=min(L.depth_of(backbone, tips(b)) for b in best) - L.depth_of(backbone, tips(t)),
                lnl_best_minus_true=round(hi - v[t], 3))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", required=True); ap.add_argument("--calib", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    runs, unconverged = read_scores(a.scores)
    rows = []
    for (x, var, rep, k), v in sorted(runs.items()):
        rows.append(dict(control=x, variant=var, rep=rep, categories=k, n_trees=len(v), **judge(v, os.path.join(a.calib, x, "trees"))))
    os.makedirs(a.out, exist_ok=True)
    keys = ["control", "variant", "rep", "categories", "n_trees", "true_sibling", "best", "n_tied", "true_is_best", "edges_from_true",
            "depth_best_minus_true", "lnl_best_minus_true", "best_sister_tips"]
    with open(os.path.join(a.out, "calibration_runs.tsv"), "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in rows:
            f.write("\t".join(str(r[k]) for k in keys) + "\n")
    agg = collections.defaultdict(list)
    for r in rows:
        agg[(r["control"], r["variant"], r["categories"])].append(r)
    lines = ["control\tvariant\tcategories\treplicates\ttrue_is_best\tmedian_edges_from_true\tmedian_depth_best_minus_true\tmost_common_best_sister"]
    for (x, var, k), rs in sorted(agg.items()):
        hits = sum(r["true_is_best"] == "yes" for r in rs)
        ds = [r["edges_from_true"] for r in rs if r["edges_from_true"] != ""]
        dd = [r["depth_best_minus_true"] for r in rs if r["depth_best_minus_true"] != ""]
        common = collections.Counter(r["best_sister_tips"] for r in rs).most_common(1)[0][0]
        common = common if len(common) < 60 else common[:57] + "..."
        lines.append(f"{x}\t{var}\t{k}\t{len(rs)}\t{hits}\t{statistics.median(ds) if ds else ''}\t"
                     f"{statistics.median(dd) if dd else ''}\t{common}")
    with open(os.path.join(a.out, "calibration_summary.tsv"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    tot = collections.defaultdict(lambda: [0, 0, 0, 0])   # (variant, k) -> runs, hits, deeper, shallower
    for r in rows:
        t = tot[(r["variant"], r["categories"])]
        t[0] += 1; t[1] += r["true_is_best"] == "yes"
        if r["depth_best_minus_true"] != "":
            t[2] += r["depth_best_minus_true"] < 0; t[3] += r["depth_best_minus_true"] > 0
    print("\n# variant\tcategories\truns\ttrue_is_best\tplaced_deeper\tplaced_shallower")
    for (var, k), (n, h, d, s) in sorted(tot.items()):
        print(f"# {var}\t{k}\t{n}\t{h}\t{d}\t{s}")
    if unconverged:
        print(f"# warning: {unconverged} fits did not converge", file=sys.stderr)


if __name__ == "__main__":
    main()
