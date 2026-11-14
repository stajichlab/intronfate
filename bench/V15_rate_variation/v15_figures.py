#!/usr/bin/env python3
"""V15 tables and figures: the likelihood model with rate classes (run_v15.sbatch) and the likelihood placement calibration
(calib_ml.sbatch). Uses V13's tree drawing and V14's helpers.

  results/gene_fits.tsv            per error set and number of classes: lnL, shape, totals, Alith stem and ancestor, its parent
  results/branches_measured.tsv    per branch: Dollo, likelihood with 1 class and with K classes (bootstrap interval, measured errors)
  results/placements.tsv           per placement: lnL, delta, RELL win, for each error set with K classes
  results/calibration_compare.tsv  per control, variant and method (Dollo V13, likelihood 1 class, K classes): runs, right, deeper,
                                   shallower
  figures/fig9_ml_classes_tree.png gene tree with K classes: expected gains/losses per branch, positions present [95 % CI]
  figures/fig10_classes_effect.png per branch, 1 class against K classes: losses, and positions present at each node
  figures/fig11_ml_classes_placements.png  backbone coloured by the lnL cost of attaching Alith on each branch (K classes)
  figures/fig12_calibration.png    where intron-poor stand-ins of known position are placed: Dollo against the likelihood
Usage: v15_figures.py --trees V9/trees --genomes V9/genomes.tsv --score score --results results --out figures \
           --v13results V13/results --k 4
"""
import argparse, collections, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "V13_antarctolithica_intron_map"))
sys.path.insert(0, os.path.join(HERE, "..", "V14_likelihood_alith"))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import v13lib as L  # noqa: E402
from v13_figures import SC_COLOR, draw_tree, layout, read_groups  # noqa: E402
import v14_figures as V14  # noqa: E402

tsv = L.read_tsv
by_tips = V14.by_tips


def write(path, rows, keys=None):
    keys = keys or list(rows[0])
    with open(path, "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")


def gene_fits(a, alith):
    rows = []
    for d in sorted(os.listdir(a.score)):
        p = os.path.join(a.score, d, "likelihood_branches.tsv")
        if not d.startswith("gene_") or not os.path.exists(p):
            continue
        es, k = d[5:].rsplit("_k", 1)
        b = by_tips(tsv(p))
        stem = b[frozenset(alith)]
        parent = min((r for ts, r in b.items() if frozenset(alith) < ts), key=lambda r: int(r["n_tips"]))
        s = {r["key"]: r["value"] for r in tsv(os.path.join(a.score, d, "summary_score.tsv"))}
        rows.append(dict(error_set=es, categories=k, log_lik=s["ml_log_likelihood"], gamma_shape=s.get("ml_gamma_shape", ""),
                         iterations=s["ml_iterations"], total_gains=s["ml_expected_gains"], total_losses=s["ml_expected_losses"],
                         unobserved=s["ml_unobserved_positions"], root_present=s["ml_root_expected_present"],
                         stem_gains=stem["expected_gains"], stem_losses=stem["expected_losses"],
                         stem_losses_ci=f"{stem.get('losses_lo', '')}-{stem.get('losses_hi', '')}" if stem.get("losses_lo") else "",
                         ancestor_present=stem["expected_present"],
                         ancestor_ci=f"{stem.get('present_lo', '')}-{stem.get('present_hi', '')}" if stem.get("present_lo") else "",
                         parent_present=parent["expected_present"]))
    write(os.path.join(a.results, "gene_fits.tsv"), rows)
    return rows


def branch_table(a):
    k1 = by_tips(tsv(os.path.join(a.score, "gene_measured_k1", "likelihood_branches.tsv")))
    kk = by_tips(tsv(os.path.join(a.score, f"gene_measured_k{a.k}", "likelihood_branches.tsv")))
    do = by_tips(tsv(os.path.join(a.score, f"gene_measured_k{a.k}", "dollo_branches.tsv")))
    rows = []
    for ts, m in sorted(kk.items(), key=lambda kv: -len(kv[0])):
        d, o = do.get(ts, {}), k1[ts]
        rows.append(dict(branch=m["branch"], n_tips=m["n_tips"], fixed=m["fixed"], dollo_gains=d.get("gains", ""),
                         dollo_losses=d.get("losses", ""), k1_gains=o["expected_gains"], k1_losses=o["expected_losses"],
                         k1_present=o["expected_present"], gains=m["expected_gains"], gains_lo=m.get("gains_lo", ""),
                         gains_hi=m.get("gains_hi", ""), losses=m["expected_losses"], losses_lo=m.get("losses_lo", ""),
                         losses_hi=m.get("losses_hi", ""), present=m["expected_present"], present_lo=m.get("present_lo", ""),
                         present_hi=m.get("present_hi", ""), tips=m["tips"]))
    write(os.path.join(a.results, "branches_measured.tsv"), rows)
    return k1, kk


def placements(a):
    idx = {r["id"]: r for r in tsv(os.path.join(a.trees, "index.tsv"))}
    out = []
    for fn in sorted(os.listdir(a.results)):
        if fn.startswith("rell_") and fn.endswith(".tsv"):
            for r in tsv(os.path.join(a.results, fn)):
                r["error_set"] = fn[5:-4]
                r["sister_of"] = idx[r["tree"]]["tips_below"] if r["tree"] in idx else "(gene tree)"
                out.append(r)
    write(os.path.join(a.results, "placements.tsv"), out, ["error_set", "tree", "lnL", "delta", "rell_win", "sister_of"])
    return out, idx


def calibration(a):
    """Dollo (V13, every position) beside the likelihood runs. Returns {(method, variant): {control: Counter}}."""
    agg = collections.defaultdict(lambda: collections.defaultdict(collections.Counter))
    def add(method, r):
        c = agg[(method, r["variant"])][r["control"]]
        c["runs"] += 1
        if r["true_is_best"] == "yes":
            c["right"] += 1
        elif r["depth_best_minus_true"] not in ("", None):
            d = int(float(r["depth_best_minus_true"]))
            c["deeper" if d < 0 else "shallower" if d > 0 else "level"] += 1
    for r in tsv(os.path.join(a.v13results, "calibration_runs.tsv")):
        if r["mode"] == "all":
            add("Dollo", r)
    for k in ("1", str(a.k)):
        p = os.path.join(a.results, f"calib_k{k}", "calibration_runs.tsv")
        if os.path.exists(p):
            for r in tsv(p):
                add(f"likelihood, {k} class{'es' if k != '1' else ''}", r)
    rows = [dict(method=m, variant=v, control=x, **{k: c[k] for k in ("runs", "right", "deeper", "shallower", "level")})
            for (m, v), per in sorted(agg.items()) for x, c in sorted(per.items())]
    write(os.path.join(a.results, "calibration_compare.tsv"), rows)
    return agg


def fig_tree(a, sc, ml):
    root = L.parse_newick(open(os.path.join(a.trees, "gene_tree.nwk")).read())
    obs = {r["tip"]: r for r in tsv(os.path.join(a.v13results, "tip_counts.tsv"))}
    xs, ys, tips, H = layout(root)
    fig, ax = plt.subplots(figsize=(15, 13))
    red = lambda x: SC_COLOR["Antarctolithica"] if all(t.startswith("Alith") for t in L.leaves(x)) else "black"
    draw_tree(ax, root, xs, ys, edge_color=red, edge_width=lambda x: 2.2 if red(x) != "black" else 1.2)
    parent_x = {id(k): xs[id(p)] for p in L.nodes(root) for k in p["kids"]}
    for x in L.nodes(root):
        r = ml[L.tipset(x)]
        if x is not root:
            lab = "fixed" if r["fixed"] in ("root", "unary") else f"+{float(r['expected_gains']):.0f} −{float(r['expected_losses']):.0f}"
            ax.text((parent_x[id(x)] + xs[id(x)]) / 2, ys[id(x)] - 0.18, lab, ha="center", va="bottom", fontsize=6.3, color="#444444")
        if x["kids"]:
            ci = f" [{float(r['present_lo']):.0f}-{float(r['present_hi']):.0f}]" if r.get("present_lo") else ""
            ax.scatter([xs[id(x)]], [ys[id(x)]], s=22, color="white", edgecolor="black", zorder=3)
            ax.text(xs[id(x)] - 0.05, ys[id(x)] + 0.12, f"{float(r['expected_present']):.0f}{ci}", fontsize=6.0, color="#1a5276",
                    va="top", ha="right", zorder=4)
    for t in tips:
        ax.text(H + 0.15, ys[id(t)], f"{t['name']}  {obs[t['name']]['carried']}", va="center", fontsize=8,
                color=SC_COLOR.get(sc.get(t["name"], ""), "black"))
    ax.set_xlim(-1.2, H + 3.2); ax.set_ylim(-1, len(tips)); ax.invert_yaxis(); ax.axis("off")
    ax.legend(handles=[plt.Line2D([], [], color=c, lw=6, label=g) for g, c in SC_COLOR.items()], loc="lower left", fontsize=8,
              frameon=False)
    ax.set_title(f"Likelihood model with {a.k} rate classes and annotation error (measured rates where known): expected gains (+) and "
                 "losses (−) per branch;\nleft of each node the expected positions present [95 % family bootstrap]; tips: positions "
                 "observed. 'fixed': branch held at no change (root rule)", fontsize=10, loc="left")
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig9_ml_classes_tree.png"), dpi=170); plt.close(fig)


def fig_effect(a, k1, kk):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    for ax, key, title in ((axes[0], "expected_losses", "expected losses on the branch"),
                           (axes[1], "expected_present", "expected positions present at the node")):
        for ts, m in kk.items():
            if m["fixed"] != "no" or (key == "expected_present" and len(ts) == 1):
                continue
            x, y = float(k1[ts][key]), float(m[key])
            is_a = all(t.startswith("Alith") for t in ts)
            col = SC_COLOR["Antarctolithica"] if is_a else ("#1f77b4" if len(ts) > 1 else "#7f7f7f")
            lo, hi = ("losses_lo", "losses_hi") if key == "expected_losses" else ("present_lo", "present_hi")
            if m.get(lo):
                ax.plot([x, x], [float(m[lo]), float(m[hi])], color=col, lw=0.8, alpha=0.6)
            ax.scatter(x, y, s=40 if is_a else 16, color=col, zorder=3)
            if is_a and len(ts) == 6:
                ax.annotate("Antarctolithica", (x, y), xytext=(5, -10), textcoords="offset points", fontsize=8)
        mx = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([0, mx], [0, mx], color="grey", lw=0.6, ls="--")
        ax.set_xlabel(f"1 class: {title}"); ax.set_ylabel(f"{a.k} classes [95 % bootstrap]")
        ax.set_title(f"{title} (grey: tips, blue: internal)", fontsize=10, loc="left")
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig10_classes_effect.png"), dpi=170); plt.close(fig)


def fig_calibration(a, agg):
    variants = ["full", "uniform", "age"]
    methods = ["Dollo", "likelihood, 1 class", f"likelihood, {a.k} classes"]
    controls = sorted({x for per in agg.values() for x in per})
    fig, axes = plt.subplots(len(variants), len(methods), figsize=(15, 11), sharey=True)
    cols = {"right": "#2e8b57", "deeper": "#e67e22", "shallower": "#8e44ad", "level": "#95a5a6"}
    for i, v in enumerate(variants):
        for j, m in enumerate(methods):
            ax = axes[i][j]
            per = agg.get((m, v), {})
            if not per:
                ax.text(0.5, 0.5, "not run", ha="center", va="center", transform=ax.transAxes, color="grey")
            left = [0.0] * len(controls)
            for kind in ("right", "level", "shallower", "deeper"):
                w = [per[x][kind] / per[x]["runs"] if x in per and per[x]["runs"] else 0 for x in controls]
                ax.barh(range(len(controls)), w, left=left, color=cols[kind], label=kind if (i, j) == (0, 0) else None)
                left = [l + x for l, x in zip(left, w)]
            tot = collections.Counter()
            for c in per.values():
                tot.update(c)
            ax.set_title(f"{m} — {v}" + (f": {tot['right']}/{tot['runs']} right, {tot['deeper']} deeper, "
                                         f"{tot['shallower']} shallower" if tot["runs"] else ""), fontsize=8.5, loc="left")
            ax.set_yticks(range(len(controls)), controls, fontsize=7.5); ax.set_xlim(0, 1)
            if i == len(variants) - 1:
                ax.set_xlabel("fraction of runs")
    axes[0][0].invert_yaxis()
    fig.legend(loc="lower center", ncol=4, frameon=False, fontsize=9,
               title="best placement relative to the true one (deeper = closer to the root)")
    fig.suptitle("Calibration with intron-poor stand-ins: full genome (no thinning), uniform thinning to Alith's intron count, "
                 "age-matched thinning (V13 design)", fontsize=10.5, x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0.05, 1, 0.97)); fig.savefig(os.path.join(a.out, "fig12_calibration.png"), dpi=150); plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("trees", "genomes", "score", "results", "out", "v13results"):
        ap.add_argument("--" + k, required=True)
    ap.add_argument("--k", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sc = read_groups(a.genomes)
    alith = sorted(g for g in sc if g.startswith("Alith"))
    gene_fits(a, alith)
    k1, kk = branch_table(a)
    pl, idx = placements(a)
    fig_tree(a, sc, kk)
    fig_effect(a, k1, kk)
    a2 = argparse.Namespace(**vars(a))
    V14.fig_place(a2, sc, pl, idx)
    os.replace(os.path.join(a.out, "fig8_ml_placements.png"), os.path.join(a.out, "fig11_ml_classes_placements.png"))
    fig_calibration(a, calibration(a))
    print("V15 tables in", a.results, "figures in", a.out)
