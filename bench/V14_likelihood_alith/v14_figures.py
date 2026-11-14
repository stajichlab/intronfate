#!/usr/bin/env python3
"""V14 tables and figures from the likelihood runs (run_v14.sbatch). Uses V13's tree drawing (v13lib, v13_figures).

  results/branches_measured.tsv   per branch: Dollo gains/losses beside the likelihood's expected gains/losses/present with
                                  bootstrap intervals (error set "measured")
  results/alith_sweep.tsv         Antarctolithica's stem and ancestor under each error set
  results/placements.tsv          per placement and error set: fitted lnL, delta to the best, RELL win fraction
  figures/fig6_ml_tree.png        gene tree: expected gains/losses per branch, expected positions present at each node [95 % CI]
  figures/fig7_dollo_vs_ml.png    per-branch losses and gains, Dollo against the likelihood model
  figures/fig8_ml_placements.png  backbone coloured by the lnL cost of attaching Antarctolithica on each branch
Usage: v14_figures.py --trees V9/trees --genomes V9/genomes.tsv --score score --place place --results results --out figures
"""
import argparse, csv, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "V13_antarctolithica_intron_map"))
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import cm, colors  # noqa: E402
import v13lib as L  # noqa: E402
from v13_figures import SC_COLOR, draw_tree, layout, read_groups  # noqa: E402


def tsv(path):
    return L.read_tsv(path)


def by_tips(rows):
    return {frozenset(r["tips"].split(",")): r for r in rows}


def branch_table(a):
    ml = by_tips(tsv(os.path.join(a.score, "gene_measured", "likelihood_branches.tsv")))
    do = by_tips(tsv(os.path.join(a.score, "gene_measured", "dollo_branches.tsv")))
    out = os.path.join(a.results, "branches_measured.tsv")
    cols = ["branch", "n_tips", "fixed", "dollo_gains", "dollo_losses", "ml_gains", "ml_gains_lo", "ml_gains_hi", "ml_losses",
            "ml_losses_lo", "ml_losses_hi", "ml_present", "ml_present_lo", "ml_present_hi", "tips"]
    with open(out, "w") as f:
        f.write("\t".join(cols) + "\n")
        for ts, m in sorted(ml.items(), key=lambda kv: -len(kv[0])):
            d = do.get(ts, {})
            f.write("\t".join([m["branch"], m["n_tips"], m["fixed"], d.get("gains", ""), d.get("losses", ""), m["expected_gains"],
                               m.get("gains_lo", ""), m.get("gains_hi", ""), m["expected_losses"], m.get("losses_lo", ""),
                               m.get("losses_hi", ""), m["expected_present"], m.get("present_lo", ""), m.get("present_hi", ""),
                               m["tips"]]) + "\n")
    return ml, do


def alith_sweep(a, alith):
    rows = []
    for d in sorted(os.listdir(a.score)):
        p = os.path.join(a.score, d, "likelihood_branches.tsv")
        if not d.startswith("gene_") or not os.path.exists(p):
            continue
        b = by_tips(tsv(p))
        stem = b[frozenset(alith)]
        summ = {r["key"]: r["value"] for r in tsv(os.path.join(a.score, d, "summary_score.tsv"))}
        tips = {r["taxon"]: r for r in tsv(os.path.join(a.score, d, "likelihood_tips.tsv"))}
        missed = sum(float(tips[t]["expected_missed"]) for t in alith)
        rows.append(dict(error_set=d[5:], log_lik=summ["ml_log_likelihood"], stem_gains=stem["expected_gains"],
                         stem_losses=stem["expected_losses"], ancestor_present=stem["expected_present"],
                         alith_expected_missed=f"{missed:.1f}", total_losses=summ["ml_expected_losses"],
                         total_gains=summ["ml_expected_gains"], unobserved=summ["ml_unobserved_positions"]))
    keys = list(rows[0])
    with open(os.path.join(a.results, "alith_sweep.tsv"), "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in rows:
            f.write("\t".join(str(r[k]) for k in keys) + "\n")
    return rows


def placements(a):
    idx = {r["id"]: r for r in tsv(os.path.join(a.trees, "index.tsv"))}
    out = []
    for fn in sorted(os.listdir(a.results)):
        if not (fn.startswith("rell_") and fn.endswith(".tsv")):
            continue
        es = fn[5:-4]
        for r in tsv(os.path.join(a.results, fn)):
            r["error_set"] = es
            r["sister_of"] = idx[r["tree"]]["tips_below"] if r["tree"] in idx else "(gene tree)"
            out.append(r)
    with open(os.path.join(a.results, "placements.tsv"), "w") as f:
        f.write("error_set\ttree\tlnL\tdelta\trell_win\tsister_of\n")
        for r in out:
            f.write(f"{r['error_set']}\t{r['tree']}\t{r['lnL']}\t{r['delta']}\t{r['rell_win']}\t{r['sister_of']}\n")
    return out, idx


def fig_tree(a, sc, ml):
    root = L.parse_newick(open(os.path.join(a.trees, "gene_tree.nwk")).read())
    obs = {r["tip"]: r for r in tsv(os.path.join(a.v13results, "tip_counts.tsv"))}
    xs, ys, tips, H = layout(root)
    fig, ax = plt.subplots(figsize=(15, 13))
    red = lambda x: SC_COLOR["Antarctolithica"] if all(t.startswith("Alith") for t in L.leaves(x)) else "black"
    draw_tree(ax, root, xs, ys, edge_color=red, edge_width=lambda x: 2.2 if red(x) != "black" else 1.2)
    parent_x = {}
    for p in L.nodes(root):
        for k in p["kids"]:
            parent_x[id(k)] = xs[id(p)]
    for x in L.nodes(root):
        r = ml[L.tipset(x)]
        g, l = float(r["expected_gains"]), float(r["expected_losses"])
        if x is not root:
            mid = (parent_x[id(x)] + xs[id(x)]) / 2
            lab = "fixed" if r["fixed"] in ("root", "unary") else f"+{g:.0f} −{l:.0f}"
            ax.text(mid, ys[id(x)] - 0.18, lab, ha="center", va="bottom", fontsize=6.3, color="#444444")
        if x["kids"]:
            ci = f" [{float(r['present_lo']):.0f}-{float(r['present_hi']):.0f}]" if r.get("present_lo") else ""
            ax.scatter([xs[id(x)]], [ys[id(x)]], s=22, color="white", edgecolor="black", zorder=3)
            ax.text(xs[id(x)] - 0.05, ys[id(x)] + 0.12, f"{float(r['expected_present']):.0f}{ci}", fontsize=6.0, color="#1a5276",
                    va="top", ha="right", zorder=4)
    for t in tips:
        ax.text(H + 0.15, ys[id(t)], f"{t['name']}  {obs[t['name']]['carried']}", va="center", fontsize=8,
                color=SC_COLOR.get(sc.get(t["name"], ""), "black"))
    ax.set_xlim(-1.2, H + 3.2); ax.set_ylim(-1, len(tips)); ax.invert_yaxis(); ax.axis("off")
    handles = [plt.Line2D([], [], color=c, lw=6, label=g) for g, c in SC_COLOR.items()]
    ax.legend(handles=handles, loc="lower left", fontsize=8, frameon=False)
    ax.set_title("Likelihood model with annotation error (measured error rates where known): expected gains (+) and losses (−) per branch;\n"
                 "left of each node the expected positions present [95 % family bootstrap]; tips: positions observed. 'fixed': branch held at "
                 "no change (root rule)", fontsize=10, loc="left")
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig6_ml_tree.png"), dpi=170); plt.close(fig)


def fig_compare(a, ml, do, sc):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    for ax, key_d, key_m, lo, hi, title in ((axes[0], "losses", "expected_losses", "losses_lo", "losses_hi", "losses"),
                                             (axes[1], "gains", "expected_gains", "gains_lo", "gains_hi", "gains")):
        for ts, m in ml.items():
            if ts not in do or m["fixed"] != "no" or len(ts) == 42:
                continue
            x, y = float(do[ts][key_d]), float(m[key_m])
            is_a = all(t.startswith("Alith") for t in ts)
            col = SC_COLOR["Antarctolithica"] if is_a else ("#1f77b4" if len(ts) > 1 else "#7f7f7f")
            if m.get(lo):
                ax.plot([x, x], [float(m[lo]), float(m[hi])], color=col, lw=0.8, alpha=0.6)
            ax.scatter(x, y, s=40 if is_a else 16, color=col, zorder=3)
            if is_a and len(ts) == 6:
                ax.annotate("Antarctolithica stem", (x, y), xytext=(5, -10), textcoords="offset points", fontsize=8)
        mx = max(ax.get_xlim()[1], ax.get_ylim()[1])
        ax.plot([0, mx], [0, mx], color="grey", lw=0.6, ls="--")
        ax.set_xlabel(f"Dollo {title}"); ax.set_ylabel(f"likelihood: expected {title} [95 % bootstrap]")
        ax.set_title(f"{title} per branch (grey: tips, blue: internal)", fontsize=10, loc="left")
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig7_dollo_vs_ml.png"), dpi=170); plt.close(fig)


def fig_place(a, sc, pl, idx):
    sets = sorted({r["error_set"] for r in pl}, key=lambda s: ["none", "measured"].index(s) if s in ("none", "measured") else 9)
    backbone = L.parse_newick(open(os.path.join(a.trees, "backbone.nwk")).read())
    xs, ys, tips, H = layout(backbone)
    fig, axes = plt.subplots(1, len(sets), figsize=(8.5 * len(sets), 10))
    axes = [axes] if len(sets) == 1 else axes
    parent_x = {}
    for p in L.nodes(backbone):
        for k in p["kids"]:
            parent_x[id(k)] = xs[id(p)]
    for ax, es in zip(axes, sets):
        rows = {r["tree"]: r for r in pl if r["error_set"] == es and r["tree"] in idx}
        cost = {frozenset(idx[t]["tips_below"].split(",")): float(rows[t]["delta"]) for t in rows}
        norm = colors.PowerNorm(0.5, vmin=0, vmax=max(cost.values()))
        cmap = matplotlib.colormaps["viridis_r"]
        draw_tree(ax, backbone, xs, ys, edge_color=lambda x: cmap(norm(cost[L.tipset(x)])) if L.tipset(x) in cost else "lightgrey",
                  edge_width=lambda x: 4.0)
        for x in L.nodes(backbone):
            ts = L.tipset(x)
            if ts in cost and x is not backbone:
                win = float(rows[next(t for t in rows if frozenset(idx[t]["tips_below"].split(",")) == ts)]["rell_win"])
                lab = f"{cost[ts]:.0f}" + (f" ({win:.2f})" if win >= 0.01 else "")
                ax.text((parent_x[id(x)] + xs[id(x)]) / 2, ys[id(x)] - 0.2, lab, ha="center", va="bottom", fontsize=6.3)
        for t in tips:
            ax.text(H + 0.15, ys[id(t)], t["name"], va="center", fontsize=8, color=SC_COLOR.get(sc.get(t["name"], ""), "black"))
        ax.set_xlim(-0.5, H + 2.5); ax.set_ylim(-1, len(tips)); ax.invert_yaxis(); ax.axis("off")
        fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.5, label="lnL below the best placement")
        ax.set_title(f"error set '{es}': Δ lnL per attachment branch (RELL win fraction)", fontsize=10, loc="left")
    fig.tight_layout(); fig.savefig(os.path.join(a.out, "fig8_ml_placements.png"), dpi=150); plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("trees", "genomes", "score", "results", "out", "v13results"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sc = read_groups(a.genomes)
    alith = sorted(g for g in sc if g.startswith("Alith"))
    ml, do = branch_table(a)
    alith_sweep(a, alith)
    pl, idx = placements(a)
    fig_tree(a, sc, ml)
    fig_compare(a, ml, do, sc)
    fig_place(a, sc, pl, idx)
    print("V14 tables in", a.results, "figures in", a.out)
