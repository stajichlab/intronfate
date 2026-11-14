#!/usr/bin/env python3
"""V13 figures (matplotlib; run with `pixi run python v13_figures.py ...`, see run_figures.sh).

  fig1_tree_gains_losses.png  gene tree; on every branch Dollo gains (+) and losses (-); at every node the estimated number of intron
                              positions present (Dollo, with Fitch MP in brackets); at every tip the positions the genome carries
  fig2_placement_cost.png     backbone without Antarctolithica; each branch coloured by the extra Dollo losses of attaching it there
                              (left: all positions; right: positions Antarctolithica carries)
  fig3_retention.png          retention of the oldest (pan-fungal) intron positions against intron richness, every lineage; and
                              Antarctolithica's carried/lost positions by the node where they arose
  fig4_calibration.png        intron-poor stand-ins: distance (branches) between the cheapest and the true placement
  fig5_dollo_weighted.png     total gains and losses as a second origin of a position gets cheaper (gain:loss cost sweep)
Usage: v13_figures.py --trees V9/trees --genomes V9/genomes.tsv --score score --results results --out figures
"""
import argparse, csv, math, os, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import cm, colors  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v13lib as L  # noqa: E402
from analyze import superclade  # noqa: E402

SC_COLOR = {"Mucoromycota": "#7f7f7f", "Basidiomycota": "#8c564b", "Taphrinomycotina": "#2ca02c", "Saccharomycotina": "#1f77b4",
            "Pezizomycotina": "#9467bd", "Antarctolithica": "#d62728"}


def layout(root, order_key=None):
    """Cladogram coordinates: tips at x = max depth (aligned), internal nodes one step left of their deepest child."""
    ys, xs = {}, {}
    tips = []

    def ladder(x):
        if x["kids"]:
            x["kids"].sort(key=lambda k: len(L.leaves(k)))
            for k in x["kids"]:
                ladder(k)
    ladder(root)

    def height(x):
        return 0 if not x["kids"] else 1 + max(height(k) for k in x["kids"])
    H = height(root)

    def place(x, depth):
        if not x["kids"]:
            tips.append(x)
            ys[id(x)] = len(tips) - 1
            xs[id(x)] = H
            return
        for k in x["kids"]:
            place(k, depth + 1)
        ys[id(x)] = sum(ys[id(k)] for k in x["kids"]) / len(x["kids"])
        xs[id(x)] = min(xs[id(k)] for k in x["kids"]) - 1
    place(root, 0)
    return xs, ys, tips, H


def draw_tree(ax, root, xs, ys, edge_color=lambda x: "black", edge_width=lambda x: 1.2):
    for x in L.nodes(root):
        for k in x["kids"]:
            ax.plot([xs[id(x)], xs[id(x)]], [ys[id(x)], ys[id(k)]], color="black", lw=0.8, zorder=1)
            ax.plot([xs[id(x)], xs[id(k)]], [ys[id(k)], ys[id(k)]], color=edge_color(k), lw=edge_width(k), zorder=2,
                    solid_capstyle="butt")


def read_groups(genomes):
    sc = {}
    for l in open(genomes):
        if l.strip() and not l.startswith("#"):
            f = l.rstrip("\n").split("\t")
            sc[f[0]] = superclade(f[2])
    return sc


def fig1(a, sc):
    root = L.parse_newick(open(os.path.join(a.trees, "gene_tree.nwk")).read())
    anc = {frozenset(r["tips"].split(",")): r for r in L.read_tsv(os.path.join(a.results, "ancestral_states.tsv"))}
    xs, ys, tips, H = layout(root)
    fig, ax = plt.subplots(figsize=(15, 13))
    red = lambda x: SC_COLOR["Antarctolithica"] if all(t.startswith("Alith") for t in L.leaves(x)) else "black"
    draw_tree(ax, root, xs, ys, edge_color=red, edge_width=lambda x: 2.2 if red(x) != "black" else 1.2)
    for x in L.nodes(root):
        r = anc[L.tipset(x)]
        if x is not root:
            par_x = next(xs[id(p)] for p in L.nodes(root) if any(k is x for k in p["kids"]))
            mid = (par_x + xs[id(x)]) / 2
            ax.text(mid, ys[id(x)] - 0.18, f"+{r['branch_gains']} −{r['branch_losses']}", ha="center", va="bottom", fontsize=6.3,
                    color="#444444")
        else:
            ax.text(xs[id(x)] - 0.15, ys[id(x)], f"root: +{r['branch_gains']}", ha="right", va="center", fontsize=7)
        if x["kids"]:
            ax.scatter([xs[id(x)]], [ys[id(x)]], s=22, color="white", edgecolor="black", zorder=3)
            ax.text(xs[id(x)] + 0.06, ys[id(x)] + 0.22, f"{r['dollo_present']} ({r['mp_present']})", fontsize=6.3, color="#1a5276",
                    va="top", zorder=4)
    for t in tips:
        r = anc[frozenset([t["name"]])]
        ax.text(H + 0.15, ys[id(t)], f"{t['name']}  {r['dollo_present']}", va="center", fontsize=8,
                color=SC_COLOR.get(sc.get(t["name"], ""), "black"))
    ax.set_xlim(-1.2, H + 3.2)
    ax.set_ylim(-1, len(tips))
    ax.invert_yaxis()
    ax.axis("off")
    handles = [plt.Line2D([], [], color=c, lw=6, label=g) for g, c in SC_COLOR.items()]
    ax.legend(handles=handles, loc="lower left", fontsize=8, frameon=False)
    ax.set_title("Dollo gains (+) and losses (−) of intron positions per branch; at each node the positions estimated present, "
                 "Dollo (Fitch MP);\nat each tip the positions carried. V9 tier A: 1,704 families, 29,085 positions, gene tree rooted on Mucor",
                 fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "fig1_tree_gains_losses.png"), dpi=170)
    plt.close(fig)


def fig2(a, sc):
    backbone = L.parse_newick(open(os.path.join(a.trees, "backbone.nwk")).read())
    idx = {r["id"]: frozenset(r["tips_below"].split(",")) for r in L.read_tsv(os.path.join(a.trees, "index.tsv"))}
    sup = {r["tree"]: r for r in L.read_tsv(os.path.join(a.results, "placement_support.tsv"))}
    xs, ys, tips, H = layout(backbone)
    fig, axes = plt.subplots(1, 2, figsize=(18, 10))
    for ax, key, title in ((axes[0], "delta_all", "all 29,085 positions"), (axes[1], "delta_carrier", "positions Antarctolithica carries")):
        by_tips = {idx[t]: int(sup[t][key]) for t in sup}
        vmax = max(by_tips.values())
        norm = colors.PowerNorm(0.5, vmin=0, vmax=vmax)
        cmap = matplotlib.colormaps["viridis_r"]
        col = lambda x: cmap(norm(by_tips[L.tipset(x)])) if L.tipset(x) in by_tips else "lightgrey"
        draw_tree(ax, backbone, xs, ys, edge_color=col, edge_width=lambda x: 4.0)
        for x in L.nodes(backbone):
            ts = L.tipset(x)
            if ts in by_tips:
                par_x = next((xs[id(p)] for p in L.nodes(backbone) if any(k is x for k in p["kids"])), xs[id(x)] - 1)
                ax.text((par_x + xs[id(x)]) / 2, ys[id(x)] - 0.2, str(by_tips[ts]), ha="center", va="bottom", fontsize=6.5)
        gene = next(t for t in sup if sup[t]["gene_tree_position"] == "yes")
        best = [t for t in sup if int(sup[t][key]) == 0]
        for t, mk, lab in [(gene, "*", "gene-tree position"), (best[0], "v", "fewest losses")]:
            node = next(x for x in L.nodes(backbone) if L.tipset(x) == idx[t])
            par_x = next((xs[id(p)] for p in L.nodes(backbone) if any(k is node for k in p["kids"])), xs[id(node)] - 1)
            ax.scatter([(par_x + xs[id(node)]) / 2], [ys[id(node)] + 0.25], marker=mk, s=120, color="red", zorder=5, label=lab)
        for t in tips:
            ax.text(H + 0.15, ys[id(t)], t["name"], va="center", fontsize=8, color=SC_COLOR.get(sc.get(t["name"], ""), "black"))
        ax.set_xlim(-0.5, H + 2.5)
        ax.set_ylim(-1, len(tips))
        ax.invert_yaxis()
        ax.axis("off")
        ax.legend(loc="lower left", fontsize=8, frameon=False)
        fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.5, label="extra Dollo losses vs the best placement")
        ax.set_title(f"Antarctolithica attached on each branch: {title}", fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "fig2_placement_cost.png"), dpi=170)
    plt.close(fig)


def fig3(a, sc):
    rows = L.read_tsv(os.path.join(a.results, "retention_by_origin.tsv"))
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={"width_ratios": [1, 1.25]})
    for depth, mk in (("0", "o"), ("2", "s")):
        pts = [r for r in rows if r["depth"] == depth and r["retention"] != "nan" and r["by_construction"] == "no"]
        for r in pts:
            is_a = r["lineage"] == "Antarctolithica"
            ax1.scatter(int(r["carried"]), float(r["retention"]), marker=mk, s=110 if is_a else 38,
                        color=SC_COLOR.get(r["superclade"], "black"), edgecolor="black" if is_a else "none", zorder=3 if is_a else 2,
                        alpha=1 if depth == "0" else 0.55)
            if is_a or r["lineage"] in ("Yarrowia", "Spombe", "Neolecta", "Saitoella", "Taphrina", "Lipomyces"):
                ax1.annotate(r["lineage"], (int(r["carried"]), float(r["retention"])), fontsize=7, xytext=(4, 3), textcoords="offset points")
    ax1.set_xscale("log")
    ax1.set_xlabel("intron positions carried (of 29,085)")
    ax1.set_ylabel("fraction retained")
    ax1.set_title("Retention of positions gained at the root (circles) and two nodes down (squares)", fontsize=9, loc="left")
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, label=g) for g, c in SC_COLOR.items()]
    ax1.legend(handles=handles, fontsize=7, frameon=False, loc="upper left")
    prof = L.read_tsv(os.path.join(a.results, "alith_origin.tsv"))
    labels = [p["origin_node"] if len(p["origin_node"]) < 55 else p["origin_node"][:52] + "..." for p in prof]
    pres = [int(p["present"]) for p in prof]; absn = [int(p["absent"]) for p in prof]; miss = [int(p["missing"]) for p in prof]
    y = range(len(prof))
    ax2.barh(y, pres, color=SC_COLOR["Antarctolithica"], label="carried")
    ax2.barh(y, absn, left=pres, color="#f4a6a6", label="observed absent (lost)")
    ax2.barh(y, miss, left=[p + q for p, q in zip(pres, absn)], color="#dddddd", label="not observable (?)")
    for i, p in enumerate(prof):
        txt = "carried by construction" if p["by_construction"] == "yes" else f"retained {float(p['retention']):.2f}"
        ax2.text(pres[i] + absn[i] + miss[i] + 50, i, txt, va="center", fontsize=7)
    ax2.set_yticks(list(y)); ax2.set_yticklabels(labels, fontsize=7)
    ax2.invert_yaxis()
    ax2.set_xlabel("intron positions gained at this node (Dollo, gene tree)")
    ax2.set_title("Antarctolithica: positions by where they arose (root at top)", fontsize=9, loc="left")
    ax2.legend(fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "fig3_retention.png"), dpi=170)
    plt.close(fig)


def fig4(a):
    p = os.path.join(a.results, "calibration_runs.tsv")
    if not os.path.exists(p):
        return
    rows = L.read_tsv(p)
    controls = sorted({r["control"] for r in rows}, key=lambda c: next(int(x["carried"]) for x in L.read_tsv(
        os.path.join(a.results, "calibration_plan.tsv")) if x["control"] == c))
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5), sharey=True)
    off = {"full": -0.25, "uniform": 0.0, "age": 0.25}
    colr = {"full": "black", "uniform": "#1f77b4", "age": "#d62728"}
    for ax, mode in zip(axes, ("all", "carrier")):
        for r in rows:
            if r["mode"] != mode or r["edges_from_true"] == "":
                continue
            i = controls.index(r["control"])
            ax.scatter(i + off[r["variant"]], int(r["edges_from_true"]), color=colr[r["variant"]], s=24, alpha=0.75)
        ax.set_xticks(range(len(controls))); ax.set_xticklabels(controls, rotation=60, ha="right", fontsize=8)
        ax.axhline(0, color="grey", lw=0.6)
        ax.set_title(f"cheapest placement vs true ({'all positions' if mode == 'all' else 'positions the stand-in carries'})",
                     fontsize=9, loc="left")
    axes[0].set_ylabel("branches between cheapest and true placement")
    handles = [plt.Line2D([], [], marker="o", ls="", color=c, label={"full": "real genome", "uniform": "thinned at random",
                                                                       "age": "thinned to Alith's age profile"}[k]) for k, c in colr.items()]
    axes[1].legend(handles=handles, fontsize=8, frameon=False)
    fig.suptitle("Calibration: known lineages thinned to Antarctolithica's intron content (controls sorted by intron richness)", fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "fig4_calibration.png"), dpi=170)
    plt.close(fig)


def fig5(a):
    rows = L.read_tsv(os.path.join(a.score, "gene_tree", "dollo_weighted.tsv"))
    g = [int(r["gain_weight"]) / max(1, int(r["loss_weight"])) for r in rows]
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.plot(g, [int(r["gains"]) for r in rows], "o-", label="gains")
    ax.plot(g, [int(r["losses"]) for r in rows], "s-", label="losses")
    ax.plot(g, [int(r["sites_several_gains"]) for r in rows], "^--", label="positions with >1 origin")
    ax.set_xscale("log"); ax.set_xlabel("cost of a gain relative to a loss"); ax.set_ylabel("count (all 29,085 positions)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Weighted Dollo on the gene tree", fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(a.out, "fig5_dollo_weighted.png"), dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("trees", "genomes", "score", "results", "out"):
        ap.add_argument("--" + k, required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sc = read_groups(a.genomes)
    fig1(a, sc); fig2(a, sc); fig3(a, sc); fig4(a); fig5(a)
    print("figures in", a.out)
