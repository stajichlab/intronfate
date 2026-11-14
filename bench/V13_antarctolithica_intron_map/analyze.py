#!/usr/bin/env python3
"""V13 tables: how *Antarctolithica* (Alith) intron positions map onto Dollo gains, losses, ancestral states and candidate placements.

Inputs (from run_v13.sbatch): the V9 tier-A combined matrix, V9's gene tree and its 39 placement trees (trees/index.tsv), and
`intronfate score` output for each tree in SCORE/<tree id>/ (gene_tree also with --reconstructions and -m dollo_weighted).

Writes into OUT/:
  ancestral_states.tsv     every node of the gene tree: label, tips, Dollo and MP number of positions present, branch gains/losses
  alith_origin.tsv         Alith-carried / lacked / unseen positions by the node where Dollo places their origin
  retention_by_origin.tsv  the same as a retention fraction for every lineage (single genomes and the Alith unit)
  placement_support.tsv    per candidate placement: Dollo losses (all sites, Alith-carrier sites), MP steps, sites that prefer it
  site_preference.tsv.gz   informative sites grouped by the set of placements that are cheapest for them
  summary.txt              the headline numbers quoted in README.md
Usage: analyze.py --matrix M --trees V9/trees --genomes V9/genomes.tsv --score score --out results
"""
import argparse, collections, csv, gzip, math, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v13lib as L  # noqa: E402

PEZ = {"Pezizomycetes", "Orbiliomycetes", "Lecanoromycetes", "Lichinomycetes", "Arthoniomycetes", "Geoglossomycetes", "Leotiomycetes",
       "Sordariomycetes", "Eurotiomycetes", "Dothideomycetes"}


def superclade(group):
    if group in PEZ:
        return "Pezizomycotina"
    return {"outgroup_Basidiomycota": "Basidiomycota", "outgroup_Mucoromycota": "Mucoromycota"}.get(group, group)


def labeller(genomes_tsv, all_tips):
    sc = {}
    for l in open(genomes_tsv):
        if l.strip() and not l.startswith("#"):
            f = l.rstrip("\n").split("\t")
            sc[f[0]] = superclade(f[2])
    members = collections.defaultdict(set)
    for t in all_tips:
        members[sc[t]].add(t)
    order = ["Mucoromycota", "Basidiomycota", "Taphrinomycotina", "Saccharomycotina", "Pezizomycotina", "Antarctolithica"]

    def label(ts):
        ts = set(ts)
        if len(ts) == 1:
            return next(iter(ts))
        if ts == set(all_tips):
            return "root (all)"
        if ts == set(all_tips) - members["Mucoromycota"]:
            return "Dikarya"
        parts = []
        for g in order:
            inside = ts & members[g]
            if not inside:
                continue
            if inside == members[g]:
                parts.append(g)
            elif len(inside) <= 3:
                parts.append(f"{g} ({','.join(sorted(inside))})")
            else:
                parts.append(f"{g} ({len(inside)} of {len(members[g])})")
        lab = " + ".join(parts)
        if ts == set(all_tips) - members["Mucoromycota"] - members["Basidiomycota"]:
            lab = "Ascomycota: " + lab
        return lab

    return label, sc


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True); ap.add_argument("--trees", required=True); ap.add_argument("--genomes", required=True)
    ap.add_argument("--score", required=True); ap.add_argument("--out", required=True); ap.add_argument("--group-prefix", default="Alith")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    W = lambda name: open(os.path.join(a.out, name), "w")
    summary = []

    root = L.parse_newick(open(os.path.join(a.trees, "gene_tree.nwk")).read())
    tips = L.leaves(root)
    label, sc = labeller(a.genomes, tips)
    alith = sorted(t for t in tips if t.startswith(a.group_prefix))
    g = os.path.join(a.score, "gene_tree")
    branches = L.read_branches(os.path.join(g, "dollo_branches.tsv"))
    bidx = L.branch_index(branches)
    sites = L.read_dollo_sites(os.path.join(g, "dollo_sites.tsv"), bidx)
    taxa, rows = L.read_matrix(a.matrix)
    anc = L.ancestral_counts(os.path.join(g, "reconstructions.tsv"), branches)
    summary.append(f"matrix: {len(rows)} positions x {len(taxa)} genomes; gene tree: {len(tips)} tips, rooted on {L.leaves(root['kids'][0])}")

    # ---- ancestral states on the gene tree
    with W("ancestral_states.tsv") as f:
        f.write("node\tlabel\tn_tips\tdollo_present\tmp_present\tbranch_gains\tbranch_losses\ttips\n")
        for b in branches:
            f.write(f"{b['node']}\t{label(b['tips'])}\t{len(b['tips'])}\t{anc.get('dollo', {}).get(b['tips'], 0)}\t"
                    f"{anc.get('mp', {}).get(b['tips'], 0)}\t{b['gains']}\t{b['losses']}\t{','.join(b['tip_list'])}\n")
    with W("tip_counts.tsv") as f:                      # observed cells, not the reconstruction (which fills in '?')
        f.write("tip\tcarried\tabsent\tmissing\n")
        for t in tips:
            c = collections.Counter(r[t] for r in rows.values())
            f.write(f"{t}\t{c['1']}\t{c['0']}\t{c['?']}\n")
    astem = frozenset(alith)
    stem = next(b for b in branches if b["tips"] == astem)
    parent = L.path_to(root, astem)[-2]
    summary.append(f"Alith stem branch ({len(alith)} strains): {stem['gains']} gains, {stem['losses']} losses (Dollo, gene tree); "
                   f"positions present: Dollo {anc['dollo'].get(astem, 0)} at the Alith ancestor, {anc['dollo'].get(L.tipset(parent), 0)} at its parent "
                   f"({label(L.tipset(parent))}); MP {anc['mp'].get(astem, 0)} and {anc['mp'].get(L.tipset(parent), 0)}")

    # ---- where Alith's positions come from
    prof = L.retention_by_origin(sites, rows, alith, root)
    with W("alith_origin.tsv") as f:
        f.write("origin_node\tn_tips\tpresent\tabsent\tmissing\tretention\tby_construction\n")
        for p in prof:
            nm = "inside Antarctolithica" if p["node"] == "inside" else label(p["node"])
            n = "-" if p["node"] == "inside" else len(p["node"])
            f.write(f"{nm}\t{n}\t{p['present']}\t{p['absent']}\t{p['missing']}\t{p['retention']:.3f}\t{'yes' if p['by_construction'] else 'no'}\n")
    carried = sum(p["present"] for p in prof)
    summary.append(f"Alith carries {carried} positions: " + "; ".join(
        f"{('inside' if p['node'] == 'inside' else label(p['node']))}: {p['present']}" + ("" if p["by_construction"] else f" (retention {p['retention']:.2f})")
        for p in prof) + "  [retention only for origins above Alith's parent; the rest are carried by construction]")

    # ---- retention by origin for every lineage (depth 0 = root, 1 = next ancestor, ...)
    units = [[t] for t in tips if not t.startswith(a.group_prefix)] + [alith]
    with W("retention_by_origin.tsv") as f:
        f.write("lineage\tsuperclade\tcarried\tdepth\torigin_node\tpresent\tabsent\tretention\tby_construction\n")
        for u in units:
            name = "Antarctolithica" if len(u) > 1 else u[0]
            rp = L.retention_by_origin(sites, rows, u, root)
            car = sum(p["present"] for p in rp)
            for d, p in enumerate(rp):
                nm = "inside" if p["node"] == "inside" else label(p["node"])
                r = "nan" if math.isnan(p["retention"]) else f"{p['retention']:.4f}"
                f.write(f"{name}\t{sc[u[0]]}\t{car}\t{d if p['node'] != 'inside' else 'inside'}\t{nm}\t{p['present']}\t{p['absent']}\t{r}\t"
                        f"{'yes' if p['by_construction'] else 'no'}\n")

    # ---- placements
    idx = {r["id"]: r for r in L.read_tsv(os.path.join(a.trees, "index.tsv"))}
    ids = ["gene_tree"] + sorted(idx)
    losses = {t: L.read_site_losses(os.path.join(a.score, t, "dollo_sites.tsv")) for t in ids if t != "gene_tree"}
    mp = {t: dict((r["key"], r["value"]) for r in L.read_tsv(os.path.join(a.score, t, "summary_score.tsv")))["mp_steps"] for t in ids}
    carrier = {s for s, r in rows.items() if L.unit_state(r, alith) == "1"}
    absent = {s for s, r in rows.items() if L.unit_state(r, alith) == "0"}
    per_tree, per_site = L.placement_support(losses, carrier)
    gt_sib = sorted(t for k in parent["kids"] if L.tipset(k) != astem for t in L.leaves(k))
    gene_place = next(t for t, r in idx.items() if sorted(r["tips_below"].split(",")) == gt_sib)
    with W("placement_support.tsv") as f:
        f.write("tree\tsister_of\tn_tips_below\tdollo_losses_all\tdelta_all\tdollo_losses_carrier\tdelta_carrier\tmp_steps\tdelta_mp\t"
                "sites_best_all\tsites_best_carrier\tsites_best_absent\tsites_unique_best\tgene_tree_position\n")
        mp_lo = min(int(mp[r["tree"]]) for r in per_tree)
        for r in per_tree:
            t = r["tree"]
            f.write(f"{t}\t{label(frozenset(idx[t]['tips_below'].split(',')))}\t{idx[t]['n_tips_below']}\t{r['losses_all']}\t{r['delta_all']}\t"
                    f"{r['losses_carrier']}\t{r['delta_carrier']}\t{mp[t]}\t{int(mp[t]) - mp_lo}\t{r['best_all']}\t{r['best_carrier']}\t"
                    f"{r['best_absent']}\t{r['unique_best_all']}\t{'yes' if t == gene_place else ''}\n")
    best_all = per_tree[0]
    best_car = min(per_tree, key=lambda r: (r["losses_carrier"], r["tree"]))
    gp = next(r for r in per_tree if r["tree"] == gene_place)
    summary.append(f"placements: {len(per_tree)}; gene-tree position {gene_place} (sister of {label(frozenset(gt_sib))}): "
                   f"+{gp['delta_all']} Dollo losses over the best (all sites), +{gp['delta_carrier']} (Alith-carrier sites)")
    summary.append(f"fewest Dollo losses, all sites: {best_all['tree']} = sister of {label(frozenset(idx[best_all['tree']]['tips_below'].split(',')))}")
    summary.append(f"fewest Dollo losses, Alith-carrier sites: {best_car['tree']} = sister of {label(frozenset(idx[best_car['tree']]['tips_below'].split(',')))}")
    summary.append(f"Alith-carrier positions {len(carrier)}, observed absent {len(absent)}; informative (placements differ): "
                   f"{sum(1 for s in per_site if s['informative'] and s['carrier'])} carrier, "
                   f"{sum(1 for s in per_site if s['informative'] and s['site'] in absent)} absent")

    # ---- which placements do the informative positions prefer?
    pat = collections.Counter()
    for s in per_site:
        if not s["informative"]:
            continue
        kind = "carrier" if s["carrier"] else ("absent" if s["site"] in absent else "missing")
        pat[(kind, s["best"])] += 1
    with gzip.open(os.path.join(a.out, "site_preference.tsv.gz"), "wt") as f:
        f.write("alith_state\tn_sites\tn_best_placements\tbest_placements\tsister_of_best\n")
        for (kind, best), n in pat.most_common():
            b = best.split(",")
            sis = "; ".join(label(frozenset(idx[t]["tips_below"].split(","))) for t in b[:6]) + (" ..." if len(b) > 6 else "")
            f.write(f"{kind}\t{n}\t{len(b)}\t{best}\t{sis}\n")
    with gzip.open(os.path.join(a.out, "site_support.tsv.gz"), "wt") as f:
        f.write("site\talith_state\tinformative\tmin_losses\tn_best\tbest\n")
        for s in per_site:
            st = L.unit_state(rows[s["site"]], alith)
            f.write(f"{s['site']}\t{st}\t{s['informative']}\t{s['min_losses']}\t{s['n_best']}\t{s['best']}\n")

    with W("summary.txt") as f:
        f.write("\n".join(summary) + "\n")
    print("\n".join(summary))


if __name__ == "__main__":
    main()
