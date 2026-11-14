#!/usr/bin/env python3
"""Reduced-intron lineages on the genus tree ROOTED on Chytridiomycota (outgroup = all Chytridiomycota tips, which form one
branch of the tree). Rooting uses clade bitmasks (no tree rebuilding of branch lengths). Sankoff parsimony with unit costs,
root state fixed to 0 (not reduced) if that is optimal; unknown tips (below core floors) take either state.
Reports gains (0->1, reduction) and reversions (1->0, re-expansion) under two extreme tie-breaks (changes as late as
possible = DELTRAN-like, as early as possible = ACCTRAN-like) and lists each gain clade. Read-only.
Run: reduced_rooted.py TREE DENSITY_TSV MARKER_TSV [OUTGROUP_PHYLUM] | --selftest"""
import sys, csv, collections, statistics as st
sys.path.insert(0, ".")
from reduced_mapping import parse_tree, masks

INF = 10 ** 9

def root_on(tips, children, O):
    full = (1 << len(tips)) - 1
    M = masks(tips, children)
    fam = {O, full ^ O}
    for x in M:
        if x == full or x == 0: continue
        if x & O == 0: fam.add(x)
        elif x & O == O: fam.add(full ^ x)
    fam.discard(0); fam.discard(full)
    fam = sorted(fam, key=lambda c: -bin(c).count("1"))
    node_kids = {full: []}
    for c in fam: node_kids[c] = []
    for c in fam:                                   # parent = smallest strict superset already placed (list is size-descending)
        par = full
        for d in fam:
            if d == c: break
            if (c & d) == c and bin(d).count("1") > bin(c).count("1"): par = d
        node_kids[par].append(("n", c))
    for i in range(len(tips)):                      # tip attaches to the smallest clade that contains it
        par = full
        for d in fam:
            if (d >> i) & 1: par = d
        node_kids[par].append(("t", i))
    return full, node_kids

def sankoff(tips, full, node_kids, state, tiebreak):
    order = sorted(node_kids, key=lambda c: bin(c).count("1"))     # small clades first = post-order
    cost = {}
    for c in order:
        cs = {0: 0, 1: 0}
        for k, j in node_kids[c]:
            for s in (0, 1):
                if k == "t":
                    st_ = state.get(tips[j]); b = 0 if st_ is None or st_ == s else 1   # a leaf may differ from its parent at the cost of one change
                    cs[s] += b
                else:
                    cs[s] += min(cost[j][s2] + (s != s2) for s2 in (0, 1))
        cost[c] = cs
    rs = 0 if cost[full][0] <= cost[full][1] else 1
    total = cost[full][rs]
    gains, revs, gain_clades = 0, 0, []
    stack = [(full, rs)]
    while stack:
        c, s = stack.pop()
        for k, j in node_kids[c]:
            if k == "t":                                   # change on a terminal branch (a single tip)
                obs = state.get(tips[j])
                if obs is not None and obs != s:
                    if s == 0: gains += 1; gain_clades.append(1 << j)
                    else: revs += 1
                continue
            opts = sorted((cost[j][s2] + (s != s2), (s2 == s) if tiebreak == "stay" else (s2 != s), s2) for s2 in (0, 1))
            best = min(o[0] for o in opts)
            cand = [o for o in opts if o[0] == best]
            cand.sort(key=lambda o: (not o[1], o[2]))   # prefer 'stay' (DELTRAN-like) or 'change' (ACCTRAN-like)
            s2 = cand[0][2]
            if s2 != s:
                if s == 0: gains += 1; gain_clades.append(j)
                else: revs += 1
            stack.append((j, s2))
    return total, rs, gains, revs, gain_clades

def selftest():
    t, ch, r = parse_tree("(((o1,o2),a),(b,(c,d)));")
    O = (1 << 0) | (1 << 1)
    full, nk = root_on(t, ch, O)
    state = {"o1": 0, "o2": 0, "a": 0, "b": 1, "c": 1, "d": 1}
    tot, rs, g, v, gc = sankoff(t, full, nk, state, "stay")
    assert (tot, rs, g, v) == (1, 0, 1, 0), (tot, rs, g, v)
    assert bin(gc[0]).count("1") == 3
    state2 = {"o1": 0, "o2": 0, "a": 1, "b": 0, "c": 1, "d": 1}     # a reduced, b not, c and d reduced
    tot, rs, g, v, gc = sankoff(t, full, nk, state2, "stay")
    assert tot == 2 and g == 2 and v == 0, (tot, g, v)
    print("selftest ok")

def main(tree, dens_tsv, mark_tsv, outgroup="Chytridiomycota"):
    tips, children, root = parse_tree(open(tree).read())
    dens = {r["LOCUSTAG"]: r for r in csv.DictReader(open(dens_tsv), delimiter="\t")}
    mark = {r["LOCUSTAG"]: int(r["n_marker_proteins"]) for r in csv.DictReader(open(mark_tsv), delimiter="\t")}
    tax = {r["LOCUSTAG"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/samples.csv"))}
    bycls = collections.defaultdict(list)
    for r in dens.values(): bycls[r["CLASS"]].append(float(r["cds_bp"]))
    med = {c: st.median(v) for c, v in bycls.items()}
    ok = {t: t in dens and float(dens[t]["cds_bp"]) >= 0.5 * med[dens[t]["CLASS"]] and mark.get(t, 0) >= 100 for t in tips}
    O = sum(1 << i for i, t in enumerate(tips) if t in tax and tax[t]["PHYLUM"] == outgroup)
    M = masks(tips, children); full = (1 << len(tips)) - 1
    print(f"outgroup {outgroup}: {bin(O).count('1')} tips; forms one branch: {O in M or (full ^ O) in M}")
    fullc, nk = root_on(tips, children, O)
    out_tips = [t for i, t in enumerate(tips) if (O >> i) & 1]
    for thr in (0.3, 0.5, 1.0):
        state = {t: (None if (t not in dens or not ok[t]) else (1 if float(dens[t]["introns_per_kb_cds"]) < thr else 0)) for t in tips}
        print(f"\n=== threshold < {thr}/kb, core tips only (unknown tips: {sum(v is None for v in state.values())}); outgroup tips reduced: {sum(state[t] == 1 for t in out_tips)} of {len(out_tips)}")
        for tb in ("stay", "change"):
            tot, rs, g, v, gc = sankoff(tips, fullc, nk, state, tb)
            print(f"  tie-break '{tb}': minimum changes {tot}; root state {rs}; reductions (0->1) {g}; re-expansions (1->0) {v}")
        if thr == 0.5:
            tot, rs, g, v, gc = sankoff(tips, fullc, nk, state, "stay")
            print("  reduction events (late-change reconstruction), clade size and classes:")
            for c in sorted(gc, key=lambda c: -bin(c).count("1")):
                cl = collections.Counter(tax[tips[i]]["CLASS"] or "(none)" for i in range(len(tips)) if (c >> i) & 1)
                ph = collections.Counter(tax[tips[i]]["PHYLUM"] for i in range(len(tips)) if (c >> i) & 1)
                nred = sum(state[tips[i]] == 1 for i in range(len(tips)) if (c >> i) & 1)
                print(f"     clade of {bin(c).count('1'):4d} tips ({nred:3d} reduced): {dict(cl.most_common(3))} {dict(ph.most_common(1))}")
if __name__ == "__main__":
    selftest() if sys.argv[1] == "--selftest" else main(sys.argv[1], sys.argv[2], sys.argv[3], *(sys.argv[4:5]))
