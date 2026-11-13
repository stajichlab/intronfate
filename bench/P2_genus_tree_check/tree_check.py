#!/usr/bin/env python3
"""Check a Newick tree against BFD taxonomy: tip count, and for each taxon group (PHYLUM, CLASS, ORDER) whether it is
monophyletic on the UNROOTED tree (a branch splits exactly the group from the rest) and its best F1 over all branches.
Tips are LOCUSTAGs, taxonomy from samples.csv. Run: tree_check.py TREE [SAMPLES_CSV] | --selftest. Read-only."""
import sys, re, csv, collections

def parse(nwk):
    """iterative Newick parse -> list of tip names, list of splits as int bitmasks over tips (each internal branch)"""
    toks = re.findall(r"\(|\)|,|;|[^(),;:]+(?::[0-9eE.+-]+)?|:[0-9eE.+-]+", nwk.strip())
    tips, stack, splits = [], [[]], []
    after_close = False
    for t in toks:
        if t == "(":
            stack.append([]); after_close = False
        elif t == ",":
            after_close = False
        elif t == ")":
            kids = stack.pop(); m = 0
            for k in kids: m |= k
            stack[-1].append(m); splits.append(m); after_close = True
        elif t == ";":
            break
        else:
            name = t.split(":")[0]
            if name == "" or after_close:    # branch length only, or internal node label / support value
                after_close = False
                continue
            tips.append(name); stack[-1].append(1 << (len(tips) - 1))
    return tips, splits

def group_stats(tips, splits, groups, minn=3):
    full = (1 << len(tips)) - 1
    sp = set(splits) | {full ^ s for s in splits}
    out = []
    for g, members in groups.items():
        idx = [i for i, t in enumerate(tips) if t in members]
        if len(idx) < minn: continue
        S = 0
        for i in idx: S |= 1 << i
        mono = S in sp or len(idx) == len(tips)
        best, f1 = None, 0.0
        if not mono:
            for s in sp:
                inter = bin(s & S).count("1")
                if inter == 0: continue
                p = inter / bin(s).count("1"); r = inter / len(idx); f = 2 * p * r / (p + r)
                if f > f1: f1, best = f, (inter, bin(s).count("1"))
        out.append((g, len(idx), mono, 1.0 if mono else f1, best))
    return out

def selftest():
    t, s = parse("((a:1,b:1):1,((c:1,d:1):1,(e:1,f:1):1):1);")
    assert t == list("abcdef")
    g = {"AB": {"a", "b"}, "CD": {"c", "d"}, "ACE": {"a", "c", "e"}, "ABCD": {"a", "b", "c", "d"}}
    r = {x[0]: x for x in group_stats(t, s, g, minn=2)}
    assert r["AB"][2] and r["CD"][2] and r["ABCD"][2] and not r["ACE"][2]
    t2, s2 = parse("((a:1,b:1)0.95:1,(c:1,d:1)0.80:1);")   # internal support labels must not become tips
    assert t2 == list("abcd"), t2
    print("selftest ok")

def main(tree, samples):
    nwk = open(tree).read()
    tips, splits = parse(nwk)
    print("tips:", len(tips), " internal branches:", len(splits))
    tax = {}
    for r in csv.DictReader(open(samples)):
        tax[r["LOCUSTAG"]] = r
    miss = [t for t in tips if t not in tax]
    print("tips without taxonomy:", len(miss), miss[:5])
    for level in ("PHYLUM", "CLASS", "ORDER"):
        groups = collections.defaultdict(set)
        for t in tips:
            if t in tax and tax[t][level]: groups[tax[t][level]].add(t)
        res = group_stats(tips, splits, groups, minn=3 if level != "PHYLUM" else 2)
        mono = [x for x in res if x[2]]
        print(f"\n{level}: {len(res)} groups with >=3 tips (>=2 for phylum); monophyletic {len(mono)}; tips in monophyletic groups {sum(x[1] for x in mono)} of {sum(x[1] for x in res)}")
        for g, n, m, f1, best in sorted(res, key=lambda x: -x[1])[: (20 if level != "ORDER" else 25)]:
            print(f"   {g[:28]:28s} n={n:4d} {'MONOPHYLETIC' if m else 'not mono, best F1=%.2f' % f1}")

if __name__ == "__main__":
    if sys.argv[1] == "--selftest": selftest()
    else: main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/samples.csv")
