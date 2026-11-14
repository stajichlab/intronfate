#!/usr/bin/env python3
"""P5: compare two Newick trees on the same tips (unrooted bipartitions): Robinson-Foulds, support-stratified agreement,
and which taxonomic groups the disagreements involve. Read-only. Run: compare_trees.py TREE_A NAME_A TREE_B NAME_B [SAMPLES_CSV] | --selftest"""
import sys, re, csv, collections

def parse(nwk):
    toks = re.findall(r"\(|\)|,|;|[^(),;:]+(?::[0-9eE.+-]+)?|:[0-9eE.+-]+", nwk.strip())
    tips, stack, nodes = [], [[]], []   # nodes: (mask, label)
    after_close = False; pending = None
    for t in toks:
        if t == "(": stack.append([]); after_close = False
        elif t == ",": after_close = False
        elif t == ")":
            kids = stack.pop(); m = 0
            for k in kids: m |= k
            stack[-1].append(m); nodes.append([m, None]); after_close = True
        elif t == ";": break
        else:
            name = t.split(":")[0]
            if t.startswith(":") or name == "":
                continue
            if after_close:
                nodes[-1][1] = name; after_close = False; continue
            tips.append(name); stack[-1].append(1 << (len(tips) - 1))
    return tips, nodes

def splits(tips, nodes, order):
    pos = {t: i for i, t in enumerate(order)}
    full = (1 << len(order)) - 1
    out = {}
    for m, lab in nodes:
        mm = 0
        for i, t in enumerate(tips):
            if (m >> i) & 1: mm |= 1 << pos[t]
        k = bin(mm).count("1")
        if k < 2 or k > len(order) - 2: continue
        if mm & 1: mm = full ^ mm                # canonical: side without tip 0
        try: v = float(lab) if lab is not None else None
        except ValueError: v = None
        out[mm] = v if mm not in out or out[mm] is None else max(out[mm], v or 0)
    return out

def selftest():
    ta, na = parse("((a,b)1,(c,d)0.5);"); tb, nb = parse("((a,c)1,(b,d)0.9);")
    o = sorted(ta); A = splits(ta, na, o); B = splits(tb, nb, o)
    assert len(A) == 1 and len(B) == 2 or len(A) == len(B) == 1, (A, B)   # (a,b|c,d) root bifurcation collapses to one split
    assert set(A) != set(B)
    tc, nc = parse("((a,b)1,(c,d)1);"); C = splits(tc, nc, o); assert set(C) == set(A)
    print("selftest ok")

def main(fa, la, fb, lb, samples="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/samples.csv"):
    ta, na = parse(open(fa).read()); tb, nb = parse(open(fb).read())
    assert sorted(ta) == sorted(tb), f"tip sets differ: {len(set(ta)^set(tb))} tips"
    order = sorted(ta); A = splits(ta, na, order); B = splits(tb, nb, order)
    inter = set(A) & set(B)
    rf = len(set(A) ^ set(B))
    print(f"tips {len(order)}; internal splits {la} {len(A)}, {lb} {len(B)}; shared {len(inter)}; RF {rf}; normalised {rf/(len(A)+len(B)):.3f}")
    for (x, X, y, Y) in ((la, A, lb, B), (lb, B, la, A)):
        print(f"\n{x} splits by support, fraction also in {y}:")
        for lo, hi in ((0, 0.5), (0.5, 0.9), (0.9, 0.99), (0.99, 1.001)):
            s = [m for m, v in X.items() if v is not None and lo <= v < hi]
            if s: print(f"  support [{lo},{hi}): n={len(s):4d} shared {sum(m in Y for m in s)/len(s):.2f}")
        print(f"  no support label: n={sum(v is None for v in X.values())}")
    tax = {r["LOCUSTAG"]: r for r in csv.DictReader(open(samples))}
    def desc(m):
        side = [order[i] for i in range(len(order)) if (m >> i) & 1]
        if len(side) > len(order) / 2: side = [order[i] for i in range(len(order)) if not (m >> i) & 1]
        c = collections.Counter((tax.get(t, {}).get("CLASS") or tax.get(t, {}).get("PHYLUM") or "?") for t in side)
        return len(side), c.most_common(3)
    for (x, X, y, Y) in ((la, A, lb, B), (lb, B, la, A)):
        bad = [(m, v) for m, v in X.items() if m not in Y and v is not None and v >= 0.95]
        print(f"\n{x} splits with support >= 0.95 absent from {y}: {len(bad)}; 25 largest-side examples (smaller side size, top classes):")
        for m, v in sorted(bad, key=lambda p: -desc(p[0])[0])[:25]:
            n, c = desc(m); print(f"   n={n:4d} support {v:.3f} {c}")
if __name__ == "__main__":
    selftest() if sys.argv[1] == "--selftest" else main(*sys.argv[1:6])
