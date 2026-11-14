#!/usr/bin/env python3
"""List out-of-place tips: for each CLASS (n>=5) that is not monophyletic on the unrooted tree, take the branch with the
best F1 and report (a) class members outside that clade ("strays": found elsewhere) and (b) non-members inside it
("intruders": a different-class genome sitting in the clade), with species and genus. Read-only.
Run: tree_outliers.py TREE [SAMPLES_CSV] | --selftest"""
import sys, csv, collections
sys.path.insert(0, ".")
from tree_check import parse

def best_split(tips, splits, S):
    full = (1 << len(tips)) - 1
    sp = set(splits) | {full ^ s for s in splits}
    best, f1 = None, -1.0
    n = bin(S).count("1")
    for s in sp:
        inter = bin(s & S).count("1")
        if not inter: continue
        p = inter / bin(s).count("1"); r = inter / n; f = 2 * p * r / (p + r)
        if f > f1: f1, best = f, s
    return best, f1

def outliers(tips, splits, members):
    S = 0
    for i, t in enumerate(tips):
        if t in members: S |= 1 << i
    b, f1 = best_split(tips, splits, S)
    strays = [tips[i] for i in range(len(tips)) if (S >> i) & 1 and not (b >> i) & 1]
    intr = [tips[i] for i in range(len(tips)) if (b >> i) & 1 and not (S >> i) & 1]
    return strays, intr, f1

def selftest():
    t, s = parse("(((a,b),c),((e,f),(g,(d,h))));")   # class {a,b,c,d}: d is a stray; best clade is {a,b,c} (F1 0.857)
    strays, intr, f1 = outliers(t, s, {"a", "b", "c", "d"})
    assert strays == ["d"] and intr == [], (strays, intr)
    t, s = parse("((a,b),(c,(d,x)),(e,f));")     # class {a,b,c,d}: x sits inside the clade of the class
    strays, intr, f1 = outliers(t, s, {"a", "b", "c", "d"})
    assert intr == ["x"] or strays == [] , (strays, intr)
    print("selftest ok")

def main(tree, samples):
    tips, splits = parse(open(tree).read())
    tax = {r["LOCUSTAG"]: r for r in csv.DictReader(open(samples))}
    by = collections.defaultdict(set)
    for t in tips: by[tax[t]["CLASS"]].add(t)
    rows = []
    for c, mem in sorted(by.items(), key=lambda kv: -len(kv[1])):
        if len(mem) < 5 or not c: continue
        strays, intr, f1 = outliers(tips, splits, mem)
        if not strays and not intr: continue
        print(f"\n{c} (n={len(mem)}, best F1={f1:.3f}): {len(strays)} stray(s) outside the main clade, {len(intr)} intruder(s) inside it")
        for t in strays: print(f"   stray    {t} {tax[t]['GENUS']:22s} {tax[t]['SPECIES'][:34]:34s} ORDER={tax[t]['ORDER']}")
        for t in intr:  print(f"   intruder {t} {tax[t]['GENUS']:22s} {tax[t]['SPECIES'][:34]:34s} CLASS={tax[t]['CLASS']} ORDER={tax[t]['ORDER']}")
if __name__ == "__main__":
    selftest() if sys.argv[1] == "--selftest" else main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/samples.csv")
