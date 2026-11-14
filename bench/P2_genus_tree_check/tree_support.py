#!/usr/bin/env python3
"""Support summary of a FastTree SH-like support tree, and Robinson-Foulds distance to another tree on the same tips.
Run: tree_support.py SUPPORT_TREE [OTHER_TREE] | --selftest. Read-only."""
import sys, re, csv, collections
sys.path.insert(0, ".")
from tree_check import parse

def support_values(nwk):
    return [float(x) for x in re.findall(r"\)([0-9]*\.?[0-9]+)(?=[:;,)])", nwk)]

def rf(tips1, sp1, tips2, sp2):
    idx = {t: i for i, t in enumerate(tips2)}
    assert set(tips1) == set(tips2), "tip sets differ"
    def canon(tips, splits, mapping):
        full = (1 << len(tips)) - 1; out = set()
        for s in splits:
            m = 0
            for i, t in enumerate(tips):
                if (s >> i) & 1: m |= 1 << mapping[t]
            c = full ^ m
            if m == full or c == 0: continue
            out.add(min(m, c))
        return out
    a = canon(tips1, sp1, idx); b = canon(tips2, sp2, idx)
    return len(a ^ b), len(a), len(b)

def selftest():
    n = "((a:1,b:1)0.95:1,(c:1,d:1)0.50:1);"
    assert support_values(n) == [0.95, 0.5], support_values(n)
    t1, s1 = parse("((a,b),(c,d));"); t2, s2 = parse("((a,c),(b,d));")
    d, n1, n2 = rf(t1, s1, t2, s2)
    assert d == 2 and rf(t1, s1, t1, s1)[0] == 0
    print("selftest ok")

def main(tree, other):
    nwk = open(tree).read(); v = support_values(nwk)
    tips, sp = parse(nwk)
    print(f"tips {len(tips)}; internal nodes with support {len(v)}")
    for thr in (0.5, 0.7, 0.9, 0.95, 0.99, 1.0):
        print(f"  support >= {thr}: {sum(x >= thr for x in v)} ({100*sum(x >= thr for x in v)/len(v):.1f}%)")
    print("  support < 0.7:", sum(x < 0.7 for x in v))
    if other:
        t2, s2 = parse(open(other).read())
        d, a, b = rf(tips, sp, t2, s2)
        print(f"RF distance to the other tree: {d} of {a + b} splits differ ({100*d/(a+b):.1f}%)")
if __name__ == "__main__":
    selftest() if sys.argv[1] == "--selftest" else main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
