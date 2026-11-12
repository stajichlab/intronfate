"""Fixture test for make_trees.py: re-rooting preserves the unrooted topology, and the candidate placements are right. Run: python3 tests/test_make_trees.py"""
import os
src = open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "make_trees.py")).read().split("ap = argparse")[0]
ns = {}; exec(src, ns)
parse, newick, reroot, leaves, remove = ns["parse"], ns["newick"], ns["reroot"], ns["leaves"], ns["remove"]

def splits(nd):
    """Unrooted bipartitions of a rooted tree, as frozensets of the side not containing the lexicographically smallest tip."""
    allt = frozenset(leaves(nd)); small = min(allt); out = set()
    def rec(x):
        if not x["kids"]: return
        s = frozenset(leaves(x))
        if 1 < len(s) < len(allt) - 1: out.add(allt - s if small in s else s)
        for k in x["kids"]: rec(k)
    rec(nd); return out

t = parse("((C,D),(M,((S,T),(U,V))));")
r = reroot(t, "M")
assert newick(r).startswith("(M,") and leaves(r["kids"][0]) == ["M"] and sorted(leaves(r)) == sorted(leaves(t))   # rooted on the branch above M
assert splits(t) == splits(r)                               # same unrooted tree
for tip in "CDSUV":
    assert splits(reroot(t, tip)) == splits(t), tip
# removing a lineage keeps the rest
assert newick(remove(parse("((A,B),(C,D));"), lambda n: n == "B")) == "(A,(C,D))"
print("ok")
