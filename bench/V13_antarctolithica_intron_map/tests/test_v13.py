"""Fixture tests for v13lib.py. Run: python3 -I tests/test_v13.py (from the V13 folder) or python3 -m pytest tests/."""
import os, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import v13lib as L  # noqa: E402

TREE = "(M,((A,B),(X1,X2)));"     # M outgroup; lineage X = {X1, X2}


def test_tree_helpers():
    r = L.parse_newick("(M:0.1,((A,B)90:0.2,(X1,X2)):0.3);")
    assert L.leaves(r) == ["M", "A", "B", "X1", "X2"]
    assert L.edge_distance(r, frozenset({"A"}), frozenset({"B"})) == 2
    assert L.edge_distance(r, frozenset({"A"}), frozenset({"X1", "X2"})) == 3
    assert [len(L.tipset(x)) for x in L.path_to(r, frozenset({"X1"}))] == [5, 4, 2, 1]
    assert L.recon_name(["A", "B", "C", "D", "E"]) == "nodeA-B-C-D…(5 tips)"
    assert L.recon_name(["A", "B"]) == "nodeA-B"


def test_unit_state():
    assert L.unit_state({"X1": "0", "X2": "1"}, ["X1", "X2"]) == "1"
    assert L.unit_state({"X1": "0", "X2": "?"}, ["X1", "X2"]) == "0"
    assert L.unit_state({"X1": "?", "X2": "?"}, ["X1", "X2"]) == "?"


def test_origin_profile_and_retention():
    root = L.parse_newick(TREE)
    allt = frozenset("M A B X1 X2".split()); ing = frozenset("A B X1 X2".split()); X = frozenset({"X1", "X2"})
    sites = {"s1": {"gained_at": allt}, "s2": {"gained_at": allt}, "s3": {"gained_at": ing}, "s4": {"gained_at": X},
             "s5": {"gained_at": frozenset({"X1"})}, "s6": {"gained_at": frozenset({"A", "B"})}, "s7": {"gained_at": ing}}
    rows = {"s1": {"X1": "1", "X2": "0"}, "s2": {"X1": "0", "X2": "0"}, "s3": {"X1": "?", "X2": "?"}, "s4": {"X1": "1", "X2": "1"},
            "s5": {"X1": "1", "X2": "0"}, "s6": {"X1": "0", "X2": "0"}, "s7": {"X1": "1", "X2": "?"}}
    p = L.origin_profile(sites, rows, ["X1", "X2"], root)
    assert [(len(x["node"]) if x["node"] != "inside" else "in", x["present"], x["absent"], x["missing"]) for x in p] == \
        [(5, 1, 1, 0), (4, 1, 0, 1), (2, 1, 0, 0), ("in", 1, 0, 0)]      # s6 is gained on A,B: not on X's path, not counted
    assert [x["by_construction"] for x in p] == [False, True, True, True]     # parent of X, X's stem, inside
    r = L.retention_by_origin(sites, rows, ["X1", "X2"], root)
    assert r[0]["retention"] == 0.5 and r[1]["retention"] == 1.0


def test_placement_support():
    by_tree = {"p1": {"a": 1, "b": 0, "c": 2}, "p2": {"a": 0, "b": 0, "c": 2}, "p3": {"a": 0, "b": 1, "c": 2}}
    rows, per_site = L.placement_support(by_tree, carrier_sites={"a"})
    r = {x["tree"]: x for x in rows}
    assert rows[0]["tree"] == "p2" and r["p2"]["delta_all"] == 0 and r["p1"]["delta_all"] == 1 and r["p3"]["delta_all"] == 1
    assert (r["p2"]["best_carrier"], r["p3"]["best_carrier"], r["p1"]["best_carrier"]) == (1, 1, 0)
    assert (r["p1"]["best_absent"], r["p2"]["best_absent"], r["p3"]["best_absent"]) == (1, 1, 0)
    assert r["p1"]["losses_carrier"] == 1 and r["p2"]["delta_carrier"] == 0
    s = {x["site"]: x for x in per_site}
    assert s["c"]["informative"] is False and s["c"]["best"] == "all" and s["a"]["best"] == "p2,p3"


def test_ancestral_counts():
    with tempfile.TemporaryDirectory() as d:
        b = os.path.join(d, "b.tsv")
        open(b, "w").write("branch\tnode\tn_tips\tgains\tlosses\ttips\nroot\t0\t5\t1\t0\tM,A,B,C,D\nA,B,C,D\t1\t4\t0\t0\tA,B,C,D\nM\t2\t1\t0\t0\tM\n")
        r = os.path.join(d, "r.tsv")
        open(r, "w").write("method\tsite\tnode\tis_tip\tstate\n"
                           "dollo\ts\tnodeM-A-B-C…(5 tips)\tfalse\tpresent\ndollo\ts\tnodeA-B-C-D\tfalse\tlost\ndollo\ts\tM\ttrue\tpresent\n"
                           "mp\ts\tnodeA-B-C-D\tfalse\tpresent\n")
        c = L.ancestral_counts(r, L.read_branches(b))
        assert c["dollo"] == {frozenset("MABCD"): 1, frozenset({"M"}): 1}
        assert c["mp"] == {frozenset("ABCD"): 1}


def test_thin_and_carrier_only():
    with tempfile.TemporaryDirectory() as d:
        src = os.path.join(d, "m.tsv")
        rows = "".join(f"s{i}\t{i}\t0\t2\t3\t11?\t1\t1\t?\n" for i in range(200))
        open(src, "w").write("# format=2\nsite\tcolumn\tphase\tn_taxa\tn_observable\tpattern\tA\tX\tB\n" + rows)
        dst = os.path.join(d, "t.tsv")
        L.thin_matrix(src, dst, "X", 0.25, seed=1)
        taxa, m = L.read_matrix(dst)
        kept = sum(r["X"] == "1" for r in m.values())
        assert taxa == ["A", "X", "B"] and 30 < kept < 70                              # about a quarter of 200
        assert all(r["A"] == "1" for r in m.values())                                  # other taxa untouched
        line = [l for l in open(dst) if l.startswith("s") and "\t0\t?\n" in l][0].split("\t")
        assert line[5] == "10?" and line[3] == "1"                                     # pattern and n_taxa recomputed
        L.thin_matrix(src, dst, "X", lambda s: 1.0 if s == "s3" else 0.0, seed=1)      # per-site probability
        _, m = L.read_matrix(dst)
        assert [s for s, r in m.items() if r["X"] == "1"] == ["s3"]
        co = os.path.join(d, "c.tsv")
        L.carrier_only(dst, co, "X")
        assert list(L.read_matrix(co)[1]) == ["s3"]


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f()
    print("ok: v13lib tests")
