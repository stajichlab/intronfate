"""Fixture tests for perturb.py: extreme Sn/Pr values give exactly predictable matrices. Run: python3 tests/test_perturb.py"""
import os, random, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import perturb as P

HDR = ["site", "column", "phase", "n_taxa", "n_observable", "pattern", "A", "B", "C"]
def rows(*pats):
    return [[f"s{i}", str(i), "0", str(sum(c == "1" for c in p)), str(sum(c in "01" for c in p)), p] + list(p) for i, p in enumerate(pats)]
def cells(out, j): return [x[6 + j] for x in out]

def test_identity_when_perfect():
    r = rows("110", "1?0", "001")
    cand = {"A": [2], "B": [2], "C": [0]}
    assert P.perturb(HDR, r, {"A": 1, "B": 1, "C": 1}, {"A": 1, "B": 1, "C": 1}, "ceiling", random.Random(1), cand) == r

def test_all_false_positives_removed_and_empty_sites_dropped():
    assert P.perturb(HDR, rows("110", "001"), {}, {"A": 0.0, "B": 0.0, "C": 0.0}, "ceiling", random.Random(1), {}) == []

def test_budget_is_pooled_over_the_whole_matrix():
    # A carries 3 introns; Sn=0.5 and Pr=1 => M = 3 misses, to be placed at A's candidate cells (rows 3 and 4), which are only 2: both filled.
    # The old per-family cap gave each family its own budget; here the budget is the taxon's whole genome.
    r = rows("100", "100", "100", "011", "011")
    out = P.perturb(HDR, r, {"A": 0.5}, {}, "ceiling_fn", random.Random(1), {"A": [3, 4]})
    assert cells(out, 0) == ["1"] * 5
    # a budget smaller than the candidates fills exactly that many: A has 1 intron, M = 1
    r = rows("100", "011", "011", "011")
    out = P.perturb(HDR, r, {"A": 0.5}, {}, "ceiling_fn", random.Random(1), {"A": [1, 2, 3]})
    assert sum(1 for c in cells(out, 0) if c == "1") == 2                    # its own intron + exactly one fill

def test_missed_count_uses_true_positives():
    # Pr = 0.5, Sn = 0.5, A has 4 introns: TP = 2, M = 2 * 0.5 / 0.5 = 2 fills (not 4)
    r = rows("100", "100", "100", "100", "011", "011", "011", "011")
    out = P.perturb(HDR, r, {"A": 0.5}, {"A": 0.5}, "ceiling_fn", random.Random(1), {"A": [4, 5, 6, 7]})
    assert sum(1 for c in cells(out, 0) if c == "1") == 4 + 2

def test_question_marks_and_non_candidates_are_never_filled():
    r = rows("100", "0?1", "011")
    out = P.perturb(HDR, r, {"A": 0.1}, {}, "ceiling_fn", random.Random(1), {"A": [2]})   # M = 9: far more than candidates; row 1 is not a candidate
    assert [x[6] for x in out] == ["1", "0", "1"] and out[1][7] == "?"

def test_fp_and_fn_modes_are_separate():
    r = rows("011", "100")
    fn = P.perturb(HDR, r, {"A": 0.5}, {"A": 1.0}, "ceiling_fn", random.Random(1), {"A": [0]})   # M = 1 * 1 * 0.5/0.5 = 1 fill, none removed
    assert cells(fn, 0) == ["1", "1"]
    fn_only = P.perturb(HDR, r, {"A": 0.5}, {"A": 0.5}, "ceiling_fn", random.Random(1), {"A": [0]})  # the "_fn" modes never remove introns
    assert cells(fn_only, 0).count("1") >= 1 and len(fn_only) == 2
    fp = P.perturb(HDR, r, {"A": 0.5}, {"A": 0.0}, "fp", random.Random(1), {"A": [0]})           # Sn ignored; every A intron removed
    assert all(x != "1" for x in cells(fp, 0)) and len(fp) == 1                                   # row 1 (A only) loses its carrier and is dropped

def test_typical_is_at_most_ceiling():
    r = rows("011", "011", "101", "110", "100")
    cand = {"A": [0, 1]}
    def n(mode): return sum(sum(1 for x in P.perturb(HDR, r, {"A": 0.5}, {}, mode, random.Random(s), cand) if x[6] == "1") for s in range(50))
    assert n("typical_fn") <= n("ceiling_fn")

def test_extra_leading_columns_survive_and_columns_are_recomputed():
    h = ["site", "column", "phase", "family", "family_column", "n_taxa", "n_observable", "pattern", "A", "B", "C"]
    r = [["f1:5:0", "5", "0", "f1", "5", "2", "3", "110"] + list("110"), ["f2:7:1", "10007", "1", "f2", "7", "1", "2", "01?"] + list("01?")]
    out = P.perturb(h, r, {"A": 0.6, "B": 0.6, "C": 0.6}, {"A": 0.7, "B": 0.7, "C": 0.7}, "ceiling", random.Random(3), {"A": [1]})
    for x in out:
        assert x[:5] in (r[0][:5], r[1][:5])                                   # site, column, phase, family, family_column untouched
        assert x[7] == "".join(x[8:]) and int(x[5]) == x[7].count("1") and int(x[6]) == sum(c in "01" for c in x[7])

def test_candidates_follow_the_dollo_output():
    with tempfile.TemporaryDirectory() as d:
        p = os.path.join(d, "dollo_sites.tsv")
        open(p, "w").write("site\tcolumn\tphase\tn_present\tn_missing\tgains\tlosses\tmonophyletic\tdollo_compatible\tgained_at\tlost_in\n"
                           "s0\t0\t0\t2\t0\t1\t2\tno\tno\troot\tA;B,C\n" "s1\t1\t0\t1\t0\t1\t0\tyes\tyes\tB\t\n" "s2\t2\t0\t2\t0\t1\t1\tno\tno\troot\tC\n")
        cand = P.candidates_from_dollo(p, {"s0": 0, "s1": 1, "s2": 2}, ["A", "B", "C"])
        assert cand == {"A": [0], "B": [], "C": [2]}                           # 'B,C' is an internal branch, not B's terminal branch

if __name__ == "__main__":
    for k, v in list(globals().items()):
        if k.startswith("test_"): v()
    print("ok")
