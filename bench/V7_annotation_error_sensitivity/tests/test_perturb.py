"""Fixture tests for perturb.py: extreme Sn/Pr values give exactly predictable matrices. Run: python3 tests/test_perturb.py"""
import os, random, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import perturb as P
TREE = P.parse_newick("((A,B),C);")

HDR = ["site", "column", "phase", "n_taxa", "n_observable", "pattern", "A", "B", "C"]
def rows(*pats):
    return [[f"{i}:0", str(i), "0", str(sum(c == "1" for c in p)), str(sum(c in "01" for c in p)), p] + list(p) for i, p in enumerate(pats)]

def test_identity_when_perfect():
    r = rows("110", "1?0", "001")
    assert P.perturb(HDR, r, {"A": 1, "B": 1, "C": 1}, {"A": 1, "B": 1, "C": 1}, "ceiling", random.Random(1), TREE) == r

def test_all_false_positives_removed_and_empty_sites_dropped():
    r = rows("110", "001")
    out = P.perturb(HDR, r, {}, {"A": 0.0, "B": 0.0, "C": 0.0}, "ceiling", random.Random(1), TREE)
    assert out == []                                    # Pr = 0 removes every intron, so no site keeps a carrier

def test_missed_introns_fill_only_terminal_loss_cells():
    # tree ((A,B),C). Row 0: carriers B and C -> MRCA is the root, A is inside it and its sibling B carries: a terminal loss on A (candidate).
    # Row 1: carrier C only -> A is outside the MRCA (gain on C's branch): not a candidate. Row 2: carrier B only: MRCA is the B tip, A outside.
    r = rows("011", "001", "010", "100")        # row 3: A's own intron, so M = 1 * 0.5 / 0.5 = 1 missed
    out = P.perturb(HDR, r, {"A": 0.5}, {}, "ceiling", random.Random(1), TREE)
    assert [x[6] for x in out] == ["1", "0", "0", "1"]  # only row 0 is filled
    assert P.terminal_loss_candidates(0, ["A", "B", "C"], [list(x[6:]) for x in r], *TREE) == [0]

def test_fp_and_fn_modes_are_separate():
    r = rows("011", "001", "010", "100")
    fn = P.perturb(HDR, r, {"A": 0.5}, {"A": 0.0}, "ceiling_fn", random.Random(1), TREE)      # Pr=0 would delete everything, but fn ignores it
    assert [x[6] for x in fn] == ["1", "0", "0", "1"]
    fp = P.perturb(HDR, r, {"A": 0.5}, {"A": 0.0}, "fp", random.Random(1), TREE)              # Sn ignored; every A intron removed
    assert all(x[6] != "1" for x in fp) and len(fp) == 3                                       # row 3 (A only) loses its carrier and is dropped

def test_question_marks_are_never_changed():
    r = rows("0?1", "01?")
    out = P.perturb(HDR, r, {"A": 0.5, "B": 0.5, "C": 0.5}, {"A": 0.9, "B": 0.9, "C": 0.9}, "ceiling", random.Random(2), TREE)
    for x in out:
        for orig in r:
            if orig[0] == x[0]: assert all(o != "?" or n == "?" for o, n in zip(orig[6:], x[6:]))

def test_typical_is_at_most_ceiling():
    r = rows("011", "011", "101", "110")
    def n(mode):
        return sum(sum(1 for x in P.perturb(HDR, r, {"A": 0.5}, {}, mode, random.Random(s), TREE) if x[6] == "1") for s in range(50))
    assert n("typical") <= n("ceiling")

def test_recomputed_columns_consistent():
    out = P.perturb(HDR, rows("110", "1?0", "001"), {"A": 0.6, "B": 0.6, "C": 0.6}, {"A": 0.7, "B": 0.7, "C": 0.7}, "ceiling", random.Random(3), TREE)
    for x in out:
        assert x[5] == "".join(x[6:]) and int(x[3]) == x[5].count("1") and int(x[4]) == sum(c in "01" for c in x[5])

if __name__ == "__main__":
    for k, v in list(globals().items()):
        if k.startswith("test_"): v()
    print("ok")
