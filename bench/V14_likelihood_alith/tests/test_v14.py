"""Fixture tests for the V14 helpers. Run: python3 -I tests/test_v14.py"""
import os, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.dirname(HERE))
import make_error_sets, rell  # noqa: E402


def test_error_sets():
    with tempfile.TemporaryDirectory() as d:
        v4 = os.path.join(d, "v4.tsv")
        open(v4, "w").write("genome\tscored\tsn_v2_genes\tpr_v2_genes\nA\tyes\t90\t80\nB\tyes\t100\t100\nC\tno\t\t\nZ\tyes\t50\t50\n")
        g = os.path.join(d, "g.tsv")
        open(g, "w").write("# label\tdir\tgroup\nA\ta\tx\nB\tb\tx\nC\tc\tx\nAlith_1\tq\ty\n")
        out = os.path.join(d, "e")
        make_error_sets.main(["--v4", v4, "--genomes", g, "--out", out, "--alith-sn", "0.7"])
        assert open(os.path.join(out, "none.tsv")).read() == "taxon\tsensitivity\tprecision\tsource\n"
        rows = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(os.path.join(out, "alith_sn0.70.tsv")) if not l.startswith("taxon")}
        assert rows["A"][1:3] == ["0.9000", "0.8000"] and rows["C"][1:3] == ["0.9500", "0.9000"]   # C: median of A, B (Z is not a V9 genome)
        assert rows["Alith_1"][1:3] == ["0.7000", "0.9000"]
        m = {l.split("\t")[0]: l.split("\t") for l in open(os.path.join(out, "measured.tsv")) if not l.startswith("taxon")}
        assert m["Alith_1"][1] == "0.9500"


def test_rell():
    per = {"t1": [("f1:1:0", -1.0), ("f1:2:0", -1.0), ("f2:1:0", -5.0)], "t2": [("f1:1:0", -1.5), ("f1:2:0", -1.5), ("f2:1:0", -4.0)]}
    res, nfam = rell.rell(per, 2000, 1)
    assert nfam == 2 and res["t1"][0] == -7.0 and res["t2"][0] == -7.0
    # family f1 favours t1 by 1.0, f2 favours t2 by 1.0: t1 wins when f1 is drawn more often, ties split
    assert abs(res["t1"][1] + res["t2"][1] - 1.0) < 1e-9 and 0.3 < res["t1"][1] < 0.7
    assert rell.family("a:b:3:0") == "a:b" and rell.family("3:0") == "3:0"


if __name__ == "__main__":
    test_error_sets(); test_rell(); print("ok: V14 helper tests")
