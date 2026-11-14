"""Fixture tests for summarize_ml_calibration.py. Run: python3 -I tests/test_v15.py (from the V15 folder)."""
import os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import summarize_ml_calibration as S  # noqa: E402

BACKBONE = "(M,((A,B),(C,D)));"          # control X was attached on these edges; truly sister to C
INDEX = [("place_01", "M"), ("place_02", "A"), ("place_03", "C"), ("place_04", "A,B,C,D"), ("place_05", "C,D")]


def make_calib(root):
    d = os.path.join(root, "X", "trees")
    os.makedirs(d)
    open(os.path.join(d, "backbone.nwk"), "w").write(BACKBONE)
    open(os.path.join(d, "true_sibling.txt"), "w").write("C\n")
    with open(os.path.join(d, "index.tsv"), "w") as f:
        f.write("id\tn_tips_below\tgroup\ttips_below\n")
        for i, t in INDEX:
            f.write(f"{i}\t{len(t.split(','))}\tg\t{t}\n")
    return d


def test_judge():
    with tempfile.TemporaryDirectory() as tmp:
        d = make_calib(tmp)
        r = S.judge({"place_01": -10.0, "place_02": -12.0, "place_03": -5.0, "place_04": -7.0, "place_05": -6.0}, d)
        assert (r["true_is_best"], r["edges_from_true"], r["depth_best_minus_true"], r["lnl_best_minus_true"]) == ("yes", 0, 0, 0.0)
        # best is the (A,B,C,D) stem: deeper than C by 2 edges, 3 edges away
        r = S.judge({"place_01": -10.0, "place_03": -9.5, "place_04": -7.0}, d)
        assert (r["true_is_best"], r["depth_best_minus_true"], r["lnl_best_minus_true"], r["best"]) == ("no", -2, 2.5, "place_04")
        assert r["edges_from_true"] == 2
        # ties within 1e-6 count as best
        r = S.judge({"place_03": -5.0, "place_05": -5.0000001}, d)
        assert r["true_is_best"] == "yes" and r["n_tied"] == 2


def test_main_end_to_end():
    with tempfile.TemporaryDirectory() as tmp:
        make_calib(os.path.join(tmp, "calib"))
        os.makedirs(os.path.join(tmp, "scores"))
        with open(os.path.join(tmp, "scores", "X_age_1_k4.tsv"), "w") as f:
            for t, l in [("place_01", -10), ("place_03", -9.5), ("place_04", -7)]:
                f.write(f"X\tage\t1\t4\t{t}\t{l}\tyes\t0.8\n")
        open(os.path.join(tmp, "scores", "X_age_1_k4.tsv.part"), "w").write("garbage\n")
        out = subprocess.run([sys.executable, "-I", os.path.join(os.path.dirname(HERE), "summarize_ml_calibration.py"),
                              "--scores", os.path.join(tmp, "scores", "*.tsv*"), "--calib", os.path.join(tmp, "calib"),
                              "--out", os.path.join(tmp, "res")], capture_output=True, text=True, check=True).stdout
        assert "X\tage\t4\t1\t0\t2\t-2\tA,B,C,D" in out, out
        assert "# age\t4\t1\t0\t1\t0" in out, out


if __name__ == "__main__":
    for n, f in sorted(globals().items()):
        if n.startswith("test_"):
            f(); print("ok", n)
