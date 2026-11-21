#!/usr/bin/env python3
"""J3b: compare 4-class fits stopped early (iteration caps, looser tolerances) with the 50,000-iteration fit (standard library only).

For each evidence type (`hard`, `soft`) the folder holds `<evidence>__<name>/` score folders (`summary_score.tsv`, `likelihood_branches.tsv`, `likelihood_sites.tsv`).
`<evidence>__cap50000` is the reference. Per fit it reports the log-likelihood gap, the iterations used, whether it converged, the change in the gamma shape, the largest change in
a branch gain or loss probability (and which one), the change in the total expected gains and losses, and how well the site origin table agrees with the reference (Spearman of
ln LR, overlap of the top 100). Run: analyse.py WORK | --selftest
"""
import csv, os, sys

REF = "cap50000"


def read_kv(path):
    d = {}
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 2:
            d[f[0]] = f[1]
    return d


def read_branches(path):
    """{label: (gain_prob, loss_prob)}; the root row is the root presence probability (loss 0). Labels are the tip lists, so they match across runs."""
    return {r["tips"] if r["branch"] != "root" else "root": (float(r["gain_prob"]), float(r["loss_prob"])) for r in csv.DictReader(open(path), delimiter="\t")}


def read_sites(path):
    return {r["site"]: float(r["ln_lr_multi"]) for r in csv.DictReader(open(path), delimiter="\t")}


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r, i = [0.0] * len(v), 0
    while i < len(order):
        j = i
        while j < len(order) and v[order[j]] == v[order[i]]:
            j += 1
        for k in range(i, j):
            r[order[k]] = (i + j - 1) / 2
        i = j
    return r


def spearman(a, b):
    ra, rb = ranks(a), ranks(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    den = (sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb)) ** 0.5
    return num / den if den else float("nan")


def compare(ref, run):
    """ref, run: dicts with 'sum', 'br', 'sites' -> dict of the comparison columns"""
    out = {"lnL_gap": float(ref["sum"]["ml_log_likelihood"]) - float(run["sum"]["ml_log_likelihood"]),
           "iterations": run["sum"]["ml_iterations"], "converged": run["sum"]["ml_converged"],
           "shape_diff": float(run["sum"]["ml_gamma_shape"]) - float(ref["sum"]["ml_gamma_shape"]),
           "gains_rel": float(run["sum"]["ml_expected_gains"]) / float(ref["sum"]["ml_expected_gains"]) - 1,
           "losses_rel": float(run["sum"]["ml_expected_losses"]) / float(ref["sum"]["ml_expected_losses"]) - 1}
    worst = (0.0, "", "")
    for lab, (g0, l0) in ref["br"].items():
        g1, l1 = run["br"][lab]
        for kind, a, b in (("gain", g0, g1), ("loss", l0, l1)):
            if abs(b - a) > worst[0]:
                worst = (abs(b - a), f"{kind}:{lab[:40]}", f"{a:.5f}->{b:.5f}")
    out["max_branch_change"], out["max_branch_param"], out["max_branch_values"] = worst
    keys = sorted(set(ref["sites"]) & set(run["sites"]))
    out["spearman_lnLR"] = spearman([ref["sites"][k] for k in keys], [run["sites"][k] for k in keys])
    top = lambda d: set(sorted(keys, key=lambda k: -d[k])[:100])
    out["top100_overlap"] = len(top(ref["sites"]) & top(run["sites"]))
    return out


def load(folder):
    return {"sum": read_kv(os.path.join(folder, "summary_score.tsv")), "br": read_branches(os.path.join(folder, "likelihood_branches.tsv")),
            "sites": read_sites(os.path.join(folder, "likelihood_sites.tsv"))}


COLS = ["lnL_gap", "iterations", "converged", "shape_diff", "gains_rel", "losses_rel", "max_branch_change", "max_branch_param", "max_branch_values", "spearman_lnLR", "top100_overlap"]


def main(work):
    names = sorted(d for d in os.listdir(work) if "__" in d and os.path.exists(os.path.join(work, d, "summary_score.tsv")))
    print("\t".join(["evidence", "fit"] + COLS))
    for ev in ("hard", "soft"):
        if f"{ev}__{REF}" not in names:
            continue
        ref = load(os.path.join(work, f"{ev}__{REF}"))
        for n in sorted((n for n in names if n.startswith(ev + "__")), key=lambda n: (n.split("__")[1].startswith("tol"), int("".join(c for c in n.split("__")[1] if c.isdigit()) or 0) if n.split("__")[1].startswith("cap") else 0, n)):
            c = compare(ref, load(os.path.join(work, n)))
            print("\t".join([ev, n.split("__")[1]] + [f"{c[k]:.6g}" if isinstance(c[k], float) else str(c[k]) for k in COLS]))


def selftest():
    assert spearman([1, 2, 3], [1, 2, 3]) == 1.0 and spearman([1, 2, 3], [3, 2, 1]) == -1.0
    base = {"sum": {"ml_log_likelihood": "-100", "ml_iterations": "50000", "ml_converged": "yes", "ml_gamma_shape": "3.3", "ml_expected_gains": "10", "ml_expected_losses": "20"},
            "br": {"root": (0.1, 0.0), "A": (0.01, 0.3)}, "sites": {"s1": 1.0, "s2": 2.0, "s3": 3.0}}
    run = {"sum": dict(base["sum"], ml_log_likelihood="-101", ml_iterations="500", ml_converged="no", ml_gamma_shape="3.5", ml_expected_gains="11"),
           "br": {"root": (0.1, 0.0), "A": (0.01, 0.25)}, "sites": {"s1": 3.0, "s2": 2.0, "s3": 1.0}}
    c = compare(base, run)
    assert abs(c["lnL_gap"] - 1) < 1e-12 and abs(c["shape_diff"] - 0.2) < 1e-12 and abs(c["gains_rel"] - 0.1) < 1e-12
    assert c["max_branch_param"].startswith("loss:A") and abs(c["max_branch_change"] - 0.05) < 1e-12 and c["spearman_lnLR"] == -1.0 and c["top100_overlap"] == 3
    print("selftest ok")


if __name__ == "__main__":
    selftest() if sys.argv[1:2] == ["--selftest"] else main(sys.argv[1])
