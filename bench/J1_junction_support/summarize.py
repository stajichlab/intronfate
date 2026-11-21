#!/usr/bin/env python3
"""J1 summaries from the tables the job leaves in a work folder (standard library only).

  summarize.py WORK > results/summary.txt     and writes results/support_by_taxon.tsv, results/sensitivity.tsv next to this script.

WORK holds `soft_<p0>_<k>.tsv` (the `intronfate junctions` tables), `junctions_<p0>_<k>.stderr`, and one folder `score_<name>/`
per likelihood run (`score_hard`, `score_<p0>_<k>`), each with `likelihood_sites.tsv` and `summary_score.tsv`.
"""
import csv, glob, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")


def read_tsv(path):
    with open(path) as fh:
        rows = [l.rstrip("\n").split("\t") for l in fh if not l.startswith("#") and l.strip()]
    head = rows[0]
    return head, [dict(zip(head, r)) for r in rows[1:]]


def support_by_taxon(soft_path):
    _, rows = read_tsv(soft_path)
    by = defaultdict(list)
    for r in rows:
        by[r["taxon"]].append(int(r["reads"]))
    out = []
    for t in sorted(by):
        v = by[t]
        n = len(v)
        out.append([t, n] + [round(sum(1 for x in v if x >= k) / n, 4) for k in (1, 3, 10)])
    return out


def lnl(score_dir):
    for l in open(os.path.join(score_dir, "summary_score.tsv")):
        f = l.rstrip("\n").split("\t")
        if f[0] == "ml_log_likelihood":
            return float(f[1])
    return float("nan")


def flagged(score_dir, thr=2.0):
    _, rows = read_tsv(os.path.join(score_dir, "likelihood_sites.tsv"))
    return {r["site"] for r in rows if float(r["ln_lr_multi"]) > thr}, len(rows)


def scores(score_dir, col):
    """site -> value of `col` (ln_lr_multi or p_multi_origin) in a score folder."""
    _, rows = read_tsv(os.path.join(score_dir, "likelihood_sites.tsv"))
    return {r["site"]: float(r[col]) for r in rows}


def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    i = 0
    while i < len(order):          # average ranks over ties
        j = i
        while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2
        i = j + 1
    return r


def spearman(a, b):
    ra, rb = ranks(a), ranks(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    den = (sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb)) ** 0.5
    return num / den if den else float("nan")


def top_overlap(a, b, n=100):
    ta = {s for s, _ in sorted(a.items(), key=lambda kv: -kv[1])[:n]}
    tb = {s for s, _ in sorted(b.items(), key=lambda kv: -kv[1])[:n]}
    return len(ta & tb) / n


def main():
    work = sys.argv[1]
    tag = sys.argv[2] if len(sys.argv) > 2 else ""
    suffix = f".{tag}" if tag else ""
    base = os.path.join(work, "soft_0.8_3.tsv")
    with open(os.path.join(RES, "support_by_taxon.tsv"), "w") as fh:  # the same for every tag (it only reads the junction tables)
        fh.write("taxon\tcells\tfrac_reads_ge1\tfrac_reads_ge3\tfrac_reads_ge10\n")
        for r in support_by_taxon(base):
            fh.write("\t".join(map(str, r)) + "\n")
    print("== support of annotated introns by RNA-seq junction reads (cells of the combined matrix; settings 0.8/3)")
    print(open(os.path.join(RES, "support_by_taxon.tsv")).read())

    hard_dir = os.path.join(work, f"score_{tag}hard")
    hard_flag, n_sites = flagged(hard_dir)
    hard_l = lnl(hard_dir)
    hard_lr = scores(hard_dir, "ln_lr_multi")
    hard_pm = scores(hard_dir, "p_multi_origin")
    sites = sorted(hard_lr)
    lines = [["setting", "ln_likelihood", "sites_lnLR_gt2", "spearman_lnLR_vs_hard", "top100_lnLR_overlap", "mean_p_multi", "sites_p_multi_gt0.5", "sites"]]
    lines.append(["hard", hard_l, len(hard_flag), 1.0, 1.0, round(sum(hard_pm.values()) / n_sites, 4), sum(1 for v in hard_pm.values() if v > 0.5), n_sites])
    for d in sorted(x for x in glob.glob(os.path.join(work, f"score_{tag}*")) if os.path.isdir(x)):
        name = os.path.basename(d)[len(f"score_{tag}"):]
        if name == "hard" or (not tag and name.startswith("rc")):
            continue
        fl, n = flagged(d)
        union = len(fl | hard_flag) or 1
        lr, pm = scores(d, "ln_lr_multi"), scores(d, "p_multi_origin")
        rho = spearman([hard_lr[s] for s in sites], [lr[s] for s in sites])
        lines.append([name, lnl(d), len(fl), round(rho, 4), top_overlap(hard_lr, lr), round(sum(pm.values()) / n, 4), sum(1 for v in pm.values() if v > 0.5), n])
    with open(os.path.join(RES, f"sensitivity{suffix}.tsv"), "w") as fh:
        for r in lines:
            fh.write("\t".join(map(str, r)) + "\n")
    print("== likelihood and convergent-insertion calls (ln LR > 2) under each setting")
    print(open(os.path.join(RES, f"sensitivity{suffix}.tsv")).read())
    print("== warnings from `junctions`")
    for f in sorted(glob.glob(os.path.join(work, "junctions_*.stderr"))):
        w = [l for l in open(f) if l.startswith("warn")]
        print(os.path.basename(f), "none" if not w else "")
        for l in w:
            print("  ", l.rstrip())


if __name__ == "__main__":
    main()
