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


def main():
    work = sys.argv[1]
    base = os.path.join(work, "soft_0.8_3.tsv")
    with open(os.path.join(RES, "support_by_taxon.tsv"), "w") as fh:
        fh.write("taxon\tcells\tfrac_reads_ge1\tfrac_reads_ge3\tfrac_reads_ge10\n")
        for r in support_by_taxon(base):
            fh.write("\t".join(map(str, r)) + "\n")
    print("== support of annotated introns by RNA-seq junction reads (cells of the combined matrix; settings 0.8/3)")
    print(open(os.path.join(RES, "support_by_taxon.tsv")).read())

    hard_dir = os.path.join(work, "score_hard")
    hard_flag, n_sites = flagged(hard_dir)
    hard_l = lnl(hard_dir)
    lines = [["setting", "ln_likelihood", "d_ln_lik_vs_hard", "sites_lnLR_gt2", "overlap_with_hard_jaccard", "sites"]]
    lines.append(["hard", hard_l, 0.0, len(hard_flag), 1.0, n_sites])
    for d in sorted(x for x in glob.glob(os.path.join(work, "score_*")) if os.path.isdir(x)):
        name = os.path.basename(d)[len("score_"):]
        if name == "hard":
            continue
        fl, n = flagged(d)
        union = len(fl | hard_flag) or 1
        lines.append([name, lnl(d), round(lnl(d) - hard_l, 3), len(fl), round(len(fl & hard_flag) / union, 4), n])
    with open(os.path.join(RES, "sensitivity.tsv"), "w") as fh:
        for r in lines:
            fh.write("\t".join(map(str, r)) + "\n")
    print("== likelihood and convergent-insertion calls (ln LR > 2) under each setting")
    print(open(os.path.join(RES, "sensitivity.tsv")).read())
    print("== warnings from `junctions`")
    for f in sorted(glob.glob(os.path.join(work, "junctions_*.stderr"))):
        w = [l for l in open(f) if l.startswith("warn")]
        print(os.path.basename(f), "none" if not w else "")
        for l in w:
            print("  ", l.rstrip())


if __name__ == "__main__":
    main()
