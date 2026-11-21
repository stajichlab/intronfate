#!/usr/bin/env python3
"""S1: score one fit against the simulated truth (standard library only).

Reads the truth table (`truth.tsv`: site, n_origins) and the fit's `likelihood_sites.tsv` (p_multi_origin, ln_lr_multi, p_no_origin) and prints one TSV row:
  * auc_multi_lnlr, auc_multi_p: how well ln LR / p(>= 2 origins) rank the sites with >= 2 true origins above the others;
  * called, precision, recall at p(>= 2 origins) > 0.5;  top100_precision: share of the 100 highest ln LR that truly have >= 2 origins;
  * auc_noorigin: how well p(no origin) ranks the sites with no true presence anywhere (every annotated presence is an error).
Run: evaluate.py TRUTH SCORE_DIR MODEL SCENARIO READS SEED | --header | --selftest
"""
import csv, os, sys

COLS = ["scenario", "reads", "seed", "model", "sites", "true_multi", "true_noorigin", "auc_multi_lnlr", "auc_multi_p", "called", "precision", "recall",
        "top100_precision", "auc_noorigin", "ml_log_likelihood", "converged"]


def auc(pos, neg):
    """P(positive scores above negative), ties half; None when a class is empty."""
    if not pos or not neg:
        return None
    allv = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    ranks, i = {}, 0
    while i < len(allv):
        j = i
        while j < len(allv) and allv[j][0] == allv[i][0]:
            j += 1
        ranks[allv[i][0]] = (i + j + 1) / 2
        i = j
    return (sum(ranks[v] for v in pos) - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def metrics(truth, sites):
    """truth {site: n_origins}; sites {site: (p_multi, ln_lr, p_none)} -> dict"""
    keys = [k for k in truth if k in sites]
    multi = {k for k in keys if truth[k] >= 2}
    none = {k for k in keys if truth[k] == 0}
    f = lambda v: None if v is None else round(v, 4)
    a_ln = auc([sites[k][1] for k in multi], [sites[k][1] for k in keys if k not in multi])
    a_p = auc([sites[k][0] for k in multi], [sites[k][0] for k in keys if k not in multi])
    called = [k for k in keys if sites[k][0] > 0.5]
    hit = sum(1 for k in called if k in multi)
    top = sorted(keys, key=lambda k: -sites[k][1])[:100]
    a_none = auc([sites[k][2] for k in none], [sites[k][2] for k in keys if k not in none])
    return {"sites": len(keys), "true_multi": len(multi), "true_noorigin": len(none), "auc_multi_lnlr": f(a_ln), "auc_multi_p": f(a_p),
            "called": len(called), "precision": f(hit / len(called)) if called else None, "recall": f(hit / len(multi)) if multi else None,
            "top100_precision": f(sum(1 for k in top if k in multi) / len(top)) if top else None, "auc_noorigin": f(a_none)}


def read_truth(path):
    return {r["site"]: int(r["n_origins"]) for r in csv.DictReader(open(path), delimiter="\t")}


def read_sites(path):
    return {r["site"]: (float(r["p_multi_origin"]), float(r["ln_lr_multi"]), float(r["p_no_origin"])) for r in csv.DictReader(open(path), delimiter="\t")}


def summary(path):
    d = {}
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        if len(f) >= 2:
            d[f[0]] = f[1]
    return d


def selftest():
    assert auc([3, 4], [1, 2]) == 1.0 and auc([1], [1]) == 0.5 and auc([], [1]) is None
    truth = {"a": 2, "b": 1, "c": 0, "d": 3}
    sites = {"a": (0.9, 3.0, 0.0), "b": (0.2, -1.0, 0.0), "c": (0.1, -2.0, 0.9), "d": (0.6, 1.0, 0.0)}
    m = metrics(truth, sites)
    assert m["true_multi"] == 2 and m["called"] == 2 and m["precision"] == 1.0 and m["recall"] == 1.0 and m["auc_multi_p"] == 1.0, m
    assert m["true_noorigin"] == 1 and m["auc_noorigin"] == 1.0 and m["top100_precision"] == 0.5, m
    print("selftest ok")


if __name__ == "__main__":
    if sys.argv[1:2] == ["--selftest"]:
        selftest()
    elif sys.argv[1:2] == ["--header"]:
        print("\t".join(COLS))
    else:
        truth_p, sd, model, scen, reads, seed = sys.argv[1:7]
        m = metrics(read_truth(truth_p), read_sites(os.path.join(sd, "likelihood_sites.tsv")))
        s = summary(os.path.join(sd, "summary_score.tsv"))
        m.update(scenario=scen, reads=reads, seed=seed, model=model, ml_log_likelihood=s.get("ml_log_likelihood", "NA"), converged=s.get("ml_converged", "NA"))
        print("\t".join("NA" if m[c] is None else str(m[c]) for c in COLS))
