#!/usr/bin/env python3
"""Compare each tier's terminal-branch Dollo losses (and gains) with what annotation error could explain (V4 method).

Loss ceiling  = sites_carried * (1 - Sn) / Sn   : introns the annotation missed, if every one sat at a site other species share
Gain ceiling  = sites_carried * (1 - Pr)        : introns in the annotation but absent from RefSeq (possible false positives)
sites_carried = intron rows of that genome in the tier's sites.tsv. Sn/Pr = V4 accuracy at the V2 BUSCO genes
(../V4_annotation_accuracy/results/busco_accuracy.tsv). The V6 gene sets overlap V2's but are not identical, so this is an approximation.
S. cerevisiae S288C has no gene-restricted score: the genome-wide Sn/Pr (73.3 / 66.0, 312 reference introns) is used and flagged.
Usage: loss_ceiling_tiers.py results_projection ../V4_annotation_accuracy/results/busco_accuracy.tsv
"""
import collections, csv, glob, os, sys
root, acc = sys.argv[1], sys.argv[2]
SN = {}; PR = {}
for r in csv.DictReader(open(acc), delimiter="\t"):
    if r["sn_v2_genes"]: SN[r["genome"]] = float(r["sn_v2_genes"]) / 100; PR[r["genome"]] = float(r["pr_v2_genes"]) / 100
SN["Scerevisiae"], PR["Scerevisiae"] = 0.733, 0.660
for tier in sorted(d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d))):
    carried = collections.Counter(); loss = collections.Counter(); gain = collections.Counter()
    for aln in glob.glob(os.path.join(root, tier, "*.aln")):
        with open(os.path.join(aln, "sites.tsv")) as f:
            for r in csv.DictReader(f, delimiter="\t"): carried[r["taxon"]] += 1
        with open(os.path.join(aln, "dollo_branches.tsv")) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["n_tips"] == "1": loss[r["tips"]] += int(r["losses"]); gain[r["tips"]] += int(r["gains"])
    print(f"== {tier}\n  {'genome':12s} {'carried':>8s} {'term.loss':>9s} {'Sn':>5s} {'ceiling':>8s} {'share':>7s} | {'term.gain':>9s} {'Pr':>5s} {'ceiling':>8s} {'share':>7s}")
    for g in sorted(carried):
        sn, pr = SN[g], PR[g]; c = carried[g]
        lc = c * (1 - sn) / sn; gc = c * (1 - pr)
        ls = f"{100*lc/loss[g]:6.1f}%" if loss[g] else "      -"
        gs = f"{100*gc/gain[g]:6.1f}%" if gain[g] else "      -"
        flag = "  (genome-wide Sn/Pr, 312 ref introns)" if g == "Scerevisiae" else ""
        print(f"  {g:12s} {c:8d} {loss[g]:9d} {100*sn:5.1f} {lc:8.0f} {ls} | {gain[g]:9d} {100*pr:5.1f} {gc:8.0f} {gs}{flag}")
