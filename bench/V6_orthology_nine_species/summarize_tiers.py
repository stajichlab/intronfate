#!/usr/bin/env python3
"""Per tier: families, intron sites, and Dollo gains/losses summed per branch (from each alignment's dollo_branches.tsv).
Usage: summarize_tiers.py results_projection > results_projection/REPORT.txt   (a branch is identified by the tips below it)"""
import collections, csv, glob, os, sys
root = sys.argv[1]
for tier in sorted(d for d in os.listdir(root) if os.path.isdir(os.path.join(root, d))):
    gains = collections.Counter(); losses = collections.Counter(); n_tips = {}; fams = 0; sites = 0
    for aln in glob.glob(os.path.join(root, tier, "*.aln")):
        fams += 1
        with open(os.path.join(aln, "sites.tsv")) as f: sites += sum(1 for l in f if l.strip() and not l.startswith(("#", "site")))
        with open(os.path.join(aln, "dollo_branches.tsv")) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                gains[r["tips"]] += int(r["gains"]); losses[r["tips"]] += int(r["losses"]); n_tips[r["tips"]] = int(r["n_tips"])
    print(f"== {tier}: {fams} families, {sites} intron rows (taxon x site, rows of sites.tsv), Dollo gains {sum(gains.values())}, losses {sum(losses.values())}")
    print(f"  {'branch (tips below)':70s} {'gains':>7s} {'losses':>7s}")
    for k in sorted(gains, key=lambda k: (-n_tips[k], k)):
        print(f"  {k[:70]:70s} {gains[k]:7d} {losses[k]:7d}")
