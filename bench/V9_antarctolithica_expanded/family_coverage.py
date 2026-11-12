#!/usr/bin/env python3
"""In how many families does each genome appear, per tier? A genome with many genes per family (large proteome) is often dropped by the "two or more genes = missing" rule,
so the lineage counts overstate how much it contributes. Usage: family_coverage.py families_dir genomes.tsv > family_coverage.tsv"""
import collections, glob, os, sys
fam, genomes = sys.argv[1:3]
labels = [l.split("\t")[0] for l in open(genomes) if l.strip() and not l.startswith("#")]
print("genome\t" + "\t".join(f"families_{t}\tfraction_{t}" for t in sorted(os.listdir(fam)) if os.path.isdir(os.path.join(fam, t))))
tiers = sorted(t for t in os.listdir(fam) if os.path.isdir(os.path.join(fam, t)))
count = {t: collections.Counter() for t in tiers}; nfam = {}
for t in tiers:
    n = len(glob.glob(os.path.join(fam, t, "families", "*.faa"))); nfam[t] = n
    for l in open(os.path.join(fam, t, "taxon_map.tsv")):
        prot, lab = l.rstrip("\n").split("\t"); count[t][lab] += 1
alith = [x for x in labels if x.startswith("Alith")]
for lab in labels:
    print(lab + "\t" + "\t".join(f"{count[t][lab]}\t{count[t][lab]/nfam[t]:.2f}" for t in tiers))
