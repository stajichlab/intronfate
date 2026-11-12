#!/usr/bin/env python3
"""Introns per protein for every genome (is Antarctolithica intron-poor?). Reads intronfate's run-level introns.tsv and the proteomes.
Usage: intron_density.py introns.tsv proteomes_dir genomes.tsv > intron_density.tsv"""
import collections, csv, glob, os, sys
introns, pdir, genomes = sys.argv[1:4]
owner = {}; n_prot = collections.Counter()
for fa in glob.glob(os.path.join(pdir, "*.faa")):
    lab = os.path.basename(fa)[:-4]
    for l in open(fa):
        if l.startswith(">"): owner[l[1:].split()[0]] = lab; n_prot[lab] += 1
n_int = collections.Counter(); with_int = collections.defaultdict(set)
with open(introns, newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        lab = owner.get(r["protein_id"])
        if lab: n_int[lab] += 1; with_int[lab].add(r["protein_id"])
grp = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2] for l in open(genomes) if l.strip() and not l.startswith("#")}
print("genome\tgroup\tproteins\tintrons\tintrons_per_protein\tfraction_proteins_with_introns")
for lab in sorted(n_prot, key=lambda x: (grp.get(x, ""), x)):
    p = n_prot[lab]; i = n_int[lab]
    print(f"{lab}\t{grp.get(lab,'?')}\t{p}\t{i}\t{i/p:.2f}\t{len(with_int[lab])/p:.2f}")
