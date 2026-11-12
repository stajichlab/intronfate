#!/usr/bin/env python3
"""Prune the BFD genus-representative tree (FastTree on fungi_odb10 markers) to the genera of genomes.tsv, root it on the outgroup tip, and print the clades.
Uses make_trees.py's parse/remove/reroot (reroot is fixture-tested to keep the unrooted topology).
Usage: bfd_backbone_check.py TREE genus_representatives.tsv genomes.tsv OUTGROUP_LABEL"""
import csv, sys
src = open(__file__.replace("bfd_backbone_check.py", "make_trees.py")).read().split("ap = argparse")[0]
ns = {}; exec(src, ns)
parse, remove, reroot, leaves = ns["parse"], ns["remove"], ns["reroot"], ns["leaves"]
tree, reps, genomes, outgroup = sys.argv[1:5]
label_genus = {}
for l in open(genomes):
    if l.strip() and not l.startswith("#") and not l.startswith("Alith"):
        lab, d, g = l.rstrip("\n").split("\t"); label_genus[lab] = d.split("_")[0]
tag2label = {}; used = set()
for r in csv.DictReader(open(reps), delimiter="\t"):
    for lab, genus in label_genus.items():
        if r["GENUS"] == genus and lab not in used and r["LOCUSTAG"] not in tag2label:
            tag2label[r["LOCUSTAG"]] = lab; used.add(lab)     # one genus representative per genus: the first label of a genus gets it
print(f"{len(tag2label)} of {len(label_genus)} genomes have a genus representative in the BFD tree; none for: {sorted(set(label_genus) - used)}")
root = remove(parse(open(tree).read()), lambda nm: nm not in tag2label)
def rename(x):
    if not x["kids"]: x["name"] = tag2label[x["name"]]
    for k in x["kids"]: rename(k)
rename(root)
rr = reroot(root, outgroup)
def show(x, d=0):
    if not x["kids"]: print("  " * d + x["name"]); return
    print("  " * d + f"+ [{len(leaves(x))} tips]")
    for k in x["kids"]: show(k, d + 1)
show(rr)
