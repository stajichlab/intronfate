#!/usr/bin/env python3
"""Candidate placements of one lineage on a backbone.

Reads OrthoFinder's rooted species tree, strips branch lengths and supports, removes the tips whose label starts with --group-prefix (the lineage under
test, kept as the subtree OrthoFinder gave it), and writes one tree per attachment edge of the remaining backbone: the lineage becomes the sister of the
node/tip on that edge. Also writes the observed gene-based tree and an index with a description of every edge.
Usage: make_trees.py SpeciesTree_rooted.txt --group-prefix Alith --genomes genomes.tsv --out trees
Outputs: trees/gene_tree.nwk, trees/backbone.nwk, trees/place_NN.nwk, trees/index.tsv (id, n_tips_below, group, tips below)"""
import argparse, os, sys

def parse(text):
    text = text.strip().rstrip(";"); i = 0; n = len(text)
    def skip():
        nonlocal i
        if i < n and text[i] == ":":
            i += 1
            while i < n and text[i] not in ",()": i += 1
    def node():
        nonlocal i
        if text[i] == "(":
            i += 1; kids = [node()]
            while text[i] == ",": i += 1; kids.append(node())
            i += 1
            while i < n and text[i] not in ",():;": i += 1       # support label
            skip(); return {"kids": kids, "name": None}
        j = i
        while j < n and text[j] not in ",():": j += 1
        nm = text[i:j]; i = j; skip(); return {"kids": [], "name": nm}
    return node()

def leaves(x): return [x["name"]] if not x["kids"] else [l for k in x["kids"] for l in leaves(k)]
def newick(x): return x["name"] if not x["kids"] else "(" + ",".join(newick(k) for k in x["kids"]) + ")"

def remove(x, pred):
    if not x["kids"]: return None if pred(x["name"]) else x
    kids = [k for k in (remove(c, pred) for c in x["kids"]) if k]
    if not kids: return None
    return kids[0] if len(kids) == 1 else {"kids": kids, "name": None}

def reroot(nd, target):
    """The same unrooted tree, rooted on the branch above leaf `target`: (target, rest)."""
    path = []
    def find(x):
        if not x["kids"]: return x["name"] == target
        for k in x["kids"]:
            if find(k): path.append((x, k)); return True
        return False
    if not find(nd): sys.exit(f"--root-on: no tip named {target}")
    cur = None
    for parent, child in reversed(path):                    # from the old root down to the target's parent: each node keeps its other children
        others = [k for k in parent["kids"] if k is not child]
        up = {"kids": others + ([cur] if cur else []), "name": None}
        cur = up["kids"][0] if len(up["kids"]) == 1 else up
    return {"kids": [{"kids": [], "name": target}, cur], "name": None}

def extract(x, pred):
    """The smallest subtree containing exactly the tips matching pred (must be a clade)."""
    L = leaves(x)
    if all(pred(l) for l in L): return x
    for k in x["kids"]:
        r = extract(k, pred)
        if r: return r
    return None

def attach_at(x, target, clade):
    if x is target: return {"kids": [x, clade], "name": None}
    return {"kids": [attach_at(k, target, clade) for k in x["kids"]], "name": x["name"]}

def edges(x, root=True):
    if not root: yield x
    for k in x["kids"]: yield from edges(k, False)

ap = argparse.ArgumentParser(); ap.add_argument("tree"); ap.add_argument("--group-prefix", required=True)
ap.add_argument("--genomes", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--root-on", default="", help="re-root the tree on the branch above this tip (OrthoFinder's own rooting is not always the outgroup)")
ap.add_argument("--drop-prefix", default="", help="remove tips with this prefix from every tree (e.g. the lineage under test, for a calibration run)")
ap.add_argument("--taph-monophyletic", action="store_true",
                help="force the Taphrinomycotina genomes into one clade, sister to the other Ascomycota (the gene tree splits them)")
a = ap.parse_args()
group_of = {}
for l in open(a.genomes):
    if l.strip() and not l.startswith("#"):
        f = l.rstrip("\n").split("\t"); group_of[f[0]] = f[2]
root = parse(open(a.tree).read())
if a.root_on: root = reroot(root, a.root_on)
if a.drop_prefix: root = remove(root, lambda nm: nm.startswith(a.drop_prefix))
is_g = lambda nm: nm.startswith(a.group_prefix)
clade = extract(root, is_g)
matching = sorted(l for l in leaves(root) if is_g(l))
if clade is None or sorted(leaves(clade)) != matching:
    sys.exit(f"the tips starting with {a.group_prefix!r} ({len(matching)}) are not one clade in the tree (the smallest clade holds "
             f"{0 if clade is None else len(leaves(clade))} of them); give the lineage's subtree explicitly")
backbone = remove(root, is_g)
def sibling_tips(x, target_clade):                         # tips of the sibling of `target_clade` in tree x
    for k in x["kids"]:
        if k is target_clade: return sorted(l for o in x["kids"] if o is not k for l in leaves(o))
        r = sibling_tips(k, target_clade)
        if r is not None: return r
    return None
os.makedirs(a.out, exist_ok=True)
_sib = sibling_tips(root, clade) or []
# when the lineage is a child of the root, its sibling is the whole rest of the tree, which is not an edge: the true position is "the root"
open(f"{a.out}/true_sibling.txt", "w").write(("ROOT" if len(_sib) == len(leaves(backbone)) else ",".join(_sib)) + "\n")
if a.taph_monophyletic:
    is_taph = lambda nm: group_of.get(nm) == "Taphrinomycotina"
    def max_clades(x):                                       # maximal clades made only of Taphrinomycotina tips
        if all(is_taph(l) for l in leaves(x)): return [x]
        return [c for k in x["kids"] for c in max_clades(k)]
    parts = max_clades(backbone)
    taph = parts[0]
    for q in parts[1:]: taph = {"kids": [taph, q], "name": None}
    rest = remove(backbone, is_taph)
    ing = {l for l in leaves(rest) if group_of.get(l, "").startswith(("Saccharomycotina", "Pezizomycetes", "Orbilio", "Lecano", "Lichino", "Arthonio", "Geoglosso", "Leotio", "Sordario", "Eurotio", "Dothideo"))}
    def smallest(x):
        for k in x["kids"]:
            r = smallest(k)
            if r: return r
        return x if ing <= set(leaves(x)) else None
    backbone = attach_at(rest, smallest(rest), taph)
os.makedirs(a.out, exist_ok=True)
open(f"{a.out}/gene_tree.nwk", "w").write(newick(root) + ";\n")
open(f"{a.out}/backbone.nwk", "w").write(newick(backbone) + ";\n")
with open(f"{a.out}/index.tsv", "w") as ix:
    ix.write("id\tn_tips_below\tgroup\ttips_below\n")
    for k, e in enumerate(edges(backbone), 1):
        tips = sorted(leaves(e)); groups = sorted({group_of.get(t, "?") for t in tips})
        gl = groups[0] if len(groups) == 1 else "mixed(" + "+".join(groups) + ")"
        open(f"{a.out}/place_{k:02d}.nwk", "w").write(newick(attach_at(backbone, e, clade)) + ";\n")
        ix.write(f"place_{k:02d}\t{len(tips)}\t{gl}\t{','.join(tips)}\n")
print(f"{len(leaves(backbone))} backbone tips, {len(leaves(clade))} strains in the clade, {k} candidate placements -> {a.out}/")
