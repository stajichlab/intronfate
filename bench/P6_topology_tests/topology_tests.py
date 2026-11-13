#!/usr/bin/env python3
"""P6 helpers for constrained topology tests. Read-only on BFD/atlas files; writes into OUTDIR.
  topology_tests.py constraints OUTDIR            -> OUTDIR/constraint_<group>.fa (FastTree -constraints: 1 = group, 0 = rest)
  topology_tests.py concordance TREES_NW OUTDIR   -> monophyly of groups in each of the per-marker gene trees
  topology_tests.py --selftest"""
import sys, re, csv, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
SAMPLES = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/samples.csv"
ALN = "/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data/phylo/results/protein/buildtree/protein-atlas_genus_reps-taxa_1654.fungi_odb10.fa"
GROUPS = {"Microsporidia": ("PHYLUM", "Microsporidia"), "Kickxellomycota": ("PHYLUM", "Kickxellomycota"),
          "Mortierellomycota": ("PHYLUM", "Mortierellomycota"), "Glomeromycota": ("PHYLUM", "Glomeromycota")}   # last two = known-monophyletic calibration

def parse_nwk(nwk):
    toks = re.findall(r"\(|\)|,|;|[^(),;:]+(?::[0-9eE.+-]+)?|:[0-9eE.+-]+", nwk.strip())
    tips, stack, masks, after = [], [[]], [], False
    for t in toks:
        if t == "(": stack.append([]); after = False
        elif t == ",": after = False
        elif t == ")":
            m = 0
            for k in stack.pop(): m |= k
            stack[-1].append(m); masks.append(m); after = True
        elif t == ";": break
        else:
            name = t.split(":")[0]
            if t.startswith(":") or name == "": continue
            if after: after = False; continue
            tips.append(name); stack[-1].append(1 << (len(tips) - 1))
    return tips, masks

def is_split(tips, masks, members):
    """True if `members` (set of tip names, all present) is one side of an internal branch (unrooted)."""
    n = len(tips); idx = {t: i for i, t in enumerate(tips)}
    G = sum(1 << idx[t] for t in members); full = (1 << n) - 1
    return any(m == G or m == full ^ G for m in masks)

def alignment_names(path):
    return [l[1:].split()[0] for l in open(path) if l.startswith(">")]

def tax():
    return {r["LOCUSTAG"]: r for r in csv.DictReader(open(SAMPLES))}

def constraints(out):
    os.makedirs(out, exist_ok=True); names = alignment_names(ALN); T = tax()
    for g, (col, val) in GROUPS.items():
        mem = {n for n in names if n in T and T[n][col] == val}
        with open(f"{out}/constraint_{g}.fa", "w") as o:
            for n in names: o.write(f">{n}\n{'1' if n in mem else '0'}\n")
        print(g, "members", len(mem), "of", len(names))

def concordance(trees, out):
    T = tax(); lines = [l for l in open(trees) if l.strip()]
    print("gene trees", len(lines))
    res = collections.defaultdict(lambda: [0, 0])
    for l in lines:
        tips, masks = parse_nwk(l)
        for g, (col, val) in GROUPS.items():
            mem = {t for t in tips if t in T and T[t][col] == val}
            if len(mem) < 3 or len(tips) - len(mem) < 2: continue
            res[g][1] += 1
            res[g][0] += is_split(tips, masks, mem)
    for g, (a, b) in res.items(): print(f"{g}: monophyletic in {a} of {b} gene trees with >=3 members present ({a/b:.2f})" if b else g)

def selftest():
    tips, masks = parse_nwk("((a,b),(c,(d,e)));")
    assert is_split(tips, masks, {"a", "b"}) and is_split(tips, masks, {"d", "e"}) and not is_split(tips, masks, {"a", "c"})
    assert is_split(tips, masks, {"c", "d", "e"})
    print("selftest ok")

if __name__ == "__main__":
    a = sys.argv[1]
    if a == "--selftest": selftest()
    elif a == "constraints": constraints(sys.argv[2])
    elif a == "concordance": concordance(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else ".")
