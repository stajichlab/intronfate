#!/usr/bin/env python3
"""Perturb one intronfate `matrix.tsv` the way annotation error would, given each taxon's measured sensitivity (Sn) and precision (Pr).

The observed matrix O is what the annotation says. We do not know the truth T, so we draw plausible T from O by undoing the errors:
  * missed introns (false negatives): a taxon with sensitivity Sn is missing about  M = present * (1 - Sn) / Sn  introns. A miss only matters
    to the Dollo score when it hides a loss: the taxon is '0' at a site whose gain is inferred above it (it lies inside the clade of the site's
    carriers) and a relative in its parent clade carries the site, so Dollo puts a loss on the taxon's own terminal branch. Those cells are the
    candidates; setting one to '1' is "the loss was really a missed intron". Placing misses on arbitrary '0' cells would instead create sharing
    with distant lineages and push gains deeper, adding losses elsewhere, which is not what a missed intron at a lost-looking site does.
    `ceiling` fills min(M, candidates) of them (the V4 worst case: every miss hides a terminal loss); `typical` scales M by f, the fraction of the
    taxon's own introns that sit at sites with >= 2 carriers (missed introns are assumed to be distributed like annotated ones).
  * false introns (false positives): every '1' is turned to '0' with probability 1 - Pr.
'?' cells are never changed. Sites left with no carrier are dropped; n_taxa, n_observable and pattern are recomputed.
"""
import random

def read_matrix(path):
    with open(path) as f:
        lines = f.read().rstrip("\n").split("\n")
    meta = lines[0] if lines[0].startswith("#") else None
    hdr = lines[1 if meta else 0].split("\t")
    rows = [l.split("\t") for l in lines[(2 if meta else 1):] if l]
    return meta, hdr, rows

def write_matrix(path, meta, hdr, rows):
    with open(path, "w") as f:
        if meta: f.write(meta + "\n")
        f.write("\t".join(hdr) + "\n")
        for r in rows: f.write("\t".join(r) + "\n")

def parse_newick(text):
    """Return (leaf_sets, parent_leaves): all clades as frozensets, and {leaf: its parent's clade}. Branch lengths and labels are ignored."""
    text = text.strip().rstrip(";"); pos = [0]; clades = []; parent = {}
    def node():
        if text[pos[0]] == "(":
            pos[0] += 1; kids = [node()]
            while text[pos[0]] == ",": pos[0] += 1; kids.append(node())
            assert text[pos[0]] == ")"; pos[0] += 1
            leaves = frozenset().union(*kids); label()
            clades.append(leaves)
            for k in kids:
                if len(k) == 1: parent[next(iter(k))] = leaves
            return leaves
        name = label(); return frozenset([name])
    def label():
        i = pos[0]
        while pos[0] < len(text) and text[pos[0]] not in ",()": pos[0] += 1
        tok = text[i:pos[0]]; return tok.split(":")[0]
    root = node(); return sorted(clades, key=len), parent

def terminal_loss_candidates(taxon_idx, taxa, cells, clades, parent):
    """Row indices where `taxa[taxon_idx]` is '0' and Dollo would put a loss on its terminal branch."""
    t = taxa[taxon_idx]; sib = parent[t] - {t}; out = []
    for i, c in enumerate(cells):
        if c[taxon_idx] != "0": continue
        carriers = {taxa[j] for j, x in enumerate(c) if x == "1"}
        if not carriers or not (carriers & sib): continue
        if len(carriers) == 1: continue                    # gain on that tip's own branch: nothing above it to lose from
        mrca = next(cl for cl in clades if carriers <= cl)
        if t in mrca: out.append(i)
    return out

def stochastic_round(x, rng):
    k = int(x)
    return k + (1 if rng.random() < x - k else 0)

def perturb(hdr, rows, sn, pr, mode, rng, tree=None):
    """sn, pr: {taxon: fraction}. mode: 'ceiling' or 'typical'. tree: (clades, parent) from parse_newick. Returns new rows."""
    first = hdr.index("pattern") + 1
    taxa = hdr[first:]
    cells = [r[first:] for r in rows]                       # mutable copy of the cell grid
    orig = [list(c) for c in cells]
    carriers = [sum(1 for x in c if x == "1") for c in orig]
    for j, t in enumerate(taxa):
        present = [i for i in range(len(rows)) if orig[i][j] == "1"]
        zeros = terminal_loss_candidates(j, taxa, orig, *tree)
        s, p = sn.get(t, 1.0), pr.get(t, 1.0)
        if present and 0 < s < 1 and mode != "fp":
            missed = len(present) * (1 - s) / s
            if mode.startswith("typical"):
                missed *= sum(1 for i in present if carriers[i] >= 2) / len(present)
            k = min(stochastic_round(missed, rng), len(zeros))
            for i in rng.sample(zeros, k): cells[i][j] = "1"
        if p < 1 and not mode.endswith("_fn"):
            for i in present:
                if rng.random() < 1 - p: cells[i][j] = "0"
    out = []
    for r, c in zip(rows, cells):
        n1 = sum(1 for x in c if x == "1")
        if n1 == 0: continue
        nobs = sum(1 for x in c if x in "01")
        out.append(r[:3] + [str(n1), str(nobs), "".join(c)] + c)    # keeps site column phase; n_taxa n_observable pattern recomputed
    return out
