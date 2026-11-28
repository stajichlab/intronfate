#!/usr/bin/env python3
"""Write a synthetic MALIN intron table and tree with T taxa and N sites, for the B2 scaling series (no third-party data).

The tree is a random binary tree (random splits of the taxon list, fixed seed). Each site has one gain at a random node and an
independent loss (probability LOSS) on every branch below it, so Dollo, parsimony and the likelihood all have something to find.
`#MAP` lines keep the phase so `convert --from-malin` works. Usage: make_scaling_data.py OUT_PREFIX TAXA SITES [SEED]
Writes OUT_PREFIX.table and OUT_PREFIX.nwk.
"""
import random
import sys

LOSS = 0.1


def build_tree(taxa, rng):
    """Nested tuples of tip names; each internal node splits its tips at a random point."""
    if len(taxa) == 1:
        return taxa[0]
    k = rng.randint(1, len(taxa) - 1)
    return (build_tree(taxa[:k], rng), build_tree(taxa[k:], rng))


def newick(t):
    return t if isinstance(t, str) else "(" + newick(t[0]) + "," + newick(t[1]) + ")"


def nodes(tree):
    """All subtrees (the places a gain can happen), preorder; the iterative walk keeps deep trees off the recursion limit."""
    out, stack = [], [tree]
    while stack:
        n = stack.pop()
        out.append(n)
        if not isinstance(n, str):
            stack.extend(n)
    return out


def tips_with_loss(node, rng, loss, present):
    """Add the tips of `node` that keep the intron to `present`, each branch losing it with probability `loss`."""
    stack = [node]
    while stack:
        n = stack.pop()
        if rng.random() < loss:
            continue
        if isinstance(n, str):
            present.add(n)
        else:
            stack.extend(n)


def build(n_taxa, n_sites, seed=1):
    rng = random.Random(seed)
    taxa = [f"T{i:04d}" for i in range(n_taxa)]
    rng.shuffle(taxa)
    tree = build_tree(taxa, rng)
    places = nodes(tree)
    rows = {t: bytearray(b"0" * n_sites) for t in taxa}
    for s in range(n_sites):
        present = set()
        tips_with_loss(rng.choice(places), rng, LOSS, present)
        for t in present:
            rows[t][s] = ord("1")
    return sorted(taxa), rows, tree


def write(prefix, n_taxa, n_sites, seed=1):
    taxa, rows, tree = build(n_taxa, n_sites, seed)
    with open(prefix + ".table", "w") as fh:
        for t in taxa:
            fh.write(f"#ORGANISM {t} {t}\n")
        for t in taxa:
            fh.write(f"{t}\t{rows[t].decode()}\n")
        fh.write("#MAP\tsite\tkog\talipos\n")
        for i in range(n_sites):
            fh.write(f"#MAP\t{i}\tK0\t{3 + 7 * i}\n")
    with open(prefix + ".nwk", "w") as fh:
        fh.write(newick(tree) + ";\n")


if __name__ == "__main__":
    a = sys.argv
    write(a[1], int(a[2]), int(a[3]), int(a[4]) if len(a) > 4 else 1)
