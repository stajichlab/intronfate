#!/usr/bin/env python3
"""Concatenate the V9 family alignments (one MAFFT alignment per orthogroup, one protein per genome) into a supermatrix for tree-likelihood tests.
Lineages = the genomes of genomes.tsv, with the Antarctolithica strains merged into one lineage `Alith` (per family the strain with the most aligned residues
is used). A family needs >= --min-lineages lineages; columns gapped in more than --max-gap of ALL lineages are dropped; a family keeping fewer than --min-cols
columns is dropped; families are ranked by (lineages, kept columns) and added until --max-sites. Missing lineages are filled with gaps.
Outputs: <out>.fasta (supermatrix), <out>.partitions (RAxML-style, one per family), <out>.genes.tsv (family, lineages, columns, offset).
Usage: build_supermatrix.py --aln-dir ALN --taxon-map MAP --genomes genomes.tsv --out OUT [--min-lineages 35 --max-gap 0.5 --min-cols 100 --max-sites 200000]"""
import argparse, glob, os, sys

def read_fasta(path):
    seqs, name = {}, None
    for l in open(path):
        l = l.rstrip("\n")
        if l.startswith(">"): name = l[1:].split()[0]; seqs[name] = []
        elif name: seqs[name].append(l)
    return {k: "".join(v).upper() for k, v in seqs.items()}

def lineages_from(genomes_tsv, strain_group="Antarctolithica", merged="Alith"):
    """label -> lineage (strains merged); lineage order = first appearance."""
    lin, order = {}, []
    for l in open(genomes_tsv):
        if not l.strip() or l.startswith("#"): continue
        f = l.rstrip("\n").split("\t"); lab = merged if f[2] == strain_group else f[0]
        lin[f[0]] = lab
        if lab not in order: order.append(lab)
    return lin, order

def family_columns(aln, prot2lin, order, max_gap):
    """One sequence per lineage (the one with most residues), kept columns, per-lineage strings; None if no usable columns."""
    best = {}
    for p, s in aln.items():
        lab = prot2lin.get(p)
        if lab is None: continue
        n = sum(c not in "-X?" for c in s)
        if lab not in best or n > best[lab][0]: best[lab] = (n, s)
    if not best: return None
    L = len(next(iter(best.values()))[1])
    keep = [i for i in range(L) if sum(best[t][1][i] in "-X?" for t in best) + (len(order) - len(best)) <= max_gap * len(order)]
    return len(best), keep, {t: "".join(best[t][1][i] for i in keep) for t in best}

def build(aln_dir, taxon_map, genomes, min_lineages=35, max_gap=0.5, min_cols=100, max_sites=200000):
    lin, order = lineages_from(genomes)
    prot2lin = {}
    for l in open(taxon_map):
        p, lab = l.rstrip("\n").split("\t")[:2]
        if lab in lin: prot2lin[p] = lin[lab]
    fams = []
    for f in sorted(glob.glob(os.path.join(aln_dir, "*.aln.fasta"))):
        r = family_columns(read_fasta(f), prot2lin, order, max_gap)
        if r and r[0] >= min_lineages and len(r[1]) >= min_cols:
            fams.append((os.path.basename(f).split(".")[0], r[0], len(r[1]), r[2]))
    fams.sort(key=lambda x: (-x[1], -x[2], x[0]))
    chosen, total = [], 0
    for fam in fams:
        if total + fam[2] > max_sites: continue
        chosen.append(fam); total += fam[2]
    chosen.sort(key=lambda x: x[0])
    seq = {t: [] for t in order}; parts, offset = [], 0
    for name, n, w, cols in chosen:
        for t in order: seq[t].append(cols.get(t, "-" * w))
        parts.append((name, n, w, offset)); offset += w
    return order, {t: "".join(v) for t, v in seq.items()}, parts

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--aln-dir", required=True); ap.add_argument("--taxon-map", required=True); ap.add_argument("--genomes", required=True)
    ap.add_argument("--out", required=True); ap.add_argument("--min-lineages", type=int, default=35); ap.add_argument("--max-gap", type=float, default=0.5)
    ap.add_argument("--min-cols", type=int, default=100); ap.add_argument("--max-sites", type=int, default=200000)
    a = ap.parse_args()
    order, seq, parts = build(a.aln_dir, a.taxon_map, a.genomes, a.min_lineages, a.max_gap, a.min_cols, a.max_sites)
    with open(a.out + ".fasta", "w") as f:
        for t in order: f.write(f">{t}\n{seq[t]}\n")
    with open(a.out + ".partitions", "w") as f:
        for name, n, w, off in parts: f.write(f"LG, {name} = {off + 1}-{off + w}\n")
    with open(a.out + ".genes.tsv", "w") as f:
        f.write("family\tlineages\tcolumns\toffset\n")
        for p in parts: f.write("\t".join(map(str, p)) + "\n")
    tot = sum(p[2] for p in parts)
    miss = sum(seq[t].count("-") for t in order) / (len(order) * tot)
    print(f"{len(order)} lineages, {len(parts)} families, {tot} columns, {miss:.1%} gaps", file=sys.stderr)
