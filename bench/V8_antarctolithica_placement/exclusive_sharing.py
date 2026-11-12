#!/usr/bin/env python3
"""Introns shared EXCLUSIVELY between two groups: sites whose carriers all lie in group A or group B and include a carrier of each.
This is the rare-genomic-change evidence: an intron present in Antarctolithica and in group X and nowhere else (among the genomes observed) is
a candidate synapomorphy of (Antarctolithica, X). Known sister groups calibrate the signal (Saccharomycotina + Pezizomycotina, Taphrinomycetes).
Usage: exclusive_sharing.py matrix.tsv > exclusive_sharing.tsv"""
import itertools, sys
GROUPS = {
    "Alith": lambda t: t.startswith("Alith"),
    "Mucor": lambda t: t == "Mmucedo",
    "Basidio": lambda t: t in ("Cneoformans", "Smaydis"),
    "Taphrinomycetes": lambda t: t in ("Taphrina", "Protomyces", "Savitreella"),
    "Spombe": lambda t: t == "Spombe",
    "Saitoella+Neolecta": lambda t: t in ("Saitoella", "Neolecta"),
    "Saccharomycotina": lambda t: t in ("Yarrowia", "Lipomyces"),
    "Pez_early(Pyronema,Orbilia)": lambda t: t in ("Pyronema", "Orbilia"),
    "Pez_core": lambda t: t in ("Botrytis", "Antarctomyces", "Ncrassa", "Geoglossum", "Lichina", "Exophiala", "Anidulans", "Cladonia", "Arthonia", "Cryomyces", "Zymoseptoria", "Hortaea"),
}
with open(sys.argv[1]) as f:
    lines = [l.rstrip("\n") for l in f if not l.startswith("#")]
hdr = lines[0].split("\t"); first = hdr.index("pattern") + 1; taxa = hdr[first:]
gid = {}
for t in taxa:
    gs = [g for g, p in GROUPS.items() if p(t)]
    if len(gs) != 1: sys.exit(f"{t} is in {gs}")
    gid[t] = gs[0]
names = list(GROUPS)
pair = {p: 0 for p in itertools.combinations(names, 2)}
single = {g: 0 for g in names}; carried = {g: 0 for g in names}; shared_any = {p: 0 for p in itertools.combinations(names, 2)}
n_sites = 0
for line in lines[1:]:
    c = line.split("\t"); cells = c[first:]
    carriers = {gid[t] for t, x in zip(taxa, cells) if x == "1"}
    if not carriers: continue
    n_sites += 1
    for g in carriers: carried[g] += 1
    if len(carriers) == 1: single[next(iter(carriers))] += 1
    elif len(carriers) == 2: pair[tuple(sorted(carriers, key=names.index))] += 1
    for p in shared_any:
        if p[0] in carriers and p[1] in carriers: shared_any[p] += 1
print(f"# {n_sites} sites with at least one carrier; groups: " + "; ".join(f"{g}={sum(1 for t in taxa if gid[t]==g)} genomes" for g in names))
print("group\tsites_carried\tsites_carried_by_this_group_only")
for g in names: print(f"{g}\t{carried[g]}\t{single[g]}")
print("\ngroup_a\tgroup_b\tsites_exclusive_to_the_two\tsites_carried_by_both_and_possibly_others")
for p, n in sorted(pair.items(), key=lambda kv: -kv[1]): print(f"{p[0]}\t{p[1]}\t{n}\t{shared_any[p]}")
