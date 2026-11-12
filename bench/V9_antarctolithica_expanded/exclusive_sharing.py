#!/usr/bin/env python3
"""Introns shared EXCLUSIVELY between two groups: sites whose carriers all lie in group A or group B and include a carrier of each.
This is the rare-genomic-change evidence: an intron present in Antarctolithica and in group X and nowhere else (among the genomes observed) is
a candidate synapomorphy of (Antarctolithica, X). Known sister groups calibrate the signal (Saccharomycotina + Pezizomycotina, Taphrinomycetes).
Usage: exclusive_sharing.py matrix.tsv genomes.tsv > exclusive_sharing.tsv"""
import itertools, sys
def groups_from(genomes_tsv):
    """Groups for the pairwise exclusive-sharing table, from genomes.tsv (its group column) plus a few labels."""
    grp = {}
    for l in open(genomes_tsv):
        if l.strip() and not l.startswith("#"):
            f = l.rstrip("\n").split("\t"); grp[f[0]] = f[2]
    lip = {"Lipomyces", "Babjevia", "Kockiozyma", "Myxozyma", "Limtongia", "Dipodascopsis", "Neoaidaea"}
    taphm = {"Taphrina", "Tdeformans", "Protomyces", "Savitreella"}
    schizo = {"Spombe", "Sjaponicus"}
    def which(t):
        g = grp.get(t, "")
        if t.startswith("Alith"): return "Alith"
        if t == "Mmucedo": return "Mucor"
        if g == "outgroup_Basidiomycota": return "Basidio"
        if t in taphm: return "Taphrinomycetes"
        if t in schizo: return "Schizosaccharomyces"
        if t in ("Saitoella", "Neolecta"): return "Saitoella+Neolecta"
        if t in lip: return "Lipomycetes"
        if g == "Saccharomycotina": return "Sacch_other"
        if t in ("Pyronema", "Orbilia"): return "Pez_early(Pyronema,Orbilia)"
        if g.startswith("outgroup") or g in ("", "Taphrinomycotina", "Saccharomycotina"):
            sys.exit(f"{t} (group {g!r}) is not assigned to any group in exclusive_sharing.py")
        return "Pez_core"
    return which

with open(sys.argv[1]) as f:
    lines = [l.rstrip("\n") for l in f if not l.startswith("#")]
hdr = lines[0].split("\t"); first = hdr.index("pattern") + 1; taxa = hdr[first:]
which = groups_from(sys.argv[2])
gid = {t: which(t) for t in taxa}
names = ["Alith", "Mucor", "Basidio", "Taphrinomycetes", "Schizosaccharomyces", "Saitoella+Neolecta", "Lipomycetes", "Sacch_other", "Pez_early(Pyronema,Orbilia)", "Pez_core"]
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
