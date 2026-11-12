#!/usr/bin/env python3
"""Clade-specific intron retention: does a lineage carry the introns that define a clade?

A clade-specific site: >= M members of the clade (the target itself excluded) carry the intron, NO genome outside the clade carries it (the target and the
lineage under test, Alith, excluded), and >= 2 outside genomes were observed without it. An ancient site: carried by >= 1 of Mucor/Basidiomycota and by >= 2
ascomycete genomes. For a target genome T (or the Alith lineage: carrier if any strain carries) the retention of a site set is
   carriers / sites where T was observed (present or absent, not `?`).
The ratio R = retention(clade-specific) / retention(ancient) is ~ 0 for a lineage outside the clade (only convergent insertions) and of order 1 for a
lineage inside it, and it is normalised by the lineage's own intron richness, which is what biased the whole-matrix scans.
Usage: clade_retention.py matrix.tsv > table.tsv"""
import sys
CLADES = {
    "Pezizomycotina": ({"Pyronema", "Orbilia", "Botrytis", "Antarctomyces", "Ncrassa", "Geoglossum", "Lichina", "Exophiala", "Anidulans", "Cladonia", "Arthonia", "Cryomyces", "Zymoseptoria", "Hortaea"}, 3),
    "Taphrinomycotina": ({"Taphrina", "Protomyces", "Savitreella", "Spombe", "Saitoella", "Neolecta"}, 3),
    "Saccharomycotina": ({"Yarrowia", "Lipomyces"}, 2),
    "Taphrinomycetes": ({"Taphrina", "Protomyces", "Savitreella"}, 2),
    # broader clades, to ask whether the lineage is inside Ascomycota, and inside or outside Saccharomycotina + Pezizomycotina
    "Ascomycota": ({"Taphrina", "Protomyces", "Savitreella", "Spombe", "Saitoella", "Neolecta", "Yarrowia", "Lipomyces", "Pyronema", "Orbilia", "Botrytis", "Antarctomyces", "Ncrassa", "Geoglossum", "Lichina", "Exophiala", "Anidulans", "Cladonia", "Arthonia", "Cryomyces", "Zymoseptoria", "Hortaea"}, 4),
    "Sacch+Pezizo": ({"Yarrowia", "Lipomyces", "Pyronema", "Orbilia", "Botrytis", "Antarctomyces", "Ncrassa", "Geoglossum", "Lichina", "Exophiala", "Anidulans", "Cladonia", "Arthonia", "Cryomyces", "Zymoseptoria", "Hortaea"}, 4),
}
OUTSIDE_ASCO_ANCIENT = {"Mmucedo", "Cneoformans", "Smaydis"}
with open(sys.argv[1]) as f:
    lines = [l.rstrip("\n") for l in f if not l.startswith("#")]
hdr = lines[0].split("\t"); first = hdr.index("pattern") + 1; taxa = hdr[first:]
alith = [i for i, t in enumerate(taxa) if t.startswith("Alith")]
genomes = [t for t in taxa if not t.startswith("Alith")]
gi = {t: i for i, t in enumerate(taxa)}
rows = [l.split("\t")[first:] for l in lines[1:]]
def state(cells, T):
    """'1', '0' or '?' for target T (a genome name or 'Alith')."""
    if T == "Alith":
        v = [cells[i] for i in alith]
        return "1" if "1" in v else ("0" if "0" in v else "?")
    return cells[gi[T]]
asco = [t for t in genomes if t not in OUTSIDE_ASCO_ANCIENT]
print("target\tclade\tclade_specific_sites\tretention_clade_specific\tancient_sites\tretention_ancient\tR")
for T in ["Alith"] + genomes:
    p_anc = None
    cnt = {}
    n_anc = c_anc = 0
    for cells in rows:
        s = state(cells, T)
        if s == "?": continue
        car = {t for t in genomes if t != T and cells[gi[t]] == "1"}
        if (car & OUTSIDE_ASCO_ANCIENT) and len(car & set(asco)) >= 2:
            n_anc += 1; c_anc += (s == "1")
        for cname, (members, m) in CLADES.items():
            mem = members - {T}
            if len(car & mem) < m: continue
            outside = [t for t in genomes if t not in members and t != T]
            if any(cells[gi[t]] == "1" for t in outside): continue
            if sum(1 for t in outside if cells[gi[t]] == "0") < 2: continue
            a = cnt.setdefault(cname, [0, 0]); a[0] += 1; a[1] += (s == "1")
    ra = c_anc / n_anc if n_anc else float("nan")
    for cname in CLADES:
        n, c = cnt.get(cname, [0, 0])
        rc = c / n if n else float("nan")
        R = rc / ra if n and ra else float("nan")
        print(f"{T}\t{cname}\t{n}\t{rc:.4f}\t{n_anc}\t{ra:.4f}\t{R:.3f}")
