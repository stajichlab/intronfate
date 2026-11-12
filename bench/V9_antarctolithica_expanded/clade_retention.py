#!/usr/bin/env python3
"""Clade-specific intron retention: does a lineage carry the introns that define a clade?

A clade-specific site: >= M members of the clade (the target itself excluded) carry the intron, NO genome outside the clade carries it (the target and the
lineage under test, Alith, excluded), and >= 2 outside genomes were observed without it. An ancient site: carried by >= 1 outgroup genome and by >= 2
ascomycete genomes. For a target genome T (or the Alith lineage: carrier if any strain carries) the retention of a site set is
   carriers / sites where T was observed (present or absent, not `?`).
R = retention(clade-specific) / retention(ancient) is ~ 0 for a lineage outside the clade (only convergent insertions) and of order 1 for a lineage inside it,
and it is normalised by the lineage's own intron richness, which is what biased the whole-matrix scans.
Clades come from genomes.tsv (its group column): Pezizomycotina = Ascomycota minus Saccharomycotina and Taphrinomycotina; the outgroups are the genomes whose
group starts with `outgroup`; Taphrinomycetes and Lipomycetes are given by label below.
Usage: clade_retention.py matrix.tsv genomes.tsv > table.tsv"""
import sys
TAPHRINOMYCETES = {"Taphrina", "Tdeformans", "Protomyces", "Savitreella"}
LIPOMYCETES = {"Lipomyces", "Babjevia", "Kockiozyma", "Myxozyma", "Limtongia", "Dipodascopsis", "Neoaidaea"}

def clades_from(genomes_tsv):
    group = {}
    for l in open(genomes_tsv):
        if l.strip() and not l.startswith("#"):
            f = l.rstrip("\n").split("\t"); group[f[0]] = f[2]
    outgroups = {t for t, g in group.items() if g.startswith("outgroup")}
    asco = {t for t, g in group.items() if not g.startswith("outgroup") and g != "Antarctolithica"}
    sacc = {t for t in asco if group[t] == "Saccharomycotina"}
    taph = {t for t in asco if group[t] == "Taphrinomycotina"}
    pez = asco - sacc - taph
    clades = {"Pezizomycotina": (pez, 3), "Taphrinomycotina": (taph, 3), "Taphrinomycetes": (TAPHRINOMYCETES & asco, 2),
              "Saccharomycotina": (sacc, 3), "Lipomycetes": (LIPOMYCETES & asco, 3),
              "Ascomycota": (asco, 4), "Sacch+Pezizo": (sacc | pez, 4)}
    return clades, outgroups, asco

if __name__ == "__main__":
    CLADES, OUTGROUPS, ASCO = clades_from(sys.argv[2])
    with open(sys.argv[1]) as f:
        lines = [l.rstrip("\n") for l in f if not l.startswith("#")]
    hdr = lines[0].split("\t"); first = hdr.index("pattern") + 1; taxa = hdr[first:]
    alith = [i for i, t in enumerate(taxa) if t.startswith("Alith")]
    genomes = [t for t in taxa if not t.startswith("Alith")]
    gi = {t: i for i, t in enumerate(taxa)}
    rows = [l.split("\t")[first:] for l in lines[1:]]
    def state(cells, T):
        if T == "Alith":
            v = [cells[i] for i in alith]
            return "1" if "1" in v else ("0" if "0" in v else "?")
        return cells[gi[T]]
    print("target\tclade\tclade_specific_sites\tretention_clade_specific\tancient_sites\tretention_ancient\tR")
    for T in ["Alith"] + genomes:
        cnt = {}; n_anc = c_anc = 0
        for cells in rows:
            s = state(cells, T)
            if s == "?": continue
            car = {t for t in genomes if t != T and cells[gi[t]] == "1"}
            if (car & OUTGROUPS) and len(car & ASCO) >= 2:
                n_anc += 1; c_anc += (s == "1")
            for cname, (members, m) in CLADES.items():
                members = members & set(genomes); mem = members - {T}
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
