#!/usr/bin/env python3
"""Constraint trees (Newick, multifurcating, for IQ-TREE -g) for the candidate placements of Alith. Each constrains only the named clades; everything else is free.
Groups come from genomes.tsv via ../V9_antarctolithica_expanded/clade_retention.py (Pezizomycotina = ascomycetes minus Saccharomycotina and Taphrinomycotina).
  H1 sister to Pezizomycotina            ((Pezizomycotina),Alith)
  H2 sister to Saccharomycotina+Pezizo   (((Saccharomycotina),(Pezizomycotina)),Alith)
  H3 sister to Taphrinomycotina          ((Taphrinomycotina),Alith)
  H4 sister to Saccharomycotina          ((Saccharomycotina),Alith)
  H5 sister to all sampled ascomycetes   ((ascomycetes),Alith)
  H6 with Saitoella+Neolecta             ((Saitoella,Neolecta),Alith)
  H7 with the Dothideomycetes black fungi (control: the placement the black-yeast assumption would give)  ((Cryomyces,Hortaea,Zymoseptoria),Alith)
  H8 Taphrinomycotina monophyletic, Alith free (control: the likelihood cost of the monophyly that H3 needs, with no claim about Alith)
Each file is the full constraint: the other lineages sit in a root polytomy.  Usage: make_hypotheses.py genomes.tsv OUTDIR"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "V9_antarctolithica_expanded"))
from clade_retention import clades_from

def hypotheses(genomes_tsv):
    C, outgroups, asco = clades_from(genomes_tsv)
    pez, sacc, taph = sorted(C["Pezizomycotina"][0]), sorted(C["Saccharomycotina"][0]), sorted(C["Taphrinomycotina"][0])
    asco = sorted(asco); A = "Alith"
    g = lambda xs: "(" + ",".join(xs) + ")"
    H = {"H1_sister_Pezizomycotina": (g(pez), [g([g(pez), A])]),
         "H2_sister_Sacch_Pezizo": (None, [g([g([g(sacc), g(pez)]), A])]),
         "H3_sister_Taphrinomycotina": (g(taph), [g([g(taph), A])]),
         "H4_sister_Saccharomycotina": (g(sacc), [g([g(sacc), A])]),
         "H5_sister_all_ascomycetes": (g(asco), [g([g(asco), A])]),
         "H6_with_Saitoella_Neolecta": (None, [g([g(["Saitoella", "Neolecta"]), A])]),
         "H7_with_Dothideomycetes_black_fungi": (None, [g([g(["Cryomyces", "Hortaea", "Zymoseptoria"]), A])]),
         "H8_Taphrinomycotina_monophyletic_Alith_free": (None, [g(taph), A])}
    out = {}
    for name, (_, clade) in H.items():
        used = set(x for c in clade for x in c.replace("(", " ").replace(")", " ").replace(",", " ").split())
        rest = [t for t in sorted(asco + list(outgroups)) if t not in used]
        out[name] = "(" + ",".join(rest + clade) + ");"
    return out

if __name__ == "__main__":
    os.makedirs(sys.argv[2], exist_ok=True)
    for k, v in hypotheses(sys.argv[1]).items():
        open(os.path.join(sys.argv[2], k + ".constraint.nwk"), "w").write(v + "\n"); print(k, v.count(","), "commas")
