"""Fixture tests for build_supermatrix.py and make_hypotheses.py. Run: python3 tests/test_topology_setup.py"""
import os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from build_supermatrix import build
from make_hypotheses import hypotheses
d = tempfile.mkdtemp()
open(f"{d}/genomes.tsv", "w").write("# h\nM\tm\toutgroup_Mucoromycota\nS1\ts\tSaccharomycotina\nS2\ts\tSaccharomycotina\nP1\tp\tSordariomycetes\nP2\tp\tEurotiomycetes\n"
                                    "T1\tt\tTaphrinomycotina\nT2\tt\tTaphrinomycotina\nA1\ta\tAntarctolithica\nA2\ta\tAntarctolithica\n")
open(f"{d}/map.tsv", "w").write("".join(f"{p}\t{l}\n" for p, l in
    [("m1", "M"), ("s1", "S1"), ("s2", "S2"), ("p1", "P1"), ("p2", "P2"), ("t1", "T1"), ("t2", "T2"), ("a1", "A1"), ("a2", "A2")]))
os.makedirs(f"{d}/aln")
def fam(name, rows): open(f"{d}/aln/{name}.aln.fasta", "w").write("".join(f">{k}\n{v}\n" for k, v in rows.items()))
L = 120
fam("FAM1", {"m1": "A" * L, "s1": "C" * L, "s2": "D" * L, "p1": "E" * L, "p2": "F" * L, "t1": "G" * L, "t2": "H" * L, "a1": "K" * 60 + "-" * 60, "a2": "L" * L})  # a2 has more residues
fam("FAM2", {"m1": "A" * L, "s1": "C" * L, "p1": "E" * L, "t1": "G" * L, "a1": "K" * L})                                                          # 5 of 6 lineages
fam("FAM3", {"m1": "A" * 50, "s1": "C" * 50, "s2": "D" * 50, "p1": "E" * 50, "p2": "F" * 50, "t1": "G" * 50, "t2": "H" * 50, "a1": "K" * 50})            # too short
fam("FAM4", {"m1": "A" * L, "s1": "-" * L, "s2": "-" * L, "p1": "-" * L, "p2": "-" * L, "t1": "G" * L, "t2": "-" * L, "a1": "K" * L})               # gap-only columns dropped
order, seq, parts = build(f"{d}/aln", f"{d}/map.tsv", f"{d}/genomes.tsv", min_lineages=6, max_gap=0.5, min_cols=100)
assert order == ["M", "S1", "S2", "P1", "P2", "T1", "T2", "Alith"], order
assert [p[0] for p in parts] == ["FAM1"], parts          # FAM2 has 5 lineages, FAM3 is too short, FAM4 keeps one column
assert parts[0][:3] == ("FAM1", 8, L)
assert seq["Alith"][:L] == "L" * L                      # the strain with more residues was used, not a1
assert {len(v) for v in seq.values()} == {sum(p[2] for p in parts)}
order, seq, parts = build(f"{d}/aln", f"{d}/map.tsv", f"{d}/genomes.tsv", min_lineages=5, max_gap=0.5, min_cols=100)
assert [p[0] for p in parts] == ["FAM1", "FAM2"], parts
assert seq["S2"][L:] == "-" * L                          # lineage missing from FAM2 filled with gaps
H = hypotheses(f"{d}/genomes.tsv")
assert set(H) >= {"H1_sister_Pezizomycotina", "H5_sister_all_ascomycetes", "H8_Taphrinomycotina_monophyletic_Alith_free"}
for name, nwk in H.items():
    if name.startswith(('H6', 'H7')): continue            # these name real genomes (Saitoella, Cryomyces ...) absent from the fixture
    tips = nwk.replace("(", " ").replace(")", " ").replace(",", " ").replace(";", " ").split()
    assert sorted(tips) == sorted(order), (name, tips)    # every lineage once, Alith included
assert "((P1,P2),Alith)" in H["H1_sister_Pezizomycotina"]
assert "(((S1,S2),(P1,P2)),Alith)" in H["H2_sister_Sacch_Pezizo"]
assert "((T1,T2),Alith)" in H["H3_sister_Taphrinomycotina"]
print("ok")
