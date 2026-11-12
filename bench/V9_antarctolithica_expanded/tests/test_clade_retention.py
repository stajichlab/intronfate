"""Fixture test for clade_retention.py (V9: clades read from genomes.tsv): a hand-built matrix with known answers. Run: python3 tests/test_clade_retention.py"""
import os, subprocess, sys, tempfile
here = os.path.dirname(os.path.abspath(__file__)); script = os.path.join(os.path.dirname(here), "clade_retention.py")
GENOMES = ["Mmucedo", "Cneoformans", "Smaydis", "Taphrina", "Protomyces", "Savitreella", "Spombe", "Saitoella", "Neolecta", "Yarrowia", "Lipomyces",
           "Pyronema", "Orbilia", "Botrytis", "Antarctomyces", "Ncrassa", "Geoglossum", "Lichina", "Exophiala", "Anidulans", "Cladonia", "Arthonia",
           "Cryomyces", "Zymoseptoria", "Hortaea"]
TAXA = GENOMES + ["Alith_1", "Alith_2"]

def row(carriers, observed_absent=None):
    cells = {t: "?" for t in TAXA}
    for t in GENOMES: cells[t] = "0"                       # everything observed absent unless listed
    for t in carriers: cells[t] = "1"
    for t, v in (observed_absent or {}).items(): cells[t] = v
    return [cells[t] for t in TAXA]

rows = []
# 5 Pezizomycotina-specific sites: carried by Orbilia, Pyronema, Botrytis; Hortaea carries 2 of them; Alith is observed absent at all 5
for i in range(5):
    carriers = ["Orbilia", "Pyronema", "Botrytis"] + (["Hortaea"] if i < 2 else [])
    rows.append(row(carriers, {"Alith_1": "0", "Alith_2": "?"}))
# 4 ancient sites: Mucor + two Pezizomycotina; Hortaea carries 2, Alith (via Alith_2) carries 1
for i in range(4):
    carriers = ["Mmucedo", "Orbilia", "Pyronema"] + (["Hortaea"] if i < 2 else [])
    rows.append(row(carriers, {"Alith_1": "?", "Alith_2": "1" if i == 0 else "0"}))
GROUPS = {t: "Taphrinomycotina" for t in ["Taphrina", "Protomyces", "Savitreella", "Spombe", "Saitoella", "Neolecta"]}
GROUPS.update({t: "Saccharomycotina" for t in ["Yarrowia", "Lipomyces"]})
GROUPS.update({t: "outgroup_x" for t in ["Mmucedo", "Cneoformans", "Smaydis"]})
hdr = ["site", "column", "phase", "n_taxa", "n_observable", "pattern"] + TAXA
with tempfile.TemporaryDirectory() as d:
    p = os.path.join(d, "m.tsv")
    with open(p, "w") as f:
        f.write("# test\n" + "\t".join(hdr) + "\n")
        for i, r in enumerate(rows): f.write("\t".join([f"s{i}", str(i), "0", "0", "0", "x"] + r) + "\n")
    g = os.path.join(d, "genomes.tsv")
    with open(g, "w") as f:
        for t in GENOMES: f.write(f"{t}\t{t}_dir\t{GROUPS.get(t, 'Pezizomycotina_class')}\n")
    out = subprocess.run([sys.executable, "-I", script, p, g], capture_output=True, text=True, check=True).stdout
res = {(l.split("\t")[0], l.split("\t")[1]): l.split("\t") for l in out.strip().split("\n")[1:]}
h = res[("Hortaea", "Pezizomycotina")]
assert h[2] == "5" and abs(float(h[3]) - 0.4) < 1e-9, h                  # 5 clade-specific sites, 2 carried
assert h[4] == "4" and abs(float(h[5]) - 0.5) < 1e-9 and abs(float(h[6]) - 0.8) < 1e-9, h   # ancient retention 0.5, R = 0.4/0.5
a = res[("Alith", "Pezizomycotina")]
assert a[2] == "5" and float(a[3]) == 0.0, a                              # an outsider-like lineage: carries none of the clade-specific introns
assert a[4] == "4" and abs(float(a[5]) - 0.25) < 1e-9, a                  # Alith counted as carrier if ANY strain carries: 1 of 4 ancient sites
t = res[("Taphrina", "Pezizomycotina")]
assert float(t[3]) == 0.0, t                                              # a non-member that carries none
print("ok")
