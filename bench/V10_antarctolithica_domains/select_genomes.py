#!/usr/bin/env python3
"""Choose the genomes for the Pfam comparison: all V9 genomes (by the locus tag in their proteome ids) + the Antarctolithica strains
+ a broad BFD sample (up to N_PER_ORDER genomes per taxonomic order, best fungal BUSCO first, fixed seed).
Reads BFD tables/{samples,busco_genome}.parquet read-only.  Usage: select_genomes.py V9_DIR OUT.tsv"""
import sys, glob, os, random
import pandas as pd
T = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/"
v9, out = sys.argv[1], sys.argv[2]
N_PER_ORDER = 3
s = pd.read_parquet(T + "samples.parquet")
b = pd.read_parquet(T + "busco_genome.parquet")
b = b[b.lineage.str.contains("fungi", na=False)].sort_values("complete_pct", ascending=False).drop_duplicates("ASMID")
s = s.merge(b[["ASMID", "complete_pct", "duplicated_pct"]], on="ASMID", how="left")
rows = []
# V9 genomes
grp = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2] for l in open(f"{v9}/genomes.tsv") if l.strip() and not l.startswith("#")}
v9tags = {}
for f in sorted(glob.glob(f"{v9}/proteomes/*.faa")):
    label = os.path.basename(f)[:-4]
    tag = open(f).readline()[1:].split("_")[0]
    v9tags[tag] = label
    g = label.split("_")[0] if label.startswith("Alith") else label
    rows.append(dict(locustag=tag, label=label, set="V9", group=grp.get(label, "Antarctolithica" if label.startswith("Alith") else "?"),
                     phylum="", subphylum="", cls="", order="", complete_pct=float("nan")))
tax = s.set_index("LOCUSTAG")
for r in rows:
    if r["locustag"] in tax.index:
        t = tax.loc[r["locustag"]]
        r.update(phylum=t.PHYLUM or "", subphylum=t.SUBPHYLUM or "", cls=t.CLASS or "", order=t.ORDER or "", complete_pct=t.complete_pct)
# broad sample
cand = s[(s.complete_pct >= 80) & (s.duplicated_pct < 10) & s.ORDER.notna() & ~s.LOCUSTAG.isin(v9tags)]
random.seed(1)
for order, d in cand.groupby("ORDER"):
    d = d.sort_values("complete_pct", ascending=False).head(15).sample(frac=1, random_state=1).head(N_PER_ORDER)
    for _, t in d.iterrows():
        rows.append(dict(locustag=t.LOCUSTAG, label=f"{t.GENUS}_{t.SPECIES}_{t.STRAIN}".replace(" ", "_")[:40] + "_" + t.LOCUSTAG, set="BFD",
                         group=(t.SUBPHYLUM or t.PHYLUM or "?"), phylum=t.PHYLUM or "", subphylum=t.SUBPHYLUM or "", cls=t.CLASS or "",
                         order=order, complete_pct=t.complete_pct))
pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
d = pd.DataFrame(rows); print(len(d), "genomes"); print(d.groupby(["set", "group"]).size().sort_values(ascending=False).head(25))
