#!/usr/bin/env python3
"""Genomic composition of Antarctolithica (Alith) vs fungi in general, and whether it is a stand-alone group.
Inputs: ../V10_antarctolithica_domains/{selected_genomes.tsv,pfam_counts.tsv}, gene_counts.tsv (this folder), BFD tables {samples,asm_stats,aa_freq,codon_freq}.parquet (read-only).
Outputs (./out/): scalars.tsv, scalar_percentiles.tsv, neighbours.tsv, isolation.tsv, rf.tsv, pcoa_*.tsv/png, composition.log (stdout).
 1. scalar traits (genome size, GC, GC3, genes, gene density, masked fraction): Alith strains vs each group (z-score and percentile).
 2. nearest genomes by amino-acid composition (CLR), codon usage (CLR) and Pfam content (Jaccard).
 3. isolation test: for every genome, the distance to the nearest genome OUTSIDE its own order (and outside its own class); Alith (strains pooled, each strain's
    distance to the nearest non-Alith genome) is placed in that distribution. A value in the upper tail means no sampled lineage is a close relative on that feature.
 4. random forest (leave-one-order-out) on amino-acid + codon + scalar features; Alith class probabilities.
Usage: composition_analysis.py"""
import numpy as np, pandas as pd, os, sys
from scipy.spatial.distance import pdist, squareform
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
T = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/"
V10 = "../V10_antarctolithica_domains/"
os.makedirs("out", exist_ok=True)
g = pd.read_csv(V10 + "selected_genomes.tsv", sep="\t", keep_default_na=False)
sam = pd.read_parquet(T + "samples.parquet")[["ASMID", "LOCUSTAG"]]
asm = pd.read_parquet(T + "asm_stats.parquet")[["ASMID", "total_length_bp", "gc_pct", "masked_pct", "contig_count", "N50_bp"]]
gc = pd.read_csv("gene_counts.tsv", sep="\t")
g = g.merge(sam, left_on="locustag", right_on="LOCUSTAG", how="left").merge(asm, on="ASMID", how="left").merge(gc, on="locustag", how="left")
aa = pd.read_parquet(T + "aa_freq.parquet"); aa = aa[aa.species_prefix.isin(g.locustag)]
cod = pd.read_parquet(T + "codon_freq.parquet"); cod = cod[cod.species_prefix.isin(g.locustag)]
AA = aa.pivot(index="species_prefix", columns="amino_acid", values="frequency").drop(columns=["*"], errors="ignore")
STOP = {"TAA", "TAG", "TGA"}
CO = cod.pivot(index="species_prefix", columns="codon", values="frequency").drop(columns=[c for c in STOP if c in set(cod.codon)], errors="ignore")
g = g[g.locustag.isin(AA.index) & g.locustag.isin(CO.index)].reset_index(drop=True)
AA = AA.loc[g.locustag]; CO = CO.loc[g.locustag]
print("genomes with a missing amino-acid / codon frequency (filled with 0):", int(AA.isna().any(axis=1).sum()), int(CO.isna().any(axis=1).sum()))
AA = AA.fillna(0); CO = CO.fillna(0)
g["gc3"] = [(CO.loc[t, [c for c in CO.columns if c[2] in "GC"]].sum() / CO.loc[t].sum()) * 100 for t in g.locustag]
g["genome_mb"] = g.total_length_bp / 1e6; g["genes_per_mb"] = g.n_genes / g.genome_mb
from importlib import import_module
sys.path.insert(0, "../V9_antarctolithica_expanded")
TAPH = {"Taphrinomycetes", "Schizosaccharomycetes", "Neolectomycetes", "Pneumocystomycetes", "Saitoellomycetes", "Taphrinomycotina"}
SACC = {"Pichiomycetes", "Saccharomycetes", "Lipomycetes", "Dipodascomycetes", "Sporopachydermiomycetes", "Trigonopsidomycetes", "Alloascoideomycetes",
        "Serinomycetes", "Phaffomycetes", "Saccharomycodomycetes", "Saccharomycotina"}
def rf_group(r):            # same grouping as V10
    if r.group == "Antarctolithica": return "Antarctolithica"
    if r.group.startswith("outgroup_"): return r.group[len("outgroup_"):]
    k = r.cls or (r.group if r.group.endswith("mycetes") or r.group.endswith("mycotina") else "")
    if r.phylum != "Ascomycota" and r.set == "BFD": return r.phylum or "unassigned"
    if k in TAPH or r.group in TAPH: return "Taphrinomycotina"
    if k in SACC or r.group in SACC: return "Saccharomycotina"
    if r.group == "Basidiomycota" or r.phylum == "Basidiomycota": return "Basidiomycota"
    if r.group == "Mucoromycota" or r.phylum == "Mucoromycota": return "Mucoromycota"
    if k in {"", "Novakomyces", "Novakomycetes"}: return "unassigned"
    return "Pezizomycotina"
g["rf_group"] = g.apply(rf_group, axis=1)
alith = (g.group == "Antarctolithica").values
print(len(g), "genomes;", alith.sum(), "Antarctolithica strains")
SC = ["genome_mb", "gc_pct", "gc3", "n_genes", "genes_per_mb", "masked_pct"]
g[["label", "rf_group", "order"] + SC + ["contig_count", "N50_bp"]].to_csv("out/scalars.tsv", sep="\t", index=False)

# 1. scalars
rows = []
for sc in SC:
    for grp in ["Pezizomycotina", "Saccharomycotina", "Taphrinomycotina", "Basidiomycota", "Mucoromycota"]:
        v = g.loc[(g.rf_group == grp) & ~alith, sc].dropna()
        if len(v) < 5: continue
        for i in np.where(alith)[0]:
            rows.append(dict(trait=sc, group=grp, strain=g.label[i], value=round(g[sc][i], 3), group_median=round(v.median(), 3), group_q05=round(v.quantile(.05), 3),
                             group_q95=round(v.quantile(.95), 3), percentile=round((v < g[sc][i]).mean() * 100, 1)))
sp = pd.DataFrame(rows); sp.to_csv("out/scalar_percentiles.tsv", sep="\t", index=False)
print("\nScalar traits: Alith strains (median over strains) and percentile within each group (0 = below all, 100 = above all):")
for sc in SC:
    d = sp[sp.trait == sc]; vals = g.loc[alith, sc]
    print(f"  {sc:13s} Alith {vals.min():.2f}-{vals.max():.2f} | " + "; ".join(f"{gr[:5]} median {x.group_median.iloc[0]:.2f} pct {x.percentile.median():.0f}" for gr, x in d.groupby("group")))

# distances
clr = lambda M: (lambda L: L - L.mean(1, keepdims=True))(np.log(np.asarray(M, float) + 1e-4))
spaces = {"aa_clr": squareform(pdist(clr(AA.values))), "codon_clr": squareform(pdist(clr(CO.values)))}
c = pd.read_csv(V10 + "pfam_counts.tsv", sep="\t"); c = c[c.locustag.isin(g.locustag)]
pf = {p: j for j, p in enumerate(sorted(c.pfam_acc.unique()))}; ix = {t: i for i, t in enumerate(g.locustag)}
P = np.zeros((len(g), len(pf)))
for r in c.itertuples(): P[ix[r.locustag], pf[r.pfam_acc]] = 1
P = P[:, P[~alith].mean(0) >= 0.02]
spaces["pfam_jaccard"] = squareform(pdist(P > 0, "jaccard"))

# 2 + 3. neighbours and isolation
nb, iso = [], []
lvl = {"order": g.order.replace("", np.nan).fillna(g.label).values,
       "class": g.cls.replace("", np.nan).fillna(g.order.replace("", np.nan)).fillna(g.label).values}
for name, D in spaces.items():
    for i in np.where(alith)[0]:
        o = [j for j in np.argsort(D[i]) if not alith[j]][:6]
        nb += [dict(space=name, strain=g.label[i], rank=r + 1, neighbour=g.label[j], group=g.rf_group[j], d=round(D[i, j], 3)) for r, j in enumerate(o)]
    for level, lab in lvl.items():
        nn = np.array([min(D[i, j] for j in range(len(g)) if not alith[j] and lab[j] != lab[i]) if not alith[i] else np.nan for i in range(len(g))])
        base = nn[~alith & ~np.isnan(nn)]
        assert len(base) > 0 and not np.isnan(D).any(), (name, level, "NaN distances")
        al = np.array([D[i, ~alith].min() for i in np.where(alith)[0]])
        within = np.array([D[i, j] for i in np.where(alith)[0] for j in np.where(alith)[0] if i < j])
        iso.append(dict(space=name, level=level, alith_nn=round(al.mean(), 3), pct_of_genomes_with_larger_nn=round((base > al.mean()).mean() * 100, 1),
                        median_nn=round(np.median(base), 3), q95_nn=round(np.quantile(base, .95), 3), alith_strain_spread=round(within.mean(), 3),
                        most_isolated=", ".join(f"{g.label[k]}({nn[k]:.2f})" for k in np.argsort(-np.nan_to_num(nn))[:4])))
pd.DataFrame(nb).to_csv("out/neighbours.tsv", sep="\t", index=False); isod = pd.DataFrame(iso); isod.to_csv("out/isolation.tsv", sep="\t", index=False)
print("\nNearest non-Alith genomes (strain Alith_52.25 if present, else first strain):")
for name in spaces:
    d = pd.DataFrame(nb); d = d[(d.space == name) & (d.strain == (d.strain[d.strain == "Alith_52.25"].iloc[0] if (d.strain == "Alith_52.25").any() else d.strain.iloc[0]))]
    print(f"  {name:13s}", [(r.neighbour[:18], r.group[:6], r.d) for r in d.itertuples()][:6])
print("\nIsolation (distance to the nearest genome outside the own order/class; Alith strains vs all other genomes):")
print(isod.drop(columns="most_isolated").to_string(index=False))
print(isod[["space", "level", "most_isolated"]].to_string(index=False))

# 4. random forest
F = np.hstack([clr(AA.values), clr(CO.values), np.nan_to_num((g[SC] - g[SC].mean()) / g[SC].std())])
lab = g.rf_group.values; train = ~alith & (lab != "unassigned")
keep = {k for k, v in pd.Series(lab[train]).value_counts().items() if v >= 5}; train &= np.isin(lab, list(keep))
idx = np.where(train)[0]; cv = []
for tr, te in GroupKFold(n_splits=8).split(F[idx], lab[idx], g.order.values[idx]):
    m = RandomForestClassifier(500, random_state=1, n_jobs=-1, class_weight="balanced").fit(F[idx[tr]], lab[idx[tr]])
    cv += [dict(label=g.label[i], true=lab[i], pred=p) for i, p in zip(idx[te], m.predict(F[idx[te]]))]
cv = pd.DataFrame(cv); cv["ok"] = cv.true == cv.pred; print("\nRF (aa + codon + scalars) leave-order-out accuracy:\n", cv.groupby("true").ok.agg(["mean", "size"]).round(3).to_string())
m = RandomForestClassifier(1000, random_state=1, n_jobs=-1, class_weight="balanced").fit(F[train], lab[train])
pa = pd.DataFrame(m.predict_proba(F[alith]), columns=m.classes_); pa.insert(0, "strain", g.label[alith].values)
pa.loc[len(pa)] = ["mean"] + pa.iloc[:, 1:].mean().tolist(); pa.round(3).to_csv("out/rf.tsv", sep="\t", index=False); print("Alith RF probabilities:\n", pa.round(3).to_string())
