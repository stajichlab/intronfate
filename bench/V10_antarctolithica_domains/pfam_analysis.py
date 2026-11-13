#!/usr/bin/env python3
"""Pfam domain composition of Antarctolithica (Alith) strains against fungi in general.
Inputs (this folder): selected_genomes.tsv, pfam_counts.tsv (locustag, pfam_acc, n_domains, n_proteins).
Outputs (./out/): pcoa_{braycurtis,jaccard}.tsv (+ png if matplotlib), neighbours.tsv, rf_cv.tsv, rf_alith.tsv, pfam_signature.tsv, core_presence.tsv.
  1. PCoA (classical MDS) of Bray-Curtis on relative domain counts and of Jaccard on presence/absence.
  2. nearest genomes to each Alith strain (both distances).
  3. Random forest (group = subphylum, else phylum) trained on all non-Alith genomes; accuracy by leave-one-ORDER-out CV (so a lineage is never
     predicted from its own order); then the Alith class probabilities. Calibration: the same CV restricted to genomes with an Alith-sized Pfam repertoire.
  4. Pfam signature test (the Pfam analogue of the intron retention test): for each group, Pfams present in >=80 % of its genomes and <=10 % of all others;
     the fraction of them each genome carries (self excluded from the definition), members vs non-members vs Alith.
  5. core presence: fraction of Pfams present in >=95 % of non-Alith genomes that each genome carries (Alith vs the distribution).
Usage: pfam_analysis.py"""
import numpy as np, pandas as pd, os, sys
from scipy.spatial.distance import pdist, squareform
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
os.makedirs("out", exist_ok=True)
g = pd.read_csv("selected_genomes.tsv", sep="\t", keep_default_na=False)
c = pd.read_csv("pfam_counts.tsv", sep="\t")
g = g[g.locustag.isin(c.locustag.unique())].reset_index(drop=True)
alith = (g.group == "Antarctolithica").values
TAPH = {"Taphrinomycetes", "Schizosaccharomycetes", "Neolectomycetes", "Pneumocystomycetes", "Saitoellomycetes", "Taphrinomycotina"}
SACC = {"Pichiomycetes", "Saccharomycetes", "Lipomycetes", "Dipodascomycetes", "Sporopachydermiomycetes", "Trigonopsidomycetes", "Alloascoideomycetes",
        "Serinomycetes", "Phaffomycetes", "Saccharomycodomycetes", "Saccharomycotina"}
UNK = {"", "Novakomycetes"}                     # class missing or uncertain placement: kept in the plots, not used for training
def rf_group(r):
    if r.group == "Antarctolithica": return "Antarctolithica"
    if r.group.startswith("outgroup_"): return r.group[len("outgroup_"):]
    k = r.cls or (r.group if r.group.endswith("mycetes") or r.group.endswith("mycotina") else "")
    if r.phylum != "Ascomycota" and r.set == "BFD": return r.phylum or "unassigned"
    if k in TAPH or r.group in TAPH: return "Taphrinomycotina"
    if k in SACC or r.group in SACC: return "Saccharomycotina"
    if r.group == "Basidiomycota" or r.phylum == "Basidiomycota": return "Basidiomycota"
    if r.group == "Mucoromycota" or r.phylum == "Mucoromycota": return "Mucoromycota"
    if k in UNK: return "unassigned"
    return "Pezizomycotina"
g["rf_group"] = g.apply(rf_group, axis=1)
print(g.groupby("rf_group").size().sort_values(ascending=False).to_string())
tags = {t: i for i, t in enumerate(g.locustag)}
pf = {p: j for j, p in enumerate(sorted(c.pfam_acc.unique()))}
X = np.zeros((len(g), len(pf)))
for r in c[c.locustag.isin(tags)].itertuples():
    X[tags[r.locustag], pf[r.pfam_acc]] = r.n_domains
P = (X > 0).astype(float)
nfam = P.sum(1); ndom = X.sum(1)
g["n_pfam_families"] = nfam; g["n_domain_hits"] = ndom
print(len(g), "genomes,", len(pf), "Pfams; Alith strains:", alith.sum(), "families/genome", nfam[alith].round().tolist())
keep = P[~alith].mean(0) >= 0.02            # drop Pfams seen in <2 % of the other genomes
Xk, Pk = X[:, keep], P[:, keep]

def pcoa(D, k=5):
    n = len(D); J = np.eye(n) - 1.0 / n
    B = -0.5 * J @ (D ** 2) @ J
    w, v = np.linalg.eigh(B); o = np.argsort(w)[::-1]; w, v = w[o], v[:, o]
    pos = w[:k] > 0
    return v[:, :k][:, pos] * np.sqrt(w[:k][pos]), w[:k] / w[w > 0].sum()
rel = Xk / Xk.sum(1, keepdims=True)
D = {"braycurtis": squareform(pdist(rel, "braycurtis")), "jaccard": squareform(pdist(Pk > 0, "jaccard"))}
nb = []
for name, d in D.items():
    co, ev = pcoa(d)
    out = g[["label", "set", "group", "order", "n_pfam_families"]].copy()
    for i in range(co.shape[1]): out[f"PC{i+1}"] = co[:, i]
    out.to_csv(f"out/pcoa_{name}.tsv", sep="\t", index=False)
    print(name, "variance explained PC1-5:", (ev * 100).round(1).tolist())
    for i in np.where(alith)[0]:
        order = np.argsort(d[i]); order = [j for j in order if j != i and not alith[j]][:10]
        for rank, j in enumerate(order, 1):
            nb.append(dict(distance=name, strain=g.label[i], rank=rank, neighbour=g.label[j], group=g.group[j], order=g.order[j], d=round(d[i, j], 4)))
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 6.5))
        cmap = plt.get_cmap("tab20")
        for k, gr in enumerate(sorted(set(g.rf_group))):
            m = (g.rf_group == gr).values & ~alith
            if m.sum(): ax.scatter(co[m, 0], co[m, 1], s=14 if gr in ("Pezizomycotina", "Saccharomycotina", "Taphrinomycotina") else 9, label=f"{gr} ({m.sum()})", alpha=.75, color=cmap(k % 20))
        ax.scatter(co[alith, 0], co[alith, 1], marker="*", s=160, c="k", label="Antarctolithica strains")
        ax.set_xlabel(f"PC1 ({ev[0]*100:.1f} %)"); ax.set_ylabel(f"PC2 ({ev[1]*100:.1f} %)"); ax.set_title(f"Pfam composition, {name}")
        ax.legend(fontsize=6, ncol=2); fig.tight_layout(); fig.savefig(f"out/pcoa_{name}.png", dpi=150)
    except Exception as e: print("no plot:", e)
pd.DataFrame(nb).to_csv("out/neighbours.tsv", sep="\t", index=False)

# random forest
feat = np.hstack([Pk, np.log1p(Xk)])
lab = g.rf_group.values; train = ~alith & (lab != "unassigned")
grp_ok = {k for k, v in pd.Series(lab[train]).value_counts().items() if v >= 5}
train &= np.isin(lab, list(grp_ok))
orders = g.order.values
cv = []
for tr, te in GroupKFold(n_splits=8).split(feat[train], lab[train], orders[train]):
    idx = np.where(train)[0]
    m = RandomForestClassifier(500, random_state=1, n_jobs=-1, class_weight="balanced").fit(feat[idx[tr]], lab[idx[tr]])
    pr = m.predict(feat[idx[te]])
    for i, p in zip(idx[te], pr): cv.append(dict(label=g.label[i], true=lab[i], pred=p, n_pfam_families=nfam[i], order=orders[i]))
cv = pd.DataFrame(cv); cv["ok"] = cv.true == cv.pred; cv.to_csv("out/rf_cv.tsv", sep="\t", index=False)
print("RF leave-order-out accuracy by group:\n", cv.groupby("true").ok.agg(["mean", "size"]).round(3).to_string())
lo, hi = np.percentile(nfam[alith], [0, 100])
sm = cv[(cv.n_pfam_families >= lo * 0.85) & (cv.n_pfam_families <= hi * 1.15)]
print("genomes with an Alith-sized repertoire (%d-%d families): n=%d accuracy %.3f" % (lo * .85, hi * 1.15, len(sm), sm.ok.mean()))
print(sm.groupby("true").ok.agg(["mean", "size"]).round(3).to_string())
m = RandomForestClassifier(1000, random_state=1, n_jobs=-1, class_weight="balanced").fit(feat[train], lab[train])
pa = pd.DataFrame(m.predict_proba(feat[alith]), columns=m.classes_); pa.insert(0, "strain", g.label[alith].values)
pa.loc[len(pa)] = ["mean"] + pa.iloc[:, 1:].mean().tolist(); pa.round(3).to_csv("out/rf_alith.tsv", sep="\t", index=False)
print("Alith RF probabilities:\n", pa.round(3).to_string())

# Pfam signature test and core presence
rows = []
for gr in sorted(grp_ok):
    mem = (lab == gr) & ~alith
    for i in np.where(~alith)[0]:                       # leave-one-out definition so a member never defines its own signature
        mm = mem.copy(); mm[i] = False
        oth = ~mem & ~alith; oth[i] = False
        sig = (Pk[mm].mean(0) >= 0.8) & (Pk[oth].mean(0) <= 0.10)
        if sig.sum() < 20: continue
        rows.append(dict(group=gr, genome=g.label[i], member=bool(mem[i]), signature_size=int(sig.sum()), carried=Pk[i, sig].mean()))
    sig = (Pk[mem].mean(0) >= 0.8) & (Pk[~mem & ~alith].mean(0) <= 0.10)
    if sig.sum() >= 20:
        for i in np.where(alith)[0]:
            rows.append(dict(group=gr, genome=g.label[i], member="Alith", signature_size=int(sig.sum()), carried=Pk[i, sig].mean()))
sg = pd.DataFrame(rows); sg.to_csv("out/pfam_signature.tsv", sep="\t", index=False)
print("Pfam signature carried fraction (members / non-members / Alith):")
for gr, d in sg.groupby("group"):
    a = d[d.member == "Alith"].carried; mm = d[d.member == True].carried; nn = d[d.member == False].carried
    print(f"  {gr:22s} sig {d.signature_size.iloc[0]:4d} | members {mm.min():.2f}-{mm.max():.2f} (n={len(mm)}, median {mm.median():.2f}) | non-members {nn.min():.2f}-{nn.max():.2f} (median {nn.median():.2f}, above Alith: {(nn >= a.mean()).sum()}/{len(nn)}) | Alith {a.mean():.2f}")
core = P[~alith].mean(0) >= 0.95
cp = pd.DataFrame(dict(label=g.label, group=g.rf_group, n_pfam_families=nfam, core_size=int(core.sum()), core_carried=P[:, core].mean(1)))
cp.to_csv("out/core_presence.tsv", sep="\t", index=False)
print("core Pfams (>=95 %% of other genomes): %d; Alith carries %s; all other genomes: median %.3f, 1st percentile %.3f" %
      (core.sum(), cp.core_carried[alith].round(3).tolist(), cp.core_carried[~alith].median(), np.percentile(cp.core_carried[~alith], 1)))
small = ~alith & (nfam >= lo * 0.85) & (nfam <= hi * 1.15)
print("  genomes with similar repertoire size: n=%d, median core carried %.3f" % (small.sum(), cp.core_carried[small].median()))
