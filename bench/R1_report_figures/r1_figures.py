#!/usr/bin/env python3
"""R1: figures for the BFD intron report. Reads existing bench result tables (read-only); writes PNGs to ./figures."""
import os, sys
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

B = os.path.expanduser("~/projects/intron_projector/bench")
OUT = sys.argv[1] if len(sys.argv) > 1 else "figures"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "figure.dpi": 150})
C1, C2, C3 = "#2b6cb0", "#c05621", "#718096"

# Fig 1: density by class
d = pd.read_csv(f"{B}/P1_intron_density/intron_density_per_genome.tsv", sep="\t")
g = d.groupby("CLASS")["introns_per_kb_cds"]
keep = g.size()[lambda s: s >= 20].index
order = g.median()[keep].sort_values().index
fig, ax = plt.subplots(figsize=(7, 6))
ax.boxplot([d.loc[d.CLASS == c, "introns_per_kb_cds"] for c in order], vert=False, tick_labels=[f"{c} ({(d.CLASS==c).sum()})" for c in order],
           showfliers=False, whis=(5, 95), medianprops={"color": C2})
ax.axvline(d.introns_per_kb_cds.median(), color=C3, ls="--", lw=.8)
ax.set_xlabel("introns per kb CDS (box = IQR, whiskers = 5-95%)"); ax.set_title("Fig 1. Intron density by class (n >= 20; dashed = overall median)", fontsize=9)
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_density_by_class.png"); plt.close(fig)

# Fig 2: marker vs genome-wide
m = pd.read_csv(f"{B}/P1_intron_density/marker_density_per_genome.tsv", sep="\t").merge(d, on="LOCUSTAG")
m = m[m.n_marker_proteins >= 100]
fig, ax = plt.subplots(figsize=(5, 4.6))
ax.hexbin(m.introns_per_kb_cds, m.marker_introns_per_kb, gridsize=60, bins="log", cmap="Blues", mincnt=1, extent=(0, 6, 0, 6))
ax.plot([0, 6], [0, 6], color=C2, lw=.8)
rho = m[["introns_per_kb_cds", "marker_introns_per_kb"]].corr(method="spearman").iloc[0, 1]
ax.set_xlabel("genome-wide introns per kb CDS"); ax.set_ylabel("marker-gene introns per kb CDS")
ax.set_title(f"Fig 2. Marker vs genome-wide density\n(n={len(m):,} genomes with >=100 markers; Spearman {rho:.2f})")
fig.tight_layout(); fig.savefig(f"{OUT}/fig2_marker_vs_genome.png"); plt.close(fig)

# Fig 3: agreement with RefSeq, all genes vs markers
a = pd.read_csv(f"{B}/P1c_reference_introns/ref_vs_bfd.tsv", sep="\t")
mk = pd.read_csv(f"{B}/P1c_reference_introns/ref_vs_bfd_markers.tsv", sep="\t")
print(mk.columns.tolist())
a = a.assign(key=a.genome.str.split().str[0].str[:6]); mk = mk.assign(key=mk.genome.str.split().str[0].str[:6])
x = a.merge(mk[["key", "precision", "recall"]], on="key", suffixes=("", "_mk"))
x = x.sort_values("recall")
fig, ax = plt.subplots(figsize=(7, 3.8))
idx = np.arange(len(x)); w = .2
for k, (col, lab, col_c) in enumerate([("precision", "precision, all genes", C1), ("recall", "recall, all genes", C2), ("precision_mk", "precision, markers", "#90cdf4"), ("recall_mk", "recall, markers", "#fbd38d")]):
    ax.bar(idx + (k - 1.5) * w, x[col], w, label=lab, color=col_c)
ax.set_xticks(idx); ax.set_xticklabels(x.genome.str.replace(" ", "\n", 1), fontsize=7); ax.set_ylim(0, 1.3)
ax.legend(ncol=4, fontsize=7, loc="upper left"); ax.set_title("Fig 3. Agreement of BFD introns with RefSeq annotation (exact coordinates)")
fig.tight_layout(); fig.savefig(f"{OUT}/fig3_reference_agreement.png"); plt.close(fig)

# Fig 4: version audit (current vs pre_rc6 recall) and recall by class
v = pd.read_csv(f"{B}/P1d_version_audit/version_audit.tsv", sep="\t")
v = v[v.cur_ref_introns >= 200]
fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
both = v.dropna(subset=["pre_rec", "cur_rec"])
axs[0].scatter(both.pre_rec, both.cur_rec, s=5, alpha=.4, color=C1)
axs[0].plot([0, 1], [0, 1], color=C2, lw=.8)
axs[0].set_xlabel("recall, pre_rc6 prediction"); axs[0].set_ylabel("recall, current prediction")
axs[0].set_title(f"Recall vs RefSeq, same genome (n={len(both)})")
cs = v.groupby("CLASS").cur_rec; ks = cs.size()[lambda s: s >= 15].index; od = cs.median()[ks].sort_values().index
axs[1].boxplot([v.loc[v.CLASS == c, "cur_rec"] for c in od], vert=False, tick_labels=[f"{c} ({(v.CLASS==c).sum()})" for c in od], showfliers=False, medianprops={"color": C2})
axs[1].set_xlabel("recall of RefSeq introns (current)"); axs[1].set_title("Recall by class (n >= 15)")
fig.suptitle("Fig 4. Annotation-version audit against RefSeq"); fig.tight_layout(); fig.savefig(f"{OUT}/fig4_version_audit.png"); plt.close(fig)

# Fig 5: RNA-seq junction support (values from adjudication_summary, plan section 16)
gen = ["A. nidulans", "Z. tritici", "B05.10", "N. crassa", "C. neoformans", "S. pombe"]
bfd_only = [88, 81, 51, 57, 61, 38]; ref_only = [6, 6, 37, 39, 35, 48]
fig, ax = plt.subplots(figsize=(6.5, 3.6)); i = np.arange(len(gen))
ax.bar(i - .2, bfd_only, .4, label="BFD-only introns", color=C1); ax.bar(i + .2, ref_only, .4, label="reference-only introns", color=C2)
ax.set_xticks(i); ax.set_xticklabels(gen, rotation=20); ax.set_ylabel("% with exact Trinity junction")
ax.set_title("Fig 5. Transcript support of disputed introns (<200 bp for A.n., Z.t.)\n(same RNA-seq BFD used; see text)", fontsize=9); ax.legend(fontsize=7)
fig.tight_layout(); fig.savefig(f"{OUT}/fig5_adjudication.png"); plt.close(fig)

# Fig 6: reductions and re-expansions by threshold (reduced_rooted.txt, 'stay' / 'change' tie-breaks)
th = ["<0.3", "<0.5", "<1.0"]
red_lo = [13, 17, 26]; red_hi = [15, 17, 34]; re_lo = [7, 6, 6]; re_hi = [9, 6, 14]
fig, ax = plt.subplots(figsize=(5.5, 3.6)); i = np.arange(3)
ax.bar(i - .2, red_hi, .4, color=C1, label="reductions (late tie-break)"); ax.bar(i - .2, red_lo, .4, color="#1a365d", label="reductions (early tie-break)")
ax.bar(i + .2, re_hi, .4, color=C2, label="re-expansions (max)"); ax.bar(i + .2, re_lo, .4, color="#7b341e", label="re-expansions (min)")
ax.set_xticks(i); ax.set_xticklabels(th); ax.set_xlabel("reduced = introns per kb CDS"); ax.set_ylabel("events on 1,643-genus tree")
ax.set_title("Fig 6. Changes of intron-density state (rooted on Chytridiomycota)", fontsize=9); ax.legend(fontsize=6)
fig.tight_layout(); fig.savefig(f"{OUT}/fig6_state_changes.png"); plt.close(fig)
print("done", os.listdir(OUT))
