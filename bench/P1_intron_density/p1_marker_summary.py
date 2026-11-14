#!/usr/bin/env python3
"""Compare marker-gene intron density with genome-wide density (P1 b2). Stdlib only."""
import csv, math, statistics as st, collections
tsv = {r["LOCUSTAG"]: r for r in csv.DictReader(open("intron_density_per_genome.tsv"), delimiter="\t")}
mk = {r["LOCUSTAG"]: r for r in csv.DictReader(open("marker_density_per_genome.tsv"), delimiter="\t")}
rows = []
for t, m in mk.items():
    if t in tsv and m["marker_introns_per_kb"] != "nan":
        r = tsv[t]
        rows.append(dict(t=t, cls=r["CLASS"], sp=r["SPECIES"], g=float(r["introns_per_kb_cds"]), m=float(m["marker_introns_per_kb"]),
                         n=int(m["n_marker_proteins"]), mbp=int(m["marker_cds_bp"]), cds=float(r["cds_bp"]), mi=int(m["marker_introns"])))
print("genomes", len(rows), " markers available (Pfams):", sum(1 for _ in open("marker_pfams.tsv")) - 1)
def ranks(x):
    o = sorted(range(len(x)), key=x.__getitem__); rk = [0.0]*len(x); i = 0
    while i < len(o):
        j = i
        while j+1 < len(o) and x[o[j+1]] == x[o[i]]: j += 1
        for k in range(i, j+1): rk[o[k]] = (i+j)/2+1
        i = j+1
    return rk
def sp(a, b):
    ra, rb = ranks(a), ranks(b); ma, mb = st.fmean(ra), st.fmean(rb)
    return sum((x-ma)*(y-mb) for x, y in zip(ra, rb))/math.sqrt(sum((x-ma)**2 for x in ra)*sum((y-mb)**2 for y in rb))
def q(v, ps): s = sorted(v); return [round(s[int(p*(len(s)-1))], 2) for p in ps]
print("markers found per genome p1/p5/p50/p95:", q([r["n"] for r in rows], (.01,.05,.5,.95)))
print("Spearman marker vs genome-wide density: %.3f" % sp([r["g"] for r in rows], [r["m"] for r in rows]))
print("marker density quantiles 1/5/50/95/99:", q([r["m"] for r in rows], (.01,.05,.5,.95,.99)))
print("genome-wide  quantiles 1/5/50/95/99:   ", q([r["g"] for r in rows], (.01,.05,.5,.95,.99)))
print("ratio marker/genome-wide median:", round(st.median([r["m"] for r in rows])/st.median([r["g"] for r in rows]), 2))
ok = [r for r in rows if r["n"] >= 50]
print("genomes with >=50 marker proteins:", len(ok), " <50:", len(rows)-len(ok), " <100:", sum(r["n"] < 100 for r in rows))
print("\ngenomes at genome-wide >=6/kb:", sum(r["g"] >= 6 for r in rows), "; of these marker density p50:", q([r["m"] for r in rows if r["g"] >= 6], (.5,)) if any(r["g"] >= 6 for r in rows) else "-")
print("genomes at marker >=6/kb (n>=50):", sum(r["m"] >= 6 for r in ok), "; >=8:", sum(r["m"] >= 8 for r in ok))
print("\nper class (n>=40): median genome-wide vs marker density, median #markers")
byc = collections.defaultdict(list)
for r in rows: byc[r["cls"]].append(r)
out = []
for c, v in byc.items():
    if len(v) >= 40:
        out.append((st.median([r["m"] for r in v]), c, len(v), st.median([r["g"] for r in v]), st.median([r["n"] for r in v]), sp([r["g"] for r in v], [r["m"] for r in v])))
for m, c, n, g, nm, rho in sorted(out)[:12] + sorted(out)[-4:]:
    print(f"  {c[:24]:24s} n={n:5d} genome-wide={g:5.2f} marker={m:5.2f} markers/genome={nm:5.0f} rho={rho:+.2f}")
print("\nclass-level Spearman of medians (n>=40 classes):", round(sp([o[3] for o in out], [o[0] for o in out]), 3))
print("\nsuspicious genomes (genome-wide >= 8 per kb or marker >= 8 with n>=50):")
for r in sorted(rows, key=lambda r: -max(r["g"], r["m"] if r["n"] >= 50 else 0))[:15]:
    print(f"  {r['t']} {r['cls'][:16]:16s} genome-wide={r['g']:5.2f} marker={r['m']:5.2f} nmark={r['n']:4d} {r['sp'][:28]}")
