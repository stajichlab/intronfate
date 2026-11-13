#!/usr/bin/env python3
"""P1 pilot (a2): do extreme intron densities track assembly quality?

Reads intron_density_per_genome.tsv (output of p1_intron_density.py). Pure stdlib.
Reports: Spearman correlation of density with log N50 / contig count / CDS bp,
overall and within the classes with n>=50 (removes the clade confound);
outliers = density > class median + 3 * scaled MAD (n>=20 classes) or > 8/kb,
compared against non-outliers on N50, contigs, CDS bp.
"""
import csv, math, statistics as st, collections, sys

f = sys.argv[1] if len(sys.argv) > 1 else "intron_density_per_genome.tsv"
rows = []
for r in csv.DictReader(open(f), delimiter="\t"):
    try:
        r["d"] = float(r["introns_per_kb_cds"])
        r["n50"] = float(r["N50_bp"])
        r["ctg"] = float(r["contig_count"])
        r["cds"] = float(r["cds_bp"])
    except ValueError:
        continue
    if r["n50"] > 0 and r["ctg"] > 0:
        rows.append(r)
print("genomes with density + assembly stats:", len(rows))


def ranks(x):
    o = sorted(range(len(x)), key=x.__getitem__)
    rk = [0.0] * len(x)
    i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and x[o[j + 1]] == x[o[i]]:
            j += 1
        for k in range(i, j + 1):
            rk[o[k]] = (i + j) / 2 + 1
        i = j + 1
    return rk


def spearman(a, b):
    if len(a) < 10:
        return float("nan")
    ra, rb = ranks(a), ranks(b)
    ma, mb = st.fmean(ra), st.fmean(rb)
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    den = math.sqrt(sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb))
    return num / den if den else float("nan")


def corr_report(label, rs):
    d = [r["d"] for r in rs]
    print(f"{label:34s} n={len(rs):5d}  rho(dens,N50)={spearman(d,[r['n50'] for r in rs]):+.2f}"
          f"  rho(dens,contigs)={spearman(d,[r['ctg'] for r in rs]):+.2f}"
          f"  rho(dens,CDSbp)={spearman(d,[r['cds'] for r in rs]):+.2f}")


print("\n== Spearman correlations ==")
corr_report("ALL", rows)
byclass = collections.defaultdict(list)
for r in rows:
    byclass[r["CLASS"]].append(r)
for c, rs in sorted(byclass.items(), key=lambda kv: -len(kv[1])):
    if len(rs) >= 50:
        corr_report("class " + c[:28], rs)

print("\n== Outliers (class median + 3 scaled MAD, or >8/kb) ==")
out = []
for c, rs in byclass.items():
    d = [r["d"] for r in rs]
    med = st.median(d)
    mad = 1.4826 * st.median([abs(x - med) for x in d]) if len(d) >= 20 else None
    for r in rs:
        hi = r["d"] > 8 or (mad is not None and mad > 0 and r["d"] > med + 3 * mad)
        lo = mad is not None and mad > 0 and r["d"] < med - 3 * mad
        r["flag"] = "high" if hi else ("low" if lo else "")
        if r["flag"]:
            out.append(r)
for flag in ("high", "low"):
    sel = [r for r in rows if r.get("flag") == flag]
    rest = [r for r in rows if r.get("flag") == ""]
    print(f"\n{flag}-density outliers: {len(sel)} of {len(rows)}")
    if sel:
        for lab, key in (("N50 bp", "n50"), ("contigs", "ctg"), ("CDS bp", "cds")):
            print(f"  {lab:8s} median outliers={st.median([r[key] for r in sel]):>12,.0f}"
                  f"  non-outliers={st.median([r[key] for r in rest]):>12,.0f}")
        print("  by class:", collections.Counter(r["CLASS"] for r in sel).most_common(6))

# Quality bins: do density spread and outlier rate shrink with N50 within the large classes?
print("\n== Outlier rate and density spread by N50 bin (all genomes) ==")
bins = [(0, 5e4), (5e4, 2e5), (2e5, 1e6), (1e6, 5e6), (5e6, 1e12)]
for lo, hi in bins:
    sel = [r for r in rows if lo <= r["n50"] < hi]
    if not sel:
        continue
    nh = sum(r["flag"] == "high" for r in sel)
    nl = sum(r["flag"] == "low" for r in sel)
    d = sorted(r["d"] for r in sel)
    print(f"  N50 {lo:>9,.0f}-{min(hi,1e9):>11,.0f}: n={len(sel):5d} high={nh/len(sel):5.1%} low={nl/len(sel):5.1%}"
          f" dens p5/p50/p95={d[int(.05*(len(d)-1))]:.2f}/{d[len(d)//2]:.2f}/{d[int(.95*(len(d)-1))]:.2f}")

print("\n== Top 25 densities (genome, class, density, N50, contigs, CDS bp) ==")
for r in sorted(rows, key=lambda r: -r["d"])[:25]:
    print(f"  {r['LOCUSTAG']} {r['CLASS'][:18]:18s} {r['d']:6.2f} N50={r['n50']:>10,.0f} ctg={r['ctg']:>7,.0f} cds={r['cds']:>11,.0f} {r['SPECIES'][:30]}")
