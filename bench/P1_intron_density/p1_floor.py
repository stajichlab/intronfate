#!/usr/bin/env python3
"""P1 (b1): how many genomes would a CDS-bp floor remove, per class? Read-only on the P1 TSV.
Floor = fraction f of the CLASS median CDS bp (classes n>=20; smaller classes use the phylum median)."""
import csv, statistics as st, collections, sys
f = sys.argv[1]
rows = [r for r in csv.DictReader(open(f), delimiter="\t")]
for r in rows:
    r["cds"] = float(r["cds_bp"]); r["d"] = float(r["introns_per_kb_cds"])
byc = collections.defaultdict(list); byp = collections.defaultdict(list)
for r in rows:
    byc[r["CLASS"]].append(r["cds"]); byp[r["PHYLUM"]].append(r["cds"])
def ref(r):
    v = byc[r["CLASS"]]
    return st.median(v) if len(v) >= 20 else st.median(byp[r["PHYLUM"]])
for r in rows: r["ref"] = ref(r)
print("genomes", len(rows))
for fr in (0.25, 0.5, 0.75):
    rem = [r for r in rows if r["cds"] < fr * r["ref"]]
    print(f"\nfloor {fr:.0%} of class median CDS bp: removes {len(rem)} ({len(rem)/len(rows):.1%})")
    c = collections.Counter(r["CLASS"] or "(none)" for r in rem)
    tot = collections.Counter(r["CLASS"] or "(none)" for r in rows)
    for k, n in c.most_common(10):
        print(f"   {k[:26]:26s} {n:5d} of {tot[k]:5d} ({n/tot[k]:.1%})")
    keep = [r["d"] for r in rows if r["cds"] >= fr * r["ref"]]
    ks = sorted(keep); print("   kept density p1/p50/p99/max:", [round(ks[int(p*(len(ks)-1))],2) for p in (.01,.5,.99)], round(ks[-1],2))
    print("   kept with density>=8:", sum(x >= 8 for x in keep), " >=6:", sum(x >= 6 for x in keep))
# Saccharomycotina specifically
SC = {"Saccharomycetes","Pichiomycetes","Dipodascomycetes","Lipomycetes","Trigonopsidomycetes","Alaninomycetes","Saccharomycodomycetes","Phaffomycetes"}
sub = [r for r in rows if r["CLASS"] in SC]  # SUBPHYLUM column is empty in species.parquet
print("SUBPHYLUM values:", collections.Counter(r["SUBPHYLUM"] for r in rows).most_common(12))
print("\nSaccharomycotina-like classes n=", len(sub), collections.Counter(r["SUBPHYLUM"] for r in sub).most_common(3))
for fr in (0.25, 0.5, 0.75):
    print(f"  floor {fr:.0%}: removes", sum(r["cds"] < fr * r["ref"] for r in sub))
cd = sorted(r["cds"] for r in sub); print("  CDS bp p1/p5/p50:", [int(cd[int(p*(len(cd)-1))]) for p in (.01,.05,.5)])
