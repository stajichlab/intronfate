import csv, statistics as st, collections
f = "/rhome/jstajich/projects/intron_projector/bench/P1_intron_density/intron_density_per_genome.tsv"
rows = list(csv.DictReader(open(f), delimiter="\t"))
d = [float(r["introns_per_kb_cds"]) for r in rows]
d_sorted = sorted(d)
q = lambda p: d_sorted[int(p * (len(d_sorted) - 1))]
print("genomes", len(rows))
print("density per kb CDS quantiles 1/5/25/50/75/95/99%:", [round(q(p), 2) for p in (.01, .05, .25, .5, .75, .95, .99)])
print("genomes <0.5:", sum(x < 0.5 for x in d), " <1:", sum(x < 1 for x in d), " <2:", sum(x < 2 for x in d), " >=4:", sum(x >= 4 for x in d))
print("annotation_status:", collections.Counter(r["annotation_status"] for r in rows).most_common(5))

def group(key, minn):
    g = collections.defaultdict(list)
    for r in rows:
        g[(r["PHYLUM"], r[key])].append(float(r["introns_per_kb_cds"]))
    return {k: v for k, v in g.items() if len(v) >= minn}

for key, minn in (("CLASS", 20), ("ORDER", 10)):
    g = group(key, minn)
    meds = sorted(((st.median(v), k, len(v), min(v), max(v)) for k, v in g.items()))
    print(f"\n{key}: {len(g)} groups with n>={minn}; lowest 12 by median, then highest 5")
    for m, k, n, lo, hi in meds[:12] + meds[-5:]:
        print(f"  {k[0][:14]:14s} {k[1][:28]:28s} n={n:5d} median={m:5.2f} range={lo:5.2f}-{hi:5.2f}")
    print(f"  groups with median <1.0: {sum(m < 1.0 for m, *_ in meds)} of {len(meds)}; <2.0: {sum(m < 2.0 for m, *_ in meds)}")

# within-group spread: orders where genomes differ a lot (candidate recent shifts)
g = group("ORDER", 10)
spread = sorted(((max(v) - min(v), k, len(v), min(v), max(v)) for k, v in g.items()), reverse=True)[:8]
print("\nOrders with largest within-order density range (n>=10):")
for s, k, n, lo, hi in spread:
    print(f"  {k[0][:14]:14s} {k[1][:28]:28s} n={n:5d} range={lo:5.2f}-{hi:5.2f}")
