#!/usr/bin/env python3
"""Report on a V2 run: does the GFF -> protein space -> alignment projection behave on real annotation?

Usage: summarize_projection.py <results/projection>
It prints observations and *expectations from the literature* side by side; it asserts nothing about
biology. The hard checks (joins resolving, no length mismatches) are listed as counts you can read.
"""
import collections, csv, glob, os, sys

root = sys.argv[1]
labels = [l.split("\t")[0] for l in open("genomes.tsv") if l.strip() and not l.startswith("#")]

def tsv(path):
    with open(path) as f:
        return list(csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"))

def kv(path):
    with open(path) as f:
        return dict(l.rstrip("\n").split("\t", 1) for l in f if "\t" in l)

print("== timing (GFF parse + protein-space projection, one process per genome; wall s / peak RSS MB)")
for r in tsv(os.path.join(root, "per_genome_timing.tsv")):
    print(f"  {r['label']:14s} {float(r['wall_s']):8.2f}  {int(r['max_rss_kb'])/1024:8.1f}")

print("\n== per genome: annotation as read (all genes, before any alignment)")
print(f"  {'genome':14s} {'proteins joined':>16s} {'introns':>8s} {'introns/gene':>13s} {'phase0':>7s} {'phase1':>7s} {'phase2':>7s}")
dens = {}
for lab in labels:
    rows = tsv(os.path.join(root, "per_genome", lab, "introns.tsv"))
    genes = {r["protein_id"] for r in rows}
    ph = collections.Counter(r["phase"] for r in rows)
    n = len(rows); g = max(len(genes), 1)
    dens[lab] = n / g
    tot = max(n, 1)
    print(f"  {lab:14s} {len(genes):16d} {n:8d} {n/g:13.2f} {ph['0']/tot:7.2f} {ph['1']/tot:7.2f} {ph['2']/tot:7.2f}")
print("  (genes with introns only; intronless genes are absent from introns.tsv, so introns/gene is an upper bound)")

print("\n== diagnostics raised while reading each genome (kind: count)")
for lab in labels:
    p = os.path.join(root, "per_genome", lab, "diagnostics.tsv")
    d = ", ".join(f"{r['kind']}: {r['count']}" for r in tsv(p)) if os.path.exists(p) else "none"
    print(f"  {lab:14s} {d}")

allroot = os.path.join(root, "all")
print("\n== ortholog families through the whole pipeline")
alns = sorted(d for d in glob.glob(os.path.join(allroot, "*/")) if os.path.exists(os.path.join(d, "matrix.tsv")))
ann = collections.Counter(); reasons = collections.Counter(); per_taxon_sites = collections.Counter()
g = l = mp = sites = 0; branch = collections.defaultdict(lambda: [0, 0])
for d in alns:
    for r in tsv(os.path.join(d, "taxa.tsv")):
        ann[r["annotated"]] += 1
        if r["annotated"] != "yes": reasons[r["reason_if_not"]] += 1
        per_taxon_sites[r["taxon"]] += int(r["n_sites"])
    s = kv(os.path.join(d, "summary_score.tsv"))
    g += int(s.get("dollo_gains", 0)); l += int(s.get("dollo_losses", 0)); mp += int(s.get("mp_steps", 0)); sites += int(s.get("characters_scored", 0))
    for r in tsv(os.path.join(d, "dollo_branches.tsv")):
        b = branch[(r["branch"], r["tips"])]; b[0] += int(r["gains"]); b[1] += int(r["losses"])
print(f"  alignments with a matrix: {len(alns)}   protein rows annotated: {ann['yes']}   not annotated: {ann['no']}  {dict(reasons)}")
print(f"  sites scored {sites}   Dollo gains {g} losses {l}   MP steps {mp}")
print("  introns placed per taxon (alignment sites carried):")
for lab in labels:
    print(f"    {lab:14s} {per_taxon_sites[lab]:7d}")
dp = os.path.join(allroot, "diagnostics.tsv")
print("  diagnostics across the run:", ", ".join(f"{r['kind']}: {r['count']}" for r in tsv(dp)) if os.path.exists(dp) else "  none")

print("\n== per-branch Dollo losses and gains over all families (top 12 by losses)")
for (name, tips), (gn, ls) in sorted(branch.items(), key=lambda kv: -kv[1][1])[:12]:
    print(f"  {name:28s} losses {ls:6d} gains {gn:6d}   {tips if len(tips) < 60 else tips[:57] + '...'}")

print("""
== expectations from the literature (NOT asserted; compare by eye)
  - Intron-poor: Saccharomyces cerevisiae (few introns/gene). Intron-rich: Cryptococcus neoformans, Spizellomyces punctatus.
  - Phase 0 is the most common intron phase and phase 2 the rarest (Fedorov 1992; Roy & Gilbert 2005).
  - Dollo losses should concentrate on the S. cerevisiae and S. pombe branches and on the branch above the Saccharomycotina/Pezizomycotina split.
  - The gene models are funannotate predictions, not curated annotation; differences from literature values can be annotation artefacts.""")
