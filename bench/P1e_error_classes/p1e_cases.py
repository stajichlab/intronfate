#!/usr/bin/env python3
"""P1e-b: (1) stratified test cases from genes.tsv (current annotation): up to N genes per (genome, class) with
coordinates, for use as a regression test set; (2) intron-length distribution of the reference vs BFD (current) per genome.
Run: p1e_cases.py [N]. Read-only."""
import csv, sys, os, random, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p1e_classify import parse_models
R = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/"
S = "/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/source/NCBI_ASM/"
N = int(sys.argv[1]) if len(sys.argv) > 1 else 25
rows = [r for r in csv.DictReader(open("genes.tsv"), delimiter="\t") if r["annotation"] == "cur" and r["class"] not in ("identical",)]
by = collections.defaultdict(list)
for r in rows: by[(r["asmid"], r["class"])].append(r)
rnd = random.Random(1)
with open("test_cases.tsv", "w") as fh:
    fh.write("asmid\tclass\tref_tid\tbfd_tid\tchrom\tstrand\tref_start\tref_end\tref_aa\tshared_aa\tmissing_aa\textra_aa\n")
    for (a, c), v in sorted(by.items()):
        for r in rnd.sample(v, min(N, len(v))):
            fh.write("\t".join(r[k] for k in ("asmid", "class", "ref_tid", "bfd_tid", "chrom", "strand", "ref_start", "ref_end", "ref_aa", "shared_aa", "missing_aa", "extra_aa")) + "\n")
outs = {r["ASMID"]: r["out"] for r in csv.DictReader(open("../P1d_version_audit/version_audit.tsv"), delimiter="\t")}
def q(v, p): v = sorted(v); return v[int(p * (len(v) - 1))]
print("intron length (bp): p50 / p90 / p99 / max, reference vs BFD current; n introns")
for a in sorted({r["asmid"] for r in csv.DictReader(open("genes.tsv"), delimiter="\t")}):
    o = outs[a]
    L = {}
    for name, p in (("ref", f"{S}{a}/{a}_genomic.gff.gz"), ("bfd", f"{R}genome_annotation/{o}/predict_results/{o}.gff3.gz")):
        if not os.path.exists(p): p = p[:-3]
        ms = parse_models(p)
        L[name] = sorted({(m["chrom"], i) for m in ms.values() for i in m["introns"]})
        L[name] = [e - s + 1 for _, (s, e) in L[name]]
    print(f"  {a[:26]:26s} ref {q(L['ref'],.5):4d}/{q(L['ref'],.9):4d}/{q(L['ref'],.99):5d}/{max(L['ref']):6d} n={len(L['ref']):6d} | BFD {q(L['bfd'],.5):4d}/{q(L['bfd'],.9):4d}/{q(L['bfd'],.99):5d}/{max(L['bfd']):6d} n={len(L['bfd']):6d}")
