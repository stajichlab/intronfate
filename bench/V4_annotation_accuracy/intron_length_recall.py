#!/usr/bin/env python3
"""Is intron recall length-dependent? For each genome with a per-intron table (colleague's gffcompare scoring), bin the distinct
reference introns by length and report the fraction matched (Sn) per bin; also the length of predicted introns that are novel.
Read-only on SCORING.  Usage: intron_length_recall.py SCORING OUT.tsv LABEL=DIR [LABEL=DIR ...]
Intron length = end - start + 1 (reference introns are distinct (seqid,start,end,strand) among source=ref rows)."""
import csv, gzip, sys, collections
scoring, out = sys.argv[1], sys.argv[2]
BINS = [(0, 50), (50, 100), (100, 200), (200, 400), (400, 800), (800, 10**9)]
def b(L):
    for lo, hi in BINS:
        if lo <= L < hi: return f"{lo}-{hi if hi < 10**9 else 'inf'}"
rows = []
for spec in sys.argv[3:]:
    label, d = spec.split("=")
    ref = {}; novel = collections.Counter(); pm = collections.Counter()
    with gzip.open(f"{scoring}/reference_introns/{d}.new.introns.tsv.gz", "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            L = int(r["end"]) - int(r["start"]) + 1
            if r["source"] == "ref":
                ref[(r["seqid"], r["start"], r["end"], r["strand"])] = (r["status"], L)
            else:
                (pm if r["status"] == "matched" else novel)[b(L)] += 1
    tot = collections.Counter(); hit = collections.Counter()
    for st, L in ref.values():
        tot[b(L)] += 1; hit[b(L)] += st == "matched"
    for lo, hi in BINS:
        k = f"{lo}-{hi if hi < 10**9 else 'inf'}"
        rows.append(dict(genome=label, bin=k, ref=tot[k], matched=hit[k], sn=round(100 * hit[k] / tot[k], 1) if tot[k] else "",
                         pred_matched=pm[k], pred_novel=novel[k]))
    Ls = sorted(L for _, L in ref.values())
    print(label, "ref introns", len(Ls), "median", Ls[len(Ls)//2] if Ls else None, "max", Ls[-1] if Ls else None, file=sys.stderr)
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
print(open(out).read())
