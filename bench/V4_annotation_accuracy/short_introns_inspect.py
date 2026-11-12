#!/usr/bin/env python3
"""List reference introns <50 bp (and predicted novel ones <50 bp) of one genome from the per-intron scoring table, with the nearest
predicted intron on the same transcript pair (offset), to see whether they are real, or an artefact of coordinates / CDS pieces.
Usage: short_introns_inspect.py SCORING DIR [MAXLEN=50]"""
import csv, gzip, sys, collections
scoring, d = sys.argv[1], sys.argv[2]; M = int(sys.argv[3]) if len(sys.argv) > 3 else 50
rows = list(csv.DictReader(gzip.open(f"{scoring}/reference_introns/{d}.new.introns.tsv.gz", "rt"), delimiter="\t"))
print("columns:", list(rows[0]))
L = lambda r: int(r["end"]) - int(r["start"]) + 1
pred = collections.defaultdict(list)
for r in rows:
    if r["source"] == "pred": pred[(r["seqid"], r["strand"])].append(r)
seen = set()
print("== reference introns <%d, status, nearest pred intron (same seqid/strand, by start)" % M)
for r in rows:
    if r["source"] == "ref" and L(r) < M:
        k = (r["seqid"], r["start"], r["end"], r["strand"])
        if k in seen: continue
        seen.add(k)
        near = min(pred[(r["seqid"], r["strand"])], key=lambda p: abs(int(p["start"]) - int(r["start"])), default=None)
        n = f'{near["start"]}-{near["end"]} len{L(near)} {near["status"]}' if near and abs(int(near["start"]) - int(r["start"])) < 3000 else "none within 3kb"
        print(r["seqid"], r["start"], r["end"], r["strand"], "len", L(r), r["status"], r["transcript"], f'{r["intron_number"]}/{r["n_introns"]}', "| near:", n)
print("== novel predicted introns <%d" % M)
for r in rows:
    if r["source"] == "pred" and r["status"] != "matched" and L(r) < M:
        print(r["seqid"], r["start"], r["end"], r["strand"], "len", L(r), r["transcript"], f'{r["intron_number"]}/{r["n_introns"]}')
