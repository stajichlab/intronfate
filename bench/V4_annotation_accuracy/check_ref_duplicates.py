#!/usr/bin/env python3
"""Audit of the per-intron reference table behind the Sn figures: can one reference intron (seqid,start,end,strand) appear in several rows with different
`status`? If so, "last row wins" (as busco_gene_accuracy.py did) makes Sn depend on row order. Usage: check_ref_duplicates.py SCORING genome_dir [genome_dir ...]"""
import collections, csv, gzip, sys
scoring = sys.argv[1]
for d in sys.argv[2:]:
    rows = collections.defaultdict(list)
    with gzip.open(f"{scoring}/reference_introns/{d}.new.introns.tsv.gz", "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["source"] != "pred": rows[(r["seqid"], r["start"], r["end"], r["strand"])].append((r["status"], r["transcript"]))
    n = len(rows); multi = sum(1 for v in rows.values() if len(v) > 1)
    conflict = sum(1 for v in rows.values() if len({s for s, _ in v}) > 1)
    last = sum(1 for v in rows.values() if v[-1][0] == "matched"); anyi = sum(1 for v in rows.values() if any(s == "matched" for s, _ in v))
    first = sum(1 for v in rows.values() if v[0][0] == "matched")
    print(f"{d}: {n} distinct ref introns; in >1 row: {multi}; with conflicting status: {conflict}; matched by last-row {last} ({100*last/n:.1f}%), first-row {first} ({100*first/n:.1f}%), any-row {anyi} ({100*anyi/n:.1f}%)")
