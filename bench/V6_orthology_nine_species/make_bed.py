#!/usr/bin/env python3
"""funannotate GFF3 -> BED for orthoSynAssign: `scaffold<TAB>start0<TAB>end<TAB>protein ids` (isoforms joined by ';').

funannotate writes `gene` rows (ID=<locus>) and `mRNA` rows (ID=<locus>-T1, Parent=<locus>); the protein FASTA id is the mRNA ID.
Usage: make_bed.py <annotation.gff3[.gz]> > species.bed
"""
import gzip, sys
from collections import OrderedDict

def attrs(s):
    out = {}
    for kv in s.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1); out[k] = v
    return out

genes = {}                       # gene id -> (seqid, start0, end)
proteins = OrderedDict()         # gene id -> [mRNA ids]
op = gzip.open if sys.argv[1].endswith(".gz") else open
with op(sys.argv[1], "rt") as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9:
            continue
        a = attrs(c[8])
        if c[2] == "gene" and "ID" in a:
            genes[a["ID"]] = (c[0], int(c[3]) - 1, int(c[4]))
        elif c[2] == "mRNA" and "ID" in a and "Parent" in a:
            proteins.setdefault(a["Parent"], []).append(a["ID"])
rows = [(*genes[g], ";".join(p)) for g, p in proteins.items() if g in genes]
rows.sort(key=lambda r: (r[0], r[1], r[2]))
for r in rows:
    print("\t".join(map(str, r)))
print(f"{len(rows)} protein-coding genes written; {sum(1 for g in proteins if g not in genes)} mRNA parents without a gene row", file=sys.stderr)
