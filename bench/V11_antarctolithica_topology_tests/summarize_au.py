#!/usr/bin/env python3
"""Table of an IQ-TREE -z/-au run: tree name (from all_trees.names, same order), logL, deltaL, bp-RELL, p-KH, p-SH, c-ELW, p-AU.
Usage: summarize_au.py au.iqtree all_trees.names [...more pairs]"""
import re, sys
def table(iq, names):
    rows, inside = [], False
    for l in open(iq):
        if l.startswith("Tree      logL"): inside = True; continue
        if inside and l.startswith("---"): continue
        if inside and re.match(r"\s*\d+ ", l):
            f = [x for x in l.split() if x not in ("+", "-")]
            rows.append(f)
        elif inside and not l.strip(): break
    nm = [x.strip() for x in open(names) if x.strip()]
    assert len(nm) == len(rows), (len(nm), len(rows))
    return [(n, r[1], r[2], r[3], r[4], r[5], r[8] if len(r) > 8 else "", r[-1]) for n, r in zip(nm, rows)]
if __name__ == "__main__":
    a = sys.argv[1:]
    for iq, names in zip(a[0::2], a[1::2]):
        print("==", iq); print("tree\tlogL\tdeltaL\tbpRELL\tpKH\tpSH\tcELW\tpAU")
        for r in table(iq, names): print("\t".join(r))
