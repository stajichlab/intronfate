#!/usr/bin/env python3
"""How much of the missing data is a short gap that column ambiguity could resolve?

For every `?` cell of a matrix made with --flank 1 from strict single-copy families (every taxon annotated, so `?` means
"not aligned around the site") find the gap run in that taxon's row at the site and classify it:
  short  : the run is at most --short columns and the taxon has residues on both sides within --reach columns
  long   : a longer gap run
  edge   : the row ends (or begins) within the window, no residues beyond
  other  : the window holds no gap (the cell is missing for another reason, e.g. a row mismatch)
Also reported: among the short ones, how many sites have at least one carrier, i.e. a presence whose absences nearby are in doubt.
Standard library only.
"""
import argparse, collections, os, sys


def read_fasta(path):
    rows, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    rows[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name is not None:
        rows[name] = "".join(buf)
    return rows


def classify(row, k, phase, flank, short):
    m = max(flank, 1)
    lo, hi = (k - flank, k + flank) if phase != 0 else (k - m, k + m - 1)
    lo0, hi0 = max(lo, 0), min(hi, len(row) - 1)
    cols = [i for i in range(lo0, hi0 + 1) if row[i] == "-"]
    if lo < 0 or hi >= len(row):
        return "edge", 0
    if not cols:
        return "other", 0
    a = cols[0]
    while a > 0 and row[a - 1] == "-":
        a -= 1
    b = cols[-1]
    while b < len(row) - 1 and row[b + 1] == "-":
        b += 1
    run = b - a + 1
    if a == 0 or b == len(row) - 1:
        return "edge", run
    return ("short" if run <= short else "long"), run


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--aln-dir", required=True)
    ap.add_argument("--aln-ext", default="aln.fasta")
    ap.add_argument("--taxon-map", required=True)
    ap.add_argument("--flank", type=int, default=1)
    ap.add_argument("--short", type=int, default=6)
    a = ap.parse_args()
    tmap = {}
    for line in open(a.taxon_map):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            tmap[p[0]] = p[1]
    cache, kinds, runs = {}, collections.Counter(), collections.Counter()
    sites_short = set()
    cells = 0
    with open(a.matrix) as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if header is None:
                header = p
                first = header.index("pattern") + 1
                taxa = header[first:]
                continue
            fam, col, phase = p[3], int(p[4]), int(p[2])
            cs = p[first:]
            fk = fam[:-4] if fam.endswith(".aln") else fam
            if fk not in cache:
                path = os.path.join(a.aln_dir, f"{fk}.{a.aln_ext}")
                cache[fk] = {tmap[k]: v for k, v in read_fasta(path).items() if k in tmap} if os.path.exists(path) else None
            rows = cache[fk]
            if rows is None:
                continue
            for i, c in enumerate(cs):
                cells += 1
                if c != "?":
                    continue
                row = rows.get(taxa[i])
                if row is None:
                    kinds["no_row"] += 1
                    continue
                kind, run = classify(row, col, phase, a.flank, a.short)
                kinds[kind] += 1
                if kind in ("short", "long"):
                    runs[min(run, 30)] += 1
                if kind == "short" and "1" in cs:
                    sites_short.add(p[0])
    miss = sum(kinds.values())
    print(f"cells {cells}  missing {miss} ({miss / cells:.3f})")
    for k in ("short", "long", "edge", "other", "no_row"):
        print(f"  {k:7s} {kinds[k]:7d}  {kinds[k] / max(miss, 1):.3f} of missing")
    print(f"sites with a carrier and at least one short-gap missing cell: {len(sites_short)}")
    print("gap-run length\tmissing cells")
    for r in sorted(runs):
        print(f"{r}\t{runs[r]}")


if __name__ == "__main__":
    main()
