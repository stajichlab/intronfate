#!/usr/bin/env python3
"""Print, per flank, the cells masked (from the matrix headers), Dollo totals, and the absence-by-gap-distance table."""
import glob, os, re, sys

d = sys.argv[1]
print("flank\tsites\tmasked_present\tdollo_losses\tdollo_gains")
for f in sorted(glob.glob(os.path.join(d, "matrix_header.flank*.txt")), key=lambda x: int(re.findall(r"flank(\d+)", x)[0])):
    fl = re.findall(r"flank(\d+)", f)[0]
    h = open(f).read()
    sites = re.search(r"sites=(\d+)", h).group(1)
    mp = re.search(r"masked_present=(\d+)", h).group(1)
    with open(os.path.join(d, f"dollo_branches.flank{fl}.tsv")) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        gi, li = hdr.index("gains"), hdr.index("losses")
        g = l = 0
        for line in fh:
            p = line.rstrip("\n").split("\t")
            g += int(p[gi]); l += int(p[li])
    print(f"{fl}\t{sites}\t{mp}\t{l}\t{g}")
print()
print(open(os.path.join(d, "absence_by_gap_distance.tsv")).read())
