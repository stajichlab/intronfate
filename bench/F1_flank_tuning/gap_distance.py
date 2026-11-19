#!/usr/bin/env python3
"""Absence rate of conserved introns as a function of the distance to the nearest alignment gap.

Input: a matrix.tsv made with --flank 0 (so every cell with a residue at the site is observable) and the alignments it came from.
A site is "conserved" when it is present in at least --min-others taxa other than the one being scored. For every observable
cell of a conserved site the distance d is the number of consecutive residue columns on the shorter side of the site in that
taxon's row (capped at --cap). Absence is `0`. If the absence rate is raised for small d and flat beyond some d0, the flank
should be d0: cells closer to a gap than that are not trustworthy absences.
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


def run_length(row, start, step, cap):
    n, i = 0, start
    while 0 <= i < len(row) and row[i] != "-" and n < cap:
        n += 1
        i += step
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--aln-dir", required=True)
    ap.add_argument("--aln-ext", default="aln.fasta")
    ap.add_argument("--taxon-map", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-others", type=int, default=2)
    ap.add_argument("--cap", type=int, default=12)
    a = ap.parse_args()

    tmap = {}
    with open(a.taxon_map) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                tmap[p[0]] = p[1]

    stats = collections.defaultdict(lambda: [0, 0])  # (phase, d) -> [cells, absent]
    by_taxon = collections.defaultdict(lambda: [0, 0])  # (taxon, distance bin) -> [cells, absent]
    examples = []
    cache = {}
    check = [0, 0]  # present phase-1/2 cells, and how many have a gap at the assumed column (must be ~0)
    with open(a.matrix) as fh:
        header = None
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if header is None:
                header = p
                first_taxon = header.index("pattern") + 1
                taxa = header[first_taxon:]
                continue
            fam, col, phase = p[3], int(p[4]), int(p[2])
            cells = p[first_taxon:]
            present = [i for i, c in enumerate(cells) if c == "1"]
            if len(present) < a.min_others:
                continue
            fam_key = fam[:-4] if fam.endswith(".aln") else fam
            if fam_key not in cache:
                path = os.path.join(a.aln_dir, f"{fam_key}.{a.aln_ext}")
                if not os.path.exists(path):
                    cache[fam_key] = None
                else:
                    rows = read_fasta(path)
                    cache[fam_key] = {tmap[k]: v for k, v in rows.items() if k in tmap}
            rows = cache[fam_key]
            if rows is None:
                continue
            for i, c in enumerate(cells):
                if c not in "01":
                    continue
                others = len(present) - (1 if c == "1" else 0)
                if others < a.min_others:
                    continue
                row = rows.get(taxa[i])
                if row is None:
                    continue
                # the matrix `column` is 0-based, like the row index
                k = col
                if c == "1" and phase != 0:
                    check[0] += 1
                    bad = k >= len(row) or row[k] == "-"
                    check[1] += bad
                    if bad and len(examples) < 8:
                        examples.append((fam, col, phase, taxa[i], len(row)))
                if phase == 0:
                    left, right = run_length(row, k - 1, -1, a.cap), run_length(row, k, 1, a.cap)
                else:
                    left, right = run_length(row, k - 1, -1, a.cap), run_length(row, k + 1, 1, a.cap)
                d = min(left, right)
                tb = by_taxon[(taxa[i], "0" if d == 0 else "1" if d == 1 else "2-3" if d < 4 else "4-7" if d < 8 else "8+")]
                tb[0] += 1
                tb[1] += c == "0"
                s = stats[(phase, d)]
                s[0] += 1
                s[1] += c == "0"
    sys.stderr.write(f"index check: {check[1]} of {check[0]} present phase-1/2 cells sit on a gap\n")
    for e in examples:
        sys.stderr.write(f"  present on a gap: family={e[0]} column={e[1]} phase={e[2]} taxon={e[3]} row_length={e[4]}\n")
    with open(a.out + ".by_taxon", "w") as bt:
        bt.write("taxon\tgap_distance\tcells\tabsent\tabsence_rate\n")
        for (tx, b) in sorted(by_taxon):
            n, z = by_taxon[(tx, b)]
            bt.write(f"{tx}\t{b}\t{n}\t{z}\t{z / n:.4f}\n")
    with open(a.out, "w") as out:
        out.write("phase\tdistance\tcells\tabsent\tabsence_rate\n")
        for (phase, d) in sorted(stats):
            n, z = stats[(phase, d)]
            out.write(f"{phase}\t{d}\t{n}\t{z}\t{z / n:.4f}\n")
        tot = collections.defaultdict(lambda: [0, 0])
        for (phase, d), (n, z) in stats.items():
            tot[d][0] += n
            tot[d][1] += z
        for d in sorted(tot):
            n, z = tot[d]
            out.write(f"all\t{d}\t{n}\t{z}\t{z / n:.4f}\n")


if __name__ == "__main__":
    main()
