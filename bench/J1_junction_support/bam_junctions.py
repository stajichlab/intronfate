#!/usr/bin/env python3
"""Splice junctions and their read counts from SAM text on stdin (`samtools view file.bam | bam_junctions.py > x.bed`).

Output is BED6 for `intronfate junctions --junction-format bed`: chrom, start (0-based), end, name, score = reads, strand.
A junction is an `N` run of the CIGAR. Its first intron base is the reference position after the preceding M/D/=/X operations
(1-based), its last is that plus N minus 1 -- the same convention as a GFF intron, so BED start = first base - 1.
Strand comes from the aligner's `XS:A:` tag, `.` when absent. Unmapped, secondary and supplementary records are skipped; each
remaining record that crosses a junction counts once (the two mates of a pair that both cross it count twice).
Standard library only.
"""
import argparse, re, sys
from collections import Counter

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def junctions(lines, min_intron=20):
    counts = Counter()
    for line in lines:
        if line.startswith("@"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 11:
            continue
        flag, chrom, pos, cigar = int(f[1]), f[2], int(f[3]), f[5]
        if flag & (0x4 | 0x100 | 0x800) or cigar == "*" or chrom == "*":
            continue
        strand = "."
        for tag in f[11:]:
            if tag.startswith("XS:A:"):
                strand = tag[5:]
        ref = pos
        for n, op in CIGAR.findall(cigar):
            n = int(n)
            if op == "N":
                if n >= min_intron:
                    counts[(chrom, ref, ref + n - 1, strand)] += 1
                ref += n
            elif op in "MD=X":
                ref += n
    return counts


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--min-intron", type=int, default=20)
    ap.add_argument("--min-reads", type=int, default=1)
    a = ap.parse_args()
    for (chrom, start, end, strand), n in sorted(junctions(sys.stdin, a.min_intron).items()):
        if n >= a.min_reads:
            print(f"{chrom}\t{start - 1}\t{end}\tj\t{n}\t{strand}")


if __name__ == "__main__":
    main()
