#!/usr/bin/env python3
"""Choose the V3 inputs from a size listing of every funannotate GFF3 in the BFD run.

Reads data/gff_sizes.tsv (bytes<TAB>genome-dir) and writes
  data/largest.tsv           the TOP_K largest GFF3s (one genome per line: label, gff, proteins, bytes)
  data/random_<N>.tsv        N genomes drawn at random with a fixed seed (for each N in G_STEPS)
Genomes whose protein FASTA is missing are skipped. Deterministic for a given listing and seed.
"""
import os, random, sys

bfd, top_k, steps, seed = sys.argv[1], int(sys.argv[2]), [int(x) for x in sys.argv[3].split()], int(sys.argv[4])
rows = []
for line in open("data/gff_sizes.tsv"):
    size, d = line.rstrip("\n").split("\t")
    gff = f"{bfd}/{d}/predict_results/{d}.gff3.gz"
    prot = f"{bfd}/{d}/predict_results/{d}.proteins.fa.gz"
    if os.path.exists(prot):
        rows.append((int(size), d, gff, prot))
rows.sort(key=lambda r: (-r[0], r[1]))
print(f"{len(rows)} genomes with GFF3 + proteins; GFF3 size MB: max {rows[0][0]/1e6:.1f}, median {rows[len(rows)//2][0]/1e6:.1f}, min {rows[-1][0]/1e6:.1f}")

def write(path, sel):
    with open(path, "w") as f:
        for size, d, gff, prot in sel:
            f.write(f"{d}\t{gff}\t{prot}\t{size}\n")

write("data/largest.tsv", rows[:top_k])
rng = random.Random(seed)
for n in steps:
    write(f"data/random_{n}.tsv", rng.sample(rows, min(n, len(rows))))
print("largest:")
for size, d, *_ in rows[:top_k]:
    print(f"  {size/1e6:8.1f} MB  {d}")
