#!/usr/bin/env python3
"""Variants of the supermatrix for robustness checks of the placement tests.
  recode     Dayhoff-6 recoding of amino acids to digits 0-5 (AGPST=0, DENQ=1, HKR=2, ILMV=3, FWY=4, C=5); gaps stay '-', ambiguous/X become '?'.
             Run with IQ-TREE `--st MORPH -m MK+G4` (compositional heterogeneity and saturation are reduced at the cost of information).
  drop-fast  remove the fastest-evolving fraction of sites using an IQ-TREE `--rate` file (per-site posterior mean rate on a fixed tree).
Usage: site_filter.py recode IN.fasta OUT.fasta
       site_filter.py drop-fast IN.fasta RATE_FILE FRACTION OUT.fasta"""
import sys
GROUPS = {"AGPST": "0", "DENQ": "1", "HKR": "2", "ILMV": "3", "FWY": "4", "C": "5"}
MAP = {a: d for g, d in GROUPS.items() for a in g}

def read_fasta(path):
    seqs, name = {}, None
    for l in open(path):
        l = l.rstrip("\n")
        if l.startswith(">"): name = l[1:].split()[0]; seqs[name] = []
        elif name: seqs[name].append(l)
    return {k: "".join(v) for k, v in seqs.items()}

def write_fasta(path, seqs):
    with open(path, "w") as f:
        for k, v in seqs.items(): f.write(f">{k}\n{v}\n")

def recode(seqs):
    return {k: "".join(MAP.get(c.upper(), "-" if c == "-" else "?") for c in v) for k, v in seqs.items()}

def read_rates(path):
    """Site -> rate from an IQ-TREE .rate file (comment lines start with '#', then a 'Site Rate Cat C_Rate' header)."""
    r = {}
    for l in open(path):
        if l.startswith("#") or l.startswith("Site") or not l.strip(): continue
        f = l.split(); r[int(f[0])] = float(f[1])
    return r

def drop_fast(seqs, rates, fraction):
    L = len(next(iter(seqs.values()))); assert sorted(rates) == list(range(1, L + 1)), "rate file does not match the alignment length"
    n_drop = int(round(fraction * L))
    order = sorted(range(1, L + 1), key=lambda i: (-rates[i], i))      # fastest first, ties by position
    dropped = set(order[:n_drop])
    return {k: "".join(c for i, c in enumerate(v, 1) if i not in dropped) for k, v in seqs.items()}, sorted(dropped)

if __name__ == "__main__":
    if sys.argv[1] == "recode":
        write_fasta(sys.argv[3], recode(read_fasta(sys.argv[2])))
    elif sys.argv[1] == "drop-fast":
        out, dropped = drop_fast(read_fasta(sys.argv[2]), read_rates(sys.argv[3]), float(sys.argv[4]))
        write_fasta(sys.argv[5], out); print(f"dropped {len(dropped)} fastest sites; {len(next(iter(out.values())))} remain", file=sys.stderr)
    else: sys.exit(__doc__)
