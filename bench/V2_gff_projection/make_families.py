#!/usr/bin/env python3
"""Turn per-genome BUSCO tables into ortholog families.

A BUSCO is used when it is `Complete` (single copy) in at least MIN_GENOMES genomes. For each, write the
unaligned proteins to families/<busco_id>.faa and add `protein_id<TAB>label` lines to taxon_map.tsv.
Usage: make_families.py <min_genomes>
"""
import collections, csv, gzip, os, sys

min_genomes = int(sys.argv[1]) if len(sys.argv) > 1 else 7
labels = [l.split("\t")[0] for l in open("genomes.tsv") if l.strip() and not l.startswith("#")]

def read_fasta(path):
    seqs, name, buf = {}, None, []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name: seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs

fam = collections.defaultdict(dict)          # busco id -> {label: protein id}
for lab in labels:
    with open(f"results/busco/{lab}/full_table.tsv") as f:
        for row in csv.reader((l for l in f if not l.startswith("#")), delimiter="\t"):
            if len(row) >= 3 and row[1] == "Complete":
                fam[row[0]][lab] = row[2]
keep = {b: m for b, m in fam.items() if len(m) >= min_genomes}
os.makedirs("families", exist_ok=True)
prot = {lab: read_fasta(f"data/{lab}/{lab}.proteins.fa.gz") for lab in labels}
# Funannotate protein ids are transcript ids (e.g. LOCUSTAG_000123-T1); BUSCO may report them with or
# without the -T1 suffix, so resolve each id against the FASTA and drop anything that does not resolve.
def resolve(lab, pid):
    for cand in (pid, pid + "-T1"):
        if cand in prot[lab]:
            return cand
    return None
n_fam = n_rows = unresolved = internal_stop = 0
with open("taxon_map.tsv", "w") as tm:
    for b, members in sorted(keep.items()):
        rows = []
        for lab, pid in sorted(members.items()):
            r = resolve(lab, pid)
            if r is None: unresolved += 1; continue
            # A premature stop codon (a pseudogene-like model) is dropped by aligners, so the aligned row
            # would no longer equal the protein and the tool would (rightly) refuse to project through it.
            if "*" in prot[lab][r].rstrip("*"): internal_stop += 1; continue
            rows.append((lab, r))
        if len(rows) < min_genomes: continue
        with open(f"families/{b}.faa", "w") as out:
            for lab, r in rows:
                out.write(f">{r}\n{prot[lab][r].rstrip('*')}\n")
                tm.write(f"{r}\t{lab}\n")
        n_fam += 1; n_rows += len(rows)
dist = collections.Counter(len(m) for m in keep.values())
print(f"BUSCOs single-copy in >= {min_genomes} genomes: {len(keep)}; families written: {n_fam}; sequences: {n_rows}; unresolved ids: {unresolved}; dropped for an internal stop codon: {internal_stop}")
print("genomes per family:", dict(sorted(dist.items())))
