#!/usr/bin/env python3
"""One protein per gene (the longest isoform) from a GFF3 and its protein FASTA, for OrthoFinder.

NCBI RefSeq/GenBank and Ensembl genomes ship every isoform in the protein FASTA. OrthoFinder wants one protein per gene, otherwise
isoforms of one gene look like paralogs and fewer orthogroups are single-copy. This keeps, for each gene, the longest protein
(ties: the smallest protein id), so the families are built from one isoform per gene.

How proteins are tied to genes (GFF3):
  CDS row  --protein_id=-->  protein id (the FASTA id; NCBI 'XP_...', Ensembl stable id)
  CDS row  --Parent=----->   transcript (mRNA); transcript --Parent=--> gene
  A CDS whose Parent is not a transcript row is taken to hang directly off its gene.

Dropped (counted on stderr): proteins whose CDS rows carry no protein_id, proteins not in the FASTA, proteins in the FASTA that no CDS
names, and proteins with an internal stop '*' (an aligner would not keep them in the row, so intronfate would refuse that row anyway).
A terminal '*' is removed. Output ids are the FASTA ids (first word), unchanged, so `intronfate` can join the alignment rows back to the
annotation with the original GFF3 and protein files.

Usage: one_protein_per_gene.py --gff X.gff.gz --proteins X.protein.faa.gz --out proteomes/X.faa
"""
import argparse, gzip, os, sys
from urllib.parse import unquote


def opener(path):
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)


def attrs(col9):
    out = {}
    for kv in col9.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = unquote(v.strip())
    return out


def read_gff(path):
    """Return {protein_id: gene_id} from the CDS rows."""
    parent_of = {}          # feature id -> its first Parent (for transcripts)
    cds_parent = {}         # protein id -> CDS Parent
    with opener(path) as f:
        for line in f:
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            a = attrs(c[8])
            if c[2] == "CDS":
                pid = a.get("protein_id")
                par = a.get("Parent", "").split(",")[0]
                if pid and par:
                    cds_parent.setdefault(pid, par)
            elif "ID" in a and "Parent" in a:
                parent_of[a["ID"]] = a["Parent"].split(",")[0]
    return {pid: parent_of.get(par, par) for pid, par in cds_parent.items()}


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with opener(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0] if line[1:].split() else "", []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True, help="GFF3 (gzip is fine)")
    ap.add_argument("--proteins", required=True, help="protein FASTA with every isoform (gzip is fine)")
    ap.add_argument("--out", required=True, help="output FASTA, one protein per gene")
    a = ap.parse_args(argv)

    gene_of = read_gff(a.gff)
    seqs = read_fasta(a.proteins)
    best = {}
    n_stop = n_missing = 0
    for pid, gene in gene_of.items():
        s = seqs.get(pid)
        if s is None:
            n_missing += 1
            continue
        s = s.rstrip("*")
        if "*" in s:
            n_stop += 1
            continue
        cur = best.get(gene)
        if cur is None or len(s) > len(cur[1]) or (len(s) == len(cur[1]) and pid < cur[0]):
            best[gene] = (pid, s)
    n_orphan = len(set(seqs) - set(gene_of))

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    tmp = f"{a.out}.tmp{os.getpid()}"
    with open(tmp, "w") as out:
        for pid, s in sorted(best.values()):
            out.write(f">{pid}\n")
            for i in range(0, len(s), 60):
                out.write(s[i:i + 60] + "\n")
    os.replace(tmp, a.out)
    label = os.path.basename(a.out).split(".")[0]
    print(f"{label}: {len(best)} genes kept from {len(seqs)} proteins "
          f"({len(gene_of) - len(best) - n_stop - n_missing} other isoforms, {n_stop} internal stop, "
          f"{n_missing} named in the GFF but not in the FASTA, {n_orphan} in the FASTA but not in the GFF)", file=sys.stderr)


if __name__ == "__main__":
    main()
