#!/usr/bin/env python3
"""J2b: J2 again, restricted to the introns that are cells of the strict-family matrix, with a relaxed match to RefSeq.

J2 (../J2_support_vs_reference) scored every BFD intron genome-wide and counted a RefSeq intron as agreeing only when chrom, first, last and
strand were identical. The likelihood only sees the introns that landed in a strict-family alignment column, and a splice site shifted by a few
bases is the same intron for it. Here:
  * the introns are the rows of the J1 work folder's per-family `sites.tsv` (the matrix cells), joined to `introns.tsv` on (protein, cds_offset)
    for their genomic coordinates;
  * a BFD intron agrees with RefSeq at level `exact` (identical), `near3` (RefSeq intron on the same sequence and strand with both ends within 3 nt:
    the same codon neighbourhood) or `near12` (both ends within 12 nt, four codons);
  * reads are those at the BFD coordinates (strand '.' counts for either), as in J2.
It still has no RefSeq projection into the alignment columns, so near3/near12 approximate "same column"; they are an upper bound on that.
Run: matrix_cells.py WORK RUN_DIR OUT_PREFIX [LABEL ...]   |   --selftest
"""
import collections, csv, glob, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "J2_support_vs_reference"))
sys.path.insert(0, os.path.join(HERE, "..", "P1f_adjudication"))
from support_vs_reference import BINS, R, S, auc, read_bed, reads_at  # noqa: E402

LEVELS = (("exact", 0), ("near3", 3), ("near12", 12))


def index_ref(rint):
    """{(chrom, strand): {first: [last, ...]}} for neighbour lookup"""
    idx = collections.defaultdict(lambda: collections.defaultdict(list))
    for c, s, e, st in rint:
        idx[(c, st)][s].append(e)
    return idx


def within(idx, key, tol):
    c, s, e, st = key
    d = idx.get((c, st))
    if not d:
        return False
    for ds in range(-tol, tol + 1):
        for le in d.get(s + ds, ()):
            if abs(le - e) <= tol:
                return True
    return False


def cell_introns(run_dir):
    """[(label, family, site, protein, cds_offset)] for every cell intron in the per-family sites.tsv files"""
    out = []
    for f in sorted(glob.glob(os.path.join(run_dir, "*", "sites.tsv"))):
        with open(f) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                out.append((r["taxon"], r["alignment"], f"{r['alignment']}:{r['column']}:{r['phase']}", r["protein_id"], int(r["cds_offset"])))
    return out


def coords(run_dir):
    """{(protein, cds_offset): (chrom, first, last, strand)} from the combined introns.tsv"""
    d = {}
    with open(os.path.join(run_dir, "introns.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["cds_offset"] not in ("", "NA", "."):
                d[(r["protein_id"], int(r["cds_offset"]))] = (r["seqid"], int(r["genomic_start"]), int(r["genomic_end"]), r["strand"])
    return d


def analyse(keys, rint, junc):
    """keys: BFD intron keys (a set). -> {level: (bins {label: [n, n_agree]}, auc)}"""
    idx = index_ref(rint)
    out = {}
    for name, tol in LEVELS:
        bins = {b[2]: [0, 0] for b in BINS}
        pos, neg = [], []
        for k in keys:
            n = reads_at(junc, k)
            lab = next(b[2] for b in BINS if b[0] <= n <= b[1])
            ok = (k in rint) if tol == 0 else within(idx, k, tol)
            bins[lab][0] += 1
            bins[lab][1] += ok
            (pos if ok else neg).append(n)
        out[name] = (bins, auc(pos, neg))
    return out


def selftest():
    rint = {("c", 100, 200, "+"), ("c", 1000, 1100, "-")}
    idx = index_ref(rint)
    assert within(idx, ("c", 100, 200, "+"), 0) and within(idx, ("c", 103, 197, "+"), 3)
    assert not within(idx, ("c", 104, 200, "+"), 3) and within(idx, ("c", 104, 200, "+"), 12)
    assert not within(idx, ("c", 100, 200, "-"), 12) and not within(idx, ("d", 100, 200, "+"), 12)
    junc = read_bed(["c\t99\t200\t.\t5\t+\n", "c\t299\t400\t.\t0\t+\n"])
    keys = {("c", 100, 200, "+"), ("c", 300, 400, "+"), ("c", 98, 201, "+")}
    res = analyse(keys, rint, junc)
    assert res["exact"][0]["3-4"] == [0, 0] and res["exact"][0]["5-9"] == [1, 1] and res["exact"][0]["0"][1] == 0
    assert res["near3"][0]["0"] == [2, 1], res["near3"][0]  # (98,201) is within 3 of RefSeq; (300,400) is not
    assert res["near12"][1] is not None
    print("selftest ok")


def main(work, run_dir, prefix, labels):
    from p1e_classify import parse_models
    outs = {r["out"]: r["ASMID"] for r in csv.DictReader(open(os.path.join(HERE, "..", "P1d_version_audit", "version_audit.tsv")), delimiter="\t")}
    genomes = {}
    for ln in open(os.path.join(HERE, "..", "V6_orthology_nine_species", "genomes.tsv")):
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.rstrip("\n").split("\t")
        genomes[f[0]] = f[1]
    xy, cells = coords(run_dir), cell_introns(run_dir)
    per = collections.defaultdict(set)
    miss = collections.Counter()
    for lab, fam, site, prot, off in cells:
        k = xy.get((prot, off))
        if k is None:
            miss[lab] += 1
        else:
            per[lab].add(k)
    fb, fd = open(prefix + "_bins.tsv", "w"), open(prefix + "_summary.tsv", "w")
    fb.write("label\tlevel\treads\tn_cell_introns\tn_agree\tfrac_agree\n")
    fd.write("label\tn_cell_introns\tunjoined\tzero_read\tfrac_zero_read\tauc_exact\tauc_near3\tauc_near12\tagree_exact\tagree_near3\tagree_near12\n")
    for lab in labels or sorted(per):
        d = genomes.get(lab)
        bed = os.path.join(work, lab + ".bed")
        if not d or d not in outs or not os.path.exists(bed) or lab not in per:
            print(lab, "skipped (no BED, assembly or cells)", file=sys.stderr)
            continue
        a = outs[d]
        ref = parse_models(f"{S}{a}/{a}_genomic.gff.gz")
        rint = {(m["chrom"], i[0], i[1], m["strand"]) for m in ref.values() for i in m["introns"]}
        junc = read_bed(open(bed))
        keys = per[lab]
        res = analyse(keys, rint, junc)
        zero = sum(1 for k in keys if reads_at(junc, k) == 0)
        for lev, (bins, _) in res.items():
            for b, (n, m) in bins.items():
                fb.write(f"{lab}\t{lev}\t{b}\t{n}\t{m}\t{(m / n if n else float('nan')):.4f}\n")
        agree = {lev: sum(v[1] for v in res[lev][0].values()) for lev in res}
        fd.write(f"{lab}\t{len(keys)}\t{miss[lab]}\t{zero}\t{zero / len(keys):.4f}\t" + "\t".join("NA" if res[l][1] is None else f"{res[l][1]:.4f}" for l in res) + "\t" + "\t".join(str(agree[l]) for l in res) + "\n")
        print(lab, len(keys), "unjoined", miss[lab], "agree", agree, flush=True)


if __name__ == "__main__":
    if sys.argv[1:2] == ["--selftest"]:
        selftest()
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:])
