#!/usr/bin/env python3
"""J2: does RNA-seq junction read support predict agreement with the RefSeq annotation?

For each genome with a J1 junction BED: every intron of the BFD (funannotate) annotation is labelled `in_ref` when RefSeq has an
intron with the same sequence, coordinates and strand (the same sets P1f uses). Its junction reads are the reads at exactly those
coordinates and strand (a `.` strand counts for either, as `intronfate junctions` does). Output, per genome:
  * reads bin x (n, fraction in RefSeq) for BFD introns,
  * AUC of read count for in_ref (the probability map p(n) is monotone in n, so its AUC is the same),
  * the fraction with >= 1 read for shared / bfd_only / ref_only introns (the adjudication, with this tool's reads).
RefSeq is not truth (V4), and the reads trained the BFD annotation (J1), so agreement here is a lower bar than independent support.

Run: support_vs_reference.py OUT_PREFIX [LABEL ...] | --selftest
"""
import collections, csv, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "P1f_adjudication"))
R = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/"
S = "/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/source/NCBI_ASM/"
BINS = [(0, 0, "0"), (1, 1, "1"), (2, 2, "2"), (3, 4, "3-4"), (5, 9, "5-9"), (10, 10**12, "10+")]


def read_bed(lines):
    """junction BED6 (start 0-based = first intron base - 1, score = reads) -> {(chrom, first, last, strand): reads}"""
    out = collections.Counter()
    for ln in lines:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 6 or ln.startswith(("#", "track")):
            continue
        out[(f[0], int(f[1]) + 1, int(f[2]), f[5])] += int(float(f[4]))
    return out


def reads_at(junc, key):
    c, s, e, st = key
    return junc.get(key, 0) + junc.get((c, s, e, "."), 0)


def auc(pos, neg):
    """P(score of a random positive > score of a random negative), ties half; None if a class is empty."""
    if not pos or not neg:
        return None
    allv = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg])
    ranks, i = {}, 0
    while i < len(allv):
        j = i
        while j < len(allv) and allv[j][0] == allv[i][0]:
            j += 1
        ranks[allv[i][0]] = (i + j + 1) / 2  # mean 1-based rank of the tie group
        i = j
    rp = sum(ranks[v] for v in pos)
    return (rp - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def analyse(bint, rint, junc):
    """-> (bins {label: [n, n_in_ref]}, auc, support {set: [n, n_ge1]})"""
    bins = {b[2]: [0, 0] for b in BINS}
    pos, neg = [], []
    for k in bint:
        n = reads_at(junc, k)
        lab = next(b[2] for b in BINS if b[0] <= n <= b[1])
        bins[lab][0] += 1
        if k in rint:
            bins[lab][1] += 1
            pos.append(n)
        else:
            neg.append(n)
    sup = {}
    for name, ks in (("shared", rint & bint), ("bfd_only", bint - rint), ("ref_only", rint - bint)):
        sup[name] = [len(ks), sum(1 for k in ks if reads_at(junc, k) >= 1)]
    return bins, auc(pos, neg), sup


def selftest():
    junc = read_bed(["c\t99\t200\t.\t12\t+\n", "c\t299\t400\t.\t0\t+\n", "c\t499\t600\t.\t3\t.\n"])
    assert junc[("c", 100, 200, "+")] == 12  # BED start 99 is the intron's first base 100
    assert reads_at(junc, ("c", 500, 600, "+")) == 3 and reads_at(junc, ("c", 500, 600, "-")) == 3  # '.' counts for either
    assert reads_at(junc, ("c", 100, 200, "-")) == 0  # wrong strand
    assert reads_at(junc, ("c", 101, 200, "+")) == 0  # shifted by one
    assert auc([5, 6], [1, 2]) == 1.0 and auc([1, 2], [5, 6]) == 0.0 and auc([3], [3]) == 0.5
    assert auc([1, 2, 3], []) is None
    b = {("c", 100, 200, "+"), ("c", 300, 400, "+"), ("c", 500, 600, "+")}
    r = {("c", 100, 200, "+"), ("c", 500, 600, "+"), ("c", 700, 800, "+")}
    bins, a, sup = analyse(b, r, junc)
    assert bins["10+"] == [1, 1] and bins["0"] == [1, 0] and bins["3-4"] == [1, 1], bins
    assert a == 1.0 - 0.0 or a == 1.0, a
    assert sup["shared"] == [2, 2] and sup["bfd_only"] == [1, 0] and sup["ref_only"] == [1, 0], sup
    print("selftest ok")


def main(prefix, labels):
    from p1e_classify import parse_models
    outs = {r["out"]: r["ASMID"] for r in csv.DictReader(open(os.path.join(HERE, "..", "P1d_version_audit", "version_audit.tsv")), delimiter="\t")}
    genomes = {}
    for ln in open(os.path.join(HERE, "..", "V6_orthology_nine_species", "genomes.tsv")):
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.rstrip("\n").split("\t")
        genomes[f[0]] = f[1]
    work = os.path.join(HERE, "..", "J1_junction_support", "work")
    fb, fs = open(prefix + "_bins.tsv", "w"), open(prefix + "_summary.tsv", "w")
    fb.write("label\treads\tn_bfd_introns\tn_in_refseq\tfrac_in_refseq\n")
    fs.write("label\tasmid\tn_bfd\tn_ref\tauc\tshared_n\tshared_ge1\tbfd_only_n\tbfd_only_ge1\tref_only_n\tref_only_ge1\n")
    for lab in labels or genomes:
        bed = os.path.join(work, lab + ".bed")
        d = genomes.get(lab)
        if not d or not os.path.exists(bed) or d not in outs:
            print(lab, "no BED or no RefSeq assembly; skipped", file=sys.stderr)
            continue
        a = outs[d]
        ref = parse_models(f"{S}{a}/{a}_genomic.gff.gz")
        bp = f"{R}genome_annotation/{d}/predict_results/{d}"
        bfd = parse_models(bp + ".gff3.gz" if os.path.exists(bp + ".gff3.gz") else bp + ".gff3")
        rint = {(m["chrom"], i[0], i[1], m["strand"]) for m in ref.values() for i in m["introns"]}
        bint = {(m["chrom"], i[0], i[1], m["strand"]) for m in bfd.values() for i in m["introns"]}
        bins, au, sup = analyse(bint, rint, read_bed(open(bed)))
        for k, (n, m) in bins.items():
            fb.write(f"{lab}\t{k}\t{n}\t{m}\t{(m / n if n else float('nan')):.4f}\n")
        fs.write(f"{lab}\t{a}\t{len(bint)}\t{len(rint)}\t{'NA' if au is None else f'{au:.4f}'}\t" + "\t".join(str(x) for s in ("shared", "bfd_only", "ref_only") for x in sup[s]) + "\n")
        print(lab, a, "auc", au, sup, flush=True)


if __name__ == "__main__":
    if sys.argv[1:2] == ["--selftest"]:
        selftest()
    else:
        main(sys.argv[1], sys.argv[2:])
