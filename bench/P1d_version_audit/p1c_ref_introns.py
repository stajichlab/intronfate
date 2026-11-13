#!/usr/bin/env python3
"""P1c: compare BFD intron calls with RefSeq reference annotations for model genomes.

Reference introns = gaps between consecutive CDS features of the same RNA (GFF Parent), unique per
(chrom,start,end,strand). BFD introns from gene_introns (LOCUSTAG = first 8 chars of intron_id).
Reports per genome: counts, exact-site precision/recall, recall within BFD-gene-overlapping reference
genes (separates missing genes from wrong structure), boundary-shift matches, non-canonical rates.
Read-only on all inputs. Run: p1c_ref_introns.py OUTDIR   |   p1c_ref_introns.py --selftest
"""
import sys, gzip, collections, re

SRC = "/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/source/NCBI_ASM/"
T = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/"
REFS = {  # LOCUSTAG: (label, source dir)
    "F197696E": ("A.nidulans FGSC A4", "GCF_000011425.1_ASM1142v1"),
    "F416552B": ("C.neoformans H99", "GCF_000149245.1_CNA3"),
    "FC69C3D3": ("N.crassa OR74A", "GCF_000182925.2_NC12"),
    "F24A1226": ("S.cerevisiae S288C", "GCF_000146045.2_R64"),
    "FD3AB542": ("S.pombe 972h-", "GCF_000002945.2_ASM294v3"),
    "F54A7E8F": ("Z.tritici IPO323", "GCF_000219625.1_MYCGR_v2.0"),
    "F9A9046A": ("B.cinerea B05.10", "GCF_000143535.2_ASM14353v4"),
    "FECD96D2": ("B.cinerea BcDW1", "GCA_000349525.1_BcDW1"),
}


def ref_introns(cds_by_rna):
    """cds_by_rna: {rna: [(chrom,start,end,strand), ...]} -> set of (chrom,start,end,strand) introns."""
    out = set()
    for rna, segs in cds_by_rna.items():
        segs = sorted(segs, key=lambda s: s[1])
        for a, b in zip(segs, segs[1:]):
            if b[1] - a[2] > 1:
                out.add((a[0], a[2] + 1, b[1] - 1, a[3]))
    return out


def parse_gff(path):
    cds = collections.defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            m = re.search(r"Parent=([^;]+)", f[8])
            if m:
                cds[m.group(1)].append((f[0], int(f[3]), int(f[4]), f[6]))
    return cds


def compare(ref, bfd, tol=3):
    """ref/bfd: sets of (chrom,start,end,strand). Returns dict of counts."""
    exact = ref & bfd
    r_by = collections.defaultdict(list)
    for c, s, e, st in ref:
        r_by[(c, st)].append((s, e))
    shifted = 0
    for c, s, e, st in bfd - ref:
        if any(abs(s - rs) <= tol and abs(e - re_) <= tol for rs, re_ in r_by.get((c, st), ())):
            shifted += 1
    return dict(ref=len(ref), bfd=len(bfd), exact=len(exact), shifted=shifted)


def selftest():
    cds = {"r1": [("c", 1, 10, "+"), ("c", 21, 30, "+"), ("c", 41, 50, "+")], "r2": [("c", 1, 10, "+"), ("c", 21, 30, "+")]}
    assert ref_introns(cds) == {("c", 11, 20, "+"), ("c", 31, 40, "+")}
    assert ref_introns({"r": [("c", 1, 10, "+"), ("c", 11, 20, "+")]}) == set()  # adjacent CDS, no intron
    r = {("c", 11, 20, "+"), ("c", 31, 40, "+")}
    b = {("c", 11, 20, "+"), ("c", 32, 40, "+"), ("c", 60, 70, "+")}
    assert compare(r, b) == dict(ref=2, bfd=3, exact=1, shifted=1)
    print("selftest ok")


def main(outdir):
    import pyarrow as pa, pyarrow.parquet as pq, pyarrow.compute as pc
    ref = {}
    for tag, (lab, d) in REFS.items():
        cds = parse_gff(f"{SRC}{d}/{d}_genomic.gff.gz")
        ref[tag] = ref_introns(cds)
        print(lab, "reference CDS introns:", len(ref[tag]), file=sys.stderr, flush=True)
    tags = pa.array(list(REFS))
    bfd = {t: {} for t in REFS}   # (chrom,start,end,strand) -> (splice5, splice3)
    pf = pq.ParquetFile(T + "gene_introns.parquet")
    cols = ["intron_id", "chrom", "start", "end", "strand", "splice_5", "splice_3"]
    sample = None
    for b in pf.iter_batches(batch_size=4_000_000, columns=cols):
        t = pa.Table.from_batches([b])
        t = t.filter(pc.is_in(pc.utf8_slice_codeunits(t["intron_id"], 0, 8), value_set=tags))
        for r in t.to_pylist():
            if sample is None:
                sample = r
            bfd[r["intron_id"][:8]][(r["chrom"], r["start"], r["end"], "+" if r["strand"] in (1, "+", "1") else "-")] = (r["splice_5"], r["splice_3"])
    print("example BFD intron row:", sample, file=sys.stderr)
    # chromosome-name overlap check
    with open(outdir + "/ref_vs_bfd.tsv", "w") as fh:
        fh.write("genome\tref_introns\tbfd_introns\texact\tprecision\trecall\tshift_le3bp\tbfd_only_noncanon_pct\tmatched_noncanon_pct\tchrom_overlap\n")
        for tag, (lab, d) in REFS.items():
            r, bset = ref[tag], set(bfd[tag])
            res = compare(r, bset)
            rc = {c for c, *_ in r}; bc = {c for c, *_ in bset}
            def noncanon(keys):
                keys = list(keys)
                if not keys:
                    return float("nan")
                bad = sum(1 for k in keys if (bfd[tag][k][0] or "")[:2].upper() != "GT" or (bfd[tag][k][1] or "")[-2:].upper() != "AG")
                return 100.0 * bad / len(keys)
            fh.write(f"{lab}\t{res['ref']}\t{res['bfd']}\t{res['exact']}\t{res['exact']/max(res['bfd'],1):.4f}\t"
                     f"{res['exact']/max(res['ref'],1):.4f}\t{res['shifted']}\t{noncanon(bset - r):.1f}\t{noncanon(bset & r):.1f}\t{len(rc & bc)}/{len(rc)}\n")
            print(lab, res, "chrom overlap", len(rc & bc), "of", len(rc), file=sys.stderr, flush=True)
            # gene-region recall: reference introns lying in a locus where BFD has any intron on same chrom/strand within span
            idx = collections.defaultdict(list)
            for c, s, e, st in bset:
                idx[(c, st)].append((s, e))
            for k in idx:
                idx[k].sort()
            import bisect
            has_near = 0
            for c, s, e, st in r:
                L = idx.get((c, st), [])
                i = bisect.bisect_left(L, (s - 2000, 0))
                if i < len(L) and L[i][0] <= e + 2000:
                    has_near += 1
            print(f"   ref introns with any BFD intron within 2 kb same strand: {has_near}/{len(r)} "
                  f"(recall inside covered loci ~ {res['exact']/max(has_near,1):.3f})", file=sys.stderr, flush=True)


if __name__ == "__main__":
    selftest() if sys.argv[1] == "--selftest" else main(sys.argv[1])
