#!/usr/bin/env python3
"""P1c-c: where do BFD intron calls disagree with the reference? Position, length, gene intron count, gene-level failure.
All genes (not only markers). Per reference RNA, introns are ordered 5'->3' along the transcript.
Missed = reference intron absent from BFD (exact site). Extra = BFD intron absent from reference.
Run: p1c_miss_profile.py OUTDIR | --selftest. Read-only."""
import sys, collections
from p1c_ref_introns import REFS, SRC, T, parse_gff

LEN_BINS = [(0, 50), (50, 70), (70, 100), (100, 200), (200, 500), (500, 10**9)]
NIN_BINS = [(1, 1), (2, 2), (3, 5), (6, 10), (11, 10**9)]

def ordered_introns(segs):
    """segs [(chrom,s,e,strand)] of one RNA -> list of introns in 5'->3' order: (key, idx, n)."""
    segs = sorted(segs, key=lambda x: x[1])
    ints = [(a[0], a[2] + 1, b[1] - 1, a[3]) for a, b in zip(segs, segs[1:]) if b[1] - a[2] > 1]
    if segs and segs[0][3] == "-":
        ints.reverse()
    return [(k, i, len(ints)) for i, k in enumerate(ints)]

def binlab(v, bins):
    for lo, hi in bins:
        if lo <= v <= hi if bins is NIN_BINS else lo <= v < hi:
            return f"{lo}-{hi if hi < 10**9 else 'inf'}"

def relpos(i, n):
    if n == 1: return "only"
    if i == 0: return "first"
    if i == n - 1: return "last"
    return "middle"

def selftest():
    segs = [("c", 1, 10, "+"), ("c", 21, 30, "+"), ("c", 41, 50, "+")]
    o = ordered_introns(segs)
    assert [x[0] for x in o] == [("c", 11, 20, "+"), ("c", 31, 40, "+")] and [x[1] for x in o] == [0, 1]
    neg = [("c", 1, 10, "-"), ("c", 21, 30, "-"), ("c", 41, 50, "-")]
    assert ordered_introns(neg)[0][0] == ("c", 31, 40, "-")   # 5' end of a minus-strand gene is the high coordinate
    assert relpos(0, 1) == "only" and relpos(0, 3) == "first" and relpos(2, 3) == "last" and relpos(1, 3) == "middle"
    assert binlab(55, LEN_BINS) == "50-70" and binlab(3, NIN_BINS) == "3-5" and binlab(1, NIN_BINS) == "1-1"
    print("selftest ok")

def main(outdir):
    import pyarrow as pa, pyarrow.parquet as pq, pyarrow.compute as pc
    tags = pa.array(list(REFS))
    bfd = {t: set() for t in REFS}
    btr = {t: collections.defaultdict(list) for t in REFS}   # transcript -> [(start,end,key,strand)]
    for b in pq.ParquetFile(T + "gene_introns.parquet").iter_batches(batch_size=4_000_000, columns=["intron_id", "transcript_id", "chrom", "start", "end", "strand"]):
        t = pa.Table.from_batches([b])
        t = t.filter(pc.is_in(pc.utf8_slice_codeunits(t["intron_id"], 0, 8), value_set=tags))
        for r in t.to_pylist():
            k = (r["chrom"], r["start"], r["end"], "+" if r["strand"] == 1 else "-")
            tg = r["intron_id"][:8]
            bfd[tg].add(k); btr[tg][r["transcript_id"]].append(k)
    out = open(outdir + "/miss_profile.tsv", "w")
    out.write("genome\tfeature\tcategory\tref_n\tmissed\trecall\n")
    gl = open(outdir + "/miss_gene_level.tsv", "w")
    gl.write("genome\tref_rnas_with_introns\tall_matched\tpartial\tnone_matched\tmissed_introns_in_none_matched_pct\tbfd_only_introns\tbfd_only_in_ref_intronless_region_pct\n")
    for tag, (lab, d) in REFS.items():
        cds = parse_gff(f"{SRC}{d}/{d}_genomic.gff.gz")
        seen = {}
        genes = []
        for rna, segs in cds.items():
            oi = ordered_introns(segs)
            if not oi: continue
            genes.append(oi)
            for k, i, n in oi:
                seen.setdefault(k, (i, n))
        cnt = collections.defaultdict(lambda: [0, 0])
        for k, (i, n) in seen.items():
            miss = k not in bfd[tag]
            ln = k[2] - k[1] + 1
            for feat, cat in (("length", binlab(ln, LEN_BINS)), ("position", relpos(i, n)), ("gene_intron_count", binlab(n, NIN_BINS))):
                cnt[(feat, cat)][0] += 1; cnt[(feat, cat)][1] += miss
        for (feat, cat), (n_, m_) in sorted(cnt.items()):
            out.write(f"{lab}\t{feat}\t{cat}\t{n_}\t{m_}\t{1 - m_/n_:.3f}\n")
        allm = part = none = 0; miss_none = miss_all = 0
        for oi in genes:
            ms = [k in bfd[tag] for k, _, _ in oi]
            nm = len(ms) - sum(ms)
            if nm == 0: allm += 1
            elif sum(ms) == 0: none += 1; miss_none += nm
            else: part += 1
            miss_all += nm
        refset = set(seen)
        extra = bfd[tag] - refset
        # extra intron lying in a locus where the reference has no CDS intron within 2 kb (same chrom/strand)
        import bisect
        idx = collections.defaultdict(list)
        for c, s, e, st in refset: idx[(c, st)].append(s)
        for k in idx: idx[k].sort()
        iso = 0
        for c, s, e, st in extra:
            L = idx.get((c, st), []); j = bisect.bisect_left(L, s - 2000)
            if not (j < len(L) and L[j] <= e + 2000): iso += 1
        gl.write(f"{lab}\t{len(genes)}\t{allm}\t{part}\t{none}\t{100*miss_none/max(miss_all,1):.1f}\t{len(extra)}\t{100*iso/max(len(extra),1):.1f}\n")
        # extra introns: length and position profile vs matched BFD introns
        pos = collections.defaultdict(lambda: [0, 0])
        for tr, ks in btr[tag].items():
            n = len(ks)
            ks2 = sorted(ks, key=lambda k: k[1], reverse=(ks[0][3] == "-"))
            for i, k in enumerate(ks2):
                ex = k not in refset
                pos[("bfd_position", relpos(i, n))][0] += 1; pos[("bfd_position", relpos(i, n))][1] += ex
                lb = binlab(k[2] - k[1] + 1, LEN_BINS)
                pos[("bfd_length", lb)][0] += 1; pos[("bfd_length", lb)][1] += ex
        for (feat, cat), (n_, e_) in sorted(pos.items()):
            out.write(f"{lab}\t{feat}\t{cat}\t{n_}\t{e_}\t{1 - e_/n_:.3f}\n")   # recall column = precision here
    out.close(); gl.close()

if __name__ == "__main__":
    selftest() if sys.argv[1] == "--selftest" else main(sys.argv[1])
