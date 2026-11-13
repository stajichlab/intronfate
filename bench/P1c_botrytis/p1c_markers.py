#!/usr/bin/env python3
"""P1c-b: BFD vs reference introns restricted to conserved single-copy marker genes (Pfam set from P1 b2).
Marker genes per reference genome = single-hit proteins for the 453 marker Pfams (protein_id == transcript_id).
Reference side = reference CDS introns lying inside the marker gene spans (same chrom/strand).
Run: p1c_markers.py MARKER_PFAMS OUTDIR | --selftest.  Read-only on BFD and source files."""
import sys, bisect, collections
from p1c_ref_introns import REFS, SRC, T, parse_gff, ref_introns

def spans_introns(spans, ref_sorted):
    """spans: [(chrom,s,e,strand)], ref_sorted: {(chrom,strand): sorted [(s,e)]} -> set of contained introns."""
    out = set()
    for c, s, e, st in spans:
        L = ref_sorted.get((c, st), [])
        i = bisect.bisect_left(L, (s, 0))
        while i < len(L) and L[i][0] <= e:
            if L[i][1] <= e:
                out.add((c, L[i][0], L[i][1], st))
            i += 1
    return out

def selftest():
    r = {("c", 11, 20, "+"), ("c", 31, 40, "+"), ("c", 101, 110, "+"), ("c", 15, 18, "-")}
    rs = collections.defaultdict(list)
    for c, s, e, st in r: rs[(c, st)].append((s, e))
    for k in rs: rs[k].sort()
    assert spans_introns([("c", 1, 50, "+")], rs) == {("c", 11, 20, "+"), ("c", 31, 40, "+")}
    assert spans_introns([("c", 1, 50, "-")], rs) == {("c", 15, 18, "-")}
    assert spans_introns([("c", 1, 30, "+")], rs) == {("c", 11, 20, "+")}  # intron beyond span end excluded
    print("selftest ok")

def main(pf_path, outdir):
    import pyarrow as pa, pyarrow.parquet as pq, pyarrow.compute as pc
    markers = pa.array([l.strip() for l in open(pf_path)][1:])
    tags = pa.array(list(REFS))
    tag = lambda c: pc.utf8_slice_codeunits(c, 0, 8)
    parts = []
    for b in pq.ParquetFile(T + "pfam.parquet").iter_batches(batch_size=4_000_000, columns=["protein_id", "pfam_acc", "domain_num_of", "domain_i_evalue"]):
        t = pa.Table.from_batches([b])
        m = pc.and_(pc.and_(pc.equal(t["domain_num_of"], 1), pc.less(t["domain_i_evalue"], 1e-5)),
                    pc.and_(pc.is_in(t["pfam_acc"], value_set=markers), pc.is_in(tag(t["protein_id"]), value_set=tags)))
        t = t.filter(m)
        if t.num_rows:
            parts.append(pa.table({"tag": tag(t["protein_id"]), "pf": t["pfam_acc"], "pid": t["protein_id"]}))
    allp = pa.concat_tables(parts)
    g = allp.group_by(["tag", "pf"]).aggregate([("pid", "count"), ("pid", "min")])
    g = g.filter(pc.equal(g["pid_count"], 1))
    mk = collections.defaultdict(set)
    for tg, p in zip(g["tag"].to_pylist(), g["pid_min"].to_pylist()):
        mk[tg].add(p)
    allm = pa.array(sorted(set().union(*mk.values())))
    print({REFS[k][0]: len(v) for k, v in mk.items()}, file=sys.stderr, flush=True)
    span = collections.defaultdict(dict)   # tag -> tid -> (chrom,s,e,strand)
    for b in pq.ParquetFile(T + "gene_transcripts.parquet").iter_batches(batch_size=4_000_000, columns=["transcript_id", "chrom", "start", "end", "strand"]):
        t = pa.Table.from_batches([b])
        t = t.filter(pc.is_in(t["transcript_id"], value_set=allm))
        for r in t.to_pylist():
            span[r["transcript_id"][:8]][r["transcript_id"]] = (r["chrom"], r["start"], r["end"], "+" if r["strand"] == 1 else "-")
    bint = collections.defaultdict(lambda: collections.defaultdict(set))  # tag -> tid -> introns
    for b in pq.ParquetFile(T + "gene_introns.parquet").iter_batches(batch_size=4_000_000, columns=["transcript_id", "chrom", "start", "end", "strand"]):
        t = pa.Table.from_batches([b])
        t = t.filter(pc.is_in(t["transcript_id"], value_set=allm))
        for r in t.to_pylist():
            bint[r["transcript_id"][:8]][r["transcript_id"]].add((r["chrom"], r["start"], r["end"], "+" if r["strand"] == 1 else "-"))
    with open(outdir + "/ref_vs_bfd_markers.tsv", "w") as fh:
        fh.write("genome\tmarker_genes\tref_introns_in_spans\tbfd_marker_introns\texact\tprecision\trecall\tgenes_identical_pct\tgenes_ref_has_introns_bfd_none\tgenes_bfd_has_introns_ref_none\n")
        for tag_, (lab, d) in REFS.items():
            ref = ref_introns(parse_gff(f"{SRC}{d}/{d}_genomic.gff.gz"))
            rs = collections.defaultdict(list)
            for c, s, e, st in ref: rs[(c, st)].append((s, e))
            for k in rs: rs[k].sort()
            sp = span[tag_]
            rin = spans_introns(list(sp.values()), rs)
            bset = set().union(*bint[tag_].values()) if bint[tag_] else set()
            ex = rin & bset
            ident = a = b2 = 0
            for tid, spn in sp.items():
                rg = spans_introns([spn], rs); bg = bint[tag_].get(tid, set())
                ident += rg == bg; a += bool(rg) and not bg; b2 += bool(bg) and not rg
            n = max(len(sp), 1)
            fh.write(f"{lab}\t{len(sp)}\t{len(rin)}\t{len(bset)}\t{len(ex)}\t{len(ex)/max(len(bset),1):.4f}\t{len(ex)/max(len(rin),1):.4f}\t{100*ident/n:.1f}\t{a}\t{b2}\n")
            print(lab, len(sp), len(rin), len(bset), len(ex), file=sys.stderr, flush=True)

if __name__ == "__main__":
    selftest() if sys.argv[1] == "--selftest" else main(sys.argv[1], sys.argv[2])
