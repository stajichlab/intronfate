#!/usr/bin/env python3
"""P1 pilot (b2): intron density on a conserved single-copy PFAM marker set (a proxy for the
BUSCO/own-HMM marker set of Phase 0; BFD has no per-gene BUSCO calls).

Steps (read-only on BFD tables):
 1. Sample genomes stratified by class (<=PER_CLASS each, seed 1) from the P1 TSV, excluding genomes
    below FLOOR x class-median CDS bp.
 2. Stream pfam.parquet restricted to the sample; count single-domain hits per (genome, Pfam).
    Marker Pfams: present in >=PRESENT of sampled genomes, exactly one copy in >=SINGLE of those.
 3. Stream pfam.parquet for the marker Pfams over all genomes; keep (genome, Pfam) with exactly one hit.
 4. Count introns and CDS bp for those transcripts (protein_id == transcript_id in BFD) per genome.
Run:  p1_marker_density.py TSV OUTDIR   |   p1_marker_density.py --selftest
"""
import sys, csv, random, collections, statistics as st

PER_CLASS, FLOOR, PRESENT, SINGLE, EVAL = 40, 0.5, 0.90, 0.90, 1e-5
T = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/"


def select_markers(counts, n_sample, present=PRESENT, single=SINGLE):
    """counts: {(tag, pfam): n_single_domain_hits}. Return sorted Pfams meeting both thresholds."""
    per = collections.defaultdict(list)
    for (tag, pf), n in counts.items():
        per[pf].append(n)
    keep = []
    for pf, ns in per.items():
        if len(ns) / n_sample >= present and sum(n == 1 for n in ns) / len(ns) >= single:
            keep.append(pf)
    return sorted(keep)


def selftest():
    c = {}
    for g in range(10):
        c[(f"G{g}", "PF_single")] = 1                 # kept
        c[(f"G{g}", "PF_multi")] = 2                  # present but multicopy -> dropped
        if g < 5:
            c[(f"G{g}", "PF_rare")] = 1               # present in 50% -> dropped
    c[(("G0"), "PF_single")] = 1
    assert select_markers(c, 10) == ["PF_single"], select_markers(c, 10)
    c[("G1", "PF_single")] = 2                        # 90% single still passes
    assert select_markers(c, 10) == ["PF_single"]
    c[("G2", "PF_single")] = 2                        # 80% single -> dropped
    assert select_markers(c, 10) == []
    print("selftest ok")


def main(tsv, outdir):
    import pyarrow as pa, pyarrow.parquet as pq, pyarrow.compute as pc
    tag = lambda col: pc.utf8_slice_codeunits(col, 0, 8)
    rows = list(csv.DictReader(open(tsv), delimiter="\t"))
    for r in rows:
        r["cds"] = float(r["cds_bp"])
    byc = collections.defaultdict(list)
    for r in rows:
        byc[r["CLASS"]].append(r["cds"])
    ok = [r for r in rows if r["cds"] >= FLOOR * st.median(byc[r["CLASS"]])]
    rnd = random.Random(1)
    sample = []
    bc = collections.defaultdict(list)
    for r in ok:
        bc[r["CLASS"]].append(r["LOCUSTAG"])
    for c in sorted(bc):
        v = sorted(bc[c]); rnd.shuffle(v); sample += v[:PER_CLASS]
    print("sampled genomes:", len(sample), "from", len(bc), "classes", file=sys.stderr, flush=True)
    sset = pa.array(sample)

    pf = pq.ParquetFile(T + "pfam.parquet")
    cols = ["protein_id", "pfam_acc", "domain_num_of", "domain_i_evalue"]
    counts = collections.Counter()
    seen = 0
    for b in pf.iter_batches(batch_size=4_000_000, columns=cols):
        t = pa.Table.from_batches([b])
        m = pc.and_(pc.and_(pc.equal(t["domain_num_of"], 1), pc.less(t["domain_i_evalue"], EVAL)),
                    pc.is_in(tag(t["protein_id"]), value_set=sset))
        t = t.filter(m)
        if t.num_rows:
            g = pa.table({"tag": tag(t["protein_id"]), "pf": t["pfam_acc"]}).group_by(["tag", "pf"]).aggregate([("pf", "count")])
            for a, p, n in zip(g["tag"].to_pylist(), g["pf"].to_pylist(), g["pf_count"].to_pylist()):
                counts[(a, p)] += n
        seen += b.num_rows
        print(f"pass1 {seen:,}", file=sys.stderr, flush=True)
    markers = select_markers(counts, len(sample))
    print("marker Pfams:", len(markers), file=sys.stderr, flush=True)
    with open(outdir + "/marker_pfams.tsv", "w") as fh:
        fh.write("pfam_acc\n" + "\n".join(markers) + "\n")
    mset = pa.array(markers)

    parts = []
    seen = 0
    for b in pf.iter_batches(batch_size=4_000_000, columns=cols):
        t = pa.Table.from_batches([b])
        m = pc.and_(pc.and_(pc.equal(t["domain_num_of"], 1), pc.less(t["domain_i_evalue"], EVAL)),
                    pc.is_in(t["pfam_acc"], value_set=mset))
        t = t.filter(m)
        if t.num_rows:
            parts.append(pa.table({"tag": tag(t["protein_id"]), "pf": t["pfam_acc"], "pid": t["protein_id"]}))
        seen += b.num_rows
        print(f"pass2 {seen:,}", file=sys.stderr, flush=True)
    allp = pa.concat_tables(parts)
    g = allp.group_by(["tag", "pf"]).aggregate([("pid", "count"), ("pid", "min")])
    g = g.filter(pc.equal(g["pid_count"], 1))
    pids = g["pid_min"]
    print("single-copy marker proteins:", len(pids), file=sys.stderr, flush=True)
    pset = pa.array(pids.to_pylist())

    def stream(name, cols, fn):
        tot = collections.Counter()
        f = pq.ParquetFile(T + name)
        for b in f.iter_batches(batch_size=4_000_000, columns=cols):
            t = pa.Table.from_batches([b])
            t = t.filter(pc.is_in(t[cols[0]], value_set=pset))
            if t.num_rows:
                a = pa.table({"tag": tag(t[cols[0]]), "v": fn(t)}).group_by("tag").aggregate([("v", "sum")])
                for k, v in zip(a["tag"].to_pylist(), a["v_sum"].to_pylist()):
                    tot[k] += v
        return tot

    nm = collections.Counter(g["tag"].to_pylist())
    tot = collections.Counter()
    f = pq.ParquetFile(T + "gene_introns.parquet")
    for b in f.iter_batches(batch_size=4_000_000, columns=["transcript_id"]):
        t = pa.Table.from_batches([b])
        t = t.filter(pc.is_in(t["transcript_id"], value_set=pset))
        for k, v in collections.Counter(tag(t["transcript_id"]).to_pylist()).items():
            tot[k] += v
    ni = tot
    cb = stream("gene_CDS.parquet", ["transcript_id", "start", "end"],
                lambda t: pc.add(pc.subtract(t["end"], t["start"]), 1))
    with open(outdir + "/marker_density_per_genome.tsv", "w") as fh:
        fh.write("LOCUSTAG\tn_marker_proteins\tmarker_introns\tmarker_cds_bp\tmarker_introns_per_kb\n")
        for k in sorted(nm):
            bp = cb.get(k, 0)
            fh.write(f"{k}\t{nm[k]}\t{ni.get(k,0)}\t{bp}\t{(ni.get(k,0)/(bp/1000) if bp else float('nan')):.4f}\n")
    print("done", file=sys.stderr)


if __name__ == "__main__":
    selftest() if sys.argv[1] == "--selftest" else main(sys.argv[1], sys.argv[2])
