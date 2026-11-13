#!/usr/bin/env python3
"""P1e: classify disagreements between BFD funannotate gene models and RefSeq models, gene by gene.

For each reference gene (representative = longest-CDS transcript) find the best BFD transcript (same strand, most
codon-frame-consistent shared CDS bases). Reports, per gene, amino acids shared / missing from BFD / extra in BFD
(frame-aware, coordinate based; no sequence needed) and a primary class:
  identical, ends_only, shifted_boundary, exon_difference, retained_intron, end_difference, different_structure,
  split (ref gene covered by >=2 BFD genes), merged (BFD gene covers >=2 ref genes), opposite_strand_only, missing_gene.
Intron-level events for reference introns absent from BFD: retained, shifted(<=50bp), intron_differs, outside_span;
BFD introns absent from the reference: splices_out_ref_coding (BFD removes reference coding sequence), extra_noncoding.
Run: p1e_classify.py OUTDIR ASMID,out,[cur|pre] ...  (see sbatch) | --selftest. Read-only.
"""
import sys, re, os, gzip, bisect, collections

TOL = 50

def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)

def parse_models(path):
    """-> {tid: dict(chrom,strand,gene,cds=[(s,e,off)] 5'->3' order with cumulative frame offset, s,e span, len)}"""
    parent, cds = {}, collections.defaultdict(list)
    with opener(path) as fh:
        for ln in fh:
            if ln.startswith("#"): continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9: continue
            a = f[8]
            i = re.search(r"(?:^|;)ID=([^;]+)", a); p = re.search(r"(?:^|;)Parent=([^;]+)", a)
            if f[2] == "CDS" and p:
                cds[p.group(1)].append((f[0], int(f[3]), int(f[4]), f[6]))
            elif f[2] in ("mRNA", "transcript") and i:
                parent[i.group(1)] = p.group(1) if p else i.group(1)
    out = {}
    for t, segs in cds.items():
        chrom, strand = segs[0][0], segs[0][3]
        out[t] = mk(chrom, strand, parent.get(t, t), [(s, e) for _, s, e, _ in segs])
    return out

def mk(chrom, strand, gene, segs):
    segs = sorted(segs)
    order = segs if strand == "+" else segs[::-1]
    off, withoff = 0, []
    for s, e in order:
        withoff.append((s, e, off)); off += e - s + 1
    withoff.sort()
    introns = [(a[1] + 1, b[0] - 1) for a, b in zip(withoff, withoff[1:]) if b[0] - a[1] > 1]
    return dict(chrom=chrom, strand=strand, gene=gene, cds=withoff, s=withoff[0][0], e=withoff[-1][1], len=off, introns=introns)

def cp(seg, p, strand):
    s, e, off = seg
    return (off + (p - s)) % 3 if strand == "+" else (off + (e - p)) % 3

def shared_inframe(R, B):
    """bases in both CDS at the same codon position (same strand assumed), and total overlapping bases."""
    i = j = inf = tot = 0
    A, C = R["cds"], B["cds"]
    while i < len(A) and j < len(C):
        a, b = max(A[i][0], C[j][0]), min(A[i][1], C[j][1])
        if a <= b:
            tot += b - a + 1
            if cp(A[i], a, R["strand"]) == cp(C[j], a, B["strand"]): inf += b - a + 1
        if A[i][1] < C[j][1]: i += 1
        else: j += 1
    return inf, tot

def cov(B, s, e):
    """fraction of [s,e] covered by B CDS"""
    n = 0
    for a, b, _ in B["cds"]:
        x, y = max(a, s), min(b, e)
        if x <= y: n += y - x + 1
    return n / (e - s + 1)

def classify(R, B):
    """-> (primary class, events list [(kind, intron_len)], shared_aa, missing_aa, extra_aa)"""
    inf, tot = shared_inframe(R, B)
    sa, ma, ea = inf / 3, (R["len"] - inf) / 3, (B["len"] - inf) / 3
    ri, bi = set(R["introns"]), set(B["introns"])
    ev = []
    unmatched_b = list(bi - ri)
    for (s, e) in sorted(ri - bi):
        L = e - s + 1
        if cov(B, s, e) >= 0.9: ev.append(("retained", L)); continue
        if any(abs(bs - s) <= TOL and abs(be - e) <= TOL for bs, be in unmatched_b): ev.append(("shifted", L)); continue
        if e < B["s"] or s > B["e"]: ev.append(("outside_span", L)); continue
        ev.append(("intron_differs", L))
    for (s, e) in sorted(bi - ri):
        if any(abs(rs - s) <= TOL and abs(re_ - e) <= TOL for rs, re_ in ri - bi): continue  # counted as shifted
        ev.append(("splices_out_ref_coding" if cov(R, s, e) >= 0.1 else "extra_noncoding", e - s + 1))
    same_cds = [(a, b) for a, b, _ in R["cds"]] == [(a, b) for a, b, _ in B["cds"]]
    if same_cds: return "identical", ev, sa, ma, ea
    if inf / max(R["len"], B["len"]) < 0.5: return "different_structure", ev, sa, ma, ea
    kinds = {k for k, _ in ev}
    if ri == bi: return "ends_only", ev, sa, ma, ea
    if "retained" in kinds: c = "retained_intron"
    elif kinds & {"splices_out_ref_coding", "intron_differs"}: c = "exon_difference"
    elif "shifted" in kinds: c = "shifted_boundary"
    elif kinds & {"outside_span"}: c = "end_difference"
    else: c = "end_difference" if abs(R["s"] - B["s"]) + abs(R["e"] - B["e"]) > 9 else "other"
    return c, ev, sa, ma, ea

def index(models):
    idx = collections.defaultdict(list)
    for t, m in models.items(): idx[(m["chrom"], m["strand"])].append((m["s"], m["e"], t))
    pm = {}
    for k, L in idx.items():
        L.sort(); m = 0; P = []
        for s, e, t in L: m = max(m, e); P.append(m)
        pm[k] = P
    return idx, pm

def overlapping(idx, pm, chrom, strand, s, e):
    L = idx.get((chrom, strand), [])
    if not L: return []
    starts = [x[0] for x in L]
    i = bisect.bisect_right(starts, e) - 1
    out = []
    P = pm[(chrom, strand)]
    while i >= 0 and P[i] >= s:
        if L[i][1] >= s: out.append(L[i][2])
        i -= 1
    return out

def run_pair(ref_models, qry_models):
    # representative reference transcript per gene
    rep = {}
    for t, m in ref_models.items():
        g = m["gene"]
        if g not in rep or m["len"] > ref_models[rep[g]]["len"]: rep[g] = t
    ridx = {t: ref_models[t] for t in rep.values()}
    qidx, qpm = index(qry_models)
    opp_idx, opp_pm = index({t: dict(m, strand="-" if m["strand"] == "+" else "+") for t, m in qry_models.items()})
    q_to_refs = collections.defaultdict(list)
    rows = []
    for t, R in ridx.items():
        cands = overlapping(qidx, qpm, R["chrom"], R["strand"], R["s"], R["e"])
        best, bi, shares = None, -1, []
        for q in cands:
            inf, tot = shared_inframe(R, qry_models[q])
            if tot: shares.append((inf, tot, q))
        if not shares:
            opp = overlapping(opp_idx, opp_pm, R["chrom"], R["strand"], R["s"], R["e"])
            rows.append((t, R, None, "opposite_strand_only" if opp else "missing_gene", [], 0, R["len"] / 3, 0))
            continue
        shares.sort(reverse=True)
        q = shares[0][2]
        for inf, tot, qq in shares:
            if inf >= 0.2 * R["len"]: q_to_refs[qq].append(t)
        c, ev, sa, ma, ea = classify(R, qry_models[q])
        second = shares[1][0] if len(shares) > 1 else 0
        rows.append((t, R, q, c, ev, sa, ma, ea, second))
    out = []
    for r in rows:
        if len(r) == 8: t, R, q, c, ev, sa, ma, ea = r; second = 0
        else: t, R, q, c, ev, sa, ma, ea, second = r
        if q and c != "missing_gene":
            if second >= 0.2 * R["len"] and c != "identical": c = "split"
            if len(q_to_refs[q]) >= 2 and c != "identical": c = "merged"
        out.append(dict(ref=t, q=q or "", cls=c, chrom=R["chrom"], strand=R["strand"], rs=R["s"], re=R["e"],
                        ref_aa=R["len"] // 3, shared_aa=int(sa), missing_aa=int(ma), extra_aa=int(ea), n_ref_introns=len(R["introns"]), events=ev))
    return out

def selftest():
    P = lambda segs, st="+": mk("c", st, "g", segs)
    R = P([(1, 90), (201, 290), (399, 488)])                  # 3 exons of 90 bases = 90 aa; introns 110 and 108 bp
    assert classify(R, P([(1, 90), (201, 290), (399, 488)]))[0] == "identical"
    c, ev, sa, ma, ea = classify(R, P([(1, 96), (207, 290), (399, 488)]))      # intron 1 moved 6 bp, in frame
    assert c == "shifted_boundary" and ("shifted", 110) in ev and ma == 2 and ea == 2, (c, ev, ma, ea)  # 6 bp moved: 2 aa replaced
    c, ev, *_ = classify(R, P([(1, 90), (201, 488)]))                          # intron 2 (108 bp) retained
    assert c == "retained_intron" and ("retained", 108) in ev, (c, ev)
    c, ev, sa, ma, ea = classify(R, P([(1, 90), (399, 488)]))                  # middle exon skipped
    assert c == "exon_difference" and ma == 30 and ea == 0, (c, ev, ma, ea)
    c, ev, sa, ma, ea = classify(R, P([(1, 90), (201, 290)]))                  # last exon missing
    assert c == "end_difference" and ma == 30 and ea == 0, (c, ev, ma, ea)
    c, *_ = classify(R, P([(2, 90), (202, 290), (400, 488)]))                  # frame shifted
    assert c == "different_structure", c
    Rm = P([(1, 90), (201, 290), (399, 488)], "-")
    assert classify(Rm, P([(1, 90), (201, 290), (399, 488)], "-"))[0] == "identical"
    c, ev, sa, ma, ea = classify(Rm, P([(201, 290), (399, 488)], "-"))         # minus strand: lowest exon = C-terminal, missing
    assert ma == 30 and ea == 0 and c == "end_difference", (c, ev, ma, ea)
    ref = {"r1": mk("c", "+", "g1", [(1, 90), (201, 290)]), "r2": mk("c", "+", "g2", [(1000, 1089), (1200, 1289)])}
    q = {"q1": mk("c", "+", "h", [(1, 90)]), "q2": mk("c", "+", "h2", [(201, 290)])}
    res = {r["ref"]: r["cls"] for r in run_pair(ref, q)}
    assert res["r1"] == "split" and res["r2"] == "missing_gene", res
    ref2 = {"r1": mk("c", "+", "g1", [(1, 90)]), "r2": mk("c", "+", "g2", [(300, 389)])}
    q2 = {"q": mk("c", "+", "h", [(1, 90), (300, 389)])}
    res = {r["ref"]: r["cls"] for r in run_pair(ref2, q2)}
    assert res["r1"] == "merged" and res["r2"] == "merged", res
    ref3 = {"r": mk("c", "+", "g", [(1, 90)])}
    opp = {"q": mk("c", "-", "h", [(1, 90)])}
    assert run_pair(ref3, opp)[0]["cls"] == "opposite_strand_only"
    print("selftest ok")

def main(outdir, specs):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    R = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/"
    S = "/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/source/NCBI_ASM/"
    def gpath(base, out):
        for ext in (".gff3.gz", ".gff3"):
            p = f"{base}/{out}/predict_results/{out}{ext}"
            if os.path.exists(p): return p
    genes = open(outdir + "/genes.tsv", "w"); events = open(outdir + "/events.tsv", "w")
    genes.write("asmid\tannotation\tref_tid\tbfd_tid\tclass\tchrom\tstrand\tref_start\tref_end\tref_aa\tshared_aa\tmissing_aa\textra_aa\tn_ref_introns\n")
    events.write("asmid\tannotation\tclass\tevent\tintron_len\n")
    import csv
    outs = {r["ASMID"]: r["out"] for r in csv.DictReader(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../P1d_version_audit/version_audit.tsv")), delimiter="\t")}
    for asm in specs:
        out = outs[asm]
        ref = parse_models(f"{S}{asm}/{asm}_genomic.gff.gz")
        for ann, base in (("cur", R + "genome_annotation"), ("pre", R + "genome_annotation.pre_rc6")):
            p = gpath(base, out)
            if not p: continue
            res = run_pair(ref, parse_models(p))
            for r in res:
                genes.write("\t".join(map(str, [asm, ann, r["ref"], r["q"], r["cls"], r["chrom"], r["strand"], r["rs"], r["re"], r["ref_aa"], r["shared_aa"], r["missing_aa"], r["extra_aa"], r["n_ref_introns"]])) + "\n")
                for k, L in r["events"]: events.write(f"{asm}\t{ann}\t{r['cls']}\t{k}\t{L}\n")
            print(asm, ann, len(res), file=sys.stderr, flush=True)
    genes.close(); events.close()

if __name__ == "__main__":
    if sys.argv[1] == "--selftest": selftest()
    else: main(sys.argv[1], sys.argv[2:])
