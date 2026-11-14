#!/usr/bin/env python3
"""P1f: adjudicate reference-vs-BFD intron disagreements with the genome's own RNA-seq evidence, and locate the stage
where the model went wrong.

Evidence (read-only): training/transcript.alignments.gff3 (Trinity-GG transcripts aligned by minimap2; junctions from
cDNA_match segments with the same ID) and training/funannotate_train.pasa.gff3 (PASA/TransDecoder gene models).
Intron-level: for introns only in the reference ("ref_only"), only in BFD ("bfd_only") and shared: exact junction support
by transcript alignments, or conflict (an overlapping transcript intron with different coordinates), or no junction.
Gene-level (ref genes whose BFD structure is not identical): verdict from RNA junctions in the gene span:
  ref_right (RNA supports a ref-only intron and no BFD-only intron), bfd_right, mixed, no_junction_support.
For verdict ref_right, where did the correct structure get lost: PASA_has_ref (PASA model carries the ref-only introns:
EVM or later steps overrode it) or PASA_lacks_ref (RNA junctions exist but PASA did not model them: evidence not used).
Run: p1f_adjudicate.py OUT_PREFIX ASMID ... | --selftest
"""
import sys, os, re, csv, bisect, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p1e_classify import parse_models, opener

R = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/"
S = "/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/source/NCBI_ASM/"
MINGAP = 30

def parse_ta(path):
    """transcript alignment introns: set of (chrom,s,e,strand), plus transcript spans [(chrom,strand,s,e)]"""
    segs = collections.defaultdict(list)
    with opener(path) as fh:
        for ln in fh:
            if ln.startswith("#"): continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in ("cDNA_match", "match_part", "expressed_sequence_match"): continue
            m = re.search(r"ID=([^;]+)", f[8])
            if m: segs[(m.group(1), f[0], f[6])].append((int(f[3]), int(f[4])))
    introns, spans = set(), []
    for (tid, c, st), L in segs.items():
        L.sort()
        spans.append((c, st, L[0][0], max(e for _, e in L)))
        for a, b in zip(L, L[1:]):
            if b[0] - a[1] - 1 >= MINGAP: introns.add((c, a[1] + 1, b[0] - 1, st))
    return introns, spans

class Idx:
    """sorted introns per (chrom,strand) with range queries"""
    def __init__(self, introns):
        d = collections.defaultdict(list)
        for c, s, e, st in introns: d[(c, st)].append((s, e))
        self.d = {k: sorted(v) for k, v in d.items()}
        self.starts = {k: [x[0] for x in v] for k, v in self.d.items()}
    def within(self, c, st, s, e):
        v = self.d.get((c, st), []); i = bisect.bisect_left(self.starts.get((c, st), []), s); out = set()
        while i < len(v) and v[i][0] <= e:
            if v[i][1] <= e: out.add((c, v[i][0], v[i][1], st))
            i += 1
        return out
    def overlapping(self, c, st, s, e):
        v = self.d.get((c, st), []); out = []
        for a, b in v[max(0, bisect.bisect_left(self.starts.get((c, st), []), s - 100000)):]:
            if a > e: break
            if b >= s: out.append((a, b))
        return out

def verdict(ref_only_supp, bfd_only_supp, ta_empty):
    if ta_empty: return "no_junction_support"
    if ref_only_supp and not bfd_only_supp: return "ref_right"
    if bfd_only_supp and not ref_only_supp: return "bfd_right"
    if ref_only_supp and bfd_only_supp: return "mixed"
    return "no_junction_support"

def selftest():
    ref = {(1, 2, "+", 0)}
    assert verdict(2, 0, False) == "ref_right" and verdict(0, 1, False) == "bfd_right"
    assert verdict(1, 1, False) == "mixed" and verdict(0, 0, False) == "no_junction_support" and verdict(3, 0, True) == "no_junction_support"
    import tempfile, gzip
    t = tempfile.mkdtemp(); p = os.path.join(t, "ta.gff3")
    open(p, "w").write("##gff-version 3\nc\tgenome\tcDNA_match\t1\t100\t100\t+\t.\tID=T1;Target=T1 1 100\nc\tgenome\tcDNA_match\t201\t300\t100\t+\t.\tID=T1;Target=T1 101 200\n"
                       "c\tgenome\tcDNA_match\t400\t500\t100\t+\t.\tID=T2;Target=T2 1 101\nc\tgenome\tcDNA_match\t510\t600\t100\t+\t.\tID=T2;Target=T2 102 192\n")
    ints, spans = parse_ta(p)
    assert ints == {("c", 101, 200, "+")}, ints          # 9-bp gap (T2) is below MINGAP and ignored
    ix = Idx(ints)
    assert ix.within("c", "+", 1, 300) == {("c", 101, 200, "+")} and ix.within("c", "+", 1, 150) == set()
    assert ix.overlapping("c", "+", 120, 130) == [(101, 200)]
    print("selftest ok")

def main(prefix, asms):
    outs = {r["ASMID"]: r["out"] for r in csv.DictReader(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../P1d_version_audit/version_audit.tsv")), delimiter="\t")}
    genes = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "../P1e_error_classes/genes.tsv")), delimiter="\t"):
        if r["annotation"] == "cur": genes[r["asmid"]].append(r)
    fi = open(prefix + "_intron_support.tsv", "w"); fg = open(prefix + "_gene_verdicts.tsv", "w")
    fi.write("asmid\tset\tlength_class\tn\tta_exact\tta_conflict\tta_none\tpasa_exact\n")
    fg.write("asmid\tclass\tverdict\tpasa_status\tn\n")
    cases = open(prefix + "_cases.tsv", "w"); cases.write("asmid\tclass\tref_tid\tbfd_tid\tverdict\tpasa_status\tchrom\tstrand\tref_start\tref_end\n")
    for a in asms:
        o = outs[a]; tr = f"{R}genome_annotation/{o}/training/"
        tap, pp = tr + "transcript.alignments.gff3", tr + "funannotate_train.pasa.gff3"
        if not os.path.exists(tap):
            print(a, "no transcript alignments; skipped", file=sys.stderr); continue
        ta, tspans = parse_ta(tap)
        pasa = parse_models(pp) if os.path.exists(pp) else {}
        pint = {(m["chrom"], i[0], i[1], m["strand"]) for m in pasa.values() for i in m["introns"]}
        ref = parse_models(f"{S}{a}/{a}_genomic.gff.gz")
        bfdp = f"{R}genome_annotation/{o}/predict_results/{o}"
        bfd = parse_models(bfdp + ".gff3.gz" if os.path.exists(bfdp + ".gff3.gz") else bfdp + ".gff3")
        rint = {(m["chrom"], i[0], i[1], m["strand"]) for m in ref.values() for i in m["introns"]}
        bint = {(m["chrom"], i[0], i[1], m["strand"]) for m in bfd.values() for i in m["introns"]}
        tai = Idx(ta); pai = Idx(pint)
        tsp = collections.defaultdict(list)
        for c, st, s, e in tspans: tsp[(c, st)].append((s, e))
        for k in tsp: tsp[k].sort()
        def covered(c, st, s, e):
            L = tsp.get((c, st), [])
            return any(x <= s and e <= y for x, y in L[:0]) or any(x <= s and y >= e for x, y in L if x <= s)
        # intron level: only introns lying inside a transcript span are assessable
        cnt = collections.defaultdict(lambda: [0, 0, 0, 0, 0])
        for name, S_ in (("shared", rint & bint), ("ref_only", rint - bint), ("bfd_only", bint - rint)):
            for k in S_:
                c, s, e, st = k
                L = tsp.get((c, st), [])
                if not any(x <= s and y >= e for x, y in L): continue
                lc = ">=200" if e - s + 1 >= 200 else "<200"
                v = cnt[(name, lc)]; v[0] += 1
                if k in ta: v[1] += 1
                elif any(not (b == e and a_ == s) for a_, b in tai.overlapping(c, st, s, e)): v[2] += 1
                else: v[3] += 1
                if k in pint: v[4] += 1
        for (name, lc), v in sorted(cnt.items()): fi.write(f"{a}\t{name}\t{lc}\t{v[0]}\t{v[1]}\t{v[2]}\t{v[3]}\t{v[4]}\n")
        # gene level
        refrep = {}
        for t, m in ref.items():
            if m["gene"] not in refrep or m["len"] > ref[refrep[m["gene"]]]["len"]: refrep[m["gene"]] = t
        vc = collections.Counter()
        for r in genes[a]:
            if r["class"] in ("identical", "ends_only", "missing_gene", "opposite_strand_only") or not r["bfd_tid"]: continue
            R_ = ref[r["ref_tid"]]; B_ = bfd[r["bfd_tid"]]
            c, st, s, e = R_["chrom"], R_["strand"], min(R_["s"], B_["s"]), max(R_["e"], B_["e"])
            Rg = {(c, x, y, st) for x, y in R_["introns"]}; Bg = {(c, x, y, st) for x, y in B_["introns"]}
            Tg = tai.within(c, st, s, e)
            ro, bo = Rg - Bg, Bg - Rg
            rs, bs = len(ro & Tg), len(bo & Tg)
            v = verdict(rs, bs, not Tg)
            ps = ""
            if v == "ref_right":
                Pg = pai.within(c, st, s, e)
                ps = "PASA_has_ref" if (ro & Tg) <= Pg else "PASA_lacks_ref"
            elif v == "bfd_right":
                Pg = pai.within(c, st, s, e)
                ps = "PASA_has_bfd" if (bo & Tg) <= Pg else "PASA_lacks_bfd"
            vc[(r["class"], v, ps)] += 1
            cases.write(f"{a}\t{r['class']}\t{r['ref_tid']}\t{r['bfd_tid']}\t{v}\t{ps}\t{c}\t{st}\t{R_['s']}\t{R_['e']}\n")
        for (cl, v, ps), n in sorted(vc.items()): fg.write(f"{a}\t{cl}\t{v}\t{ps}\t{n}\n")
        print(a, "done: TA introns", len(ta), "PASA introns", len(pint), file=sys.stderr, flush=True)
    fi.close(); fg.close(); cases.close()

if __name__ == "__main__":
    selftest() if sys.argv[1] == "--selftest" else main(sys.argv[1], sys.argv[2:])
