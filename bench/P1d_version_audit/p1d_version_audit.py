#!/usr/bin/env python3
"""P1d: agreement of BFD funannotate GFF3 predictions with RefSeq annotations, by funannotate version,
and before/after the rc.6 re-prediction (genome_annotation.pre_rc6 vs current), for all RefSeq (GCF_) genomes in BFD.
Everything is read from GFF3 files (read-only). Per genome and annotation: intron-site precision/recall, recall of
introns >=200 bp, fraction of reference RNAs with introns that have no matched intron (gene-level failure),
mRNA-count ratio. Run: p1d_version_audit.py OUT.tsv [NPROC] | --selftest"""
import sys, csv, re, os, gzip, collections, multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p1c_ref_introns import parse_gff, ref_introns

R = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/"
S = "/bigdata/stajichlab/shared/projects/1KFG/2026/NCBI_fungi/source/NCBI_ASM/"

def tag(species, strain):
    s = strain.split(";")[0].replace(":", " ").strip().strip("*").replace("*", "-")
    t = (species + " " + s).replace('"', "").replace("'", "")
    return re.sub(r"[\s/#\[\]\?\{\}]+", "_", t.strip())

def metrics(ref_cds, qry_cds):
    """ref_cds/qry_cds: {rna: [(chrom,s,e,strand)]}. Returns dict of agreement metrics."""
    ri, qi = ref_introns(ref_cds), ref_introns(qry_cds)
    ex = ri & qi
    long_ref = {k for k in ri if k[2] - k[1] + 1 >= 200}
    ref_rnas = [ref_introns({r: v}) for r, v in ref_cds.items()]
    ref_rnas = [x for x in ref_rnas if x]
    none = sum(1 for x in ref_rnas if not (x & qi))
    return dict(ref_introns=len(ri), q_introns=len(qi), exact=len(ex),
                prec=len(ex) / max(len(qi), 1), rec=len(ex) / max(len(ri), 1),
                rec_ge200=len(long_ref & qi) / max(len(long_ref), 1), n_long=len(long_ref),
                none_frac=none / max(len(ref_rnas), 1), mrna_ratio=len(qry_cds) / max(len(ref_cds), 1))

def selftest():
    ref = {"r1": [("c", 1, 10, "+"), ("c", 21, 30, "+")], "r2": [("c", 101, 110, "+"), ("c", 401, 410, "+")]}
    q = {"q1": [("c", 1, 10, "+"), ("c", 21, 30, "+")], "q2": [("c", 101, 110, "+"), ("c", 411, 420, "+")]}
    m = metrics(ref, q)
    assert (m["ref_introns"], m["q_introns"], m["exact"]) == (2, 2, 1)
    assert m["exact"] == 1 and abs(m["rec"] - 0.5) < 1e-9 and m["n_long"] == 1 and m["rec_ge200"] == 0.0
    assert abs(m["none_frac"] - 0.5) < 1e-9 and m["mrna_ratio"] == 1.0
    assert tag("Abortiporus biennis", "LGAM 436") == "Abortiporus_biennis_LGAM_436"
    print("selftest ok")

def gff_path(base, out):
    for ext in (".gff3.gz", ".gff3"):
        p = f"{base}/{out}/predict_results/{out}{ext}"
        if os.path.exists(p):
            return p

def parse_any(p):
    if p.endswith(".gz"):
        return parse_gff(p)
    cds = collections.defaultdict(list)
    for ln in open(p):
        if ln.startswith("#"): continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        m = re.search(r"Parent=([^;]+)", f[8])
        if m: cds[m.group(1)].append((f[0], int(f[3]), int(f[4]), f[6]))
    return cds

def pre_version(out):
    p = f"{R}genome_annotation.pre_rc6/{out}/logfiles/funannotate-predict.log"
    try:
        for i, ln in enumerate(open(p, errors="replace")):
            m = re.search(r"Running funannotate (\S+)", ln)
            if m: return m.group(1)
            if i > 30: break
    except OSError:
        pass
    return ""

def work(a):
    asm, out, cls, ver = a
    refp = f"{S}{asm}/{asm}_genomic.gff.gz"
    row = dict(ASMID=asm, out=out, CLASS=cls, version=ver)
    try:
        ref = parse_gff(refp)
        cur = gff_path(R + "genome_annotation", out)
        if cur:
            row.update({"cur_" + k: v for k, v in metrics(ref, parse_any(cur)).items()})
        pre = gff_path(R + "genome_annotation.pre_rc6", out)
        if pre:
            row.update({"pre_" + k: v for k, v in metrics(ref, parse_any(pre)).items()})
            row["pre_version"] = pre_version(out)
        row["has_pasa_gff3"] = int(os.path.exists(f"{R}genome_annotation/{out}/training/funannotate_train.pasa.gff3"))
    except Exception as e:  # keep going; record failure
        row["error"] = repr(e)[:120]
    return row

def main(outp, nproc):
    ver = {r["out"]: r["funannotate_version"] for r in csv.DictReader(open(R + "do_annotation_wave1/funannotate_version_by_genome_20261021.tsv"), delimiter="\t")}
    jobs = []
    for r in csv.DictReader(open(R + "samples.csv")):
        if r["ASMID"].startswith("GCF_"):
            o = tag(r["SPECIES"], r["STRAIN"])
            if o in ver:
                jobs.append((r["ASMID"], o, r["CLASS"], ver[o]))
    print("genomes:", len(jobs), file=sys.stderr, flush=True)
    with mp.Pool(nproc) as p:
        rows = p.map(work, jobs, chunksize=4)
    keys = sorted({k for r in rows for k in r}, key=lambda k: (not k in ("ASMID", "out", "CLASS", "version"), k))
    with open(outp, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print("done; errors:", sum("error" in r for r in rows), file=sys.stderr)

if __name__ == "__main__":
    selftest() if sys.argv[1] == "--selftest" else main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 4)
