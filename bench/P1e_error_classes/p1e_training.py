#!/usr/bin/env python3
"""P1e-c: what training source did each RefSeq genome's ab initio predictors use (funannotate-predict.log), and does it
explain agreement with the RefSeq annotation (version_audit.tsv)? Read-only. Reads the first 400 lines of each log."""
import csv, re, os, collections, statistics as st
R = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation/"
rows = list(csv.DictReader(open("../P1d_version_audit/version_audit.tsv"), delimiter="\t"))
out = []
for r in rows:
    p = f"{R}{r['out']}/logfiles/funannotate-predict.log"
    d = dict(ASMID=r["ASMID"], out=r["out"], CLASS=r["CLASS"], version=r["version"])
    for k in ("cur_rec", "cur_prec", "cur_rec_ge200", "cur_none_frac", "cur_ref_introns", "cur_mrna_ratio", "has_pasa_gff3"): d[k] = r.get(k, "")
    try:
        txt = []
        with open(p, errors="replace") as fh:
            for i, ln in enumerate(fh):
                txt.append(ln)
                if i > 400: break
        t = "".join(txt)
        for tool in ("augustus", "snap", "glimmerhmm"):
            m = re.search(rf"TRAINING-DECISION[^\n]*{tool} training source = (\w+)", t)
            d[tool + "_src"] = m.group(1) if m else "NA"
        m = re.search(r"augustus --species=(\S+)", t); d["aug_species"] = m.group(1) if m else "NA"
        d["found_training_files"] = int("Found training files" in t)
        d["genemark_skipped"] = int("Will skip GeneMark" in t)
        d["glimmer_failed"] = int("GlimmerHMM training failed" in t)
        m = re.search(r"Summary of gene models: (\{[^}]*\})", t)
        d["summary"] = m.group(1) if m else ""
    except OSError:
        d["augustus_src"] = "no_log"
    out.append(d)
with open("training_by_genome.tsv", "w") as fh:
    keys = list(out[0].keys())
    for d in out:
        for k in d:
            if k not in keys: keys.append(k)
    w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t"); w.writeheader(); w.writerows(out)
f = lambda d, k: float(d[k]) if d.get(k) not in (None, "") else None
ok = [d for d in out if f(d, "cur_rec") is not None and (f(d, "cur_ref_introns") or 0) >= 200]
print("genomes", len(ok))
def table(keyf, title):
    print("\n" + title)
    g = collections.defaultdict(list)
    for d in ok: g[keyf(d)].append(d)
    for k, v in sorted(g.items(), key=lambda kv: -len(kv[1])):
        if len(v) >= 5:
            med = lambda m: st.median([f(d, m) for d in v if f(d, m) is not None])
            print(f"  {str(k)[:44]:44s} n={len(v):4d} rec={med('cur_rec'):.2f} prec={med('cur_prec'):.2f} rec200={med('cur_rec_ge200'):.2f} none={med('cur_none_frac'):.2f}")
table(lambda d: d["augustus_src"], "by Augustus training source")
table(lambda d: (d["augustus_src"], d["snap_src"], d["glimmerhmm_src"]), "by (augustus, snap, glimmerhmm) training source")
table(lambda d: d["aug_species"], "by Augustus species used (n>=5)")
table(lambda d: (d["augustus_src"], d["version"]), "by (augustus source, version)")
table(lambda d: (d["augustus_src"], d["found_training_files"]), "by (augustus source, RNA-seq training files found)")
table(lambda d: d["glimmer_failed"], "GlimmerHMM training failed flag")
table(lambda d: (d["augustus_src"], d["CLASS"] if d["CLASS"] in ("Saccharomycetes","Pichiomycetes","Dipodascomycetes") else "other"), "yeast-like classes vs others, by augustus source")
