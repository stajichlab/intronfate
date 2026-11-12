#!/usr/bin/env python3
"""Intron accuracy (funannotate rc.6 vs RefSeq, from a colleague's gffcompare run) measured over ONLY the genes used in the V2
alignments, not the genome-wide rate. Read-only on the shared scoring folder.

Usage: busco_gene_accuracy.py SCORING V2 [--genomes G] [--taxon-map M] [--projection DIR] [--out DIR]
  (the options let the same scoring run over another gene set, e.g. the V6 tiers; the S288C row needs its BFD folder in G)
Inputs
  SCORING   BFD/Fungi_BFD_runs/analysis/rc6_vs_pre_rc6/   (reference_introns/*.new.introns.tsv.gz, reference_gffcompare/*.new.cds.tmap.gz)
  V2        bench/V2_gff_projection/   (taxon_map.tsv: V2 protein id -> genome label; genomes.tsv; results/projection/all/*/{taxa,dollo_branches}.tsv)
Output (./results/)
  busco_accuracy.tsv           one row per genome: genome-wide vs V2-gene intron Sn / Pr
  busco_accuracy_by_class.tsv  the same split by intron phase and by position in the gene
  loss_ceiling.tsv             terminal-branch Dollo losses vs the number annotation misses could explain, using the measured Sn

Definitions (from the per-intron table): pred intron = a row with source=pred; status matched/novel. Ref intron = a distinct
(seqid,start,end,strand) among rows with source=ref; status matched/missed. Pr = matched / (matched + novel) over pred introns of
the V2 transcripts. Sn = matched / (matched + missed) over distinct ref introns of the ref transcripts that the V2 transcripts
map to in the gffcompare .tmap (class codes in ACCEPT, i.e. the reference locus is the same gene).
"""
import collections, csv, glob, gzip, os, sys

import argparse
ap = argparse.ArgumentParser()
ap.add_argument("scoring"); ap.add_argument("v2", help="folder holding genomes.tsv, taxon_map.tsv and results/projection/all/ (the V2 layout)")
ap.add_argument("--genomes", help="genomes.tsv (default: <v2>/genomes.tsv)")
ap.add_argument("--taxon-map", help="protein<TAB>label map of the gene set to score (default: <v2>/taxon_map.tsv)")
ap.add_argument("--projection", help="folder of per-alignment result folders holding taxa.tsv and dollo_branches.tsv "
                                      "(default: <v2>/results/projection/all)")
ap.add_argument("--out", default="results")
ap.add_argument("--min-ref-intron", type=int, default=10, help="ignore reference introns shorter than this (bp): S288C RefSeq has 47 one-base CDS gaps (start==end, frameshifts) that the scoring table lists as introns")
args = ap.parse_args()
SCORING = args.scoring
V2 = args.v2
OUT = args.out
GENOMES = args.genomes or os.path.join(V2, "genomes.tsv")
TAXON_MAP = args.taxon_map or os.path.join(V2, "taxon_map.tsv")
PROJECTION = args.projection or os.path.join(V2, "results", "projection", "all")
ACCEPT = set("=cjkmneo")        # reference transcript shares the locus (exact chain, contained, new junction, ...)
os.makedirs(OUT, exist_ok=True)

genomes = {}
for l in open(GENOMES):
    if l.strip() and not l.startswith("#"):
        label, d, _ = l.rstrip("\n").split("\t")
        genomes[label] = d
tx_of = collections.defaultdict(set)
for l in open(TAXON_MAP):
    pid, lab = l.rstrip("\n").split("\t")
    tx_of[lab].add(pid)

def rows(path):
    with gzip.open(path, "rt") as f:
        r = csv.DictReader(f, delimiter="\t")
        yield from r

def position(row):
    n, k = int(row["n_introns"]), int(row["intron_number"])
    if n == 1: return "only"
    return "first" if k == 1 else ("last" if k == n else "interior")

summary, by_class = [], []
for label, d in genomes.items():
    ipath = f"{SCORING}/reference_introns/{d}.new.introns.tsv.gz"
    tpath = f"{SCORING}/reference_gffcompare/{d}.new.cds.tmap.gz"
    if not (os.path.exists(ipath) and os.path.exists(tpath)):
        summary.append(dict(genome=label, scored="no (no RefSeq comparison for this BFD genome)"))
        continue
    S = tx_of[label]
    # ref transcripts the V2 transcripts map to
    R, classes = set(), collections.Counter()
    found = 0
    with gzip.open(tpath, "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["qry_id"] in S:
                found += 1
                classes[r["class_code"]] += 1
                if r["class_code"] in ACCEPT and r["ref_id"] not in ("", "-"):
                    R.add(r["ref_id"])
    pred_all = collections.Counter(); pred_v2 = collections.Counter()
    ref_all = {}; ref_v2 = {}
    pc = collections.defaultdict(collections.Counter)      # (kind,value) -> {matched/novel}
    rc = collections.defaultdict(lambda: collections.defaultdict(str))
    for r in rows(ipath):
        if r["source"] == "pred":
            pred_all[r["status"]] += 1
            if r["transcript"] in S:
                pred_v2[r["status"]] += 1
                pc[("phase", r["phase"])][r["status"]] += 1
                pc[("position", position(r))][r["status"]] += 1
        else:
            if int(r["end"]) - int(r["start"]) + 1 < args.min_ref_intron: continue
            key = (r["seqid"], r["start"], r["end"], r["strand"])
            ref_all[key] = r["status"]
            if r["transcript"] in R:
                ref_v2[key] = r["status"]
                rc[("phase", r["phase"])][key] = r["status"]
                rc[("position", position(r))][key] = r["status"]
    def sn(dct):
        m = sum(1 for v in dct.values() if v == "matched"); t = len(dct)
        return (m / t * 100 if t else float("nan")), t
    def pr(c):
        t = c["matched"] + c["novel"]
        return (c["matched"] / t * 100 if t else float("nan")), t
    sn_all, n_ref_all = sn(ref_all); sn_v2, n_ref_v2 = sn(ref_v2)
    pr_all, n_pred_all = pr(pred_all); pr_v2, n_pred_v2 = pr(pred_v2)
    summary.append(dict(genome=label, scored="yes", v2_transcripts=len(S), in_tmap=found, with_ref_locus=len(R),
                        class_codes=" ".join(f"{k}:{v}" for k, v in classes.most_common()),
                        sn_genome_wide=round(sn_all, 3), ref_introns_genome=n_ref_all,
                        sn_v2_genes=round(sn_v2, 3), ref_introns_v2=n_ref_v2,
                        pr_genome_wide=round(pr_all, 3), pred_introns_genome=n_pred_all,
                        pr_v2_genes=round(pr_v2, 3), pred_introns_v2=n_pred_v2))
    for kind, value in sorted(set(k for k in list(pc) + list(rc))):
        s, ns = sn(rc[(kind, value)]); p, npred = pr(pc[(kind, value)])
        by_class.append(dict(genome=label, kind=kind, value=value, sn=round(s, 1), n_ref=ns, pr=round(p, 1), n_pred=npred))

def write(path, rs):
    keys = []
    for r in rs:
        for k in r:
            if k not in keys: keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", restval="")
        w.writeheader(); w.writerows(rs)
write(f"{OUT}/busco_accuracy.tsv", summary)
write(f"{OUT}/busco_accuracy_by_class.tsv", by_class)

# ---- loss ceiling with the MEASURED Sn ----
carried, losses = collections.Counter(), collections.Counter()
for d in glob.glob(os.path.join(PROJECTION, "*", "")):
    try:
        for r in csv.DictReader(open(d + "taxa.tsv"), delimiter="\t"):
            carried[r["taxon"]] += int(r["n_sites"])
        for r in csv.DictReader(open(d + "dollo_branches.tsv"), delimiter="\t"):
            if "," not in r["tips"]:
                losses[r["tips"]] += int(r["losses"])
    except FileNotFoundError:
        pass
ceil = []
for r in summary:
    if r.get("scored") != "yes": continue
    g = r["genome"]; s = r["sn_v2_genes"] / 100; sa = r["sn_genome_wide"] / 100
    c = carried[g]; l = losses[g]
    if s != s or sa != sa:        # no reference introns at this gene set (Sn is NaN): there is nothing to bound
        ceil.append(dict(genome=g, sites_carried=c, terminal_losses=l, sn_genome_wide=r["sn_genome_wide"], sn_v2_genes="nan")); continue
    miss_v2 = c * (1 - s) / s if s else float("nan"); miss_all = c * (1 - sa) / sa if sa else float("nan")
    ceil.append(dict(genome=g, sites_carried=c, terminal_losses=l,
                     sn_genome_wide=r["sn_genome_wide"], missed_if_genome_wide=round(miss_all), share_of_losses_genome_wide=round(miss_all / l * 100, 1) if l else "",
                     sn_v2_genes=r["sn_v2_genes"], missed_if_v2_genes=round(miss_v2), share_of_losses_v2_genes=round(miss_v2 / l * 100, 1) if l else ""))
write(f"{OUT}/loss_ceiling.tsv", ceil)
for name in ("busco_accuracy", "loss_ceiling"):
    print(f"== {name}.tsv"); print(open(f"{OUT}/{name}.tsv").read())
