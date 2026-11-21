#!/usr/bin/env python3
"""J2b step 3: how many matrix cells can test the soft evidence against RefSeq (power check before the likelihood test).

A matrix cell is one (site, taxon) with the intron present in the BFD annotation. It is a *disagreement* cell when RefSeq has no intron within 12 nt on
the same sequence and strand (the `near12` level of matrix_cells.py), and a *zero-read* cell when no junction read lies at its BFD coordinates.
Counted per genome, over all cells and over *informative* sites (present in two or more taxa, the only ones that carry gain/loss/origin information), and
for the hard model's convergent-insertion candidates (`--candidates` = the likelihood_sites.tsv of the hard fit: p_multi_origin > 0.5, and the top 100 by ln LR),
compared with all other sites. Reads: `matrix_cells.reads_at`; RefSeq: the same intron sets as J2.
Run: cell_power.py WORK RUN_DIR HARD_SCORE_DIR OUT_PREFIX | --selftest
"""
import collections, csv, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from matrix_cells import cell_introns, coords, index_ref, within, read_bed, reads_at, S, csv as _csv  # noqa: E402,F401


def cell_table(cells, xy, flag):
    """{(site, taxon): {'agree': bool, 'zero': bool}} from `cells` (label, family, site, protein, off); flag(label, key) -> (agree, zero)"""
    out = {}
    for lab, fam, site, prot, off in cells:
        k = xy.get((prot, off))
        if k is None:
            continue
        a, z = flag(lab, k)
        c = out.setdefault((site, lab), {"agree": False, "zero": True})
        c["agree"] |= a      # one agreeing intron makes the cell agree
        c["zero"] &= z       # the cell has reads when any of its introns does
    return out


def tally(table, sites_n, subset=None):
    """-> {label: [cells, disagree, disagree_with_reads, zero_read]} over the sites in `subset` (all when None)"""
    res = collections.defaultdict(lambda: [0, 0, 0, 0])
    for (site, lab), c in table.items():
        if subset is not None and site not in subset:
            continue
        r = res[lab]
        r[0] += 1
        r[1] += not c["agree"]
        r[2] += (not c["agree"]) and not c["zero"]
        r[3] += c["zero"]
    return res


def selftest():
    cells = [("A", "f", "f:1:0", "p1", 5), ("A", "f", "f:2:0", "p1", 9), ("B", "f", "f:1:0", "p2", 5)]
    xy = {("p1", 5): ("c", 1, 2, "+"), ("p1", 9): ("c", 5, 6, "+"), ("p2", 5): ("c", 1, 2, "+")}
    flag = lambda lab, k: (k[1] == 1, lab == "B")          # agree when first == 1; B has no reads
    t = cell_table(cells, xy, flag)
    assert t[("f:1:0", "A")] == {"agree": True, "zero": False} and t[("f:2:0", "A")]["agree"] is False
    r = tally(t, None)
    assert r["A"] == [2, 1, 1, 0] and r["B"] == [1, 0, 0, 1], dict(r)
    assert tally(t, None, {"f:2:0"})["A"] == [1, 1, 1, 0] and "B" not in tally(t, None, {"f:2:0"})
    print("selftest ok")


def main(work, run_dir, hard_dir, prefix):
    from p1e_classify import parse_models
    outs = {r["out"]: r["ASMID"] for r in csv.DictReader(open(os.path.join(HERE, "..", "P1d_version_audit", "version_audit.tsv")), delimiter="\t")}
    genomes = {}
    for ln in open(os.path.join(HERE, "..", "V6_orthology_nine_species", "genomes.tsv")):
        if not ln.startswith("#") and ln.strip():
            f = ln.rstrip("\n").split("\t")
            genomes[f[0]] = f[1]
    xy, cells = coords(run_dir), cell_introns(run_dir)
    labs = sorted({c[0] for c in cells} & set(genomes))
    ref, junc = {}, {}
    for lab in labs:
        d = genomes[lab]
        bed = os.path.join(work, lab + ".bed")
        if d not in outs or not os.path.exists(bed):
            continue
        a = outs[d]
        m = parse_models(f"{S}{a}/{a}_genomic.gff.gz")
        ref[lab] = index_ref({(x["chrom"], i[0], i[1], x["strand"]) for x in m.values() for i in x["introns"]})
        junc[lab] = read_bed(open(bed))
    cells = [c for c in cells if c[0] in ref]

    def flag(lab, k):
        return within(ref[lab], k, 12), reads_at(junc[lab], k) == 0

    table = cell_table(cells, xy, flag)
    ntaxa = collections.Counter(site for (site, lab) in table)       # taxa with a BFD intron at the site (RefSeq-genome taxa only)
    informative = {s for s, n in ntaxa.items() if n >= 2}
    sites = []
    with open(os.path.join(hard_dir, "likelihood_sites.tsv")) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sites.append((r["site"], float(r["ln_lr_multi"]), float(r["p_multi_origin"])))
    multi = {s for s, _, p in sites if p > 0.5}
    top100 = {s for s, _, _ in sorted(sites, key=lambda x: -x[1])[:100]}
    sets = (("all cells", None), ("informative sites (>=2 taxa)", informative), ("hard candidates p(multi)>0.5", multi), ("hard top 100 by ln LR", top100),
            ("informative non-candidates", informative - multi))
    with open(prefix + "_power.tsv", "w") as fo:
        fo.write("set\tlabel\tcells\tdisagree\tdisagree_with_reads\tzero_read\tfrac_disagree\n")
        for name, sub in sets:
            t = tally(table, None, sub)
            tot = [0, 0, 0, 0]
            for lab in sorted(t):
                v = t[lab]
                tot = [a + b for a, b in zip(tot, v)]
                fo.write(f"{name}\t{lab}\t{v[0]}\t{v[1]}\t{v[2]}\t{v[3]}\t{v[1] / v[0]:.4f}\n")
            if tot[0]:
                fo.write(f"{name}\tALL\t{tot[0]}\t{tot[1]}\t{tot[2]}\t{tot[3]}\t{tot[1] / tot[0]:.4f}\n")
    print(open(prefix + "_power.tsv").read())


if __name__ == "__main__":
    if sys.argv[1:2] == ["--selftest"]:
        selftest()
    else:
        main(*sys.argv[1:5])
