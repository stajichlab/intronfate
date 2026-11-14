#!/usr/bin/env python3
"""P4: inside the yeast (Saccharomycotina) reduced clade. Read-only. Run from bench/P3_reduced_mapping (imports reduced_rooted).
Run: yeast_clade.py TREE DENSITY_TSV MARKER_TSV OUTDIR | --selftest"""
import sys, csv, collections, random, statistics as st, math, os
sys.path.insert(0, ".")
from reduced_mapping import parse_tree, masks
from reduced_rooted import root_on, sankoff
YEAST = {"Pichiomycetes", "Saccharomycetes", "Dipodascomycetes", "Lipomycetes", "Trigonopsidomycetes", "Alaromycetes", "Phaffomycetes", "Saccharomycotina"}

def summarize(tips, fullc, nk, state, cls):
    tot, rs, g, v, gc = sankoff(tips, fullc, nk, state, "stay")
    yg = [c for c in gc if sum(1 for i in range(len(tips)) if (c >> i) & 1 and cls[tips[i]] in YEAST) * 2 > bin(c).count("1")]
    big = max(gc, key=lambda c: bin(c).count("1")) if gc else 0
    return tot, g, v, len(yg), bin(big).count("1")

def selftest():
    t, ch, r = parse_tree("(((o1,o2),a),(b,(c,d)));")
    O = 3; full, nk = root_on(t, ch, O)
    cls = {x: "Pichiomycetes" for x in t}
    s = {"o1": 0, "o2": 0, "a": 0, "b": 1, "c": 1, "d": 1}
    assert summarize(t, full, nk, s, cls) == (1, 1, 0, 1, 3)
    print("selftest ok")

def main(tree, dens_tsv, mark_tsv, outdir):
    os.makedirs(outdir, exist_ok=True)
    tips, children, root = parse_tree(open(tree).read())
    dens = {r["LOCUSTAG"]: r for r in csv.DictReader(open(dens_tsv), delimiter="\t")}
    mark = {r["LOCUSTAG"]: int(r["n_marker_proteins"]) for r in csv.DictReader(open(mark_tsv), delimiter="\t")}
    tax = {r["LOCUSTAG"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/samples.csv"))}
    bycls = collections.defaultdict(list)
    for r in dens.values(): bycls[r["CLASS"]].append(float(r["cds_bp"]))
    med = {c: st.median(v) for c, v in bycls.items()}
    ok = {t: t in dens and float(dens[t]["cds_bp"]) >= 0.5 * med[dens[t]["CLASS"]] and mark.get(t, 0) >= 100 for t in tips}
    cls = {t: (tax[t]["CLASS"] if t in tax else "") for t in tips}
    O = sum(1 << i for i, t in enumerate(tips) if t in tax and tax[t]["PHYLUM"] == "Chytridiomycota")
    fullc, nk = root_on(tips, children, O)
    def mkstate(thr, f=lambda t, d: d):
        return {t: (None if (t not in dens or not ok[t]) else (1 if f(t, float(dens[t]["introns_per_kb_cds"])) < thr else 0)) for t in tips}
    # 1. interrupters of the main yeast clade at 0.5
    state = mkstate(0.5)
    tot, rs, g, v, gc = sankoff(tips, fullc, nk, state, "stay")
    big = max(gc, key=lambda c: bin(c).count("1"))
    idx = [i for i in range(len(tips)) if (big >> i) & 1]
    with open(f"{outdir}/yeast_clade_tips.tsv", "w") as o:
        o.write("LOCUSTAG\tCLASS\tORDER\tGENUS\tSPECIES\tdensity\tcds_bp\tmarkers\tN50\tcontigs\tcore\tstate\n")
        for i in idx:
            t = tips[i]; d = dens.get(t, {})
            o.write("\t".join(str(x) for x in [t, d.get("CLASS", ""), d.get("ORDER", ""), d.get("GENUS", ""), d.get("SPECIES", ""), d.get("introns_per_kb_cds", ""), d.get("cds_bp", ""), mark.get(t, ""), d.get("N50_bp", ""), d.get("contig_count", ""), ok[t], state[t]]) + "\n")
    inter = [tips[i] for i in idx if state[tips[i]] == 0]
    print(f"main clade: {len(idx)} tips; reduced {sum(state[tips[i]]==1 for i in idx)}; unknown {sum(state[tips[i]] is None for i in idx)}; interrupters (core, not reduced) {len(inter)}")
    for t in sorted(inter, key=lambda t: float(dens[t]["introns_per_kb_cds"])):
        d = dens[t]; print(f"   {d['CLASS'][:16]:16s} {d['ORDER'][:18]:18s} {d['GENUS'][:22]:22s} dens {float(d['introns_per_kb_cds']):.2f} CDS {float(d['cds_bp'])/1e6:5.1f}Mb markers {mark.get(t,0):3d} N50 {float(d['N50_bp'])/1e6:5.2f}Mb")
    dv = [float(dens[tips[i]]["introns_per_kb_cds"]) for i in idx if ok[tips[i]]]
    qs = sorted(dv); q = lambda p: qs[int(p * (len(qs) - 1))]
    print(f"density in main clade (core tips, n={len(dv)}): p5 {q(.05):.2f} p25 {q(.25):.2f} p50 {q(.5):.2f} p75 {q(.75):.2f} p95 {q(.95):.2f} max {qs[-1]:.2f}")
    byc = collections.defaultdict(list)
    for i in idx:
        if ok[tips[i]]: byc[cls[tips[i]]].append(float(dens[tips[i]]["introns_per_kb_cds"]))
    for c, v in sorted(byc.items(), key=lambda kv: -len(kv[1])): print(f"   {c:20s} n={len(v):3d} median {st.median(v):.3f} max {max(v):.2f}")
    # 2. threshold sweep
    print("\nthreshold sweep (yeast-majority gain clades; size of largest gain clade):")
    for thr in (0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5):
        s = mkstate(thr); tot, g_, v_, ny, bg = summarize(tips, fullc, nk, s, cls)
        print(f"  <{thr}: min changes {tot}; reductions {g_}; re-expansions {v_}; yeast gain clades {ny}; largest clade {bg}")
    # 3. error sensitivity: yeast tips scaled by factor (miss -> true density higher) and lognormal noise on all tips
    print("\nscaling density of yeast-class tips (x factor) at 0.5/kb:")
    for fac in (0.5, 0.75, 1.0, 1.5, 2.0, 3.0):
        s = mkstate(0.5, lambda t, d: d * fac if cls[t] in YEAST else d); tot, g_, v_, ny, bg = summarize(tips, fullc, nk, s, cls)
        print(f"  x{fac}: min changes {tot}; reductions {g_}; re-expansions {v_}; yeast gain clades {ny}; largest clade {bg}")
    print("\nlognormal noise on every tip's density (200 replicates, seed 1), threshold 0.5:")
    random.seed(1)
    for sig in (0.2, 0.4, 0.7):
        res = []
        for _ in range(200):
            noise = {t: math.exp(random.gauss(0, sig)) for t in tips}
            s = mkstate(0.5, lambda t, d: d * noise[t]); res.append(summarize(tips, fullc, nk, s, cls))
        def rng(k): x = sorted(r[k] for r in res); return f"{st.median(x)} [{x[5]}-{x[194]}]"
        print(f"  sigma {sig}: min changes {rng(0)}; reductions {rng(1)}; re-expansions {rng(2)}; yeast gain clades {rng(3)}; largest clade {rng(4)}")
if __name__ == "__main__":
    selftest() if sys.argv[1] == "--selftest" else main(*sys.argv[1:5])
