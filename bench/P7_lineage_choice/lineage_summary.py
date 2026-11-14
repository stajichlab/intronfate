#!/usr/bin/env python3
"""P7: per-candidate-lineage numbers for choosing the site-level loss analysis sets.

Reads only P1 TSVs (genome-wide density, marker density) and the P4 yeast tip table.
Writes lineage_summary.tsv. Self-test: python3 lineage_summary.py --selftest
"""
import csv, statistics as st, sys, os

P1 = os.path.expanduser("~/projects/intron_projector/bench/P1_intron_density")
P4 = os.path.expanduser("~/projects/intron_projector/bench/P3_reduced_mapping/P4_out/yeast_clade_tips.tsv")

def med(x):
    return st.median(x) if x else float("nan")

def pct(x, q):
    if not x:
        return float("nan")
    s = sorted(x)
    return s[min(len(s) - 1, int(q * len(s)))]

def summarize(rows):
    """rows: list of dicts with genus, markers, marker_introns, density. Returns dict."""
    genera = {r["genus"] for r in rows if r["genus"]}
    ok = [r for r in rows if r["markers"] >= 100]
    mi = [r["marker_introns"] for r in ok]
    return {
        "genomes": len(rows), "genera": len(genera), "genomes_ge100_markers": len(ok),
        "median_markers": med([r["markers"] for r in rows]),
        "median_marker_introns": med(mi), "p90_marker_introns": pct(mi, 0.9),
        "median_density": med([r["density"] for r in rows]),
    }

def selftest():
    rows = [dict(genus="A", markers=450, marker_introns=10, density=0.1),
            dict(genus="A", markers=450, marker_introns=30, density=0.3),
            dict(genus="B", markers=50, marker_introns=1, density=0.0)]
    s = summarize(rows)
    assert s["genomes"] == 3 and s["genera"] == 2 and s["genomes_ge100_markers"] == 2
    assert s["median_marker_introns"] == 20
    print("selftest ok")

def load():
    md = {}
    with open(f"{P1}/marker_density_per_genome.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            md[r["LOCUSTAG"]] = (int(r["n_marker_proteins"]), int(r["marker_introns"]))
    out = []
    with open(f"{P1}/intron_density_per_genome.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["LOCUSTAG"] not in md or not r["introns_per_kb_cds"]:
                continue
            m, mi = md[r["LOCUSTAG"]]
            out.append(dict(tag=r["LOCUSTAG"], phylum=r["PHYLUM"], cls=r["CLASS"], order=r["ORDER"],
                            genus=r["GENUS"], markers=m, marker_introns=mi,
                            density=float(r["introns_per_kb_cds"])))
    return out

def main():
    allr = load()
    groups = {}
    for r in allr:
        groups.setdefault(("class", r["cls"] or "(none) " + r["phylum"]), []).append(r)
        groups.setdefault(("phylum", r["phylum"]), []).append(r)
    yeast = set()
    with open(P4) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            yeast.add(r["GENUS"])
    groups[("clade", "yeast clade genera (132 tips, genus-level match)")] = [r for r in allr if r["genus"] in yeast and r["phylum"] == "Ascomycota"]
    wanted = ["Saccharomycetes", "Pichiomycetes", "Dipodascomycetes", "Trigonopsidomycetes", "Ustilaginomycetes",
              "Malasseziomycetes", "Exobasidiomycetes", "Microbotryomycetes", "Taphrinomycetes",
              "Schizosaccharomycetes", "Kickxellomycetes", "Dimargaritomycetes", "Agaricomycetes",
              "Tremellomycetes", "Pucciniomycetes", "Dothideomycetes", "Sordariomycetes", "Eurotiomycetes",
              "Leotiomycetes", "Pezizomycetes", "Orbiliomycetes", "Neolectomycetes", "Pneumocystidomycetes"]
    keys = [("clade", "yeast clade genera (132 tips, genus-level match)")] + [("class", w) for w in wanted] + \
           [("phylum", p) for p in ("Microsporidia", "Kickxellomycota", "Chytridiomycota", "Mucoromycota", "Zoopagomycota")]
    cols = ["level", "group", "genomes", "genera", "genomes_ge100_markers", "median_markers",
            "median_marker_introns", "p90_marker_introns", "median_density"]
    with open("lineage_summary.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for k in keys:
            rows = groups.get(k, [])
            if not rows:
                continue
            s = summarize(rows)
            o.write("\t".join([k[0], k[1]] + [str(round(s[c], 3)) if isinstance(s[c], float) else str(s[c]) for c in cols[2:]]) + "\n")
    print(open("lineage_summary.tsv").read())

if __name__ == "__main__":
    selftest() if "--selftest" in sys.argv else main()
