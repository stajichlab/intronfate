#!/usr/bin/env python3
"""Error-rate files for V14 (`intronfate --error-rates`): measured where V4 measured them, a default elsewhere, and a sweep for Alith.

Measured: V4's gffcompare-based intron sensitivity and precision restricted to the V6 tier-A genes (`sn_v2_genes`, `pr_v2_genes` in
bench/V4_annotation_accuracy/results_v6/A_strict/busco_accuracy.tsv), for the V9 genomes that are the same assemblies (V6 and V9 use
the same BFD folders for these six). Default for the other genomes: the median of the measured values. Antarctolithica has no
reference annotation, so its sensitivity is swept (precision kept at the default).

Writes OUT/none.tsv (header only: every genome error-free), OUT/measured.tsv (Alith at the default), OUT/alith_sn<X>.tsv.
Usage: make_error_sets.py --v4 busco_accuracy.tsv --genomes V9/genomes.tsv --out errors [--alith-sn 1.0,0.9,0.8,0.7]
"""
import argparse, csv, os, statistics


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--v4", required=True); ap.add_argument("--genomes", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--alith-sn", default="1.0,0.9,0.8,0.7"); ap.add_argument("--group-prefix", default="Alith")
    a = ap.parse_args(argv)
    labels = [l.split("\t")[0] for l in open(a.genomes) if l.strip() and not l.startswith("#")]
    meas = {}
    for r in csv.DictReader(open(a.v4), delimiter="\t"):
        if r.get("scored") == "yes" and r["genome"] in labels and r.get("sn_v2_genes") and r.get("pr_v2_genes"):
            meas[r["genome"]] = (float(r["sn_v2_genes"]) / 100, float(r["pr_v2_genes"]) / 100)
    sn0 = statistics.median(v[0] for v in meas.values())
    pr0 = statistics.median(v[1] for v in meas.values())
    os.makedirs(a.out, exist_ok=True)

    def write(name, alith_sn):
        with open(os.path.join(a.out, name), "w") as f:
            f.write("taxon\tsensitivity\tprecision\tsource\n")
            if alith_sn is None:
                return
            for g in labels:
                if g in meas:
                    f.write(f"{g}\t{meas[g][0]:.4f}\t{meas[g][1]:.4f}\tV4 measured\n")
                elif g.startswith(a.group_prefix):
                    f.write(f"{g}\t{alith_sn:.4f}\t{pr0:.4f}\tsweep\n")
                else:
                    f.write(f"{g}\t{sn0:.4f}\t{pr0:.4f}\tdefault (median of measured)\n")

    write("none.tsv", None)
    write("measured.tsv", sn0)
    for x in a.alith_sn.split(","):
        write(f"alith_sn{float(x):.2f}.tsv", float(x))
    print(f"measured {len(meas)} genomes: {', '.join(sorted(meas))}; default sensitivity {sn0:.4f}, precision {pr0:.4f}")


if __name__ == "__main__":
    main()
