#!/usr/bin/env python3
"""Ortholog families for a set with near-identical strains, from OrthoFinder root orthogroups (Orthogroups.tsv).

For each orthogroup, a genome with exactly one gene contributes that gene; a genome with two or more genes (ambiguous paralogs) or none is MISSING in
the family (no row, so intronfate reads `?`, never a loss). Proteins with an internal stop are dropped too (an aligner would not keep them).
Strains of one lineage (genome labels starting with --group-prefix) count as ONE lineage when the family is classified: the lineage is present when any
strain contributes, and the family must contain it. With S lineages in all:
  tier A : at least A_MIN lineages present        tier B : B_MIN <= lineages < A_MIN
Writes <out>/<tier>/families/<OG>.faa (proteins without the terminal '*'), <out>/<tier>/taxon_map.tsv (protein<TAB>label), <out>/summary.tsv.
Usage: make_families_pruned.py --of results/plain --data data --out families --group-prefix Alith --tiers A=22,B=18
"""
import argparse, csv, gzip, os, sys

def read_fasta(path):
    seqs, name, buf = {}, None, []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name: seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--of", required=True); ap.add_argument("--data", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--group-prefix", default="", help="genome labels with this prefix are strains of one lineage")
    ap.add_argument("--tiers", default="A=22,B=18", help="A=<min lineages>,B=<min lineages>")
    a = ap.parse_args(argv)
    t = dict(x.split("=") for x in a.tiers.split(","))
    a_min, b_min = int(t["A"]), int(t["B"])
    if not b_min <= a_min: sys.exit("B must not exceed A")

    with open(os.path.join(a.of, "Orthogroups", "Orthogroups.tsv"), newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        names = r.fieldnames[1:]
        rows = [(row["Orthogroup"], {n: [g for g in row[n].replace(", ", ",").split(",") if g] if row.get(n) else [] for n in names}) for row in r]
    group = [n for n in names if a.group_prefix and n.startswith(a.group_prefix)]
    others = [n for n in names if n not in group]
    n_lineages = len(others) + (1 if group else 0)
    prot = {n: read_fasta(os.path.join(a.data, n, f"{n}.proteins.fa.gz")) for n in names}

    tiers = {"A": {}, "B": {}}; dropped_stop = 0
    for og, c in rows:
        pick = {}
        for n in names:
            if len(c[n]) != 1: continue
            s = prot[n].get(c[n][0])
            if s is None: sys.exit(f"{c[n][0]} ({n}) not in the proteome")
            if "*" in s.rstrip("*"): dropped_stop += 1; continue
            pick[n] = (c[n][0], s.rstrip("*"))
        if group and not any(n in pick for n in group): continue
        lineages = sum(1 for n in others if n in pick) + (1 if any(n in pick for n in group) else 0)
        if lineages >= a_min: tiers["A"][og] = pick
        elif lineages >= b_min: tiers["B"][og] = pick

    summary = []
    for tier, table in tiers.items():
        tdir = os.path.join(a.out, tier); os.makedirs(os.path.join(tdir, "families"), exist_ok=True)
        with open(os.path.join(tdir, "taxon_map.tsv"), "w") as tm:
            for og, pick in sorted(table.items()):
                with open(os.path.join(tdir, "families", f"{og}.faa"), "w") as out:
                    for n, (g, s) in sorted(pick.items()):
                        out.write(f">{g}\n{s}\n"); tm.write(f"{g}\t{n}\n")
        summary.append((tier, len(table), sum(len(p) for p in table.values())))
    with open(os.path.join(a.out, "summary.tsv"), "w") as f:
        f.write("tier\tfamilies\tsequences\n")
        for s in summary: f.write("\t".join(map(str, s)) + "\n")
    print(open(os.path.join(a.out, "summary.tsv")).read(), end="")
    print(f"{n_lineages} lineages (group {a.group_prefix!r}: {len(group)} strains); proteins dropped for an internal stop: {dropped_stop}")

if __name__ == "__main__":
    main()
