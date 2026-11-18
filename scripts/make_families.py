#!/usr/bin/env python3
"""Ortholog families for intronfate from OrthoFinder results, in tiers, with at most one protein per species.

User-facing version of bench/V6_orthology_nine_species/make_families.py (which stays as the record of the V6 run). The difference:
proteins are read from the folder OrthoFinder was given (--proteomes, <label>.faa), so no particular data layout is needed.

Tiers (S = number of species, MIN = least species a family may span; a species with no protein in a family is absent from it, which
intronfate reads as missing data and never as an intron loss):
  A_strict      orthogroups (Orthogroups/Orthogroups.tsv) single-copy in all S species
  B_partial     orthogroups single-copy in MIN..S-1 species
  C_hog         from multi-copy orthogroups: hierarchical orthogroups (Phylogenetic_Hierarchical_Orthogroups/N<k>.tsv) single-copy in
                >= MIN species; only the shallowest node that yields such HOGs is used for a given orthogroup. A HOG file only covers the
                species below its node, so with a high MIN only N1 can qualify
  D_synteny     orthoSynAssign synteny-refined groups (--refined) single-copy in >= MIN species, if given
Paralog sets that do not resolve to single copy are left out, never guessed at.

Writes <out>/<tier>/families/<id>.faa, <out>/<tier>/taxon_map.tsv (protein<TAB>label, for `intronfate --taxon-map`) and
<out>/summary.tsv. A protein with an internal stop is dropped from its family (the aligned row would no longer equal the protein, and
intronfate would refuse it); a family left with fewer than MIN proteins is dropped.

Usage: make_families.py --of orthofinder/Results_x --proteomes proteomes --min 4 --out families [--refined Refined_SOGs.tsv]
"""
import argparse, csv, glob, gzip, os, sys

EXTS = (".faa", ".fa", ".fasta", ".faa.gz", ".fa.gz", ".fasta.gz")


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as f:
        for line in f:
            if line.startswith(">"):
                if name: seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs


def proteome(folder, label):
    for e in EXTS:
        p = os.path.join(folder, label + e)
        if os.path.exists(p):
            return p
    sys.exit(f"no proteome for {label} in {folder} (looked for {label}{{{','.join(EXTS)}}})")


def split(cell):
    return [g for g in cell.replace(", ", ",").split(",") if g] if cell else []


def groups(path, first_species_col, id_col=0, extra_col=None):
    """Yield (group id, parent OG id or None, {species: [genes]})."""
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        names = r.fieldnames[first_species_col:]
        for row in r:
            yield row[r.fieldnames[id_col]], (row[extra_col] if extra_col else None), {n: split(row[n]) for n in names if row.get(n)}


def single_copy(m):
    return all(len(v) == 1 for v in m.values())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--of", required=True, help="OrthoFinder Results_<name> folder")
    ap.add_argument("--proteomes", required=True, help="folder with <label>.faa (the files given to OrthoFinder)")
    ap.add_argument("--min", type=int, required=True, help="least species a family may span (tiers B and C)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--refined", help="orthoSynAssign Refined_SOGs.tsv (tier D)")
    a = ap.parse_args(argv)

    og_table = os.path.join(a.of, "Orthogroups", "Orthogroups.tsv")
    with open(og_table) as f:
        S = len(f.readline().rstrip("\n").split("\t")) - 1
    fams = {"A_strict": {}, "B_partial": {}, "C_hog": {}, "D_synteny": {}}
    og_multi = set()
    for og, _, m in groups(og_table, 1):
        if single_copy(m):
            if len(m) == S: fams["A_strict"][og] = m
            elif len(m) >= a.min: fams["B_partial"][og] = m
        else:
            og_multi.add(og)
    hogs = sorted(glob.glob(os.path.join(a.of, "Phylogenetic_Hierarchical_Orthogroups", "N*.tsv")), key=lambda p: int(os.path.basename(p)[1:-4]))
    done_og = set()
    for p in hogs:
        found = set()
        for hog, og, m in groups(p, 3, 0, "OG"):
            if og in og_multi and og not in done_og and len(m) >= a.min and single_copy(m):
                fams["C_hog"][hog] = m; found.add(og)
        done_og |= found
    if a.refined:
        for sog, _, m in groups(a.refined, 1):
            if len(m) >= a.min and single_copy(m): fams["D_synteny"][sog] = m

    labels = sorted({l for t in fams.values() for m in t.values() for l in m})
    prot = {l: read_fasta(proteome(a.proteomes, l)) for l in labels}
    rows = []; dropped_stop = 0
    for tier, table in fams.items():
        tdir = os.path.join(a.out, tier); os.makedirs(os.path.join(tdir, "families"), exist_ok=True)
        n_fam = n_seq = 0; per = {}
        with open(os.path.join(tdir, "taxon_map.tsv"), "w") as tm:
            for fid, m in sorted(table.items()):
                keep = []
                for lab, genes in sorted(m.items()):
                    s = prot[lab].get(genes[0])
                    if s is None: sys.exit(f"{genes[0]} ({lab}) is not in {proteome(a.proteomes, lab)}")
                    if "*" in s.rstrip("*"): dropped_stop += 1; continue
                    keep.append((lab, genes[0], s.rstrip("*")))
                if len(keep) < a.min: continue
                with open(os.path.join(tdir, "families", f"{fid}.faa"), "w") as out:
                    for lab, g, s in keep:
                        out.write(f">{g}\n{s}\n"); tm.write(f"{g}\t{lab}\n")
                n_fam += 1; n_seq += len(keep); per[len(keep)] = per.get(len(keep), 0) + 1
        rows.append((tier, n_fam, n_seq, ";".join(f"{k}:{v}" for k, v in sorted(per.items()))))
    with open(os.path.join(a.out, "summary.tsv"), "w") as f:
        f.write("tier\tfamilies\tsequences\tspecies_per_family\n")
        for r in rows: f.write("\t".join(map(str, r)) + "\n")
    print(open(os.path.join(a.out, "summary.tsv")).read(), end="")
    print(f"S={S} MIN={a.min}; proteins dropped for an internal stop: {dropped_stop}")


if __name__ == "__main__":
    main()
