#!/usr/bin/env python3
"""How single-copy are the OrthoFinder results? Classifies every orthogroup (OG, from Orthogroups.tsv) and every root-level
hierarchical orthogroup (HOG, Phylogenetic_Hierarchical_Orthogroups/N<k>.tsv, one file per internal node of the species tree) by
how many species it spans and whether any species has more than one gene in it.

The root node N0 has no HOG file of its own: OrthoFinder's Orthogroups.tsv is the root-level grouping, so it is reported as the
first row. An N<k> file only covers the species below that node (N1 here excludes the outgroup Spunctatus), so S for a HOG file is the
number of species that appear in it, not the number of columns.

Usage: orthogroup_summary.py <results/variant> [Refined_SOGs.tsv ...]
Classes (S = number of species in the data):
  strict single-copy      present in all S species, exactly one gene each
  single-copy, >= k of S  at most one gene per species, present in at least k species (k = 7, 5)
  multi-copy              at least one species has >= 2 genes (paralogs / duplications)
"""
import csv, sys, os, glob

res = sys.argv[1]
KS = [7, 5]

def cells(row, names):
    return [len([g for g in row[n].replace(", ", ",").split(",") if g]) if row.get(n, "") else 0 for n in names]

def classify(counts, S):
    present = sum(1 for c in counts if c > 0)
    if present == 0:
        return "empty", 0
    if max(counts) >= 2:
        return "multi-copy", present
    if present == S:
        return "strict single-copy", present
    return "single-copy, partial", present

def summarise(path, label, first_species_col, clade_only=False):
    with open(path, newline="") as f:
        r = csv.DictReader(f, delimiter="\t")
        names = r.fieldnames[first_species_col:]
        rows = list(r)
        if clade_only:      # species below this node = columns that have a gene somewhere in the file
            names = [n for n in names if any(row.get(n) for row in rows)]
        S = len(names)
        tally = {"strict single-copy": 0, "multi-copy": 0}
        partial = {k: 0 for k in KS}
        genes = 0; multi_genes = 0; total = 0; max_copies = []
        for row in rows:
            c = cells(row, names); total += 1; genes += sum(c)
            cls, present = classify(c, S)
            if cls == "strict single-copy": tally[cls] += 1
            elif cls == "multi-copy":
                tally[cls] += 1; multi_genes += sum(c)
                max_copies.append(max(c))
            elif cls == "single-copy, partial":
                for k in KS:
                    if present >= k: partial[k] += 1
        print(f"  {label:30s} groups {total:6d} genes {genes:6d} | strict single-copy {tally['strict single-copy']:5d} | "
              + " ".join(f"single-copy >= {k}/{S}: {tally['strict single-copy'] + partial[k]:5d}" for k in KS)
              + f" | multi-copy {tally['multi-copy']:5d} ({100*multi_genes/max(genes,1):.0f}% of grouped genes; median max copies "
              + (f"{sorted(max_copies)[len(max_copies)//2]}" if max_copies else "-") + ")")

print(f"== {res}")
og = os.path.join(res, "Orthogroups", "Orthogroups.tsv")
summarise(og, "orthogroups (OG, root)", 1)
hog_dir = os.path.join(res, "Phylogenetic_Hierarchical_Orthogroups")
for n in sorted(glob.glob(os.path.join(hog_dir, "N*.tsv")), key=lambda p: int(os.path.basename(p)[1:-4])):
    summarise(n, "HOG " + os.path.basename(n)[:-4], 3, clade_only=True)
for extra in sys.argv[2:]:                      # refined files from orthosynassign: SOG_ID then one column per species
    summarise(extra, "refined " + os.path.basename(extra), 1)
un = os.path.join(res, "Orthogroups", "Orthogroups_UnassignedGenes.tsv")
if os.path.exists(un):
    with open(un) as f:
        n = sum(1 for _ in f) - 1
    print(f"  genes in no orthogroup (unassigned singletons): {n}")
