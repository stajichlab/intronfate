#!/usr/bin/env python3
"""Summarise a V1 validation run and compare it with the figures quoted in docs/SPEC.md.

Usage: summarize.py <results-dir> <binary-table> <map-table>
Writes <results-dir>/{rogozin_branches.tsv,chromo23_branches.tsv}; prints a report. Never asserts the
spec's numbers: it prints observed next to claimed, because the point is to find out which is right.
"""
import collections
import csv
import glob
import os
import sys

out, binary_table, map_table = sys.argv[1:4]

# Figures quoted in docs/SPEC.md before they had been re-derived (all produced locally by an earlier session).
CLAIMED = {
    "rogozin_sites": 7221, "rogozin_dollo_losses": 2770, "rogozin_incompatible": 1619, "rogozin_mp_steps": 9007,
    "chromo23_sites": 9256, "alphabet": {"0": 160418, "?": 33102, "3": 10493, "2": 4568, "1": 4307},
}
# MALIN guide p. 27, per-branch Dollo losses on the Rogozin table (visible rows only; they sum to 2,719).
GUIDE_P27 = {"Dm": 156, "Ag": 194, "Ce": 411, "Hs": 112, "Sp": 796, "At": 46, "Pf": 0, "Diptera": 436, "Ecdysozoa": 568}


def kv(path):
    with open(path) as f:
        return dict(line.rstrip("\n").split("\t", 1) for line in f if "\t" in line)


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"))


def branches(score_dirs):
    tot = collections.defaultdict(lambda: [0, 0, ""])
    for d in score_dirs:
        for r in read_tsv(os.path.join(d, "dollo_branches.tsv")):
            t = tot[r["tips"]]
            t[0] += int(r["gains"]); t[1] += int(r["losses"]); t[2] = r["branch"]
    return tot


def write_branches(path, tot):
    with open(path, "w") as f:
        f.write("branch\tgains\tlosses\ttips\n")
        for tips, (g, l, name) in sorted(tot.items(), key=lambda kv: -kv[1][1]):
            f.write(f"{name}\t{g}\t{l}\t{tips}\n")


print("== timings (s)")
for line in open(os.path.join(out, "timings.tsv")).read().splitlines()[1:]:
    print("  " + line.replace("\t", "  "))

# ---- Rogozin ----
r = kv(glob.glob(os.path.join(out, "rogozin/score/summary_score.tsv"))[0])
print("\n== Rogozin 2003 table (binary, no #MAP)")
for label, key, claim in [("sites", "sites", "rogozin_sites"), ("Dollo losses", "dollo_losses", "rogozin_dollo_losses"),
                          ("sites needing a loss", "dollo_needing_losses", "rogozin_incompatible"), ("MP steps", "mp_steps", "rogozin_mp_steps")]:
    print(f"  {label:22s} observed {r[key]:>8}   spec claimed {CLAIMED[claim]:>8}   {'OK' if int(r[key]) == CLAIMED[claim] else 'DIFFERS'}")
rb = branches([os.path.join(out, "rogozin/score")])
write_branches(os.path.join(out, "rogozin_branches.tsv"), rb)
print("  per-branch losses (top 12) vs MALIN guide p.27 where the branch name matches a tip or label:")
for tips, (g, l, name) in sorted(rb.items(), key=lambda kv: -kv[1][1])[:12]:
    ref = GUIDE_P27.get(name)
    print(f"    {name:24s} {l:6d}   guide: {ref if ref is not None else '-'}")
print(f"  sum over all branches: {sum(v[1] for v in rb.values())}   (guide's visible rows sum to {sum(GUIDE_P27.values())})")

# ---- Chromo23 ----
print("\n== Chromo23 (phase codes + #MAP)")
rows = [l.split() for l in open(map_table, errors="replace") if l.strip() and not l.startswith("#")]
alpha = collections.Counter(c for r_ in rows for c in r_[1])
print(f"  taxa {len(rows)}  sites {len(rows[0][1])}  (spec claimed {CLAIMED['chromo23_sites']} sites)")
print("  alphabet  observed:", dict(sorted(alpha.items())))
print("  alphabet  claimed :", CLAIMED["alphabet"], "->", "OK" if dict(alpha) == CLAIMED["alphabet"] else "DIFFERS")
score_dirs = sorted(glob.glob(os.path.join(out, "chromo23/score/*/")))
tot_g = tot_l = tot_mp = tot_sites = 0
for d in score_dirs:
    s = kv(os.path.join(d, "summary_score.tsv"))
    tot_g += int(s.get("dollo_gains", 0)); tot_l += int(s.get("dollo_losses", 0))
    tot_mp += int(s.get("mp_steps", 0)); tot_sites += int(s.get("characters_scored", 0))
print(f"  alignments {len(score_dirs)}  sites scored {tot_sites}  Dollo gains {tot_g} losses {tot_l}  MP steps {tot_mp}")
write_branches(os.path.join(out, "chromo23_branches.tsv"), branches(score_dirs))

# ---- round trip: canonical (alignment, alipos) -> row string, order independent ----
def canonical(path):
    taxa, cells, maps = [], [], {}
    for l in open(path, errors="replace"):
        if l.startswith("#MAP"):
            p = l.split()
            if len(p) >= 4 and p[1].isdigit():
                maps[int(p[1])] = (p[2], int(p[3]))
        elif l.strip() and not l.startswith("#"):
            taxa.append(l.split()[0]); cells.append(l.split()[1])
    return {maps[s]: tuple(c[s] for c in cells) for s in range(len(cells[0]))} if cells else {}, taxa

a, ta = canonical(map_table)
b, tb = canonical(os.path.join(out, "chromo23/roundtrip.table"))
print("\n== round trip (import -> per-alignment matrices -> export)")
print(f"  taxa identical: {ta == tb}   sites: {len(a)} -> {len(b)}   every (alignment, alipos) row identical: {a == b}")
if a != b:
    sys.exit("ROUND TRIP DIFFERS")
