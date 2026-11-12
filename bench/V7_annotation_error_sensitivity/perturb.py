#!/usr/bin/env python3
"""Perturb an intronfate `matrix.tsv` the way annotation error would, given each taxon's measured sensitivity (Sn) and precision (Pr).

Work on ONE matrix that holds every family (the `combined/matrix.tsv` of `analyze --combine`), because a missed-intron budget belongs to the genome, not to
a family: a family where a taxon has no candidate cell must not use up, or discard, the budget of another.

The observed matrix O is what the annotation says; we do not know the truth T, so we draw plausible T from O by undoing the errors:
  * missed introns (false negatives): a taxon with P introns in the matrix, precision Pr and sensitivity Sn is missing about  M = P*Pr*(1-Sn)/Sn  introns
    (true positives P*Pr, and FN = TP*(1-Sn)/Sn). A miss matters to the Dollo score only where it hides a loss, i.e. where Dollo puts a loss on the taxon's own
    terminal branch. The candidate cells are taken from the Dollo scorer's own output (`lost_in` of dollo_sites.tsv), so they follow its exact rule, including
    the treatment of `?` cells and ties. Setting a candidate to '1' is "the loss was really a missed intron".
    `ceiling` fills min(M, candidates) of them (every miss hides a terminal loss: the V4 worst case); `typical` scales M by f, the fraction of the taxon's own
    introns at sites with >= 2 carriers (missed introns are assumed to be distributed like annotated ones). Cells are drawn uniformly from the pooled candidates.
  * false introns (false positives): every '1' becomes '0' with probability 1 - Pr.
Modes: ceiling_fn / typical_fn (misses only), fp (false introns only), ceiling (both). '?' cells are never changed. Sites left with no carrier are dropped;
n_taxa, n_observable and pattern are recomputed and every column before `n_taxa` (site, column, phase, family, ...) is kept.
"""

def read_matrix(path):
    with open(path) as f:
        lines = f.read().rstrip("\n").split("\n")
    meta = lines[0] if lines[0].startswith("#") else None
    hdr = lines[1 if meta else 0].split("\t")
    rows = [l.split("\t") for l in lines[(2 if meta else 1):] if l]
    return meta, hdr, rows

def write_matrix(path, meta, hdr, rows):
    with open(path, "w") as f:
        if meta: f.write(meta + "\n")
        f.write("\t".join(hdr) + "\n")
        for r in rows: f.write("\t".join(r) + "\n")

def stochastic_round(x, rng):
    k = int(x)
    return k + (1 if rng.random() < x - k else 0)

def candidates_from_dollo(dollo_sites_path, site_index, taxa):
    """{taxon: [row index]} of the sites where the Dollo scorer put a loss on that taxon's own terminal branch.
    `lost_in` lists the branches that lost the site, separated by ';', each branch as the comma-joined tips below it (a terminal branch is one tip name)."""
    cand = {t: [] for t in taxa}
    with open(dollo_sites_path) as f:
        head = f.readline().rstrip("\n").split("\t"); si, li = head.index("site"), head.index("lost_in")
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) <= li or not c[li]: continue
            i = site_index.get(c[si])
            if i is None: continue
            for br in c[li].split(";"):
                if br in cand: cand[br].append(i)
    return cand

def perturb(hdr, rows, sn, pr, mode, rng, candidates):
    """sn, pr: {taxon: fraction}. candidates: {taxon: [row indices]} (see candidates_from_dollo). Returns new rows; the input is not modified."""
    first = hdr.index("pattern") + 1
    taxa = hdr[first:]
    cells = [list(r[first:]) for r in rows]
    carriers = [sum(1 for x in c if x == "1") for c in cells]
    orig = [list(c) for c in cells]
    for j, t in enumerate(taxa):
        present = [i for i in range(len(rows)) if orig[i][j] == "1"]
        s, p = sn.get(t, 1.0), pr.get(t, 1.0)
        if present and 0 < s < 1 and mode != "fp":
            missed = len(present) * p * (1 - s) / s
            if mode.startswith("typical"):
                missed *= sum(1 for i in present if carriers[i] >= 2) / len(present)
            pool = [i for i in candidates.get(t, []) if orig[i][j] == "0"]
            for i in rng.sample(pool, min(stochastic_round(missed, rng), len(pool))): cells[i][j] = "1"
        if p < 1 and not mode.endswith("_fn"):
            for i in present:
                if rng.random() < 1 - p: cells[i][j] = "0"
    out = []
    for r, c in zip(rows, cells):
        n1 = sum(1 for x in c if x == "1")
        if n1 == 0: continue
        out.append(r[:first - 3] + [str(n1), str(sum(1 for x in c if x in "01")), "".join(c)] + c)
    return out
