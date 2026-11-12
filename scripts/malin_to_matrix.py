#!/usr/bin/env python3
"""Convert a MALIN intron table into this project's `matrix.tsv` (one file per alignment).

MALIN table layout (see the MALIN guide, section 2.4, and `Malin/test/Chromo23.table`):

* one row per taxon, one character per intron SITE (not per alignment column);
* characters: `0` no intron, `?` not observable, and for an intron present: `3` = phase 0 (the intron
  is at a codon boundary), `1` = phase 1, `2` = phase 2. The guide calls boundary introns "phase 3";
* `#MAP <site> <kog> <alipos>` lines place each site: `alipos` is a nucleotide position in the
  protein alignment named `<kog>`, with `alipos = 3 * codon_index0 + phase` (3 is between codon 1 and
  codon 2). That is this project's `(column, phase) = (alipos // 3, alipos % 3)`.

Because a site's column is only unique within its alignment, a table with `#MAP` lines is written as
one `<kog>.matrix.tsv` per alignment, into the output directory.

Tables without `#MAP` lines (for example `rogozin_introns_table`, which uses 0/1 only) carry no phase
or alignment: every site becomes `column=<site index>, phase=0` in a single file marked
`phase_mode=ignore`, and nothing in it should be read as phase information.

Usage:  scripts/malin_to_matrix.py <malin.table> <out> [alignment-name]
        <out> is a directory when the table has #MAP lines, a file otherwise.

`Malin/` is third-party data and stays out of version control; reference it by path.
"""
from __future__ import annotations

import os
import sys
from collections import defaultdict

CODE_TO_PHASE = {"3": 0, "1": 1, "2": 2}


def read_table(path: str):
    taxa: list[str] = []
    rows: list[str] = []
    site_map: dict[int, tuple[str, int]] = {}
    with open(path, errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("#MAP"):
                parts = line.split("\t") if "\t" in line else line.split()
                if len(parts) >= 4 and parts[1].isdigit() and parts[3].split()[0].isdigit():
                    site_map[int(parts[1])] = (parts[2], int(parts[3].split()[0]))
                continue
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) >= 2:
                taxa.append(fields[0])
                rows.append(fields[1])
    return taxa, rows, site_map


def write_matrix(path: str, name: str, taxa: list[str], sites, phase_mode: str, n_cols: int) -> None:
    """`sites`: list of (column, phase, per-taxon cell string over `1`, `0`, `?`)."""
    unannotated = [t for i, t in enumerate(taxa) if sites and all(s[2][i] == "?" for s in sites)]
    with open(path, "w") as out:
        out.write(
            f"# format=2 alignment={name} alen={n_cols} phase_mode={phase_mode} flank=0 "
            f"masked_present=0 taxa={len(taxa)} sites={len(sites)} "
            f"unannotated={','.join(unannotated) if unannotated else '-'} residues=-\n"
        )
        out.write("site\tcolumn\tphase\tn_taxa\tn_observable\tpattern\t" + "\t".join(taxa) + "\n")
        for col, ph, cells in sites:
            out.write(
                f"{col}:{ph}\t{col}\t{ph}\t{cells.count('1')}\t{len(cells) - cells.count('?')}\t{cells}\t"
                + "\t".join(cells)
                + "\n"
            )


def main(path: str, out: str, name: str = "") -> int:
    taxa, rows, site_map = read_table(path)
    if not taxa:
        print(f"{path}: no data rows found", file=sys.stderr)
        return 1
    if len({len(r) for r in rows}) != 1:
        print(f"{path}: rows differ in length", file=sys.stderr)
        return 1
    n_sites = len(rows[0])
    bad = {c for r in rows for c in r} - set("0123?")
    if bad:
        print(f"{path}: unexpected characters {sorted(bad)}", file=sys.stderr)
        return 1

    def cells_for(site: int, code_to_one) -> str:
        return "".join(code_to_one(r[site]) for r in rows)

    if site_map:
        missing = [s for s in range(n_sites) if s not in site_map]
        if missing:
            print(f"{path}: {len(missing)} sites have no #MAP line (first: {missing[0]})", file=sys.stderr)
            return 1
        per_aln: dict[str, list] = defaultdict(list)
        for site in range(n_sites):
            aln, alipos = site_map[site]
            col, phase = divmod(alipos, 3)
            # The digit repeats the phase; a disagreement means the table or the map is corrupt.
            for r in rows:
                c = r[site]
                if c in CODE_TO_PHASE and CODE_TO_PHASE[c] != phase:
                    print(
                        f"{path}: site {site} has code {c} (phase {CODE_TO_PHASE[c]}) but alipos "
                        f"{alipos} implies phase {phase}",
                        file=sys.stderr,
                    )
                    return 1
            cells = cells_for(site, lambda c: "1" if c in CODE_TO_PHASE else c)
            per_aln[aln].append((col, phase, cells))
        os.makedirs(out, exist_ok=True)
        for aln, sites in per_aln.items():
            sites.sort(key=lambda s: (s[0], s[1]))
            write_matrix(os.path.join(out, f"{aln}.matrix.tsv"), aln, taxa, sites, "strict", 0)
        n_q = sum(r.count("?") for r in rows)
        print(
            f"{len(taxa)} taxa x {n_sites} sites -> {len(per_aln)} alignments in {out}; "
            f"{n_q} per-cell '?'"
        )
        return 0

    # No #MAP: binary table, phase unknown.
    if any(c in "23" for r in rows for c in r):
        print(
            f"{path}: contains phase codes 2/3 but no #MAP lines, so positions cannot be recovered",
            file=sys.stderr,
        )
        return 1
    sites = [(s, 0, cells_for(s, lambda c: c)) for s in range(n_sites)]
    sites = [s for s in sites if "1" in s[2]]
    write_matrix(out, name or "malin", taxa, sites, "ignore", n_sites)
    print(f"{len(taxa)} taxa x {n_sites} sites -> {len(sites)} sites with at least one intron; phase unknown")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(*sys.argv[1:4]))
