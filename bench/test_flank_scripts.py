#!/usr/bin/env python3
"""Fixture test for F1/gap_distance.py and F2/missing_cells.py: tiny alignments whose gap positions were worked out by hand.

The matrix `column` is 0-based (the same as the row index); a script that treats it as 1-based fails every assertion here.
Run: python3 -I bench/test_flank_scripts.py     Standard library only.
"""
import os, subprocess, sys, tempfile, unittest

HERE = os.path.dirname(os.path.abspath(__file__))
GAP = os.path.join(HERE, "F1_flank_tuning", "gap_distance.py")
MISS = os.path.join(HERE, "F2_column_ambiguity", "missing_cells.py")

#            0123456789
ROWS = {"px": "MKVLAAGTSS", "py": "MK--AAGTSS", "pz": "MK--------"}
TAXA = {"px": "X", "py": "Y", "pz": "Z"}


def write_fixture(d, taxa, sites):
    os.makedirs(os.path.join(d, "aln"))
    with open(os.path.join(d, "aln", "FAM.aln.fasta"), "w") as fh:
        for k, v in ROWS.items():
            fh.write(f">{k}\n{v}\n")
    with open(os.path.join(d, "taxon_map.tsv"), "w") as fh:
        for k, v in TAXA.items():
            fh.write(f"{k}\t{v}\n")
    with open(os.path.join(d, "matrix.tsv"), "w") as fh:
        fh.write("# format=2 alignment=combined flank=1\n")
        fh.write("site\tcolumn\tphase\tfamily\tfamily_column\tn_taxa\tn_observable\tpattern\t" + "\t".join(taxa) + "\n")
        for col, phase, cells in sites:
            fh.write(f"FAM.aln:{col}:{phase}\t{col}\t{phase}\tFAM.aln\t{col}\t{len(taxa)}\t{len(taxa)}\t{''.join(cells)}\t" + "\t".join(cells) + "\n")


def run(script, *args):
    r = subprocess.run([sys.executable, "-I", script, *args], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r


class GapDistance(unittest.TestCase):
    def test_distances_and_absences(self):
        with tempfile.TemporaryDirectory() as d:
            # site A (phase 1, column 5): X and Y carry it.  Y has residues at 4 and 6.. but a gap at 3: left run 1, right run 4 -> d 1;
            # X is gap-free: left run 5 (columns 4..0), right run 4 (6..9) -> d 4.
            # (the index check counts the two carriers of site A; X at site B is skipped, its only other taxon being absent)
            # site B (phase 1, column 6): X carries it, Y is absent: left run from 5 is 2 (5, 4), right run from 7 is 3 -> d 2.
            write_fixture(d, ["X", "Y"], [(5, 1, ["1", "1"]), (6, 1, ["1", "0"])])
            out = os.path.join(d, "out.tsv")
            r = run(GAP, "--matrix", os.path.join(d, "matrix.tsv"), "--aln-dir", os.path.join(d, "aln"),
                    "--taxon-map", os.path.join(d, "taxon_map.tsv"), "--min-others", "1", "--out", out)
            self.assertIn("index check: 0 of 2 present phase-1/2 cells sit on a gap", r.stderr)
            rows = [l.split("\t") for l in open(out).read().splitlines()[1:]]
            got = {(p, int(dist)): (int(n), int(z)) for p, dist, n, z, _ in rows}
            self.assertEqual(got[("1", 1)], (1, 0))
            self.assertEqual(got[("1", 2)], (1, 1))
            self.assertEqual(got[("1", 4)], (1, 0))
            self.assertEqual(sum(n for (p, _), (n, _) in got.items() if p == "1"), 3)

    def test_a_one_based_reading_would_be_caught(self):
        with tempfile.TemporaryDirectory() as d:
            # Y is a gap at column 3: a present phase-1/2 cell there is reported by the index check
            write_fixture(d, ["X", "Y"], [(3, 1, ["1", "1"])])
            r = run(GAP, "--matrix", os.path.join(d, "matrix.tsv"), "--aln-dir", os.path.join(d, "aln"),
                    "--taxon-map", os.path.join(d, "taxon_map.tsv"), "--min-others", "1", "--out", os.path.join(d, "o.tsv"))
            self.assertIn("index check: 1 of 2", r.stderr)


class MissingCells(unittest.TestCase):
    def test_classes(self):
        with tempfile.TemporaryDirectory() as d:
            # site C (phase 1, column 3): Y gap run 2..3 with residues each side -> short; Z gap run 2..9 reaches the row end -> edge.
            # site D (phase 1, column 6): Z gap again, edge.
            write_fixture(d, ["X", "Y", "Z"], [(3, 1, ["1", "?", "?"]), (6, 1, ["1", "1", "?"])])
            r = run(MISS, "--matrix", os.path.join(d, "matrix.tsv"), "--aln-dir", os.path.join(d, "aln"),
                    "--taxon-map", os.path.join(d, "taxon_map.tsv"), "--flank", "1")
            lines = r.stdout.splitlines()
            self.assertTrue(lines[0].startswith("cells 6  missing 3"), lines[0])
            kinds = {l.split()[0]: int(l.split()[1]) for l in lines[1:6]}
            self.assertEqual(kinds, {"short": 1, "long": 0, "edge": 2, "other": 0, "no_row": 0})
            self.assertIn("sites with a carrier and at least one short-gap missing cell: 1", r.stdout)

    def test_a_missing_cell_without_a_gap_is_other(self):
        with tempfile.TemporaryDirectory() as d:
            write_fixture(d, ["X", "Y"], [(7, 1, ["1", "?"])])  # Y has residues at 6..8: not a gap
            r = run(MISS, "--matrix", os.path.join(d, "matrix.tsv"), "--aln-dir", os.path.join(d, "aln"),
                    "--taxon-map", os.path.join(d, "taxon_map.tsv"), "--flank", "1")
            self.assertIn("other", r.stdout)
            self.assertRegex(r.stdout, r"other\s+1 ")


if __name__ == "__main__":
    unittest.main()
