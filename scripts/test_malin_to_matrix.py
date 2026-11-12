#!/usr/bin/env python3
"""Fixture tests for malin_to_matrix.py (synthetic tables in tests/data, no third-party data)."""
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "malin_to_matrix.py")
DATA = os.path.join(HERE, "..", "tests", "data")


def run(*args):
    return subprocess.run([sys.executable, "-I", SCRIPT, *args], capture_output=True, text=True)


def rows(path):
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f]
    return lines[0], lines[1], [l.split("\t") for l in lines[2:]]


class Converter(unittest.TestCase):
    def test_phase_codes_and_map_give_column_and_phase(self):
        with tempfile.TemporaryDirectory() as d:
            r = run(os.path.join(DATA, "malin_tiny.table"), d)
            self.assertEqual(r.returncode, 0, r.stderr)
            meta, header, k1 = rows(os.path.join(d, "K1.matrix.tsv"))
            self.assertIn("format=2", meta)
            self.assertEqual(header.split("\t")[:6], ["site", "column", "phase", "n_taxa", "n_observable", "pattern"])
            # alipos 3 -> column 1, phase 0 (code 3); alipos 7 -> column 2, phase 1 (code 1)
            self.assertEqual([(x[1], x[2], x[5]) for x in k1], [("1", "0", "11"), ("2", "1", "10")])
            _, _, k2 = rows(os.path.join(d, "K2.matrix.tsv"))
            # alipos 8 -> column 2, phase 2 (code 2); A is '?', B present
            self.assertEqual([(x[1], x[2], x[5]) for x in k2], [("2", "2", "?1")])

    def test_phase_code_that_contradicts_alipos_is_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            bad = os.path.join(d, "bad.table")
            with open(bad, "w") as f:
                f.write("A\t1\nB\t0\n#MAP\t0\tK1\t3\n")  # code 1 = phase 1, but alipos 3 is phase 0
            r = run(bad, os.path.join(d, "out"))
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("implies phase", r.stderr)

    def test_binary_table_has_unknown_phase(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "m.tsv")
            r = run(os.path.join(DATA, "malin_binary.table"), out)
            self.assertEqual(r.returncode, 0, r.stderr)
            meta, _, body = rows(out)
            self.assertIn("phase_mode=ignore", meta)
            self.assertEqual([x[5] for x in body], ["10", "10", "01"])


if __name__ == "__main__":
    unittest.main()
