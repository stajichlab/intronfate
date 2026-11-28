#!/usr/bin/env python3
"""Fixture test for make_smoke_table.py: deterministic, 6 taxa, one row per taxon, one #MAP line per site."""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from make_smoke_table import write  # noqa: E402


class SmokeTable(unittest.TestCase):
    def test_shape_and_determinism(self):
        with tempfile.TemporaryDirectory() as d:
            write(os.path.join(d, "a"), 30, 7)
            write(os.path.join(d, "b"), 30, 7)
            a, b = (open(os.path.join(d, x + ".table")).read() for x in "ab")
            self.assertEqual(a, b)
            lines = a.splitlines()
            self.assertEqual(sum(l.startswith("#MAP\t") and l.split("\t")[1] != "site" for l in lines), 30)
            rows = [l for l in lines if l and l[0] in "ABCDEF"]
            self.assertEqual(len(rows), 6)
            self.assertTrue(all(len(r.split("\t")[1]) == 30 for r in rows))


if __name__ == "__main__":
    unittest.main()
