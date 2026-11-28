#!/usr/bin/env python3
"""Fixture test for make_scaling_data.py. Run: python3 -I test_make_scaling_data.py"""
import os
import re
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_scaling_data as m


def test_shape_and_determinism():
    with tempfile.TemporaryDirectory() as d:
        a, b = os.path.join(d, "a"), os.path.join(d, "b")
        m.write(a, 12, 50, 3)
        m.write(b, 12, 50, 3)
        assert open(a + ".table").read() == open(b + ".table").read()
        lines = open(a + ".table").read().splitlines()
        assert sum(l.startswith("#ORGANISM") for l in lines) == 12
        data = [l for l in lines if not l.startswith("#")]
        assert len(data) == 12 and all(len(l.split("\t")[1]) == 50 for l in data)
        assert sum(l.startswith("#MAP\t") and l[5].isdigit() for l in lines) == 50
        nwk = open(a + ".nwk").read().strip()
        assert nwk.endswith(");") and sorted(re.findall(r"T\d{4}", nwk)) == [f"T{i:04d}" for i in range(12)]
        assert nwk.count("(") == 11  # a rooted binary tree on 12 tips


def test_sites_are_informative():
    taxa, rows, _ = m.build(30, 400, 1)
    present = [sum(rows[t][s] == ord("1") for t in taxa) for s in range(400)]
    assert max(present) > 1 and min(present) >= 0 and sum(p > 0 for p in present) > 300


def test_two_taxa_and_deep_tree():
    m.build(2, 5, 1)
    m.build(2000, 3, 1)  # a deep random tree must not hit the recursion limit in the walks


if __name__ == "__main__":
    test_shape_and_determinism(); test_sites_are_informative(); test_two_taxa_and_deep_tree(); print("ok")
