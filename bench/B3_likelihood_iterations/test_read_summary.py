#!/usr/bin/env python3
"""Fixture test for read_summary.py. Run: python3 -I test_read_summary.py"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import read_summary as r

with tempfile.TemporaryDirectory() as d:
    open(os.path.join(d, "summary_score.tsv"), "w").write("key\tvalue\nml_iterations\t33\nml_converged\tyes\n")
    assert r.read(d, "ml_iterations") == "33" and r.read(d, "ml_converged") == "yes" and r.read(d, "nope") is None
print("ok")
