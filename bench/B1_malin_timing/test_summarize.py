#!/usr/bin/env python3
"""Fixture test for summarize.py: medians, peak memory, and failed runs kept out of the statistics."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from summarize import summarise  # noqa: E402


def row(tool, task, wall, rss, ex="0"):
    return {"tool": tool, "task": task, "rep": "1", "wall_s": wall, "max_rss_kb": rss, "exit": ex}


class Summarise(unittest.TestCase):
    def test_median_and_peak(self):
        lines, fails = summarise([row("a", "dollo", "1.0", "1024"), row("a", "dollo", "3.0", "2048"), row("a", "dollo", "2.0", "512")])
        self.assertEqual(fails, [])
        self.assertIn("| a | dollo | 3 | 2.00 | 1.00-3.00 | 2 |", lines)

    def test_failures_are_excluded_and_reported(self):
        lines, fails = summarise([row("m", "likelihood", "9.0", "100"), row("m", "likelihood", "0.1", "1", ex="1")])
        self.assertEqual(fails, [("m", "likelihood")])
        self.assertIn("| m | likelihood | 1 | 9.00 | 9.00-9.00 | 0 |", lines)


if __name__ == "__main__":
    unittest.main()
