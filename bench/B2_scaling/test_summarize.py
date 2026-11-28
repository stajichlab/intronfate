#!/usr/bin/env python3
"""Fixture test for summarize.py. Run: python3 -I test_summarize.py"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import summarize as s


def row(task, t, n, w, rss, ex="0"):
    return dict(task=task, taxa=str(t), sites=str(n), rep="1", wall_s=str(w), max_rss_kb=str(rss), exit=ex)


def test_median_and_failures():
    rows = [row("dollo", 10, 100, 1.0, 2048), row("dollo", 10, 100, 3.0, 4096), row("dollo", 10, 100, 2.0, 1024),
            row("dollo", 100, 100, 9.0, 1024, "124"), row("likelihood", 10, 100, 5.0, 1024)]
    c = s.cells(rows)
    assert c["dollo"][(10, 100)] == "2.00 s (4 MB)"
    assert c["dollo"][(100, 100)] == "fail124"
    assert c["likelihood"][(10, 100)] == "5.00 s (1 MB)"
    text = "\n".join(s.tables(rows))
    assert "### dollo" in text and "| 100 | fail124 |" in text


if __name__ == "__main__":
    test_median_and_failures(); print("ok")
