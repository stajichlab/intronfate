#!/usr/bin/env python3
"""Summarise bench.sh's timings.tsv: median wall time and peak memory per (tool, task), as a Markdown table.

Rows whose exit status is not 0 are reported, never averaged in. Usage: summarize.py timings.tsv
"""
import csv
import statistics
import sys


def summarise(rows):
    """rows: dicts with tool, task, wall_s, max_rss_kb, exit. Returns (table lines, failures)."""
    groups, failures = {}, []
    for r in rows:
        key = (r["tool"], r["task"])
        if r["exit"] != "0":
            failures.append(key)
            continue
        groups.setdefault(key, []).append((float(r["wall_s"]), int(r["max_rss_kb"])))
    lines = ["| tool | task | n | median wall (s) | min-max wall (s) | peak RSS (MB) |", "|---|---|---|---|---|---|"]
    for (tool, task), v in sorted(groups.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        walls = [w for w, _ in v]
        rss = max(m for _, m in v) / 1024
        lines.append(f"| {tool} | {task} | {len(v)} | {statistics.median(walls):.2f} | {min(walls):.2f}-{max(walls):.2f} | {rss:.0f} |")
    return lines, failures


def main(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    lines, failures = summarise(rows)
    print("\n".join(lines))
    if failures:
        print("\nFAILED runs (excluded): " + ", ".join(f"{t}/{k}" for t, k in sorted(set(failures))))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
