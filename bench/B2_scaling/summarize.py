#!/usr/bin/env python3
"""Summarise scaling.sh's timings.tsv: one Markdown table per task, rows = taxa, columns = sites, cells = median wall s (peak RSS MB).

Runs whose exit status is not 0 are shown as `fail<code>` (124 = timed out), never averaged in. Usage: summarize.py timings.tsv
"""
import csv
import statistics
import sys


def cells(rows):
    """rows: dicts with task, taxa, sites, wall_s, max_rss_kb, exit -> {task: {(taxa, sites): text}}."""
    groups = {}
    for r in rows:
        groups.setdefault((r["task"], int(r["taxa"]), int(r["sites"])), []).append(r)
    out = {}
    for (task, t, s), g in groups.items():
        bad = [r["exit"] for r in g if r["exit"] != "0"]
        if bad:
            text = "fail" + bad[0]
        else:
            w = statistics.median(float(r["wall_s"]) for r in g)
            m = max(int(r["max_rss_kb"]) for r in g) / 1024
            text = f"{w:.2f} s ({m:.0f} MB)"
        out.setdefault(task, {})[(t, s)] = text
    return out


def tables(rows):
    lines = []
    for task, c in sorted(cells(rows).items()):
        taxa = sorted({t for t, _ in c}); sites = sorted({s for _, s in c})
        lines += [f"### {task}: median wall time (peak RSS), taxa down, sites across", "",
                  "| taxa \\ sites | " + " | ".join(str(s) for s in sites) + " |", "|---|" + "---|" * len(sites)]
        lines += [f"| {t} | " + " | ".join(c.get((t, s), "-") for s in sites) + " |" for t in taxa]
        lines.append("")
    return lines


if __name__ == "__main__":
    with open(sys.argv[1], newline="") as fh:
        print("\n".join(tables(list(csv.DictReader(fh, delimiter="\t")))))
