#!/usr/bin/env python3
"""Print one value from <out-dir>/summary_score.tsv (`key<TAB>value` lines). Usage: read_summary.py OUT_DIR KEY. Exit 1 if absent."""
import os
import sys


def read(out_dir, key):
    with open(os.path.join(out_dir, "summary_score.tsv")) as fh:
        for line in fh:
            k, _, v = line.rstrip("\n").partition("\t")
            if k == key:
                return v
    return None


if __name__ == "__main__":
    v = read(sys.argv[1], sys.argv[2])
    if v is None:
        sys.exit(1)
    print(v)
