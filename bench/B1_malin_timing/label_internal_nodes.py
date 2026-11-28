#!/usr/bin/env python3
"""Copy a Newick tree, naming every unnamed internal node n1, n2, ... (left to right, post-order).

MALIN's `.rates` output names nodes by their Newick label, and the malin_rates.rs comparison finds nodes by name, so the
tree both tools read in the B1 likelihood check must label its internal nodes. Usage: label_internal_nodes.py IN OUT
"""
import re
import sys


def label(newick):
    count = iter(range(1, 10**6))
    return re.sub(r"\)(?=[,);:])", lambda _m: f")n{next(count)}", newick)


if __name__ == "__main__":
    src = open(sys.argv[1]).read().strip()
    open(sys.argv[2], "w").write(label(src) + "\n")
