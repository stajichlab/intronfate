#!/usr/bin/env python3
"""Turn the ```sh blocks of a Markdown how-to into one bash script, so the documentation can be tested exactly as a reader runs it.

Blocks are taken in order, starting at the first heading that begins with --from (the setup section before it is usually done by
the test harness itself). --sub OLD NEW replaces a line that equals OLD (after stripping) with NEW, e.g. to pick the other example
or to use a container; every OLD must be found, so a test cannot silently drift from the document.

Only blocks fenced exactly ```sh are taken; alternatives the reader may choose instead (fenced ```bash, ```text, ...) are left out.
--trace prints the time before each block to stderr, for per-step timings.

Usage: doc_blocks.py docs/HOWTO.md --from "## 1." [--sub 'OLD LINE' 'NEW LINE' ...] [--trace] > howto.sh
"""
import argparse, sys


def blocks(text, start):
    on = start is None
    inside = False
    cur, out = [], []
    for line in text.splitlines():
        if not inside and line.startswith("#") and start is not None and line.startswith(start):
            on = True
        if line.startswith("```"):
            if inside:
                inside = False
                if keep:
                    out.append(cur)
                cur = []
            else:
                inside = True
                keep = on and line.strip() == "```sh"
            continue
        if inside:
            cur.append(line)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("doc")
    ap.add_argument("--from", dest="start", help="heading prefix to start at, e.g. '## 1.'")
    ap.add_argument("--sub", nargs=2, action="append", default=[], metavar=("OLD", "NEW"))
    ap.add_argument("--trace", action="store_true", help="print the time before each block to stderr")
    a = ap.parse_args(argv)
    bl = blocks(open(a.doc).read(), a.start)
    if not bl:
        sys.exit(f"no ```sh blocks found in {a.doc}" + (f" after a heading starting {a.start!r}" if a.start else ""))
    used = set()
    print("set -euo pipefail")
    for i, b in enumerate(bl, 1):
        print(f"\n# ---- block {i}")
        if a.trace:
            print(f'echo "[doc block {i}/{len(bl)}] $(date +%s) $(date +%T)" >&2')
        for line in b:
            for old, new in a.sub:
                if line.strip() == old.strip():
                    line = new; used.add(old)
            print(line)
    missing = [o for o, _ in a.sub if o not in used]
    if missing:
        sys.exit(f"--sub line(s) not found in {a.doc}: {missing}")


if __name__ == "__main__":
    main()
