#!/usr/bin/env python3
"""Mutation check for scripts/test_doc_samples.py: proves the doc checks still fail when a document is wrong.

`test_doc_samples.py` compares the prose pages against the binary, and a green run of it is only evidence that the two agree —
not that the comparison is capable of noticing when they don't. Golden files have the same weakness (docs/TESTING.md says so:
"Golden files only catch *changes*"), so this applies the same idea the 2026-11-11 suite review applied to the Rust tests:
break one documented fact at a time in a throwaway copy of the tree, run the checks against the copy, and require the specific
check that should notice to notice — and, on the unmutated copy, that nothing fails.

Each mutation is a drift that either actually happened on main or is the twin of one that did: the `demo` table showing the
pre-fix `ambiguous` value, the sample dropping the `# format=2` header or the self-check line, an OUTPUTS.md row for a file
the tool writes deleted, a `summary_score.tsv` key renamed away, a §8 troubleshooting row deleted, an option name mangled,
and a §4 row for a flag the binary has never had.

Run: python3 -I scripts/mutate_doc_samples.py [path/to/intronfate]      Standard library only; needs a built binary.
Exit status: 0 if every mutation was caught by the intended check and the clean copy passed, 1 otherwise.

The copy is a temporary directory holding README.md, docs/ and the script under test; the repository is never modified.
"""
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
CHECKS = "test_doc_samples.py"

# (label, file to mutate or None, fragment to replace, replacement, test that must fail, replace every occurrence?)
MUTATIONS = [
    ("baseline: nothing mutated", None, None, None, None, False),
    ("README: the pre-fix `ambiguous yes`", "README.md",
     "4:0   4       0      4          0          1      no           no",
     "4:0   4       0      4          0          1      no           yes", "test_demo_block_matches_the_binary", False),
    ("README: `# format=2` header deleted", "README.md", "# format=2 alignment=demo", "# deleted",
     "test_the_sample_shows_the_run_identity", False),
    ("README: self-check line deleted", "README.md", "expect 3 sites, 3 Dollo gains", "dropped",
     "test_the_sample_shows_the_run_identity", False),
    ("OUTPUTS: row for a written file deleted", "docs/OUTPUTS.md", "| `characters.nex` |", "| `deleted` |",
     "test_every_written_file_has_a_row", False),
    ("OUTPUTS: a summary_score key deleted", "docs/OUTPUTS.md", "tips_without_data", "removed_key",
     "test_every_summary_score_key_is_named", False),
    ("USAGE §8: a reachable diagnostic row deleted", "docs/USAGE.md", "| `row_sequence_mismatch` |", "| `gone` |",
     "test_every_emitted_code_has_a_troubleshooting_row", False),
    ("USAGE: an option name mangled", "docs/USAGE.md", "--keep-noncoding", "keep-noncoding",
     "test_every_option_is_documented", True),
    ("USAGE §4: a flag the binary does not have", "docs/USAGE.md", "| `--strict` |", "| `--does-not-exist` |",
     "test_the_options_table_describes_real_options", False),
]


def find_binary():
    """`argv[1]`, else $INTRONFATE, else this workspace's build. An explicit path must exist."""
    arg = sys.argv[1] if len(sys.argv) > 1 else ""
    explicit = arg if (os.sep in arg or arg.endswith("intronfate")) else os.environ.get("INTRONFATE")
    if explicit:
        if not os.path.exists(explicit):
            sys.exit(f"error: no intronfate binary at {explicit!r} (build it, or fix the path/INTRONFATE)")
        return os.path.abspath(explicit)
    for profile in ("release", "debug"):
        p = os.path.join(REPO, "target", profile, "intronfate")
        if os.path.exists(p):
            return p
    sys.exit("error: no intronfate binary in target/{release,debug}; cargo build --release first or pass its path")


def repo_files():
    """The files the checks read: README.md, docs/*.md and the check script itself."""
    return ["README.md", os.path.join("scripts", CHECKS)] + [
        os.path.join("docs", n) for n in sorted(os.listdir(os.path.join(REPO, "docs"))) if n.endswith(".md")
    ]


def main():
    binary = find_binary()
    work = tempfile.mkdtemp(prefix="doc-mutation-")
    for rel in repo_files():
        dest = os.path.join(work, rel)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        shutil.copy(os.path.join(REPO, rel), dest)
    script = os.path.join(work, "scripts", CHECKS)

    print("mutation".ljust(46), "result")
    print("-" * 100)
    problems = 0
    try:
        for label, rel, old, new, expect, every in MUTATIONS:
            for f in repo_files():                       # reset the copy before each mutation
                shutil.copy(os.path.join(REPO, f), os.path.join(work, f))
            if rel:
                path = os.path.join(work, rel)
                text = open(path, encoding="utf-8").read()
                if old not in text:
                    print(label.ljust(46), f"HARNESS BUG: {old!r} is not in {rel} (the document changed; update this list)")
                    problems += 1
                    continue
                open(path, "w", encoding="utf-8").write(text.replace(old, new) if every else text.replace(old, new, 1))

            r = subprocess.run([sys.executable, "-I", script], capture_output=True, text=True,
                               cwd=work, env=dict(os.environ, INTRONFATE=binary))
            failed = set(re.findall(r"(?:FAIL|ERROR): (\w+)", r.stdout + r.stderr))
            hit = next((t for t in sorted(failed) if t.startswith(expect)), None) if expect else None
            if expect is None:                            # the clean copy has to pass, not merely fail elsewhere
                ok = r.returncode == 0
                verdict = "passed, as it must" if ok else "FAILED on an unmutated tree: " + ", ".join(sorted(failed))
            elif hit:
                ok, verdict = True, f"caught by {hit}"
            else:
                ok = False
                verdict = ("MISSED (nothing failed)" if not failed
                           else "MISSED (a different check failed: " + ", ".join(sorted(failed)) + ")")
            problems += 0 if ok else 1
            print(label.ljust(46), verdict)
    finally:
        shutil.rmtree(work, ignore_errors=True)

    if problems:
        print(f"\n{problems} problem(s): a documented fact can be wrong without the right check saying so.")
        return 1
    print(f"\nall {len(MUTATIONS) - 1} mutations caught by the intended check, and the clean copy passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
