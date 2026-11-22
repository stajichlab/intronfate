#!/usr/bin/env python3
"""Checks that the prose pages still describe the binary, so a sample cannot drift silently.

The examples have byte-for-byte golden files and the CLI has end-to-end tests; nothing covered the README and docs pages
themselves, and four drifts had accumulated on main by 2026-11-21 (issues #16, #17): a `demo` table still printed the
`ambiguous` value from before that flag was fixed, "Not yet" listed three shipped features, codes that a real
`diagnostics.tsv` emits were absent from the section-8 table, and `summary_score.tsv` had keys OUTPUTS.md did not name.

Every check compares a document against the tool as built now, never against another document, so the two cannot agree by
being wrong together. Two honest limits: this cannot tell whether a sample was right when it was written (only that it is
right now), and the diagnostic-code check sees only the codes its fixture reaches — `no_protein_id` in the transcript->protein
fallback path, for example, is documented but is not one of them, which is why the fixture also asserts its own breadth.

Run: python3 -I scripts/test_doc_samples.py [path/to/intronfate]     Standard library only; skipped without a built binary.
     INTRONFATE=/path/to/intronfate python3 -I scripts/test_doc_samples.py
"""
import os
import re
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
DOC_FILES = ["README.md"] + [
    os.path.join("docs", n) for n in
    ("HOWTO.md", "TUTORIAL.md", "USAGE.md", "OUTPUTS.md", "SPEC.md", "TESTING.md", "HANDOFF.md")
]
COMMANDS = ["project", "matrix", "score", "analyze", "convert", "junctions", "demo"]
OUTPUT_SUFFIXES = (".tsv", ".nex", ".nwk")


def read(*rel):
    """Read a repository file as text."""
    with open(os.path.join(REPO, *rel), encoding="utf-8") as h:
        return h.read()


DOCS = "\n".join(read(rel) for rel in DOC_FILES if os.path.exists(os.path.join(REPO, rel)))


def find_binary():
    """`argv[1]`, else $INTRONFATE, else the release/debug build in this workspace.

    An explicit path that does not exist is an error, never a skip: a typo in $INTRONFATE must not turn the whole file into
    eight tests that pass by running nothing. Only the auto-detected default may be absent.
    """
    arg = sys.argv[1] if len(sys.argv) > 1 else ""
    explicit = arg if (os.sep in arg or arg.endswith("intronfate")) else os.environ.get("INTRONFATE")
    if explicit:
        if not os.path.exists(explicit):
            sys.exit(f"error: no intronfate binary at {explicit!r} (build it, or fix the path/INTRONFATE)")
        # Absolute: several checks run the tool in a temporary working directory, where a relative path would not resolve.
        return os.path.abspath(explicit)
    for profile in ("release", "debug"):
        p = os.path.join(REPO, "target", profile, "intronfate")
        if os.path.exists(p):
            return p
    return os.path.join(REPO, "target", "release", "intronfate")


BIN = find_binary()
if len(sys.argv) > 1 and (os.sep in sys.argv[1] or sys.argv[1].endswith("intronfate")):
    sys.argv.pop(1)                       # the binary path is ours to read; leave unittest its own arguments


def run(args, cwd=None):
    """One intronfate invocation; the caller decides whether the exit status matters."""
    return subprocess.run([BIN] + args, capture_output=True, text=True, cwd=cwd, check=False)


def must_run(args, cwd=None):
    r = run(args, cwd=cwd)
    if r.returncode != 0:
        raise RuntimeError(f"intronfate {' '.join(args[:2])} failed ({r.returncode}): {r.stderr.strip()[:400]}")
    return r


def help_text(subcommand=None):
    return run([subcommand, "--help"] if subcommand else ["--help"]).stdout


def long_flags_with_aliases(text):
    """Every long flag in a `--help` dump, with its short alias where clap prints one (`-m, --methods`)."""
    pairs = {long: short for short, long in re.findall(r"-(\w), --([a-z][a-z0-9-]+)", text)}
    for long in re.findall(r"--([a-z][a-z0-9-]+)", text):
        pairs.setdefault(long, None)
    return pairs


def fenced(text, fence, containing=None):
    """Bodies of ```fence blocks, optionally only those containing a marker string."""
    return [b for b in re.findall(rf"```{fence}\n(.*?)```", text, re.S)
            if containing is None or containing in b]


def table_rows(block):
    """A sample as comparable rows: tabs when present (the tool's own output), otherwise runs of spaces."""
    rows = []
    for line in block.splitlines():
        if not line.strip():
            continue
        line = re.sub(r"written to .*", "written to", line)      # "-o DIR" echoes the directory it was given
        cells = line.split("\t") if "\t" in line else line.split()
        while len(cells) > 1 and cells[-1] == "":
            cells.pop()                                          # the space-padded rendering cannot show a trailing empty field
        rows.append([c.strip() for c in cells])
    return rows


def documented_output_names():
    """File names OUTPUTS.md has a table row for — a row, not any mention, so prose cannot stand in for documentation."""
    outputs = read("docs", "OUTPUTS.md")
    return {n for n in re.findall(r"^\| `([^`]+)`", outputs, re.M) if n.endswith(OUTPUT_SUFFIXES)}


@unittest.skipUnless(os.path.exists(BIN), f"no intronfate binary at {BIN} (cargo build --release, or pass its path)")
class TestReadmeDemoSample(unittest.TestCase):
    """The first thing a newcomer runs and diffs against."""

    def block(self):
        readme = read("README.md")
        blocks = fenced(readme, "text", containing="n_observable")
        self.assertEqual(len(blocks), 1, "expected exactly one demo-output block in README.md")
        return blocks[0]

    def test_demo_block_matches_the_binary(self):
        with tempfile.TemporaryDirectory() as tmp:                # demo writes its worked example wherever it is run
            actual = must_run(["demo", "-o", "intron_demo"], cwd=tmp).stdout
        documented, real = table_rows(self.block()), table_rows(actual)
        self.assertEqual(len(documented), len(real), f"README block has {len(documented)} rows, `demo` prints {len(real)}")
        for i, (a, b) in enumerate(zip(documented, real)):
            self.assertEqual(a, b, f"README line {i + 1}\n  documented: {a}\n  actual:     {b}")

    def test_the_sample_shows_the_run_identity_and_the_self_check(self):
        """Both lines the README dropped the first time: the matrix header and the line where demo grades itself."""
        block = self.block()
        self.assertRegex(block, r"(?m)^# format=2 alignment=demo", "the sample drops the `# format=2 …` header `demo` prints first")
        self.assertRegex(block, r"(?m)^expect \d+ sites, .*: got ", "the sample drops the demo self-check line")


@unittest.skipUnless(os.path.exists(BIN), f"no intronfate binary at {BIN} (cargo build --release, or pass its path)")
class TestDocumentedFlags(unittest.TestCase):
    """Every option the binary offers is described somewhere, and the reference table offers nothing the binary lacks."""

    def flags(self):
        flags = long_flags_with_aliases(help_text())
        for cmd in COMMANDS:
            flags.update(long_flags_with_aliases(help_text(cmd)))
        return flags

    def test_every_option_is_documented(self):
        # A long form counts as documented when a page names it, or names the short flag clap shows beside it.
        missing = sorted(f for f, short in self.flags().items()
                         if f"--{f}" not in DOCS and not (short and f"-{short}" in DOCS))
        self.assertEqual(missing, [], "options in --help that no page mentions: " + ", ".join("--" + f for f in missing))

    def test_the_options_table_describes_real_options(self):
        """USAGE §4's table is the reference list; a row for a flag the binary dropped is worse than no row."""
        usage = read("docs", "USAGE.md")
        self.assertIn("### Options that change the result", usage, "USAGE.md lost the section the table test reads")
        table = usage.split("### Options that change the result")[1].split("\n## ")[0]
        known = self.flags()
        claimed = {row.lstrip("-") for row in re.findall(r"^\| `(-{1,2}[a-z][a-z0-9-]+)[^`]*`", table, re.M)}
        short_forms = {"m", "t", "p", "o", "v"}                   # documented by their long form
        stale = sorted(c for c in claimed if c not in known and c not in short_forms)
        self.assertTrue(claimed, "the options table yielded no rows to check")
        self.assertEqual(stale, [], "USAGE §4 lists options the binary does not offer: " + ", ".join("--" + c for c in stale))


@unittest.skipUnless(os.path.exists(BIN), f"no intronfate binary at {BIN} (cargo build --release, or pass its path)")
class TestOutputsMatchTheTool(unittest.TestCase):
    """OUTPUTS.md promises 'every output file and column' (README twice, SPEC); the tool has to agree with it."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tmp = tempfile.TemporaryDirectory()
        must_run(["demo", "-o", os.path.join(tmp.name, "intron_demo")])
        d = os.path.join(tmp.name, "intron_demo")
        gff, faa, aln, nwk = (os.path.join(d, n) for n in ("demo.gff3", "demo.faa", "demo.aln.fasta", "demo.nwk"))
        # All five methods + combining + reconstructions, so every table the tool can write exists.
        must_run(["analyze", "--gff", gff, "-p", faa, "--msa", aln, "-t", nwk, "--outgroup", "outgroup",
                  "-m", "dollo,mp,counts,dollo_weighted,likelihood", "--combine", "--reconstructions",
                  "-o", os.path.join(tmp.name, "all")])
        must_run(["project", "--gff", gff, "-p", faa, "-o", os.path.join(tmp.name, "proj")])
        # The two likelihood forms that add their own summary keys; free levels and gamma classes refuse each other, so two runs.
        # --ml-max-iter is a test-only speed knob (the rate-class fit creeps for 45 s on 3 sites); which keys a form writes
        # does not depend on it converging.
        mx = os.path.join(tmp.name, "all", "combined", "matrix.tsv")
        for i, extra in enumerate((["--ml-rate-categories", "2"], ["--ml-loss-levels", "2", "--ml-origins"])):
            must_run(["score", "--matrix", mx, "-t", nwk, "--outgroup", "outgroup", "-m", "likelihood"] + extra
                     + ["--ml-max-iter", "300", "-o", os.path.join(tmp.name, f"ml{i}")])

        cls.written, cls.summary_keys = set(), set()
        for run_name in ("all", "proj", "ml0", "ml1"):           # the demo folder holds inputs, not our outputs
            for dirpath, _, names in os.walk(os.path.join(tmp.name, run_name)):
                for n in names:
                    if not n.endswith(OUTPUT_SUFFIXES):
                        continue
                    cls.written.add(n)
                    if n == "summary_score.tsv":
                        with open(os.path.join(dirpath, n), encoding="utf-8") as h:
                            cls.summary_keys |= {l.split("\t")[0] for l in h.read().splitlines()[1:] if "\t" in l}

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_every_written_file_has_a_row(self):
        self.assertTrue(len(self.written) >= 20, f"the runs wrote only {len(self.written)} files; the fixture has shrunk")
        undocumented = sorted(f for f in self.written if f not in documented_output_names())
        self.assertEqual(undocumented, [], "files the tool writes that OUTPUTS.md has no row for: " + ", ".join(undocumented))

    def test_every_summary_score_key_is_named(self):
        """The summary is the file every page points at first; an unnamed key is an unexplained number."""
        outputs = read("docs", "OUTPUTS.md")
        row = [l for l in outputs.splitlines() if l.startswith("| `summary_score.tsv`")]
        self.assertEqual(len(row), 1, "OUTPUTS.md should have exactly one summary_score.tsv row")
        self.assertTrue(len(self.summary_keys) >= 25, f"only {len(self.summary_keys)} keys collected; the fixture has shrunk")
        missing = [k for k in sorted(self.summary_keys) if k and k not in row[0]]
        self.assertEqual(missing, [], "summary_score.tsv writes keys OUTPUTS.md does not name: " + ", ".join(missing))


@unittest.skipUnless(os.path.exists(BIN), f"no intronfate binary at {BIN} (cargo build --release, or pass its path)")
class TestDiagnosticCodesAreTroubleshot(unittest.TestCase):
    """USAGE §8 is where a reader lands with an unknown code in diagnostics.tsv; a missing row is a dead end.

    The fixture is one messy run chosen to reach the codes people actually hit: no protein behind a row, a trimmed row, a
    protein whose length disagrees with its CDS, a transcript with no CDS, an empty file, a protein no transcript names.
    """

    GFF = {
        # no protein_id= and the transcript id is not a sequence either
        "a_noprotein.gff3": "##gff-version 3\nchr1\tsrc\tCDS\t1\t30\t.\t+\t0\tID=c1;Parent=t1\n",
        # a transcript with no CDS
        "b_nocds.gff3": "##gff-version 3\nchr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=t2;gene=g2\n",
        # names a protein that is not in the FASTA
        "c_ghost.gff3": "##gff-version 3\nchr1\tsrc\tCDS\t1\t30\t.\t+\t0\tID=c3;Parent=t3;protein_id=GHOST\n",
        # the same transcript id again, in a second file: two records, reported rather than merged
        "d_collision.gff3": "##gff-version 3\nchr1\tsrc\tCDS\t1\t30\t.\t+\t0\tID=c3b;Parent=t3;protein_id=GHOST\n",
        # nothing parseable in it
        "e_empty.gff3": "##gff-version 3\n",
        # 30 coding nt against a 30-residue protein
        "f_length.gff3": "##gff-version 3\nchr1\tsrc\tCDS\t1\t30\t.\t+\t0\tID=c5;Parent=t5;protein_id=P5\n",
        # P6 is annotated, so its W-filled alignment row is checked and refused (a mismatch is only a
        # diagnostic for a row that had annotation to lose)
        "g_p6.gff3": "##gff-version 3\nchr1\tsrc\tCDS\t1\t30\t.\t+\t0\tID=c6;Parent=t6;protein_id=P6\n",
        # one file, six independent malformed-model defects: no transcript key, an unknown strand, overlapping CDS
        # and exon rows, a CDS gap that is not an exon gap, and a phase that disagrees with the derived frame
        "h_malformed.gff3": (
            "##gff-version 3\n"
            "chr1\tsrc\tCDS\t1\t30\t.\t+\t0\tID=nok;gene=g1\n"   # no Parent, no transcript_id: nothing to group it by
            "chr1\tsrc\tCDS\t1\t30\t.\t.\t0\tID=su1;Parent=su_t;protein_id=PSU\n"
            "chr1\tsrc\tCDS\t1\t40\t.\t+\t0\tID=ov1;Parent=ov_t;protein_id=POV\n"
            "chr1\tsrc\tCDS\t20\t60\t.\t+\t0\tID=ov2;Parent=ov_t;protein_id=POV\n"
            "chr1\tsrc\texon\t1\t30\t.\t+\t.\tID=ex1;Parent=ex_t\n"
            "chr1\tsrc\texon\t20\t60\t.\t+\t.\tID=ex2;Parent=ex_t\n"
            "chr1\tsrc\tCDS\t1\t30\t.\t+\t0\tID=cg1;Parent=cg_t;protein_id=PCG\n"
            "chr1\tsrc\tCDS\t60\t90\t.\t+\t0\tID=cg2;Parent=cg_t;protein_id=PCG\n"
            "chr1\tsrc\texon\t1\t90\t.\t+\t.\tID=cge;Parent=cg_t\n"
            "chr1\tsrc\tCDS\t1\t30\t.\t+\t0\tID=fr1;Parent=fr_t;protein_id=PFR\n"
            "chr1\tsrc\tCDS\t40\t60\t.\t+\t2\tID=fr2;Parent=fr_t;protein_id=PFR\n"
        ),
    }
    PROTEINS = (">P5\n" + "A" * 30 + "\n>P6\nMKVLGAWEHT\n>Lonely\nMM\n"
                + "".join(f">{n}\n{'A' * 20}\n" for n in ("PSU", "POV", "PCG", "PFR")))
    # W, not X: row_matches_protein treats X as a wildcard on either side, so a row of X's would match anything.
    ALIGNMENT = ">P6\nWWWWWWWWWW\n>NotInAnySet\nMKVLGAWEHT\n"   # a row that is not its protein, and a row with no protein behind it

    @classmethod
    def setUpClass(cls):
        cls.tmp = tmp = tempfile.TemporaryDirectory()
        w = lambda n, s: open(os.path.join(tmp.name, n), "w").write(s)
        for n, body in cls.GFF.items():
            w(n, body)
        w("proteins.faa", cls.PROTEINS); w("fam.aln", cls.ALIGNMENT)
        # `matrix`, not `analyze`: every code here is raised before scoring, and with no scored taxa left
        # there would be nothing to score once each row is unannotated.
        args = ["matrix"] + sum((["--gff", os.path.join(tmp.name, n)] for n in sorted(cls.GFF)), [])
        must_run(args + ["-p", os.path.join(tmp.name, "proteins.faa"), "--msa", os.path.join(tmp.name, "fam.aln"),
                         "-o", os.path.join(tmp.name, "out")])
        cls.codes = set()
        for name, column in (("diagnostics.tsv", "kind"), ("unmatched.tsv", "reason")):
            p = os.path.join(tmp.name, "out", name)
            if not os.path.exists(p):
                continue
            with open(p, encoding="utf-8") as h:
                header = h.readline().rstrip("\n").split("\t")
                if column not in header:
                    raise RuntimeError(f"{name} no longer has a {column} column: {header}")
                i = header.index(column)
                cls.codes |= {l.split("\t")[i].strip() for l in h if l.strip()}

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_the_fixture_still_reaches_the_common_codes(self):
        """Without this, the check below could pass by the fixture quietly emitting nothing.

        17 of the 28 codes §8 documents. The rest need a second annotation form (GTF, `--derive-cds-from-exons`),
        an alignment with gaps in it, or bytes no text editor writes (`invalid_utf8`); `no_protein_id` is not among
        them, so the coverage claim here is deliberately a floor, not "all of §8".
        """
        expected = {"row_sequence_mismatch", "protein_not_in_sequence_set", "protein_not_in_fasta",
                    "length_mismatch", "no_cds", "empty_file", "strand_unknown", "overlapping_cds",
                    "overlapping_exons", "cds_gap_not_exon_gap", "frame_mismatch", "no_transcript_key",
                    "transcript_id_collision", "protein_id_fallback_transcript", "protein_without_annotation",
                    "proteins_not_in_alignment"}
        absent = sorted(expected - self.codes)
        self.assertEqual(absent, [], f"the fixture no longer triggers {absent}; widen it before trusting the next check")

    def test_every_emitted_code_has_a_troubleshooting_row(self):
        usage = read("docs", "USAGE.md")
        self.assertIn("## 8. Troubleshooting", usage, "USAGE.md lost its troubleshooting section")
        section = usage.split("## 8. Troubleshooting")[1].split("\n## 9.")[0]
        missing = sorted(c for c in self.codes if c and f"`{c}`" not in section)
        self.assertEqual(missing, [], "codes a real run emits that USAGE §8 does not explain: " + ", ".join(missing))


if __name__ == "__main__":
    if not os.path.exists(BIN):
        print(f"note: no binary at {BIN}; build it (cargo build --release) or pass its path to run these checks",
              file=sys.stderr)
    unittest.main(verbosity=2)
