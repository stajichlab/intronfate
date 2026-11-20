#!/usr/bin/env python3
"""Fixture test for J1_junction_support/bam_junctions.py, including the coordinate convention end to end.

The last test builds `intronfate` output for the quickstart example and checks that a spliced read crossing the annotated intron
(113..200 on chrA, GFF coordinates) is exact support: if the script or the tool is off by one the intron is "unsupported".
Run: python3 -I bench/test_junction_scripts.py [path/to/intronfate]     Standard library only; the last test is skipped without a binary.
"""
import os, subprocess, sys, tempfile, unittest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SCRIPT = os.path.join(HERE, "J1_junction_support", "bam_junctions.py")
BIN = (sys.argv[1] if len(sys.argv) > 1 and os.path.exists(sys.argv[1]) else os.path.join(REPO, "target", "debug", "intronfate"))


def sam(*recs):
    return "@HD\tVN:1.6\n" + "".join("\t".join(map(str, r)) + "\n" for r in recs)


def rec(flag, chrom, pos, cigar, *tags):
    return ("r", flag, chrom, pos, 60, cigar, "*", 0, 0, "*", "*", *tags)


def run(text, *args):
    return subprocess.run([sys.executable, "-I", SCRIPT, *args], input=text, capture_output=True, text=True, check=True).stdout


class Junctions(unittest.TestCase):
    def test_first_and_last_intron_base_and_bed_start(self):
        # 50M (63..112) then 88N covers 113..200 then 50M: first intron base 113, last 200; BED start is 0-based
        out = run(sam(rec(0, "chrA", 63, "50M88N50M", "XS:A:+")))
        self.assertEqual(out, "chrA\t112\t200\tj\t1\t+\n")

    def test_counts_strand_and_skipped_records(self):
        text = sam(
            rec(0, "chrA", 63, "50M88N50M", "XS:A:+"),
            rec(16, "chrA", 63, "50M88N50M", "XS:A:+"),
            rec(0x100, "chrA", 63, "50M88N50M", "XS:A:+"),      # secondary: skipped
            rec(0x800, "chrA", 63, "50M88N50M", "XS:A:+"),      # supplementary: skipped
            rec(4, "*", 0, "*"),                                 # unmapped: skipped
            rec(0, "chrA", 10, "30M10N30M", "XS:A:-"),          # 10N is shorter than --min-intron 20: not an intron
            rec(0, "chrB", 5, "20M300N20M"),                    # no XS tag: strand unknown
        )
        out = run(text).splitlines()
        self.assertEqual(out, ["chrA\t112\t200\tj\t2\t+", "chrB\t24\t324\tj\t1\t."])

    def test_deletions_advance_the_reference_but_insertions_and_soft_clips_do_not(self):
        # 10S 20M 5D 5I 15M 100N 10M from pos 1000: ref after 20M=1020, after 5D=1025, 5I and 10S add none, after 15M=1040
        out = run(sam(rec(0, "c", 1000, "10S20M5D5I15M100N10M", "XS:A:+")))
        self.assertEqual(out, "c\t1039\t1139\tj\t1\t+\n")   # first intron base 1040 -> BED 1039, last 1139

    def test_min_reads_filter(self):
        text = sam(rec(0, "c", 1, "10M100N10M"), rec(0, "c", 1, "10M100N10M"), rec(0, "c", 1, "10M200N10M"))
        self.assertEqual(len(run(text).splitlines()), 2)
        self.assertEqual(run(text, "--min-reads", "2").splitlines(), ["c\t10\t110\tj\t2\t."])


@unittest.skipUnless(os.path.exists(BIN), "no intronfate binary (cargo build first, or pass its path)")
class EndToEnd(unittest.TestCase):
    def test_a_spliced_read_over_the_annotated_intron_is_exact_support(self):
        q = os.path.join(REPO, "examples", "01_quickstart")
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "o")
            subprocess.run([BIN, "analyze", "--gff", f"{q}/family.gff3", "-p", f"{q}/proteins.faa", "--msa", f"{q}/family.aln.fasta",
                            "-t", f"{q}/family.nwk", "-m", "counts", "-o", out], check=True, capture_output=True)
            # taxon A's first intron is chrA 113..200; 12 reads cross it exactly, and 4 reads one base short (not support)
            bed = run(sam(*[rec(0, "chrA", 63, "50M88N50M", "XS:A:+")] * 12, *[rec(0, "chrA", 63, "50M87N50M", "XS:A:+")] * 4))
            bedfile = os.path.join(d, "A.bed")
            with open(bedfile, "w") as fh:
                fh.write(bed)
            r = subprocess.run([BIN, "junctions", "--introns", f"{out}/introns.tsv", "--sites", f"{out}/family.aln/sites.tsv",
                                "--matrix", f"{out}/family.aln/matrix.tsv", "--junctions", f"A={bedfile}", "--junction-format", "bed", "-o", "-"], capture_output=True, text=True, check=True)
            rows = {tuple(l.split("\t")[:2]): l.split("\t") for l in r.stdout.splitlines() if not l.startswith("#")}
            first = rows[("4:0", "A")]
            self.assertEqual(first[3], "12", r.stdout)            # exactly the 12 reads at 113..200
            self.assertAlmostEqual(float(first[2]), 1 - 0.2 * 3 / 15, places=5)
            self.assertEqual(rows[("7:0", "A")][3], "0", r.stdout)


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
