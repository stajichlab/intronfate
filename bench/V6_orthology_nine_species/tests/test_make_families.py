"""Fixture test for make_families.py: 3 species, MIN=2. Run: python3 -m pytest tests/ (or python3 tests/test_make_families.py)."""
import os, sys, tempfile
here = os.path.dirname(os.path.abspath(__file__)); fx = os.path.join(here, "fixture")
sys.path.insert(0, os.path.dirname(here))
import make_families

def fams(out, tier):
    return sorted(f[:-4] for f in os.listdir(os.path.join(out, tier, "families")))

def seqs(out, tier, fam):
    return [l[1:].strip() for l in open(os.path.join(out, tier, "families", fam + ".faa")) if l.startswith(">")]

def test_tiers():
    with tempfile.TemporaryDirectory() as out:
        make_families.main(["--of", f"{fx}/of", "--data", f"{fx}/data", "--min", "2", "--out", out, "--refined", f"{fx}/syn_groups.tsv"])
        assert fams(out, "A_strict") == ["OG1"]                 # a1 b1 c1
        assert fams(out, "B_partial") == ["OG2"]                # a2 b2, species C absent; OG3 (one species) is too small
        assert fams(out, "C_hog") == ["N1.HOG0"]                # from multi-copy OG4; OG5's HOG still has 2 copies in A; OG1's HOG is not multi-copy
        assert seqs(out, "C_hog", "N1.HOG0") == ["a4", "b4"]    # c4 has an internal stop and is dropped
        assert fams(out, "D_synteny") == ["SOG1.OG5"]
        tm = open(os.path.join(out, "A_strict", "taxon_map.tsv")).read().split("\n")
        assert tm[:3] == ["a1\tA", "b1\tB", "c1\tC"]
        assert not open(os.path.join(out, "A_strict", "families", "OG1.faa")).read().count("*")

if __name__ == "__main__":
    test_tiers(); print("ok")
