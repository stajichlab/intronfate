"""Fixture test for make_families_pruned.py (3 genomes + 2 strains of one lineage). Run: python3 tests/test_make_families_pruned.py"""
import os, sys, tempfile
here = os.path.dirname(os.path.abspath(__file__)); fx = os.path.join(here, "fixture")
sys.path.insert(0, os.path.dirname(here))
import make_families_pruned as M

def fams(out, t): return sorted(f[:-4] for f in os.listdir(os.path.join(out, t, "families")))
def rows(out, t, fam): return [l[1:].strip() for l in open(os.path.join(out, t, "families", fam + ".faa")) if l.startswith(">")]

def test_tiers():
    with tempfile.TemporaryDirectory() as out:
        M.main(["--of", f"{fx}/of", "--data", f"{fx}/data", "--out", out, "--group-prefix", "G", "--tiers", "A=4,B=3"])
        assert fams(out, "A") == ["OG1", "OG2", "OG6"]
        assert fams(out, "B") == ["OG3"]                       # Y has two genes: missing; X, Z and the group remain
        assert rows(out, "A", "OG2") == ["g2b", "x2", "y2", "z2"] or sorted(rows(out, "A", "OG2")) == ["g2b", "x2", "y2", "z2"]
        assert sorted(rows(out, "A", "OG6")) == ["g2f", "x6", "y6", "z6"]      # G1's protein has an internal stop: dropped
        assert sorted(rows(out, "B", "OG3")) == ["g1c", "x3", "z3"]            # the multi-copy genome contributes nothing
        tm = open(os.path.join(out, "B", "taxon_map.tsv")).read()
        assert "y3b" not in tm and "g1c\tG1" in tm
        assert not any(f in fams(out, "A") + fams(out, "B") for f in ["OG4", "OG5"])   # too few lineages / group absent

if __name__ == "__main__":
    test_tiers(); print("ok")
