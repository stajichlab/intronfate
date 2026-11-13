"""Fixture test for site_filter.py. Run: python3 tests/test_site_filter.py"""
import os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from site_filter import recode, read_rates, drop_fast
s = {"a": "AGDEHKILFWC-XB", "b": "PSNQRVMYTc--?Z"}
r = recode(s)
assert r["a"] == "00112233445-??", r["a"]
assert r["b"] == "0011233405--??", r["b"]          # lowercase c recoded like C; ? kept as ?
d = tempfile.mkdtemp(); p = f"{d}/x.rate"
open(p, "w").write("# comment\n# more\nSite\tRate\tCat\tC_Rate\n1\t0.1\t1\t0.2\n2\t3.0\t4\t3.0\n3\t0.5\t2\t0.5\n4\t2.0\t3\t2.0\n5\t0.2\t1\t0.2\n")
rates = read_rates(p); assert rates == {1: 0.1, 2: 3.0, 3: 0.5, 4: 2.0, 5: 0.2}, rates
out, dropped = drop_fast({"a": "ABCDE", "b": "FGHIK"}, rates, 0.4)
assert dropped == [2, 4] and out == {"a": "ACE", "b": "FHK"}, (dropped, out)       # the two fastest sites (3.0, 2.0) removed
out, dropped = drop_fast({"a": "ABCDE"}, rates, 0.0); assert out["a"] == "ABCDE" and dropped == []
try: drop_fast({"a": "ABC"}, rates, 0.5); assert False
except AssertionError as e: assert "rate file" in str(e)
print("ok")
