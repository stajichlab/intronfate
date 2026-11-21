#!/usr/bin/env python3
"""Tests for J2_support_vs_reference/support_vs_reference.py (coordinates, strand rule, AUC, set bookkeeping)."""
import os, subprocess, sys, unittest

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "J2_support_vs_reference", "support_vs_reference.py")


class SelfTest(unittest.TestCase):
    def test_selftest_passes(self):
        r = subprocess.run([sys.executable, "-I", SCRIPT, "--selftest"], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("selftest ok", r.stdout)


if __name__ == "__main__":
    unittest.main()
