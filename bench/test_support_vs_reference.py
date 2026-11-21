#!/usr/bin/env python3
"""Tests for J2_support_vs_reference/support_vs_reference.py (coordinates, strand rule, AUC, set bookkeeping)."""
import os, subprocess, sys, unittest

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "J2_support_vs_reference", "support_vs_reference.py")


class SelfTest(unittest.TestCase):
    def test_selftest_passes(self):
        r = subprocess.run([sys.executable, "-I", SCRIPT, "--selftest"], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("selftest ok", r.stdout)

    def test_matrix_cells_selftest_passes(self):
        script = os.path.join(os.path.dirname(SCRIPT), "..", "J2b_matrix_cells", "matrix_cells.py")
        r = subprocess.run([sys.executable, "-I", script, "--selftest"], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("selftest ok", r.stdout)

    def test_cell_power_selftest_passes(self):
        script = os.path.join(os.path.dirname(SCRIPT), "..", "J2b_matrix_cells", "cell_power.py")
        r = subprocess.run([sys.executable, "-I", script, "--selftest"], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("selftest ok", r.stdout)

    def test_creep_selftest_passes(self):
        script = os.path.join(os.path.dirname(SCRIPT), "..", "J3b_creeping_parameters", "analyse.py")
        r = subprocess.run([sys.executable, "-I", script, "--selftest"], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("selftest ok", r.stdout)

    def test_simulation_selftests_pass(self):
        for name in ("simulate.py", "evaluate.py"):
            script = os.path.join(os.path.dirname(SCRIPT), "..", "S1_soft_evidence_simulation", name)
            r = subprocess.run([sys.executable, "-I", script, "--selftest"], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("selftest ok", r.stdout)


if __name__ == "__main__":
    unittest.main()
