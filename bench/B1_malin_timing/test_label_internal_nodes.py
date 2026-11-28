#!/usr/bin/env python3
"""Fixture test for label_internal_nodes.py."""
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from label_internal_nodes import label  # noqa: E402


class Label(unittest.TestCase):
    def test_unnamed_internals_are_named_in_postorder(self):
        self.assertEqual(label("((A,B),(C,D));"), "((A,B)n1,(C,D)n2)n3;")

    def test_branch_lengths_and_existing_labels_are_kept(self):
        self.assertEqual(label("((A:1,B:2):0.5,C:3)root;"), "((A:1,B:2)n1:0.5,C:3)root;")


if __name__ == "__main__":
    unittest.main()
