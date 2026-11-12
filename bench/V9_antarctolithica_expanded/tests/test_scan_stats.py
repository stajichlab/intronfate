"""Fixture tests for scan_stats.py. Run: python3 tests/test_scan_stats.py"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import scan_stats as S

# 4 families x 4 trees; trees 0 and 1 are identical and best, tree 2 is worse everywhere, tree 3 is best in family 0 only
X = np.array([[1, 1, 3, 0], [2, 2, 5, 9], [1, 1, 2, 9], [0, 0, 4, 9]])
W = np.eye(4, dtype=int) * 4 + 0                                     # four replicates, each resampling one family 4 times
f = S.win_fractions(W, X)
assert f[0] == f[1], f                                                # twins share one fraction
assert f[2] == 0.0 and abs(f[0] - 0.75) < 1e-9 and abs(f[3] - 0.25) < 1e-9, f   # twins win 3 of 4 replicates (not 3/8 each); tree 3 wins the first
# without the tie handling the twins would split their wins, so the best tree could never exceed 0.5
assert S.twins(X, ["a", "b", "c", "d"]) == [["b"], ["a"], [], []]
# ties between DIFFERENT groups split the replicate equally
Y = np.array([[1, 1, 2], [1, 3, 1]]); V = np.array([[1, 0]])        # weights pick family 0: totals 1, 1, 2 -> trees 0 and 1 (different groups) tie
assert np.allclose(S.win_fractions(V, Y), [0.5, 0.5, 0.0])
Z = np.array([[1, 2], [2, 1]]); U = np.array([[1, 1]])               # totals 3 vs 3: a genuine tie between two groups
assert np.allclose(S.win_fractions(U, Z), [0.5, 0.5])
print("ok")
