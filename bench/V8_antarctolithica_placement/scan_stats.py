"""Statistics for the placement scans (pure numpy, fixture-tested in tests/test_scan_stats.py)."""
import numpy as np

def groups(X):
    """(unique column vectors, group id per tree): trees whose per-family score vectors are identical share a group (a tie group)."""
    uniq, inv = np.unique(X, axis=1, return_inverse=True)
    return uniq, np.asarray(inv).reshape(-1)

def win_fractions(W, X):
    """Fraction of bootstrap replicates in which each tree is (jointly) best. W: replicates x families resampling counts; X: families x trees scores (lower is
    better). Trees with identical score vectors form a tie group that wins or loses together and every member reports the group's fraction; ties between
    different groups split the win equally."""
    uniq, inv = groups(X)
    tot = W @ uniq
    best = tot == tot.min(axis=1, keepdims=True)
    frac = (best / best.sum(axis=1, keepdims=True)).sum(axis=0) / W.shape[0]
    return frac[inv]

def twins(X, names):
    """For each tree, the names of the other trees with an identical score vector."""
    _, inv = groups(X)
    return [[names[u] for u in range(len(names)) if inv[u] == inv[t] and u != t] for t in range(len(names))]
