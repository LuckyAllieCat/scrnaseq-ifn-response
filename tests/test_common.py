"""Unit tests for the pure helper functions in scripts/common.py."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from common import choose_n_pcs, mad_bounds, neighbor_mixing, pseudobulk  # noqa: E402


def test_mad_bounds_symmetric_and_robust():
    vals = np.r_[np.full(99, 10.0), 1000.0]  # one huge outlier barely moves the bounds
    lo, hi = mad_bounds(vals, 3)
    assert lo == hi == 10.0
    lo, hi = mad_bounds([1, 2, 3, 4, 5], 2)
    assert np.isclose(10 - hi, lo - (-4))  # symmetric around the median (3)


def test_choose_n_pcs_finds_elbow():
    var = np.r_[[10, 6, 4, 3], np.full(20, 0.5)]
    assert choose_n_pcs(var, min_pcs=1) in (4, 5)  # where the curve flattens
    assert choose_n_pcs(var, min_pcs=8) == 8


def test_neighbor_mixing_extremes():
    labels = np.array(["a"] * 5 + ["b"] * 5)
    separated = sparse.block_diag([np.ones((5, 5)), np.ones((5, 5))])
    mixed = sparse.csr_matrix(np.ones((10, 10)))
    assert neighbor_mixing(separated, labels) < 0.01
    assert neighbor_mixing(mixed, labels) > 0.99


def test_pseudobulk_sums_counts_per_group():
    counts = sparse.csr_matrix(np.array([[1, 0], [2, 1], [0, 5], [3, 3]]))
    obs = pd.DataFrame({"patient": ["p1", "p1", "p2", "p2"], "condition": ["c", "c", "c", "s"]})
    summed, meta = pseudobulk(counts, obs, ["patient", "condition"])
    assert summed.loc["p1|c"].tolist() == [3, 1]
    assert summed.loc["p2|s"].tolist() == [3, 3]
    assert meta.loc["p1|c", "n_cells"] == 2
    assert summed.values.sum() == counts.sum()  # no counts lost
