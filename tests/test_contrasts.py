"""Matched contrasts and the block bootstrap recover known answers (X-132)."""
import numpy as np
import pandas as pd
import pytest

from insar_wetlands import contrasts as C


def test_match_nearest_caliper_order_ties_and_no_replacement():
    t = pd.DataFrame({"depth": [10.0, 12.0, 100.0], "dp6": [0.0, 50.0, 10.0]}, index=["t0", "t1", "t2"])
    c = pd.DataFrame({"depth": [11.0, 11.0, 60.0], "dp6": [40.0, 5.0, 0.0]}, index=["c0", "c1", "c2"])
    m = C.match_nearest(t, c, on="depth", caliper=40.0, tie="dp6", order="dp6")
    # order by dp6: t0 (0), t2 (10), t1 (50). t0: c0 and c1 tie on |Δdepth| = 1 → c1 (|Δdp6| 5 < 40).
    # t2: |Δdepth| to c2 = 40 (inside the caliper). t1: c1 used → c0.
    assert m.set_index("treated").control.to_dict() == {"t0": "c1", "t2": "c2", "t1": "c0"}
    m2 = C.match_nearest(t, c, on="depth", caliper=39.0, tie="dp6", order="dp6")
    assert "t2" not in set(m2.treated)                         # beyond the caliper: unmatched, not forced
    assert m.control.is_unique


def test_weighted_median_is_the_median_of_the_repeated_values():
    rng = np.random.default_rng(1)
    for _ in range(50):
        x = rng.normal(size=rng.integers(1, 12))
        w = rng.integers(0, 4, len(x))
        if w.sum() == 0:
            assert np.isnan(C.weighted_median(x, w))
            continue
        assert C.weighted_median(x, w) == pytest.approx(np.median(np.repeat(x, w)))


def test_block_counts_and_weights():
    ub, counts = C.block_counts([3, 1, 3, 7], n_boot=500, seed=0)
    assert ub.tolist() == [1, 3, 7] and counts.shape == (500, 3)
    assert (counts.sum(1) == 3).all()                           # as many draws as distinct blocks
    assert counts.mean() == pytest.approx(1.0, abs=0.1)
    w = C.unit_weights([3, 1, 3, 7], ub, counts)
    assert (w[:, 0] == w[:, 2]).all() and (w[:, 0] == counts[:, 1]).all()
    with pytest.raises(ValueError):
        C.unit_weights([2], ub, counts)


def test_block_bootstrap_covers_a_known_shift_and_keeps_shared_draws():
    rng = np.random.default_rng(0)
    blocks = np.repeat(np.arange(40), 5)
    block_effect = rng.normal(0, 0.05, 40)[blocks]               # spatially correlated noise
    on = 0.10 + block_effect + rng.normal(0, 0.02, 200)          # true contrast 0.10
    off = 0.02 + block_effect + rng.normal(0, 0.02, 200)         # true contrast 0.02 → interaction 0.08
    ub, counts = C.block_counts(blocks, 2000, seed=0)
    W = C.unit_weights(blocks, ub, counts)
    d_on = np.array([C.weighted_median(on, w) for w in W])
    d_off = np.array([C.weighted_median(off, w) for w in W])
    lo, hi = C.boot_interval(d_on)
    assert lo < 0.10 < hi
    ilo, ihi = C.boot_interval(d_on - d_off)                      # the shared block noise cancels in the difference
    assert ilo < 0.08 < ihi and (ihi - ilo) < 0.75 * (hi - lo)


def test_ols_from_block_sums_equals_a_direct_fit():
    rng = np.random.default_rng(2)
    n = 3000
    X = np.c_[np.ones(n), rng.normal(size=(n, 3))]
    y = X @ np.array([0.5, 1.0, -0.3, 0.0]) + rng.normal(0, 0.5, n)
    blocks = rng.integers(0, 30, n)
    S = C.block_sums(X, y, blocks)
    one = np.ones(len(S["blocks"]))
    beta, r2 = C.ols_from_sums(np.tensordot(one, S["xtx"], 1), one @ S["xty"], one @ S["yty"], one @ S["sy"], one @ S["n"])
    ref = np.linalg.lstsq(X, y, rcond=None)[0]
    assert beta == pytest.approx(ref, abs=1e-10)
    r2_ref = 1 - np.sum((y - X @ ref) ** 2) / np.sum((y - y.mean()) ** 2)
    assert r2 == pytest.approx(r2_ref, abs=1e-10)
    # leaving column 1 out equals a refit without it
    keep = [0, 2, 3]
    b2, r2b = C.ols_from_sums(np.tensordot(one, S["xtx"], 1), one @ S["xty"], one @ S["yty"], one @ S["sy"], one @ S["n"], keep)
    ref2 = np.linalg.lstsq(X[:, keep], y, rcond=None)[0]
    assert b2 == pytest.approx(ref2, abs=1e-10)
    assert r2b == pytest.approx(1 - np.sum((y - X[:, keep] @ ref2) ** 2) / np.sum((y - y.mean()) ** 2), abs=1e-10)


def test_ols_block_bootstrap_drops_and_intervals():
    rng = np.random.default_rng(3)
    n = 6000
    X = np.c_[np.ones(n), rng.normal(size=(n, 2))]
    y = 2.0 * X[:, 1] + rng.normal(0, 1, n)                       # factor a matters, factor b does not
    S = C.block_sums(X, y, rng.integers(0, 60, n))
    out = C.ols_block_bootstrap(S, {"a": [1], "b": [2]}, n_boot=300, seed=0)
    assert out["drop"]["a"] == pytest.approx(out["r2"], abs=0.01)  # a carries all the explained variance
    assert abs(out["drop"]["b"]) < 1e-3
    lo, hi = C.boot_interval(out["beta_draws"][:, 1])
    assert lo < 2.0 < hi
    lo, hi = C.boot_interval(out["beta_draws"][:, 2])
    assert lo < 0.0 < hi
