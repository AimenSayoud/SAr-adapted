"""The footprint routines recover known answers (X-133)."""
import numpy as np
import pytest

from insar_wetlands import footprint as F


def test_pearson_and_slope_cols_match_direct_fits_with_gaps():
    rng = np.random.default_rng(0)
    x = rng.normal(size=60)
    Y = np.column_stack([2.0 * x + rng.normal(0, 0.5, 60), -x + rng.normal(0, 2, 60), rng.normal(size=60)])
    Y[[3, 10, 40], 0] = np.nan                              # gaps in one cell only
    r = F.pearson_cols(x, Y)
    s = F.slope_origin_cols(x, Y)
    for j in range(3):
        m = np.isfinite(Y[:, j])
        assert r[j] == pytest.approx(np.corrcoef(x[m], Y[m, j])[0, 1], abs=1e-12)
        assert s[j] == pytest.approx(np.sum(x[m] * Y[m, j]) / np.sum(x[m] ** 2), abs=1e-12)
    assert s[0] == pytest.approx(2.0, abs=0.1)
    few = np.full((60, 1), np.nan)
    few[:4, 0] = 1.0
    assert np.isnan(F.pearson_cols(x, few)[0])                # under min_n pairs: no value


def test_permuted_r_keeps_groups_and_breaks_the_link():
    rng = np.random.default_rng(1)
    groups = np.repeat(np.arange(6), 10)
    season = np.repeat(rng.normal(0, 3, 6), 10)               # a shared seasonal level per group
    x = season + rng.normal(size=60)
    Y = np.column_stack([season + rng.normal(size=60),       # follows only the season
                         x + rng.normal(0, 0.3, 60)])        # follows x pair by pair
    P = F.permuted_r(x, Y, groups, n_perm=300, seed=0)
    obs = F.pearson_cols(x, Y)
    p = (1 + (P >= obs).sum(0)) / (1 + len(P))
    assert p[0] > 0.05                                       # the season alone is not a pair-by-pair match
    assert p[1] < 0.01
    # within-group permutation keeps each group's values
    xp_like = F.permuted_r(np.arange(60.0), np.arange(60.0)[:, None], groups, 1, 0)
    assert np.isfinite(xp_like).all()


def test_longest_run():
    d = np.datetime64("2021-06-01")
    t1 = np.array([d, d + 6, d + 12, d + 30, d + 36, d + 42, d + 48])
    t2 = np.array([d + 6, d + 12, d + 18, d + 36, d + 42, d + 48, d + 54])
    assert F.longest_run(t1, t2).tolist() == [3, 4, 5, 6]     # a missing acquisition breaks the first run
    perm = np.array([6, 0, 4, 2, 5, 1, 3])
    assert sorted(perm[F.longest_run(t1[perm], t2[perm])].tolist()) == [3, 4, 5, 6]
    assert F.longest_run(t1[:3], t2[:3]).tolist() == [0, 1, 2]


def test_resample_index():
    counts = np.array([2, 0, 1])                             # block 0 twice, block 1 not, block 2 once
    unit_pos = np.array([0, 0, 1, 2])
    assert F.resample_index(counts, unit_pos).tolist() == [0, 0, 1, 1, 3]
