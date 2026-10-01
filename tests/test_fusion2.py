"""fusion2 (X-066): each piece recovers a known answer from simulated data.

OU precision = inverse OU covariance; the smoother recovers a buoyant surface from noisy short-pair
increments better than chaining them; a 2π jump is found and undone; two tracks seen through different
incidences recover one vertical motion better than either; spatial coupling beats independent pixels on a
smooth field; σ and τ of an OU series are recovered.
"""
from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from insar_wetlands import fusion2 as f2

rng = np.random.default_rng(42)


def ou_path(t, tau, sigma, n=1):
    x = np.empty((n, len(t)))
    x[:, 0] = rng.normal(0, sigma, n)
    for k in range(1, len(t)):
        a = np.exp(-(t[k] - t[k - 1]) / tau)
        x[:, k] = a * x[:, k - 1] + rng.normal(0, sigma * np.sqrt(1 - a ** 2), n)
    return x


def pairs_for(t_idx, A, h, sd, px=0, jump=None):
    i, j = t_idx[:-1], t_idx[1:]
    y = A * (h[j] - h[i]) + rng.normal(0, sd, len(i))
    if jump is not None:
        y[jump] += f2.CYCLE_LOS_MM
    n = len(i)
    return f2.Pairs(np.full(n, px), i, j, np.full(n, A), y, np.full(n, sd))


def test_ou_precision_is_the_inverse_covariance():
    t = np.array([0, 6, 12, 30, 31, 50.0])
    tau, s = 20.0, 3.0
    C = s ** 2 * np.exp(-np.abs(t[:, None] - t[None, :]) / tau)
    np.testing.assert_allclose(f2.ou_precision(t, tau, s).toarray(), np.linalg.inv(C), rtol=1e-8, atol=1e-10)


def test_smoother_beats_chained_increments_and_undoes_a_cycle_jump():
    t = np.arange(0, 6 * 120, 6.0)
    W = 8 * np.sin(2 * np.pi * t / 365)                         # water-table anomaly, cm
    g = 2.0
    h = g * W + ou_path(t, 25, 6)[0]                             # vertical mm
    A, sd = 0.9 * np.cos(np.radians(32)), 3.0
    p = pairs_for(np.arange(len(t)), A, h, sd, jump=40)
    out = f2.smooth(p, g * W, f2.ou_precision(t, 25, 6), len(t), 1)
    est = out["h"][0]
    chained = np.r_[0, np.cumsum((p.y - np.where(np.arange(len(p.y)) == 40, f2.CYCLE_LOS_MM, 0)) / A)]  # even with the jump removed
    rmse = lambda a: np.sqrt(np.mean((a - a.mean() - (h - h.mean())) ** 2))  # noqa: E731
    assert out["k"][40] == -1 and np.count_nonzero(out["k"]) == 1
    assert rmse(est) < 0.6 * rmse(chained)


def test_two_tracks_recover_one_vertical_motion_better_than_either():
    t = np.arange(0, 400, 1.0)
    prior = 2.0 * 6 * np.sin(2 * np.pi * t / 365)
    Q = f2.ou_precision(t, 20, 5)
    asc, desc = np.arange(0, 400, 6), np.arange(2, 400, 6)
    Aa, Ad, sd = 0.85 * np.cos(np.radians(32.26)), 0.85 * np.cos(np.radians(39.17)), 3.0
    err = {"asc": [], "desc": [], "both": []}
    for _ in range(10):
        h = prior + ou_path(t, 20, 5)[0]
        pa, pd_ = pairs_for(asc, Aa, h, sd), pairs_for(desc, Ad, h, sd)
        both = f2.Pairs(*(np.r_[getattr(pa, k), getattr(pd_, k)] for k in ("px", "i", "j", "A", "y", "sd")))
        for name, p in (("asc", pa), ("desc", pd_), ("both", both)):
            d = f2.smooth(p, prior, Q, len(t), 1, ambiguity_iters=0)["h"][0] - h
            err[name].append(np.sqrt(np.mean((d - d.mean()) ** 2)))
    m = {k: np.mean(v) for k, v in err.items()}
    assert m["both"] < 0.95 * min(m["asc"], m["desc"])


def mat_case(sd, seed):
    rng_ = np.random.default_rng(seed)
    globals()["rng"] = rng_
    ny, nx, nt = 8, 8, 40
    t = np.arange(nt) * 6.0
    L, idx = f2.grid_laplacian(np.ones((ny, nx), bool))
    common = ou_path(t, 30, 6)[0]
    yy, xx = np.mgrid[0:ny, 0:nx]
    gain = (1 + 0.3 * np.sin(xx / 3.0) * np.cos(yy / 4.0)).ravel()
    h = gain[:, None] * common[None, :]                           # one motion, smoothly varying amplitude
    rows = [pairs_for(np.arange(nt), 0.9, h[k], sd, px=k) for k in range(len(idx))]
    p = f2.Pairs(*(np.concatenate([getattr(r, f) for r in rows]) for f in ("px", "i", "j", "A", "y", "sd")))
    return t, L, idx, h, p, nt


def test_shared_motion_plus_smooth_departures_beats_independent_pixels():
    errs = {"independent": [], "two-level": []}
    for seed in range(3):
        t, L, idx, h, p, nt = mat_case(6.0, seed)
        n = len(idx)
        for name, kw in (("independent", {"Qprior": f2.space_time_precision(f2.ou_precision(t, 30, 6), n)}),
                         ("two-level", {"Qprior": f2.space_time_precision(f2.ou_precision(t, 30, 1.5), n, L, 16.0),
                                        "Qcommon": f2.ou_precision(t, 30, 6)})):
            d = f2.smooth(p, np.zeros((n, nt)), n_t=nt, n_px=n, ambiguity_iters=0, **kw)["h"] - h
            d = d - d.mean(axis=1, keepdims=True)                 # increments fix shapes, not levels
            errs[name].append(np.sqrt(np.mean(d ** 2)))
    assert np.mean(errs["two-level"]) < 0.6 * np.mean(errs["independent"])


def test_shared_series_equals_one_pixel_with_pooled_noise():
    """With rigid departures (σ_d → 0) the shared series is a single pixel seen with noise sd/√n_px and
    its OWN prior — not a prior n_px times tighter (the flaw of tying independent pixel priors together)."""
    t, L, idx, h, p, nt = mat_case(6.0, 7)
    n = len(idx)
    p = f2.Pairs(p.px, p.i, p.j, p.A, np.tile(p.y.reshape(n, -1)[0], n), p.sd)     # every pixel sees pixel 0's data
    two = f2.smooth(p, np.zeros((n, nt)), f2.space_time_precision(f2.ou_precision(t, 30, 1e-3), n), nt, n,
                    ambiguity_iters=0, Qcommon=f2.ou_precision(t, 30, 6))
    one = f2.Pairs(np.zeros(nt - 1, int), np.arange(nt - 1), np.arange(1, nt), np.full(nt - 1, 0.9),
                   p.y[:nt - 1], np.full(nt - 1, 6.0 / np.sqrt(n)))
    single = f2.smooth(one, np.zeros(nt), f2.ou_precision(t, 30, 6), nt, 1, ambiguity_iters=0)["h"][0]
    np.testing.assert_allclose(two["h"].mean(axis=0), single, atol=1e-2)


def test_grid_laplacian_rows_sum_to_zero_and_count_neighbours():
    m = np.array([[1, 1, 0], [1, 1, 1]], bool)
    L, idx = f2.grid_laplacian(m)
    assert sp.issparse(L) and len(idx) == 5
    np.testing.assert_allclose(np.asarray(L.sum(axis=1)).ravel(), 0)
    assert sorted(np.diag(L.toarray()).astype(int).tolist()) == [1, 2, 2, 2, 3]


def test_ou_from_series_recovers_sigma_and_tau():
    t = np.arange(0, 3000, 1.0)
    taus, sigmas = [], []
    for _ in range(5):
        fit = f2.ou_from_series(t, ou_path(t, 15, 4)[0], max_lag_days=40)
        taus.append(fit["tau_days"])
        sigmas.append(fit["sigma"])
    assert np.mean(sigmas) == pytest.approx(4, rel=0.15)
    assert np.mean(taus) == pytest.approx(15, rel=0.25)


def test_posterior_sd_shrinks_where_observed():
    t = np.arange(0, 120, 6.0)
    Q = f2.ou_precision(t, 25, 6)
    p = pairs_for(np.arange(10), 0.9, np.zeros(len(t)), 1.0)       # pairs only over the first 10 dates
    out = f2.smooth(p, np.zeros(len(t)), Q, len(t), 1, ambiguity_iters=0)
    sd = f2.posterior_sd_pixel(out["P"], len(t), 0, out["n_common"])
    assert sd[5] < sd[-1] <= 6.0 + 1e-6


def test_conjugate_gradients_match_the_direct_solve():
    t, L, idx, h, p, nt = mat_case(4.0, 3)
    n = len(idx)
    kw = {"Qprior": f2.space_time_precision(f2.ou_precision(t, 30, 1.5), n, L, 16.0), "Qcommon": f2.ou_precision(t, 30, 6),
          "n_t": nt, "n_px": n, "ambiguity_iters": 0}
    a = f2.smooth(p, np.zeros((n, nt)), **kw)["h"]
    b = f2.smooth(p, np.zeros((n, nt)), solver="cg", **kw)["h"]
    np.testing.assert_allclose(a, b, atol=1e-5)


def test_ambiguities_decided_per_pixel_find_a_planted_cycle():
    t = np.arange(0, 6 * 80, 6.0)
    h = ou_path(t, 25, 6, n=3)
    rows = [pairs_for(np.arange(len(t)), 0.8, h[x], 2.5, px=x, jump=20 + 10 * x) for x in range(3)]
    p = f2.Pairs(*(np.concatenate([getattr(r, f) for r in rows]) for f in ("px", "i", "j", "A", "y", "sd")))
    k = f2.resolve_ambiguities_per_pixel(p, np.zeros((3, len(t))), f2.ou_precision(t, 25, 6), len(t))
    planted = np.concatenate([np.arange(len(t) - 1) == 20 + 10 * x for x in range(3)])
    assert np.all(k[planted] == -1) and np.count_nonzero(k) == 3


def test_pixelwise_sd_equals_the_exact_single_pixel_posterior():
    t = np.arange(0, 6 * 30, 6.0)
    Q = f2.ou_precision(t, 25, 6)
    rows = [pairs_for(np.arange(len(t)), 0.8, np.zeros(len(t)), 2.0 + x, px=x) for x in range(2)]
    p = f2.Pairs(*(np.concatenate([getattr(r, f) for r in rows]) for f in ("px", "i", "j", "A", "y", "sd")))
    sd = f2.pixelwise_sd(p, Q, len(t), 2)
    for x in range(2):
        q = p.subset(p.px == x)
        q = f2.Pairs(np.zeros(len(q.y), int), q.i, q.j, q.A, q.y, q.sd)
        P = f2.smooth(q, np.zeros(len(t)), Q, len(t), 1, ambiguity_iters=0)["P"]
        np.testing.assert_allclose(sd[x], f2.posterior_sd_pixel(P, len(t), 0), rtol=1e-8)
    assert np.all(sd[1] > sd[0])                                   # noisier pixel, wider uncertainty
