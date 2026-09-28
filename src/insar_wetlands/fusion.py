"""Hybrid fusion of field data and Sentinel-1 over the floating mat (X-062, branch fusion-x062).

Five blocks, each small enough to test on synthetic ground truth:

A. **Water model** — a daily water-balance ("bucket") model per plot,
   ΔWTD = a·P − b·PET − c·(WTD − w0), fitted by least squares on the logger series; Hamon PET
   from air temperature and day length.
B. **Zones** — k-means on per-pixel behaviour features (k by silhouette) and a dissimilarity
   index to the plots (after Meyer & Pebesma's area of applicability): where the plots are
   representative and where nothing we measured applies.
C. **Spatial transfer** — Gaussian-process regression of plot traits on pixel features,
   judged leave-one-plot-out; predictions carry their own standard deviation.
D. **Phase physics** — the 12-day LOS change as motion + moisture:
   s1 = m·laser_LOS + e·ΔWTD + c, coherence-weighted Bayesian regression (phase noise from the
   Cramér–Rao bound for an L-look interferogram), with a fitted noise scale.
E. **Fusion** — per pixel and pair, the surface increment the water forcing predicts and the one
   the radar observes (moisture removed, scaled by m) are combined by inverse-variance weighting:
   the Kalman update of a random walk observed through its increments, with innovation gating
   against unwrapping errors. Increments are chained, never inverted over a network.

Nothing here holds data; the caller reads the field data (hub) and the interferograms (Drive).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

WAVELENGTH_MM = 55.465763


# ------------------------------------------------------------------ A. water model

def day_length_h(doy, lat_deg: float) -> np.ndarray:
    """Astronomical day length (hours) for day of year and latitude."""
    doy = np.asarray(doy, float)
    decl = 0.409 * np.sin(2 * np.pi * doy / 365 - 1.39)
    x = np.clip(-np.tan(np.radians(lat_deg)) * np.tan(decl), -1, 1)
    return 24 / np.pi * np.arccos(x)


def hamon_pet(temp_c, doy, lat_deg: float) -> np.ndarray:
    """Hamon potential evapotranspiration (mm/day) from daily mean air temperature (°C)."""
    t = np.asarray(temp_c, float)
    es = 6.108 * np.exp(17.26939 * t / (t + 237.3))            # saturation vapour pressure, hPa
    rho = 216.7 * es / (t + 273.3)                             # saturated vapour density, g/m³
    return np.maximum(0.1651 * (day_length_h(doy, lat_deg) / 12) * rho * 1.2, 0)


def fit_bucket(wtd, rain, pet, valid=None) -> dict:
    """Least-squares fit of ΔWTD_t = a·rain_t − b·pet_t − c·(WTD_t − w0) (daily, WTD in cm,
    rain and PET in mm). ``valid``: days usable as a start (e.g. not at the well floor)."""
    w, p, e = (np.asarray(v, float) for v in (wtd, rain, pet))
    dw = w[1:] - w[:-1]
    ok = np.isfinite(dw) & np.isfinite(p[:-1]) & np.isfinite(e[:-1]) & np.isfinite(w[:-1])
    if valid is not None:
        v = np.asarray(valid, bool)
        ok &= v[:-1] & v[1:]
    A = np.column_stack([p[:-1], -e[:-1], -w[:-1], np.ones(len(dw))])[ok]
    coef, *_ = np.linalg.lstsq(A, dw[ok], rcond=None)
    a, b, c, k = coef
    return {"a": float(a), "b": float(b), "c": float(c), "w0": float(k / c) if c != 0 else np.nan, "n": int(ok.sum())}


def simulate_bucket(par: dict, w_start: float, rain, pet, floor: float | None = None) -> np.ndarray:
    """Run the fitted bucket forward from ``w_start``; optionally hold it above a well floor."""
    p, e = np.asarray(rain, float), np.asarray(pet, float)
    w = np.empty(len(p))
    w[0] = w_start
    for t in range(len(p) - 1):
        w[t + 1] = w[t] + par["a"] * p[t] - par["b"] * e[t] - par["c"] * (w[t] - par["w0"])
        if floor is not None:
            w[t + 1] = max(w[t + 1], floor)
    return w


def nse(obs, sim) -> float:
    """Nash–Sutcliffe efficiency (1 = perfect, 0 = as good as the mean, < 0 worse)."""
    o, s = np.asarray(obs, float), np.asarray(sim, float)
    ok = np.isfinite(o) & np.isfinite(s)
    return float(1 - np.sum((o[ok] - s[ok]) ** 2) / np.sum((o[ok] - o[ok].mean()) ** 2))


# ------------------------------------------------------------------ B. zones

def standardize(X: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mu, sd = np.nanmean(X, axis=0), np.nanstd(X, axis=0)
    sd = np.where(sd > 0, sd, 1.0)
    return (X - mu) / sd, mu, sd


def choose_kmeans(Z: np.ndarray, ks=range(3, 9), sample: int = 3000, seed: int = 0) -> tuple[np.ndarray, int, dict]:
    """k-means on standardised features, k chosen by the silhouette score on a sample."""
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(Z), min(sample, len(Z)), replace=False)
    scores, fits = {}, {}
    for k in ks:
        km = KMeans(n_clusters=k, n_init=10, random_state=seed).fit(Z)
        fits[k] = km
        scores[k] = float(silhouette_score(Z[idx], km.labels_[idx]))
    best = max(scores, key=scores.get)
    return fits[best].labels_, best, scores


def dissimilarity_index(Z: np.ndarray, train: np.ndarray) -> tuple[np.ndarray, float]:
    """Distance of every row of Z to its nearest training row, divided by the mean distance
    between training rows; and the threshold beyond which a pixel is outside the area the
    training data represent (upper whisker of the training rows' leave-one-out index)."""
    T = Z[train]
    d_tt = np.sqrt(((T[:, None, :] - T[None, :, :]) ** 2).sum(-1))
    n = len(T)
    mean_d = d_tt[~np.eye(n, dtype=bool)].mean()
    d_all = np.sqrt(((Z[:, None, :] - T[None, :, :]) ** 2).sum(-1)).min(axis=1)
    di = d_all / mean_d
    np.fill_diagonal(d_tt, np.inf)
    loo = d_tt.min(axis=1) / mean_d
    q1, q3 = np.percentile(loo, [25, 75])
    return di, float(q3 + 1.5 * (q3 - q1))


# ------------------------------------------------------------------ C. spatial transfer

def gp_transfer(Xtrain: np.ndarray, y: np.ndarray, Xall: np.ndarray) -> dict:
    """Gaussian-process regression (RBF + white noise, standardised inputs) of a plot trait on
    pixel features: leave-one-out predictions and skill against the mean, and the prediction
    (mean, sd) for every pixel from the fit on all training points."""
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel

    def model():
        k = ConstantKernel(1.0, (1e-2, 1e2)) * RBF(np.ones(Xtrain.shape[1]), (0.3, 30.0)) + WhiteKernel(0.1, (1e-3, 10.0))
        return GaussianProcessRegressor(k, normalize_y=True, n_restarts_optimizer=3, random_state=0)

    Z, mu, sd = standardize(Xtrain)
    Zall = (Xall - mu) / sd
    n = len(y)
    loo = np.empty(n)
    for i in range(n):
        m = np.arange(n) != i
        loo[i] = model().fit(Z[m], y[m]).predict(Z[i:i + 1])[0]
    base = np.array([y[np.arange(n) != i].mean() for i in range(n)])
    rmse, rmse0 = float(np.sqrt(np.mean((loo - y) ** 2))), float(np.sqrt(np.mean((base - y) ** 2)))
    full = model().fit(Z, y)
    mean, std = full.predict(Zall, return_std=True)
    return {"loo_pred": loo, "rmse_loo": rmse, "rmse_mean_baseline": rmse0,
            "skill": 1 - rmse ** 2 / rmse0 ** 2, "mean": mean, "sd": std, "kernel": str(full.kernel_)}


# ------------------------------------------------------------------ D. phase physics

def phase_sigma_mm(coh, looks: int = 20) -> np.ndarray:
    """1-σ LOS noise (mm) of a multilooked interferogram from its coherence (Cramér–Rao bound)."""
    g = np.clip(np.asarray(coh, float), 0.05, 0.99)
    sig_rad = np.sqrt(1 - g ** 2) / (g * np.sqrt(2 * looks))
    return sig_rad * WAVELENGTH_MM / (4 * np.pi)


def bayes_linear(X: np.ndarray, y: np.ndarray, sigma: np.ndarray, prior_sd: float = 10.0,
                 iters: int = 20) -> dict:
    """Bayesian linear regression y = X·β + ε, ε_i ~ N(0, (s·σ_i)²), β ~ N(0, prior_sd²).
    The noise scale s is fitted (σ from coherence sets the relative weights, s the level).
    Returns posterior mean and covariance of β, and s."""
    X, y, sigma = np.asarray(X, float), np.asarray(y, float), np.asarray(sigma, float)
    s2 = 1.0
    for _ in range(iters):
        W = 1 / (s2 * sigma ** 2)
        P = X.T @ (W[:, None] * X) + np.eye(X.shape[1]) / prior_sd ** 2
        cov = np.linalg.inv(P)
        beta = cov @ (X.T @ (W * y))
        r = y - X @ beta
        s2_new = max(float(np.mean(r ** 2 / sigma ** 2)), 1e-6)
        if abs(s2_new - s2) < 1e-6 * s2:
            s2 = s2_new
            break
        s2 = s2_new
    return {"beta": beta, "cov": cov, "sd": np.sqrt(np.diag(cov)), "noise_scale": float(np.sqrt(s2))}


# ------------------------------------------------------------------ E. fusion

def fuse_increments(u, Q, z, R, gate: float = 3.0) -> dict:
    """Inverse-variance fusion of a forcing increment u (variance Q) and an observed increment
    z (variance R), per step — the Kalman update of a random-walk surface observed through its
    increments. Steps where the innovation |z − u| exceeds ``gate``·√(Q+R) are treated as
    radar outliers (e.g. an unwrapping error) and keep the forcing alone. Missing z: forcing.
    Arrays of shape (steps, ...) broadcast; returns fused increment, its variance, the gate
    mask, and the cumulative position with its variance (first date = 0)."""
    u, Q, z, R = (np.asarray(v, float) for v in (u, Q, z, R))
    has = np.isfinite(z) & np.isfinite(R)
    innov = np.where(has, z - u, 0.0)
    rejected = has & (np.abs(innov) > gate * np.sqrt(Q + np.where(has, R, 0)))
    use = has & ~rejected
    K = np.where(use, Q / (Q + np.where(use, R, 1.0)), 0.0)
    inc = u + K * innov
    var = (1 - K) * Q
    cum = np.concatenate([np.zeros((1,) + inc.shape[1:]), np.cumsum(inc, axis=0)])
    cvar = np.concatenate([np.zeros((1,) + var.shape[1:]), np.cumsum(var, axis=0)])
    return {"inc": inc, "var": var, "rejected": rejected, "used": use, "cum": cum, "cum_var": cvar}


def harmonic_amplitude_map(t_years: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """Annual semi-amplitude of every column of Y (dates × pixels) after a linear trend."""
    M = np.column_stack([np.cos(2 * np.pi * t_years), np.sin(2 * np.pi * t_years), t_years, np.ones_like(t_years)])
    out = np.full(Y.shape[1], np.nan)
    ok_cols = np.isfinite(Y).sum(axis=0) >= 8
    for j in np.flatnonzero(ok_cols):
        m = np.isfinite(Y[:, j])
        b, *_ = np.linalg.lstsq(M[m], Y[m, j], rcond=None)
        out[j] = np.hypot(b[0], b[1])
    return out


def interval_series(dates, values: pd.Series) -> np.ndarray:
    """Change of a series between consecutive ``dates`` (value at the later minus the earlier)."""
    v = values.reindex(pd.to_datetime(dates)).to_numpy(float)
    return v[1:] - v[:-1]
