"""Fusion v2 (X-066, branch fusion-x062): the mat's vertical surface height, every pixel and every date of both
tracks, as the exact posterior of one linear-Gaussian model — solved as a sparse system, not a filter.

Unknowns: h(x, t), vertical height (mm) at every acquisition date t of either track, for every pixel x of a mask.

Prior (physics):   h(x, t) = g · W(t) + c(t) + d(x, t)
    W   water-table anomaly forcing (cm; field loggers), g buoyancy (mm per cm, fitted on the laser);
    c   the motion the whole mat shares beyond the forcing: Ornstein–Uhlenbeck in time (τ, σ_c) —
        mean-reverting, so no random-walk drift;
    d   each pixel's departure from it: OU in time (σ_d, smaller) and smooth in space (graph Laplacian, κ).
    The two levels matter: tying independent per-pixel priors together instead would give the shared motion a
    prior n_px times too tight and pull it toward zero (the test suite checks the equivalence).
Observations (radar), one per short pair (i, j) of one track at pixel x:
    s = m_p · cos θ_p · (h_j − h_i) + e · (W_j − W_i) + n,   n ~ N(0, σ_p²)
    m_p the fraction of the motion the pair records (depends on revisit and track; X-065), σ_p from coherence.
Solve: the posterior mean of r is P⁻¹ b with P = prior precision + Dᵀ W D (D: pair differences) — the exact
Kalman/RTS smoother for this model. Phase ambiguities (multiples of λ/2 in line of sight, |k| ≤ 1) are resolved by
model comparison: a cycle is kept only when it lowers the posterior objective clearly (``smooth``).

No network inversion: each pixel's heights are tied by its own short pairs and by the prior.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

WAVELENGTH_MM = 55.465763
CYCLE_LOS_MM = WAVELENGTH_MM / 2           # one 2π of interferometric phase, line of sight


# ------------------------------------------------------------------------------ priors

def ou_precision(t_days, tau_days: float, sigma: float) -> sp.csr_matrix:
    """Precision matrix of a stationary Ornstein–Uhlenbeck process sampled at increasing times ``t_days``
    (variance σ², correlation exp(−Δt/τ)). Tridiagonal: the Markov property of the OU process."""
    t = np.asarray(t_days, float)
    n = len(t)
    a = np.exp(-np.diff(t) / tau_days)                       # AR coefficients between samples
    q = sigma ** 2 * (1 - a ** 2)                            # innovation variances
    main = np.zeros(n)
    main[0] += 1 / sigma ** 2
    main[1:] += 1 / q
    main[:-1] += a ** 2 / q
    off = -a / q
    return sp.diags([main, off, off], [0, 1, -1], format="csr")


def grid_laplacian(mask: np.ndarray) -> tuple[sp.csr_matrix, np.ndarray]:
    """4-neighbour graph Laplacian over the True cells of ``mask``; returns (L, flat indices of the cells)."""
    idx = np.flatnonzero(mask.ravel())
    pos = -np.ones(mask.size, int)
    pos[idx] = np.arange(len(idx))
    ny, nx = mask.shape
    rows, cols = [], []
    for k, f in enumerate(idx):
        r, c = divmod(f, nx)
        for rr, cc in ((r + 1, c), (r, c + 1)):
            if rr < ny and cc < nx and mask[rr, cc]:
                rows.append(k)
                cols.append(pos[rr * nx + cc])
    n = len(idx)
    A = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    A = (A + A.T).tocsr()
    return (sp.diags(np.asarray(A.sum(axis=1)).ravel()) - A).tocsr(), idx


def space_time_precision(Qt: sp.spmatrix, n_px: int, L: sp.spmatrix | None = None, kappa: float = 0.0) -> sp.csr_matrix:
    """Separable prior precision for unknowns ordered pixel-major (pixel · n_t + t): (I + κL) ⊗ Q_t."""
    Qs = sp.identity(n_px, format="csr") if L is None or kappa == 0 else (sp.identity(n_px) + kappa * L).tocsr()
    return sp.kron(Qs, Qt, format="csr")


# ------------------------------------------------------------------------------ observations

@dataclass
class Pairs:
    """Short-pair observations. ``px`` pixel index (0..n_px−1), ``i``/``j`` date indices into the joint date
    axis, ``A`` = m·cos θ (LOS mm per mm of vertical change), ``y`` observed LOS change (mm) minus the moisture
    term e·ΔW, ``sd`` its noise (mm)."""
    px: np.ndarray
    i: np.ndarray
    j: np.ndarray
    A: np.ndarray
    y: np.ndarray
    sd: np.ndarray

    def subset(self, keep: np.ndarray) -> Pairs:
        return Pairs(*(getattr(self, f)[keep] for f in ("px", "i", "j", "A", "y", "sd")))


def difference_operator(p: Pairs, n_t: int, n_px: int) -> sp.csr_matrix:
    """D (n_obs × n_px·n_t): row k = A_k·(e_{px,j} − e_{px,i})."""
    k = np.arange(len(p.y))
    rows = np.r_[k, k]
    cols = np.r_[p.px * n_t + p.j, p.px * n_t + p.i]
    vals = np.r_[p.A, -p.A]
    return sp.csr_matrix((vals, (rows, cols)), shape=(len(p.y), n_px * n_t))


# ------------------------------------------------------------------------------ the solve

def block_preconditioner(P: sp.spmatrix, n_c: int, n_t: int):
    """Solve with the block-diagonal part of P — the shared series, then one n_t × n_t time block per pixel.
    The blocks do not fill one another, so this factorises fast where P itself does not (space × time)."""
    P = P.tocoo()
    blk = lambda i: np.where(i < n_c, -1, (i - n_c) // n_t)  # noqa: E731
    keep = blk(P.row) == blk(P.col)
    B = sp.csc_matrix((P.data[keep], (P.row[keep], P.col[keep])), shape=P.shape)
    return spla.factorized(B)


def cg_solver(P: sp.spmatrix, n_c: int, n_t: int, rtol: float = 1e-9):
    """``solve(b)`` by preconditioned conjugate gradients (P is symmetric positive definite)."""
    pre = block_preconditioner(P, n_c, n_t)
    M = spla.LinearOperator(P.shape, matvec=pre)
    P = P.tocsr()

    def solve(b):
        x, info = spla.cg(P, b, rtol=rtol, maxiter=5000, M=M)
        if info != 0:
            raise RuntimeError(f"conjugate gradients did not converge (info {info})")
        return x
    return solve


def resolve_ambiguities_per_pixel(p: Pairs, prior_mean: np.ndarray, Qsingle: sp.spmatrix, n_t: int,
                                  min_gain: float = 9.0) -> np.ndarray:
    """Phase ambiguities decided pixel by pixel (``smooth`` on each pixel's own pairs, with a single-pixel
    prior of the combined variance) — small direct solves, so the model comparison stays cheap; the mat-wide
    solve then takes them as given. Returns k aligned with ``p``."""
    k = np.zeros(len(p.y))
    mu = np.asarray(prior_mean, float)
    for x in np.unique(p.px):
        sel = np.flatnonzero(p.px == x)
        q = p.subset(sel)
        q = Pairs(np.zeros(len(sel), int), q.i, q.j, q.A, q.y, q.sd)
        k[sel] = smooth(q, mu[x], Qsingle, n_t, 1, min_gain=min_gain)["k"]
    return k


def smooth(p: Pairs, prior_mean: np.ndarray, Qprior: sp.spmatrix, n_t: int, n_px: int,
           ambiguity_iters: int = 3, min_gain: float = 9.0, max_flips: int = 50,
           Qcommon: sp.spmatrix | None = None, solver: str = "direct") -> dict:
    """Posterior mean of h (n_px × n_t) given the pairs, a prior mean (n_px × n_t) and the prior precision of
    the residual, with phase ambiguities resolved by model comparison.

    For ambiguities k (one per pair, |k| ≤ 1, in cycles of λ/2 LOS) the best residual has the closed-form
    objective J(k) = y'ᵀWy' − bᵀP⁻¹b, y' = y + k·cycle − Dμ, b = DᵀWy' (the data misfit plus the prior penalty,
    minimised over h). Starting from k = 0, the pairs with the largest residuals are tried one by one and a
    cycle is kept only if it lowers J by more than ``min_gain`` (in χ² units: 9 = a 3-σ improvement), so a
    real large motion that the prior and the neighbouring pairs support is restored, and noise is not.
    ``Qprior`` is the precision of the pixel departures d (n_px·n_t, pixel-major); with ``Qcommon`` (n_t × n_t)
    a shared series c is estimated too (unknowns [c, d], h = μ + c + d). ``solver="cg"`` for the mat-wide
    system (a direct factorisation of space × time fills in badly); use it with ``ambiguity_iters=0`` after
    ``resolve_ambiguities_per_pixel``.
    Returns h, the shared c (or None), k, the final pair residuals, and P with its layout (for variances)."""
    mu = np.asarray(prior_mean, float).reshape(-1)
    D = difference_operator(p, n_t, n_px)
    w = 1.0 / p.sd ** 2
    base = p.y - D @ mu                                       # data minus the prior mean (physics) term
    n_c = 0 if Qcommon is None else n_t
    if n_c:
        M = sp.kron(sp.csr_matrix(np.ones((n_px, 1))), sp.identity(n_t), format="csr")   # c → every pixel
        D = sp.hstack([D @ M, D], format="csr")
        Qprior = sp.block_diag([Qcommon, Qprior], format="csr")
    P = (Qprior + D.T @ sp.diags(w) @ D).tocsc()
    solve = spla.factorized(P) if solver == "direct" else cg_solver(P, n_c, n_t)

    def fit(k):
        yk = base + k * CYCLE_LOS_MM
        b = D.T @ (w * yk)
        z = solve(b)
        return z, float(yk @ (w * yk) - b @ z), yk - D @ z

    k = np.zeros(len(p.y))
    r, J, res = fit(k)
    for _ in range(max(ambiguity_iters, 0)):
        changed = False
        for c in np.argsort(-np.abs(res / p.sd))[:max_flips]:
            if abs(res[c]) < 3 * p.sd[c] and abs(res[c]) < CYCLE_LOS_MM / 4:
                break
            trial = k.copy()
            trial[c] = np.clip(k[c] - np.sign(res[c]), -1, 1)
            if trial[c] == k[c]:
                continue
            r2, J2, res2 = fit(trial)
            if J - J2 > min_gain:
                k, r, J, res, changed = trial, r2, J2, res2, True
        if not changed:
            break
    c = r[:n_c] if n_c else None
    d = r[n_c:]
    h = (mu + d + (np.tile(c, n_px) if n_c else 0)).reshape(n_px, n_t)
    return {"h": h, "c": c, "d": d.reshape(n_px, n_t), "k": k, "residual": res, "objective": J, "P": P, "n_common": n_c}


def posterior_sd_pixel(P: sp.spmatrix, n_t: int, px: int, n_common: int = 0) -> np.ndarray:
    """Posterior standard deviation of pixel ``px``'s n_t heights (h = c + d when a shared series is
    estimated): the quadratic form eᵀP⁻¹e for each date, by n_t solves with one factorisation."""
    solve = spla.factorized(P.tocsc())
    sd = np.empty(n_t)
    for t in range(n_t):
        e = np.zeros(P.shape[0])
        e[n_common + px * n_t + t] = 1.0
        if n_common:
            e[t] = 1.0
        sd[t] = np.sqrt(max(float(e @ solve(e)), 0.0))
    return sd


def pixelwise_sd(p: Pairs, Qsingle: sp.spmatrix, n_t: int, n_px: int) -> np.ndarray:
    """Posterior SD of every pixel's heights from its own pairs and a single-pixel prior (n_px × n_t): the
    diagonal of each pixel's small n_t × n_t inverse precision. It leaves out what neighbouring pixels add
    through the spatial coupling, so it is an upper bound on the mat model's uncertainty (conservative)."""
    Qd = Qsingle.toarray()
    out = np.empty((n_px, n_t))
    order = np.argsort(p.px, kind="stable")
    bounds = np.searchsorted(p.px[order], np.arange(n_px + 1))
    for x in range(n_px):
        sel = order[bounds[x]:bounds[x + 1]]
        P = Qd.copy()
        w = 1.0 / p.sd[sel] ** 2
        a, i, j = p.A[sel], p.i[sel], p.j[sel]
        np.add.at(P, (j, j), w * a * a)
        np.add.at(P, (i, i), w * a * a)
        np.add.at(P, (i, j), -w * a * a)
        np.add.at(P, (j, i), -w * a * a)
        out[x] = np.sqrt(np.clip(np.diag(np.linalg.inv(P)), 0, None))
    return out


# ------------------------------------------------------------------------------ fitting from the laser

def ou_from_series(t_days, x, max_lag_days: float = 60.0) -> dict:
    """σ and τ of an OU process from an irregular series (e.g. the daily laser residual): σ from the variance,
    τ from the decay of the lag-binned autocorrelation (least squares on log ρ for 0 < ρ). Pairwise in the
    length of the series — pass daily values, not hourly."""
    t, x = np.asarray(t_days, float), np.asarray(x, float)
    ok = np.isfinite(x)
    t, x = t[ok], x[ok] - np.nanmean(x[ok])
    var = float(np.mean(x ** 2))
    dt = np.abs(t[:, None] - t[None, :])
    prod = x[:, None] * x[None, :]
    bins = np.arange(1, max_lag_days + 2, 2.0)
    lags, rho = [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        m = (dt >= lo) & (dt < hi)
        if m.sum() > 20:
            lags.append(dt[m].mean())
            rho.append(prod[m].mean() / var)
    lags, rho = np.array(lags), np.array(rho)
    use = rho > 0.2                                          # the weak long-lag tail biases τ upward
    tau = float(-1 / np.polyfit(lags[use], np.log(rho[use]), 1)[0]) if use.sum() >= 3 else np.nan
    return {"sigma": float(np.sqrt(var)), "tau_days": tau, "lags": lags, "rho": rho}
