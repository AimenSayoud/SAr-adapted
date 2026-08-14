"""Tests written in response to two referee reports.

Each function here answers one specific objection. The objection is named in the
docstring, because a test whose motivation is lost is a test nobody can judge
later.

The reports converge on three things, and those drive most of this module:

- **The null and the control do not share a geometry.** `null_distribution`
  builds both halves of its null as compact patches of stable ground, while the
  real reference (zone C) is 398 scattered pixels with a measured N_eff of 5. If
  the fragmented grassland carries more variance than a compact patch of equal
  pixel count, the null is too quiet and the empirical *p* is anti-conservative.
  `null_against_real_reference` re-runs the null keeping the real C.
- **The lake control may not be independent.** Zone B is 65 pixels wholly
  enclosed by the mat; at ~40 m spacing, border pixels necessarily mix the two.
  If B is contaminated by A then "B oscillates like A" is tautological.
  `erode_zone` lets the test be repeated on an eroded lake.
- **Aggregation of a decorrelated stack is biased toward zero.** A recovered
  amplitude is a lower bound on the true one, so a small measured amplitude
  cannot be read as an absence of motion. `attenuation_bias` measures the
  recovery ratio at the site's own coherence.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

from .aggregate import (_compact_blob, aggregate_unwrapped, empirical_pvalue,
                        invert_aggregate, seasonal_amplitude)
from .inversion.isbas import PHASE_TO_MM


# --------------------------------------------------------------- mask surgery
def erode_zone(mask: xr.DataArray, iterations: int = 1) -> xr.DataArray:
    """Shrink a zone by `iterations` pixels of border.

    Answers the contamination objection to the lake control: zone B is 65 pixels
    entirely enclosed by the mat, so at ~40 m spacing its border pixels contain
    mat vegetation through the resolution cell and sidelobes. If the B − C
    seasonal amplitude survives erosion, the control is independent; if it
    collapses toward the mat's value, "the lake oscillates too" was partly a
    measurement of the mat."""
    from scipy.ndimage import binary_erosion
    out = binary_erosion(mask.values, iterations=iterations,
                         border_value=0)
    return xr.DataArray(out, coords=mask.coords, dims=mask.dims)


def subdivide_zone(mask: xr.DataArray, n_px: int, seed: int = 0,
                   compact: bool = True) -> xr.DataArray | None:
    """A sub-patch of `n_px` pixels drawn from `mask`.

    Answers the block-deformation objection. Aggregation assumes the mat behaves
    as one hydrological unit; that assumption is load-bearing for the signal it
    produced, so it must be tested now rather than deferred. Recomputing the
    amplitude at 499 -> 250 -> 125 -> 60 pixels separates two cases: an amplitude
    that holds while the null floor rises as 1/sqrt(N) means the signal is
    spatially coherent across the mat; an amplitude that collapses faster than
    the null means the aggregate is not a single unit.

    `compact=True` draws a contiguous blob, which is the geometry the null uses;
    `compact=False` draws scattered pixels, which is the geometry zone C has."""
    yx = np.argwhere(mask.values)
    if len(yx) < n_px:
        return None
    rng = np.random.default_rng(seed)
    if compact:
        idx = _compact_blob(yx, int(rng.integers(len(yx))), n_px)
    else:
        idx = rng.choice(len(yx), n_px, replace=False)
    a = np.zeros_like(mask.values)
    a[yx[idx, 0], yx[idx, 1]] = True
    return xr.DataArray(a, coords=mask.coords, dims=mask.dims)


def amplitude_vs_size(unw: xr.DataArray, corr: xr.DataArray, zones: dict,
                      sizes=(499, 250, 125, 60), reference: str = "C",
                      target: str = "A", n_draws: int = 12,
                      seed: int = 0, subdivide: str = "target") -> pd.DataFrame:
    """Seasonal amplitude of `target` - `reference` as one of them shrinks.

    Returns one row per (size, draw). Aggregate noise falls as 1/sqrt(N), so the
    expectation under "no coherent signal" is that amplitude RISES as the patch
    shrinks, tracking the noise floor. A flat amplitude is the signature of a
    real, spatially coherent signal.

    `subdivide` selects WHICH side shrinks, and running both is what makes the
    test conclusive. Subdividing only the target is blind to the one alternative
    that still matters: a seasonal term carried by the REFERENCE zone, or by any
    regionally common signal, is equally flat under subdivision of the target,
    because every sub-patch keeps the same reference. Subdividing the reference
    instead, at fixed target, separates the two."""
    if subdivide not in ("target", "reference"):
        raise ValueError("subdivide must be 'target' or 'reference'")
    shrink = target if subdivide == "target" else reference
    rows = []
    for n in sizes:
        for d in range(n_draws):
            sub = subdivide_zone(zones[shrink], n, seed=seed + 1000 * d)
            if sub is None:
                continue
            z = {**zones, shrink: sub}
            try:
                dd = aggregate_unwrapped(unw, corr, z, target, reference)
                if len(dd) < 10:
                    continue
                s = seasonal_amplitude(invert_aggregate(dd))
                rows.append({"n_px": n, "draw": d, "subdivided": shrink,
                             "amplitude_mm": s["amplitude_mm"],
                             "phase_doy": s.get("phase_doy"),
                             "r2_seasonal": s["r2_seasonal"]})
            except Exception:
                continue
    return pd.DataFrame(rows, columns=["n_px", "draw", "subdivided",
                                       "amplitude_mm", "phase_doy",
                                       "r2_seasonal"])


# ------------------------------------------------------------ null, revisited
def null_against_real_reference(unw: xr.DataArray, corr: xr.DataArray,
                                zones: dict, template: xr.DataArray,
                                n_target: int, reference: str = "C",
                                ref_pool: str = "D", n_trials: int = 100,
                                seed: int = 0) -> pd.DataFrame:
    """Null distribution that keeps the REAL reference zone.

    The published null draws BOTH halves from stable ground, so it never
    contains the real zone C. That matters because C is fragmented (398
    scattered pixels, measured N_eff = 5): if most of the variance of the A − C
    series comes from the grassland, a null built without it is too quiet and
    the empirical *p* is anti-conservative.

    Here only the target patch is randomised, drawn from `ref_pool` and matched
    in pixel count to the real target, while `reference` stays the actual zone.
    Comparing this distribution with the published one answers the question
    directly: if they agree, the objection is answered; if this one is wider,
    the reference term dominates and the *p*-values need restating."""
    pool = np.argwhere(zones[ref_pool].values)
    if len(pool) < n_target:
        return pd.DataFrame(columns=["trial", "amplitude_mm", "r2_seasonal"])
    rows = []
    rng = np.random.default_rng(seed)
    for t in range(n_trials):
        idx = _compact_blob(pool, int(rng.integers(len(pool))), n_target)
        a = np.zeros_like(zones[ref_pool].values)
        a[pool[idx, 0], pool[idx, 1]] = True
        z = {**zones, "_null": xr.DataArray(a, coords=template.coords,
                                            dims=template.dims)}
        try:
            dd = aggregate_unwrapped(unw, corr, z, "_null", reference)
            if len(dd) < 10:
                continue
            s = seasonal_amplitude(invert_aggregate(dd))
            rows.append({"trial": t, "amplitude_mm": s["amplitude_mm"],
                         "r2_seasonal": s["r2_seasonal"]})
        except Exception:
            continue
    return pd.DataFrame(rows, columns=["trial", "amplitude_mm", "r2_seasonal"])


def multi_control_patches(unw: xr.DataArray, corr: xr.DataArray, zones: dict,
                          template: xr.DataArray, n_px: int,
                          target: str = "A", pool: str = "D",
                          n_controls: int = 8, seed: int = 0) -> pd.DataFrame:
    """Repeat the headline test against several independent compact controls.

    Answers two objections at once: that zone C is a single control, and that
    its fragmented geometry makes its N_eff uninterpretable. Each control here
    is compact and contiguous, so its N_eff is meaningful and its construction
    matches the null's. A seasonal amplitude stable across controls is far
    stronger evidence than one measured against a single fragmented patch."""
    cand = np.argwhere(zones[pool].values)
    rows = []
    rng = np.random.default_rng(seed)
    for c in range(n_controls):
        if len(cand) < n_px:
            break
        idx = _compact_blob(cand, int(rng.integers(len(cand))), n_px)
        a = np.zeros_like(zones[pool].values)
        a[cand[idx, 0], cand[idx, 1]] = True
        z = {**zones, "_ctrl": xr.DataArray(a, coords=template.coords,
                                            dims=template.dims)}
        try:
            dd = aggregate_unwrapped(unw, corr, z, target, "_ctrl")
            if len(dd) < 10:
                continue
            s = seasonal_amplitude(invert_aggregate(dd))
            rows.append({"control": c, "n_px": n_px,
                         "amplitude_mm": s["amplitude_mm"],
                         "phase_doy": s.get("phase_doy"),
                         "r2_seasonal": s["r2_seasonal"]})
        except Exception:
            continue
    return pd.DataFrame(rows, columns=["control", "n_px", "amplitude_mm",
                                       "phase_doy", "r2_seasonal"])


# ------------------------------------------------------------- uncertainty
def bootstrap_amplitude_ci(dd: pd.DataFrame, n_boot: int = 2000,
                           seed: int = 0, level: float = 95.0) -> dict:
    """Confidence interval on the seasonal amplitude, resampling ACQUISITIONS.

    The headline amplitude is reported with an empirical *p*-value but no
    uncertainty. Resampling is over acquisition dates rather than over pairs,
    for the same reason the coherence test uses a date-jackknife: the ~356 pairs
    share ~90 dates, so resampling pairs would treat correlated observations as
    independent and produce an interval that is too narrow."""
    if "pair" not in dd.columns or dd.empty:
        return {"n_boot": 0}
    dates = sorted({d for p in dd["pair"] for d in (p[:8], p[9:17])})
    rng = np.random.default_rng(seed)
    amps = []
    for _ in range(n_boot):
        keep = set(rng.choice(dates, len(dates), replace=True))
        sub = dd[dd["pair"].map(lambda p: p[:8] in keep and p[9:17] in keep)]
        if len(sub) < 10:
            continue
        try:
            amps.append(seasonal_amplitude(invert_aggregate(sub))["amplitude_mm"])
        except Exception:
            continue
    if not amps:
        return {"n_boot": 0}
    lo, hi = np.percentile(amps, [(100 - level) / 2, 100 - (100 - level) / 2])
    return {"n_boot": len(amps), "mean": float(np.mean(amps)),
            "median": float(np.median(amps)), "ci_low": float(lo),
            "ci_high": float(hi), "level": level,
            "std": float(np.std(amps, ddof=1))}


def detectable_closure_bias(se: float, n_triplets: int,
                            sigma: float = 2.0) -> dict:
    """Smallest closure bias the network could have resolved.

    A dielectric mechanism is expected to bias closure phase; none was detected.
    That non-detection only argues against the mechanism if the test had the
    power to see the expected bias, so the detection limit must be stated
    alongside it. `se` is the standard error already reported per zone."""
    limit = sigma * se
    return {"n_triplets": int(n_triplets), "se_rad": float(se),
            "sigma": sigma, "detectable_bias_rad": float(limit),
            "detectable_bias_mm": float(abs(limit * PHASE_TO_MM))}


# --------------------------------------------------------------- simulation
def noise_floor_simulation(n_dates: int = 90, n_pairs: int = 356,
                           n_trials: int = 200, seed: int = 0,
                           coherence: float = 0.4,
                           idx: np.ndarray | None = None) -> dict:
    """Temporal coherence returned by a FULLY DECORRELATED pixel.

    The 0.55 floor carries the whole H1 reading — A at 0.604 as near-noise, B at
    0.584 as at-floor, C at 0.734 as signal — but was stated as a single number
    with no design, no spread and no sensitivity.

    The floor is high (~0.55, not ~0.05) because temporal coherence measures the
    agreement between a phase history that was ESTIMATED FROM THE SAME DATA and
    the interferograms it was fitted to. With `n_dates - 1` free parameters
    fitted to `n_pairs` observations the estimator explains part of pure noise,
    and the floor is set by that redundancy. Simulating it therefore requires
    running the real EVD, not just averaging random phasors — which is why this
    reproduces the whole per-pixel inversion on random input.

    Pass `idx` (the real ``(i, j)`` date-index pairs — the SECOND element
    returned by ``_pair_date_index(pairs)``, not the first, which is the date
    array itself) whenever it is available. The floor is strongly
    topology-dependent — roughly 0.69 at redundancy 2 against 0.36 at redundancy
    8 — so a floor quoted without its network is not reproducible. The synthetic
    fallback below only fixes the redundancy, not the baseline distribution."""
    from .inversion.phaselinking import evd_pixel

    rng = np.random.default_rng(seed)
    if idx is None:
        # fallback: a network with the requested redundancy only
        built = []
        span = 1
        while len(built) < n_pairs and span < n_dates:
            for i in range(n_dates - span):
                built.append((i, i + span))
                if len(built) >= n_pairs:
                    break
            span += 1
        idx = np.array(built[:n_pairs])
    else:
        idx = np.asarray(idx)
        # `idx` must hold small non-negative date INDICES, not dates or any
        # other large-valued array. The most likely way to get this wrong is
        # unpacking `_pair_date_index(pairs)` in the wrong order — its first
        # return value is the date array, not the index pairs — which yields
        # nanosecond-epoch integers here and would otherwise try to allocate a
        # matrix with tens of trillions of entries. Fail with a clear message
        # instead of letting that allocation raise deep inside numpy.
        if idx.dtype.kind not in "iu" or idx.min() < 0 or idx.max() > 100_000:
            raise ValueError(
                "idx must be small non-negative integer (i, j) date indices, "
                f"got dtype={idx.dtype}, range=[{idx.min()}, {idx.max()}]. "
                "Did you unpack _pair_date_index(pairs) as `idx, _ = ...` "
                "instead of `_, idx = ...`? It returns (dates, idx).")
        n_dates = max(n_dates, int(idx.max()) + 1)

    vals = []
    for _ in range(n_trials):
        phi = rng.uniform(-np.pi, np.pi, len(idx))       # fully decorrelated
        coh = np.full(len(idx), coherence)
        _, tc = evd_pixel(phi, coh, idx, n_dates)
        if np.isfinite(tc):
            vals.append(float(tc))
    v = np.asarray(vals)
    if not v.size:
        return {"n_dates": n_dates, "n_pairs": n_pairs, "n_trials": 0}
    return {"n_dates": n_dates, "n_pairs": len(idx), "n_trials": len(v),
            "redundancy": len(idx) / max(n_dates - 1, 1),
            "median": float(np.median(v)), "mean": float(v.mean()),
            "p05": float(np.percentile(v, 5)),
            "p95": float(np.percentile(v, 95)),
            "std": float(v.std(ddof=1)), "values": v}


def estimator_noise_bias(true_amp_mm: float = 3.3, coherence: float = 0.4,
                         n_eff: int = 31, n_dates: int = 90,
                         n_trials: int = 400, seed: int = 0) -> pd.DataFrame:
    """Bias of the seasonal-amplitude estimator under phase noise.

    Note the SIGN, which is the opposite of what "attenuation" suggests.
    Amplitude is ``sqrt(a^2 + b^2)``, a strictly positive function of two noisy
    coefficients, so noise INFLATES it (a Rice bias): a measured 3.3 mm on a
    noisy series corresponds to a slightly smaller truth, not a larger one.

    This is a distinct effect from the loss of signal described by
    `coupling_scaled_bound`, and the two push in opposite directions. Reporting
    only one of them would misstate the uncertainty, so both are measured."""
    rng = np.random.default_rng(seed)
    dates = pd.date_range("2022-01-01", periods=n_dates, freq="12D")
    t = (dates - dates[0]).days.values / 365.25
    truth = true_amp_mm * np.cos(2 * np.pi * t)
    sigma_phi = np.sqrt((1 - coherence ** 2) / (2 * coherence ** 2))
    sigma_mm = abs(sigma_phi * PHASE_TO_MM) / np.sqrt(max(n_eff, 1))
    rows = []
    for k in range(n_trials):
        s = seasonal_amplitude(pd.DataFrame(
            {"date": dates, "disp_mm": truth + rng.normal(0, sigma_mm, n_dates)}))
        rows.append({"trial": k, "recovered_mm": s["amplitude_mm"],
                     "true_mm": true_amp_mm,
                     "ratio": s["amplitude_mm"] / true_amp_mm})
    df = pd.DataFrame(rows)
    df.attrs["sigma_mm_aggregate"] = float(sigma_mm)
    df.attrs["sigma_phi_per_pixel_rad"] = float(sigma_phi)
    return df


def coupling_scaled_bound(observed_bound_mm: float = 3.9,
                          couplings=(1.0, 0.75, 0.5, 0.25, 0.1)) -> pd.DataFrame:
    """What the measured bound implies about MAT motion, per coupling fraction.

    This is the manuscript's central internal tension made quantitative. The
    paper concludes that the phase centre sits in a saturated canopy volume
    decoupled from the substrate and tracks moisture rather than the peat. If
    that is right, the phase does not observe the mat directly: writing `f` for
    the fraction of the phase that responds to mat motion, the observable is
    ``f x true_motion``, so a bound of `observed_bound_mm` on the observable
    implies a bound of ``observed_bound_mm / f`` on the mat.

    At f = 1 (phase centre rigidly attached to the surface) the published bound
    stands. As f falls the implied bound diverges, and f is not constrained by
    these data. The honest statement is therefore a bound on APPARENT
    phase-centre displacement, plus an explicit coupling assumption — which is
    what this table is for."""
    rows = [{"coupling_f": float(f),
             "implied_mat_bound_mm": float(observed_bound_mm / f) if f > 0
             else float("inf"),
             "interpretation": ("phase centre rigid with the surface"
                                if f >= 0.999 else
                                f"only {100 * f:.0f} % of the phase follows the mat")}
            for f in couplings]
    df = pd.DataFrame(rows)
    df.attrs["observed_bound_mm"] = float(observed_bound_mm)
    return df


def synthetic_seasonal_recovery(true_amp_mm: float = 3.3,
                                coherence: float = 0.4,
                                n_px_aggregate: int = 499,
                                n_eff: int = 31, n_dates: int = 90,
                                n_trials: int = 200, seed: int = 0) -> dict:
    """Per-pixel versus aggregated recovery of a SEASONAL cycle.

    The published end-to-end validation recovers a VELOCITY, which the same
    manuscript argues has no power on a periodic signal. So the observable
    carrying the headline result had no synthetic validation at all.

    What the comparison shows is DISPERSION, not a shift in the median: fitting
    an annual cycle to 90 dates averages noise down even per pixel, so the
    median per-pixel estimate is not far off. The difference is that any single
    per-pixel estimate is unusable — its interquartile range is several times
    the signal — while the aggregate is tight. Reporting this as "per-pixel
    fails to recover the signal" would overstate it; the accurate claim is that
    per-pixel cannot resolve it."""
    rng = np.random.default_rng(seed)
    dates = pd.date_range("2022-01-01", periods=n_dates, freq="12D")
    t = (dates - dates[0]).days.values / 365.25
    truth = true_amp_mm * np.cos(2 * np.pi * t)
    sigma_phi = np.sqrt((1 - coherence ** 2) / (2 * coherence ** 2))
    sigma_px = abs(sigma_phi * PHASE_TO_MM)
    per_px, agg = [], []
    for _ in range(n_trials):
        per_px.append(seasonal_amplitude(pd.DataFrame(
            {"date": dates, "disp_mm": truth + rng.normal(0, sigma_px, n_dates)}
        ))["amplitude_mm"])
        agg.append(seasonal_amplitude(pd.DataFrame(
            {"date": dates,
             "disp_mm": truth + rng.normal(0, sigma_px / np.sqrt(n_eff), n_dates)}
        ))["amplitude_mm"])
    per_px, agg = np.asarray(per_px), np.asarray(agg)

    def iqr(a):
        return float(np.percentile(a, 75) - np.percentile(a, 25))

    return {"true_mm": true_amp_mm, "coherence": coherence, "n_eff": n_eff,
            "sigma_px_mm": float(sigma_px),
            "per_pixel_median": float(np.median(per_px)),
            "per_pixel_iqr": iqr(per_px),
            "per_pixel_p95": float(np.percentile(per_px, 95)),
            "aggregate_median": float(np.median(agg)),
            "aggregate_iqr": iqr(agg),
            "aggregate_p05": float(np.percentile(agg, 5)),
            "aggregate_p95": float(np.percentile(agg, 95)),
            "iqr_ratio": iqr(per_px) / iqr(agg) if iqr(agg) else float("nan"),
            "per_pixel": per_px, "aggregate": agg}


# ------------------------------------------------------------- multiplicity
def family_corrected_pvalue(observed_best: float, null_family: np.ndarray
                            ) -> dict:
    """Empirical *p* for the BEST forcing at its best lag.

    The lag sweep is already handled by replaying it on the null. The family of
    forcings is not: seven forcings times sixteen lags is a large selection
    space, and the reported result is the maximum over it. The null must
    therefore undergo the same maximisation — the same rule, applied one level
    up. `null_family` holds, per realisation, the best |r| over the whole
    family."""
    return empirical_pvalue(abs(observed_best), np.abs(np.asarray(null_family)))


def matrix_fill(n_dates: int, n_pairs: int) -> dict:
    """How much of the N x N coherence matrix the network actually populates.

    Phase linking is described as the maximum-likelihood estimator over the full
    coherence matrix. With 356 pairs across ~90 dates the matrix is ~9 % filled,
    which bounds how much that description can carry."""
    possible = n_dates * (n_dates - 1) // 2
    return {"n_dates": int(n_dates), "n_pairs": int(n_pairs),
            "possible_pairs": int(possible),
            "fill_fraction": float(n_pairs / possible) if possible else float("nan"),
            "redundancy": float(n_pairs / max(n_dates - 1, 1))}


# ------------------------------------------- second-round referee additions
def count_closed_triplets(pairs: list[str]) -> dict:
    """Closed triangles in the interferogram network, two independent ways.

    The count is a deterministic property of the pair list, so leaving it
    "unresolved" is indefensible: trace(A^3)/6 on the adjacency matrix and a
    direct enumeration must agree. Reported alongside the count a network of the
    same density would give if its pairs were drawn at random, because HyP3
    pairs follow constrained baselines and therefore close far fewer triangles
    than an unstructured graph of equal density."""
    dates = sorted({d for p in pairs for d in (str(p)[:8], str(p)[9:17])})
    di = {d: i for i, d in enumerate(dates)}
    n = len(dates)
    A = np.zeros((n, n), dtype=np.int64)
    for p in pairs:
        i, j = di[str(p)[:8]], di[str(p)[9:17]]
        A[i, j] = A[j, i] = 1
    trace_count = int(np.trace(A @ A @ A) // 6)

    enumerated = 0
    for i in range(n):
        for j in range(i + 1, n):
            if not A[i, j]:
                continue
            enumerated += int(np.sum(A[j, j + 1:] & A[i, j + 1:]))

    m = int(A.sum() // 2)
    p_edge = 2 * m / (n * (n - 1)) if n > 1 else 0.0
    random_expect = p_edge ** 3 * n * (n - 1) * (n - 2) / 6
    return {"n_dates": n, "n_pairs": m, "triplets_trace": trace_count,
            "triplets_enumerated": enumerated,
            "agree": trace_count == enumerated,
            "random_graph_expectation": float(random_expect),
            "ratio_to_random": float(trace_count / random_expect)
            if random_expect else float("nan")}


def excess_above_floor(values: dict, floor: float) -> pd.DataFrame:
    """Coherence expressed as excess over the network noise floor.

    Replaces a threshold-based dichotomy with a statement that needs no
    threshold at all: what fraction of the reference zone's excess above the
    floor does the target retain? Robust to the floor moving, which it does —
    it is topology-dependent."""
    rows = [{"zone": z, "value": float(v), "excess": float(v) - floor}
            for z, v in values.items()]
    df = pd.DataFrame(rows)
    ref = df["excess"].max()
    df["frac_of_max_excess"] = df["excess"] / ref if ref else np.nan
    df.attrs["floor"] = float(floor)
    return df


def toroidal_permutation_test(subset: xr.DataArray, zone: xr.DataArray,
                              field: xr.DataArray, n_trials: int = 2000,
                              seed: int = 0, stat=np.median) -> dict:
    """Test a spatial subset's statistic while PRESERVING its clustering.

    A rank test on clustered pixels treats neighbours as independent
    observations, which they are not at a correlation length of several pixels;
    it therefore returns a *p*-value that is far too small. This shifts the
    subset rigidly and toroidally inside the zone instead, so every null
    realisation has the same size, the same shape and the same internal
    autocorrelation as the observed one, and only its position changes.

    Returns the observed statistic, the null distribution, and an empirical
    *p*-value that is two-sided by default in the sense of |deviation|."""
    sm, zm = subset.values, zone.values
    obs_vals = field.values[sm & np.isfinite(field.values)]
    if not obs_vals.size:
        return {"n_subset": 0}
    observed = float(stat(obs_vals))

    yx = np.argwhere(sm)
    y0, x0 = yx.min(0)
    rel = yx - (y0, x0)                     # shape preserved exactly
    ny, nx = sm.shape
    rng = np.random.default_rng(seed)
    nulls = []
    for _ in range(n_trials):
        dy, dx = rng.integers(0, ny), rng.integers(0, nx)
        yy = (rel[:, 0] + dy) % ny
        xx = (rel[:, 1] + dx) % nx
        if not zm[yy, xx].all():            # must land wholly inside the zone
            continue
        v = field.values[yy, xx]
        v = v[np.isfinite(v)]
        if v.size:
            nulls.append(float(stat(v)))
    nulls = np.asarray(nulls)
    if not nulls.size:
        return {"n_subset": int(sm.sum()), "observed": observed, "n_null": 0,
                "note": "no shift placed the subset wholly inside the zone"}
    zone_vals = field.values[zm & np.isfinite(field.values)]
    k = int(np.sum(np.abs(nulls - np.median(nulls))
                   >= abs(observed - np.median(nulls))))
    return {"n_subset": int(sm.sum()), "observed": observed,
            "zone_statistic": float(stat(zone_vals[np.isfinite(zone_vals)])),
            "n_null": int(nulls.size), "null_median": float(np.median(nulls)),
            "null_p05": float(np.percentile(nulls, 5)),
            "null_p95": float(np.percentile(nulls, 95)),
            "p_value": float((1 + k) / (1 + nulls.size)),
            "p_floor": float(1 / (1 + nulls.size)), "nulls": nulls}


def wrapped_seasonal_amplitude(dd: pd.DataFrame, date_col: str = "pair",
                               value_col: str = "ddphase_rad") -> dict:
    """Fit the annual cycle directly on per-pair WRAPPED phase differences.

    Closes a hole opened by conceding that immunity to unwrapping errors does
    not survive the network inversion to date-referenced values — because the
    published amplitude comes from exactly that inversion.

    Each pair contributes ``phi_j - phi_i``, so a cycle ``a·cos + b·sin`` on the
    dates predicts the pair difference directly and no inversion is needed. Safe
    here only because the amplitude is small: 3.3 mm is ~0.75 rad two-way, well
    inside ±pi, so wrapping does not alias. The check would be invalid for a
    signal approaching half a wavelength."""
    from .inversion.isbas import PHASE_TO_MM

    d = dd.dropna(subset=[value_col]).copy()
    ref = pd.to_datetime([str(p)[:8] for p in d[date_col]])
    sec = pd.to_datetime([str(p)[9:17] for p in d[date_col]])
    t0 = min(ref.min(), sec.min())
    tr = (ref - t0).days.values / 365.25
    ts = (sec - t0).days.values / 365.25
    # design: difference of the harmonic between the two dates, plus a trend
    X = np.column_stack([
        ts - tr,
        np.cos(2 * np.pi * ts) - np.cos(2 * np.pi * tr),
        np.sin(2 * np.pi * ts) - np.sin(2 * np.pi * tr),
    ])
    y = d[value_col].values * PHASE_TO_MM        # rad -> mm, no inversion
    w = d["weight"].values if "weight" in d else np.ones(len(d))
    W = np.sqrt(w)[:, None]
    coef, *_ = np.linalg.lstsq(X * W, y * np.sqrt(w), rcond=None)
    trend, a, b = coef
    resid = y - X @ coef
    ss = float(np.sum((y - y.mean()) ** 2))
    return {"amplitude_mm": float(np.hypot(a, b)),
            "phase_doy": float((np.degrees(np.arctan2(b, a)) % 360) * 365.25 / 360),
            "trend_mm_yr": float(trend), "n_pairs": int(len(d)),
            "r2": float(1 - np.sum(resid ** 2) / ss) if ss else np.nan,
            "rms_residual_mm": float(np.sqrt(np.mean(resid ** 2)))}


def matched_cover_pool(zones: dict, worldcover: xr.DataArray,
                       dominant_class: int | None = None,
                       exclude_reference: str | None = "C") -> xr.DataArray:
    """Pixels of the SAME land-cover class as the mat, for independent controls.

    The multi-control test drew its patches from zone D, and that was wrong in a
    way the zone definitions make structural: D is built as
    ``outside & ~water & flat & ~C``, i.e. as the COMPLEMENT of the matched
    reference. Controls drawn from it are guaranteed *not* to be land-cover
    matched, so their spread measures how much the cover varies, not how much
    the mat's amplitude depends on the choice of a comparable control.

    The right pool is the same-class terrain: ``(C | D)`` recovers
    ``outside & ~water & flat`` exactly, and intersecting it with the mat's
    dominant WorldCover class gives matched candidates. The published reference
    is excluded by default so the controls are independent of it.

    This is the pool for the test that decides whether the seasonal result is a
    property of the mat or of the reference chosen to measure it against."""
    if dominant_class is None:
        from .stratify import dominant_class as _dom
        dominant_class = _dom(worldcover, zones["A"])
    eligible = zones["C"] | zones["D"]            # outside & ~water & flat
    pool = eligible & (worldcover == dominant_class)
    if exclude_reference and exclude_reference in zones:
        pool = pool & ~zones[exclude_reference]
    pool = pool.astype(bool)
    pool.attrs["dominant_class"] = int(dominant_class)
    pool.attrs["n_px"] = int(pool.values.sum())
    return pool


def matched_null_pairs(unw: xr.DataArray, corr: xr.DataArray, zones: dict,
                       template: xr.DataArray, pool: str, n_target: int,
                       n_reference: int, n_trials: int = 300,
                       seed: int = 0) -> pd.DataFrame:
    """Null from two INDEPENDENT compact patches of the same matched pool.

    `null_distribution` cuts a single contiguous blob of ``n_target +
    n_reference`` pixels into two adjacent halves. That has two problems here.

    It is **geometrically unlike the real observable**: the mat and its
    reference lie about a kilometre apart, not side by side, so a null built
    from adjacent halves shares far more of the atmospheric screen than the
    quantity it is supposed to mimic.

    And it is **infeasible on a matched pool**. Requiring 897 pixels of the
    mat's own land-cover class in one contiguous patch outside the site is a
    much stronger demand than drawing two patches of 499 and 398 anywhere in
    that class; when it cannot be met, `null_distribution` returns an empty
    frame and every statistic computed from it fails.

    Drawing the two patches independently, disjointly, and from the same
    matched pool fixes both. If even the disjoint requirement cannot be met the
    two sizes are scaled down by a common factor rather than failing: fewer
    pixels means more aggregate noise, so the null becomes WIDER and the test
    more conservative, which is the safe direction to err in. The factor is
    recorded in ``.attrs['size_scale']`` so it is reported rather than hidden."""
    cand = np.argwhere(zones[pool].values)
    scale = 1.0
    nt, nr = int(n_target), int(n_reference)
    if len(cand) < nt + nr:
        scale = len(cand) / (nt + nr)
        nt, nr = max(int(nt * scale), 20), max(int(nr * scale), 20)
    if len(cand) < nt + nr:
        out = pd.DataFrame(columns=["trial", "amplitude_mm", "r2_seasonal"])
        out.attrs.update(size_scale=0.0, n_target=0, n_reference=0,
                         pool_px=int(len(cand)),
                         note="pool too small even after scaling")
        return out

    rng = np.random.default_rng(seed)
    rows = []
    for t in range(n_trials):
        first = _compact_blob(cand, int(rng.integers(len(cand))), nt)
        keep = np.ones(len(cand), bool)
        keep[first] = False                      # disjoint: no shared pixels
        rest = cand[keep]
        if len(rest) < nr:
            continue
        second = _compact_blob(rest, int(rng.integers(len(rest))), nr)
        z = dict(zones)
        for name, pts in (("_n1", cand[first]), ("_n2", rest[second])):
            a = np.zeros_like(zones[pool].values)
            a[pts[:, 0], pts[:, 1]] = True
            z[name] = xr.DataArray(a, coords=template.coords, dims=template.dims)
        try:
            dd = aggregate_unwrapped(unw, corr, z, "_n1", "_n2")
            if len(dd) < 10:
                continue
            s = seasonal_amplitude(invert_aggregate(dd))
            rows.append({"trial": t, "amplitude_mm": s["amplitude_mm"],
                         "r2_seasonal": s["r2_seasonal"]})
        except Exception:
            continue
    out = pd.DataFrame(rows, columns=["trial", "amplitude_mm", "r2_seasonal"])
    out.attrs.update(size_scale=float(scale), n_target=nt, n_reference=nr,
                     pool_px=int(len(cand)),
                     note="" if scale == 1.0 else
                     f"patch sizes scaled by {scale:.2f} to fit the pool; "
                     "the null is therefore conservative (wider)")
    return out
