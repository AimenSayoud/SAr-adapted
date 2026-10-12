"""One scorer for every test against a surface truth (C-061, the InSAR full test, hub
``08_deliverables/insar_full_test/PLAN.md`` §2).

Every test is scored here, with the same views, so that no test invents its own statistics:

- **increments** — per pair: Pearson r, RMSE, slope through the origin, share within ±λ/4 (the ambiguity bound);
- **levels** — the estimate chained along consecutive pairs (a missing pair starts a new segment), each segment
  re-zeroed at its first scored epoch, then compared with the truth after one linear trend is removed (Pearson r,
  RMSE, shares within ±5 / ±10 mm, amplitude ratio — the Hrysiewicz et al. 2024 statistics, X-082):
  over the full span, per truth datum, and per season of each year (re-zeroed at the season start);
- **winters** — the level step across each winter gap (last epoch of the season → first of the next spring) against
  the truth's step: the across-winter loss of the level (SETTLED §4);
- **seasons** — increments by meteorological season;
- **years** — increments per year with the sign and slope of the relation, to see whether the relation is stable
  between years (Mleczko et al. 2025 found its sign can flip).

Nothing site-specific lives in this module (the repository is public): truth datum breaks, season windows and the
lockbox years come from a ``Protocol``, which the hub keeps in ``08_deliverables/insar_full_test/protocol.yaml``.
The **lockbox** guard refuses to score lockbox years unless the call says ``lockbox=True`` — that call is made once,
by the frozen method (X-107).

Pre-registration: ``assert_preregistered(prereg_path, outputs)`` refuses a run whose pre-registration file is
missing, invalid, uncommitted, or committed after the run's first output.
"""
from __future__ import annotations

import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from .cube.core import WAVELENGTH_MM

QUARTER_WAVE_MM = WAVELENGTH_MM / 4.0
SEASON_OF_MONTH = {12: "winter", 1: "winter", 2: "winter", 3: "spring", 4: "spring", 5: "spring",
                   6: "summer", 7: "summer", 8: "summer", 9: "autumn", 10: "autumn", 11: "autumn"}


class LockboxError(RuntimeError):
    """A lockbox year was about to be scored outside the single lockbox run."""


@dataclass
class Protocol:
    """The site-specific part of the scoring, kept outside this public module."""
    datum_breaks: list[tuple[pd.Timestamp, pd.Timestamp]] = field(default_factory=list)  # truth gaps between datums
    season_months: tuple[int, int] = (5, 11)            # the snow-free season scored "within season"
    spring_window: tuple[int, int] = (3, 6)             # months in which the first post-winter epoch is looked for
    min_epochs: int = 8
    lockbox_years: tuple[int, ...] = ()

    @classmethod
    def load(cls, path: str | Path) -> "Protocol":
        d = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
        breaks = [(pd.Timestamp(a, tz="UTC"), pd.Timestamp(b, tz="UTC")) for a, b in d.get("datum_breaks", [])]
        return cls(datum_breaks=breaks, season_months=tuple(d.get("season_months", (5, 11))),
                   spring_window=tuple(d.get("spring_window", (3, 6))), min_epochs=int(d.get("min_epochs", 8)),
                   lockbox_years=tuple(int(y) for y in d.get("lockbox_years", [])))

    def datum(self, times: pd.Series) -> np.ndarray:
        """Datum label per time: 'A' before the first break, 'B' after it, 'C' after the second…"""
        t = pd.to_datetime(times, utc=True)
        lab = np.zeros(len(t), dtype=int)
        for a, b in self.datum_breaks:
            lab += (t >= b).to_numpy().astype(int)
        return np.array([chr(ord("A") + k) for k in lab])


# --------------------------------------------------------------------------------------------- guards
def guard_lockbox(times, protocol: Protocol, lockbox: bool = False) -> None:
    if lockbox or not protocol.lockbox_years:
        return
    years = set(pd.to_datetime(pd.Series(times), utc=True).dt.year.unique())
    hit = sorted(years & set(protocol.lockbox_years))
    if hit:
        raise LockboxError(f"lockbox years {hit} in the data: drop them, or score once with lockbox=True (X-107)")


def _git_commit_time(path: Path) -> float | None:
    try:
        r = subprocess.run(["git", "-C", str(path.parent), "log", "-1", "--format=%ct", "--", path.name],
                           capture_output=True, text=True, timeout=20)
        dirty = subprocess.run(["git", "-C", str(path.parent), "status", "--porcelain", "--", path.name],
                               capture_output=True, text=True, timeout=20).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None
    return float(r.stdout.strip()) if r.stdout.strip() and not dirty else None


PREREG_REQUIRED = ("ticket", "question", "primary", "comparison", "years", "better_if")


def load_preregistration(path: str | Path) -> dict:
    d = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    missing = [k for k in PREREG_REQUIRED if not d.get(k)]
    if missing:
        raise ValueError(f"{path}: pre-registration lacks {missing}")
    return d


def assert_preregistered(path: str | Path, outputs: list[str | Path] = ()) -> dict:
    """The pre-registration exists, is complete, is committed (clean), and was committed before any output exists."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"no pre-registration at {path}: write and commit it before the run")
    d = load_preregistration(path)
    t = _git_commit_time(path)
    if t is None:
        raise RuntimeError(f"{path} is not committed (or has uncommitted changes): commit it before the run")
    early = [Path(o) for o in outputs if Path(o).exists() and Path(o).stat().st_mtime < t]
    if early:
        raise RuntimeError(f"outputs older than the pre-registration commit: {[str(o) for o in early[:3]]}")
    return d


# ----------------------------------------------------------------------------------------- statistics
def _r(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.corrcoef(a, b)[0, 1]) if len(a) > 2 and np.std(a) > 0 and np.std(b) > 0 else float("nan")


def detrended_stats(time, truth, est) -> dict:
    """Hrysiewicz et al. (2024) / X-082: one linear trend removed from each, then r, RMSE, shares, amplitude ratio."""
    t = pd.to_datetime(pd.Series(time), utc=True)
    tt = ((t - t.iloc[0]).dt.total_seconds() / 86400).to_numpy()
    tr, es = np.asarray(truth, float), np.asarray(est, float)
    dl = tr - np.polyval(np.polyfit(tt, tr, 1), tt)
    ds = es - np.polyval(np.polyfit(tt, es, 1), tt)
    dev = ds - dl
    return {"epochs": len(tr), "pearson_detrended": _r(dl, ds), "rmse_detrended_mm": float(np.sqrt(np.mean(dev ** 2))),
            "laser_sd_mm": float(np.std(dl)), "within_5mm": float(np.mean(np.abs(dev) < 5)),
            "within_10mm": float(np.mean(np.abs(dev) < 10)),
            "amplitude_ratio": float(np.std(ds) / np.std(dl)) if np.std(dl) > 0 else float("nan")}


def increment_stats(truth_inc, est_inc) -> dict:
    tr, es = np.asarray(truth_inc, float), np.asarray(est_inc, float)
    m = np.isfinite(tr) & np.isfinite(es)
    tr, es = tr[m], es[m]
    dev = es - tr
    return {"pairs": int(m.sum()), "pearson_increments": _r(tr, es),
            "rmse_increments_mm": float(np.sqrt(np.mean(dev ** 2))) if len(dev) else float("nan"),
            "slope_through_origin": float(np.sum(tr * es) / np.sum(tr * tr)) if np.sum(tr * tr) > 0 else float("nan"),
            "within_quarter_wave": float(np.mean(np.abs(dev) < QUARTER_WAVE_MM)) if len(dev) else float("nan")}


# ----------------------------------------------------------------------------------------- chains
def chain(pairs: pd.DataFrame, est_col: str = "est_mm") -> pd.DataFrame:
    """Cumulate ``est_col`` along consecutive pairs (sorted by t1); a gap (t1 ≠ previous t2) starts a new segment."""
    rows, seg, level, prev = [], 0, 0.0, None
    for r in pairs.sort_values("t1").itertuples():
        if prev is None or r.t1 != prev:
            if prev is not None:
                seg += 1
            level = 0.0
            rows.append({"time": r.t1, "seg": seg, "est": 0.0})
        level += getattr(r, est_col)
        rows.append({"time": r.t2, "seg": seg, "est": level})
        prev = r.t2
    return pd.DataFrame(rows)


def level_views(c: pd.DataFrame, protocol: Protocol) -> list[tuple[str, pd.DataFrame]]:
    """``c`` has time, seg, est, truth (levels at the epochs). Each segment is re-zeroed at its first scored epoch;
    views: full span, per datum, per season of each year and datum (re-zeroed at the season start)."""
    c = c.dropna(subset=["truth"]).copy()
    for _, cs in c.groupby("seg"):
        c.loc[cs.index, "truth"] -= cs.truth.iloc[0]
        c.loc[cs.index, "est"] -= cs.est.iloc[0]
    c["datum"] = protocol.datum(c.time)
    views = [("full span", c)]
    views += [(f"datum {k}", c[c.datum == k]) for k in sorted(set(c.datum))]
    m0, m1 = protocol.season_months
    for y in sorted(c.time.dt.year.unique()):
        for k in sorted(set(c.datum)):
            s = c[(c.time.dt.year == y) & c.time.dt.month.between(m0, m1) & (c.datum == k)].copy()
            if len(s):
                s["truth"] -= s.truth.iloc[0]
                s["est"] -= s.est.iloc[0]
                views.append((f"season {y} {k}", s))
    return views


def winter_steps(c: pd.DataFrame, protocol: Protocol) -> pd.DataFrame:
    """Level step across each winter: last epoch in the season of year y → first epoch in the spring window of y+1,
    on the same datum (``c`` already re-zeroed per segment, as ``level_views`` does)."""
    c = c.dropna(subset=["truth"]).copy()
    c["datum"] = protocol.datum(c.time)
    m0, m1 = protocol.season_months
    s0, s1 = protocol.spring_window
    out = []
    for y in sorted(c.time.dt.year.unique()):
        a = c[(c.time.dt.year == y) & c.time.dt.month.between(m0, m1)]
        b = c[(c.time.dt.year == y + 1) & c.time.dt.month.between(s0, s1)]
        if not len(a) or not len(b):
            continue
        a, b = a.iloc[-1], b.iloc[0]
        if a.datum != b.datum:
            continue
        out.append({"winter": f"{y}/{str(y + 1)[2:]}", "from": a.time.date(), "to": b.time.date(),
                    "laser_change_mm": b.truth - a.truth, "insar_change_mm": b.est - a.est,
                    "difference_mm": (b.est - a.est) - (b.truth - a.truth), "same_segment": a.seg == b.seg})
    return pd.DataFrame(out)


# ----------------------------------------------------------------------------------------- the score
def score(pairs: pd.DataFrame, truth_level: pd.Series, protocol: Protocol, est_col: str = "est_mm",
          truth_inc_col: str | None = "truth_inc_mm", lockbox: bool = False) -> dict[str, pd.DataFrame]:
    """Score one estimate on one track.

    ``pairs``: t1, t2 (UTC), ``est_col`` (the estimate's increment, mm, same geometry as the truth) and, optionally,
    ``truth_inc_col`` (the truth's increment over the pair). ``truth_level``: the truth's level (mm), indexed by UTC
    time, read at each epoch (nearest hour). Returns tables: levels, winters, increments, seasons, years."""
    pairs = pairs.copy()
    for col in ("t1", "t2"):
        pairs[col] = pd.to_datetime(pairs[col], utc=True)
    guard_lockbox(pd.concat([pairs.t1, pairs.t2]), protocol, lockbox)
    c = chain(pairs, est_col)
    c["truth"] = truth_level.reindex(pd.DatetimeIndex(c.time).round("h")).to_numpy()
    lev = []
    for view, d in level_views(c, protocol):
        if len(d) >= protocol.min_epochs:
            lev.append({"view": view, "segments": d.seg.nunique(), **detrended_stats(d.time, d.truth, d.est)})
    z = c.dropna(subset=["truth"]).copy()
    for _, cs in z.groupby("seg"):
        z.loc[cs.index, "truth"] -= cs.truth.iloc[0]
        z.loc[cs.index, "est"] -= cs.est.iloc[0]
    out = {"levels": pd.DataFrame(lev), "winters": winter_steps(z, protocol)}
    if truth_inc_col and truth_inc_col in pairs:
        p = pairs.dropna(subset=[truth_inc_col])
        out["increments"] = pd.DataFrame([increment_stats(p[truth_inc_col], p[est_col])])
        p = p.assign(season=p.t1.dt.month.map(SEASON_OF_MONTH), year=p.t1.dt.year)
        out["seasons"] = pd.DataFrame([{"season": s, **increment_stats(g[truth_inc_col], g[est_col])}
                                       for s, g in p.groupby("season")])
        out["years"] = pd.DataFrame([{"year": int(y), **increment_stats(g[truth_inc_col], g[est_col]),
                                      "sign": int(np.sign(_r(g[truth_inc_col], g[est_col])))
                                      if np.isfinite(_r(g[truth_inc_col], g[est_col])) else 0}
                                     for y, g in p.groupby("year")])
    return out
