"""X-055: what a rain event does — to the water table, the mat surface and the radar.

Superposed epochs: every rain event (≥ 10 mm within 24 h, snow-free season April–October, at
least 5 days after the previous event) is aligned at its onset (t = 0, the first rainy hour of the
24-h window), and the hourly response is stacked from −2 to +10 days:
- WTD at each plot and the plot median, relative to t = 0 (dry-well hours out);
- the laser surface at P6 (snow-free), relative to t = 0.
Median and quartiles across events; the lag to the peak rise.

Radar: do short pairs (≤ 24 d) that span an event have lower season-cleaned coherence on the mat
than pairs that do not? Null: the event onsets shifted together by a random number of days (keeps
their number and spacing, breaks their link to the pairs). Zone-mean coherence per pair comes from
the X-053 arc table. Outputs only to the private hub.

    PYTHONPATH=src python scripts/field_events_x055.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from insar_wetlands import field  # noqa: E402
from insar_wetlands import field_link as fl  # noqa: E402

RAIN_MM, WINDOW_H, SEP_D = 10.0, 24, 5
BEFORE_H, AFTER_H = 48, 240
SEASON = range(4, 11)       # April–October: liquid rain, no snow on the laser


def rain_events(rain: pd.Series) -> pd.DatetimeIndex:
    """Onsets: first rainy hour of a 24-h window reaching RAIN_MM, ≥ SEP_D days after the previous event."""
    tot = rain.rolling(f"{WINDOW_H}h").sum()
    onsets, last = [], None
    for t in tot.index[tot >= RAIN_MM]:
        if t.month not in SEASON or (last is not None and t - last < pd.Timedelta(days=SEP_D)):
            continue
        win = rain[t - pd.Timedelta(hours=WINDOW_H - 1): t]
        onset = win[win > 0].index[0]
        onsets.append(onset)
        last = onset
    return pd.DatetimeIndex(onsets)


def epochs(series: pd.Series, onsets: pd.DatetimeIndex) -> pd.DataFrame:
    """Rows: lag (h); columns: events; values relative to t = 0."""
    lags = np.arange(-BEFORE_H, AFTER_H + 1)
    cols = {}
    for t0 in onsets:
        seg = series.reindex(t0 + pd.to_timedelta(lags, unit="h"))
        if seg.notna().mean() >= 0.8 and np.isfinite(seg.iloc[BEFORE_H]):
            cols[t0] = seg.to_numpy() - seg.iloc[BEFORE_H]
    return pd.DataFrame(cols, index=lags)


def summarize(ep: pd.DataFrame, name: str) -> pd.DataFrame:
    q = ep.quantile([0.25, 0.5, 0.75], axis=1).T
    q.columns = ["q1", "median", "q3"]
    return q.assign(series=name, n_events=ep.shape[1]).rename_axis("lag_h").reset_index()


def spans_event(onset_days: pd.DatetimeIndex, d1: np.ndarray, d2: np.ndarray) -> np.ndarray:
    """True for pairs whose interval (d1, d2] contains an event onset."""
    e = onset_days.to_numpy()
    return np.array([((e > a) & (e <= b)).any() for a, b in zip(d1, d2)])


def split_diff(values: np.ndarray, hit: np.ndarray) -> float:
    """Mean of values without an event minus with one."""
    return float(values[~hit].mean() - values[hit].mean()) if hit.any() and (~hit).any() else np.nan


def main(argv=None):
    hub = field.field_root().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--arcs", default=str(hub / "08_deliverables" / "field_pixel_x053" / "arc_pairs.csv"))
    ap.add_argument("--out", default=str(hub / "08_deliverables" / "field_events_x055"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    wtd = field.load_wtd_hourly()
    laser = field.load_laser()
    surface = laser.surface_cm.where(~field.snow_mask(wtd["Air_2m"], laser.index)).reindex(wtd.index)
    cens = pd.DataFrame({p: field.censored_flag(wtd[p]) for p in field.PLOTS})
    w = wtd[field.PLOTS].where(~cens)
    onsets = rain_events(wtd["Rain_mm_Tot"])
    pd.DataFrame({"onset_utc": onsets, "rain_24h_mm": [float(wtd["Rain_mm_Tot"][t: t + pd.Timedelta(hours=WINDOW_H - 1)].sum()) for t in onsets]}
                 ).to_csv(out / "events.csv", index=False)

    parts, peaks = [], []
    for name, s in [("WTD plot median", w.median(axis=1)), *[(f"WTD {p}", w[p]) for p in field.PLOTS], ("laser surface P6", surface)]:
        ep = epochs(s, onsets)
        if ep.shape[1] < 5:
            continue
        sm = summarize(ep, name)
        parts.append(sm)
        after = sm[sm.lag_h >= 0]
        k = after["median"].idxmax()
        peaks.append({"series": name, "n_events": ep.shape[1], "peak_rise_cm": float(after.loc[k, "median"]),
                      "lag_to_peak_h": int(after.loc[k, "lag_h"]),
                      "rise_after_24h_cm": float(sm.loc[sm.lag_h == 24, "median"].iloc[0])})
    resp = pd.concat(parts, ignore_index=True)
    resp.to_csv(out / "epoch_response.csv", index=False)
    pk = pd.DataFrame(peaks)
    pk.to_csv(out / "peak_response.csv", index=False)

    # buoyancy seen through events: surface rise per cm of water-table rise at P6 (median events)
    radar = []
    arcs_path = Path(a.arcs)
    if arcs_path.exists():
        arcs = pd.read_csv(arcs_path, parse_dates=["date1", "date2"])
        days = pd.DatetimeIndex(onsets.tz_localize(None).normalize())
        for track, g in arcs.groupby("track"):
            g = g[(g.dt_days <= 24) & ~g.frozen_any].dropna(subset=["coh_A"]).copy()
            mid = g.date1 + (g.date2 - g.date1) / 2
            res = fl.season_residual(g.coh_A.to_numpy(), mid, g.dt_days.to_numpy())
            d1, d2 = g.date1.to_numpy(), g.date2.to_numpy()
            hit = spans_event(days, d1, d2)
            d0 = split_diff(res, hit)
            rng = np.random.default_rng(0)
            null = np.array([split_diff(res, spans_event(days + pd.Timedelta(days=int(k)), d1, d2))
                             for k in rng.integers(30, 1800, 1000)])
            null = null[np.isfinite(null)]
            radar.append({"track": track, "pairs": len(g), "pairs_spanning_event": int(hit.sum()),
                          "coh_A_no_event_minus_event": d0, "p_shifted_events": float((np.sum(np.abs(null) >= abs(d0)) + 1) / (len(null) + 1))})
        pd.DataFrame(radar).to_csv(out / "radar_event_test.csv", index=False)

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
    for name, c in (("WTD plot median", "C0"), ("laser surface P6", "k")):
        s = resp[resp.series == name]
        if s.empty:
            continue
        ax[0].plot(s.lag_h / 24, s["median"], color=c, label=f"{name} (n={int(s.n_events.iloc[0])})")
        ax[0].fill_between(s.lag_h / 24, s.q1, s.q3, color=c, alpha=0.15)
    ax[0].axvline(0, color="k", lw=0.6); ax[0].set_xlabel("days from rain onset"); ax[0].set_ylabel("change from t = 0 (cm)")
    ax[0].legend(fontsize=8); ax[0].set_title(f"Response to rain events (≥ {RAIN_MM:.0f} mm / 24 h, Apr–Oct)", fontsize=9)
    for p in field.PLOTS:
        s = resp[resp.series == f"WTD {p}"]
        if not s.empty:
            ax[1].plot(s.lag_h / 24, s["median"], lw=1, label=p)
    ax[1].axvline(0, color="k", lw=0.6); ax[1].set_xlabel("days from rain onset"); ax[1].legend(fontsize=7, ncol=3)
    ax[1].set_title("WTD per plot (median over events)", fontsize=9)
    fig.tight_layout(); fig.savefig(out / "fig_rain_events.png", dpi=140); plt.close(fig)
    pd.set_option("display.width", 200)
    print(len(onsets), "events"); print(pk.round(2).to_string(index=False)); print(pd.DataFrame(radar).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
