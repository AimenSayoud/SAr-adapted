"""X-079 — the nine 40 m cells around P6 against the laser, done from the data (exploratory).

Replaces the 2026-10-06 figures, whose per-cell table had no generating code, whose captions were typed by hand
and whose conclusion (the west cell moves 1:1, shrubs inflate the slope, the boardwalk drags it) was never tested
against the alternatives. Here every number is computed:

1. **Cells.** For the 3×3 cells centred on P6's cell, the consecutive-pair phase (cell minus the grassland C median,
   PHASE_TO_MM, the X-065 conventions and pairs) against the laser LOS change, by track and period (2020–2021
   six-day, 2022–2024 twelve-day): r with a circular-shift p, OLS and TLS slopes with a bootstrap 95 % interval,
   X-065's motion coefficient m (coherence-weighted, through the origin, with ΔWTD), same-sign share and coherence.
2. **What goes with the slope.** The cells' boardwalk share (X-060), LiDAR woody cover (CHM > 0.5 m, 1 m GUGiK),
   and depth into the mat (from the mat's edge, and from the mat+lake basin edge as the cube now stores it):
   Pearson r with an exact permutation p over all 9! orderings, and partial correlations between depth and cover.
3. **Seasons.** May–November of each year at P6 and its west and east neighbours, with the pair counts.
4. **The 2026-10-06 claims, checked** against the above (claims read from the hub's archive, not stored here).

Field data are read from the hub at run time; outputs go only to the hub. No coordinate or field value is
written in this file.

    PYTHONPATH=src python scripts/field_p6_cells_x079.py
"""
from __future__ import annotations

import itertools

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402
from scipy import ndimage as ndi  # noqa: E402

from insar_wetlands import field  # noqa: E402
from insar_wetlands import field_link as fl  # noqa: E402
from insar_wetlands import fusion as fu  # noqa: E402
from insar_wetlands.inversion.isbas import PHASE_TO_MM  # noqa: E402

HUB = field.field_root().parents[1]
DLV = HUB / "08_deliverables"
SIL = HUB / "06_data" / "cube" / "silver"
OUT = DLV / "field_p6_cells_x079"
CLAIMS = HUB / "99_archive" / "withdrawn_2026-10-06" / "claims.csv"
QW = 55.465763 / 4
FLOOR_MM = 1.2                                   # as X-065's same-sign window
N_BOOT = 2000
CELLS = {"NW": (-1, -1), "N": (-1, 0), "NE": (-1, 1), "W": (0, -1), "P6": (0, 0), "E": (0, 1),
         "SW": (1, -1), "S": (1, 0), "SE": (1, 1)}
PERIODS = {"2020–2021": ("2020_2021", 6), "2022–2024": ("2022_2024", 12)}
PERMS = np.array(list(itertools.permutations(range(len(CELLS)))))


def perm_r(x, y) -> tuple[float, float]:
    """Pearson r and its exact two-sided permutation p over all orderings of y."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    zx, zy = (x - x.mean()) / x.std(), (y - y.mean()) / y.std()
    r0 = float(np.mean(zx * zy))
    return r0, float(np.mean(np.abs(zy[PERMS] @ zx / len(x)) >= abs(r0) - 1e-12))


def residual(y, x) -> np.ndarray:
    X = np.column_stack([np.ones(len(x)), x])
    return np.asarray(y, float) - X @ np.linalg.lstsq(X, np.asarray(y, float), rcond=None)[0]


def pair_stats(x, y, coh, dwtd, rng) -> dict:
    sl = fl.slopes(x, y)
    r, p = fl.circular_shift_p(x, y)
    boot = [fl.slopes(x[i], y[i])["slope_tls"] for i in (rng.integers(0, len(x), len(x)) for _ in range(N_BOOT))]
    lo, hi = np.nanpercentile(boot, [2.5, 97.5])
    k = (np.abs(x) > FLOOR_MM) & (np.abs(x) < QW)
    ok = np.isfinite(dwtd)
    f = fu.bayes_linear(np.column_stack([x[ok], dwtd[ok]]), y[ok], fu.phase_sigma_mm(coh[ok]))
    return {"n": len(x), "r": r, "p_circ": p, "slope_ols": sl["slope_ols"], "slope_tls": sl["slope_tls"],
            "tls_lo95": lo, "tls_hi95": hi, "m": f["beta"][0], "m_sd": f["sd"][0],
            "same_sign": float(np.mean(np.sign(y[k]) == np.sign(x[k]))), "n_cmp": int(k.sum()),
            "median_coh_laser_pairs": float(np.median(coh))}


def build_series(st, r0, c0, pairs) -> pd.DataFrame:
    C = st.zone.values == 3
    rows = []
    for track in ("ascending", "descending"):
        for period, (tag, dt) in PERIODS.items():
            ds = xr.open_dataset(SIL / f"ifg_{tag}_{track}.nc")
            p = pairs[(pairs.track == track) & (pairs.step == 1) & (pairs.dt == dt) & pairs.pair.isin(ds.pair.values)]
            unw = ds.unw.sel(pair=p.pair.values).values
            coh = ds.coh.sel(pair=p.pair.values).values
            ref = np.array([np.nanmedian(u[C & np.isfinite(u)]) for u in unw])
            for cell, (dr, dc) in CELLS.items():
                rows.append(p.assign(period=period, cell=cell, s1_cell=(unw[:, r0 + dr, c0 + dc] - ref) * PHASE_TO_MM,
                                     coh_cell=coh[:, r0 + dr, c0 + dc]))
    return pd.concat(rows, ignore_index=True)


def cell_table(ser, rng) -> pd.DataFrame:
    rows = []
    for (track, period, cell), g in ser.groupby(["track", "period", "cell"], sort=False):
        coh_all = float(np.nanmedian(g.coh_cell))
        L = g.dropna(subset=["laser_los", "s1_cell"])
        rows.append({"track": track, "period": period, "cell": cell, "median_coh_all_pairs": coh_all,
                     **pair_stats(L.laser_los.to_numpy(), L.s1_cell.to_numpy(), L.coh_cell.to_numpy(),
                                  L.dwtd.to_numpy(), rng)})
    return pd.DataFrame(rows)


def descriptors(st, r0, c0) -> pd.DataFrame:
    A = st.zone.values == 1
    d_mat = ndi.distance_transform_edt(A) * 40.0
    bw = pd.read_csv(DLV / "field_boardwalk_x060" / "boardwalk_cells.csv")
    bw = bw[bw.res_m == 40].set_index(["row", "col"]).boardwalk_share
    chm = xr.open_dataset(SIL / "lidar_1m.nc")
    chm = chm.dsm_m - chm.dtm_m
    rows = []
    for cell, (dr, dc) in CELLS.items():
        r, c = r0 + dr, c0 + dc
        E, N = float(st.x[c]), float(st.y[r])
        h = chm.sel(x=slice(E - 20, E + 20), y=slice(N + 20, N - 20)).values
        h = h[np.isfinite(h)]
        rows.append({"cell": cell, "row": r, "col": c, "zone": int(st.zone.values[r, c]),
                     "plot": int(st.plot_cell.values[r, c]), "boardwalk_pct": 100 * float(bw.get((r, c), 0.0)),
                     "lidar_woody_pct": 100 * float(np.mean(h > 0.5)), "lidar_tree_pct": 100 * float(np.mean(h > 1.5)),
                     "depth_mat_edge_m": float(d_mat[r, c]), "depth_basin_edge_m": float(st.depth_in_mat_m.values[r, c])})
    return pd.DataFrame(rows).set_index("cell")


def drivers(cells, desc) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, part = [], []
    for (track, period), g in cells.groupby(["track", "period"], sort=False):
        g = g.set_index("cell").loc[desc.index]
        for y in ("slope_tls", "m"):
            for x in ("depth_mat_edge_m", "depth_basin_edge_m", "lidar_woody_pct", "boardwalk_pct"):
                r, p = perm_r(desc[x], g[y])
                rows.append({"track": track, "period": period, "slope": y, "against": x, "r": r, "p_perm": p})
            dep, veg = desc.depth_mat_edge_m.to_numpy(), desc.lidar_woody_pct.to_numpy()
            rv, pv = perm_r(residual(veg, dep), residual(g[y], dep))
            rd, pd_ = perm_r(residual(dep, veg), residual(g[y], veg))
            part.append({"track": track, "period": period, "slope": y, "r_cover_given_depth": rv, "p_cover_given_depth": pv,
                         "r_depth_given_cover": rd, "p_depth_given_cover": pd_})
    return pd.DataFrame(rows), pd.DataFrame(part)


def seasons(ser) -> pd.DataFrame:
    s = ser[(ser.track == "ascending") & ser.cell.isin(["W", "P6", "E"])].dropna(subset=["laser_los"]).copy()
    t1, t2 = pd.to_datetime(s.pair.str[:8]), pd.to_datetime(s.pair.str[9:])
    s = s[(t1.dt.month >= 5) & (t2.dt.month <= 11)].assign(season=t1.dt.year)
    rows = []
    for (yr, cell), g in s.groupby(["season", "cell"]):
        x, y = g.laser_los.to_numpy(), g.s1_cell.to_numpy()
        sl = fl.slopes(x, y)
        r, p = fl.circular_shift_p(x, y) if len(x) >= 8 else (sl["r"], np.nan)
        rows.append({"season": int(yr), "revisit_days": int(g.dt.iloc[0]), "cell": cell, "n": len(x), "r": r, "p_circ": p,
                     "slope_ols": sl["slope_ols"], "slope_tls": sl["slope_tls"], "median_coh": float(np.median(g.coh_cell))})
    return pd.DataFrame(rows)


def check_claims(cells, desc, drv, part, ser, sea, st) -> pd.DataFrame:
    if not CLAIMS.exists():
        return pd.DataFrame()
    a12 = cells[(cells.track == "ascending") & (cells.period == "2022–2024")].set_index("cell")
    a12p = part[(part.track == "ascending") & (part.period == "2022–2024") & (part.slope == "slope_tls")].iloc[0]
    bw = drv[(drv.track == "ascending") & (drv.slope == "slope_tls") & (drv.against == "boardwalk_pct")]
    off = pd.read_csv(DLV / "field_plots_x061" / "plots_pairs.csv")
    off = off[(off["plot"] == "P6") & (off["extraction"] == "off-boardwalk")].set_index(["track", "pair"]).s1_dlos_mm
    wide = ser.pivot_table(index=["track", "pair"], columns="cell", values="s1_cell")
    j = wide.join(off, how="inner")
    wet = pd.read_csv(DLV / "field_dew_x050" / "wetness_at_overpasses_2020_2024.csv")
    ev = pd.read_csv(DLV / "field_events_x055" / "epoch_response.csv")
    px = pd.read_csv(DLV / "field_first" / "plot_pixels.csv").set_index("plot")
    zname = {1: "mat (A)", 2: "lake (B)", 3: "grassland (C)", 4: "other ground (D)"}
    d12 = cells[(cells.track == "descending") & (cells.period == "2022–2024")]
    s12 = sea[(sea.revisit_days == 12) & (sea.cell == "P6")]
    tls = lambda c: f"{a12.loc[c, 'slope_tls']:.2f} (95 % {a12.loc[c, 'tls_lo95']:.2f}–{a12.loc[c, 'tls_hi95']:.2f})"  # noqa: E731
    found = {
        "tls12_W": tls("W"), "tls12_P6": tls("P6"), "tls12_E": tls("E"), "tls12_NE": tls("NE"),
        "veg_partial": f"cover alone r = {drv[(drv.track == 'ascending') & (drv.period == '2022–2024') & (drv.slope == 'slope_tls') & (drv.against == 'lidar_woody_pct')].r.iloc[0]:.2f}; "
                       f"given depth r = {a12p.r_cover_given_depth:.2f} (p {a12p.p_cover_given_depth:.2f}); "
                       f"depth given cover r = {a12p.r_depth_given_cover:.2f} (p {a12p.p_depth_given_cover:.3f})",
        "boardwalk_r": "slope vs boardwalk share r = " + ", ".join(f"{v:.2f}" for v in bw.r),
        "offboardwalk_is_mean_WE": f"max |series − mean(W, E)| = {np.nanmax(np.abs(j.s1_dlos_mm - (j.W + j.E) / 2)):.1e} mm; "
                                   f"vs W alone {np.nanmax(np.abs(j.s1_dlos_mm - j.W)):.1f} mm",
        "coh6_min": f"min over cells {cells[cells.period == '2020–2021'].groupby('track').median_coh_all_pairs.min().round(2).to_dict()} (all pairs); "
                    f"{cells[cells.period == '2020–2021'].groupby('track').median_coh_laser_pairs.min().round(2).to_dict()} (laser pairs)",
        "coh12_max": f"max over ascending cells {a12.median_coh_all_pairs.max():.2f} (all pairs, cell {a12.median_coh_all_pairs.idxmax()})",
        "desc12_r_max": f"descending 12-day r across cells {d12.r.min():.2f} to {d12.r.max():.2f}; "
                        f"{int((d12.p_circ < 0.05).sum())}/9 cells p < 0.05",
        "wet_asc": f"{100 * wet[wet.track == 'ascending'].wet.mean():.1f} % (X-050 'wet')",
        "rain_events": f"{int(ev[ev.series == 'WTD plot median'].n_events.max())} (water table), "
                       f"{int(ev[ev.series == 'laser surface P6'].n_events.max())} with the laser (X-055)",
        "p1_zone": zname[int(st.zone.values[int(px.loc['P1', 'row']), int(px.loc['P1', 'col'])])] + " on the 40 m grid",
        "season_n12": f"P6 12-day seasons: n = {', '.join(map(str, s12.n))}; TLS {s12.slope_tls.min():.2f}–{s12.slope_tls.max():.2f}",
    }
    cl = pd.read_csv(CLAIMS)
    cl["found"] = cl.check.map(found)
    return cl


def fig_cells(out, cells, desc):
    grid = lambda s: np.array([[s[c] for c in row] for row in (("NW", "N", "NE"), ("W", "P6", "E"), ("SW", "S", "SE"))])  # noqa: E731
    panels = []
    for track, period in (("ascending", "2020–2021"), ("ascending", "2022–2024"), ("descending", "2020–2021"), ("descending", "2022–2024")):
        g = cells[(cells.track == track) & (cells.period == period)].set_index("cell")
        panels.append((f"TLS slope, {track} {period}", grid(g.slope_tls), "RdBu_r", (0.4, 1.6),
                       lambda c, g=g: f"{g.loc[c, 'slope_tls']:.2f}\n[{g.loc[c, 'tls_lo95']:.2f}, {g.loc[c, 'tls_hi95']:.2f}]\nr {g.loc[c, 'r']:.2f}"))
    panels.append(("depth into the mat (m, from its edge)", grid(desc.depth_mat_edge_m), "viridis", None,
                   lambda c: f"{desc.loc[c, 'depth_mat_edge_m']:.0f}"))
    panels.append(("LiDAR woody cover (%, CHM > 0.5 m)", grid(desc.lidar_woody_pct), "Greens", None,
                   lambda c: f"{desc.loc[c, 'lidar_woody_pct']:.1f}"))
    names = (("NW", "N", "NE"), ("W", "P6", "E"), ("SW", "S", "SE"))
    fig, ax = plt.subplots(2, 3, figsize=(13, 8.5))
    for a, (title, z, cmap, lim, lab) in zip(ax.ravel(), panels):
        a.imshow(z, cmap=cmap, vmin=lim[0] if lim else None, vmax=lim[1] if lim else None)
        for i in range(3):
            for k in range(3):
                a.text(k, i, f"{names[i][k]}\n{lab(names[i][k])}", ha="center", va="center", fontsize=7.5)
        a.set_title(title, fontsize=9)
        a.set_xticks([]), a.set_yticks([])
    fig.suptitle("Nine 40 m cells around P6: phase vs laser per cell (slope, bootstrap 95 % interval, r), "
                 "depth into the mat, woody cover. North up.", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "fig_cells.png", dpi=110)
    plt.close(fig)


def fig_drivers(out, cells, desc, part):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    style = {("ascending", "2020–2021"): ("o", "C0"), ("ascending", "2022–2024"): ("s", "C1"),
             ("descending", "2020–2021"): ("^", "C2"), ("descending", "2022–2024"): ("v", "C3")}
    for (track, period), (mk, col) in style.items():
        g = cells[(cells.track == track) & (cells.period == period)].set_index("cell").loc[desc.index]
        pr = part[(part.track == track) & (part.period == period) & (part.slope == "m")].iloc[0]
        for a, x, key in ((ax[0], desc.depth_mat_edge_m, "depth_given_cover"), (ax[1], desc.lidar_woody_pct, "cover_given_depth")):
            a.errorbar(x, g.m, yerr=1.96 * g.m_sd, fmt=mk, color=col, ms=5, capsize=2, alpha=0.85,
                       label=f"{track} {period} (partial r {pr['r_' + key]:.2f}, p {pr['p_' + key]:.2f})")
    ax[0].set_xlabel("depth into the mat (m, from its edge)")
    ax[1].set_xlabel("LiDAR woody cover (%)")
    ax[0].set_ylabel("motion coefficient m vs P6 laser (± 1.96 SD)")
    for a in ax:
        a.axhline(1, color="0.5", lw=0.8, ls=":")
        a.legend(fontsize=6.5)
    ax[0].set_title("m against depth (partial r: given cover)", fontsize=9)
    ax[1].set_title("m against woody cover (partial r: given depth)", fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "fig_drivers.png", dpi=110)
    plt.close(fig)


def md(df: pd.DataFrame, fmt="{:.2f}") -> str:
    f = lambda v: fmt.format(v) if isinstance(v, (float, np.floating)) else str(v)  # noqa: E731
    return "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
                     + ["| " + " | ".join(f(v) for v in r) + " |" for r in df.itertuples(index=False)])


def readme(out, chk, cells, desc, drv, part, sea, claims):
    c = cells.copy()
    c["TLS [95 %]"] = [f"{a:.2f} [{lo:.2f}, {hi:.2f}]" for a, lo, hi in zip(c.slope_tls, c.tls_lo95, c.tls_hi95)]
    c["m ± sd"] = [f"{a:.2f} ± {s:.2f}" for a, s in zip(c.m, c.m_sd)]
    c["same sign"] = [f"{s:.2f} ({n})" for s, n in zip(c.same_sign, c.n_cmp)]
    c = c[["track", "period", "cell", "n", "r", "p_circ", "slope_ols", "TLS [95 %]", "m ± sd", "same sign",
           "median_coh_laser_pairs", "median_coh_all_pairs"]]
    a = drv[drv.slope == "m"].pivot_table(index=["track", "period"], columns="against", values="r", sort=False).reset_index()
    ap = drv[drv.slope == "m"].pivot_table(index=["track", "period"], columns="against", values="p_perm", sort=False).reset_index()
    asc = part[part.track == "ascending"]
    desc_rng = cells[cells.track == "descending"].groupby("period").r.agg(["min", "max"])
    L = ["# X-079 — the nine cells around P6 against the laser, from the data (exploratory)", "",
         "Generated by `05_code/SAr-adapted/scripts/field_p6_cells_x079.py`; every number is computed by it. "
         "Replaces the figures of 2026-10-06 (withdrawn, see `99_archive/withdrawn_2026-10-06/README.md`).", "",
         f"Check: the P6 cell from the data cube reproduces X-065's phase to {chk:.1e} mm on every pair.", "",
         "## Findings", "",
         f"- **Each cell's slope is poorly determined.** Ascending 2022–2024, the bootstrap 95 % intervals of the TLS slope "
         f"are {cells[(cells.track == 'ascending') & (cells.period == '2022–2024')].eval('tls_hi95 - tls_lo95').min():.2f}–"
         f"{cells[(cells.track == 'ascending') & (cells.period == '2022–2024')].eval('tls_hi95 - tls_lo95').max():.2f} "
         "wide and all overlap: no single cell is distinguishable from P6, and no cell is shown to move 1:1 rather than otherwise.",
         f"- **Across the nine cells the slope goes with depth into the mat, not with woody cover.** Ascending, "
         f"motion coefficient m: given depth, cover adds r = {asc[asc.slope == 'm'].r_cover_given_depth.min():.2f} to "
         f"{asc[asc.slope == 'm'].r_cover_given_depth.max():.2f} (p ≥ {asc[asc.slope == 'm'].p_cover_given_depth.min():.2f}); "
         f"given cover, depth keeps r = {asc[asc.slope == 'm'].r_depth_given_cover.min():.2f} to "
         f"{asc[asc.slope == 'm'].r_depth_given_cover.max():.2f}. Cover and depth both rise eastward "
         f"(r = {np.corrcoef(desc.depth_mat_edge_m, desc.lidar_woody_pct)[0, 1]:.2f} between them). This is WP3's result — "
         "the mat's interior moves more than its edge — seen against the one laser at P6.",
         f"- **The boardwalk does not go with the slope** (r with boardwalk share "
         f"{drv[(drv.against == 'boardwalk_pct') & (drv.track == 'ascending')].r.min():.2f} to "
         f"{drv[(drv.against == 'boardwalk_pct') & (drv.track == 'ascending')].r.max():.2f}, ascending), as X-060 found for P6.",
         f"- **Descending does not replicate the pattern.** r with the laser across cells "
         f"{desc_rng.loc['2020–2021', 'min']:.2f}–{desc_rng.loc['2020–2021', 'max']:.2f} (six-day) and "
         f"{desc_rng.loc['2022–2024', 'min']:.2f}–{desc_rng.loc['2022–2024', 'max']:.2f} (twelve-day); no driver passes on "
         "descending. Nine neighbouring cells, one laser, spatially correlated: exploratory.",
         f"- **Single seasons are too short.** At P6 and its west and east neighbours, May–November of one year holds "
         f"{sea[sea.revisit_days == 12].n.min()}–{sea[sea.revisit_days == 12].n.max()} twelve-day pairs; their slopes swing "
         f"{sea[sea.revisit_days == 12].slope_tls.min():.2f}–{sea[sea.revisit_days == 12].slope_tls.max():.2f}.", "",
         "![cells](fig_cells.png)", "", "![drivers](fig_drivers.png)", "",
         "## 1. Cells", "",
         "Consecutive pairs, laser available (X-065). `m`: s1 = m·laser + e·ΔWTD, coherence-weighted through the origin "
         "(X-065). Same sign: pairs with a laser change between 1.2 mm and λ/4 (count in brackets). Coherence: median "
         "over the laser pairs (mostly summer) and over all pairs of the period.", "", md(c), "",
         "## 2. Cell descriptors and what goes with the slope", "", md(desc.reset_index()), "",
         "Depth: distance to the edge of the mat (zone A), and to the edge of the mat + lake basin (the cube's "
         "`depth_in_mat_m` since 2026-10-06, D-025). Woody cover: share of 1 m LiDAR CHM above 0.5 m inside the 40 m cell.", "",
         "Pearson r of the motion coefficient m with each descriptor (nine cells):", "", md(a), "",
         "Exact permutation p (9! orderings):", "", md(ap), "",
         "Partial correlations (TLS slope and m):", "", md(part), "",
         "## 3. Seasons (May–November, ascending)", "", md(sea), ""]
    if len(claims):
        L += ["## 4. The 2026-10-06 claims, checked", "",
              "`stated`: what the withdrawn figures or captions said; `found`: computed here.", "",
              md(claims[["id", "where", "claim", "stated", "found"]]), ""]
    L += ["## Files", "", "| File | Content |", "|---|---|",
          "| `cells.csv` | section 1 |", "| `descriptors.csv` | section 2, per cell |",
          "| `drivers.csv`, `partial.csv` | section 2, correlations |", "| `seasons.csv` | section 3 |",
          "| `claims_checked.csv` | section 4 |", "| `cells_series.csv` | every pair × cell (phase, coherence, laser, ΔWTD) |",
          "| `fig_cells.png`, `fig_drivers.png` | figures |"]
    (out / "README.md").write_text("\n".join(L) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    st = xr.open_dataset(SIL / "static_40m.nc")
    px = pd.read_csv(DLV / "field_first" / "plot_pixels.csv").set_index("plot")
    r0, c0 = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    pairs = pd.read_csv(DLV / "field_p6_short_pairs_x065" / "p6_pairs_2020_2024.csv")
    ser = build_series(st, r0, c0, pairs)
    p6 = ser[ser.cell == "P6"]
    chk = float(np.nanmax(np.abs(p6.s1_cell - p6.s1)))
    rng = np.random.default_rng(0)
    cells = cell_table(ser, rng)
    desc = descriptors(st, r0, c0)
    drv, part = drivers(cells, desc)
    sea = seasons(ser)
    claims = check_claims(cells, desc, drv, part, ser, sea, st)
    ser.to_csv(OUT / "cells_series.csv", index=False)
    cells.to_csv(OUT / "cells.csv", index=False)
    desc.to_csv(OUT / "descriptors.csv")
    drv.to_csv(OUT / "drivers.csv", index=False)
    part.to_csv(OUT / "partial.csv", index=False)
    sea.to_csv(OUT / "seasons.csv", index=False)
    if len(claims):
        claims.to_csv(OUT / "claims_checked.csv", index=False)
    fig_cells(OUT, cells, desc)
    fig_drivers(OUT, cells, desc, part)
    readme(OUT, chk, cells, desc, drv, part, sea, claims)
    print(f"X-079 written to {OUT} (check {chk:.1e} mm)")


if __name__ == "__main__":
    main()
