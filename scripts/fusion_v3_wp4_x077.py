"""X-077 — fusion v3 WP4: phase cycles and the L-band anchor (exploratory, Mac-safe, reads the cube)

E4.1 Cycle choice. For each P6 pair, the wrapped phase relative to the hard targets (phasor) admits the values
     y_k = y_wrapped + k·λ/2 (line of sight). SNAPHU picked one (the unwrapped phase). A prior from the driver
     (buoyancy × water-table change, the bucket model where no logger; WP1) picks the k nearest the predicted change.
     Truth: the laser — a choice is right if it lands within λ/4 of the laser change. Compared for all pairs and for
     those whose motion exceeds λ/4 (where SNAPHU is expected to fail).
E4.1b Unwrapper check of WP3's edge gradient: the same depth-vs-raft test on WRAPPED departures (each pixel's wrapped
     phase relative to the mat's phasor mean; departures are a few mm, far below λ/2, so no unwrapping is needed).
E4.2 L-band. NISAR pairs at P6 (hard-target reference) against the laser where both exist (late 2025), and against
     the Sentinel-1 short-pair chain over the same interval; where C-band and L-band disagree by whole C-band cycles.

    PYTHONPATH=src python scripts/fusion_v3_wp4_x077.py
Outputs only to the hub: ``08_deliverables/fusion_v3/wp4/``.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy import stats

from insar_wetlands.cube.core import WAVELENGTH_MM, rad_to_los_mm

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fusion_v3_wp2_x075 as W2  # noqa: E402
import fusion_v3_wp3_x076 as W3  # noqa: E402

HUB = Path(__file__).resolve().parents[3]
G, S = HUB / "06_data" / "cube" / "gold", HUB / "06_data" / "cube" / "silver"
OUT = HUB / "08_deliverables" / "fusion_v3" / "wp4"
CYCLE = WAVELENGTH_MM / 2
RNG = np.random.default_rng(77)


def e41() -> tuple[pd.DataFrame, pd.DataFrame]:
    p = W2.p6_pairs()
    p["driver_dh_mm"] = W2.driver_prediction(p)
    q = p[p.usable & (p.revisit_days <= 24)].copy()
    # wrapped P6 phase relative to the hard targets' phasor, in mm (−φ·λ/4π), wrapped to ±λ/4
    w = np.angle(np.exp(1j * (q.wrapped_1x1 - q.ref_hard_wrapped)))
    q["y_wrapped"] = -rad_to_los_mm(w)
    q["y_snaphu"] = q.y
    k = np.round((q.driver_dh_mm.fillna(0) - q.y_wrapped) / CYCLE)
    q["y_prior"] = q.y_wrapped + k * CYCLE
    k0 = np.round((0 - q.y_wrapped) / CYCLE)
    q["y_zero"] = q.y_wrapped + k0 * CYCLE                  # the cycle nearest zero (no information)
    big = q.d_laser_los_mm.abs() > WAVELENGTH_MM / 4
    rows = []
    for (rv_lab, sel) in (("all ≤ 24 d", np.ones(len(q), bool)), ("6 d", q.revisit_days == 6),
                          ("12 d", q.revisit_days == 12), ("18–24 d", q.revisit_days >= 18)):
        for mlab, m in (("all", np.ones(len(q), bool)), ("motion > λ/4", big.to_numpy())):
            g = q[np.asarray(sel) & m]
            if len(g) < 5:
                continue
            row = {"pairs": rv_lab, "subset": mlab, "n": len(g)}
            for name in ("y_snaphu", "y_prior", "y_zero"):
                err = (g[name] - g.d_laser_los_mm).abs()
                row[f"right_cycle_{name[2:]}"] = float(np.mean(err < WAVELENGTH_MM / 4))
                row[f"rmse_{name[2:]}_mm"] = float(np.sqrt(np.mean(err ** 2)))
            rows.append(row)
    return pd.DataFrame(rows), q


def e41b() -> pd.DataFrame:
    rows = []
    for period, rv in (("2020_2021", 6), ("2025", 6), ("2022_2024", 12)):
        for track in ("ascending", "descending"):
            ds = xr.open_dataset(G / f"mat_pairs_{period}_{track}.nc")
            keep = (ds.revisit_days.values == rv) & np.array(["S1D" not in x for x in ds.platforms.values.astype(str)])
            Wr = ds.wrapped.values[:, keep].T.astype(float)
            U = ds.unw.values[:, keep].T.astype(float)
            zone = ds.st_zone.values
            hard = ds.st_hard_target.values == 1
            mat = zone == 1
            ref = np.nanmedian(U[:, hard], 1)
            raft = np.nanmedian(-rad_to_los_mm(U[:, mat] - ref[:, None]), 1)
            Z = np.exp(1j * Wr[:, mat])
            C = np.nan_to_num(ds.coh.values[:, keep].T[:, mat], nan=0)
            mean = np.angle(np.nansum(np.where(np.isfinite(Z), Z * C, 0), 1))
            dep = -rad_to_los_mm(np.angle(Z * np.exp(-1j * mean[:, None])))   # wrapped departures, mm
            rr, cc = ds.row.values[mat], ds.col.values[mat]
            depth = ds.st_depth_in_mat_m.values[mat]
            cohm = np.nanmedian(ds.coh.values[:, keep], 1)[mat]
            Sh = W3.shapes_for(rr, cc, depth, ds.st_lidar_microrelief_sd_m.values[mat], ds.st_lidar_vegh_mean_m.values[mat], cohm)
            for lab, sel in (("all mat pixels", np.ones(mat.sum(), bool)), ("interior only (≥ 80 m)", depth >= 80)):
                A = W3.coefficients(dep[:, sel], W3.shapes_for(rr[sel], cc[sel], depth[sel],
                                                               ds.st_lidar_microrelief_sd_m.values[mat][sel],
                                                               ds.st_lidar_vegh_mean_m.values[mat][sel], cohm[sel])) \
                    if lab != "all mat pixels" else W3.coefficients(dep, Sh)
                rows.append({"period": period, "track": track, "revisit_days": rv, "subset": lab, "phase": "wrapped",
                             **W3.slope_test(A[:, 0], raft)})
    return pd.DataFrame(rows)


def e42(q: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    f = G / "plots_pair_nisar.csv.gz"
    if not f.exists():
        return pd.DataFrame(), pd.DataFrame()
    n = pd.read_csv(f, low_memory=False, parse_dates=["t1", "t2"])
    lam = float(xr.open_dataset(S / "nisar_lband.nc").attrs["wavelength_m"]) * 1000
    p6 = n[n["plot"] == "P6"].copy()
    p6["y_l_mm"] = -(p6.unw_3x3 - p6.ref_hard_unw) * lam / (4 * np.pi)
    p6["y_l_stable_mm"] = -(p6.unw_3x3 - p6.ref_stable_unw) * lam / (4 * np.pi)
    las = p6.dropna(subset=["d_laser_vertical_mm", "y_l_mm"])
    t1 = las[["granule", "track_no", "direction", "t1", "t2", "revisit_days", "coh_3x3", "y_l_mm", "y_l_stable_mm",
              "d_laser_vertical_mm"]].copy()
    # Sentinel-1 6-day chain (hard-target reference) over the NISAR interval at P6, ascending
    s1 = W2.p6_pairs()
    s1 = s1[(s1.revisit_days == 6)]
    rows = []
    for r in p6.itertuples():
        for track in ("ascending", "descending"):
            g = s1[(s1.track == track) & (s1.t1 >= r.t1 - pd.Timedelta(days=3)) & (s1.t2 <= r.t2 + pd.Timedelta(days=3))]
            if len(g) and np.isfinite(r.y_l_mm):
                c_sum = g.y.sum()
                rows.append({"granule": r.granule, "t1": r.t1, "t2": r.t2, "track_s1": track, "s1_pairs": len(g),
                             "lband_los_mm": r.y_l_mm, "cband_chain_los_mm": c_sum,
                             "difference_mm": r.y_l_mm - c_sum,
                             "difference_in_cband_cycles": (r.y_l_mm - c_sum) / CYCLE})
    return t1, pd.DataFrame(rows)


def md(df: pd.DataFrame, nd: int = 2) -> str:
    if df is None or not len(df):
        return "(none)"
    dd = df.copy()
    for c in dd.columns:
        if dd[c].dtype.kind == "f":
            dd[c] = dd[c].round(nd)
    lines = ["| " + " | ".join(map(str, dd.columns)) + " |", "|" + "---|" * len(dd.columns)]
    for r in dd.itertuples(index=False):
        lines.append("| " + " | ".join("" if (isinstance(v, float) and np.isnan(v)) else str(v) for v in r) + " |")
    return "\n".join(lines)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t41, q = e41()
    t41b = e41b()
    t42a, t42b = e42(q)
    for n, t in (("e41_cycle_choice", t41), ("e41b_wrapped_edge_gradient", t41b), ("e42_nisar_vs_laser", t42a),
                 ("e42_nisar_vs_sentinel1", t42b)):
        t.to_csv(OUT / f"{n}.csv", index=False)
    q[["track", "pair", "t1", "t2", "revisit_days", "d_laser_los_mm", "driver_dh_mm", "y_wrapped", "y_snaphu", "y_prior"]] \
        .to_csv(OUT / "e41_pairs.csv", index=False)
    rl = ""
    if len(t42a) >= 3:
        r = stats.pearsonr(t42a.y_l_mm, t42a.d_laser_vertical_mm)[0]
        rl = f"L-band vs laser (vertical), {len(t42a)} pairs: r = {r:.2f}."
    txt = ["# X-077 — fusion v3 WP4: phase cycles and the L-band anchor (exploratory)", "",
           "Generated by `05_code/SAr-adapted/scripts/fusion_v3_wp4_x077.py` from `06_data/cube`; every number is computed.", "",
           "## E4.1 Which cycle? SNAPHU, the driver prior, or the cycle nearest zero", "",
           "Share of usable P6 pairs whose chosen value lands within λ/4 of the laser change (`right_cycle_*`), and the "
           "error. The driver prior = buoyancy × water-table change (bucket model where no logger), projected to the line "
           "of sight; it never sees the laser.", "", md(t41), "",
           "## E4.1b Is WP3's edge gradient an unwrapping artefact?", "",
           "The same depth-vs-raft test on wrapped departures (no unwrapping involved):", "", md(t41b), "",
           "## E4.2 L-band", "", "NISAR at P6 (hard-target reference) where the laser exists:", "", md(t42a), "", rl, "",
           "L-band against the Sentinel-1 6-day chain over the same interval (P6, hard-target reference):", "", md(t42b), "",
           "Read: differences near whole numbers of C-band cycles (λ/2 ≈ 27.7 mm) are C-band cycle slips that one L-band "
           "pair would fix; small differences mean the two agree.", ""]
    (OUT / "README.md").write_text("\n".join(txt))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
