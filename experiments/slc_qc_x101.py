"""X-101 QC (Colab): do our coregistered CSLC pairs, looked and filtered like HyP3, agree with HyP3 on hard targets?

Without this nothing built on the CSLCs (X-102 onward) is comparable with the HyP3 results. For each consecutive pair
of 2021–2022-02 with two CSLCs (``slc_cslc_x101``): our interferogram s1·conj(s2) summed over the 5 m × 10 m samples in
each 40 m cube cell (32 looks), its coherence, and a Goldstein–Werner filtered copy (alpha 0.5, HyP3's dampening factor;
32 × 32-cell windows on the 40 m grid — an approximation of ISCE's filter, which runs before geocoding). The phase
difference to the HyP3 pair of the same dates (cube silver, wrapped) is measured on hard targets and stable cells and
compared with the Cramér–Rao noise expected from coherence. Rule fixed in the hub's pre-registration
(``08_deliverables/insar_full_test/preregistration/X-101_qc_v2.yaml``; v1: ``X-101.yaml``; the commit is passed as
``--prereg-commit``, the hub is not on the VM).

Writes to Drive ``insar_rzecin/slc_pilot/2021/qc/`` (tables + README); the hub copy is fetched with ``colab download``.

    (Colab) python -m insar_wetlands.register run slc_qc_x101 -- --prereg-commit <sha>
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from insar_wetlands.cube import cube_root
from insar_wetlands.paths import resolve_drive

CELL = 40.0
L_HYP3, L_OURS = 20, 32
ALPHA, WIN, STEP = 0.5, 32, 8
OFFSET_MAX = 0.3
PLACE_TOL_M = 20.0
DATE_RE = re.compile(r"(20\d{6})")


def cslc(path: Path, pol: str = "VV") -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    import h5py  # noqa: PLC0415
    with h5py.File(path, "r") as h:
        d = h["data"]
        return d[pol][()].astype("complex64"), d["x_coordinates"][()], d["y_coordinates"][()]


def axis_index(x: np.ndarray, y: np.ndarray, gx: np.ndarray, gy: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Cube row of each CSLC row (from y) and cube column of each CSLC column (from x); -1 outside the grid.
    Cube x ascending and y descending, values at cell centres."""
    row = np.floor(((gy[0] + CELL / 2) - y) / CELL).astype(int)
    col = np.floor((x - (gx[0] - CELL / 2)) / CELL).astype(int)
    row[(row < 0) | (row >= len(gy))] = -1
    col[(col < 0) | (col >= len(gx))] = -1
    return row, col


def look(s1: np.ndarray, s2: np.ndarray, rows: np.ndarray, cols: np.ndarray, shape) -> tuple[np.ndarray, np.ndarray]:
    """Complex interferogram and coherence on the cube grid from the CSLC samples of each cell."""
    rr, cc = np.meshgrid(rows, cols, indexing="ij")
    ok = (rr >= 0) & (cc >= 0) & np.isfinite(s1) & np.isfinite(s2)
    idx = (rr * shape[1] + cc)[ok]
    n = shape[0] * shape[1]
    prod = (s1 * np.conj(s2))[ok]
    re_ = np.bincount(idx, prod.real, n)
    im_ = np.bincount(idx, prod.imag, n)
    p1 = np.bincount(idx, (np.abs(s1) ** 2)[ok], n)
    p2 = np.bincount(idx, (np.abs(s2) ** 2)[ok], n)
    cnt = np.bincount(idx, minlength=n)
    z = (re_ + 1j * im_).reshape(shape)
    coh = (np.abs(z) / np.sqrt(p1 * p2).reshape(shape))
    z[cnt.reshape(shape) < L_OURS // 2] = np.nan
    coh[cnt.reshape(shape) < L_OURS // 2] = np.nan
    return z, coh


def goldstein(z: np.ndarray, alpha: float = ALPHA, win: int = WIN, step: int = STEP) -> np.ndarray:
    """Goldstein–Werner filter: in overlapping windows, weight the spectrum by its smoothed magnitude ** alpha."""
    from scipy.ndimage import uniform_filter  # noqa: PLC0415
    a = np.where(np.isfinite(z), z / np.maximum(np.abs(z), 1e-12), 0)
    ny, nx = a.shape
    pad_y, pad_x = (-ny) % step + win, (-nx) % step + win
    a = np.pad(a, ((win // 2, pad_y), (win // 2, pad_x)))
    out = np.zeros_like(a)
    wsum = np.zeros(a.shape)
    taper = np.outer(np.hanning(win), np.hanning(win)) + 1e-6
    for i in range(0, a.shape[0] - win + 1, step):
        for j in range(0, a.shape[1] - win + 1, step):
            f = np.fft.fft2(a[i:i + win, j:j + win])
            s = uniform_filter(np.abs(f), 3)
            h = (s / max(s.max(), 1e-12)) ** alpha
            out[i:i + win, j:j + win] += np.fft.ifft2(f * h) * taper
            wsum[i:i + win, j:j + win] += taper
    out = (out / np.maximum(wsum, 1e-12))[win // 2:win // 2 + ny, win // 2:win // 2 + nx]
    return np.where(np.isfinite(z), out, np.nan)


def crb(g: np.ndarray, looks: int) -> np.ndarray:
    g = np.clip(g, 0.05, 0.999)
    return np.sqrt((1 - g ** 2) / (2 * looks * g ** 2))


def circ_stats(d: np.ndarray) -> tuple[float, float]:
    """(circular mean, circular std) of phases (rad)."""
    d = d[np.isfinite(d)]
    if len(d) < 3:
        return np.nan, np.nan
    m = np.mean(np.exp(1j * d))
    return float(np.angle(m)), float(np.sqrt(-2 * np.log(max(abs(m), 1e-12))))


def hyp3_pairs(track: str) -> dict[tuple[str, str], tuple[np.ndarray, np.ndarray]]:
    out = {}
    for s in ("2020_2021", "2022_2024"):
        f = cube_root() / "silver" / f"ifg_{s}_{track}.nc"
        if not f.exists():
            continue
        with xr.open_dataset(f) as d:
            t1 = pd.to_datetime(d.t1.values).strftime("%Y%m%d")
            t2 = pd.to_datetime(d.t2.values).strftime("%Y%m%d")
            for i, (a, b) in enumerate(zip(t1, t2)):
                out[(a, b)] = (d.wrapped.isel(pair=i).values, d.coh.isel(pair=i).values)
    return out


def qc_track(track: str, st: xr.Dataset) -> pd.DataFrame:
    folder = resolve_drive() / "slc_pilot" / "2021" / "cslc" / track
    files = {DATE_RE.search(f.name).group(1): f for f in sorted(folder.glob("*.h5"))
             if "static" not in f.name.lower() and DATE_RE.search(f.name)}
    dates = sorted(files)
    hy = hyp3_pairs(track)
    gx, gy = st.x.values, st.y.values
    shape = (len(gy), len(gx))
    hard = st.hard_target.values.astype(bool)
    stable = st.stable.values.astype(bool) & ~hard if "stable" in st else np.zeros(shape, bool)
    masks = {"hard_targets": hard, "stable": stable, "mat": st.zone.values == 1, "grassland": st.zone.values == 3}
    rows_out, cache = [], {}
    for a, b in zip(dates[:-1], dates[1:]):
        if (a, b) not in hy:
            continue
        for d_ in (a, b):
            if d_ not in cache:
                cache[d_] = cslc(files[d_])
        (s1, x, y), (s2, _, _) = cache[a], cache[b]
        rr, cc = axis_index(x, y, gx, gy)
        z, coh = look(s1, s2, rr, cc, shape)
        zf = goldstein(z)
        w_h, coh_h = hy[(a, b)]
        rec = {"track": track, "t1": a, "t2": b, "days": (pd.Timestamp(b) - pd.Timestamp(a)).days}
        for name, ours in (("filtered", np.angle(zf)), ("unfiltered", np.angle(z))):
            diff = np.angle(np.exp(1j * (ours - w_h)))
            m_h, s_h = circ_stats(diff[hard])
            rec[f"{name}_spread_hard_rad"] = s_h
            for k in ("stable", "mat", "grassland"):
                m_k, s_k = circ_stats(np.angle(np.exp(1j * (diff[masks[k]] - m_h))))
                rec[f"{name}_offset_{k}_rad"] = m_k
                rec[f"{name}_spread_{k}_rad"] = s_k
        g = coh_h[hard]       # v1 yardstick: HyP3's coherence (biased high: computed on its filtered interferogram)
        rec["expected_spread_hard_rad_v1"] = float(np.nanmedian(np.sqrt(crb(g, L_HYP3) ** 2 + crb(g, L_OURS) ** 2)))
        g = coh[hard]         # v2 yardstick: our unfiltered coherence
        rec["expected_spread_hard_rad"] = float(np.nanmedian(np.sqrt(crb(g, L_HYP3) ** 2 + crb(g, L_OURS) ** 2)))
        for k, m in masks.items():
            rec[f"coh_hyp3_{k}"] = float(np.nanmedian(coh_h[m]))
            rec[f"coh_ours_{k}"] = float(np.nanmedian(coh[m]))
        rows_out.append(rec)
    return pd.DataFrame(rows_out)


def amp40(path: Path, gx: np.ndarray, gy: np.ndarray) -> np.ndarray:
    s, x, y = cslc(path)
    shape = (len(gy), len(gx))
    r, c = axis_index(x, y, gx, gy)
    rr, cc = np.meshgrid(r, c, indexing="ij")
    ok = (rr >= 0) & (cc >= 0) & np.isfinite(s)
    idx = (rr * shape[1] + cc)[ok]
    n = shape[0] * shape[1]
    pw = np.bincount(idx, (np.abs(s) ** 2)[ok], n)
    k = np.bincount(idx, minlength=n)
    a = (pw / np.maximum(k, 1)).reshape(shape)
    a[k.reshape(shape) < L_OURS // 2] = np.nan
    return 10 * np.log10(a)


def best_shift(a: np.ndarray, b: np.ndarray, rng: float = 3.0, step: float = 0.125) -> tuple[float, float, float]:
    """(r, dy, dx) in cells: the shift of ``a`` that best correlates with ``b``."""
    from scipy.ndimage import shift as ndshift  # noqa: PLC0415
    m = np.isfinite(a) & np.isfinite(b)
    fill = np.nan_to_num(a, nan=np.nanmean(a))
    out = (-2.0, 0.0, 0.0)
    for dy in np.arange(-rng, rng + 1e-9, step):
        for dx in np.arange(-rng, rng + 1e-9, step):
            r = np.corrcoef(ndshift(fill, (dy, dx), order=1, mode="nearest")[m], b[m])[0, 1]
            if r > out[0]:
                out = (float(r), float(dy), float(dx))
    return out


def placement(st: xr.Dataset, tracks: list[str]) -> pd.DataFrame:
    """Mean June–August amplitude against the airborne LiDAR canopy: shift in metres (east +, north +)."""
    rows = []
    for t in tracks:
        folder = resolve_drive() / "slc_pilot" / "2021" / "cslc" / t
        fs = [f for f in sorted(folder.glob("*.h5")) if "static" not in f.name.lower() and DATE_RE.search(f.name)
              and DATE_RE.search(f.name).group(1)[4:6] in ("06", "07", "08")]
        img = np.nanmean([amp40(f, st.x.values, st.y.values) for f in fs], axis=0)
        r, dy, dx = best_shift(img, st.lidar_vegh_p90_m.values.astype(float))
        rows.append({"track": t, "dates": len(fs), "r": r, "dx_m": dx * CELL, "dy_m": -dy * CELL})
    return pd.DataFrame(rows)


def placement_ok(pl: pd.DataFrame) -> bool:
    d = pl.set_index("track")
    if not {"ascending", "descending"} <= set(d.index):
        return False
    a, b = d.loc["ascending"], d.loc["descending"]
    return bool(np.sign(a.dx_m) != np.sign(b.dx_m) and abs(a.dx_m + b.dx_m) <= PLACE_TOL_M
                and abs(a.dy_m) <= PLACE_TOL_M and abs(b.dy_m) <= PLACE_TOL_M)


def verdict(t: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for track, g in t.groupby("track"):
        obs = g.filtered_spread_hard_rad.median()
        exp = g.expected_spread_hard_rad.median()
        off = float(np.angle(np.mean(np.exp(1j * g.filtered_offset_stable_rad.dropna())))) if len(g) else np.nan
        rows.append({"track": track, "pairs": len(g), "observed_spread_hard_rad": obs, "expected_spread_hard_rad": exp,
                     "stable_offset_rad": off,
                     "agrees": bool(obs <= exp and abs(off) <= OFFSET_MAX)})
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--prereg-commit", required=True, help="hub commit of preregistration/X-101.yaml")
    ap.add_argument("--tracks", nargs="+", default=["ascending", "descending"])
    a = ap.parse_args(argv)
    st = xr.open_dataset(cube_root() / "silver" / "static_40m.nc")
    t = pd.concat([qc_track(tr, st) for tr in a.tracks], ignore_index=True)
    v = verdict(t)
    v["expected_spread_hard_rad_v1"] = [t[t.track == tr].expected_spread_hard_rad_v1.median() for tr in v.track]
    pl = placement(st, a.tracks)
    ok = placement_ok(pl)
    v["placement_ok"] = ok
    v["agrees"] = v.agrees & ok
    out = resolve_drive() / "slc_pilot" / "2021" / "qc"
    out.mkdir(parents=True, exist_ok=True)
    t.to_csv(out / "pairs.csv", index=False)
    v.to_csv(out / "verdict.csv", index=False)
    pl.to_csv(out / "placement.csv", index=False)
    (out / "README.md").write_text(readme(t, v, a.prereg_commit, pl), encoding="utf-8")
    print(pl.to_string(index=False))
    print(v.to_string(index=False))
    print(f"wrote {out}")


def _md(df: pd.DataFrame, nd: int = 3) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(t: pd.DataFrame, v: pd.DataFrame, prereg: str, pl: pd.DataFrame | None = None) -> str:
    cols = ["track", "t1", "t2", "filtered_spread_hard_rad", "expected_spread_hard_rad", "unfiltered_spread_hard_rad",
            "filtered_offset_stable_rad", "coh_hyp3_hard_targets", "coh_ours_hard_targets", "coh_hyp3_mat", "coh_ours_mat"]
    coh = t.groupby("track")[[c for c in t.columns if c.startswith("coh_")]].median().reset_index()
    return f"""# X-101 QC — our CSLC pairs vs HyP3 (generated)

Generated by `experiments/slc_qc_x101.py` on Colab; do not edit. Pre-registration: hub
`08_deliverables/insar_full_test/preregistration/X-101.yaml` (commit {prereg}). Ours: COMPASS CSLCs (5 m x 10 m), 32 looks
per 40 m cube cell, Goldstein–Werner alpha {ALPHA} (an approximation of HyP3's ISCE filter); HyP3: the cube's wrapped
pairs. Spread = circular std over hard targets of (ours − HyP3) after removing its mean there; expected = Cramér–Rao
from OUR unfiltered coherence (v2; v1 used HyP3's, biased high by its filter — shown as `_v1`), L = {L_HYP3} and
{L_OURS}, independent noise (a generous bound: both come from the same SLCs).

## Verdict

{_md(v)}

## Placement against the airborne LiDAR canopy (mean June–August amplitude; correct = opposite east-west signs,
|dx asc + dx desc| <= {PLACE_TOL_M:.0f} m, |dy| <= {PLACE_TOL_M:.0f} m)

{_md(pl) if pl is not None else "—"}

## Coherence, median over pairs (ours unfiltered at 32 looks, HyP3 as delivered)

{_md(coh)}

## Every pair

{_md(t[[c for c in cols if c in t.columns]])}
"""


if __name__ == "__main__":
    main()
