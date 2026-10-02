"""X-064 — the lake control with a pure open-water mask.

H4 argues "the lake oscillates too (B − C), and open water cannot breathe". X-063 found that at
40 m the lake pixels are about as coherent as the mat and not dark like open water: a mixed
target. This builds open-water masks from the native 10 m RTC (``rtc_10m_window.py``, summer
2022–2024, both tracks: a 10 m cell is water when its median σ0 VV and VH are below the
thresholds on **both** tracks; a 40 m cell's water share is the share of its 16 sub-cells), then
puts each mask through the same chain as the existing lake series — coherence-weighted
unwrapped phase, mask minus grassland C per pair (``aggregate_unwrapped``), the network inverted
as one super-pixel (``invert_aggregate``), the seasonal model (``seasonal_amplitude``), as X-054 —
on both tracks and both stacks (2020–2021, 2022–2024).

Masks, from the current to the strictest: zone B; phaseDter's clean lake (X-063); 40 m cells
≥ 75 % open water; 40 m cells 100 % open water. A small mask has a larger noise floor: the
amplitudes here are descriptive (95 % intervals from the fit). The size-matched null the project
requires for an inferential p is a Colab job (never run locally) and is listed, not run.

No field data; outputs in the hub next to X-063.

    INSAR_DRIVE_ROOT=../local/drive_pristine PYTHONPATH=src python scripts/lake_open_water_x064.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from field_series_x054 import EPOCH, amp_ci, invert  # noqa: E402

from insar_wetlands.aggregate import aggregate_unwrapped, seasonal_amplitude  # noqa: E402
from insar_wetlands.bootstrap import start  # noqa: E402
from insar_wetlands.paths import make_paths  # noqa: E402
from insar_wetlands.stack import list_pairs, load_layer  # noqa: E402
from insar_wetlands.stratify import clean_lake_mask, load_worldcover  # noqa: E402

HUB = Path(__file__).resolve().parents[3]
VV_MAX_DB = -15.0          # summer-median σ0 VV below this on both tracks …
VH_MAX_DB = -20.0          # … and VH below this: open water (calm-to-light wind, no emergent vegetation)
SHARE_MIXED = 0.75
MASKS = ["zone B (current lake)", "clean lake (X-063)", "open water ≥ 75 %", "open water 100 %"]


def water_share(rtc_root: Path, shape40: tuple[int, int]) -> tuple[np.ndarray, dict]:
    med = {}
    w = None
    for t in ("ascending", "descending"):
        ds = xr.open_dataset(rtc_root / f"rtc10m_summer_2022_2024_{t}.nc")
        vv, vh = np.nanmedian(ds.gamma0_vv_db.values, 0), np.nanmedian(ds.gamma0_vh_db.values, 0)
        med[t] = (vv, vh)
        wt = (vv < VV_MAX_DB) & (vh < VH_MAX_DB)
        w = wt if w is None else w & wt
    H, W = shape40
    return w.reshape(H, 4, W, 4).mean((1, 3)), med


def masks(ctx, share: np.ndarray) -> dict[str, xr.DataArray]:
    tpl = ctx.template
    B = ctx.zones["B"]
    wc = load_worldcover(tpl, ctx.cfg, cache_dir=ctx.paths.cache)
    clean = clean_lake_mask(tpl, ctx.cfg, worldcover=wc, s2=xr.load_dataset(ctx.paths.drive_file("s2_stack.nc")))
    as_da = lambda m: xr.DataArray(m, coords=B.coords, dims=B.dims)  # noqa: E731
    return {MASKS[0]: B.astype(bool), MASKS[1]: as_da(np.asarray(clean, bool)),
            MASKS[2]: as_da(share >= SHARE_MIXED), MASKS[3]: as_da(share >= 1.0)}


def stacks(ctx) -> list[tuple[str, str, Path]]:
    out = []
    for track in ("ascending", "descending"):
        out.append((track, "2022–2024", make_paths(cfg=ctx.cfg, track=track).cropped))
        ext = HUB / "05_code" / "local" / "s1_2020_2021" / f"hyp3_cropped_{track}"
        if ext.exists():
            out.append((track, "2020–2021", ext))
    return out


def run(ctx, mk: dict, rtc_med: dict, share: np.ndarray):
    C = ctx.zones["C"]
    A = ctx.zones["A"]
    desc, pairs_rows, fits = [], [], []
    for track, stack, root in stacks(ctx):
        pairs = list_pairs(root)
        unw, corr = load_layer(root, "unw_phase", pairs), load_layer(root, "corr", pairs)
        dt = np.array([(pd.Timestamp(p.split("_")[1]) - pd.Timestamp(p.split("_")[0])).days for p in pairs])
        short = dt <= 12
        for name, m in {**mk, "mat (zone A)": A.astype(bool), "grassland (zone C)": C.astype(bool)}.items():
            mv = m.values
            valid = np.isfinite(unw.values[:, mv]).sum(1) if mv.any() else np.zeros(len(pairs))
            coh = np.nanmedian(corr.values[:, mv], 1) if mv.any() else np.full(len(pairs), np.nan)
            row = {"track": track, "stack": stack, "mask": name, "n_px": int(mv.sum()),
                   "water_share_median": float(np.median(share[mv])) if mv.any() else np.nan,
                   "median_valid_px_per_pair": float(np.median(valid)),
                   "coh_short_pairs_median": float(np.nanmedian(coh[short])) if short.any() else np.nan,
                   "coh_24_48d_median": float(np.nanmedian(coh[(dt > 12) & (dt <= 48)])) if ((dt > 12) & (dt <= 48)).any() else np.nan,
                   "coh_floor_gt96d_median": float(np.nanmedian(coh[dt > 96])) if (dt > 96).any() else np.nan}
            if stack == "2022–2024" and mv.any():
                vv10, vh10 = rtc_med[track]
                big = np.kron(mv, np.ones((4, 4), bool))
                row.update({"vv10m_summer_median_db": float(np.nanmedian(vv10[big])),
                            "vh10m_summer_median_db": float(np.nanmedian(vh10[big]))})
            desc.append(row)
            if name not in mk or mv.sum() < 5:
                continue
            for label, z, tgt, ref in ((f"{name} − C", {"T": m, "C": C}, "T", "C"),
                                       (f"A − {name}", {"A": A, "T": m}, "A", "T")):
                dd = aggregate_unwrapped(unw, corr, z, target=tgt, reference=ref)
                pairs_rows.append(dd.assign(track=track, stack=stack, mask=name, zones=label))
                if len(dd) < 20:
                    continue
                s = invert(dd)
                f = seasonal_amplitude(s, epoch=EPOCH)
                lo, hi = amp_ci(f)
                fits.append({"track": track, "stack": stack, "mask": name, "zones": label, "n_pairs": len(dd),
                             "n_dates": len(s), "amplitude_mm": f["amplitude_mm"], "ci95_low": lo, "ci95_high": hi,
                             "phase_doy": f["phase_doy"], "r2_seasonal": f["r2_seasonal"]})
        # A − C on the same stack, the reference the lake is compared with
        dd = aggregate_unwrapped(unw, corr, {"A": A, "C": C}, target="A", reference="C")
        s = invert(dd)
        f = seasonal_amplitude(s, epoch=EPOCH)
        lo, hi = amp_ci(f)
        fits.append({"track": track, "stack": stack, "mask": "mat (zone A)", "zones": "A − C", "n_pairs": len(dd),
                     "n_dates": len(s), "amplitude_mm": f["amplitude_mm"], "ci95_low": lo, "ci95_high": hi,
                     "phase_doy": f["phase_doy"], "r2_seasonal": f["r2_seasonal"]})
        del unw, corr
    return pd.DataFrame(desc), pd.concat(pairs_rows, ignore_index=True), pd.DataFrame(fits)


def check_t07(fits: pd.DataFrame) -> pd.DataFrame:
    """Zone B and the mat through this chain, ascending 2022–2024, must reproduce T07 (the committed amplitudes)."""
    t07 = pd.read_csv(Path(__file__).resolve().parents[1] / "results" / "tables" / "T07_seasonal_amplitudes.csv").set_index("series")
    g = fits[(fits.track == "ascending") & (fits["stack"] == "2022–2024")].set_index("zones")
    rows = [("B−C", "zone B (current lake) − C"), ("A−C", "A − C"), ("A−B", "A − zone B (current lake)")]
    return pd.DataFrame([{"series": k, "T07_amplitude_mm": float(t07.loc[k, "amplitude_mm"]),
                          "here_amplitude_mm": float(g.loc[z, "amplitude_mm"]),
                          "abs_diff_mm": abs(float(t07.loc[k, "amplitude_mm"]) - float(g.loc[z, "amplitude_mm"]))} for k, z in rows])


def fig(out: Path, ctx, mk, share, desc, fits):
    fig, ax = plt.subplots(1, 3, figsize=(17, 5.2))
    B = ctx.zones["B"].values.astype(bool)
    rows, cols = np.nonzero(B | (share > 0))
    r0, r1, c0, c1 = rows.min() - 3, rows.max() + 4, cols.min() - 3, cols.max() + 4
    im = ax[0].imshow(share[r0:r1, c0:c1], cmap="Blues", vmin=0, vmax=1)
    for name, col in ((MASKS[0], "orange"), (MASKS[1], "magenta"), (MASKS[3], "red")):
        ax[0].contour(mk[name].values[r0:r1, c0:c1].astype(float), levels=[0.5], colors=col, linewidths=1.2)
    ax[0].set_title("open-water share of each 40 m cell (10 m RTC)\norange: zone B · magenta: clean lake · red: 100 % water", fontsize=8)
    plt.colorbar(im, ax=ax[0], shrink=0.8)
    d = desc[desc["stack"] == "2022–2024"]
    for k, (t, c) in enumerate((("ascending", "tab:orange"), ("descending", "tab:blue"))):
        g = d[d.track == t].set_index("mask").reindex(MASKS)
        ax[1].bar(np.arange(len(MASKS)) + (k - 0.5) * 0.35, g.coh_short_pairs_median, 0.35, color=c, label=t)
    ax[1].set_xticks(range(len(MASKS)), [m.replace(" (", "\n(") for m in MASKS], fontsize=7)
    ax[1].set_ylabel("median coherence, pairs ≤ 12 d (2022–2024)")
    ax[1].legend(fontsize=7)
    f = fits[fits.zones.str.endswith("− C")]
    labels = list(dict.fromkeys(f["mask"]))
    for k, (t, c) in enumerate((("ascending", "tab:orange"), ("descending", "tab:blue"))):
        for j, st in enumerate(("2020–2021", "2022–2024")):
            g = f[(f.track == t) & (f["stack"] == st)].set_index("mask").reindex(labels)
            x = np.arange(len(labels)) + (2 * k + j - 1.5) * 0.18
            ax[2].errorbar(x, g.amplitude_mm, yerr=[g.amplitude_mm - g.ci95_low, g.ci95_high - g.amplitude_mm], fmt="o" if j else "s",
                           color=c, ms=4, lw=0.8, label=f"{t} {st}")
    ax[2].set_xticks(range(len(labels)), [m.replace(" (", "\n(") for m in labels], fontsize=7)
    ax[2].set_ylabel("seasonal amplitude vs grassland C (mm, 95 % CI)")
    ax[2].legend(fontsize=6)
    fig.tight_layout()
    fig.savefig(out / "fig_lake_open_water.png", dpi=110)
    plt.close(fig)


def md(df: pd.DataFrame) -> str:
    lines = ["| " + " | ".join(df.columns) + " |", "|" + "---|" * len(df.columns)]
    for r in df.itertuples(index=False):
        lines.append("| " + " | ".join(f"{v:.2f}" if isinstance(v, float) else str(v) for v in r) + " |")
    return "\n".join(lines)


def readme(out: Path, desc, fits, chk):
    d = desc.drop(columns=[c for c in ("water_share_median",) if c in desc.columns])
    find = []
    for t in ("ascending", "descending"):
        g = desc[(desc.track == t) & (desc["stack"] == "2022–2024")].set_index("mask")
        find.append(f"- **{t}**: short-pair coherence {g.loc[MASKS[0], 'coh_short_pairs_median']:.2f} (zone B), "
                    f"{g.loc[MASKS[1], 'coh_short_pairs_median']:.2f} (clean lake), {g.loc[MASKS[2], 'coh_short_pairs_median']:.2f} "
                    f"(≥ 75 % water), {g.loc[MASKS[3], 'coh_short_pairs_median']:.2f} (100 % water, {int(g.loc[MASKS[3], 'n_px'])} cells).")
    for t in ("ascending", "descending"):
        for st in ("2020–2021", "2022–2024"):
            g = fits[(fits.track == t) & (fits["stack"] == st)].set_index("zones")
            parts = [f"{z.replace(' − C', '')} {g.loc[z, 'amplitude_mm']:.2f} [{g.loc[z, 'ci95_low']:.2f}, {g.loc[z, 'ci95_high']:.2f}]"
                     for z in g.index if z.endswith("− C")]
            find.append(f"- *{t} {st}* seasonal amplitude vs C, mm [95 %]: " + "; ".join(parts) + ".")
    L = ["# X-064 — the lake control with a pure open-water mask (exploratory)", "",
         "Generated by `05_code/SAr-adapted/scripts/lake_open_water_x064.py`; every number is computed by it. Why: H4 uses "
         "\"the lake oscillates too (B − C); open water cannot breathe\", and X-063 found the 40 m lake pixels about as "
         "coherent as the mat — a mixed target.", "",
         "## Findings", "", *find, "",
         "## Masks", "",
         f"Open water at 10 m: summer-median σ0 VV < {VV_MAX_DB:.0f} dB **and** VH < {VH_MAX_DB:.0f} dB on both tracks "
         f"(RTC 10 m, May–September 2022–2024). A 40 m cell's water share = share of its 16 sub-cells. Masks: zone B; "
         f"phaseDter's clean lake (X-063); cells ≥ {SHARE_MIXED:.0%} water; cells 100 % water.", "",
         md(d), "",
         "`median_valid_px_per_pair`: mask pixels with an unwrapped phase in a typical pair. `coh_floor_gt96d_median`: "
         "coherence of pairs longer than 96 days, where every surface here has decorrelated — the estimator's floor in "
         "these products; a target is coherent only to the extent it sits above it (2020–2021 has no such pairs).", "",
         "## Seasonal amplitude, mask − grassland C and mat − mask (same chain as X-054)", "",
         md(fits), "",
         "Reproduction: zone B and the mat through this chain (ascending, 2022–2024) against the committed T07:", "", md(chk), "",
         "## How to read this", "",
         "- The mask − C amplitude and its timing barely change from zone B to the purest open water, and the same timing "
         "appears in A − C: a cycle common to every target measured against C sits at least partly in the grassland "
         "reference itself (as X-050 found for the water-linked phase), not in a lake that breathes.",
         "- Open water should be dark and incoherent. Where the strict mask keeps low coherence, its unwrapped phase is "
         "mostly noise: an \"oscillation\" there would be noise or unwrapping, not a dielectric cycle.",
         "- If the seasonal B − C amplitude shrinks or loses its interval as the mask gets purer, the lake signal came from "
         "the vegetated and shore pixels — then \"the lake oscillates too\" does not show that open water oscillates.",
         "- Small masks have larger noise floors; the size-matched null needed for a p-value is a Colab job (not run here).", "",
         "## Files", "", "| File | Content |", "|---|---|",
         "| `masks_summary.csv` | per mask × track × stack: size, water share, valid pixels, coherence, 10 m σ0 |",
         "| `zone_pairs.csv` | per pair: mask − C and A − mask aggregated phase |",
         "| `seasonal_fits.csv` | amplitudes and intervals |",
         "| `open_water_share_40m.npy` | water share of every 40 m cell (template grid) |",
         "| `fig_lake_open_water.png` | the masks, coherence, amplitudes |"]
    (out / "README.md").write_text("\n".join(L) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(HUB / "08_deliverables" / "lake_x064"))
    ap.add_argument("--rtc", default=str(HUB / "05_code" / "local" / "drive_mirror"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ctx = start("lake_open_water_x064", mount=False, git=False)
    share, rtc_med = water_share(Path(a.rtc), ctx.template.shape)
    np.save(out / "open_water_share_40m.npy", share)
    mk = masks(ctx, share)
    desc, pairs, fits = run(ctx, mk, rtc_med, share)
    desc.to_csv(out / "masks_summary.csv", index=False)
    pairs.to_csv(out / "zone_pairs.csv", index=False)
    fits.to_csv(out / "seasonal_fits.csv", index=False)
    chk = check_t07(fits)
    fig(out, ctx, mk, share, desc, fits)
    readme(out, desc, fits, chk)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
