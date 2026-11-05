"""Build public presentation plots from bundled GIS and the synthetic demo only.

No station data, network access or private workspace is used. The conceptual
overview is an editable SVG maintained directly in docs/assets/.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
from examples.synthetic_scoring import synthetic_series  # noqa: E402
from insar_wetlands.scoring import increment_stats  # noqa: E402

ASSETS = ROOT / "docs/assets"
INK, MUTED, GREEN, AMBER = "#193c36", "#64766c", "#2f7d66", "#b67b28"
BG = "#f8faf7"


def setup() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "text.color": INK,
        "axes.labelcolor": INK, "axes.edgecolor": "#cdd8d0", "axes.titlecolor": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "figure.facecolor": BG,
        "axes.facecolor": BG, "savefig.facecolor": BG, "axes.spines.top": False,
        "axes.spines.right": False, "grid.color": "#dee7df", "grid.linewidth": 0.7,
    })


def synthetic_plot() -> None:
    epochs, truth = synthetic_series()
    days = (epochs - epochs[0]).days.to_numpy()
    increments = np.diff(truth)
    half = increment_stats(increments, increments / 2)
    perfect = increment_stats(increments, increments)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={"width_ratios": [1.3, 1]})
    fig.subplots_adjust(left=0.075, right=0.96, top=0.72, bottom=0.27, wspace=0.26)
    fig.text(0.075, 0.92, "SYNTHETIC SCORING / NO FIELD OBSERVATIONS", size=10, weight="bold", color=GREEN)
    fig.text(0.075, 0.83, "A shared shape does not guarantee the right amplitude.", size=20, weight="bold")
    ax = axes[0]
    ax.plot(days, truth, color=GREEN, lw=2.8, label="Synthetic truth / perfect estimate")
    ax.plot(days, truth / 2, color=AMBER, lw=2.4, ls="--", label="Half-sized estimate")
    ax.axhline(0, color="#a3b3a8", lw=0.7)
    ax.set(xlabel="Days since synthetic start", ylabel="Displacement (mm)")
    ax.grid(axis="y")
    fig.legend(*ax.get_legend_handles_labels(), loc="lower left", bbox_to_anchor=(0.065, 0.07),
               frameon=False, fontsize=9)
    ax = axes[1]
    limit = float(np.max(np.abs(increments)) * 1.18)
    ax.plot([-limit, limit], [-limit, limit], color=GREEN, lw=1.5, label="Perfect estimate")
    ax.scatter(increments, increments / 2, color=AMBER, s=25, alpha=0.9, zorder=3, label="Half-sized estimate")
    ax.set(xlim=(-limit, limit), ylim=(-limit, limit), xlabel="Truth increment (mm)", ylabel="Estimated increment (mm)")
    ax.set_aspect("equal", adjustable="box")
    ax.grid()
    ax.text(0.04, 0.96, f"Half-sized slope: {half['slope_through_origin']:.3f}\n"
            f"RMSE: {half['rmse_increments_mm']:.3f} mm", transform=ax.transAxes, va="top", fontsize=9,
            bbox={"boxstyle": "round,pad=0.5", "facecolor": BG, "edgecolor": "none"})
    fig.text(0.075, 0.025, f"Deterministic demo · {perfect['pairs']} increments · "
             f"perfect-estimate RMSE = {perfect['rmse_increments_mm']:.3f} mm · not measured Rzecin performance",
             size=9, color=MUTED)
    fig.savefig(ASSETS / "synthetic-scoring.png", dpi=180, metadata={"Description": "Synthetic demo; no field observations."})
    plt.close(fig)


def gis_plot() -> None:
    gis = ROOT / "data/gis"
    with rasterio.open(gis / "zones_classified_raster.tif") as src:
        zones = src.read(1)
        dx, dy = abs(src.transform.a) / 1000, abs(src.transform.e) / 1000
        crs, transform = src.crs, src.transform
    with rasterio.open(gis / "mean_coherence_spatial.tif") as src:
        coherence = np.ma.masked_invalid(src.read(1))
        if src.crs != crs or src.transform != transform or coherence.shape != zones.shape:
            raise ValueError("Bundled GIS layers must share their grid; do not resample presentation data.")
    row, col = np.where((zones == 1) | (zones == 2))
    if not len(row):
        raise ValueError("Bundled zone raster contains no mat/lake region.")
    # Crop around the mat/lake, keeping surrounding comparison context; values are unchanged.
    r0, r1 = max(0, int(row.min()) - 15), min(zones.shape[0], int(row.max()) + 16)
    c0, c1 = max(0, int(col.min()) - 15), min(zones.shape[1], int(col.max()) + 16)
    zones, coherence = zones[r0:r1, c0:c1], coherence[r0:r1, c0:c1]
    extent = [0, zones.shape[1] * dx, 0, zones.shape[0] * dy]
    fig, axes = plt.subplots(1, 2, figsize=(12, 6.5))
    fig.subplots_adjust(left=0.07, right=0.96, bottom=0.21, top=0.72, wspace=0.2)
    fig.text(0.07, 0.94, "PUBLIC GIS PREVIEW / FIRST-STUDY CONTEXT", size=10, weight="bold", color=GREEN)
    fig.text(0.07, 0.865, "One landscape. Two views of the same grid.", size=22, weight="bold")
    fig.text(0.07, 0.81, "Rzecin · bundled 40 m layers · cropped for display, with no resampling", size=11, color=MUTED)
    colors = ["#edf0eb", "#4e8561", "#7eb6c5", "#dfbf71", "#d6dcd1"]
    axes[0].imshow(zones, cmap=ListedColormap(colors), norm=BoundaryNorm(np.arange(-0.5, 5.5), 5),
                   extent=extent, interpolation="nearest", origin="upper")
    image = axes[1].imshow(coherence, cmap="cividis", vmin=0.3, vmax=0.8, extent=extent,
                           interpolation="nearest", origin="upper")
    axes[0].set_title("Comparison zones", loc="left", fontsize=13, weight="bold", pad=12)
    axes[1].set_title("Mean spatial coherence", loc="left", fontsize=13, weight="bold", pad=12)
    for i, ax in enumerate(axes):
        foreground = INK if i == 0 else "#ffffff"
        ax.set(xlabel="Distance east from crop edge (km)", ylabel="Distance north from crop edge (km)")
        ax.annotate("N", xy=(0.94, 0.94), xytext=(0.94, 0.77), xycoords="axes fraction",
                    ha="center", color=foreground, arrowprops={"arrowstyle": "-|>", "color": foreground})
        x, y = 0.15, 0.17
        ax.plot([x, x + 0.5], [y, y], color=foreground, lw=3)
        ax.text(x + 0.25, y + 0.06, "500 m", ha="center", fontsize=8, color=foreground)
    labels = ["Mat (A)", "Lake (B)", "Matched grassland (C)", "Other comparison area (D)"]
    fig.legend([Patch(facecolor=color) for color in colors[1:]], labels, loc="lower left",
               bbox_to_anchor=(0.06, 0.095), ncol=4, frameon=False, fontsize=9)
    bar = fig.colorbar(image, ax=axes[1], fraction=0.046, pad=0.035)
    bar.set_label("Mean coherence")
    fig.text(0.07, 0.035, "Source: data/gis · first-study 2022–2024 products · descriptive context, not displacement or validation",
             size=9, color=MUTED)
    fig.savefig(ASSETS / "public-gis-preview.png", dpi=180,
                metadata={"Description": "Bundled first-study public GIS layers; no station records or station locations."})
    plt.close(fig)


if __name__ == "__main__":
    ASSETS.mkdir(parents=True, exist_ok=True)
    setup()
    synthetic_plot()
    gis_plot()
    print("Built synthetic-scoring.png and public-gis-preview.png using public/synthetic inputs only.")
