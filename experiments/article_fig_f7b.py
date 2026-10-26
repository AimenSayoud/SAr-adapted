"""F7b — each laser winter left out in turn (X-135): how far the plain chain and the re-levelled chain miss the laser's
step across that winter, per track.

No laser years outside the design years exist (hub DECISIONS 2026-10-26), so Fig. 7's second panel shows X-135's nested
leave-one-winter-out instead of a held-out score. For each winter: the error of the step (estimate − laser, vertical mm)
of the plain chain, and of the chain re-levelled by the weather model as X-119 chose it (its fit leaving the winter's two
years out). Repeated without the winter, X-119's rule chose the plain chain in every fold, so the procedure's error is
the grey points. The laser's own step is written under each winter. The anchor's leave-one-year-out spread is not drawn
here: four of its five fits include the winter's own years, so it is not an interval on these points (it is Fig. 7a's
band). Read from X-135's table only (no cube, no refit).

Writes only to the hub: ``article_plan/figures/final/F7b_nested_winters.*``. No field value in this file.

    make run E=article_fig_f7b
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402

HUB = hub_root()
SRC = HUB / "08_deliverables" / "insar_full_test" / "nested_winter_x135"
OUT = HUB / "08_deliverables" / "article_plan" / "figures" / "final"
COL = {"chain": "#6e7781", "X-119's A": fs.CLASS["mat"]}
LABEL = {"chain": "plain chain (the procedure's choice in every fold)",
         "X-119's A": "re-levelled by the weather model (fit without the winter)"}


def main() -> None:
    n = pd.read_csv(SRC / "nested_steps.csv")
    v = pd.read_csv(SRC / "verdict.csv").iloc[0]
    winters = sorted(n.outer_winter.unique())
    x = np.arange(len(winters))
    with fs.style():
        fig, axes = plt.subplots(1, 2, figsize=fs.size("full", 75), sharey=True,
                                 gridspec_kw={"wspace": 0.06, "left": 0.08, "right": 0.99, "top": 0.84, "bottom": 0.20})
        record = {}
        for ax, track in zip(axes, ("ascending", "descending")):
            ax.set_facecolor("white")
            ax.grid(False)
            ax.axhline(0, color=fs.INK, lw=0.6, zorder=1)
            t = n[n.track == track]
            for role, dx in (("chain", -0.12), ("X-119's A", 0.12)):
                r = t[t.role == role].set_index("outer_winter").reindex(winters)
                err = r.difference_mm.to_numpy()
                ax.plot(x + dx, err, ls="none", marker="o" if role == "chain" else "s", ms=5, color=COL[role],
                        zorder=3, label=LABEL[role] if track == "ascending" else None)
            laser = t[t.role == "chain"].set_index("outer_winter").reindex(winters).laser_change_mm
            ax.set_xticks(x, [f"{w}\nlaser {s:+.0f} mm" for w, s in zip(winters, laser)])
            ax.set_xlim(-0.6, len(winters) - 0.4)
            ax.text(0.02, 0.97, f"{track} track", transform=ax.transAxes, ha="left", va="top", fontsize=8,
                    fontweight="bold")
            record[track] = {"winters": len(winters)}
        axes[0].set_ylabel("step error: estimate − laser (mm, vertical)")
        fig.text(0.01, 0.955, "(b) each laser winter left out in turn", fontsize=8, fontweight="bold", ha="left",
                 va="top")
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.27, 0.995), ncol=2, fontsize=7)
        rec = fs.export(fig, OUT / "F7b_nested_winters", provenance={
            "inputs": ["hub 08_deliverables/insar_full_test/nested_winter_x135/nested_steps.csv (X-135)",
                       "hub 08_deliverables/insar_full_test/nested_winter_x135/verdict.csv"],
            "reading": str(v.reading), "choices": str(v.choices), "laser_rule": str(v.laser_rule),
            "per_track": record})
        plt.close(fig)
    print("F7b →", OUT, rec["width_mm"], "×", rec["height_mm"], "mm")


if __name__ == "__main__":
    main()
