"""F6, final version — the systematic test as a verdict matrix: every InSAR route tried, grouped by what it varies, the
product it ran on and its outcome against the laser at P6 (or its own criterion).

Rows and outcomes come from ``article_plan/f6_routes.yaml``; every row's tries must exist in the tries map
(``independent_review/01_tries/tries.csv``), as for Table 1. The tries stay in the provenance record, not on the figure (the article carries no internal codes). The outcome is drawn twice — a symbol and a colour — so the
matrix reads in grey and with colour-vision deficiency. Writes only to the hub
``08_deliverables/article_plan/figures/final/``, and Table 1 (``article_plan/tables/table1_systematic_test.csv``) from
the same rows, with the radar years of each row's tries.

    make run E=article_fig_f6
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402

HUB = hub_root()
PLAN = HUB / "08_deliverables" / "article_plan"
OUT = PLAN / "figures" / "final"
TRIES = HUB / "08_deliverables" / "independent_review" / "01_tries" / "tries.csv"
# outcome → (symbol, fill, edge, words)
OUTCOME = {"holds": ("o", fs.CLASS["mat"], fs.CLASS["mat"], "follows the motion"),
           "helps": ("^", fs.CLASS["mat"], fs.CLASS["mat"], "helps (design years)"),
           "partly": ("o", "white", fs.CLASS["mat"], "partly: ascending only"),
           "mixed": ("X", "white", fs.CLASS["laser"], "one track better, one worse: fails by its rule"),
           "lost": ("X", fs.CLASS["laser"], fs.CLASS["laser"], "loses it"),
           "not_motion": ("D", "#CC79A7", "#CC79A7", "a signal, not the motion"),
           "not_applicable": ("o", "white", fs.MUTED, "fails its own conditions here"),
           "no_effect": ("_", fs.MUTED, fs.MUTED, "no gain, no change"),
           "inconclusive": ("s", "white", fs.MUTED, "inconclusive"),
           "not_tried": (None, None, None, "not tried here")}
PRODUCTS = (("HyP3", "HyP3\n40 m"), ("SLC", "own SLC\n2021"), ("sim", "simulation"))
PRODUCT_WORDS = {"HyP3": "HyP3", "SLC": "own SLC", "both": "HyP3 + own SLC", "sim": "simulation", "none": "—"}


def table1() -> None:
    """Table 1: the figure's rows, each with its product, its outcome in words and the radar years of its tries."""
    spec = yaml.safe_load((PLAN / "f6_routes.yaml").read_text(encoding="utf-8"))
    tries = pd.read_csv(TRIES).set_index("id")
    pd.DataFrame([{"group": g["name"], "route": r["route"], "product": PRODUCT_WORDS[r["product"]],
                   "outcome": r.get("words") or OUTCOME[r["outcome"]][3],
                   "years": "; ".join(sorted({str(tries.loc[i, "period"]) for i in r.get("tries", [])
                                              if i in tries.index})) if r["outcome"] != "not_tried" else "",
                   "tries": ", ".join(r.get("tries", []))}
                  for g in spec["groups"] for r in g["rows"]]).to_csv(PLAN / "tables" / "table1_systematic_test.csv",
                                                                      index=False)


def rows() -> list[tuple]:
    spec = yaml.safe_load((PLAN / "f6_routes.yaml").read_text(encoding="utf-8"))
    known = set(pd.read_csv(TRIES)["id"])
    out = []
    for g in spec["groups"]:
        out.append(("group", g["name"], None, None, None, None))
        for r in g["rows"]:
            missing = set(r.get("tries", [])) - known
            if missing:
                raise SystemExit(f"F6: {r['route']} cites unknown tries {missing}")
            out.append(("row", r["route"], r["product"], r["outcome"], ", ".join(r.get("tries", [])), r.get("words")))
    return out


def main() -> None:
    rr = rows()
    n = len(rr)
    with fs.style():
        fig = plt.figure(figsize=fs.size("full", 7 + 4.15 * n))
        ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
        ax.set_xlim(0, 190), ax.set_ylim(n + 2.2, -1.25)
        ax.axis("off")
        xp = {"HyP3": 91, "SLC": 104, "sim": 117}
        xo, xw, xt = 128, 130, 188
        for p, lab in PRODUCTS:
            ax.text(xp[p], -0.35, lab, ha="center", va="center", fontsize=7, color=fs.MUTED, linespacing=1.1)
        ax.text(xo, -0.35, "outcome: against the laser at P6,\nor by the route's own test", ha="left", va="center",
                fontsize=7, color=fs.MUTED, linespacing=1.1)
        y = 0.6
        for kind, name, prod, out, tries, words_override in rr:
            if kind == "group":
                y += 0.35
                ax.text(2, y, name, ha="left", va="center", fontsize=7.5, fontweight="bold")
                ax.plot([2, 188], [y + 0.5, y + 0.5], color=fs.GRID, lw=0.5)
                y += 1.0
                continue
            tried = out != "not_tried"
            ax.text(5, y, name, ha="left", va="center", fontsize=7, color=fs.INK if tried else fs.MUTED)
            for p in xp:
                ax.plot(xp[p], y, "o", ms=2.2, color=fs.GRID)
            used = () if not tried else ("HyP3", "SLC") if prod == "both" else (prod,)
            mk, fc, ec, words = OUTCOME[out]
            for p in used:
                ax.plot(xp[p], y, mk, ms=6.2 if mk != "_" else 8, mfc=fc, mec=ec, mew=1.1, color=ec)
            ax.text(xw, y, words_override or words, ha="left", va="center", fontsize=7,
                    color=fs.INK if tried else fs.MUTED, style="normal" if tried else "italic")
            y += 1.0
        ax.set_ylim(y + 0.3, -1.25)
        rec = fs.export(fig, OUT / "F6_systematic_test", provenance={
            "inputs": ["hub 08_deliverables/article_plan/f6_routes.yaml",
                       "hub 08_deliverables/independent_review/01_tries/tries.csv (every row's tries checked)"],
            "rows": sum(1 for r in rr if r[0] == "row"),
            "tries_by_route": {r[1]: r[4] for r in rr if r[0] == "row"},
            "brief": "08_deliverables/article_plan/figures/briefs/F6.md"})
        plt.close(fig)
    table1()
    print("F6 →", OUT, rec["width_mm"], "×", rec["height_mm"], "mm")


if __name__ == "__main__":
    main()
