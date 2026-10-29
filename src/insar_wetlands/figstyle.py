"""The article figures' common style and export: Elsevier sizes, one type scale, one palette per meaning, and an export
that writes a vector PDF, a PNG preview and a provenance record side by side.

Why one module. Each figure script used to set its own sizes, fonts and colours, so a quantity could change colour
between figures and a figure could be drawn at a width the journal does not print. Here the widths are the publisher's
(Elsevier artwork sizing: 90, 140, 190 mm), text is at least 7 pt at printed size, and a class (the mat, the
grassland, the hard targets, open water) keeps the same colour in every figure.

    from insar_wetlands import figstyle as fs
    with fs.style():
        fig = plt.figure(figsize=fs.size("full", 120))
        ...
        fs.export(fig, out_dir / "F1_site", provenance={...})
"""
from __future__ import annotations

import contextlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl

MM = 1 / 25.4
WIDTH_MM = {"single": 90.0, "one_half": 140.0, "full": 190.0}        # Elsevier artwork sizing
MIN_FONT_PT = 7.0                                                    # 6 pt allowed only for sub/superscripts
RASTER_DPI = 600                                                     # rasterised layers inside the vector PDF

INK, MUTED, GRID = "#1f2328", "#57606a", "#d0d7de"
# One colour per class, in every figure (Okabe–Ito hues, distinguishable with colour-vision deficiency).
CLASS = {"mat": "#009E73", "grassland": "#E69F00", "hard targets": "#0072B2", "open water": "#56B4E9",
         "laser": "#D55E00", "reference wells": "#57606a"}
# Darker companion hues (>= 3.5:1 contrast on white) for coloured text labels or thin un-outlined strokes.
CLASS_TEXT = {"mat": "#00805D", "grassland": "#B37700", "hard targets": "#0072B2", "open water": "#0277BD",
              "laser": "#D55E00", "reference wells": "#57606a"}
TRACK = {"ascending": "#0072B2", "descending": "#D55E00"}

RC = {
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7.5, "axes.labelsize": 7.5, "axes.titlesize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "legend.fontsize": 7, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 3, "ytick.major.size": 3, "xtick.direction": "out", "ytick.direction": "out",
    "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "lines.linewidth": 1.2,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "pdf.fonttype": 42, "ps.fonttype": 42,                  # TrueType embedded: text stays text, nothing substituted
}


def size(width: str | float, height_mm: float) -> tuple[float, float]:
    """(width, height) in inches for a figure ``width`` wide ('single', 'one_half', 'full' or millimetres)."""
    w = WIDTH_MM[width] if isinstance(width, str) else float(width)
    return w * MM, height_mm * MM


@contextlib.contextmanager
def style():
    with mpl.rc_context(RC):
        yield


def panel_label(ax, letter: str, x: float = -0.02, y: float = 1.02) -> None:
    """A bold panel letter at the axes' top-left corner, outside the data."""
    ax.text(x, y, f"({letter})", transform=ax.transAxes, fontsize=8, fontweight="bold", ha="right", va="bottom")


def _git_sha(path: Path) -> str:
    try:
        return subprocess.run(["git", "-C", str(path), "rev-parse", "--short", "HEAD"], capture_output=True,
                              text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def export(fig, stem: Path, provenance: dict, png_dpi: int = 600) -> dict:
    """Write ``stem``.pdf (vector, fonts embedded, rasters at their own resolution), ``stem``.png (a preview at
    ``png_dpi``) and ``stem``.json (what the figure was made from), at the figure's own page size — never cropped,
    so the printed width is the width it was drawn at. Returns the provenance record."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    # dpi applies only to rasterised artists (dense scatters) and resampled images; `interpolation="none"` images keep
    # their native resolution. 600 dpi meets Elsevier's 500 dpi for combination art.
    fig.savefig(stem.with_suffix(".pdf"), dpi=RASTER_DPI, metadata={"Creator": "insar_wetlands.figstyle"})
    fig.savefig(stem.with_suffix(".png"), dpi=png_dpi)
    w, h = fig.get_size_inches()
    rec = {"figure": stem.name, "width_mm": round(w / MM, 1), "height_mm": round(h / MM, 1),
           "made": datetime.now(timezone.utc).isoformat(timespec="seconds"),
           "code_sha": _git_sha(Path(__file__).parent), **provenance}
    stem.with_suffix(".json").write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
    return rec
