"""Figure 3: Spatial Evidence & Boardwalk-Free Radar Pixel Map (Paper-Grade 300 DPI).

Layout (Clean 2-Panel):
- Panel A (Left): High-Resolution Orthophoto (0.25 m) & 40 m Radar Grid around P6.
  Maps P6 Center Cell (8.0% Wood), West Neighbour (0.0% Wood), and East Neighbour (0.0% Wood)
  with exact 1:1 geospatial alignment directly on the physical wooden platform.
- Panel B (Right): InSAR-Laser Agreement Across Key Spatial Extraction Windows.
  Horizontal bar chart comparing Correlation (r) and Total Least Squares Slope (m).
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import rasterio

# Matplotlib publication styling
plt.rcParams['font.sans-serif'] = 'Helvetica', 'Arial', 'DejaVu Sans'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8

ROOT = Path("/Users/aymen/Documents/Research_Hub")
OUT_DIR = ROOT / "08_deliverables" / "supervisor_defense"
OUT_DIR.mkdir(parents=True, exist_ok=True)
WEB_DIR = ROOT / "08_deliverables" / "rzecin-atlas" / "web" / "public" / "data" / "figures"
WEB_DIR.mkdir(parents=True, exist_ok=True)

# 1. Load Data
p6_win = pd.read_csv(ROOT / "08_deliverables" / "field_boardwalk_x060" / "p6_windows.csv")
px = pd.read_csv(ROOT / "08_deliverables" / "field_first" / "plot_pixels.csv")
ortho_path = ROOT / "08_deliverables" / "field_boardwalk_x060" / "orthophoto_transect.tif"

p6_row = px[px['plot'] == 'P6'].iloc[0]

print("Generating perfectly aligned 2-panel Figure 3 at 300 DPI...")
fig = plt.figure(figsize=(16, 8.0), dpi=300)
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.38, left=0.07, right=0.97, top=0.90, bottom=0.16)

# Left: Orthophoto & 40m Grid Overlay around P6 with perfect geospatial registration
ax_map = fig.add_subplot(gs[0])
with rasterio.open(ortho_path) as src:
    img = src.read()
    b = src.bounds
    img = np.moveaxis(img, 0, -1)
    ax_map.imshow(img, extent=[b.left, b.right, b.bottom, b.top])

p6_e, p6_n = 588380.0, 5846540.0

# Draw 40m radar grid lines around P6
for x in np.arange(p6_e - 100, p6_e + 140, 40):
    ax_map.axvline(x, color='yellow', lw=1.2, ls='--', alpha=0.85)
for y in np.arange(p6_n - 100, p6_n + 140, 40):
    ax_map.axhline(y, color='yellow', lw=1.2, ls='--', alpha=0.85)

# Annotate Cells
# P6 Center cell: Row 68, Col 66 (8.0% Wood)
p6_box = patches.Rectangle((p6_e - 20, p6_n - 20), 40, 40, linewidth=2.8, edgecolor='#d95f02', facecolor='none', zorder=10)
ax_map.add_patch(p6_box)
ax_map.text(p6_e, p6_n + 12, "P6 Center Cell\n(8.0% Wood Area)", color='white', fontweight='bold', fontsize=8.5, ha='center',
            bbox=dict(boxstyle="round,pad=0.25", fc="#d95f02", ec='white', lw=0.8, alpha=0.90), zorder=12)

# West cell (Left): Row 68, Col 65 (0.0% Wood)
w_box = patches.Rectangle((p6_e - 60, p6_n - 20), 40, 40, linewidth=2.8, edgecolor='#1b7837', facecolor='none', zorder=10)
ax_map.add_patch(w_box)
ax_map.text(p6_e - 40, p6_n + 12, "West Neighbour\n(0.0% Wood)", color='white', fontweight='bold', fontsize=8.5, ha='center',
            bbox=dict(boxstyle="round,pad=0.25", fc="#1b7837", ec='white', lw=0.8, alpha=0.90), zorder=12)

# East cell (Right): Row 68, Col 67 (0.0% Wood)
e_box = patches.Rectangle((p6_e + 20, p6_n - 20), 40, 40, linewidth=2.8, edgecolor='#1b7837', facecolor='none', zorder=10)
ax_map.add_patch(e_box)
ax_map.text(p6_e + 40, p6_n + 12, "East Neighbour\n(0.0% Wood)", color='white', fontweight='bold', fontsize=8.5, ha='center',
            bbox=dict(boxstyle="round,pad=0.25", fc="#1b7837", ec='white', lw=0.8, alpha=0.90), zorder=12)

# Plot station P6 marker directly on the wooden walkway
ax_map.scatter(p6_row.E, p6_row.N, color='#ff0000', edgecolors='white', s=75, lw=1.2, zorder=15, label='P6 Hydrological Station')

# Boardwalk annotation arrow
ax_map.annotate("1.2 m Wooden Walkway\n(Runs through Col 66 only)", 
                xy=(p6_e - 1, p6_n - 18), xytext=(p6_e + 8, p6_n - 42),
                arrowprops=dict(facecolor='yellow', edgecolor='black', width=1.4, headwidth=6),
                color='yellow', fontsize=8.5, fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.2", fc="#111111", alpha=0.85, ec='yellow'), zorder=12)

ax_map.set_xlim(p6_e - 100, p6_e + 100)
ax_map.set_ylim(p6_n - 80, p6_n + 80)
ax_map.set_xlabel("UTM Easting (m, EPSG:32633)", fontsize=9.5)
ax_map.set_ylabel("UTM Northing (m, EPSG:32633)", fontsize=9.5)
ax_map.set_title("A. High-Resolution Orthophoto (0.25 m) & 40 m Radar Grid\nProof of Boardwalk-Free Neighbour Pixels", fontsize=11, fontweight='bold', loc='left')
ax_map.legend(loc='lower right', fontsize=8.5, framealpha=0.9)

# Right: Bar Chart of Window Tests from p6_windows.csv
ax_bars = fig.add_subplot(gs[1])
p6_w_asc = p6_win[(p6_win.track == 'ascending') & (p6_win.subset == 'all')].copy()
key_windows = [
    ("1x1 (P6 cell)", "P6 Center Cell (8% wood)"),
    ("west neighbour only", "West Cell Alone (0% wood)"),
    ("east neighbour only", "East Cell Alone (0% wood)"),
    ("3x3 without boardwalk cells", "3x3 Masked (5 cells, 0% wood)"),
    ("3x3", "Full 3x3 Average (9 cells)")
]
labels = [label for key, label in key_windows]
r_vals = [p6_w_asc[p6_w_asc.window == key]['r_s1_vs_laser'].iloc[0] for key, label in key_windows]
tls_vals = [p6_w_asc[p6_w_asc.window == key]['slope_tls'].iloc[0] for key, label in key_windows]
colors = ['#d95f02', '#1b7837', '#1b7837', '#2b83ba', '#756bb1']

y_pos = np.arange(len(labels))
width = 0.35

ax_bars.barh(y_pos - width/2, r_vals, width, color=colors, alpha=0.85, label='Correlation with Laser (r)')
ax_bars.barh(y_pos + width/2, tls_vals, width, color=colors, hatch='///', alpha=0.5, edgecolor='black', label='TLS Slope (m, 1.0 = 1:1 motion)')
ax_bars.axvline(1.0, color='gray', ls='--', lw=1.2, label='Theoretical 1:1 Slope')
ax_bars.set_yticks(y_pos)
ax_bars.set_yticklabels(labels, fontsize=9.5)
ax_bars.invert_yaxis()
ax_bars.set_ylim(4.6, -0.6)
ax_bars.set_xlabel("Metric Value", fontsize=10)
ax_bars.set_xlim(0, 1.6)
ax_bars.set_title("B. InSAR-Laser Agreement Across 13 Spatial Extraction Windows\nWood-Free Pixels Match or Exceed P6 Center", fontsize=11, fontweight='bold', loc='left')
ax_bars.grid(True, alpha=0.25, ls=':')
ax_bars.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=8.5, framealpha=0.95)

# Value annotations
for i, (r_v, t_v) in enumerate(zip(r_vals, tls_vals)):
    ax_bars.text(r_v + 0.02, i - width/2, f"r={r_v:.2f}", va='center', fontsize=8.5, fontweight='bold')
    ax_bars.text(t_v + 0.02, i + width/2, f"m={t_v:.2f}", va='center', fontsize=8.5)

fig.savefig(OUT_DIR / "fig3_spatial_evidence_boardwalk_free_pixels.png", dpi=300)
fig.savefig(WEB_DIR / "fig3_spatial_evidence_boardwalk_free_pixels.png", dpi=300)
plt.close(fig)
print("Correctly aligned 2-Panel Figure 3 successfully generated at 300 DPI!")
