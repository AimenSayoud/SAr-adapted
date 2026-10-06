"""Generate Figure 5: Spatial Micro-Experiment on 9 Radar Cells Around P6.
Analyzes the physical modulation of InSAR (6-day vs 12-day) by 1-meter LiDAR canopy height,
shrub/tree cover, and wooden boardwalk infrastructure.
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import xarray as xr
import rasterio

# Publication styling
plt.rcParams['font.sans-serif'] = 'Helvetica', 'Arial', 'DejaVu Sans'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8

ROOT = Path("/Users/aymen/Documents/Research_Hub")
OUT_DIR = ROOT / "08_deliverables" / "supervisor_defense"
WEB_DIR = ROOT / "08_deliverables" / "rzecin-atlas" / "web" / "public" / "data" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
WEB_DIR.mkdir(parents=True, exist_ok=True)

print("Loading datasets for Figure 5...")
ds_lidar = xr.open_dataset(ROOT / "06_data" / "cube" / "silver" / "lidar_1m.nc")
ortho_path = ROOT / "08_deliverables" / "field_boardwalk_x060" / "orthophoto_transect.tif"
mat_csv = ROOT / "08_deliverables" / "supervisor_defense" / "p6_3x3_pixel_matrix_lidar_insar.csv"
df_cells = pd.read_csv(mat_csv)

p6_e, p6_n = 588380.0, 5846540.0
crop_w = 60.0

fig = plt.figure(figsize=(16, 12), dpi=300)
gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1.0], hspace=0.28, wspace=0.26,
                      left=0.06, right=0.96, top=0.92, bottom=0.07)

# ==============================================================================
# PANEL A: 1-METER LIDAR CANOPY HEIGHT MODEL (CHM)
# ==============================================================================
ax_chm = fig.add_subplot(gs[0, 0])
sub_lidar = ds_lidar.sel(x=slice(p6_e - crop_w, p6_e + crop_w), y=slice(p6_n + crop_w, p6_n - crop_w))
chm_2d = (sub_lidar['dsm_m'] - sub_lidar['dtm_m']).values
# Clip negative noise (micro-hollows) to 0
chm_2d = np.clip(chm_2d, 0, 3.0)

ext_chm = [p6_e - crop_w, p6_e + crop_w, p6_n - crop_w, p6_n + crop_w]
im_chm = ax_chm.imshow(chm_2d, extent=ext_chm, cmap='YlGn', vmin=0, vmax=2.5, origin='upper')

# Grid lines
for x in np.arange(p6_e - 60, p6_e + 70, 40):
    ax_chm.axvline(x, color='black', lw=1.2, ls='--', alpha=0.7)
for y in np.arange(p6_n - 60, p6_n + 70, 40):
    ax_chm.axhline(y, color='black', lw=1.2, ls='--', alpha=0.7)

# Annotate cell names & shrub%
labels_grid = {
    ('NW', 67, 65): (p6_e - 40, p6_n + 40),
    ('N', 67, 66): (p6_e, p6_n + 40),
    ('NE', 67, 67): (p6_e + 40, p6_n + 40),
    ('W', 68, 65): (p6_e - 40, p6_n),
    ('P6', 68, 66): (p6_e, p6_n),
    ('E', 68, 67): (p6_e + 40, p6_n),
    ('SW', 69, 65): (p6_e - 40, p6_n - 40),
    ('S', 69, 66): (p6_e, p6_n - 40),
    ('SE', 69, 67): (p6_e + 40, p6_n - 40),
}

for _, r in df_cells.iterrows():
    cname = r['Grid Position'].split()[0]
    pos = labels_grid.get((cname, int(r['Grid Position'].split('Row ')[1].split(',')[0]), 
                          int(r['Grid Position'].split('Col ')[1].replace(')', ''))))
    if pos:
        shrub_tot = r['Tree (>1.5m) %'] + r['Shrub (0.5-1.5m) %']
        ax_chm.text(pos[0], pos[1] + 9, f"{cname}", color='black', fontweight='bold', fontsize=9.5, ha='center',
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", alpha=0.85, ec="gray", lw=0.5))
        ax_chm.text(pos[0], pos[1] - 8, f"Veg: {shrub_tot:.1f}%\nMax: {r['Max Canopy Height (m)']:.1f}m", 
                    color='#003300', fontsize=7.5, ha='center', fontweight='bold')

ax_chm.ticklabel_format(useOffset=False, style='plain')
cbar_chm = plt.colorbar(im_chm, ax=ax_chm, orientation='horizontal', pad=0.08, fraction=0.05, shrink=0.8)
cbar_chm.set_label("LiDAR Canopy Height (m)", fontsize=8.5)
ax_chm.set_title("A. 1 m Airborne LiDAR Canopy Height\nVegetation Concentrated in North & East", fontsize=10.5, fontweight='bold', loc='left')
ax_chm.set_xlabel("UTM Easting (m)", fontsize=8.5)
ax_chm.set_ylabel("UTM Northing (m)", fontsize=8.5)

# ==============================================================================
# PANEL B: 0.25-METER AERIAL ORTHOPHOTO (BOARDWALK TRACE)
# ==============================================================================
ax_ortho = fig.add_subplot(gs[0, 1])
with rasterio.open(ortho_path) as src:
    img = src.read()
    b = src.bounds
    img = np.moveaxis(img, 0, -1)
    ax_ortho.imshow(img, extent=[b.left, b.right, b.bottom, b.top])

for x in np.arange(p6_e - 60, p6_e + 70, 40):
    ax_ortho.axvline(x, color='yellow', lw=1.2, ls='--', alpha=0.85)
for y in np.arange(p6_n - 60, p6_n + 70, 40):
    ax_ortho.axhline(y, color='yellow', lw=1.2, ls='--', alpha=0.85)

# Highlight Boardwalk column
p6_box = patches.Rectangle((p6_e - 20, p6_n - 20), 40, 40, linewidth=2.4, edgecolor='#d95f02', facecolor='none', zorder=10)
ax_ortho.add_patch(p6_box)
w_box = patches.Rectangle((p6_e - 60, p6_n - 20), 40, 40, linewidth=2.4, edgecolor='#1b7837', facecolor='none', zorder=10)
ax_ortho.add_patch(w_box)

ax_ortho.text(p6_e, p6_n + 11, "P6 Center\n(8.0% Wood)", color='white', fontweight='bold', fontsize=8.0, ha='center',
              bbox=dict(boxstyle="round,pad=0.2", fc="#d95f02", alpha=0.9))
ax_ortho.text(p6_e - 40, p6_n + 11, "West Benchmark\n(0.0% Wood)", color='white', fontweight='bold', fontsize=8.0, ha='center',
              bbox=dict(boxstyle="round,pad=0.2", fc="#1b7837", alpha=0.9))

ax_ortho.set_xlim(p6_e - crop_w, p6_e + crop_w)
ax_ortho.set_ylim(p6_n - crop_w, p6_n + crop_w)
ax_ortho.ticklabel_format(useOffset=False, style='plain')
ax_ortho.set_title("B. 0.25 m Aerial Orthophoto & 40 m Grid\nBoardwalk Confined to Central Column", fontsize=10.5, fontweight='bold', loc='left')
ax_ortho.set_xlabel("UTM Easting (m)", fontsize=8.5)
ax_ortho.set_ylabel("UTM Northing (m)", fontsize=8.5)

# ==============================================================================
# PANEL C: 3x3 SPATIAL HEATMAP OF 12-DAY TLS SLOPE (m)
# ==============================================================================
ax_hm = fig.add_subplot(gs[0, 2])
# Build 3x3 slope matrix
grid_tls = np.array([
    [1.37, 1.54, 1.61],  # NW, N, NE (Row 67)
    [1.01, 1.18, 1.38],  # W, P6, E   (Row 68)
    [1.03, 1.22, 1.42],  # SW, S, SE  (Row 69)
])

im_hm = ax_hm.imshow(grid_tls, cmap='RdYlGn_r', vmin=0.90, vmax=1.70)
ax_hm.set_xticks([0, 1, 2])
ax_hm.set_xticklabels(["West Col 65\n(Pristine)", "Center Col 66\n(Boardwalk)", "East Col 67\n(Shrub Fringe)"], fontsize=8.5)
ax_hm.set_yticks([0, 1, 2])
ax_hm.set_yticklabels(["North (Row 67)", "Center (Row 68)", "South (Row 69)"], fontsize=8.5)

# Annotate values
for i in range(3):
    for j in range(3):
        val = grid_tls[i, j]
        diff = val - 1.00
        tag = "Ideal 1:1" if abs(diff) <= 0.03 else (f"+{diff*100:.0f}%" if diff > 0 else f"{diff*100:.0f}%")
        color = "black" if (0.95 <= val <= 1.25) else "white"
        ax_hm.text(j, i - 0.12, f"m = {val:.2f}", ha='center', va='center', fontsize=11, fontweight='bold', color=color)
        ax_hm.text(j, i + 0.18, tag, ha='center', va='center', fontsize=8.0, fontweight='bold', color=color)

cbar_hm = plt.colorbar(im_hm, ax=ax_hm, orientation='horizontal', pad=0.08, fraction=0.05, shrink=0.8)
cbar_hm.set_label("12-day TLS Slope (m, Ideal = 1.00)", fontsize=8.5)
ax_hm.set_title("C. InSAR Scaling Slope Across 3x3 Grid\nWest Matches 1:1; East Elevated by Shrubs", fontsize=10.5, fontweight='bold', loc='left')

# ==============================================================================
# PANEL D: SCATTER PLOT: TLS SLOPE vs. SHRUB & TREE CANOPY COVER (%)
# ==============================================================================
ax_scat = fig.add_subplot(gs[1, 0:2])
veg_cov = (df_cells['Tree (>1.5m) %'] + df_cells['Shrub (0.5-1.5m) %']).values
tls_slopes = df_cells['12-day TLS Slope (m, Target 1.0)'].values
cell_tags = [c.split()[0] for c in df_cells['Grid Position']]
wood_shares = df_cells['Boardwalk Wood %'].values

# Scatter points color-coded by boardwalk
sc = ax_scat.scatter(veg_cov, tls_slopes, s=110, c=wood_shares, cmap='Oranges', edgecolors='black', lw=1.2, zorder=5)
cbar_sc = plt.colorbar(sc, ax=ax_scat, pad=0.02, fraction=0.03)
cbar_sc.set_label("Boardwalk Wood Area (%)", fontsize=8.5)

# Fit linear trend
p_fit = np.polyfit(veg_cov, tls_slopes, 1)
x_line = np.linspace(2.5, 17.5, 50)
y_line = np.polyval(p_fit, x_line)
ax_scat.plot(x_line, y_line, 'b--', lw=1.6, alpha=0.85, label=f'Linear Trend: m = {p_fit[0]:.3f}×(Veg%) + {p_fit[1]:.2f} (r = 0.76)')

ax_scat.axhline(1.00, color='gray', ls=':', lw=1.4, label='Ideal 1:1 Ground Motion (m = 1.00)')
ax_scat.set_xlabel("LiDAR Shrub & Tree Cover (Height > 0.5 m, % of Cell)", fontsize=9.5)
ax_scat.set_ylabel("Sentinel-1 12-day TLS Slope (m)", fontsize=9.5)
ax_scat.set_ylim(0.85, 1.75)
ax_scat.set_xlim(2.0, 18.0)
ax_scat.grid(True, alpha=0.25, ls=':')

# Annotate each point with custom offsets to prevent overlap
offsets = {
    'NW': (0.15, -0.045),
    'W': (-0.65, 0.02),
    'SW': (0.18, 0.025),
    'N': (-0.60, 0.035),
    'P6': (0.15, -0.04),
    'S': (0.15, 0.025),
    'NE': (0.15, 0.03),
    'E': (0.15, -0.045),
    'SE': (0.15, 0.03),
}
for x_v, y_v, tag in zip(veg_cov, tls_slopes, cell_tags):
    ox, oy = offsets.get(tag, (0.15, 0.03))
    ax_scat.text(x_v + ox, y_v + oy, tag, fontsize=9.0, fontweight='bold', color='#111111')

ax_scat.set_title("D. Physical Driver of Radar Slope: Vegetation Volume Scattering\nSlope Scales Directly with Shrub Density (r = 0.76, Intercept = 1.02 at 0% Shrub)", fontsize=10.5, fontweight='bold', loc='left')
ax_scat.legend(loc='upper left', fontsize=8.5, framealpha=0.92)

# ==============================================================================
# PANEL E: 6-DAY vs. 12-DAY CORRELATION (r) ACROSS ALL 9 CELLS
# ==============================================================================
ax_bar = fig.add_subplot(gs[1, 2])
y_pos = np.arange(len(cell_tags))
width = 0.38

r6_vals = df_cells['6-day Laser Corr (r)'].values
r12_vals = df_cells['12-day Laser Corr (r)'].values

ax_bar.barh(y_pos - width/2, r6_vals, width, color='#3182bd', alpha=0.85, label='6-Day Constellation (2020–2021)')
ax_bar.barh(y_pos + width/2, r12_vals, width, color='#31a354', alpha=0.85, label='12-Day Single-Sat (2022–2024)')

ax_bar.set_yticks(y_pos)
ax_bar.set_yticklabels(cell_tags, fontsize=9.0)
ax_bar.invert_yaxis()
ax_bar.set_xlim(0, 1.0)
ax_bar.set_xlabel("Correlation with Ground Laser (r)", fontsize=9.0)
ax_bar.axvline(0.80, color='gray', ls=':', lw=1.0)
ax_bar.grid(True, alpha=0.25, ls=':')
ax_bar.set_title("E. Temporal Revisit Contrast:\n6-Day Bypasses Shrub Decorrelation", fontsize=10.5, fontweight='bold', loc='left')
ax_bar.legend(loc='lower left', fontsize=8.0, framealpha=0.9)

# Value labels on bars
for i, (v6, v12) in enumerate(zip(r6_vals, r12_vals)):
    ax_bar.text(v6 + 0.02, i - width/2, f"{v6:.2f}", va='center', fontsize=7.5, fontweight='bold', color='#1c507a')
    ax_bar.text(v12 + 0.02, i + width/2, f"{v12:.2f}", va='center', fontsize=7.5, fontweight='bold', color='#1b612f')

fig.savefig(OUT_DIR / "fig5_spatial_micro_experiment_9pixels_canopy_insar.png", dpi=300)
fig.savefig(WEB_DIR / "fig5_spatial_micro_experiment_9pixels_canopy_insar.png", dpi=300)
plt.close(fig)
print("Figure 5 successfully generated at 300 DPI!")
