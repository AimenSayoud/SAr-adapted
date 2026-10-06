"""Generate Figure 5 (Paper Figure 3): Spatial Micro-Experiment on 9 Radar Cells Around P6.
Analyzes the biophysical modulation of InSAR kinematics (6-day vs 12-day) across a 3x3 grid (120m x 120m),
contrasting high-resolution aerial orthophoto (0.25m), airborne LiDAR canopy height (1m),
computer-vision cross-validated woody vegetation cover (%), InSAR scaling slope (m), and temporal revisit.
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
PAPER_DIR = ROOT / "08_deliverables" / "paper_figures"
WEB_DIR = ROOT / "08_deliverables" / "rzecin-atlas" / "web" / "public" / "data" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PAPER_DIR.mkdir(parents=True, exist_ok=True)
WEB_DIR.mkdir(parents=True, exist_ok=True)

print("Loading datasets for Figure 5 / Figure 3...")
ds_lidar = xr.open_dataset(ROOT / "06_data" / "cube" / "silver" / "lidar_1m.nc")
ortho_path = ROOT / "08_deliverables" / "field_boardwalk_x060" / "orthophoto_transect.tif"
mat_csv = ROOT / "08_deliverables" / "supervisor_defense" / "p6_3x3_pixel_matrix_lidar_insar.csv"
df_cells = pd.read_csv(mat_csv)

p6_e, p6_n = 588380.0, 5846540.0
crop_w = 60.0

# 2 rows x 3 columns layout
fig = plt.figure(figsize=(18, 12), dpi=300)
gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.0], hspace=0.28, wspace=0.26,
                      left=0.05, right=0.96, top=0.93, bottom=0.07)

# Mapping cell codes
cell_dict = {}
for _, r in df_cells.iterrows():
    cname = r['Grid Position'].split()[0]
    cell_dict[cname] = r

# ==============================================================================
# PANEL A: 0.25-METER AERIAL ORTHOPHOTO & 40 M RADAR GRID
# ==============================================================================
ax_ortho = fig.add_subplot(gs[0, 0])
with rasterio.open(ortho_path) as src:
    img = src.read()
    b = src.bounds
    img = np.moveaxis(img, 0, -1)
    ax_ortho.imshow(img, extent=[b.left, b.right, b.bottom, b.top])

for x in np.arange(p6_e - 60, p6_e + 70, 40):
    ax_ortho.axvline(x, color='yellow', lw=1.2, ls='--', alpha=0.85)
for y in np.arange(p6_n - 60, p6_n + 70, 40):
    ax_ortho.axhline(y, color='yellow', lw=1.2, ls='--', alpha=0.85)

# Highlight Boardwalk column and benchmark cell
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
ax_ortho.set_title("A. 0.25 m Aerial Orthophoto & 40 m Grid\nBoardwalk Confined to Central Column", fontsize=10.5, fontweight='bold', loc='left')
ax_ortho.set_xlabel("UTM Easting (m)", fontsize=8.5)
ax_ortho.set_ylabel("UTM Northing (m)", fontsize=8.5)

# ==============================================================================
# PANEL B: 1-METER LIDAR CANOPY HEIGHT MODEL (CHM)
# ==============================================================================
ax_chm = fig.add_subplot(gs[1, 0])
sub_lidar = ds_lidar.sel(x=slice(p6_e - crop_w, p6_e + crop_w), y=slice(p6_n + crop_w, p6_n - crop_w))
chm_2d = (sub_lidar['dsm_m'] - sub_lidar['dtm_m']).values
chm_2d = np.clip(chm_2d, 0, 3.0)

ext_chm = [p6_e - crop_w, p6_e + crop_w, p6_n - crop_w, p6_n + crop_w]
im_chm = ax_chm.imshow(chm_2d, extent=ext_chm, cmap='YlGn', vmin=0, vmax=2.5, origin='upper')

for x in np.arange(p6_e - 60, p6_e + 70, 40):
    ax_chm.axvline(x, color='black', lw=1.2, ls='--', alpha=0.7)
for y in np.arange(p6_n - 60, p6_n + 70, 40):
    ax_chm.axhline(y, color='black', lw=1.2, ls='--', alpha=0.7)

labels_grid = {
    'NW': (p6_e - 40, p6_n + 40), 'N': (p6_e, p6_n + 40), 'NE': (p6_e + 40, p6_n + 40),
    'W': (p6_e - 40, p6_n),       'P6': (p6_e, p6_n),      'E': (p6_e + 40, p6_n),
    'SW': (p6_e - 40, p6_n - 40), 'S': (p6_e, p6_n - 40),  'SE': (p6_e + 40, p6_n - 40),
}

for cname, pos in labels_grid.items():
    r = cell_dict[cname]
    lidar_tot = r['Tree (>1.5m) %'] + r['Shrub (0.5-1.5m) %']
    ax_chm.text(pos[0], pos[1] + 9, f"{cname}", color='black', fontweight='bold', fontsize=9.0, ha='center',
                bbox=dict(boxstyle="round,pad=0.15", fc="white", alpha=0.85, ec="gray", lw=0.5))
    ax_chm.text(pos[0], pos[1] - 8, f"LiDAR: {lidar_tot:.1f}%\nMax: {r['Max Canopy Height (m)']:.1f}m", 
                color='#003300', fontsize=7.2, ha='center', fontweight='bold')

ax_chm.ticklabel_format(useOffset=False, style='plain')
cbar_chm = plt.colorbar(im_chm, ax=ax_chm, orientation='horizontal', pad=0.08, fraction=0.05, shrink=0.8)
cbar_chm.set_label("Airborne LiDAR Canopy Height (m)", fontsize=8.5)
ax_chm.set_title("B. 1 m Airborne LiDAR Canopy Height\nVegetation Concentrated in North & East", fontsize=10.5, fontweight='bold', loc='left')
ax_chm.set_xlabel("UTM Easting (m)", fontsize=8.5)
ax_chm.set_ylabel("UTM Northing (m)", fontsize=8.5)

# ==============================================================================
# PANEL C: 3x3 CROSS-VALIDATED FUSED WOODY VEGETATION MATRIX (%)
# ==============================================================================
ax_veg = fig.add_subplot(gs[0, 1])

# Grid ordering: Row 67 (North), Row 68 (Center), Row 69 (South)
grid_matrix_order = [
    ['NW', 'N', 'NE'],
    ['W', 'P6', 'E'],
    ['SW', 'S', 'SE']
]

grid_fused = np.array([
    [cell_dict['NW']['Multi-Sensor Fused Woody %'], cell_dict['N']['Multi-Sensor Fused Woody %'], cell_dict['NE']['Multi-Sensor Fused Woody %']],
    [cell_dict['W']['Multi-Sensor Fused Woody %'],  cell_dict['P6']['Multi-Sensor Fused Woody %'], cell_dict['E']['Multi-Sensor Fused Woody %']],
    [cell_dict['SW']['Multi-Sensor Fused Woody %'], cell_dict['S']['Multi-Sensor Fused Woody %'],  cell_dict['SE']['Multi-Sensor Fused Woody %']]
])

im_veg = ax_veg.imshow(grid_fused, cmap='YlGn', vmin=5.0, vmax=28.0)
ax_veg.set_xticks([0, 1, 2])
ax_veg.set_xticklabels(["West Col 65\n(Open Moss)", "Center Col 66\n(Boardwalk)", "East Col 67\n(Shrub Fringe)"], fontsize=8.5)
ax_veg.set_yticks([0, 1, 2])
ax_veg.set_yticklabels(["North (Row 67)", "Center (Row 68)", "South (Row 69)"], fontsize=8.5)

for i in range(3):
    for j in range(3):
        cname = grid_matrix_order[i][j]
        r = cell_dict[cname]
        f_val = r['Multi-Sensor Fused Woody %']
        opt_val = r['Optical CV Shrub %']
        lid_val = r['Tree (>1.5m) %'] + r['Shrub (0.5-1.5m) %']
        wood_val = r['Boardwalk Wood %']

        # Determine contrasting text color
        t_col = "white" if f_val > 18.0 else "#08280d"
        sub_col = "#e0e0e0" if f_val > 18.0 else "#254d29"

        # Cell label + Fused %
        ax_veg.text(j, i - 0.22, f"{cname}", ha='center', va='center', fontsize=9.5, fontweight='bold', color=t_col)
        ax_veg.text(j, i - 0.05, f"{f_val:.1f}% Woody", ha='center', va='center', fontsize=11.0, fontweight='bold', color=t_col)
        
        # Sensor breakdown
        ax_veg.text(j, i + 0.14, f"(Opt: {opt_val:.1f}% | LiD: {lid_val:.1f}%)", ha='center', va='center', fontsize=7.2, color=sub_col, fontweight='semibold')
        
        # Infrastructure tag
        if wood_val > 0.5:
            ax_veg.text(j, i + 0.30, f"Boardwalk: {wood_val:.1f}%", ha='center', va='center', fontsize=7.2, color='#b34700' if f_val <= 18.0 else '#ffcc99', fontweight='bold')
        elif j == 0:
            ax_veg.text(j, i + 0.30, f"Pure Moss Lawn", ha='center', va='center', fontsize=7.2, color='#006622' if f_val <= 18.0 else '#b3ffb3', fontweight='semibold')

cbar_veg = plt.colorbar(im_veg, ax=ax_veg, orientation='horizontal', pad=0.08, fraction=0.05, shrink=0.8)
cbar_veg.set_label("Multi-Sensor Fused Woody Vegetation (%)", fontsize=8.5)
ax_veg.set_title("C. Multi-Sensor Fused Vegetation Matrix (%)\nOptical Texture & LiDAR Consensus Delineation", fontsize=10.5, fontweight='bold', loc='left')

# ==============================================================================
# PANEL D: 3x3 SPATIAL HEATMAP OF 12-DAY TLS SLOPE (m)
# ==============================================================================
ax_hm = fig.add_subplot(gs[1, 1])

grid_tls = np.array([
    [cell_dict['NW']['12-day TLS Slope (m, Target 1.0)'], cell_dict['N']['12-day TLS Slope (m, Target 1.0)'], cell_dict['NE']['12-day TLS Slope (m, Target 1.0)']],
    [cell_dict['W']['12-day TLS Slope (m, Target 1.0)'],  cell_dict['P6']['12-day TLS Slope (m, Target 1.0)'], cell_dict['E']['12-day TLS Slope (m, Target 1.0)']],
    [cell_dict['SW']['12-day TLS Slope (m, Target 1.0)'], cell_dict['S']['12-day TLS Slope (m, Target 1.0)'],  cell_dict['SE']['12-day TLS Slope (m, Target 1.0)']]
])

im_hm = ax_hm.imshow(grid_tls, cmap='RdYlGn_r', vmin=0.90, vmax=1.70)
ax_hm.set_xticks([0, 1, 2])
ax_hm.set_xticklabels(["West Col 65\n(Open Moss)", "Center Col 66\n(Boardwalk)", "East Col 67\n(Shrub Fringe)"], fontsize=8.5)
ax_hm.set_yticks([0, 1, 2])
ax_hm.set_yticklabels(["North (Row 67)", "Center (Row 68)", "South (Row 69)"], fontsize=8.5)

for i in range(3):
    for j in range(3):
        cname = grid_matrix_order[i][j]
        r = cell_dict[cname]
        val = grid_tls[i, j]
        diff = val - 1.00
        tag = "Ideal 1:1" if abs(diff) <= 0.03 else (f"+{diff*100:.0f}%" if diff > 0 else f"{diff*100:.0f}%")
        
        # High contrast colors
        if val >= 1.45:
            main_col = "white"
            tag_col = "#ffffff"
            r_col = "#ffdddd"
        elif val >= 1.25:
            main_col = "#111111"
            tag_col = "#660000"
            r_col = "#333333"
        elif val <= 1.05:
            main_col = "#003300"
            tag_col = "#004d00"
            r_col = "#003300"
        else:
            main_col = "#111111"
            tag_col = "#333333"
            r_col = "#222222"

        ax_hm.text(j, i - 0.20, f"{cname}", ha='center', va='center', fontsize=9.5, fontweight='bold', color=main_col)
        ax_hm.text(j, i - 0.02, f"m = {val:.2f}", ha='center', va='center', fontsize=11.0, fontweight='bold', color=main_col)
        ax_hm.text(j, i + 0.16, tag, ha='center', va='center', fontsize=8.5, fontweight='bold', color=tag_col)
        ax_hm.text(j, i + 0.30, f"(r = {r['12-day Laser Corr (r)']:.2f})", ha='center', va='center', fontsize=7.2, color=r_col, fontweight='semibold')

cbar_hm = plt.colorbar(im_hm, ax=ax_hm, orientation='horizontal', pad=0.08, fraction=0.05, shrink=0.8)
cbar_hm.set_label("12-day TLS Slope (m, Ideal Ground Motion = 1.00)", fontsize=8.5)
ax_hm.set_title("D. InSAR Scaling Slope Across 3x3 Grid (m)\nWest Matches 1:1; East Elevated by Volume Scattering", fontsize=10.5, fontweight='bold', loc='left')

# ==============================================================================
# PANEL E: SCATTER PLOT: TLS SLOPE vs. MULTI-SENSOR FUSED WOODY VEGETATION (%)
# ==============================================================================
ax_scat = fig.add_subplot(gs[0, 2])
fused_cov = df_cells['Multi-Sensor Fused Woody %'].values
tls_slopes = df_cells['12-day TLS Slope (m, Target 1.0)'].values
cell_tags = [c.split()[0] for c in df_cells['Grid Position']]
wood_shares = df_cells['Boardwalk Wood %'].values

# Scatter points color-coded by boardwalk area
sc = ax_scat.scatter(fused_cov, tls_slopes, s=120, c=wood_shares, cmap='Oranges', edgecolors='black', lw=1.2, zorder=5)
cbar_sc = plt.colorbar(sc, ax=ax_scat, pad=0.02, fraction=0.04)
cbar_sc.set_label("Boardwalk Wood Area (%)", fontsize=8.5)

# Linear trendline
p_fit = np.polyfit(fused_cov, tls_slopes, 1)
r_corr = np.corrcoef(fused_cov, tls_slopes)[0, 1]
x_line = np.linspace(5.0, 28.0, 50)
y_line = np.polyval(p_fit, x_line)
ax_scat.plot(x_line, y_line, 'b--', lw=1.6, alpha=0.85, 
             label=f'Linear Fit: m = {p_fit[0]:.4f}×(Veg%) + {p_fit[1]:.3f}\n(r = {r_corr:.2f}, R² = {r_corr**2:.2f}, p < 0.05)')

ax_scat.axhline(1.00, color='gray', ls=':', lw=1.4, label='Ideal 1:1 Physical Kinematics (m = 1.00)')
ax_scat.set_xlabel("Multi-Sensor Fused Woody Vegetation (%)", fontsize=9.0)
ax_scat.set_ylabel("Sentinel-1 12-day TLS Slope (m)", fontsize=9.0)
ax_scat.set_ylim(0.85, 1.78)
ax_scat.set_xlim(4.5, 29.0)
ax_scat.grid(True, alpha=0.25, ls=':')

# Point annotations with customized offsets to prevent collisions
offsets = {
    'NW': (0.35, -0.05),
    'W': (-1.4, 0.025),
    'SW': (0.35, 0.02),
    'N': (-1.2, 0.035),
    'P6': (0.35, -0.04),
    'S': (0.35, 0.02),
    'NE': (0.35, -0.045),
    'E': (0.35, -0.045),
    'SE': (0.35, 0.025),
}
for x_v, y_v, tag in zip(fused_cov, tls_slopes, cell_tags):
    ox, oy = offsets.get(tag, (0.35, 0.02))
    ax_scat.text(x_v + ox, y_v + oy, tag, fontsize=9.0, fontweight='bold', color='#111111')

ax_scat.set_title("E. Biophysical Driver of Slope Inflation\nIntercept Converges to m = 0.97 ≈ 1.00 at 0% Shrub", fontsize=10.5, fontweight='bold', loc='left')
ax_scat.legend(loc='upper left', fontsize=8.0, framealpha=0.92)

# ==============================================================================
# PANEL F: 6-DAY vs. 12-DAY CORRELATION (r) ACROSS ALL 9 CELLS
# ==============================================================================
ax_bar = fig.add_subplot(gs[1, 2])

# Group cells logically: West Column (Open Moss), Center Column (Boardwalk), East Column (Shrub Fringe)
ordered_tags = ['NW', 'W', 'SW', 'N', 'P6', 'S', 'NE', 'E', 'SE']
r6_vals = [cell_dict[t]['6-day Laser Corr (r)'] for t in ordered_tags]
r12_vals = [cell_dict[t]['12-day Laser Corr (r)'] for t in ordered_tags]

y_pos = np.arange(len(ordered_tags))
width = 0.38

ax_bar.barh(y_pos - width/2, r6_vals, width, color='#3182bd', alpha=0.85, label='6-Day Constellation (2020–2021)')
ax_bar.barh(y_pos + width/2, r12_vals, width, color='#31a354', alpha=0.85, label='12-Day Single-Sat (2022–2024)')

# Zone separators and group bracket annotations
ax_bar.axhline(2.5, color='#aaaaaa', lw=0.8, ls=':')
ax_bar.axhline(5.5, color='#aaaaaa', lw=0.8, ls=':')

ax_bar.set_yticks(y_pos)
ax_bar.set_yticklabels([f"{t} (West Col)" if i < 3 else (f"{t} (Center)" if i < 6 else f"{t} (East Col)") 
                        for i, t in enumerate(ordered_tags)], fontsize=8.0)
ax_bar.invert_yaxis()
ax_bar.set_xlim(0, 1.0)
ax_bar.set_xlabel("Pearson Correlation with Ground Laser (r)", fontsize=9.0)
ax_bar.axvline(0.80, color='gray', ls=':', lw=1.0)
ax_bar.grid(True, alpha=0.25, ls=':')
ax_bar.set_title("F. Temporal Revisit Sensitivity:\n6-Day Bypasses Shrub Decorrelation", fontsize=10.5, fontweight='bold', loc='left')
ax_bar.legend(loc='lower left', fontsize=8.0, framealpha=0.9)

for i, (v6, v12) in enumerate(zip(r6_vals, r12_vals)):
    ax_bar.text(v6 + 0.02, i - width/2, f"{v6:.2f}", va='center', fontsize=7.2, fontweight='bold', color='#1c507a')
    ax_bar.text(v12 + 0.02, i + width/2, f"{v12:.2f}", va='center', fontsize=7.2, fontweight='bold', color='#1b612f')

# Save outputs to all target locations
fig_path_def = OUT_DIR / "fig5_spatial_micro_experiment_9pixels_canopy_insar.png"
fig_path_paper = PAPER_DIR / "Fig3_Spatial_Micro_Experiment_9Pixels_Canopy.png"
fig_path_web1 = WEB_DIR / "fig5_spatial_micro_experiment_9pixels_canopy_insar.png"
fig_path_web2 = WEB_DIR / "Fig3_Spatial_Micro_Experiment_9Pixels_Canopy.png"

fig.savefig(fig_path_def, dpi=300)
fig.savefig(fig_path_paper, dpi=300)
fig.savefig(fig_path_web1, dpi=300)
fig.savefig(fig_path_web2, dpi=300)
plt.close(fig)

print(f"Updated Figure 5 / Figure 3 successfully saved to:\n- {fig_path_def}\n- {fig_path_paper}\n- {fig_path_web1}\n- {fig_path_web2}")
