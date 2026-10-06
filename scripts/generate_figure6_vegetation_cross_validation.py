"""Generate Figure 6: Multi-Sensor Vegetation Detection & Cross-Validation.
Resolves physical limitations of airborne LiDAR vs high-resolution optical imagery
across the 9 radar cells (230,400 pixels at 0.25 m resolution).
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from scipy.ndimage import uniform_filter, zoom

# Publication styling
plt.rcParams['font.sans-serif'] = 'Helvetica', 'Arial', 'DejaVu Sans'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8

ROOT = Path("/Users/aymen/Documents/Research_Hub")
OUT_DIR = ROOT / "08_deliverables" / "paper_figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

p6_e, p6_n = 588380.0, 5846540.0
crop_w = 60.0 # 120x120 m around P6 (480x480 pixels at 0.25 m)

# 1. Load Orthophoto (0.25 m)
ortho_path = ROOT / "08_deliverables" / "field_boardwalk_x060" / "orthophoto_transect.tif"
with rasterio.open(ortho_path) as src:
    win = rasterio.windows.from_bounds(p6_e - crop_w, p6_n - crop_w, p6_e + crop_w, p6_n + crop_w, src.transform)
    rgb = src.read(window=win)
    ext_m = [p6_e - crop_w, p6_e + crop_w, p6_n - crop_w, p6_n + crop_w]

# 2. Load Boardwalk mask
with rasterio.open(ROOT / "08_deliverables" / "field_boardwalk_x060" / "boardwalk_mask.tif") as src_bw:
    win_bw = rasterio.windows.from_bounds(p6_e - crop_w, p6_n - crop_w, p6_e + crop_w, p6_n + crop_w, src_bw.transform)
    bw_mask = src_bw.read(1, window=win_bw) > 0

r, g, b = rgb[0].astype(float), rgb[1].astype(float), rgb[2].astype(float)
lum = 0.299 * r + 0.587 * g + 0.114 * b

# Local texture: standard deviation in 5x5 window (1.25 m footprint)
mean_lum = uniform_filter(lum, size=5)
sq_mean_lum = uniform_filter(lum**2, size=5)
tex_std = np.sqrt(np.maximum(0, sq_mean_lum - mean_lum**2))

# Optical shrub signature (high local texture / shadows, excluding boardwalk)
opt_shrub = (tex_std > 9.5) & (lum < 130) & (~bw_mask)

# 3. Load LiDAR CHM
ds_lidar = xr.open_dataset(ROOT / "06_data" / "cube" / "silver" / "lidar_1m.nc")
sub_lidar = ds_lidar.sel(x=slice(p6_e - crop_w, p6_e + crop_w), y=slice(p6_n + crop_w, p6_n - crop_w))
chm_1m = np.clip((sub_lidar['dsm_m'] - sub_lidar['dtm_m']).values, 0, 5.0)
chm_25cm = zoom(chm_1m, (480 / chm_1m.shape[0], 480 / chm_1m.shape[1]), order=1)
lidar_shrub = (chm_25cm > 0.5) & (~bw_mask)

# 4. Multi-Sensor Cross-Validation Categorization
# Class 0: Consensus Pure Sphagnum Lawn
# Class 1: Consensus Shrub (Both agree)
# Class 2: LiDAR Missed (Optical Shrub, LiDAR Ground - laser penetration)
# Class 3: LiDAR Ghost (LiDAR > 0.5m, Optical Smooth Moss - hummock artifact)
# Class 4: Boardwalk Structure
discrim_map = np.zeros((480, 480), dtype=int)
discrim_map[opt_shrub & lidar_shrub] = 1
discrim_map[opt_shrub & (~lidar_shrub)] = 2
discrim_map[(~opt_shrub) & lidar_shrub] = 3
discrim_map[bw_mask] = 4

# Figure setup: 5 panels (Gridspec 2 rows, 3 cols)
fig = plt.figure(figsize=(16, 11), dpi=300)
gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.95], hspace=0.26, wspace=0.22,
                      left=0.06, right=0.95, top=0.92, bottom=0.07)

# PANEL A: AERIAL ORTHOPHOTO (0.25 m)
ax_a = fig.add_subplot(gs[0, 0])
rgb_disp = np.moveaxis(rgb, 0, -1)
ax_a.imshow(rgb_disp, extent=ext_m)
for x in np.arange(p6_e - 60, p6_e + 70, 40):
    ax_a.axvline(x, color='yellow', lw=1.0, ls='--', alpha=0.85)
for y in np.arange(p6_n - 60, p6_n + 70, 40):
    ax_a.axhline(y, color='yellow', lw=1.0, ls='--', alpha=0.85)
ax_a.ticklabel_format(useOffset=False, style='plain')
ax_a.set_title("A. 0.25 m Aerial Orthophoto (RGB)\nVisual Shrub Crowns & Boardwalk Trace", fontsize=10.0, fontweight='bold', loc='left')
ax_a.set_xlabel("UTM Easting (m)", fontsize=8.5)
ax_a.set_ylabel("UTM Northing (m)", fontsize=8.5)

# PANEL B: 1 m AIRBORNE LIDAR CHM
ax_b = fig.add_subplot(gs[0, 1])
im_b = ax_b.imshow(chm_25cm, extent=ext_m, cmap='YlGn', vmin=0, vmax=2.5, origin='upper')
for x in np.arange(p6_e - 60, p6_e + 70, 40):
    ax_b.axvline(x, color='black', lw=1.0, ls='--', alpha=0.7)
for y in np.arange(p6_n - 60, p6_n + 70, 40):
    ax_b.axhline(y, color='black', lw=1.0, ls='--', alpha=0.7)
ax_b.ticklabel_format(useOffset=False, style='plain')
cbar_b = plt.colorbar(im_b, ax=ax_b, orientation='horizontal', pad=0.08, fraction=0.05, shrink=0.8)
cbar_b.set_label("LiDAR Canopy Height (m)", fontsize=8.0)
ax_b.set_title("B. 1 m Airborne LiDAR CHM\nMisses Wiry Scrub; Flags Hummocks", fontsize=10.0, fontweight='bold', loc='left')
ax_b.set_xlabel("UTM Easting (m)", fontsize=8.5)
ax_b.set_ylabel("UTM Northing (m)", fontsize=8.5)

# PANEL C: OPTICAL COMPUTER VISION TEXTURE & SHADOW
ax_c = fig.add_subplot(gs[0, 2])
im_c = ax_c.imshow(tex_std, extent=ext_m, cmap='magma', vmin=3, vmax=18, origin='upper')
for x in np.arange(p6_e - 60, p6_e + 70, 40):
    ax_c.axvline(x, color='cyan', lw=1.0, ls='--', alpha=0.7)
for y in np.arange(p6_n - 60, p6_n + 70, 40):
    ax_c.axhline(y, color='cyan', lw=1.0, ls='--', alpha=0.7)
ax_c.ticklabel_format(useOffset=False, style='plain')
cbar_c = plt.colorbar(im_c, ax=ax_c, orientation='horizontal', pad=0.08, fraction=0.05, shrink=0.8)
cbar_c.set_label("Local Luminance Texture (σ)", fontsize=8.0)
ax_c.set_title("C. Computer Vision Crown Texture\nIdentifies Micro-Canopies & Shadows", fontsize=10.0, fontweight='bold', loc='left')
ax_c.set_xlabel("UTM Easting (m)", fontsize=8.5)
ax_c.set_ylabel("UTM Northing (m)", fontsize=8.5)

# PANEL D: CROSS-VALIDATION DISCREPANCY MAP
ax_d = fig.add_subplot(gs[1, 0:2])
# Custom colormap for categories
from matplotlib.colors import ListedColormap
cmap_disc = ListedColormap(['#f0f0f0', '#2ca02c', '#d62728', '#ff7f0e', '#7f7f7f'])
# 0: Lawn (grey), 1: Consensus (green), 2: LiDAR Missed (red), 3: LiDAR Ghost (orange), 4: Boardwalk (dark grey)
im_d = ax_d.imshow(discrim_map, extent=ext_m, cmap=cmap_disc, vmin=0, vmax=4, origin='upper')
for x in np.arange(p6_e - 60, p6_e + 70, 40):
    ax_d.axvline(x, color='black', lw=1.2, ls='--', alpha=0.8)
for y in np.arange(p6_n - 60, p6_n + 70, 40):
    ax_d.axhline(y, color='black', lw=1.2, ls='--', alpha=0.8)

# Add cell labels
grid_tags = [
    ('NW', 0, 0), ('N', 0, 1), ('NE', 0, 2),
    ('W',  1, 0), ('P6', 1, 1), ('E',  1, 2),
    ('SW', 2, 0), ('S',  2, 1), ('SE', 2, 2)
]
for tag, ri, ci in grid_tags:
    cx = p6_e - 40 + ci * 40
    cy = p6_n + 40 - ri * 40
    ax_d.text(cx, cy, tag, color='black', fontweight='bold', fontsize=9.0, ha='center', va='center',
              bbox=dict(boxstyle="round,pad=0.15", fc="white", alpha=0.85, ec="gray", lw=0.5))

ax_d.ticklabel_format(useOffset=False, style='plain')
ax_d.set_title("D. Sensor Discrepancy & Cross-Validation Map (230,400 Pixels)\nResolving Physical Laser Penetration vs Optical Shadows", fontsize=10.5, fontweight='bold', loc='left')
ax_d.set_xlabel("UTM Easting (m)", fontsize=8.5)
ax_d.set_ylabel("UTM Northing (m)", fontsize=8.5)

# Custom legend for Panel D
legend_elements = [
    patches.Patch(facecolor='#2ca02c', edgecolor='black', label='Consensus Shrub (Both Agree, 1.7%)'),
    patches.Patch(facecolor='#d62728', edgecolor='black', label='LiDAR Missed (Optical Shrub, Laser Penetrated, 10.8%)'),
    patches.Patch(facecolor='#ff7f0e', edgecolor='black', label='LiDAR Ghost (Hummock/Ground Artifact, 5.5%)'),
    patches.Patch(facecolor='#7f7f7f', edgecolor='black', label='Engineered Boardwalk Structure (2.4%)'),
    patches.Patch(facecolor='#f0f0f0', edgecolor='gray', label='Consensus Pure Sphagnum Lawn (79.6%)')
]
ax_d.legend(handles=legend_elements, loc='upper left', fontsize=8.0, framealpha=0.92)

# PANEL E: PER-CELL COMPARISON (BAR CHART)
ax_e = fig.add_subplot(gs[1, 2])
cell_names = [t[0] for t in grid_tags]
opt_vals, lid_vals, missed_vals, ghost_vals = [], [], [], []

for tag, ri, ci in grid_tags:
    sub_map = discrim_map[ri*160:(ri+1)*160, ci*160:(ci+1)*160]
    sub_lid = lidar_shrub[ri*160:(ri+1)*160, ci*160:(ci+1)*160]
    sub_opt = opt_shrub[ri*160:(ri+1)*160, ci*160:(ci+1)*160]
    opt_vals.append(sub_opt.mean() * 100)
    lid_vals.append(sub_lid.mean() * 100)
    missed_vals.append((sub_map == 2).mean() * 100)
    ghost_vals.append((sub_map == 3).mean() * 100)

y_p = np.arange(len(cell_names))
bar_w = 0.38

ax_e.barh(y_p - bar_w/2, lid_vals, height=bar_w, color='#31a354', alpha=0.85, label='LiDAR CHM > 0.5 m')
ax_e.barh(y_p + bar_w/2, opt_vals, height=bar_w, color='#d62728', alpha=0.85, label='Optical Texture & Shadow')

ax_e.set_yticks(y_p)
ax_e.set_yticklabels(cell_names, fontsize=8.5, fontweight='bold')
ax_e.invert_yaxis()
ax_e.set_xlim(0, 25)
ax_e.set_xlabel("Vegetation Area Share (%)", fontsize=8.5)
ax_e.grid(True, alpha=0.25, ls=':')
ax_e.legend(loc='lower right', fontsize=8.0, framealpha=0.92)
ax_e.set_title("E. Sensor Discrepancy by Cell\nLiDAR Underestimates East & North", fontsize=10.0, fontweight='bold', loc='left')

for idx, (lv, ov) in enumerate(zip(lid_vals, opt_vals)):
    ax_e.text(lv + 0.5, idx - bar_w/2, f"{lv:.1f}%", va='center', fontsize=7.0, color='#1b612f', fontweight='bold')
    ax_e.text(ov + 0.5, idx + bar_w/2, f"{ov:.1f}%", va='center', fontsize=7.0, color='#a51d24', fontweight='bold')

fig.suptitle("Cross-Validation of Peatland Vegetation Detection: Airborne LiDAR vs High-Resolution Orthophoto", 
             fontsize=12.5, fontweight='bold', y=0.98)

# Save
fig.savefig(OUT_DIR / "Fig6_Vegetation_Cross_Validation_LiDAR_vs_Optical.png", dpi=300)
plt.close(fig)
print("Figure 6 successfully generated!")
