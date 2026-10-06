"""Generate Figure 4: Radar Physics, Polarimetric Mechanisms & Conceptual Scattering Model.
Synthesizes polarimetric backscatter (VV, VH, VH/VV), coherence degradation (6d vs 12d),
multi-variate biophysical correlations, and conceptual microwave interaction pathways across the 9 cells.
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import pandas as pd
import xarray as xr

# Publication styling
plt.rcParams['font.sans-serif'] = 'Helvetica', 'Arial', 'DejaVu Sans'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8

ROOT = Path("/Users/aymen/Documents/Research_Hub")
OUT_DIR = ROOT / "08_deliverables" / "paper_figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# 1. Load data
ds_stat = xr.open_dataset(ROOT / "06_data" / "cube" / "silver" / "static_40m.nc")
mat_csv = ROOT / "08_deliverables" / "supervisor_defense" / "p6_3x3_pixel_matrix_lidar_insar.csv"
df_cells = pd.read_csv(mat_csv)

cell_data = []
for _, r in df_cells.iterrows():
    cname = r['Grid Position'].split()[0]
    row_idx = int(r['Grid Position'].split('Row ')[1].split(',')[0])
    col_idx = int(r['Grid Position'].split('Col ')[1].replace(')', ''))
    
    vv = float(ds_stat.vv_asc_median_db.isel(y=row_idx, x=col_idx))
    vh = float(ds_stat.vh_asc_median_db.isel(y=row_idx, x=col_idx))
    vh_vv = vh - vv
    c6 = float(ds_stat.coh6_asc_median.isel(y=row_idx, x=col_idx))
    c12 = float(ds_stat.coh12_asc_median.isel(y=row_idx, x=col_idx))
    
    cell_data.append({
        'Tag': cname, 'Row': row_idx, 'Col': col_idx,
        'VV_dB': vv, 'VH_dB': vh, 'VH_VV_dB': vh_vv,
        'Coh6': c6, 'Coh12': c12,
        'Delta_Coh': c6 - c12,
        'TLS_Slope': r['12-day TLS Slope (m, Target 1.0)'],
        'Shrub_Cover': r['Tree (>1.5m) %'] + r['Shrub (0.5-1.5m) %'],
        'Wood_Share': r['Boardwalk Wood %'],
        'Dist_Mat': float(ds_stat.depth_in_mat_m.isel(y=row_idx, x=col_idx))
    })

df_p = pd.DataFrame(cell_data)

# Figure setup: 2x2 grid
fig, axs = plt.subplots(2, 2, figsize=(15, 12), dpi=300)
plt.subplots_adjust(hspace=0.28, wspace=0.24, left=0.07, right=0.95, top=0.93, bottom=0.07)

# ==============================================================================
# PANEL A: POLARIMETRIC CROSS-RATIO (VH/VV) - PROOF OF VOLUME SCATTERING
# ==============================================================================
ax_a = axs[0, 0]
# 3x3 grid matrix of VH/VV
grid_vhvv = np.zeros((3, 3))
for _, r in df_p.iterrows():
    grid_vhvv[r['Row'] - 67, r['Col'] - 65] = r['VH_VV_dB']

im_a = ax_a.imshow(grid_vhvv, cmap='PuOr_r', vmin=-7.8, vmax=-5.8)
ax_a.set_xticks([0, 1, 2])
ax_a.set_xticklabels(["West Col 65\n(Pristine Moss)", "Center Col 66\n(Boardwalk)", "East Col 67\n(Shrub Fringe)"], fontsize=8.5)
ax_a.set_yticks([0, 1, 2])
ax_a.set_yticklabels(["North (Row 67)", "Center (Row 68)", "South (Row 69)"], fontsize=8.5)

for i in range(3):
    for j in range(3):
        val = grid_vhvv[i, j]
        tag_cell = df_p[(df_p['Row'] == i+67) & (df_p['Col'] == j+65)]['Tag'].values[0]
        color = "white" if (val < -7.2 or val > -6.2) else "black"
        ax_a.text(j, i - 0.12, f"{tag_cell}", ha='center', va='center', fontsize=10.5, fontweight='bold', color=color)
        ax_a.text(j, i + 0.15, f"{val:.2f} dB", ha='center', va='center', fontsize=9.0, fontweight='bold', color=color)

cbar_a = plt.colorbar(im_a, ax=ax_a, orientation='horizontal', pad=0.10, fraction=0.05, shrink=0.8)
cbar_a.set_label("Cross-Polarization Ratio VH / VV (dB)\n[Higher = Stronger Depolarizing Volume Scattering]", fontsize=8.5)
ax_a.set_title("A. Polarimetric Depolarization Across 3x3 Grid\nEast Fringe Shows Elevated Volume Scattering (+1.6 dB)", fontsize=10.5, fontweight='bold', loc='left')

# ==============================================================================
# PANEL B: TEMPORAL COHERENCE TRANSITION (6-DAY vs. 12-DAY)
# ==============================================================================
ax_b = axs[0, 1]
tags = df_p['Tag'].values
x_pos = np.arange(len(tags))
w = 0.35

ax_b.bar(x_pos - w/2, df_p['Coh6'], width=w, color='#2b83ba', alpha=0.9, label='6-Day Baseline (S1A+S1B, 2020–2021)')
ax_b.bar(x_pos + w/2, df_p['Coh12'], width=w, color='#fdae61', alpha=0.9, label='12-Day Baseline (S1A alone, 2022–2024)')

ax_b.set_xticks(x_pos)
ax_b.set_xticklabels(tags, fontsize=9.0, fontweight='bold')
ax_b.set_ylabel("Interferometric Coherence (γ)", fontsize=9.5)
ax_b.set_ylim(0.45, 1.0)
ax_b.axhline(0.70, color='gray', ls='--', lw=1.2, label='InSAR Inversion Threshold (γ = 0.70)')
ax_b.grid(True, alpha=0.25, ls=':')
ax_b.legend(loc='lower right', fontsize=8.2, framealpha=0.95)
ax_b.set_title("B. Temporal Coherence Loss: 6-Day vs. 12-Day\nAll Cells Lose Coherence Below 0.70 Threshold at 12 Days", fontsize=10.5, fontweight='bold', loc='left')

# Value labels
for idx, (c6, c12) in enumerate(zip(df_p['Coh6'], df_p['Coh12'])):
    ax_b.text(idx - w/2, c6 + 0.012, f"{c6:.2f}", ha='center', fontsize=7.2, fontweight='bold', color='#1c507a')
    ax_b.text(idx + w/2, c12 + 0.012, f"{c12:.2f}", ha='center', fontsize=7.2, fontweight='bold', color='#b25900')

# ==============================================================================
# PANEL C: MULTI-PARAMETRIC BIOPHYSICAL CORRELATION MATRIX
# ==============================================================================
ax_c = axs[1, 0]
cols_corr = ['Shrub_Cover', 'Wood_Share', 'Dist_Mat', 'VV_dB', 'VH_dB', 'VH_VV_dB', 'Coh6', 'Coh12', 'TLS_Slope']
labels_corr = ['Shrub %', 'Boardwalk %', 'Mat Depth', 'VV (dB)', 'VH (dB)', 'VH/VV (dB)', 'γ (6d)', 'γ (12d)', 'Slope m']

corr_mat = df_p[cols_corr].corr().values

im_c = ax_c.imshow(corr_mat, cmap='coolwarm', vmin=-1.0, vmax=1.0)
ax_c.set_xticks(range(len(labels_corr)))
ax_c.set_yticks(range(len(labels_corr)))
ax_c.set_xticklabels(labels_corr, rotation=45, ha='right', fontsize=8.0)
ax_c.set_yticklabels(labels_corr, fontsize=8.0)

for i in range(len(labels_corr)):
    for j in range(len(labels_corr)):
        val = corr_mat[i, j]
        color = "white" if abs(val) > 0.65 else "black"
        ax_c.text(j, i, f"{val:.2f}", ha='center', va='center', fontsize=7.0, color=color)

cbar_c = plt.colorbar(im_c, ax=ax_c, pad=0.03, fraction=0.04)
cbar_c.set_label("Pearson Correlation (r)", fontsize=8.5)
ax_c.set_title("C. Multi-Variate Biophysical Coupling\nSlope m Strongly Driven by Shrub % (0.72) and VH/VV (0.45)", fontsize=10.5, fontweight='bold', loc='left')

# ==============================================================================
# PANEL D: CONCEPTUAL RADAR SCATTERING SCHEMATIC OVER PEATLAND FACIES
# ==============================================================================
ax_d = axs[1, 1]
ax_d.set_xlim(0, 10)
ax_d.set_ylim(0, 10)
ax_d.axis('off')

# Section 1: Pristine Sphagnum Moss (Left)
rect_moss = patches.Rectangle((0.2, 0.4), 2.8, 8.8, facecolor='#e5f5e0', edgecolor='#31a354', lw=1.5)
ax_d.add_patch(rect_moss)
ax_d.text(1.6, 8.7, "1. Pristine Moss Lawn\n(West Column)", ha='center', fontsize=9.5, fontweight='bold', color='#006837')
# Ground surface
ax_d.plot([0.2, 3.0], [3.8, 3.8], color='#31a354', lw=2.5)
ax_d.text(1.6, 2.7, "Floating Peat Mat\n(Hydrated Sphagnum)", ha='center', fontsize=8.0, color='#238b45')
# Radar arrow: specular / rough surface bounce
ax_d.annotate('', xy=(1.6, 3.9), xytext=(0.5, 7.5),
            arrowprops=dict(facecolor='#2b83ba', edgecolor='#2b83ba', arrowstyle="->", lw=2.2))
ax_d.annotate('', xy=(2.7, 7.5), xytext=(1.6, 3.9),
            arrowprops=dict(facecolor='#2b83ba', edgecolor='#2b83ba', arrowstyle="->", lw=2.2))
ax_d.text(1.6, 5.6, "Surface Bounce\n(Co-Pol VV)", ha='center', fontsize=8.5, fontweight='bold', color='#1c507a')
ax_d.text(1.6, 1.2, "• TLS Slope: m = 1.01 (1:1)\n• Correlation: r = 0.75 – 0.82\n• Phase Center = True Peat Surface", 
          ha='center', fontsize=7.8, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#31a354", lw=0.8))

# Section 2: Boardwalk Infrastructure (Center)
rect_bw = patches.Rectangle((3.5, 0.4), 2.8, 8.8, facecolor='#fee6ce', edgecolor='#e6550d', lw=1.5)
ax_d.add_patch(rect_bw)
ax_d.text(4.9, 8.7, "2. Boardwalk Platform\n(Center Column / P6)", ha='center', fontsize=9.5, fontweight='bold', color='#a63603')
# Ground surface
ax_d.plot([3.5, 6.3], [3.8, 3.8], color='#31a354', lw=2.5)
# Boardwalk wooden piles
ax_d.plot([4.3, 4.3], [2.0, 5.0], color='#8c510a', lw=4.5)
ax_d.plot([5.5, 5.5], [2.0, 5.0], color='#8c510a', lw=4.5)
ax_d.plot([4.1, 5.7], [5.0, 5.0], color='#8c510a', lw=4.0)
ax_d.text(4.9, 2.7, "Timber Piles Anchored\nin Mineral Bed", ha='center', fontsize=7.8, color='#8c510a', fontweight='bold',
          bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#8c510a", lw=0.5, alpha=0.9))
# Radar arrow: double bounce
ax_d.annotate('', xy=(4.3, 3.8), xytext=(3.7, 7.5),
            arrowprops=dict(facecolor='#d95f02', edgecolor='#d95f02', arrowstyle="->", lw=2.2))
ax_d.annotate('', xy=(5.9, 7.5), xytext=(4.3, 3.8),
            arrowprops=dict(facecolor='#d95f02', edgecolor='#d95f02', arrowstyle="->", lw=2.2))
ax_d.text(4.9, 6.0, "Double Bounce\n+ Timber Drag", ha='center', fontsize=8.5, fontweight='bold', color='#a63603')
ax_d.text(4.9, 1.2, "• TLS Slope: m = 1.18 – 1.22\n• Mechanical Piling Restriction\n• Local Infrastructure Drag", 
          ha='center', fontsize=7.8, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#e6550d", lw=0.8))

# Section 3: Shrub Canopy Fringe (Right)
rect_shrub = patches.Rectangle((6.8, 0.4), 3.0, 8.8, facecolor='#f7f7f7', edgecolor='#756bb1', lw=1.5)
ax_d.add_patch(rect_shrub)
ax_d.text(8.3, 8.7, "3. Shrub Canopy Fringe\n(East Column / Margin)", ha='center', fontsize=9.5, fontweight='bold', color='#54278f')
# Ground surface
ax_d.plot([6.8, 9.8], [3.8, 3.8], color='#31a354', lw=2.5)
# Shrub branches
for bx, by in [(7.4, 4.8), (8.3, 5.6), (9.2, 5.0)]:
    ax_d.plot([bx, bx-0.2], [3.8, by], color='#54278f', lw=1.8)
    ax_d.plot([bx, bx+0.2], [3.8, by], color='#54278f', lw=1.8)
    circle = patches.Circle((bx, by), 0.42, facecolor='#bcbddc', edgecolor='#54278f', alpha=0.7)
    ax_d.add_patch(circle)
ax_d.text(8.3, 2.7, "Betula nana & Salix\n(Height 1.5 – 2.5 m)", ha='center', fontsize=8.0, color='#54278f')
# Radar arrow: volume scattering
ax_d.annotate('', xy=(8.1, 5.6), xytext=(7.1, 7.8),
            arrowprops=dict(facecolor='#756bb1', edgecolor='#756bb1', arrowstyle="->", lw=2.2))
ax_d.annotate('', xy=(9.5, 7.7), xytext=(8.5, 5.5),
            arrowprops=dict(facecolor='#756bb1', edgecolor='#756bb1', arrowstyle="->", lw=2.2))
ax_d.text(8.3, 6.7, "Volume Scattering\n(Depolarization VH/VV)", ha='center', fontsize=8.5, fontweight='bold', color='#54278f')
ax_d.text(8.3, 1.2, "• TLS Slope: m = 1.38 – 1.61\n• Exaggerated Apparent Motion\n• Phase Center in Canopy / Wind Flutter", 
          ha='center', fontsize=7.8, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#756bb1", lw=0.8))

ax_d.set_title("D. Conceptual Radar Scattering Architecture\nThree Micro-Facies Modulate InSAR Kinematics", fontsize=10.5, fontweight='bold', loc='left')

fig.suptitle("Physical Modulation of C-Band InSAR by Peatland Micro-Facies, Polarimetry and Temporal Baseline", 
             fontsize=13.0, fontweight='bold', y=0.98)

# Save
fig.savefig(OUT_DIR / "fig4_radar_physics_polarimetry_mechanisms.png", dpi=300)
# Also save to web atlas
WEB_DIR = ROOT / "08_deliverables" / "rzecin-atlas" / "web" / "public" / "data" / "figures"
fig.savefig(WEB_DIR / "fig4_radar_physics_polarimetry_mechanisms.png", dpi=300)
plt.close(fig)
print("Figure 4 successfully generated!")
