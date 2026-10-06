"""Generate publication-grade figures for supervisor defense and paper manuscript:
Figure 1: 5-Year Synchronized Environmental & InSAR Timeline (2020–2024).
Figure 2: The 1:1 Motion Verification Grid (Ascending & Descending × 6-Day & 12-Day).
Figure 3: Spatial Evidence & Boardwalk-Free Radar Pixel Map (Orthophoto + 40m Grid).
Figure 4: Physical Drivers: Rain Event Buoyancy Step & Dawn Dew Decorrelation Mechanism.
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
p6_pairs = pd.read_csv(ROOT / "08_deliverables" / "field_p6_short_pairs_x065" / "p6_pairs_2020_2024.csv")
plots = pd.read_csv(ROOT / "08_deliverables" / "field_plots_x061" / "plots_pairs.csv")
p6_web = pd.read_csv(ROOT / "08_deliverables" / "field_p6" / "p6_web_series.csv")
wtd_meteo = pd.read_csv(ROOT / "06_data" / "field" / "raw" / "2026-09_abdallah_delivery" / "WTD_hourly_2020-2024_9plots_filled_meteo.csv")
p6_win = pd.read_csv(ROOT / "08_deliverables" / "field_boardwalk_x060" / "p6_windows.csv")
rain_ep = pd.read_csv(ROOT / "08_deliverables" / "field_events_x055" / "epoch_response.csv")
wetness = pd.read_csv(ROOT / "08_deliverables" / "field_dew_x050" / "wetness_at_overpasses_2020_2024.csv")

# Clean strings
p6_pairs['track'] = p6_pairs['track'].astype(str).str.strip()
p6_pairs['pair'] = p6_pairs['pair'].astype(str).str.strip()
plots['track'] = plots['track'].astype(str).str.strip()
plots['pair'] = plots['pair'].astype(str).str.strip()

p6_off = plots[(plots['plot'].astype(str) == 'P6') & (plots['extraction'].astype(str) == 'off-boardwalk')]
m = p6_pairs.merge(p6_off[['track', 'pair', 's1_dlos_mm', 'coh', 'dwtd_cm']], on=['track', 'pair'], suffixes=('', '_off'))
m_laser = m[(m.step == 1) & m.laser_los.notna()].copy()

print("Data loaded successfully!")

# ==============================================================================
# FIGURE 1: 5-YEAR SYNCHRONIZED ENVIRONMENTAL & InSAR TIMELINE (2020–2024)
# ==============================================================================
print("Generating Figure 1...")
fig, axes = plt.subplots(4, 1, figsize=(15, 13), sharex=True, 
                         gridspec_kw={'height_ratios': [1.1, 1.4, 1.0, 1.1]})

wtd_meteo['dt'] = pd.to_datetime(wtd_meteo['TIMESTAMP'])
daily_rain = wtd_meteo.set_index('dt')['Rain_mm_Tot'].resample('D').sum()
daily_wtd_p6 = wtd_meteo.set_index('dt')['WTD_P6 (CR)'].resample('D').median()
daily_wtd_p1 = wtd_meteo.set_index('dt')['WTD_P1'].resample('D').median()

# Panel A: Rain & WTD
ax0_rain = axes[0].twinx()
ax0_rain.bar(daily_rain.index, daily_rain.values, color='#7ec8e3', alpha=0.45, width=1.2, label='Daily Rain (mm)')
ax0_rain.set_ylabel('Rain (mm/day)', color='#006699', fontsize=9.5)
ax0_rain.set_ylim(0, 55)
ax0_rain.tick_params(axis='y', labelcolor='#006699')

axes[0].plot(daily_wtd_p6.index, daily_wtd_p6.values, color='#08306b', lw=1.5, label='WTD P6 (Floating Mat Central, cm)')
axes[0].plot(daily_wtd_p1.index, daily_wtd_p1.values, color='#d95f02', lw=1.2, ls='--', label='WTD P1 (Fixed Mineral Edge, cm)')
axes[0].set_ylabel('Water Table Depth (cm)', fontsize=10)
axes[0].set_ylim(-75, 5)
axes[0].set_title('A. Hydrological Forcing: Daily Precipitation and Water Table Dynamics across Transect (2020–2024)', fontsize=11, fontweight='bold', loc='left')
axes[0].grid(True, alpha=0.25, ls=':')
lines_0, labels_0 = axes[0].get_legend_handles_labels()
lines_0r, labels_0r = ax0_rain.get_legend_handles_labels()
axes[0].legend(lines_0 + lines_0r, labels_0 + labels_0r, loc='lower left', fontsize=8.5, ncol=3, framealpha=0.9)

# Panel B: Peat Surface Motion (Laser vs S1 P6 vs S1 Neighbour)
p6_asc = p6_web[(p6_web.track == 'ascending') & (p6_web.window == '1x1')].copy()
p6_asc['dt'] = pd.to_datetime(p6_asc['date'])

axes[1].plot(p6_asc['dt'], p6_asc['laser_los_matched_mm'], 'o-', color='#1a1a1a', ms=4, lw=1.6, label='P6 Ground Truth Laser Surface (LOS mm, snow-masked)')
axes[1].plot(p6_asc['dt'], p6_asc['s1_chain_mm'], 'd--', color='#e65c00', ms=4.5, lw=1.3, label='Sentinel-1 InSAR: P6 Cell (with Wooden Platform, LOS mm)')

m_asc_12 = m[(m.track == 'ascending') & (m.dt == 12)].copy()
m_asc_12['dt'] = pd.to_datetime(m_asc_12['pair'].str.split('_').str[1])
m_asc_12 = m_asc_12.sort_values('dt')
m_asc_12['s1_off_cum'] = m_asc_12['s1_dlos_mm'].cumsum()
shift_off = p6_asc['laser_los_matched_mm'].dropna().mean() - m_asc_12['s1_off_cum'].mean()
axes[1].plot(m_asc_12['dt'], m_asc_12['s1_off_cum'] + shift_off, '^-.', color='#1b7837', ms=4, lw=1.3,
             label='Sentinel-1 InSAR: Boardwalk-Free Neighbour (0% Wood, LOS mm)')

axes[1].axhspan(-14, 14, color='#e0ecf4', alpha=0.6, label='Unambiguous Phase Boundary (±λ/4 = ±14 mm LOS)')
axes[1].set_ylabel('Displacement (mm LOS)', fontsize=10)
axes[1].set_title('B. Peat Surface Dynamics: InSAR Displacement Tracking vs. Verified Laser Ground Truth', fontsize=11, fontweight='bold', loc='left')
axes[1].grid(True, alpha=0.25, ls=':')
axes[1].legend(loc='lower left', fontsize=8.5, ncol=2, framealpha=0.9)

# Panel C: Consecutive Pair Displacements
m_consec = m[(m.track == 'ascending') & (m.step == 1)].copy()
m_consec['dt'] = pd.to_datetime(m_consec['pair'].str.split('_').str[1])
m_consec = m_consec.sort_values('dt')

axes[2].bar(m_consec['dt'] - pd.Timedelta(days=2), m_consec['s1'], width=3, color='#fdb863', alpha=0.9, label='ΔLOS: P6 Cell (with wood, mm)')
axes[2].bar(m_consec['dt'] + pd.Timedelta(days=2), m_consec['s1_dlos_mm'], width=3, color='#80cdc1', alpha=0.9, label='ΔLOS: Boardwalk-Free Neighbour (0% wood, mm)')
axes[2].axhline(0, color='gray', lw=0.8)
axes[2].set_ylabel('ΔLOS per Pair (mm)', fontsize=10)
axes[2].set_ylim(-25, 25)
axes[2].set_title('C. Step-by-Step Consecutive Pair Displacements: P6 Cell vs. Wood-Free Neighbour', fontsize=11, fontweight='bold', loc='left')
axes[2].grid(True, alpha=0.25, ls=':')
axes[2].legend(loc='upper right', fontsize=8.5, ncol=2, framealpha=0.9)

# Panel D: Coherence across eras
m_sorted = m[(m.track == 'ascending') & (m.step == 1)].copy()
m_sorted['dt'] = pd.to_datetime(m_sorted['pair'].str.split('_').str[1])
m_sorted = m_sorted.sort_values('dt')

axes[3].plot(m_sorted['dt'], m_sorted['coh'], 'o-', color='#2b83ba', ms=3.5, lw=1.1, label='P6 Cell Coherence (γ)')
axes[3].plot(m_sorted['dt'], m_sorted['coh_off'], 's--', color='#2ca25f', ms=3.5, lw=1.1, label='Boardwalk-Free Neighbour Coherence (γ)')
axes[3].axhline(0.40, color='#d7191c', ls=':', lw=1.2, label='Coherence Quality Floor (γ = 0.40)')

# Shading for 6-day era vs 12-day era
axes[3].axvspan(pd.Timestamp('2020-01-01'), pd.Timestamp('2021-12-31'), color='#ffffbf', alpha=0.25, label='6-Day Revisit Era (Sentinel-1A + 1B)')
axes[3].axvspan(pd.Timestamp('2022-01-01'), pd.Timestamp('2024-12-31'), color='#fee090', alpha=0.15, label='12-Day Revisit Era (Sentinel-1A only)')

axes[3].set_ylabel('Coherence (γ)', fontsize=10)
axes[3].set_ylim(0.15, 1.02)
axes[3].set_title('D. Sentinel-1 Interferometric Coherence: High 6-Day Era vs. Vegetative Decorrelation in 12-Day Era', fontsize=11, fontweight='bold', loc='left')
axes[3].grid(True, alpha=0.25, ls=':')
axes[3].legend(loc='lower left', fontsize=8.5, ncol=3, framealpha=0.9)

# Align x limits
axes[3].set_xlim(pd.Timestamp('2020-01-01'), pd.Timestamp('2024-12-31'))

fig.tight_layout()
fig.savefig(OUT_DIR / "fig1_annual_timeseries_wtd_laser_insar_coh.png", dpi=300)
fig.savefig(WEB_DIR / "fig1_annual_timeseries_wtd_laser_insar_coh.png", dpi=300)
plt.close(fig)
print("Figure 1 generated at 300 DPI!")

# ==============================================================================
# FIGURE 2: THE 1:1 MOTION VERIFICATION GRID (Ascending/Descending × 6-day/12-day)
# ==============================================================================
print("Generating Figure 2...")
fig, axes = plt.subplots(2, 2, figsize=(14, 12), sharex=True, sharey=True)

lim = 28
configs = [
    (0, 0, "ascending", 6, "A. Ascending (Dusk) — 6-Day Pairs (2020–2021, S1A+S1B)"),
    (0, 1, "ascending", 12, "B. Ascending (Dusk) — 12-Day Pairs (2022–2024, S1A only)"),
    (1, 0, "descending", 6, "C. Descending (Dawn) — 6-Day Pairs (2020–2021, S1A+S1B)"),
    (1, 1, "descending", 12, "D. Descending (Dawn) — 12-Day Pairs (2022–2024, S1A only)")
]

for r_idx, c_idx, track, dt_target, title in configs:
    ax = axes[r_idx, c_idx]
    sub = m_laser[(m_laser.dt == dt_target) & (m_laser.track == track)].copy()
    
    # 1:1 Reference Line & Shading
    ax.plot([-lim, lim], [-lim, lim], 'k--', lw=1.6, label='1 : 1 Physical Motion (y = x)')
    ax.axvspan(-3.5, 3.5, color='#e0e0e0', alpha=0.7, label='Laser Noise Floor (±3.5 mm)')
    ax.axhspan(-14, 14, color='#e0ecf4', alpha=0.45, label='±λ/4 InSAR Limit (±14 mm)')
    ax.axvline(0, color='gray', lw=0.6, ls=':')
    ax.axhline(0, color='gray', lw=0.6, ls=':')
    
    # P6 Center Cell
    r_p6 = np.corrcoef(sub.laser_los, sub.s1)[0, 1]
    poly_p6 = np.polyfit(sub.laser_los, sub.s1, 1)
    sign_p6 = (np.sign(sub.laser_los) == np.sign(sub.s1)).mean() * 100
    ax.scatter(sub.laser_los, sub.s1, color='#d95f02', s=55, alpha=0.85, edgecolors='#4d1e00', zorder=4,
               label=f'P6 Cell (8% Wood): r = {r_p6:.2f}, slope = {poly_p6[0]:.2f}, sign = {sign_p6:.0f}%')
    x_vals = np.linspace(sub.laser_los.min(), sub.laser_los.max(), 50)
    ax.plot(x_vals, np.polyval(poly_p6, x_vals), color='#d95f02', lw=1.8, ls='-')
    
    # Wood-Free Neighbour
    r_off = np.corrcoef(sub.laser_los, sub.s1_dlos_mm)[0, 1]
    poly_off = np.polyfit(sub.laser_los, sub.s1_dlos_mm, 1)
    sign_off = (np.sign(sub.laser_los) == np.sign(sub.s1_dlos_mm)).mean() * 100
    ax.scatter(sub.laser_los, sub.s1_dlos_mm, color='#1b7837', s=50, marker='^', alpha=0.85, edgecolors='#0a3314', zorder=5,
               label=f'Wood-Free Neighbour (0% Wood): r = {r_off:.2f}, slope = {poly_off[0]:.2f}, sign = {sign_off:.0f}%')
    ax.plot(x_vals, np.polyval(poly_off, x_vals), color='#1b7837', lw=1.8, ls='--')
    
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_title(f"{title} (n = {len(sub)})", fontsize=10.5, fontweight='bold', loc='left')
    ax.grid(True, alpha=0.25, ls=':')
    ax.legend(loc='upper left', fontsize=8.0, framealpha=0.92)
    
    if r_idx == 1:
        ax.set_xlabel('Laser Displacement in LOS (mm, Ground Truth)', fontsize=9.5)
    if c_idx == 0:
        ax.set_ylabel('Sentinel-1 InSAR Displacement (mm LOS)', fontsize=9.5)

fig.suptitle('Validation of Sentinel-1 InSAR Peat Mat Motion Against In-Situ Laser (2020–2024)\nTesting Spatial Confounding by Wooden Boardwalk Structures',
             fontsize=12, fontweight='bold', y=0.99)
fig.tight_layout()
fig.savefig(OUT_DIR / "fig2_laser_1to1_regression_p6_vs_neighbours.png", dpi=300)
fig.savefig(WEB_DIR / "fig2_laser_1to1_regression_p6_vs_neighbours.png", dpi=300)
plt.close(fig)
print("Figure 2 generated at 300 DPI!")

# ==============================================================================
# FIGURE 3: SPATIAL EVIDENCE & BOARDWALK-FREE RADAR PIXEL MAP
# ==============================================================================
print("Generating Figure 3...")
fig = plt.figure(figsize=(15, 8.5))
gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1.0])

# Left: Orthophoto & 40m Grid Overlay around P6
ax_map = fig.add_subplot(gs[0])
with rasterio.open(ROOT / "08_deliverables" / "field_boardwalk_x060" / "orthophoto_transect.tif") as src:
    img = src.read()
    # Crop around P6
    p6_e, p6_n = 588380.0, 5846540.0
    crop_w = 140
    ext = (p6_e - crop_w, p6_e + crop_w, p6_n - crop_w, p6_n + crop_w)
    
    # Read subwindow
    window = rasterio.windows.from_bounds(*ext, src.transform)
    sub_img = src.read(window=window)
    sub_img = np.moveaxis(sub_img, 0, -1)
    ax_map.imshow(sub_img, extent=[ext[0], ext[1], ext[2], ext[3]])

# Draw 40m radar grid cells around P6
for x in np.arange(p6_e - 100, p6_e + 140, 40):
    ax_map.axvline(x, color='yellow', lw=1.2, ls='--', alpha=0.85)
for y in np.arange(p6_n - 100, p6_n + 140, 40):
    ax_map.axhline(y, color='yellow', lw=1.2, ls='--', alpha=0.85)

# Annotate Cells
# P6 cell
p6_box = patches.Rectangle((p6_e - 20, p6_n - 20), 40, 40, linewidth=2.5, edgecolor='#d95f02', facecolor='none', zorder=10)
ax_map.add_patch(p6_box)
ax_map.text(p6_e, p6_n + 12, "P6 Center Cell\n(8.0% Wood)", color='white', fontweight='bold', fontsize=8.5, ha='center',
            bbox=dict(boxstyle="round,pad=0.2", fc="#d95f02", alpha=0.85))

# West cell (Left)
w_box = patches.Rectangle((p6_e - 60, p6_n - 20), 40, 40, linewidth=2.5, edgecolor='#1b7837', facecolor='none', zorder=10)
ax_map.add_patch(w_box)
ax_map.text(p6_e - 40, p6_n + 12, "West Neighbour\n(0.0% Wood)", color='white', fontweight='bold', fontsize=8.5, ha='center',
            bbox=dict(boxstyle="round,pad=0.2", fc="#1b7837", alpha=0.85))

# East cell (Right)
e_box = patches.Rectangle((p6_e + 20, p6_n - 20), 40, 40, linewidth=2.5, edgecolor='#1b7837', facecolor='none', zorder=10)
ax_map.add_patch(e_box)
ax_map.text(p6_e + 40, p6_n + 12, "East Neighbour\n(0.0% Wood)", color='white', fontweight='bold', fontsize=8.5, ha='center',
            bbox=dict(boxstyle="round,pad=0.2", fc="#1b7837", alpha=0.85))

ax_map.set_xlim(p6_e - 100, p6_e + 100)
ax_map.set_ylim(p6_n - 80, p6_n + 80)
ax_map.set_xlabel("UTM Easting (m, EPSG:32633)", fontsize=9.5)
ax_map.set_ylabel("UTM Northing (m, EPSG:32633)", fontsize=9.5)
ax_map.set_title("A. High-Resolution Orthophoto (0.25 m) & 40 m Radar Grid\nProof of Boardwalk-Free Neighbour Pixels", fontsize=11, fontweight='bold', loc='left')

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
ax_bars.set_xlabel("Metric Value", fontsize=10)
ax_bars.set_xlim(0, 1.6)
ax_bars.set_title("B. InSAR-Laser Agreement Across 13 Spatial Extraction Windows\nWood-Free Pixels Match or Exceed P6 Center", fontsize=11, fontweight='bold', loc='left')
ax_bars.grid(True, alpha=0.25, ls=':')
ax_bars.legend(loc='lower right', fontsize=8.5)

# Value annotations
for i, (r_v, t_v) in enumerate(zip(r_vals, tls_vals)):
    ax_bars.text(r_v + 0.02, i - width/2, f"r={r_v:.2f}", va='center', fontsize=8.5, fontweight='bold')
    ax_bars.text(t_v + 0.02, i + width/2, f"m={t_v:.2f}", va='center', fontsize=8.5)

fig.tight_layout()
fig.savefig(OUT_DIR / "fig3_spatial_evidence_boardwalk_free_pixels.png", dpi=300)
fig.savefig(WEB_DIR / "fig3_spatial_evidence_boardwalk_free_pixels.png", dpi=300)
plt.close(fig)
print("Figure 3 generated at 300 DPI!")

# ==============================================================================
# FIGURE 4: PHYSICAL DRIVERS: RAIN STEP-RESPONSE & DEW DECORRELATION
# ==============================================================================
print("Generating Figure 4...")
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: Rain Event Superposed Epoch (44 Storms)
ax_rain = axes[0]
wtd_ep = rain_ep[rain_ep.series == 'WTD plot median']
laser_ep = rain_ep[rain_ep.series == 'laser surface P6']
p1_ep = rain_ep[rain_ep.series == 'WTD P1']

ax_rain.plot(wtd_ep['lag_h'], wtd_ep['median'], color='#08306b', lw=2.0, label='Median Water Table Rise (cm)')
ax_rain.fill_between(wtd_ep['lag_h'], wtd_ep['q1'], wtd_ep['q3'], color='#08306b', alpha=0.15)

ax_rain.plot(p1_ep['lag_h'], p1_ep['median'], color='#d95f02', lw=1.5, ls='--', label='P1 Mineral Edge WTD Rise (cm)')

ax_laser = ax_rain.twinx()
ax_laser.plot(laser_ep['lag_h'], laser_ep['median'] * 10, color='#1b7837', lw=2.2, label='P6 Laser Surface Rise (mm)')
ax_laser.fill_between(laser_ep['lag_h'], laser_ep['q1'] * 10, laser_ep['q3'] * 10, color='#1b7837', alpha=0.2)
ax_laser.set_ylabel('Laser Surface Displacement (mm)', color='#1b7837', fontsize=10)
ax_laser.tick_params(axis='y', labelcolor='#1b7837')
ax_laser.set_ylim(-2, 10)

ax_rain.axvline(0, color='gray', ls=':', lw=1.2, label='Rain Event Onset (t = 0)')
ax_rain.set_xlabel('Hours Relative to Rain Event Onset (h)', fontsize=9.5)
ax_rain.set_ylabel('Water Table Depth Rise (cm)', color='#08306b', fontsize=10)
ax_rain.tick_params(axis='y', labelcolor='#08306b')
ax_rain.set_xlim(-48, 240)
ax_rain.set_title('A. 44 Rain Storm Events (≥10 mm / 24 h):\nWater Table Peak vs. Floating Mat Buoyancy Step', fontsize=11, fontweight='bold', loc='left')
ax_rain.grid(True, alpha=0.25, ls=':')

lines_r1, labels_r1 = ax_rain.get_legend_handles_labels()
lines_r2, labels_r2 = ax_laser.get_legend_handles_labels()
ax_rain.legend(lines_r1 + lines_r2, labels_r1 + labels_r2, loc='lower right', fontsize=8.5)

# Right: Dawn vs Dusk Relative Humidity & Dew Contrast
ax_dew = axes[1]
wet_asc = wetness[wetness.track == 'ascending']['rh'].dropna()
wet_desc = wetness[wetness.track == 'descending']['rh'].dropna()

parts = ax_dew.violinplot([wet_asc, wet_desc], positions=[1, 2], showmeans=True, showmedians=True)
parts['bodies'][0].set_facecolor('#d95f02')
parts['bodies'][0].set_edgecolor('#d95f02')
parts['bodies'][0].set_alpha(0.65)
parts['bodies'][1].set_facecolor('#2b83ba')
parts['bodies'][1].set_edgecolor('#2b83ba')
parts['bodies'][1].set_alpha(0.65)

ax_dew.axhline(95, color='#d7191c', ls='--', lw=1.5, label='Liquid Dew / Canopy Saturation Threshold (RH ≥ 95%)')
ax_dew.set_xticks([1, 2])
ax_dew.set_xticklabels(['Ascending Track (Dusk)\n16:36 UTC (Dry Canopy)\n21% Wet Acquisitions', 
                        'Descending Track (Dawn)\n05:09 UTC (Dew Formed)\n46% Wet Acquisitions'], fontsize=9.5)
ax_dew.set_ylabel('Station Relative Humidity at Overpass (%)', fontsize=10)
ax_dew.set_ylim(30, 105)
ax_dew.set_title('B. Overpass Atmospheric Humidity:\nWhy 12-Day Dawn Track Experiences Severe Decorrelation', fontsize=11, fontweight='bold', loc='left')
ax_dew.grid(True, alpha=0.25, ls=':')
ax_dew.legend(loc='lower left', fontsize=8.5)

fig.tight_layout()
fig.savefig(OUT_DIR / "fig4_physical_drivers_rain_and_dew.png", dpi=300)
fig.savefig(WEB_DIR / "fig4_physical_drivers_rain_and_dew.png", dpi=300)
plt.close(fig)
print("Figure 4 generated at 300 DPI!")
print("All 4 paper-grade figures successfully generated!")
