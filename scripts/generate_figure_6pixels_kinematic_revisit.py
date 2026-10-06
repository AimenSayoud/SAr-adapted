"""Generate Definitive 6-Pixel Publication Figure: 6-Day vs. 12-Day Kinematic Tracking Against Laser.

Compares Sentinel-1 InSAR vs. Campbell SDMS40 ground truth laser across 6 contrasting biophysical pixels:
1. West Benchmark (0% Wood, 7.5% Veg) -> Pure Sphagnum Moss Lawn
2. South-West Replicate (0.6% Wood, 9.0% Veg) -> Open Moss Lawn Corridor
3. Station P6 Center (8.0% Wood, 13.9% Veg) -> EC Tower Platform & Pilings
4. South Walkway (7.0% Wood, 13.4% Veg) -> Timber Boardwalk Transect
5. East Margin (0% Wood, 21.0% Veg) -> Dwarf Birch Shrub Margin
6. North / P7 Canopy (2.6% Wood, 25.9% Veg) -> Dense Birch Canopy + Pilings

Each panel presents both 6-day (2021) and 12-day (2022-2024) seasonal pairs under snow-free conditions.
"""
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

# Publication styling
plt.rcParams['font.sans-serif'] = 'Helvetica', 'Arial', 'DejaVu Sans'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8

ROOT = Path("/Users/aymen/Documents/Research_Hub")
OUT_PAPER = ROOT / "08_deliverables" / "paper_figures"
OUT_DEFENSE = ROOT / "08_deliverables" / "supervisor_defense"
OUT_WEB = ROOT / "08_deliverables" / "rzecin-atlas" / "web" / "public" / "data" / "figures"

for d in [OUT_PAPER, OUT_DEFENSE, OUT_WEB]:
    d.mkdir(parents=True, exist_ok=True)

# 1. Load Data
p6_pairs = pd.read_csv(ROOT / "08_deliverables" / "field_p6_short_pairs_x065" / "p6_pairs_2020_2024.csv")
ds_21 = xr.open_dataset(ROOT / "06_data" / "cube" / "silver" / "ifg_2020_2021_ascending.nc")
ds_22 = xr.open_dataset(ROOT / "06_data" / "cube" / "silver" / "ifg_2022_2024_ascending.nc")

p6_x, p6_y = 588380.0, 5846540.0
lam = 55.465763

# Filter 6-day (2021 May-Nov) and 12-day (2022-2024 May-Nov)
sub_21 = p6_pairs[(p6_pairs.track == 'ascending') & (p6_pairs.dt == 6) & 
                  p6_pairs.laser_los.notna() & (p6_pairs.year == 2021)].copy()
pairs_21 = [p for p in sub_21.pair.values if p in ds_21.pair.values]
sub_21 = sub_21[sub_21.pair.isin(pairs_21)].copy()

sub_22 = p6_pairs[(p6_pairs.track == 'ascending') & (p6_pairs.dt == 12) & 
                  p6_pairs.laser_los.notna() & (p6_pairs.year >= 2022)].copy()
sub_22['d1'] = pd.to_datetime(sub_22['pair'].str.split('_').str[0])
sub_22['d2'] = pd.to_datetime(sub_22['pair'].str.split('_').str[1])
sub_22 = sub_22[(sub_22.d1.dt.month >= 5) & (sub_22.d2.dt.month <= 11)].copy()
pairs_22 = [p for p in sub_22.pair.values if p in ds_22.pair.values]
sub_22 = sub_22[sub_22.pair.isin(pairs_22)].copy()

# The 6 Selected Pixels
pixel_configs = [
    {
        'title': 'A. West Benchmark (Row 68, Col 65)\nPure Sphagnum Lawn (0% Wood, 7.5% Veg)',
        'dx': -40, 'dy': 0, 'wood': '0.0%', 'veg': '7.5%'
    },
    {
        'title': 'B. South-West Replicate (Row 69, Col 65)\nOpen Moss Mat Corridor (0.6% Wood, 9.0% Veg)',
        'dx': -40, 'dy': -40, 'wood': '0.6%', 'veg': '9.0%'
    },
    {
        'title': 'C. Station P6 Center (Row 68, Col 66)\nEC Platform & Pilings (8.0% Wood, 13.9% Veg)',
        'dx': 0, 'dy': 0, 'wood': '8.0%', 'veg': '13.9%'
    },
    {
        'title': 'D. South Walkway (Row 69, Col 66)\nBoardwalk Walkway (7.0% Wood, 13.4% Veg)',
        'dx': 0, 'dy': -40, 'wood': '7.0%', 'veg': '13.4%'
    },
    {
        'title': 'E. East Margin (Row 68, Col 67)\nDwarf Birch Shrub (0% Wood, 21.0% Veg)',
        'dx': 40, 'dy': 0, 'wood': '0.0%', 'veg': '21.0%'
    },
    {
        'title': 'F. North / P7 Canopy (Row 67, Col 66)\nDense Birch Canopy (2.6% Wood, 25.9% Veg)',
        'dx': 0, 'dy': 40, 'wood': '2.6%', 'veg': '25.9%'
    }
]

def extract_disp(sub_df, ds_obj, dx, dy):
    s1_vals = []
    for p in sub_df.pair.values:
        sub_row = sub_df[sub_df.pair == p].iloc[0]
        p_da = ds_obj.sel(pair=p)
        u_p6 = float(p_da.unw.sel(x=p6_x, y=p6_y, method='nearest').values)
        u_pix = float(p_da.unw.sel(x=p6_x + dx, y=p6_y + dy, method='nearest').values)
        u_ref = u_p6 + sub_row.s1 * (4 * np.pi) / lam
        s1_pix = -(u_pix - u_ref) * lam / (4 * np.pi)
        s1_vals.append(s1_pix)
    return np.array(s1_vals)

def tls_regression(x, y):
    r = np.corrcoef(x, y)[0, 1]
    p_ols = np.polyfit(x, y, 1)
    
    xc = x - x.mean()
    yc = y - y.mean()
    cov = np.cov(xc, yc)
    eigvals, eigvecs = np.linalg.eigh(cov)
    tls_m = eigvecs[1, 1] / eigvecs[0, 1]
    tls_b = y.mean() - tls_m * x.mean()
    sign_agree = (np.sign(x) == np.sign(y)).mean() * 100
    return r, p_ols[0], tls_m, tls_b, sign_agree

# ==============================================================================
# FIGURE CREATION: 2 ROWS x 3 COLUMNS (6 PANELS)
# ==============================================================================
fig, axes = plt.subplots(2, 3, figsize=(18, 12), sharex=True, sharey=True, dpi=300)
lim = 24

for idx, conf in enumerate(pixel_configs):
    r_idx = idx // 3
    c_idx = idx % 3
    ax = axes[r_idx, c_idx]
    dx, dy = conf['dx'], conf['dy']
    
    # Extract 6-day (2021) and 12-day (2022-2024) displacements
    s1_21 = extract_disp(sub_21, ds_21, dx, dy)
    s1_22 = extract_disp(sub_22, ds_22, dx, dy)
    
    r_21, ols_21, tls_m_21, tls_b_21, sign_21 = tls_regression(sub_21.laser_los.values, s1_21)
    r_22, ols_22, tls_m_22, tls_b_22, sign_22 = tls_regression(sub_22.laser_los.values, s1_22)
    
    coh_21 = sub_21.coh.median()
    coh_22 = sub_22.coh.median()
    
    # 1:1 Line & Boundaries
    ax.plot([-lim, lim], [-lim, lim], 'k--', lw=1.4, zorder=1, label='1 : 1 Physical Motion')
    ax.axhspan(-14, 14, color='#e0ecf4', alpha=0.35, zorder=0, label='±λ/4 InSAR Limit')
    ax.axvspan(-3.5, 3.5, color='#f0f0f0', alpha=0.6, zorder=0, label='Laser Noise Floor')
    ax.axvline(0, color='gray', lw=0.6, ls=':')
    ax.axhline(0, color='gray', lw=0.6, ls=':')
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.grid(True, ls=':', color='#cccccc', alpha=0.6)
    
    # Scatter 6-Day 2021 (Blue/Teal Squares)
    ax.scatter(sub_21.laser_los, s1_21, color='#2b83ba', marker='s', s=65, alpha=0.85,
               edgecolors='#08306b', zorder=4,
               label=f"6-Day (2021): r = {r_21:.2f}, TLS m = {tls_m_21:.2f}")
    
    # Scatter 12-Day 2022-2024 (Coral/Orange Circles)
    ax.scatter(sub_22.laser_los, s1_22, color='#d95f02', marker='o', s=70, alpha=0.9,
               edgecolors='#4d1e00', zorder=5,
               label=f"12-Day (2022–24): r = {r_22:.2f}, TLS m = {tls_m_22:.2f}")
    
    # Regression Lines
    x_line_21 = np.linspace(-15, 12, 50)
    ax.plot(x_line_21, tls_m_21 * x_line_21 + tls_b_21, color='#2b83ba', lw=1.8, ls='--', zorder=3)
    
    x_line_22 = np.linspace(-16, 12, 50)
    ax.plot(x_line_22, tls_m_22 * x_line_22 + tls_b_22, color='#d95f02', lw=2.0, zorder=3)
    
    ax.set_title(conf['title'], fontsize=11, fontweight='bold', pad=8)
    ax.legend(loc='upper left', fontsize=8.5, framealpha=0.92)
    
    # Bottom-right Metrics Box (empty quadrant)
    box_txt = (
        f"6-DAY REVISIT (2021, n=20):\n"
        f"  γ = {coh_21:.2f} | Sign: {sign_21:.0f}%\n"
        f"  TLS Slope: m = {tls_m_21:.2f}\n"
        f"12-DAY REVISIT (2022–24, n=26):\n"
        f"  γ = {coh_22:.2f} | Sign: {sign_22:.0f}%\n"
        f"  TLS Slope: m = {tls_m_22:.2f}"
    )
    ax.text(0.50, 0.04, box_txt, transform=ax.transAxes, fontsize=8.0,
            verticalalignment='bottom', bbox=dict(boxstyle='round,pad=0.4', facecolor='#f7f7f7', edgecolor='#666666', alpha=0.95))

for ax in axes[1, :]:
    ax.set_xlabel('In-Situ Laser Ground Truth in LOS (mm)', fontsize=11, fontweight='bold')
for ax in axes[:, 0]:
    ax.set_ylabel('Sentinel-1 InSAR Displacement (mm LOS)', fontsize=11, fontweight='bold')

fig.suptitle("Direct Biophysical Modulation of Radar Kinematics Across 6 Micro-Experiment Pixels\nComparing 6-Day Constellation (2021) vs. 12-Day Revisit (2022–2024) Against In-Situ Ground Laser",
             fontsize=13.5, fontweight='bold', y=0.99)
fig.tight_layout()

f_name = "Fig_6Pixels_Kinematic_Revisit_Comparison.png"
fig.savefig(OUT_PAPER / f_name, dpi=300)
fig.savefig(OUT_DEFENSE / f_name, dpi=300)
fig.savefig(OUT_WEB / f_name, dpi=300)
plt.close(fig)
print(f"6-Pixel Figure saved to {f_name}")
