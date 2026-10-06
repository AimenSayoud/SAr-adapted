"""Generate 4-Pixel 1:1 Laser Kinematic Regression Comparison.
Compares Sentinel-1 12-day InSAR displacement against in-situ laser ground truth
across 4 contrasting biophysical pixels around Station P6:
1. West Benchmark (0% wood, 7.5% veg) -> Pure Moss Lawn (m = 1.01)
2. P6 Center (8% wood, 13.9% veg) -> Platform & Station (m = 1.18)
3. East Neighbour (0% wood, 21.0% veg) -> Shrub Fringe (m = 1.38)
4. North / P7 (2.6% wood, 25.9% veg) -> Dense Shrub Canopy (m = 1.54)
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
OUT_DIR = ROOT / "08_deliverables" / "supervisor_defense"
PAPER_DIR = ROOT / "08_deliverables" / "paper_figures"
WEB_DIR = ROOT / "08_deliverables" / "rzecin-atlas" / "web" / "public" / "data" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)
PAPER_DIR.mkdir(parents=True, exist_ok=True)
WEB_DIR.mkdir(parents=True, exist_ok=True)

# 1. Load data
p6_pairs = pd.read_csv(ROOT / "08_deliverables" / "field_p6_short_pairs_x065" / "p6_pairs_2020_2024.csv")
plots = pd.read_csv(ROOT / "08_deliverables" / "field_plots_x061" / "plots_pairs.csv")
ds = xr.open_dataset(ROOT / "06_data" / "cube" / "silver" / "ifg_2022_2024_ascending.nc")

p6_x, p6_y = 588380.0, 5846540.0
lam = 55.465763

valid = p6_pairs[(p6_pairs.track == 'ascending') & (p6_pairs.dt == 12) & p6_pairs.laser_los.notna() & (p6_pairs.step == 1)].copy()
pairs_in_ds = [p for p in valid.pair.values if p in ds.pair.values]
valid = valid[valid.pair.isin(pairs_in_ds)].copy()

# Extract West, Center, East, North
p6_off = plots[(plots['plot'].astype(str) == 'P6') & (plots['extraction'].astype(str) == 'off-boardwalk')]
p7_cell = plots[(plots['plot'].astype(str) == 'P7') & (plots['extraction'].astype(str) == 'cell')]

valid = valid.merge(p6_off[['pair', 's1_dlos_mm']], on='pair', suffixes=('', '_west'))
valid = valid.merge(p7_cell[['pair', 's1_dlos_mm']], on='pair', suffixes=('', '_north'))

# Extract East neighbour directly from ifg
east_s1 = []
for p in valid.pair.values:
    sub_row = valid[valid.pair == p].iloc[0]
    p_da = ds.sel(pair=p)
    u_p6 = float(p_da.unw.sel(x=p6_x, y=p6_y, method='nearest').values)
    u_east = float(p_da.unw.sel(x=p6_x + 40, y=p6_y, method='nearest').values)
    u_ref = u_p6 + sub_row.s1 * (4 * np.pi) / lam
    s1_e = -(u_east - u_ref) * lam / (4 * np.pi)
    east_s1.append(s1_e)
valid['s1_east'] = east_s1

# TLS slope function
def calc_metrics(x, y):
    r = np.corrcoef(x, y)[0, 1]
    p_ols = np.polyfit(x, y, 1)
    
    # PCA / TLS
    xc = x - x.mean()
    yc = y - y.mean()
    cov = np.cov(xc, yc)
    eigvals, eigvecs = np.linalg.eigh(cov)
    tls_m = eigvecs[1, 1] / eigvecs[0, 1]
    tls_b = y.mean() - tls_m * x.mean()
    
    sign_agree = (np.sign(x) == np.sign(y)).mean() * 100
    return r, p_ols[0], p_ols[1], tls_m, tls_b, sign_agree

# 2. Setup Figure
fig, axes = plt.subplots(2, 2, figsize=(14, 12), sharex=True, sharey=True, dpi=300)
lim = 28

pixel_configs = [
    {
        'ax': axes[0, 0],
        'title': 'A. West Benchmark Cell (Row 68, Col 65)\nPure Pristine Sphagnum Moss Lawn (0% Wood, 7.5% Veg)',
        'y_data': valid['s1_dlos_mm'].values,
        'color': '#1b7837',
        'edge': '#0a3314',
        'label_tag': 'West Moss Benchmark'
    },
    {
        'ax': axes[0, 1],
        'title': 'B. Station P6 Center Cell (Row 68, Col 66)\nFlux Platform & Pilings (8.0% Wood, 13.9% Veg)',
        'y_data': valid['s1'].values,
        'color': '#d95f02',
        'edge': '#4d1e00',
        'label_tag': 'P6 Center (Boardwalk)'
    },
    {
        'ax': axes[1, 0],
        'title': 'C. East Neighbour Cell (Row 68, Col 67)\nDwarf Birch Shrub Margin (0% Wood, 21.0% Veg)',
        'y_data': valid['s1_east'].values,
        'color': '#b2182b',
        'edge': '#4a0811',
        'label_tag': 'East Shrub Margin'
    },
    {
        'ax': axes[1, 1],
        'title': 'D. North Cell / Plot P7 (Row 67, Col 66)\nDense Birch Canopy + Walkway (2.6% Wood, 25.9% Veg)',
        'y_data': valid['s1_dlos_mm_north'].values,
        'color': '#762a83',
        'edge': '#300f38',
        'label_tag': 'North Dense Canopy'
    }
]

x_laser = valid['laser_los'].values

for cfg in pixel_configs:
    ax = cfg['ax']
    y_val = cfg['y_data']
    col = cfg['color']
    edg = cfg['edge']
    
    # 1:1 Reference Line & Bands
    ax.plot([-lim, lim], [-lim, lim], 'k--', lw=1.6, label='1 : 1 Physical Motion Line (y = x)', zorder=2)
    ax.axvspan(-3.5, 3.5, color='#e0e0e0', alpha=0.6, label='Laser Noise Floor (±3.5 mm)', zorder=1)
    ax.axhspan(-14, 14, color='#e0ecf4', alpha=0.45, label='±λ/4 InSAR Limit (±14 mm)', zorder=1)
    ax.axvline(0, color='gray', lw=0.6, ls=':', zorder=1)
    ax.axhline(0, color='gray', lw=0.6, ls=':', zorder=1)
    
    r_val, ols_m, ols_b, tls_m, tls_b, sign_agree = calc_metrics(x_laser, y_val)
    
    # Scatter points
    ax.scatter(x_laser, y_val, color=col, s=65, alpha=0.88, edgecolors=edg, lw=1.1, zorder=5,
               label=f'{cfg["label_tag"]} (n = {len(x_laser)})\n'
                     f'• r = {r_val:.2f} (Pearson)\n'
                     f'• TLS Slope: m = {tls_m:.2f} (Target 1.00)\n'
                     f'• OLS Slope: m = {ols_m:.2f}\n'
                     f'• Sign Agreement: {sign_agree:.0f}%')
    
    # Regression line
    x_grid = np.linspace(-18, 18, 50)
    ax.plot(x_grid, tls_m * x_grid + tls_b, color=col, lw=2.0, ls='-', zorder=4)
    
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_title(cfg['title'], fontsize=10.5, fontweight='bold', loc='left')
    ax.grid(True, alpha=0.25, ls=':')
    ax.legend(loc='upper left', fontsize=8.2, framealpha=0.92)

# Axis labels
axes[1, 0].set_xlabel('Ground Laser Displacement in LOS (mm, Ground Truth)', fontsize=10, fontweight='bold')
axes[1, 1].set_xlabel('Ground Laser Displacement in LOS (mm, Ground Truth)', fontsize=10, fontweight='bold')
axes[0, 0].set_ylabel('Sentinel-1 InSAR Displacement (mm LOS)', fontsize=10, fontweight='bold')
axes[1, 0].set_ylabel('Sentinel-1 InSAR Displacement (mm LOS)', fontsize=10, fontweight='bold')

fig.suptitle('Empirical 1:1 Kinematic Validation of Sentinel-1 12-Day InSAR Across 4 Contrasting Radar Pixels\n'
             'Demonstrating Exact 1:1 Moss Breathing vs. Boardwalk Anchor Drag and Canopy Volume Scattering Inflation',
             fontsize=12.5, fontweight='bold', y=0.99)

plt.tight_layout()

# Save
out_p1 = OUT_DIR / "fig2_laser_1to1_regression_4pixels_comparison.png"
out_p2 = PAPER_DIR / "Fig2_Laser_1to1_Kinematic_Regressions.png"
out_p3 = WEB_DIR / "fig2_laser_1to1_regression_4pixels_comparison.png"

fig.savefig(out_p1, dpi=300)
fig.savefig(out_p2, dpi=300)
fig.savefig(out_p3, dpi=300)
plt.close(fig)

print(f"Successfully generated 4-pixel comparison figure:\n- {out_p1}\n- {out_p2}\n- {out_p3}")
