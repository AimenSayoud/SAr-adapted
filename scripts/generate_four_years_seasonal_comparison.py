"""Generate Comprehensive 4-Year Seasonal Comparison (May-Nov 2021, 2022, 2023, 2024).

Compares Sentinel-1 InSAR against Campbell SDMS40 ground laser across 4 consecutive years
during the snow-free summer/autumn growth season (May 1 - November 30).
"""
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import pandas as pd

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
plots = pd.read_csv(ROOT / "08_deliverables" / "field_plots_x061" / "plots_pairs.csv")

p6_off = plots[(plots['plot'].astype(str) == 'P6') & (plots['extraction'].astype(str) == 'off-boardwalk')][['pair', 's1_dlos_mm']]
df = p6_pairs.merge(p6_off, on='pair', suffixes=('', '_off'))

df['d1'] = pd.to_datetime(df['pair'].str.split('_').str[0])
df['d2'] = pd.to_datetime(df['pair'].str.split('_').str[1])
df['mid_date'] = df['d1'] + (df['d2'] - df['d1']) / 2

def calc_stats(sub, y_col='s1'):
    r = np.corrcoef(sub.laser_los, sub[y_col])[0, 1]
    p_ols = np.polyfit(sub.laser_los, sub[y_col], 1)
    
    # PCA / TLS
    xc = sub.laser_los - sub.laser_los.mean()
    yc = sub[y_col] - sub[y_col].mean()
    cov = np.cov(xc, yc)
    eigvals, eigvecs = np.linalg.eigh(cov)
    tls_m = eigvecs[1, 1] / eigvecs[0, 1]
    tls_b = sub[y_col].mean() - tls_m * sub.laser_los.mean()
    
    sign_agree = (np.sign(sub.laser_los) == np.sign(sub[y_col])).mean() * 100
    med_coh = sub.coh.median()
    return {
        'n': len(sub),
        'r': r,
        'ols_m': p_ols[0],
        'ols_b': p_ols[1],
        'tls_m': tls_m,
        'tls_b': tls_b,
        'sign_agree': sign_agree,
        'med_coh': med_coh
    }

years = [2021, 2022, 2023, 2024]
revisit_dict = {2021: 6, 2022: 12, 2023: 12, 2024: 12}
subsets = {}
stats_p6 = {}
stats_off = {}

for yr in years:
    dt_target = revisit_dict[yr]
    sub = df[(df.d1.dt.year == yr) & (df.d1.dt.month >= 5) & (df.d2.dt.month <= 11) & 
             (df.dt == dt_target) & (df.track == 'ascending') & df.laser_los.notna()].copy()
    subsets[yr] = sub
    stats_p6[yr] = calc_stats(sub, 's1')
    stats_off[yr] = calc_stats(sub, 's1_dlos_mm')
    print(f"Year {yr} ({dt_target}d): n={len(sub)}, P6 TLS m={stats_p6[yr]['tls_m']:.2f}, Lawn TLS m={stats_off[yr]['tls_m']:.2f}, med coh={stats_p6[yr]['med_coh']:.2f}")

# ==============================================================================
# FIGURE 1: 4-PANEL (2x2) 1:1 KINEMATIC REGRESSIONS ACROSS 4 YEARS
# ==============================================================================
fig, axes = plt.subplots(2, 2, figsize=(15, 13), sharex=True, sharey=True, dpi=300)
lim = 24

panel_labels = ['A', 'B', 'C', 'D']
for idx, yr in enumerate(years):
    ax = axes[idx // 2, idx % 2]
    sub = subsets[yr]
    sp6 = stats_p6[yr]
    soff = stats_off[yr]
    dt_val = revisit_dict[yr]
    
    # Reference Lines
    ax.plot([-lim, lim], [-lim, lim], 'k--', lw=1.5, zorder=1, label='1 : 1 Identity Line')
    ax.axhspan(-14, 14, color='#e0ecf4', alpha=0.35, zorder=0, label='±λ/4 InSAR Range (±14 mm)')
    ax.axvspan(-3.5, 3.5, color='#f0f0f0', alpha=0.6, zorder=0, label='Laser Noise Floor (±3.5 mm)')
    ax.axvline(0, color='gray', lw=0.6, ls=':')
    ax.axhline(0, color='gray', lw=0.6, ls=':')
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.grid(True, ls=':', color='#cccccc', alpha=0.6)
    
    # Scatter points
    ax.scatter(sub.laser_los, sub.s1_dlos_mm, color='#1b7837', marker='^', s=80, alpha=0.85,
               edgecolors='#0a3314', zorder=4,
               label=f"Wood-Free Lawn: r = {soff['r']:.2f}, TLS m = {soff['tls_m']:.2f}")
    ax.scatter(sub.laser_los, sub.s1, color='#d95f02', marker='o', s=85, alpha=0.9,
               edgecolors='#4d1e00', zorder=5,
               label=f"Station P6 (8% Wood): r = {sp6['r']:.2f}, TLS m = {sp6['tls_m']:.2f}")
    
    # Fit lines
    x_line = np.linspace(sub.laser_los.min() - 2, sub.laser_los.max() + 2, 50)
    ax.plot(x_line, sp6['tls_m'] * x_line + sp6['tls_b'], color='#d95f02', lw=2.0, zorder=3)
    ax.plot(x_line, soff['tls_m'] * x_line + soff['tls_b'], color='#1b7837', lw=1.8, ls='--', zorder=3)
    
    title_era = "6-Day Constellation (S1A+S1B)" if yr == 2021 else "12-Day Revisit (S1A only)"
    ax.set_title(f"{panel_labels[idx]}. Season {yr} ({title_era})\nMedian Coherence γ = {sp6['med_coh']:.2f} (n = {sp6['n']} pairs)",
                 fontsize=11.5, fontweight='bold', pad=10)
    ax.legend(loc='upper left', fontsize=9, framealpha=0.92)
    
    # Metrics info box
    box_txt = (
        f"YEAR {yr} (MAY–NOV)\n"
        f"• Cadence: {dt_val} days\n"
        f"• Median Coherence: γ = {sp6['med_coh']:.2f}\n"
        f"• Sign Agreement: {sp6['sign_agree']:.0f}%\n"
        f"• P6 TLS Slope: m = {sp6['tls_m']:.2f}\n"
        f"• Lawn TLS Slope: m = {soff['tls_m']:.2f}"
    )
    box_color = '#f7fcf5' if yr == 2021 else '#fff5eb'
    border_color = '#1b7837' if yr == 2021 else '#d95f02'
    ax.text(0.55, 0.04, box_txt, transform=ax.transAxes, fontsize=8.2,
            verticalalignment='bottom', bbox=dict(boxstyle='round,pad=0.5', facecolor=box_color, edgecolor=border_color, alpha=0.95))

for ax in axes[1, :]:
    ax.set_xlabel('Ground Truth Laser Displacement in LOS (mm)', fontsize=11, fontweight='bold')
for ax in axes[:, 0]:
    ax.set_ylabel('Sentinel-1 InSAR Displacement (mm LOS)', fontsize=11, fontweight='bold')

fig.suptitle("Multi-Year Seasonal Kinematic Evaluation (May–November: 2021, 2022, 2023, 2024)\nGround-Truth Laser vs. InSAR: Distinguishing Cadence Decorrelation from Structural Boardwalk Damping",
             fontsize=13.5, fontweight='bold', y=0.99)
fig.tight_layout()

f1_name = "Fig_4Year_Seasonal_2021_2024_1to1.png"
fig.savefig(OUT_PAPER / f1_name, dpi=300)
fig.savefig(OUT_DEFENSE / f1_name, dpi=300)
fig.savefig(OUT_WEB / f1_name, dpi=300)
plt.close(fig)
print(f"Figure 1 saved to {f1_name}")

# ==============================================================================
# FIGURE 2: 4-PANEL STACKED TIME SERIES (MAY-NOV FOR EACH YEAR)
# ==============================================================================
fig, axes = plt.subplots(4, 1, figsize=(14, 16), sharey=True, dpi=300)

for idx, yr in enumerate(years):
    ax = axes[idx]
    sub = subsets[yr]
    dt_val = revisit_dict[yr]
    era_txt = "6-Day (S1A+S1B)" if yr == 2021 else "12-Day (S1A)"
    
    ax.axhline(0, color='gray', lw=0.7, ls=':')
    
    # Logger gap in 2022
    if yr == 2022:
        ax.axvspan(pd.to_datetime('2022-07-05'), pd.to_datetime('2022-08-08'), color='#fee8c8', alpha=0.6,
                   label='SDMS40 Logger Outage / QC Exclusion')
    
    # Motion tracks
    ax.plot(sub.mid_date, sub.laser_los, color='#2b83ba', lw=2.0, marker='s', ms=5.5,
            label='In-Situ Laser Ground Truth (LOS mm)', zorder=4)
    ax.plot(sub.mid_date, sub.s1, color='#d95f02', lw=2.0, marker='o', ms=5.5,
            label='Sentinel-1 InSAR (P6 Center Cell, 8% Wood)', zorder=5)
    ax.plot(sub.mid_date, sub.s1_dlos_mm, color='#1b7837', lw=1.8, ls='--', marker='^', ms=5.0,
            label='Sentinel-1 InSAR (Wood-Free Lawn, 0% Wood)', zorder=4)
    
    # Coherence secondary axis
    ax_coh = ax.twinx()
    bar_width = 3 if dt_val == 6 else 6
    ax_coh.bar(sub.mid_date, sub.coh, width=bar_width, color='#cccccc', alpha=0.35, label='Coherence (γ)', zorder=1)
    ax_coh.set_ylim(0, 1.05)
    ax_coh.set_ylabel('Coherence (γ)', color='#666666', fontsize=9.5)
    
    ax.set_xlim(pd.to_datetime(f"{yr}-05-01"), pd.to_datetime(f"{yr}-11-30"))
    ax.set_ylim(-18, 18)
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=14))
    ax.set_ylabel('Displacement Step (mm)', fontsize=10, fontweight='bold')
    
    ax.set_title(f"{panel_labels[idx]}. Season {yr} ({era_txt}) — Median Coherence γ = {stats_p6[yr]['med_coh']:.2f} (n = {len(sub)} pairs)",
                 fontsize=11, fontweight='bold', pad=8)
    ax.grid(True, ls=':', color='#cccccc', alpha=0.6)
    ax.legend(loc='lower left', fontsize=8.5, framealpha=0.92)

axes[3].set_xlabel('Observation Date', fontsize=11, fontweight='bold')
fig.suptitle("Four-Year Pair-by-Pair Kinematic Tracking Time Series (May–November: 2021, 2022, 2023, 2024)\nIn-Situ Laser vs. Sentinel-1 InSAR Under Identical Summer–Autumn Seasonality",
             fontsize=13, fontweight='bold', y=0.99)
fig.tight_layout()

f2_name = "Fig_4Year_Seasonal_2021_2024_TimeSeries.png"
fig.savefig(OUT_PAPER / f2_name, dpi=300)
fig.savefig(OUT_DEFENSE / f2_name, dpi=300)
fig.savefig(OUT_WEB / f2_name, dpi=300)
plt.close(fig)
print(f"Figure 2 saved to {f2_name}")
