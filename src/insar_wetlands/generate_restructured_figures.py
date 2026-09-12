"""Generate restructured figure set for the manuscript.

Restructures the 16 main figures into:
- 7 main-text figures (F01 to F07)
- 14 supplementary figures (S01 to S14)
Per referee consolidated report §9 and tickets P-037, P-038, P-039.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as patches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402


def create_observable_hierarchy_and_mechanism_figure(out_path: Path) -> None:
    """Generate Figure 7: (a) Observable Hierarchy, (b) Physical mechanism causal chain."""
    fig = plt.figure(figsize=(13.5, 9.2), dpi=300)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15], wspace=0.18)

    # -------------------------------------------------------------
    # Panel (a): InSAR Observable Hierarchy
    # -------------------------------------------------------------
    ax_a = fig.add_subplot(gs[0])
    ax_a.axis("off")
    ax_a.set_title("(a) Four-Level Observable Hierarchy for Wetland InSAR",
                    fontsize=12, fontweight="bold", pad=14, loc="left")

    hierarchy_levels = [
        {
            "level": "LEVEL 1: PHYSICAL QUANTITY",
            "desc": "Ground-truth bio-physical state & kinematics",
            "items": "• Water-table depth (WTD)\n• Peat surface deformation (elastic breathing, compaction)\n• Soil & canopy moisture content\n• Emergent vegetation structure & biomass",
            "color": "#d0e1fd",
            "border": "#2b5c8f",
        },
        {
            "level": "LEVEL 2: RADAR INTERACTION",
            "desc": "Electromagnetic wave propagation & scattering",
            "items": "• Penetration depth (δ_p ≈ 3.6 mm at C-band)\n• Volumetric scattering vs surface reflection\n• Phase-centre elevation & vertical displacement\n• Temporal & volumetric decorrelation (γ_temp, γ_vol)",
            "color": "#e8d8f8",
            "border": "#6a3d9a",
        },
        {
            "level": "LEVEL 3: InSAR OBSERVABLE",
            "desc": "Direct interferometric phase & coherence measurements",
            "items": "• Multi-looked wrapped phase (ϕ_ij)\n• Interferometric coherence (γ_ij)\n• Non-zero closure phase (ψ_ijk)\n• Spatial aggregate differential phase (ΔΦ_AC, ΔΦ_AB)",
            "color": "#fff2cc",
            "border": "#b27b00",
        },
        {
            "level": "LEVEL 4: INFERRED QUANTITY",
            "desc": "Derived geophysical interpretation (requires model assumptions)",
            "items": "• Apparent line-of-sight (LOS) displacement (d_LOS)\n• Vertical apparent displacement bound (d_z ≤ 8.7 mm)\n• Surface wetness / dielectric anomaly association (r ≈ 0.45)\n• Mechanical peat motion (unconstrained without laser ground truth)",
            "color": "#d5e8d4",
            "border": "#274e13",
        },
    ]

    y_top = 0.94
    box_h = 0.175
    box_gap = 0.055

    for idx, lvl in enumerate(hierarchy_levels):
        y = y_top - idx * (box_h + box_gap)
        rect = patches.FancyBboxPatch(
            (0.03, y - box_h), 0.94, box_h,
            boxstyle="round,pad=0.015,rounding_size=0.02",
            fc=lvl["color"], ec=lvl["border"], lw=1.5, alpha=0.95
        )
        ax_a.add_patch(rect)

        ax_a.text(0.06, y - 0.032, lvl["level"],
                  fontsize=9.5, fontweight="bold", color=lvl["border"], va="top")
        ax_a.text(0.06, y - 0.060, lvl["desc"],
                  fontsize=8.0, style="italic", color="#333333", va="top")

        ax_a.text(0.08, y - 0.088, lvl["items"],
                  fontsize=8.0, color="#111111", va="top", linespacing=1.35)

        if idx < len(hierarchy_levels) - 1:
            arrow_y_start = y - box_h
            arrow_y_end = y - box_h - box_gap
            ax_a.annotate(
                "", xy=(0.5, arrow_y_end + 0.008), xytext=(0.5, arrow_y_start - 0.005),
                arrowprops=dict(arrowstyle="-|>", lw=2.0, color="#444444", mutation_scale=15)
            )
            ax_a.text(0.53, (arrow_y_start + arrow_y_end) / 2, "maps via forward model",
                      fontsize=7.2, color="#555555", va="center", style="italic")

    ax_a.set_xlim(0, 1)
    ax_a.set_ylim(0.0, 1.0)

    # -------------------------------------------------------------
    # Panel (b): Mechanistic Causal Chain at Rzecin
    # -------------------------------------------------------------
    ax_b = fig.add_subplot(gs[1])
    ax_b.axis("off")
    ax_b.set_title("(b) Mechanistic Causal Chain at Rzecin Floating Fen",
                    fontsize=12, fontweight="bold", pad=14, loc="left")

    causal_steps = [
        ("Near-surface water table & floating root mat", "#cfe0f3", "Level 1: Physical"),
        ("Saturated Sphagnum canopy & variable moisture (Δm_v ≈ 0.25)", "#cfe0f3", "Level 1: Physical"),
        ("High dielectric constant (ε_w = 73 − 25j) & shallow skin depth (δ_p ≈ 3.6 mm)", "#d9c7ef", "Level 2: Interaction"),
        ("Variable penetration depth + dominant volume scattering (RVI ↑, σ⁰ ↑)", "#d9c7ef", "Level 2: Interaction"),
        ("Scattering phase-centre elevation unstable between satellite passes", "#ffe6b3", "Level 2: Interaction"),
        ("Severe temporal decorrelation (coherence collapse: mean Δ = −0.081)", "#ffe6b3", "Level 3: Observable"),
        ("Per-pixel phase unwrapping & inversion fail across all 6 tested methods", "#f4c7c3", "Level 3: Observable"),
        ("Spatial aggregation over mat: coherent averaging divides noise by √N_eff", "#c9e2c9", "Level 3: Observable"),
        ("Detectable seasonal phase oscillation: amplitude = 3.29 mm (p = 0.014)", "#c9e2c9", "Level 3: Observable"),
        ("Correlates with independent optical wetness anomalies at near-zero lag (r ≈ 0.45)", "#c9e2c9", "Level 4: Inferred"),
        ("Lake oscillates with matched phase/amplitude; mat − lake cancels;\nsignal consistent with dielectric perturbation; mechanical motion unconstrained", "#b7d7b7", "Level 4: Inferred"),
    ]

    n_b = len(causal_steps)
    h_b = 0.068
    gap_b = 0.015
    y_top_b = 0.94

    for i, (text, col, tag) in enumerate(causal_steps):
        y = y_top_b - i * (h_b + gap_b)
        rect = patches.FancyBboxPatch(
            (0.04, y - h_b), 0.92, h_b,
            boxstyle="round,pad=0.01,rounding_size=0.015",
            fc=col, ec="#333333", lw=1.1, alpha=0.95
        )
        ax_b.add_patch(rect)

        ax_b.text(0.94, y - 0.015, tag, fontsize=6.5, color="#555555",
                  ha="right", va="top", style="italic")

        ax_b.text(0.07, y - h_b / 2 - 0.003, text, ha="left", va="center",
                  fontsize=7.8, linespacing=1.25, color="#111111")

        if i < n_b - 1:
            ax_b.annotate(
                "", xy=(0.5, y - h_b - gap_b + 0.003), xytext=(0.5, y - h_b),
                arrowprops=dict(arrowstyle="-|>", lw=1.4, color="#333333", mutation_scale=11)
            )

    ax_b.set_xlim(0, 1)
    ax_b.set_ylim(0.0, 1.0)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Created {out_path}")


def composite_temporal_coherence(f05_path: Path, f06_path: Path, out_path: Path) -> None:
    """Composite Figure 3: Left F05 (distributions + multithreshold), Right F06 (map)."""
    im_f05 = Image.open(f05_path).convert("RGBA")
    im_f06 = Image.open(f06_path).convert("RGBA")

    target_h = max(im_f05.height, im_f06.height)
    w05 = int(im_f05.width * (target_h / im_f05.height))
    im_f05_scaled = im_f05.resize((w05, target_h), Image.Resampling.LANCZOS)

    w06 = int(im_f06.width * (target_h / im_f06.height))
    im_f06_scaled = im_f06.resize((w06, target_h), Image.Resampling.LANCZOS)

    gap = 40
    total_w = w05 + gap + w06

    composite = Image.new("RGBA", (total_w, target_h), (255, 255, 255, 255))
    composite.paste(im_f05_scaled, (0, 0), im_f05_scaled)
    composite.paste(im_f06_scaled, (w05 + gap, 0), im_f06_scaled)

    draw = ImageDraw.Draw(composite)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 44)
    except Exception:
        font = ImageFont.load_default()

    draw.text((w05 + gap + 30, 25), "(c) Spatial map (EVD)", fill=(0, 0, 0, 255), font=font)

    composite.convert("RGB").save(out_path, "PNG", dpi=(300, 300))
    print(f"Created {out_path}")


def composite_stacked(top_path: Path, bot_path: Path, out_path: Path) -> None:
    """Composite two figures vertically."""
    im_top = Image.open(top_path).convert("RGBA")
    im_bot = Image.open(bot_path).convert("RGBA")

    target_w = max(im_top.width, im_bot.width)

    h_top = int(im_top.height * (target_w / im_top.width))
    im_top_scaled = im_top.resize((target_w, h_top), Image.Resampling.LANCZOS)

    h_bot = int(im_bot.height * (target_w / im_bot.width))
    im_bot_scaled = im_bot.resize((target_w, h_bot), Image.Resampling.LANCZOS)

    gap = 30
    total_h = h_top + gap + h_bot

    composite = Image.new("RGBA", (target_w, total_h), (255, 255, 255, 255))
    composite.paste(im_top_scaled, (0, 0), im_top_scaled)
    composite.paste(im_bot_scaled, (0, h_top + gap), im_bot_scaled)

    composite.convert("RGB").save(out_path, "PNG", dpi=(300, 300))
    print(f"Created {out_path}")


def build_all_restructured_figures(fig_dir: Path) -> dict:
    """Generate all 7 main figures and 14 supplementary figures."""
    report = {}

    # Main Figure 1: Study area (copy of F02_study_area.png)
    shutil.copyfile(fig_dir / "F02_study_area.png", fig_dir / "F01_study_area.png")
    report["F01_study_area.png"] = "copied from F02_study_area.png"

    # Main Figure 2: Hypotheses / sequential falsification (copy of F01_hypotheses.png)
    shutil.copyfile(fig_dir / "F01_hypotheses.png", fig_dir / "F02_hypotheses.png")
    report["F02_hypotheses.png"] = "copied from F01_hypotheses.png"

    # Main Figure 3: Phase linking temporal coherence composite (F05 + F06)
    composite_temporal_coherence(
        fig_dir / "F05_temporal_coherence.png",
        fig_dir / "F06_tcoh_map.png",
        fig_dir / "F03_temporal_coherence.png"
    )
    report["F03_temporal_coherence.png"] = "composite of F05 and F06"

    # Main Figure 4: Paired test (copy of F07_paired_test.png)
    shutil.copyfile(fig_dir / "F07_paired_test.png", fig_dir / "F04_paired_test.png")
    report["F04_paired_test.png"] = "copied from F07_paired_test.png"

    # Main Figure 5: Aggregated seasonal series + null significance (F11 + F12)
    composite_stacked(
        fig_dir / "F11_aggregate_series.png",
        fig_dir / "F12_significance.png",
        fig_dir / "F05_aggregate_and_significance.png"
    )
    report["F05_aggregate_and_significance.png"] = "composite of F11 and F12"

    # Main Figure 6: InSAR phase vs optical wetness anomaly (F14 + F15)
    composite_stacked(
        fig_dir / "F14_phase_vs_wetness.png",
        fig_dir / "F15_seasonal_vs_anomalies.png",
        fig_dir / "F06_wetness_anomaly_composite.png"
    )
    report["F06_wetness_anomaly_composite.png"] = "composite of F14 and F15"

    # Main Figure 7: Conceptual framework (Observable Hierarchy + Causal Chain)
    create_observable_hierarchy_and_mechanism_figure(
        fig_dir / "F07_conceptual_framework.png"
    )
    report["F07_conceptual_framework.png"] = "Observable Hierarchy + Causal Chain composite"

    # -------------------------------------------------------------
    # Supplementary Figures S01 to S14
    # -------------------------------------------------------------
    supp_mapping = [
        ("S01_network.png", "F03_network.png"),
        ("S02_rgb_composite.png", "S01_rgb_composite.png"),
        ("S03_flooded_fraction.png", "S02_flooded_fraction.png"),
        ("S04_protocol.png", "F04_protocol.png"),
        ("S05_synthetic_validation.png", "S03_synthetic_validation.png"),
        ("S06_coherence_decay.png", "S04_coherence_decay.png"),
        ("S07_zone_distributions.png", "F08_zone_distributions.png"),
        ("S08_radial_profiles.png", "F09_radial_profiles.png"),
        ("S09_predictors.png", "F10_predictors.png"),
        ("S10_amplitude_dispersion.png", "S05_amplitude_dispersion.png"),
        ("S11_hydrology_freeze.png", "S06_hydrology_freeze.png"),
        ("S12_closure_phase.png", "F13_closure_phase.png"),
        ("S13_aggregation_gain.png", "F17_aggregation_gain.png"),
        ("S14_literature_context.png", "F16_literature_context.png"),
    ]

    for new_name, src_name in supp_mapping:
        src = fig_dir / src_name
        dst = fig_dir / new_name
        shutil.copyfile(src, dst)
        report[new_name] = f"copied from {src_name}"

    return report


if __name__ == "__main__":
    hub_figs = Path("/Users/aymen/Documents/Research_Hub/05_code/SAr-adapted/docs/paper/figures")
    rep = build_all_restructured_figures(hub_figs)
    print(f"Generated {len(rep)} figure files.")
