# Figures and tables — inventory and captions

Figures are produced by `notebooks/export_figures_en.ipynb` into
`docs/paper/figures/` as **300 dpi PNG**, named `FXX_slug.png`.
Zone colours are constant throughout: **A** red `#d62728` · **B** blue
`#1f77b4` · **C** green `#2ca02c` · **D** grey `#9e9e9e`.

---

## Main figures

| # | File | Caption (draft) |
|---|---|---|
| **1** | `F01_hypotheses.png` | Study design: four competing hypotheses, the test applied to each, and the verdict reached. |
| **2** | `F02_study_area.png` | Study area. (a) Rzecin peatland outline with the Sentinel-1 burst footprint; (b) zone stratification A–D on the radar grid; (c) zone outlines over mean interferometric coherence. |
| **3** | `F03_network.png` | Interferometric network. (a) 356 pairs over ~90 acquisitions, 2022–2024; (b) temporal-baseline distribution with the 60-day robustness filter marked. |
| **4** | `F04_protocol.png` | Processing protocol: per-pixel versus aggregated observable, and the size-matched null test yielding an empirical *p*-value. |
| **5** | `F05_temporal_coherence.png` | Phase-linking result. (a) Temporal-coherence distributions by zone with the 0.55 noise floor and the 0.7 reliability threshold; (b) full multi-threshold curve. |
| **6** | `F06_tcoh_map.png` | Map of temporal coherence from EVD phase linking, with zone outlines. |
| **7** | `F07_zone_distributions.png` | Per-zone distributions across five independent sensors: coherence, σ⁰ VV, RVI, Sentinel-2 wetness, temporal coherence. |
| **8** | `F08_radial_profiles.png` | Radial profiles of coherence, σ⁰ VV and RVI against signed distance to the peatland boundary; the step at distance zero marks a physical edge. |
| **9** | `F09_paired_test.png` | Matched-cover paired test. (a) Distribution of coh(A) − coh(C) per interferogram; (b) scatter of A against C, with points below the 1:1 line indicating lower mat coherence. |
| **10** | `F10_predictors.png` | Within-mat predictive model. (a) Standardised coefficients of the collinearity-cleaned model; (b) Spearman correlations in mat versus grassland, showing sign reversal. |
| **11** | `F11_aggregate_series.png` | Aggregated series. (a) A−C, B−C, A−B and the size-matched null; (b) seasonal amplitudes — the lake oscillates like the mat, and A−B cancels. |
| **12** | `F12_significance.png` | Significance of the seasonal amplitude against size-matched null distributions. (a) Full network; (b) winter pairs excluded. |
| **13** | `F13_closure_phase.png` | Closure phase by zone. (a) Mean bias with 2σ error bars — none significant; (b) median \|closure\| dispersion against the π/2 random reference. |
| **14** | `F14_phase_vs_wetness.png` | Aggregated InSAR phase and Sentinel-2 optical wetness — two independent sensors. |
| **15** | `F15_seasonal_vs_anomalies.png` | Correlation with each forcing, raw versus deseasonalised: temperature collapses while wetness survives. |
| **16** | `F16_causal_chain.png` | Proposed mechanism as a causal chain, from water table to the measured dielectric signal. |
| **17** | `F17_literature_context.png` | Our bound in context: raised-bog breathing, drained-fen subsidence, expected free flotation, and the value measured here. |

## Supplementary figures

| # | File | Caption |
|---|---|---|
| **S1** | `S01_rgb_composite.png` | False-colour composite (R = σ⁰ VV, G = coherence, B = S2 wetness) with zone outlines and a zoom on the peatland. |
| **S2** | `S02_flooded_fraction.png` | Inundated-time fraction from the water mask, with zones A and B outlined. |
| **S3** | `S03_synthetic_validation.png` | Synthetic validation: per-pixel inversion versus aggregation on identical data against a known ground truth. |
| **S4** | `S04_amplitude_dispersion.png` | Amplitude dispersion D_A: map and per-zone distributions against the 0.25 persistent-scatterer threshold. |
| **S5** | `S05_coherence_decay.png` | Coherence decay with temporal baseline, with fitted decorrelation times per zone. |
| **S6** | `S06_hydrology_freeze.png` | Coherence sensitivity to water-table proxy, and coherence gain on freezing, by zone. |

## Tables

| # | File | Content |
|---|---|---|
| **1** | `T01_zones.csv` | Zone pixel counts, areas, and departure from the documented 89.7 ha |
| **2** | `T02_temporal_coherence.csv` | Temporal coherence statistics and fraction ≥ 0.7 by zone |
| **3** | `T03_multithreshold.csv` | Full multi-threshold curve |
| **4** | `T04_zone_signature.csv` | Multi-sensor median signature by zone |
| **5** | `T05_paired_test.csv` | Paired A−C test: Wilcoxon, date-jackknife, confidence interval |
| **6** | `T06_predictors.csv` | Within-mat predictors and the A-versus-C sign reversal |
| **7** | `T07_seasonal_amplitudes.csv` | Seasonal amplitudes of the four aggregated series |
| **8** | `T08_closure_phase.csv` | Closure-phase bias and dispersion by zone |
| **9** | `T09_forcings.csv` | Forcing correlations, raw versus deseasonalised |
| **10** | `T10_hydrology_freeze.csv` | Hydrological coupling and freeze test by zone |

---

## Production notes

- **One figure requires an external basemap** (the regional locator inset of
  Figure 2a) and is drawn schematically from the AOI polygon coordinates; replace
  with a proper basemap before submission if the journal requires one.
- Expensive steps — per-pair coherence (~5 min), null distributions (~15 min),
  closure triplets (~3 min) — are **cached** in `figures_cache/`, so re-runs are
  immediate. Cells are marked `[SLOW]`.
- Font sizes, dpi and zone colours are set once through `mpl.rcParams` and
  `zone_viz.ZONE_COLORS`, so all figures are visually consistent.
