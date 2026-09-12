# Figures and tables — inventory and captions

> **Internal record. Not part of the manuscript.** Excluded from
> `SECTION_ORDER`; the captions that reach the article live in the section files
> next to each `![…]` reference.

Figures are produced by `notebooks/06_manuscript/export_figures_en.ipynb` into
`docs/paper/figures/` as **300 dpi PNG**, named `FXX_slug.png`.
Zone colours are constant throughout: **A** red `#d62728` · **B** blue
`#1f77b4` · **C** green `#2ca02c` · **D** grey `#9e9e9e`.

---

## Main figures

| # | File | Caption |
|---|---|---|
| **1** | `F01_study_area.png` | Study area. (a) Rzecin peatland outline with the Sentinel-1 burst footprint; (b) zone stratification A–D on the radar grid; (c) zone outlines over mean interferometric coherence. |
| **2** | `F02_hypotheses.png` | Study design: sequential falsification chain H1–H4, testing protocols applied, and verdicts reached. |
| **3** | `F03_temporal_coherence.png` | Phase-linking temporal coherence. (a) Temporal-coherence distributions by zone with the 0.55 noise floor and 0.7 reliability threshold; (b) multi-threshold usable fraction curve; (c) spatial map of temporal coherence from EVD phase linking with zone outlines. |
| **4** | `F04_paired_test.png` | Matched-cover paired comparison. (a) Distribution of coh(A) − coh(C) per interferogram; (b) scatter of A against C, points below the 1:1 line indicating lower mat coherence. |
| **5** | `F05_aggregate_and_significance.png` | Aggregated seasonal phase series and empirical significance against size-matched null distributions. (a) Time series for A−C, B−C, A−B and the size-matched null; (b) seasonal amplitudes across series; (c) null distribution and significance for the full network; (d) null distribution with winter pairs excluded. |
| **6** | `F06_wetness_anomaly_composite.png` | InSAR phase and optical surface wetness relationship. (a) Aggregated InSAR phase and Sentinel-2 optical wetness anomaly series; (b) cross-correlation with hydro-climatic forcings, comparing raw seasonal cycles with deseasonalised anomalies. |
| **7** | `F07_conceptual_framework.png` | Conceptual framework for wetland InSAR. (a) Four-level observable hierarchy distinguishing physical quantities, radar interactions, raw InSAR observables, and inferred quantities; (b) mechanistic causal chain at Rzecin, demonstrating how hydrologically driven phase-centre variability produces apparent displacement without requiring true mechanical peat breathing. |

## Supplementary figures

| # | File | Caption |
|---|---|---|
| **S1** | `S01_network.png` | Interferometric network. (a) 356 pairs over ~90 acquisitions, 2022–2024; (b) temporal-baseline distribution, with the 60-day robustness filter marked. |
| **S2** | `S02_rgb_composite.png` | False-colour composite (R = σ⁰ VV, G = coherence, B = Sentinel-2 wetness) with zone outlines, and a zoom on the peatland. |
| **S3** | `S03_flooded_fraction.png` | Inundated-time fraction from the water mask, with zones A and B outlined. |
| **S4** | `S04_protocol.png` | Processing protocol: per-pixel versus aggregated observable, and the size-matched null test yielding an empirical p-value. |
| **S5** | `S05_synthetic_validation.png` | Synthetic validation. On identical simulated data, per-pixel inversion returns −13.7 mm yr⁻¹ (36 % usable pixels) whereas aggregation returns −19.8 mm yr⁻¹ against a ground truth of −20. |
| **S6** | `S06_coherence_decay.png` | Coherence decay with temporal baseline, with fitted decorrelation times per zone. |
| **S7** | `S07_zone_distributions.png` | Per-zone distributions across five independent sensors: coherence, σ⁰ VV, RVI, Sentinel-2 wetness and temporal coherence. |
| **S8** | `S08_radial_profiles.png` | Radial profiles of coherence, σ⁰ VV and RVI against signed distance to the peatland boundary. The step at distance zero marks a physical edge. |
| **S9** | `S09_predictors.png` | Within-mat predictive model. (a) Standardised coefficients of the collinearity-cleaned model; (b) Spearman correlations in mat versus grassland, showing active environmental sensitivity in the mat against its absence in grassland. |
| **S10** | `S10_amplitude_dispersion.png` | Amplitude dispersion D_A: map and per-zone distributions against the 0.25 persistent-scatterer threshold. |
| **S11** | `S11_hydrology_freeze.png` | Coherence sensitivity to the water-table proxy, and coherence gain on freezing, by zone. |
| **S12** | `S12_closure_phase.png` | Closure phase by zone. (a) Mean bias with 2σ error bars — none significant; (b) median |closure| dispersion against the π/2 random reference. |
| **S13** | `S13_aggregation_gain.png` | Aggregation gain curve. Empirical phase standard deviation as a function of aggregated pixel count $N$ across Zone A, overlaid with theoretical independent $1/\sqrt{N}$ and autocorrelated $1/\sqrt{N_{\text{eff}}}$ scaling ($N_{\text{eff}} \approx 31$). |
| **S14** | `S14_literature_context.png` | Our bound in context: raised-bog breathing, drained-fen subsidence, expected free flotation, and the value measured here. |

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
