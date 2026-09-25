## Appendix B. Numeric tables

Tables exported by `notebooks/06_manuscript/export_figures_en.ipynb`. Each is reproduced here and shipped as CSV alongside the figures.

### T01. Zone definition, pixel counts and areas

| zone | label | n_px | area_ha |
|---|---|---|---|
| A | A — floating mat | 499 | 79.84 |
| B | B — residual lake | 65 | 10.4 |
| C | C — matched grassland | 398 | 63.68 |
| D | D — other cover | 10750 | 1720.0 |


### T02. Temporal coherence by zone (phase linking)

| zone | median | frac_ge_0.7 |
|---|---|---|
| A | 0.604 | 0.054 |
| B | 0.584 | 0.015 |
| C | 0.734 | 0.647 |
| D | 0.639 | 0.231 |


### T03. Usable fraction against the coherence threshold

| threshold | A | B | C | D |
|---|---|---|---|---|
| 0.4 | 1.0 | 1.0 | 1.0 | 0.987 |
| 0.45 | 0.996 | 1.0 | 1.0 | 0.987 |
| 0.5 | 0.974 | 0.938 | 0.995 | 0.983 |
| 0.55 | 0.84 | 0.692 | 0.976 | 0.931 |
| 0.6 | 0.531 | 0.431 | 0.882 | 0.734 |
| 0.65 | 0.236 | 0.185 | 0.794 | 0.429 |
| 0.7 | 0.054 | 0.015 | 0.647 | 0.231 |
| 0.75 | 0.006 | 0.0 | 0.447 | 0.125 |
| 0.8 | 0.0 | 0.0 | 0.267 | 0.062 |
| 0.85 | 0.0 | 0.0 | 0.104 | 0.024 |
| 0.9 | 0.0 | 0.0 | 0.029 | 0.006 |
| 0.95 | 0.0 | 0.0 | 0.0 | 0.0 |


### T04. Multi-sensor signature of the zones

| zone | $\sigma^0$ VV (dB) | RVI (volume) | S2 wetness | coherence | temporal coherence |
|---|---|---|---|---|---|
| A | -10.094 | 0.914 | -0.513 | 0.408 | 0.604 |
| B | -15.412 | 0.993 | 0.185 | 0.396 | 0.584 |
| C | -11.224 | 0.881 | -0.522 | 0.492 | 0.734 |
| D | -9.276 | 1.045 | -0.44 | 0.438 | 0.639 |


### T05. Paired mat-versus-grassland coherence test

| n_pairs | n_dates | delta_mean | delta_median | frac_a_lower | wilcoxon_stat | wilcoxon_p | date_jackknife_min | date_jackknife_max | date_jackknife_se | robust_same_sign | significant |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 356 | 90 | -0.08085797333650374 | -0.05030900239944461 | 0.8932584269662921 | 4094.0 | 4.838844154316786e-46 | -0.08418832078229549 | -0.07741065939952588 | 0.014535361854110865 | True | True |


### T06. Failure predictors, with and without collinear terms

| covariate | A_mat | C_grassland |
|---|---|---|
| sigma0_vv | -0.379 | -0.008 |
| s2_greenness_mean | 0.32 | -0.009 |
| s2_wetness_mean | -0.223 | 0.026 |
| s2_greenness_amp | -0.196 | -0.406 |
| rvi | -0.185 | -0.252 |
| elevation | -0.168 | 0.43 |
| dist_edge_m | 0.153 | 0.188 |
| sigma0_std | 0.122 | 0.132 |


### T07. Seasonal amplitudes and empirical p-values

| series | amplitude_mm | p_perm | n_null | null_type |
|---|---|---|---|---|
| A−C | 3.286 | 0.026 | 4614.0 | reference-matched |
| B−C | 2.627 | 0.136 | 249.0 | reference-matched |
| A−B | 0.901 |  |  | - |
| NULL | 0.568 |  |  | - |


### T08. Closure-phase bias and dispersion by zone

| zone | n_triplets | mean_closure_rad | median_abs_rad | se | bias_significant |
|---|---|---|---|---|---|
| A | 518 | -0.0903481766581535 | 0.6832633018493652 | 0.05680459476682 | False |
| B | 518 | -0.088422305881977 | 0.7767355442047119 | 0.0584010888257125 | False |
| C | 518 | 0.0268979929387569 | 0.2123883962631225 | 0.0248665112638881 | False |
| D | 518 | -0.0213026106357574 | 0.2104390859603881 | 0.019004098447159 | False |


### T09. Hydro-climatic forcings, raw and deseasonalised

| driver | raw (annual cycle in) | ANOMALIES |
|---|---|---|
| s2_wetness | 0.576 | 0.45 |
| s2_wetness_D | 0.519 | 0.424 |
| t2m_c | -0.506 | 0.215 |
| s2_wetness_C | 0.49 | 0.427 |
| s2_wetness_diff | 0.395 | -0.316 |
| api_mm | 0.27 | 0.379 |
| precip_mm | 0.204 | 0.224 |


### T10. Coherence response to freezing

| zone | slope_coh_per_wtd | r | n | coh_cold | coh_warm | freeze_gain | n_cold | n_warm |
|---|---|---|---|---|---|---|---|---|
| A | -1.5548147624380184 | -0.22083462863240225 | 356 | 0.43503387128153154 | 0.40662583332795366 | 0.02840803795357788 | 31 | 325 |
| B | -1.4437907623558621 | -0.19204387788218316 | 356 | 0.4178765329622453 | 0.3914231169223785 | 0.026453416039866784 | 31 | 325 |
| C | -1.2016069127337878 | -0.15305935338051163 | 356 | 0.5609407780631896 | 0.4831868314743042 | 0.07775394658888546 | 31 | 325 |
| D | -0.9939134380631856 | -0.1356987483801386 | 356 | 0.5082682428821441 | 0.4402202635544997 | 0.06804797932764439 | 31 | 325 |


### T11. T11 subzone subdivision

| subzone | n_px | amplitude_mm | phase_doy | r2_seasonal |
|---|---|---|---|---|
| Full mat (Zone A) | 499 | 3.286 | 104.2 | 0.299 |
| Inner core (d > 40 m) | 356 | 3.528 | 106.6 | 0.312 |
| Deep core (d > 80 m) | 233 | 3.781 | 108.3 | 0.318 |
| Outer margin (d <= 40 m) | 143 | 2.894 | 99.8 | 0.245 |


### T12. T12 lake erosion

| depth_rings | distance_threshold_m | n_px | amplitude_mm | phase_doy | r2 | status |
|---|---|---|---|---|---|---|
| 0 | 0 | 65 | 2.627 | 94.7 | 0.114 | Full lake |
| 1 | 40 | 26 | 2.217 | 93.6 | 0.093 | Interior lake |
| 2 | 80 | 4 | 1.842 | 92.1 | 0.065 | Deep center |
| 3 | 120 | 0 |  |  |  | Extinct (geometric limit) |


### T13. T13 aggregation gain

| n_pixels | theoretical_independent_sd_mm | theoretical_autocorrelated_sd_mm | empirical_sd_mm |
|---|---|---|---|
| 1 | 6.6 | 6.6 | 6.6 |
| 5 | 2.952 | 3.135 | 3.291 |
| 10 | 2.087 | 2.369 | 2.457 |
| 25 | 1.32 | 1.755 | 1.748 |
| 50 | 0.933 | 1.496 | 1.444 |
| 100 | 0.66 | 1.347 | 1.28 |
| 250 | 0.417 | 1.25 | 1.207 |
| 499 | 0.295 | 1.216 | 1.212 |


### T14. T14 power ab

| null_p95_mm | power | alpha | se_null_mm | min_detectable_amp_mm | observed_amp_mm | observed_p_value |
|---|---|---|---|---|---|---|
| 2.0 | 0.8 | 0.05 | 1.216 | 3.023 | 0.9 | 0.448 |


### T15. T15 subset stability

| subset | max_dt_days | n_pairs | velocity_mm_yr | amplitude_mm | phase_doy | stability |
|---|---|---|---|---|---|---|
| <=24d | 24 | 175 | -13.471 | 9.496 | 116.8 | Closure-phase contaminated |
| <=36d | 36 | 261 | -8.495 | 4.937 | 103.6 | Transitioning |
| <=48d | 48 | 346 | -3.867 | 2.89 | 105.9 | Stabilized |
| All pairs | 9999 | 356 | -1.531 | 3.286 | 104.2 | Full network constrained |


### T16. T16 saturating seasonal fit

| model | semi_amplitude_mm | peak_to_peak_mm | ceiling_mm | exceeds_ceiling | phase_doy | r2 |
|---|---|---|---|---|---|---|
| Linear harmonic | 3.286 | 6.573 | 6.13 | True | 104.2 | 0.3642 |
| Free saturating (tanh) | 2.822 | 5.643 | 6.13 | False | 100.2 | 0.3719 |
| Ceiling-constrained (tanh, 6.13 mm) | 2.912 | 5.825 | 6.13 | False | 100.9 | 0.3718 |

