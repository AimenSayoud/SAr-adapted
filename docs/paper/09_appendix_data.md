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
| D | 0.639 | 0.232 |


### T03. Usable fraction against the coherence threshold

| threshold | A | B | C | D |
|---|---|---|---|---|
| 0.4 | 1.0 | 1.0 | 1.0 | 0.987 |
| 0.45 | 0.996 | 1.0 | 1.0 | 0.987 |
| 0.5 | 0.974 | 0.938 | 0.995 | 0.983 |
| 0.55 | 0.84 | 0.692 | 0.976 | 0.931 |
| 0.6 | 0.531 | 0.431 | 0.882 | 0.734 |
| 0.65 | 0.236 | 0.185 | 0.794 | 0.429 |
| 0.7 | 0.054 | 0.015 | 0.647 | 0.232 |
| 0.75 | 0.006 | 0.0 | 0.446 | 0.125 |
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

| series | amplitude_mm |
|---|---|
| A−C | 3.286 |
| B−C | 2.627 |
| A−B | 0.901 |
| NULL | 0.568 |


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
| t2m_c | -0.509 | 0.224 |
| s2_wetness_C | 0.49 | 0.427 |
| s2_wetness_diff | 0.395 | -0.316 |
| api_mm | 0.23 | 0.293 |
| precip_mm | 0.19 | -0.225 |


### T10. Coherence response to freezing

| zone | slope_coh_per_wtd | r | n | coh_cold | coh_warm | freeze_gain | n_cold | n_warm |
|---|---|---|---|---|---|---|---|---|
| A | -10.69738714860094 | -0.2619408567691676 | 356 | 0.43503387128153154 | 0.40662583332795366 | 0.02840803795357788 | 31 | 325 |
| B | -10.052890114607944 | -0.23052825023619802 | 356 | 0.4178765329622453 | 0.3914231169223785 | 0.026453416039866784 | 31 | 325 |
| C | -9.49856925824952 | -0.20858960425208958 | 356 | 0.5609407780631896 | 0.4831868314743042 | 0.07775394658888546 | 31 | 325 |
| D | -8.040306375196218 | -0.1892504355196324 | 356 | 0.5082682428821441 | 0.4402202635544997 | 0.06804797932764439 | 31 | 325 |

