# Traceability — number → source

> **Internal record. Not part of the manuscript.** This file is deliberately
> excluded from `SECTION_ORDER`, so it is never assembled into the article or
> the `.docx`. It exists to map every published number back to the code that
> produced it — which is why processing-stage names appear here, and nowhere in
> the manuscript itself.

Every number in the manuscript must be findable here. No figure appears in the
text without a corresponding row.

## Processing stages → manuscript sections

| Source | Contribution | Section |
|---|---|---|
| `phase01`–`phase07` | acquisition, network, masks, S2 fusion | §2 |
| `phase08_sbas_mintpy` | SBAS — failure | §4.1 |
| `phase09_isbas` | ISBAS — failure | §4.1 |
| `phase15_annual_pairs` | annual pairs — failure | §4.1 |
| `phaseA_hybrid_insar_test` | hybrid network — 0 reliable pixel over AOI | §4.1 |
| `phaseE2_evd_phaselinking` | **phase linking — decisive test** | §4.1 |
| `phaseD_inside_vs_outside` | paired Δcoherence A vs C | §4.2 |
| `phaseDbis_mechanism_spatial` | hydrology, freeze, radial profile, residuals | §4.2 |
| `phaseDter_scattering_scatterers` | σ⁰, VH/VV, D_A, lake, phenology | §4.2 |
| `phaseH_predict_failure` | within-mat model, driver reversal | §4.2 |
| `phaseG_aggregation` | **aggregation, seasonal, lake, closure** | §4.3 |
| `phaseI_hydro_sensor` | zone validation + hydrological link | §2, §4.4 |
| `phaseJ_falsification` | N_eff, incidence, snow/frost, *p* floor | §4.5, App. A |
| `phaseF_nisar_lband` | L-band chain (prospective) | §5.4 |
| `phaseB_egms_check`, `phaseC1_licsbas` | external controls | §5 |

## Site and zones (§2)

| Quantity | Value | Source |
|---|---|---|
| Documented area | 89.7 ha | site literature |
| A / B / C / D | 499 / 65 / 398 / 10 750 px | `phaseI` |
| Areas | 79.84 / 10.40 / 63.68 / 1 720 ha | `phaseI` |
| **A+B vs documented** | **90.24 ha, +0.6 %** | `phaseI` |
| Interferograms | 356 pairs, ~90 dates | `phase03` |
| **Incidence (measured)** | **32.26°** → LOS→vertical **1.183** | `phaseJ` J2 |

## H1 — algorithm (§4.1)

| Quantity | Value | Source |
|---|---|---|
| Temporal coherence, median A/B/C/D | 0.604 / 0.584 / 0.734 / 0.639 | `phaseE2` |
| % ≥ 0.7 | 5.4 / 1.5 / 64.7 / 23.1 | `phaseE2` |
| Noise floor (redundancy ≈ 4) | ≈ 0.55 | `tests/test_synthetic_phaselinking.py` |
| WLS residual A / C | 2.46 / 1.92 rad | `phaseDbis` |
| Hybrid network | 14 238 "resolved", 0 reliable over AOI | `phaseA` |

## H2 — target (§4.2)

| Quantity | Value | Source |
|---|---|---|
| **Paired Δ coh(A) − coh(C)** | **−0.069** | `phaseD` |
| Wilcoxon | *p* = 2.2 × 10⁻⁴⁹ | `phaseD` |
| Date-jackknife | [−0.0705, −0.0652], sign stable | `phaseD` |
| Coherence A/B/C/D | 0.408 / 0.396 / 0.492 / 0.438 | `phaseI` |
| σ⁰ VV A/B/C/D | −10.09 / −15.41 / −11.22 / −9.28 dB | `phaseI` |
| RVI A/B/C/D | 0.914 / 0.993 / 0.881 / 1.045 | `phaseI` |
| S2 wetness A/B/C/D | −0.513 / +0.185 / −0.522 / −0.440 | `phaseI` |
| Decorrelation time τ | A 21 d / C 32 d | `phaseD` |
| Closure dispersion | A 0.683 / C 0.212 (**×3.2**) | `phaseG` |
| D_A (% < 0.25) | A 59 %, C 69 %, B 3 % | `phaseDter` |
| Cross-validated R² within mat | **0.239 ± 0.022** (RF 0.326) | `phaseH` |
| Cleaned coefficients | σ⁰ −0.275, RVI −0.251, wetness −0.240 | `phaseH` |
| VIF before cleaning | vh_vv_db 240.8, rvi 235.7 | `phaseH` |
| Driver reversal | σ⁰ −0.379 vs −0.008; elevation −0.168 vs +0.430 | `phaseH` |

## H3 — mechanism (§4.3)

| Quantity | Value | Source |
|---|---|---|
| Synthetic validation | per-pixel −13.7 vs aggregate −19.8 (truth −20) | `tests/test_aggregate.py` |
| \|R\| A/B/C/D | 0.234 / 0.395 / 0.569 / 0.426 | `phaseG` |
| L_corr / N_eff measured | A 160 m / 31 · C 360 m / 5 · D 280 m / 219 | `phaseJ` J1 |
| Velocity A−C vs null | −1.53 vs −1.50 mm yr⁻¹ → not significant | `phaseG` |
| Velocity detection floor | ≈ 1.5–5 mm yr⁻¹ | `phaseG` |
| **Amplitude A−C** | **3.29 mm**, DOY 104, ***p* = 0.014** (280 nulls, floor 0.0036) | `phaseJ` J4 |
| Amplitude A−C **winter excluded** | **3.282 mm** (−0.1 %), R² 0.309, *p* = 0.022 | `phaseJ` J3 |
| Matched null | median 0.86, p95 1.87 mm | `phaseJ` |
| **Amplitude B−C (lake)** | **2.63 mm**, DOY 95, *p* = 0.036 | `phaseG` |
| **Amplitude A−B** | **0.90 mm**, DOY 146, *p* = 0.448 | `phaseG` |
| **ROBUST ceiling (level 1)** | **≤ 3.9 mm vertical** (3.29 LOS ÷ cos 32.26°, all attributed to motion, NO lake assumption) | §4.3.7 |
| Refined bound (level 2) | < 2 mm LOS / **2.4 mm** vertical — **assumes stable lake** | §4.3.7 |
| Closure bias A | −0.090 rad, 1.6σ, *n* = 518 | `phaseG` |

## H4 — hydrology (§4.4)

| Quantity | Value | Source |
|---|---|---|
| NDWI(A) anomalies | *r* = +0.450, lag 12 d, *p* ≤ 0.011 | `phaseI` |
| NDWI(C) anomalies | *r* = +0.427, lag 12 d, *p* = 0.022 | `phaseI` |
| NDWI(D) anomalies | *r* = +0.424, lag 42 d, *p* ≤ 0.011 | `phaseI` |
| Differential A−C | *r* = −0.316, *p* = 0.150 | `phaseI` |
| **Temperature collapse** | −0.509 → **0.224**, *p* = 0.581 | `phaseI` |
| API / precipitation | not significant | `phaseI` |
| Null p95 (NDWI_A) | 0.404 | `phaseI` |

## External references used quantitatively

| Reference | Use |
|---|---|
| Hrysiewicz et al. 2024, *RSE* 291 | 10–40 mm breathing, C-band raised bogs — §5.1 |
| Patil et al. 2026, *RSASE* 41 | 0.48–1.40 cm yr⁻¹, drained Great Fen — §5.2 |
| Juszczak et al. 2013 | Rzecin hydrology (water table 0–30 cm) |
| De Zan 2015; Ansari 2021 | closure-phase bias / moisture |
| Ferretti et al. 2011; Fornaro et al. 2015 | distributed scatterers / phase linking |
| Kellndorfer et al. 2022 | global S1 coherence (external context) |

## Pre-submission checklist

- [x] Raise null draws to escape *p*-value floors *(Phase J: 280 draws)*
- [x] Snow/frost test *(Phase J3)*
- [x] Measure N_eff empirically *(Phase J1)*
- [x] Quantify Δ incidence A vs C *(Phase J2)*
- [ ] Export all figures at 300 dpi (`export_figures_en.ipynb`)
- [ ] Verify every number in the text appears in this file
- [ ] Obtain in-situ WTD series if available → strengthens §4.4
- [ ] Re-verify every reference against the publisher record
- [ ] Replace the schematic locator inset with a proper basemap if required
