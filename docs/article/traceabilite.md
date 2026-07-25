# Traçabilité — phase → résultat → section → figure

Chaque chiffre de l'article doit être retrouvable ici. Aucun nombre ne doit
figurer dans le texte sans une ligne correspondante.

---

## Correspondance phases → hypothèses

| Phase (notebook) | Ce qu'elle apporte | Section |
|---|---|---|
| `phase01`-`phase07` | acquisition, réseau, masques, fusion S2 | §2 |
| `phase08_sbas_mintpy` | SBAS — échec | §4 (H1) |
| `phase09_isbas` | ISBAS — échec | §4 (H1) |
| `phase15_annual_pairs` | paires annuelles — échec | §4 (H1) |
| `phaseA_hybrid_insar_test` | réseau hybride — 0 pixel fiable sur l'AOI | §4 (H1) |
| `phaseE2_evd_phaselinking` | **phase-linking EVD — test décisif** | §4 (H1) |
| `phaseD_inside_vs_outside` | Δcohérence apparié A vs C | §5 (H2) |
| `phaseDbis_mechanism_spatial` | hydro, gel, profil radial, résidus | §5 (H2) |
| `phaseDter_scattering_scatterers` | σ0, VH/VV, D_A, lac, phénologie | §5 (H2) |
| `phaseH_predict_failure` | modèle prédictif intra-tapis, inversion des moteurs | §5 (H2) |
| `phaseG_aggregation` | **agrégation, saisonnier, lac, fermeture** | §6 (H3) |
| `phaseI_hydro_sensor` | validation des zones + lien hydrologique | §2, §7 (H4) |
| `phaseF_nisar_lband` | chaîne bande L (prospectif) | §8 |
| `phaseB_egms_check`, `phaseC1_licsbas` | contrôles externes | §8 |
| **`phaseJ_falsification`** | **N_eff mesuré, incidence, neige/gel, plancher p** | **§10** |

---

## Chiffres clés

### Site et zones (§2)

| Chiffre | Valeur | Source |
|---|---|---|
| Surface documentée | 89.7 ha | littérature site |
| A / B / C / D | 499 / 65 / 398 / 10 750 px | `phaseI` §1 |
| Surfaces | 79.84 / 10.40 / 63.68 / 1 720 ha | `phaseI` §2 vue 0 |
| **A+B vs documenté** | **90.24 ha, +0.6 %** | `phaseI` §2 vue 0 |
| Interférogrammes | 356 paires, ~90 dates | `phase03` |

### H1 — algorithme (§4)

| Chiffre | Valeur | Source |
|---|---|---|
| tcoh médiane A / B / C / D | 0.604 / 0.584 / 0.734 / 0.639 | `phaseE2` |
| % ≥ 0.7 | 5.4 / 1.5 / 64.7 / 23.1 | `phaseE2` |
| Plancher de bruit (redondance ~4) | ≈ 0.55 | `tests/test_synthetic_phaselinking.py` |
| Résidu WLS A / C | 2.46 / 1.92 rad | `phaseDbis` |
| Réseau hybride | 14 238 « résolus », 0 fiable AOI | `phaseA` |
| Courbe multi-seuils | table §4.3 | `phaseH` §2 |

### H2 — cible (§5)

| Chiffre | Valeur | Source |
|---|---|---|
| **Δ apparié coh_A − coh_C** | **−0.069** | `phaseD` |
| Wilcoxon | p = 2.2 × 10⁻⁴⁹ | `phaseD` |
| Jackknife par date | [−0.0705, −0.0652], signe stable | `phaseD` |
| Cohérence A / B / C / D | 0.408 / 0.396 / 0.492 / 0.438 | `phaseI` §2 |
| σ0 VV A / B / C / D | −10.09 / −15.41 / −11.22 / −9.28 dB | `phaseI` §2 |
| RVI A / B / C / D | 0.914 / 0.993 / 0.881 / 1.045 | `phaseI` §2 |
| Humidité S2 A / B / C / D | −0.513 / +0.185 / −0.522 / −0.440 | `phaseI` §2 |
| Dispersion fermeture | A 0.683 / C 0.212 (**×3.2**) | `phaseG` G3 |
| D_A (% < 0.25) | A 59 %, C 69 %, B 3 % | `phaseDter` |
| R² validé croisé intra-tapis | **0.239 ± 0.022** (RF 0.326) | `phaseH` §4 |
| Coefficients nettoyés | σ0 −0.275, RVI −0.251, humidité −0.240 | `phaseH` §4 |
| VIF avant nettoyage | vh_vv_db 240.8, rvi 235.7 | `phaseH` §4 |
| Inversion des moteurs | σ0 −0.379 vs −0.008 ; altitude −0.168 vs +0.430 | `phaseH` §4 |

### H3 — mécanisme (§6)

| Chiffre | Valeur | Source |
|---|---|---|
| Validation synthétique | pixel −13.7 vs agrégat −19.8 (vérité −20) | `tests/test_aggregate.py` |
| \|R\| A / B / C / D | 0.234 / 0.395 / 0.569 / 0.426 | `phaseG` G1 |
| Vitesse A−C vs nul | −1.53 vs −1.50 mm/an | `phaseG` G2 |
| **Amplitude A−C** | **3.29 mm**, jour 104, p ≤ 0.011 (92 nuls) | `phaseG` G2-ter |
| Nul apparié | médiane 0.82, p95 1.93 mm | `phaseG` G2-ter |
| **Amplitude B−C (lac)** | **2.63 mm**, jour 95, p = 0.036 | `phaseG` G2-quater |
| **Amplitude A−B** | **0.90 mm**, p = 0.448 | `phaseG` G2-quinquies |
| **Plafond ROBUSTE (niveau 1)** | **≤ 4.2 mm vertical** (3.29 LOS ÷ cos39°, tout attribué au mouvement, SANS hypothèse sur le lac) | §6.8 |
| Borne affinée (niveau 2) | < 2 mm LOS / 2.6 mm vertical — **suppose le lac stable** | §6.8 |
| Borne convertie en vertical | < 2.6 mm (÷cos 39°, ×1.29) | `geometry.los_to_vertical` |
| Vitesse A−C vs nul | −1.53 vs −1.50 mm/an → **non significatif** | `phaseG` G2 |
| Plancher de détection en vitesse | ~1.5-5 mm/an selon réseau | `phaseG` G2 |
| Biais de fermeture A | −0.090 rad, 1.6 σ, n = 518 | `phaseG` G3 |

### H4 — hydrologie (§7)

| Chiffre | Valeur | Source |
|---|---|---|
| NDWI(A) anomalies | r = +0.450, lag 12 j, p ≤ 0.011 | `phaseI` §3.3 |
| NDWI(C) anomalies | r = +0.427, lag 12 j, p = 0.022 | `phaseI` §3.5 |
| NDWI(D) anomalies | r = +0.424, lag 42 j, p ≤ 0.011 | `phaseI` §3.5 |
| Différentiel A−C | r = −0.316, p = 0.150 | `phaseI` §3.4 |
| **t2m (effondrement)** | −0.509 → **0.224**, p = 0.581 | `phaseI` §3.3 |
| API / précip | non significatifs | `phaseI` §3.3 |
| Nul p95 (NDWI_A) | 0.404 | `phaseI` §3.3 |

---

## Références externes citées

| Référence | Usage |
|---|---|
| Hrysiewicz et al. 2024, *RSE* 291 | respiration 10-40 mm, C-band sur tourbières hautes — comparaison §8.1 |
| **Patil et al. 2026, *RSASE* 41, 101919** | subsidence 0.48-1.40 cm/an, Great Fen (RU) **drainé** — comparaison §8.1-bis |
| Juszczak et al. 2013 | hydrologie du site de Rzecin |
| De Zan et al. 2015 ; Ansari et al. 2021 | biais de phase de fermeture / humidité |
| Fornaro et al. | estimateurs de phase-linking (EMI/EVD) |
| Ferretti et al. 2011 | SqueeSAR / diffuseurs distribués |
| Kellndorfer et al. 2022, *Sci. Data* | cohérence globale S1 (contrôle externe) |
| OPERA DISP-S1 (NASA/JPL) | positionnement méthodologique §4.4, §8 |
| NISAR / ASF | perspective bande L §8.5 |

---

## Vérifications à refaire avant soumission

- [ ] Augmenter les tirages nuls à 200-500 (sortir des planchers de p-value)
- [ ] Exporter toutes les figures en 300 dpi (voir `figures.md`)
- [ ] Vérifier que chaque nombre du texte figure dans ce fichier
- [ ] Obtenir la série **WTD in situ** si disponible → renforce §7
- [ ] Relire les étiquettes [FAIT]/[INT]/[HYP] et les retirer de la version finale
- [ ] **Exécuter le test neige/gel** : `filter_pairs(..., exclude_months=(12,1,2))`
- [ ] **Mesurer N_eff empiriquement** : `effective_looks(mask, field=coh_mean)`
- [ ] **Quantifier Δincidence A vs C** depuis la couche `lv_theta`
