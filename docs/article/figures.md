# Inventaire des figures

**Statut :** export automatisé — `notebooks/phaseZ_export_figures.ipynb`

Convention : **PNG 300 dpi** dans **`docs/article/figures/`** sous le nom
`FXX_slug.png`, via `zone_viz.save_figure()`. Couleurs de zone constantes :
**A** rouge `#d62728` · **B** bleu `#1f77b4` · **C** vert `#2ca02c` ·
**D** gris `#9e9e9e` (`src/insar_wetlands/zone_viz.py`).

> **Pourquoi `docs/article/figures/` et non `outputs/`** : `outputs/` est ignoré
> par git (convention des branches `outputs/phaseXX`). Les figures d'article
> doivent vivre **à côté du texte qui les cite**, versionnées avec lui — sans
> quoi une figure et le chiffre qu'elle illustre peuvent diverger.

## Régénération

`notebooks/phaseZ_export_figures.ipynb` reconstruit les figures depuis les
**produits sauvegardés sur Drive** (`phaseD_coh_mean.nc`, `phaseE2_evd.nc`,
`s2_stack.nc`, `rtc_dualpol_stack.nc`…) plutôt qu'en tout recalculant, et met en
cache les séries agrégées dans `figures_cache/`. Les figures exigeant des
distributions nulles (F10, F11, F13) restent dans leurs notebooks d'origine :
elles sont trop longues pour un export de routine.

---

## Figures principales (article)

| # | Titre | Notebook source | Statut |
|---|---|---|---|
| **F1** | Localisation du site + emprise du burst S1 | `phase01_acquisition` | à produire |
| **F2** | Zones A/B/C/D — carte catégorielle **et** contours sur la cohérence | **`phaseZ`** | ✅ automatisé |
| **F3** | Profils radiaux (cohérence, σ0, RVI) autour du bord | **`phaseZ`** | ✅ automatisé |
| **F4** | Schéma du protocole : observable par pixel vs agrégée + nul apparié | **`phaseZ`** | ✅ **dessiné** |
| **F5** | Cohérence temporelle par zone — histogrammes + courbe multi-seuils | **`phaseZ`** | ✅ automatisé |
| **F6** | Carte de cohérence temporelle (phase-linking EVD) | **`phaseZ`** | ✅ automatisé |
| **F7** | Distributions par zone des 5 champs | **`phaseZ`** | ✅ automatisé |
| **F8** | Prédicteurs intra-tapis : coefficients + inversion A vs C | `phaseH` §4 | cellule prête |
| **F9** | Série agrégée A−C, A−B, B−C **et contrôle nul** | **`phaseZ`** | ✅ automatisé |
| **F10** | Test de significativité : amplitude saisonnière vs distribution nulle | `phaseG` (G2-ter) | cellule prête |
| **F11** | Fermeture de phase : biais et dispersion par zone | `phaseG` (G3) | cellule prête |
| **F12** | Phase agrégée vs humidité S2 (série temporelle) | **`phaseZ`** | ✅ automatisé |
| **F13** | Corrélation saisonnière **vs** anomalies, par forçage | `phaseI` §3.3 | **à composer** |

## Figures supplémentaires (matériel additionnel)

| # | Titre | Source |
|---|---|---|
| S1 | Composite RVB (σ0 / cohérence / humidité) + zoom | **`phaseZ`** ✅ |
| S2 | Réseau d'interférogrammes + distribution des baselines | **`phaseZ`** ✅ |
| S3 | Validation synthétique : par pixel vs agrégé | `tests/test_aggregate.py` |
| S4 | Dispersion d'amplitude D_A par zone | `phaseDter` |

---

## Tableaux principaux

| # | Titre | Section |
|---|---|---|
| **T1** | Jeux de données et zones (effectifs, surfaces, validation +0.6 %) | §2 |
| **T2** | Cohérence temporelle et % ≥ 0.7 par zone (H1) | §4 |
| **T3** | Signature multi-capteurs par zone (cohérence, σ0, RVI, humidité) | §2 / §5 |
| **T4** | Amplitudes saisonnières A−C, B−C, A−B + nuls (H3) | §6 |
| **T5** | Corrélations saisonnier vs anomalies par forçage (H4) | §7 |
| **T6** | Erreurs corrigées et prédictions falsifiées | §8 |

---

## Reste à produire hors de `phaseZ`

| # | Figure | Pourquoi hors export |
|---|---|---|
| F1 | Localisation du site | nécessite un fond cartographique externe |
| F8 | Prédicteurs Phase H | relancer `phaseH` puis exporter |
| F10 | Significativité (nuls) | distribution nulle longue → `phaseG`/`phaseJ` |
| F11 | Fermeture de phase | énumération des triplets → `phaseG` |
| F13 | Saisonnier vs anomalies | balayage de décalage → `phaseI` |

La police et les marges sont déjà homogénéisées via `mpl.rcParams` en tête de
`phaseZ` (9 pt, 300 dpi, fond blanc).
