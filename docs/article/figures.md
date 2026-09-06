# Inventaire des figures

**Statut :** export automatisé — `notebooks/06_manuscript/phaseZ_export_figures.ipynb`

Convention : **PNG 300 dpi** dans **`docs/article/figures/`** sous le nom
`FXX_slug.png`, via `zone_viz.save_figure()`. Couleurs de zone constantes :
**A** rouge `#d62728` · **B** bleu `#1f77b4` · **C** vert `#2ca02c` ·
**D** gris `#9e9e9e` (`src/insar_wetlands/zone_viz.py`).

> **Pourquoi `docs/article/figures/` et non `outputs/`** : `outputs/` est ignoré
> par git (convention des branches `outputs/phaseXX`). Les figures d'article
> doivent vivre **à côté du texte qui les cite**, versionnées avec lui — sans
> quoi une figure et le chiffre qu'elle illustre peuvent diverger.

## Régénération

`notebooks/06_manuscript/phaseZ_export_figures.ipynb` reconstruit les figures depuis les
**produits sauvegardés sur Drive** (`phaseD_coh_mean.nc`, `phaseE2_evd.nc`,
`s2_stack.nc`, `rtc_dualpol_stack.nc`…) plutôt qu'en tout recalculant, et met en
cache les séries agrégées dans `figures_cache/`. Les figures exigeant des
distributions nulles (F10, F11, F13) restent dans leurs notebooks d'origine :
elles sont trop longues pour un export de routine.

---

## Figures principales (article)

| # | Titre | Section | Statut |
|---|---|---|---|
| **F01** | Réseau interférométrique + distribution des baselines | §2 | ✅ `phaseZ` |
| **F02** | Zones A/B/C/D — carte + contours sur la cohérence | §2 | ✅ `phaseZ` |
| **F03** | Profils radiaux (cohérence, σ0, RVI) au bord | §5 | ✅ `phaseZ` |
| **F04** | Schéma du protocole (observable + nul apparié) | §3 | ✅ `phaseZ` *(dessiné)* |
| **F05** | Cohérence temporelle par zone + multi-seuils | §4 | ✅ `phaseZ` |
| **F06** | Carte de cohérence temporelle (EVD) | §4 | ✅ `phaseZ` |
| **F07** | Distributions par zone — 5 capteurs | §5 | ✅ `phaseZ` |
| **F08** | Prédicteurs intra-tapis + **inversion A vs C** | §5 | ✅ `phaseZ` |
| **F09** | Séries agrégées A−C, B−C, A−B + nul | §6 | ✅ `phaseZ` |
| **F10** | Significativité + robustesse hiver | §6 | ✅ `phaseZ` *[LONG]* |
| **F11** | Fermeture de phase : biais et dispersion | §6 | ✅ `phaseZ` *[LONG]* |
| **F12** | Phase agrégée vs humidité S2 | §7 | ✅ `phaseZ` |
| **F13** | Saisonnier **vs** anomalies par forçage | §7 | ✅ `phaseZ` |

## Figures supplémentaires

| # | Titre | Section | Statut |
|---|---|---|---|
| **S01** | Composite RVB + zoom tourbière | §2 | ✅ `phaseZ` |
| **S02** | Fraction inondée (masque eau) | §2 | ✅ `phaseZ` |
| **S03** | Validation synthétique : par pixel vs agrégé | §3 | ✅ `phaseZ` |
| **S04** | Dispersion d'amplitude D_A (carte + zones) | §5 | ✅ `phaseZ` |
| **S05** | Décroissance de cohérence + distributions | §4/5 | ✅ `phaseZ` *[LONG]* |
| **S06** | **Test apparié A−C** (histogramme + nuage) | §5 | ✅ `phaseZ` |
| **S07** | Couplage hydrologique et test de gel | §5 | ✅ `phaseZ` |

## Tableaux exportés (CSV, à côté des figures)

| # | Contenu |
|---|---|
| **T1** | Zones : effectifs, surfaces, écart aux 89.7 ha |
| **T2** | Courbe multi-seuils de cohérence temporelle |
| **T3** | Signature multi-capteurs par zone (médianes) |
| **T4** | Test apparié A−C (Wilcoxon, jackknife) |
| **T5** | Couplage hydrologique / gel par zone |
| **T6** | Prédicteurs A vs C (inversion des moteurs) |
| **T7** | Amplitudes saisonnières des quatre séries |
| **T8** | Fermeture de phase par zone |
| **T9** | Forçages : saisonnier vs anomalies |

---

## Reste hors `phaseZ`

| # | Figure | Pourquoi |
|---|---|---|
| — | Carte de localisation du site | nécessite un fond cartographique externe |

Tout le reste est automatisé. Les cellules `[LONG]` (cohérence par paire ~5 min,
nuls ~15 min, triplets ~3 min) sont **mises en cache** dans `figures_cache/` :
la seconde exécution est immédiate.

La police, le dpi et les couleurs de zone sont fixés une seule fois en tête de
`phaseZ` (`mpl.rcParams` + `ZONE_COLORS`), donc toutes les figures sont
homogènes.
