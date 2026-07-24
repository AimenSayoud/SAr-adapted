# Inventaire des figures

**Statut :** à produire (les cellules existent, les fichiers doivent être exportés)

Convention : exporter en **PNG 300 dpi** dans `outputs/article/` sous le nom
`FXX_slug.png`. Couleurs de zone constantes : **A** rouge `#d62728` ·
**B** bleu `#1f77b4` · **C** vert `#2ca02c` · **D** gris `#9e9e9e`
(`src/insar_wetlands/zone_viz.py`).

---

## Figures principales (article)

| # | Titre | Notebook source | Statut |
|---|---|---|---|
| **F1** | Localisation du site + emprise du burst S1 | `phase01_acquisition` | à produire |
| **F2** | Zones A/B/C/D — carte catégorielle **et** contours sur la cohérence | `phaseI` §2 (vues 1-2) | cellule prête |
| **F3** | Profils radiaux (cohérence, σ0, RVI) autour du bord | `phaseI` §2 (vue 5) | cellule prête |
| **F4** | Schéma du protocole : observable par pixel vs agrégée + nul apparié | — | **à dessiner** |
| **F5** | Cohérence temporelle par zone — histogrammes + courbe multi-seuils | `phaseE2`, `phaseH` §2 | cellule prête |
| **F6** | Carte de cohérence temporelle (phase-linking EVD) | `phaseE2` | cellule prête |
| **F7** | Distributions par zone des 5 champs | `phaseI` §2 (vue 4) | cellule prête |
| **F8** | Prédicteurs intra-tapis : coefficients + inversion A vs C | `phaseH` §4 | cellule prête |
| **F9** | Série agrégée A−C, A−B, B−C **et contrôle nul** | `phaseG` | cellule prête |
| **F10** | Test de significativité : amplitude saisonnière vs distribution nulle | `phaseG` (G2-ter) | cellule prête |
| **F11** | Fermeture de phase : biais et dispersion par zone | `phaseG` (G3) | cellule prête |
| **F12** | Phase agrégée vs humidité S2 (série temporelle) | `phaseI` §3 | cellule prête |
| **F13** | Corrélation saisonnière **vs** anomalies, par forçage | `phaseI` §3.3 | **à composer** |

## Figures supplémentaires (matériel additionnel)

| # | Titre | Source |
|---|---|---|
| S1 | Composite RVB (σ0 / cohérence / humidité) + zoom | `phaseI` §2 (vue 6) |
| S2 | Réseau d'interférogrammes (baselines temporelles/perpendiculaires) | `phase03_network_qc` |
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

## À faire

1. Ajouter une cellule d'export `plt.savefig(OUT/'FXX_slug.png', dpi=300,
   bbox_inches='tight')` dans chaque notebook concerné.
2. Dessiner **F4** (schéma conceptuel) — probablement le seul élément à produire
   hors notebook.
3. Composer **F13** (comparaison saisonnier/anomalies) à partir des deux tables
   de `phaseI`.
4. Uniformiser tailles de police et dimensions (colonne simple ~90 mm, double
   ~190 mm).
