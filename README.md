# Pipeline InSAR pour zones humides — Tourbière de Rzecin (Pologne)

Suivi de la déformation verticale (« respiration » de la tourbière) par
interférométrie radar Sentinel-1, avec comparaison méthodologique
**SBAS vs ISBAS** et **correction atmosphérique globale vs calibration locale**.

- **Site :** tourbière de Rzecin (~90 ha), 52.763°N / 16.310°E
- **Données :** Sentinel-1 SLC bursts (via HyP3), Sentinel-2 L2A, ERA5
- **Exécution :** Google Colab (notebooks) + ce repo (code versionné)
- **Sorties :** poussées vers des branches `outputs/phaseXX`

## Utilisation dans Colab

1. Stocker les secrets dans Colab (icône 🔑) : `EARTHDATA_USERNAME`,
   `EARTHDATA_PASSWORD`, `CDSAPI_KEY`, `GITHUB_TOKEN`, `GIT_USER_EMAIL`.
2. Ouvrir le notebook de la phase courante (`notebooks/phaseXX_*.ipynb`).
3. La première cellule clone/pull le repo et installe l'environnement
   (`environment/colab_setup.sh`).
4. En fin de phase, les produits légers (`outputs/phaseXX/`) sont poussés
   vers la branche `outputs/phaseXX` ; les données lourdes restent sur
   Google Drive (`/content/drive/MyDrive/insar_rzecin`).

## Structure

```
config/config.yaml        # UN seul fichier de vérité (AOI, dates, seuils)
environment/colab_setup.sh
data/aoi/                 # polygone Rzecin (GeoJSON)
notebooks/                # un notebook mince par phase
src/insar_wetlands/       # tout le code réel (importé par les notebooks)
outputs/                  # produits par phase (non versionnés ici)
```

## Phases du pipeline

**Cette table est generee** depuis `config/phases.yaml` par `make readme`.
Ne pas l'editer a la main : la version ecrite a la main documentait 14 des 38
notebooks, et la moitie non documentee est celle qui porte l'argument de
l'article (D pour H2, G pour H3, J la falsification, L le gate round-2).

`make phases` valide la declaration : notebook non declare, dependance vers une
phase superseded, cycle, fichier absent.

<!-- PHASES:START -->
| Phase | Question | Status | Paper |
|---|---|---|---|
| [`phase01`](notebooks/phase01_acquisition.ipynb) — Multi-source acquisition inventory | What Sentinel-1, Sentinel-2 and ERA5 data exist over the AOI and period? | ✅ | §2.2, Table 2 |
| [`phaseB`](notebooks/phaseB_egms_check.ipynb) — External check against EGMS | Does an independent ground-motion service see anything here? | 🔍 exploratory | - |
| [`phaseF`](notebooks/phaseF_nisar_lband.ipynb) — L-band test (NISAR) | Does the mat become coherent at 24 cm? | 🔍 exploratory | §5.4 — prospective, the archive does not cover 2022-2024 |
| [`build_manuscript_docx`](notebooks/build_manuscript_docx.ipynb) — Manuscript assembly to .docx | — | 🔧 tooling | the manuscript itself |
| [`export_figures_en`](notebooks/export_figures_en.ipynb) — Manuscript figure export (English), 300 dpi | — | 🔧 tooling | all figures and T*.csv |
| [`phase01b`](notebooks/phase01b_selection.ipynb) — Pair selection — the paper's method, improved | Which interferometric pairs should the network contain? | ✅ | §2.2, Fig. 3 |
| [`phase15`](notebooks/phase15_annual_pairs.ipynb) — Annual "seasonal-annual" pairs | Does matching the seasonal state bypass temporal decorrelation? | ✅ | §4.1, Table 7 (estimator 3) |
| [`phase02`](notebooks/phase02_hyp3_jobs.ipynb) — HyP3 burst interferograms — SBAS network | Can the planned network be produced as burst interferograms in the cloud? | ✅ | §2.2, Table 2 |
| [`phase02b`](notebooks/phase02b_hyp3.ipynb) — HyP3 submission and download | Same, for the alternative pair set. | ✅ | - |
| [`phase03`](notebooks/phase03_network_qc.ipynb) — Network quality control and topology | Is C-band coherence over the mat sufficient to proceed at all? | ✅ | §2.2, Fig. 3 |
| [`phase03b`](notebooks/phase03b_pair_rates.ipynb) — Per-pair vertical rate, deramp and quality control | What vertical rate does each pair imply on its own? | ✅ | - |
| [`phase04`](notebooks/phase04_s2_fusion.ipynb) — Sentinel-2 fusion onto the HyP3 UTM grid | How do the optical observations align with the radar grid? | ✅ | §2.3, Table 3 |
| [`phase12`](notebooks/phase12_los_vertical.ipynb) — LOS to vertical decomposition | What is the measured incidence angle, and the conversion factor? | ✅ | §A.5 — the measured 32.26 deg, factor 1.183 |
| [`phaseA`](notebooks/phaseA_hybrid_insar_test.ipynb) — Hybrid network — the decisive InSAR test | Does combining short and long baselines recover usable pixels over the AOI? | ✅ | §4.1, Table 7 (estimator 4) |
| [`phaseC1`](notebooks/phaseC1_licsbas.ipynb) — LiCSBAS — an independent SBAS pipeline | Does a completely separate implementation close the loop? | 🔍 exploratory | - |
| [`phaseC2`](notebooks/phaseC2_phaselinking.ipynb) — Phase linking / DS-InSAR — first attempt | Does distributed-scatterer processing rescue the mat? | ⛔ superseded | superseded by phaseE2 |
| [`phaseE`](notebooks/phaseE_phaselinking_miaplpy.ipynb) — Phase linking with MiaplPy | Does the standard DS-InSAR toolchain run on these products? | ⛔ superseded | superseded by phaseE2 — MiaplPy needs an SLC archive we do not have |
| [`phaseE2`](notebooks/phaseE2_evd_phaselinking.ipynb) — Lightweight EVD phase linking — the direct test of H1 | Does the maximum-likelihood estimator, run on the delivered interferograms, recover the mat? | ✅ | §4.1.2, H1 — Table 8, Fig. 5, Fig. 6 |
| [`phase04b`](notebooks/phase04b_products.ipynb) — Ensemble, uncertainty and final products | What is the per-pixel uncertainty of the pairwise rates? | ✅ | - |
| [`phase05`](notebooks/phase05_water_mask.ipynb) — Dynamic water mask | Which pixels are inundated, and how often? | ✅ | §2.4, Fig. S2 |
| [`phase13`](notebooks/phase13_hydrology.ipynb) — Hydrological analysis and dielectric separation | Does the signal track water rather than motion? | ⛔ superseded | superseded by phaseI |
| [`phase14`](notebooks/phase14_products.ipynb) — Uncertainties and final products | What are the final maps and their uncertainties? | ✅ | - |
| [`phaseH`](notebooks/phaseH_predict_failure.ipynb) — Predicting where Sentinel-1 fails | Which covariates predict temporal coherence inside the mat? | ✅ | §4.2.5, Table 11, Fig. 10 |
| [`phase06`](notebooks/phase06_classification.ipynb) — Spatio-temporal classification, 5 classes | Do the pixels fall into distinct behavioural classes? | ✅ | - |
| [`phase07`](notebooks/phase07_quality_index.ipynb) — Quality index W_i | Can coherence and non-water fraction be combined into one weight? | ✅ | §3.1 |
| [`phaseD`](notebooks/phaseD_inside_vs_outside.ipynb) — Inside versus outside — is the floating mat a distinct target? | At matched land cover, does the mat decorrelate more than the same vegetation on stable ground? | ✅ | §4.2, H2 — Table 9, Fig. 7 |
| [`phase08`](notebooks/phase08_sbas_mintpy.ipynb) — SBAS inversion (MintPy) | Does a standard small-baseline inversion recover displacement? | ✅ | §4.1, Table 7 (estimator 1) |
| [`phase09`](notebooks/phase09_isbas.ipynb) — ISBAS inversion, with closure control | Does tolerating intermittent pixels rescue the inversion? | ✅ | §4.1, Table 7 (estimator 2) |
| [`phaseDbis`](notebooks/phaseDbis_mechanism_spatial.ipynb) — Mechanism and spatial comparison — interior, edge, exterior, lake | Is the boundary sharp, and where exactly does coherence step? | ✅ | §4.2.3, Fig. 9 |
| [`phaseDter`](notebooks/phaseDter_scattering_scatterers.ipynb) — Scattering mechanism and DS-InSAR feasibility | Is the mat radar-dark, or bright with unstable phase? | ✅ | §4.2.4, Table 10 |
| [`phaseG`](notebooks/phaseG_aggregation.ipynb) — Spatial aggregation — the change of observable | Does aggregating 499 pixels reveal a signal invisible pixel by pixel? | ✅ | §4.3, H3 — Table 14, Fig. 11 |
| [`phase10`](notebooks/phase10_comparison.ipynb) — SBAS versus ISBAS, quantitative | Do the two inversions disagree anywhere that matters? | ✅ | §4.1 |
| [`phase11`](notebooks/phase11_atmosphere.ipynb) — Comparative atmospheric correction — double-blind test | Does ERA5 correction beat a local reference zone? | ✅ | §A.4 |
| [`phaseI`](notebooks/phaseI_hydro_sensor.ipynb) — Is the dielectric signal a hydrological sensor? | On deseasonalised anomalies, does the aggregated phase track surface wetness? | ✅ | §4.4, H4 — Table 18, Fig. 14, Fig. 15 |
| [`phaseJ`](notebooks/phaseJ_falsification.ipynb) — Falsification — trying to demolish our own conclusion | Which of the eight alternative explanations survive? | ✅ | Appendix A — Table A1 |
| [`phaseZ`](notebooks/phaseZ_export_figures.ipynb) — Full figure and table export | — | 🔧 tooling | all figures and T*.csv |
| [`phaseK`](notebooks/phaseK_referee_response.ipynb) — Referee response | Do the referees' objections change any reported number? | ✅ | docs/paper/response_to_referees.md |
| [`phaseL`](notebooks/phaseL_gate.ipynb) — The gate, and the tests that matter either way | Do the nine round-2 tests hold? | ✅ | docs/paper/round2_report.md |
<!-- PHASES:END -->

Le cœur scientifique (matrice de design, inversion ISBAS, récupération des
pixels intermittents, *phase closure*) est validé par un test synthétique :
`python tests/test_synthetic_inversion.py`.

### Flux de données entre phases

- **Drive** (`/content/drive/MyDrive/insar_rzecin/`) : données lourdes et
  persistantes — produits HyP3 croppés (`hyp3_cropped/`), stacks netCDF
  (S2, RTC, masque d'eau, séries temporelles), artefacts inter-phases
  (`artifacts/` : sélection de paires, pixel de référence, inventaires).
- **Branches `outputs/phaseXX`** : produits légers (CSV, PNG, GeoTIFF finaux,
  `manifest.json`) poussés à la fin de chaque notebook.

## Corrections apportées au pipeline initial

1. **Phase 1-2 :** pas de téléchargement SLC — interférométrie par burst dans
   le cloud HyP3 (`INSAR_ISCE_BURST`), indispensable pour les quotas de
   crédits et les limites de Colab.
2. **Phase 4 :** les produits HyP3 sont **déjà géocodés (UTM)**, pas en
   géométrie radar ; la fusion S1/S2 est un rééchantillonnage nearest de S2
   vers la grille HyP3.
3. **Phase 3 :** ajout d'un go/no-go explicite sur la cohérence C-band
   (échantillon de paires test avant de soumettre tout le réseau).
4. **Phase 7 :** W_i = cohérence **conditionnelle hors-eau** × fraction
   hors-eau (évite la double pénalisation des pixels inondés).
5. **Phase 9-10 :** contrôle des sauts de phase des pixels « sauvés » par
   l'ISBAS via les résidus d'inversion et la fermeture des triplets
   (*phase closure*), qui sert aussi au biais diélectrique de la Phase 13.
