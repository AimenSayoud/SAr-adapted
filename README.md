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

| # | Phase | Notebook | Statut |
|---|-------|----------|--------|
| 1 | Acquisition multi-source (inventaires S1/S2 + ERA5) | `phase01_acquisition.ipynb` | ✅ prêt |
| 2 | Interférogrammes HyP3 (burst InSAR, lot test + réseau SBAS) | `phase02_hyp3_jobs.ipynb` | ✅ prêt |
| 3 | Contrôle qualité & topologie du réseau (go/no-go cohérence) | `phase03_network_qc.ipynb` | ✅ prêt |
| 4 | Fusion multi-capteurs (S2 → grille UTM HyP3, nearest) | `phase04_s2_fusion.ipynb` | ✅ prêt |
| 5 | Masque dynamique de l'eau (NDWI/MNDWI/σ⁰ RTC) | `phase05_water_mask.ipynb` | ✅ prêt |
| 6 | Classification spatio-temporelle (5 classes + réf. Classe A) | `phase06_classification.ipynb` | ✅ prêt |
| 7 | Indice de qualité W_i & matrice de pondération | `phase07_quality_index.ipynb` | ✅ prêt |
| 8 | Inversion SBAS (MintPy) | `phase08_sbas_mintpy.ipynb` | ✅ prêt |
| 9 | Inversion ISBAS customisée (+ contrôle *phase closure*) | `phase09_isbas.ipynb` | ✅ prêt |
| 10 | Comparaison quantitative SBAS vs ISBAS | `phase10_comparison.ipynb` | ✅ prêt |
| 11 | Correction atmosphérique comparative (ERA5 vs réf. locale) | `phase11_atmosphere.ipynb` | ✅ prêt |
| 12 | Décomposition LOS → vertical (d_vert = d_LOS / cos θ) | `phase12_los_vertical.ipynb` | ✅ prêt |
| 13 | Analyse hydrologique & séparation signal diélectrique | `phase13_hydrology.ipynb` | ✅ prêt |
| 14 | Incertitudes & produits finaux (mm/an, RMSE, cartes) | `phase14_products.ipynb` | ✅ prêt |

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
