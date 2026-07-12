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
| 2 | Interférogrammes HyP3 (burst InSAR, réseau SBAS) | — | à venir |
| 3 | Contrôle qualité & topologie du réseau (go/no-go cohérence) | — | à venir |
| 4 | Fusion multi-capteurs (S2 → grille UTM HyP3, nearest) | — | à venir |
| 5 | Masque dynamique de l'eau (NDWI/MNDWI/σ⁰) | — | à venir |
| 6 | Classification spatio-temporelle (5 classes) | — | à venir |
| 7 | Indice de qualité W_i & matrice de pondération | — | à venir |
| 8 | Inversion SBAS (MintPy) | — | à venir |
| 9 | Inversion ISBAS customisée (+ contrôle *phase closure*) | — | à venir |
| 10 | Comparaison quantitative SBAS vs ISBAS | — | à venir |
| 11 | Correction atmosphérique comparative (ERA5/GACOS vs réf. locale) | — | à venir |
| 12 | Décomposition LOS → vertical (d_vert = d_LOS / cos θ) | — | à venir |
| 13 | Analyse hydrologique & séparation signal diélectrique | — | à venir |
| 14 | Incertitudes & produits finaux (mm/an, RMSE, cartes) | — | à venir |

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
