# Rzecin InSAR Visual Validation Package (QGIS)

> **Grid:** every raster here is on the 40 m HyP3 radar grid (138 × 129 px, EPSG:32633).
> `mean_coherence_spatial.tif`, `temporal_coherence_evd.tif`, `water_flooded_fraction.tif` and
> `elevation_glo30_dem.tif` were verified identical to `phaseD_coh_mean.nc`, `phaseE2_evd.nc`,
> the pipeline flooded fraction and the HyP3 DEM (2026-09-24). The script that originally
> produced this package was never committed; `scripts/build_web_atlas.py` now regenerates the
> same layers from the products with provenance.

This directory contains a complete, self-contained GIS validation package for the Rzecin peatland study (*What does C-band InSAR measure over a floating peatland?*). It allows researchers and reviewers to visually inspect all spatial domains, sensitivity masks, ground stations, topographic models, and radar coherence metrics in QGIS.

---

## 🚀 Quick Start in QGIS

1. **Open the Project File**:
   Double-click or open `rzecin_validation.qgs` directly in **QGIS** (version 3.22+ recommended).
   All layers use relative paths (`./<filename>`) and will load without missing datasource warnings.

2. **Add Basemap (Optional)**:
   - In QGIS, navigate to the **Browser** panel -> **XYZ Tiles**.
   - Add **OpenStreetMap** (`https://tile.openstreetmap.org/{z}/{x}/{y}.png`) or **Google Satellite** (`https://mt1.google.com/vt/lyrs=s&x={x}&y={y}&z={z}`) to see high-resolution optical imagery under the layers.

3. **Coordinate Reference System (CRS)**:
   - All spatial layers are projected in **WGS 84 / UTM zone 33N (`EPSG:32633`)**, which matches the native processing grid of the Sentinel-1 interferometric stack.

---

## 📂 Layer Catalog & Grouping

The QGIS project organizes layers into logical groups:

### 1. Vector Boundaries & Zones (`Vector Layers`)
- **`zones_polygons.geojson`**:
  - **Zone A**: floating mat — AOI minus water (499 px of 40 m = 79.84 ha, `T01_zones.csv`).
  - **Zone B**: residual lake — AOI water (65 px).
  - **Zone C**: matched grassland control outside the AOI — same dominant WorldCover class as A,
    S2 features within A's p10–p90, slope ≤ 5° (398 px).
  - **Zone D**: all other non-water, flat cover outside the AOI (10,750 px).
  - *Corrected 2026-09-24 (web-atlas audit): this entry previously said Zone A was
    "120 pixels, ~1.08 ha" and swapped the meaning of C and D. Counts are those of
    `T01_zones.csv`, reproduced exactly by `stratify.define_zones`.*
- **`lake_erosion_rings.geojson`** (*Experiment X-010*):
  - Sensitivity zones for open-water radar leakage:
    - **Border ring**: $d \le 40$ m from Zone A (39 pixels).
    - **Interior ring**: $d > 40$ m from Zone A (26 pixels).
    - **Deep center**: $d > 80$ m from Zone A (4 pixels).
  - Validates that the 2.22–2.90 mm seasonal signal observed over the lake is phase-locked ($\Delta = 1.5$ days) to the mat, ruling out edge contamination as the sole driver.
- **`marginal_27_pixels.geojson`** (*Experiment X-011 / Phase L7*):
  - The 27 surviving edge pixels under high-threshold spatial masking, used in the toroidal permutation test.
- **`field_monitoring_stations.geojson`**:
  - Exact coordinates of field infrastructure:
    - **Flux Tower**: ICOS Eddy Covariance Station (PL-Wet, Poznań University of Life Sciences).
    - **P1 & P2**: Automated groundwater level piezometers.
    - **Core-1**: Peat coring location. **Provenance unverified:** this README cites Milecka et al.
      2017, the GeoJSON cites Lamentowicz 2008, and no source records where any of the four
      station coordinates came from. Treat positions as approximate until confirmed.
- **`aoi_rzecin_boundary.geojson`**:
  - Overall Area of Interest (89.7 ha) enclosing the peatland-lake complex.

---

### 2. Radar Interferometry & Coherence (`Radar Coherence`)
- **`mean_coherence_spatial.tif`**:
  - Average spatial interferometric coherence ($\bar{\gamma}$) across all 356 Sentinel-1 interferometric pairs.
  - Clearly delineates the severe decorrelation over open water and floating vegetation compared to stable mineral ground.
- **`temporal_coherence_evd.tif`**:
  - Temporal coherence derived from the eigenvalue decomposition (EVD) of the multi-temporal covariance matrix.
  - Shows the loss of phase closure over the floating mat.

---

### 3. Topography & Hydrology (`Topography & Surface`)
- **`elevation_glo30_dem.tif`**:
  - Copernicus GLO-30 as delivered inside the HyP3 products, on the 40 m radar grid
    (identical to the HyP3 `*_dem.tif`). Values 93–141 m are WGS84-ellipsoidal heights
    (HyP3 convention), not a.s.l. — *corrected 2026-09-24; previously "10 m, ~50–65 m a.s.l."*
  - Shows the topographic bowl holding the peatland complex.
- **`slope_topography.tif`**:
  - Topographic slope in degrees, from the 40 m DEM. Zones A and B are flat to the DEM's
    resolution: slope ≤ 1.6° everywhere in A (median 0°), ≤ 1.0° in B. The DEM is stored in
    whole metres, so slope is quantised. *(Previously claimed "< 0.5°".)*
- **`distance_to_mat_boundary.tif`**:
  - Continuous Euclidean distance raster (in meters) from the boundary of the floating mat (Zone A).
  - Used directly to parameterize the ring erosion experiments in X-010.
- **`water_flooded_fraction.tif`**:
  - Sentinel-1 / optical-derived surface water frequency index, tracking permanent open water vs seasonally inundated *Sphagnum* lawns.
- **`esa_worldcover_rzecin.tif`**:
  - ESA WorldCover 2021 (10 m product) resampled to the 40 m radar grid.
- **`zones_classified_raster.tif`**:
  - Rasterized 40 m mask of Zones A (1), B (2), C (3), and D (4); counts match `T01_zones.csv`.

---

## 📦 Sharing & Portability

The file `rzecin_qgis_package.zip` (209 KB) contains all the GeoTIFFs, GeoJSONs, and the `.qgs` project file. You can share this zip file directly with collaborators or reviewers—extracting it anywhere retains full functionality.
