# Rzecin first-study GIS context package

A dated visual context for the first Sentinel-1 study (2022–2024), rather than
the current field-validation results or the complete multi-sensor cube.
Rasters use the 40 m HyP3 grid in EPSG:32633. Vector CRS declarations are stored
in each GeoJSON; the AOI outline is supplied in geographic coordinates.

[Preview the public GIS layers](https://github.com/AimenSayoud/SAr-adapted/blob/main/docs/PRESENTATION.md#the-public-gis-layers).
The preview uses the bundled raster values without resampling and shows no
station points. It is descriptive first-study context, not displacement.

Open `rzecin_validation.qgs` in QGIS, or extract `rzecin_qgis_package.zip` to a
new folder. The project uses relative paths. The ZIP mirrors the supplied GIS
files and documents. See [layer provenance](PROVENANCE.md) before interpreting
any point or raster as independent ground truth.

## Contents

| Layers | Purpose |
|---|---|
| AOI, zones and classified zone raster | Study extent and first-study comparison zones |
| Mean spatial coherence and EVD temporal coherence | Descriptive radar-quality context |
| GLO-30 elevation and slope | Topographic context on the radar grid |
| ESA WorldCover and flooded fraction | Land-cover/water context |
| Distance, lake erosion rings and marginal cells | First-study sensitivity regions |
| Station reference points | Approximate orientation only; coordinate provenance is unverified |

The station file contains the tower and a coring reference point. Their exact
coordinate sources have not been established. A publication associated with a
site does not establish the origin or precision of these coordinates. These
points are unsuitable for field navigation, validation or precise extraction.
The public subset omits the unverified piezometer locations.

The original package-generation script was not preserved. The current
`pipelines/build_web_atlas.py` exports related layers with provenance; that is
not a claim that it recreates every byte of this historical QGIS package.
For current study entry points, see [reproduction guidance](../../docs/REPRODUCIBILITY.md).

The software licence applies to code. Source-provider attribution and terms
remain relevant to bundled/derived data; this package does not assert a new
blanket licence over third-party source products.
