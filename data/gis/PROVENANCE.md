# GIS layer provenance and limitations

This inventory documents the source identity available for the retained
first-study context package. Exact acquisition dates or original production
commands that were not recorded are labelled as such, rather than inferred.

| Files | Recorded source / processing | Temporal or provenance limit |
|---|---|---|
| `aoi_rzecin_boundary.geojson` | Study AOI outline | Boundary survey/source precision not established here; CRS declared in file |
| `zones_polygons.geojson`, `zones_classified_raster.tif` | First-study radar-grid zones; mat/lake mask and matched outside comparison | Describes the first-study zone definition, not later cube revisions |
| `mean_coherence_spatial.tif` | Mean HyP3 interferometric spatial coherence, 2022–2024 study stack | Processing summary, not a motion/accuracy measurement |
| `temporal_coherence_evd.tif` | First-study EVD temporal coherence | Historical inversion diagnostic, not field validation |
| `elevation_glo30_dem.tif` | Copernicus GLO-30 DEM supplied with HyP3, on the radar grid | HyP3 ellipsoidal-height convention; not a recent mat-surface survey |
| `slope_topography.tif` | Slope derived from that gridded DEM | Limited by DEM/grid resolution and height quantisation |
| `esa_worldcover_rzecin.tif` | ESA WorldCover 2021 resampled to radar grid | Class map epoch differs from later field campaigns |
| `water_flooded_fraction.tif` | First-study Sentinel-1/optical water-mask product | Composite frequency; not surveyed open-water extent at one date |
| `distance_to_mat_boundary.tif` | Distance from first-study mat boundary | Follows that mask, not later basin-edge definitions |
| `lake_erosion_rings.geojson`, `marginal_27_pixels.geojson` | First-study sensitivity selections | Analysis selections, not independent reference observations |
| `field_monitoring_stations.geojson` | Approximate tower/coring labels retained for orientation | Exact coordinate sources and precision unverified; not suitable as truth |
| `rzecin_validation.qgs`, `rzecin_qgis_package.zip` | Visual project/package assembled for the first study | Original generation command not preserved; ZIP is a convenience copy |

Provider names identify the origin of source products; consult the relevant
provider's terms before redistributing or reusing its data. The MIT software
licence does not replace those terms. Unverified station coordinates are
explicitly labelled in the GeoJSON properties as well as this documentation.
