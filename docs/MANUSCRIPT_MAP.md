# Manuscript reproduction map

Study figures and tables map to entries in `config/register.yaml`. These scripts require the inputs described in [the reproduction guide](REPRODUCIBILITY.md); the code alone does not reproduce the private study results.

| Manuscript Element | Generating Script(s) & Core Library Module(s) |
|---|---|
| **Data Cube (`raw → silver → gold`)** | `pipelines/cube_build.py`, `src/insar_wetlands/cube/{sources,build,gold,schema}.py` |
| **Sub-cell grid placement (`HyP3` & `RTC`)** | `src/insar_wetlands/placement.py`, `experiments/hyp3_half_cell_refine_x120.py`, `experiments/rtc_placement_x122.py`, `experiments/hyp3_smoothing_check_x124.py` |
| **Laser QC, snow-evidence mask & datum steps** | `src/insar_wetlands/field.py`, `experiments/field_p6_laser_qc.py`, `experiments/laser_frost_rule_x129.py`, `experiments/reread_laser_rule_x130.py`, `experiments/laser_qc_x147.py` |
| **Graphical Abstract** | `experiments/article_graphical_abstract.py` |
| **Fig. 1** (Site, CIR orthophoto, transect & radar geometry) | `experiments/article_fig_f1.py` |
| **Fig. 2** (Floating-mat hydro-mechanical motion & woody cover) | `experiments/article_fig_f2.py`, `experiments/buoyancy_x145.py`, `experiments/landcover_coherence_x132.py` |
| **Fig. 3** (Physical mechanism schematic across seasons) | `experiments/article_fig_f3.py` |
| **Fig. 4 & Table 2** (Consecutive 6-d vs 12-d pairs against `P6` laser) | `experiments/article_fig_f4.py`, `experiments/article_figures_c067.py`, `experiments/x065_hard_targets_x126.py`, `experiments/winter_snow_x146.py` |
| **Fig. 5** (Open-water coherence floor & wetness change `ΔVV`) | `experiments/article_fig_f5.py`, `experiments/lake_open_water_x064.py`, `experiments/coherence_moisture_x082.py`, `experiments/interior_edge_x128.py` |
| **Fig. 6** (Spatial footprint of 6-day mat motion) | `experiments/figures_landcover_footprint_c073.py`, `src/insar_wetlands/footprint.py`, `experiments/p6_motion_footprint_x133.py` |
| **Fig. 7 & Table 1** (Synthesis across pair lengths & tested routes) | `experiments/article_fig_f6.py`, `experiments/article_figures_c067.py`, `experiments/skip_pairs_x127.py`, `experiments/closure_x143.py`, `experiments/contrast_x144.py` |
| **Fig. 8** (Within-season chains vs nested winter re-levelling `SPAMS`) | `experiments/article_fig_f7.py`, `experiments/article_fig_f7b.py`, `src/insar_wetlands/{v4,spams}.py`, `experiments/nested_winter_x135.py`, `experiments/cross_track_x136.py` |
| **Multi-year season totals & basin drying (2017–2025)** | `src/insar_wetlands/drying.py`, `experiments/season_total_x137.py`, `experiments/drying_x138.py`, `experiments/drying_x139.py`, `experiments/drying_x140.py`, `experiments/drying_x141.py` |
| **Supplementary Figs (`S15`–`S20`)** | `experiments/fig_vegetation_set.py` (`S15`, `S16`), `experiments/drying_x141.py` (`S17`), `experiments/figures_landcover_footprint_c073.py` (`S19`, `S20`) |
