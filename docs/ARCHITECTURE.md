# Architecture & Reproducibility Guide

`SAr-adapted` implements a modular, test-first scientific pipeline for evaluating C-band Sentinel-1 InSAR over the **Rzecin poor-fen floating mat** (*Schwingmoor*, western Poland, 2017–2026).

---

## 1. Layered Design

1. **Multi-Sensor Data Cube (`src/insar_wetlands/cube/`, `pipelines/cube_build.py`)**
   - **`raw`**: Native-format acquisitions (ASF HyP3 Sentinel-1 burst interferograms on relative orbits 175 ascending and 22 descending, OPERA COMPASS coregistered SLCs, Planetary Computer `sentinel-1-rtc`, Sentinel-2 L2A, Landsat 8/9 C2 L2, GUGiK 1 m LiDAR & CIR orthophotos, Open-Meteo ERA5/ERA5-Land, IMGW synoptic stations, NISAR L2 GUNW).
   - **`silver`**: Sensor-level harmonised stacks on the common 40 m UTM Zone 33N (`EPSG:32633`) grid with sub-cell geometric relocation applied before resampling.
   - **`gold`**: Analysis-ready cubes joining radar observables, optical indices, dated canopy structure, meteorology, and quality flags.

2. **Sub-Cell Grid Relocation (`src/insar_wetlands/placement.py`)**
   - Corrects the half-cell pixel-centre vs upper-left tag convention in `hyp3_isce2` (`[-20 m, +20 m]`) and the half-pixel tag shift in Planetary Computer `sentinel-1-rtc` (`[+5 m, -5 m]`), verified against 1 m LiDAR canopy edges (`X-109`, `X-120`, `X-122`).

3. **Common Ground-Truth Scorer (`src/insar_wetlands/scoring.py`, `config/protocol.yaml`)**
   - Every candidate estimator is evaluated by the same protocol (`Protocol.load()`):
     - **Within-season skill**: RMS error, correlation, and slope against continuous `P6` surface-elevation laser measurements inside snow-free seasons (`May–November`) within single sensor-datum segments.
     - **Across-winter step**: Post-winter (`March–June`) level error across winters not interrupted by sensor re-levelling.
     - **Pairwise revisit contrast**: Block-bootstrapped correlation and phase-to-laser slope as a function of temporal baseline (`6 d`, `12 d`, `18–365 d`).

4. **Tested Processing Routes (`src/insar_wetlands/`)**
   - Multi-looked network inversions (`inversion/sbas_mintpy.py`, `inversion/isbas.py`, `inversion/phaselinking.py`).
   - Coregistered single-look complex (`CSLC`) persistent scatterers, SHP adaptive multilooking, and EVD phase linking (`slc_stats.py`, `experiments/slc_*`).
   - Moisture-phase / closure bias models (`phase_bias.py`, `phase_sim.py`, `experiments/bias_correction_x085.py`, `experiments/closure_x143.py`).
   - Consecutive 6-day chains and nested weather-driven winter re-levelling (`spams.py`, `v4.py`, `experiments/nested_winter_x135.py`).

5. **Execution Register (`config/register.yaml`, `src/insar_wetlands/register.py`)**
   - Declares every pipeline, experiment, input, output, and git code reference.
   - Validated automatically via `make check` (`python -m insar_wetlands.register check`).
