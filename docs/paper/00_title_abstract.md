# C-band InSAR over a floating peatland: target decorrelation, controlled contrast, and a non-unique hydrological phase response

**Authors.** Aymen Sayoud¹*, Supervisor Name¹, Co-authors¹

¹ Institutional Affiliation, Department, University, City, Country

*Corresponding author: aimen.sayoud.polska@gmail.com

---

## Abstract

Peatland surfaces oscillate vertically with the water table, and this "bog
breathing" is a useful proxy for hydrological condition and, indirectly, for
carbon status. Sentinel-1 interferometry offers free, continuous monitoring of
such motion, and C-band retrievals have been demonstrated on raised bogs.
Whether this transfers to floating peatlands (*Schwingmoor*) — wetter, with
a dense low canopy resting on a near-surface water table — has not been tested
directly, even though these are the sites where the expected geomorphological
signal is largest.

We assess C-band InSAR over the Rzecin floating fen (Poland, 89.7 ha) using 356
Sentinel-1 interferograms (2022–2024), structuring the analysis around four
falsifiable hypotheses.

**H1 (algorithmic failure — rejected).** Under identical processing, standard
burst interferometry recovers 64.7 % of pixels in land-cover-matched mineral
grassland but collapses to 5.4 % in the floating mat. Six inversion estimators
spanning pairwise multi-looking to eigenvalue-decomposition phase linking fail
identically, demonstrating that per-pixel failure is a target property rather
than an algorithmic artifact.

**H2 (distinct radar target — supported).** At matched land cover and
phenology, the mat exhibits lower coherence (mean Δ = −0.081, lower in
89 % of pairs, date-jackknife leave-one-out range [−0.0842, −0.0774] demonstrating
sign invariance, SE-based 95 % CI [−0.109, −0.052]), a sharp boundary, 3.2×
larger closure-phase dispersion, and environmental sensitivity unobserved in
matched grassland (where effective sample size is qualified by $N_{\text{eff}} \approx 5$).

**H3 (surface motion — rejected).** Spatial aggregation over 499 pixels recovers
a seasonal point estimate of 3.29 mm LOS (95 % upper bound 7.32 mm LOS / 8.66 mm
vertical; permutation p = 0.026 with Monte Carlo 95 % CI [0.021, 0.031] against
4,614 reference-matched null realisations) invisible pixel by pixel. Forward
dielectric modeling shows this amplitude can arise entirely from seasonal water-content
variations in the upper peat moss without mechanical displacement. The adjacent
open-water lake control is inconclusive due to spatial filter leakage
($L_{\text{corr}} \approx 160\text{ m}$) across its small basin and an elevated
minimum detectable amplitude ($2.86\text{--}3.02\text{ mm}$). Attributing the entire
differential signal to motion under pure-vertical attribution yields an upper bound
on differential apparent phase-centre displacement between mat and matched grassland
of ≤ 8.7 mm vertical; because the phase centre is not rigidly coupled to the mat,
this does not exclude published raised-bog breathing amplitudes.

**H4 (hydrological state — consistent with surface wetness).** On deseasonalised
anomalies, the aggregated phase co-varies with Sentinel-2 optical wetness
(lag-0 r = 0.39–0.42, p ≤ 0.022) within the Sentinel-1 12-day sampling resolution,
whereas air temperature does not survive deseasonalisation (−0.509 → 0.224, p = 0.58).
Because buoyant hydrostatic coupling is itself rapid, this 12-day association
is consistent with a rapid dielectric response while not excluding a rapid
mechanical response. Coupling operates through a sensitivity contrast between
saturated peat and mineral grassland rather than through a moisture contrast
between them.

The displacement-retrieval limitation persists across all tested standard
burst-product processing approaches and is consistent with strong target
non-stationarity at C-band. We report two transferable contributions: a
change of observable (aggregate rather than map, because the signal lies below
the per-pixel noise floor) and a weak-signal test protocol (size-matched nulls
plus empirical p-values) that invalidated four of our own intermediate
conclusions and falsified two explicit predictions.

**Keywords:** InSAR; Sentinel-1; C-band; floating peatland; *Schwingmoor*;
decorrelation; phase linking; distributed scatterers; surface moisture; spatial
aggregation.

---

## Highlights

- Controlled contrast under identical processing: standard burst InSAR recovers 64.7 % of pixels on stable grassland but collapses to 5.4 % on the floating peat mat, demonstrating target-specific decorrelation.
- Spatial aggregation recovers a 3.3 mm seasonal signal invisible pixel by
  pixel; the response is consistent with dielectric variation and cannot be
  uniquely interpreted as mechanical displacement.
- Differential apparent phase-centre displacement between mat and grassland is
  constrained to ≤ 8.7 mm vertical under pure-vertical attribution; mat motion
  itself is not bounded, because radar coupling is unmeasured.
- Aggregated phase responds to surface-moisture anomalies within the 12-day
  sampling cycle through a sensitivity contrast between surfaces.
- A reproducible weak-signal protocol — size-matched nulls and empirical
  p-values — is demonstrated for decorrelated terrain.
