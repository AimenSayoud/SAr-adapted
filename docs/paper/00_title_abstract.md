# What does C-band InSAR measure over a floating peatland? Multi-method evidence for a dielectric-dominated signal and the limits of satellite displacement retrieval

**Authors.** [Aymen Sayoud]¹, [Supervisor Name]¹, [Co-authors]

¹ [Institutional Affiliation, Department, University, City, Country]

**Corresponding author.** [author@institution.edu]

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
competing hypotheses.

**(H1) The failure is not algorithmic.** Six estimators with distinct
mathematical assumptions — SBAS, ISBAS, annual pairs, a hybrid network,
weighted least squares, and eigenvalue-decomposition phase linking — fail
identically over standard burst interferograms. Phase linking, while
theoretically optimal under an unconstrained covariance matrix, recovers only
5.4 % of mat pixels at temporal coherence ≥ 0.7, against 64.7 % for
land-cover-matched vegetation on stable ground using the same network.

**(H2) The mat is a distinct radar target.** At matched land cover and
phenology, its coherence is significantly lower (mean Δ = −0.081, lower in
89 % of pairs, date-jackknife leave-one-out range [−0.0842, −0.0774] demonstrating
sign invariance, SE-based 95 % CI [−0.109, −0.052]), its boundary is
sharp, its closure-phase dispersion is 3.2× larger, and environmental
predictors that actively modulate coherence inside the floating mat have no
detectable effect in matched grassland.

**(H3) The detectable seasonal signal is dielectric, not mechanical.** Spatial
aggregation over 499 pixels recovers a seasonal point estimate of 3.29 mm LOS
(95 % upper bound 7.32 mm LOS / 8.66 mm vertical; permutation p = 0.026 with
Monte Carlo 95 % CI [0.021, 0.031] against ≈4 600 reference-matched null
realisations) that six per-pixel inversions could not see. However, the
residual open-water lake exhibits a consistent seasonal trajectory
(2.63 mm, same phase) although it cannot breathe mechanically; the
mat-minus-lake difference cancels (0.90 mm, p = 0.45). Attributing the entire
signal to motion and propagating its 95 % interval yields an upper bound on
apparent phase-centre displacement of ≤ 8.7 mm vertical; because the phase
centre is not rigidly coupled to the mat, this does not exclude peat motion
of the amplitude published for raised bogs.

**(H4) What the sensor does track is surface wetness.** On deseasonalised
anomalies, the aggregated phase co-varies with Sentinel-2 optical wetness
(lag-0 r = 0.42–0.45 depending on reference zone, p ≤ 0.022) at a near-zero lag
of 12 days (one revisit), whereas air temperature does not survive
deseasonalisation (−0.509 → 0.224, p = 0.58). Coupling operates through a
sensitivity contrast between saturated peat and mineral grassland rather
than through a moisture contrast between them — a model prediction confirmed by
the expected failure of a differential forcing.

We conclude that the limitation is physical rather than methodological, and
report two transferable contributions: a change of observable (aggregate
rather than map, because the signal lies below the per-pixel noise floor) and
a weak-signal test protocol (size-matched nulls plus empirical p-values)
that invalidated two of our own intermediate conclusions and falsified two
explicit predictions.

**Keywords:** InSAR; Sentinel-1; C-band; floating peatland; *Schwingmoor*;
decorrelation; phase linking; distributed scatterers; surface moisture; spatial
aggregation.

---

## Highlights

- Six InSAR inversion strategies fail identically over a floating peatland: the
  limitation is physical, not algorithmic.
- Spatial aggregation recovers a 3.3 mm seasonal signal invisible pixel by
  pixel, but three independent tests indicate it is dielectric.
- Apparent phase-centre displacement is constrained to ≤ 8.7 mm vertical;
  mat motion itself is not bounded, because radar coupling is unmeasured.
- Aggregated phase responds to surface-moisture anomalies at near-zero lag
  through a sensitivity contrast between surfaces.
- A reproducible weak-signal protocol — size-matched nulls and empirical
  p-values — is demonstrated for decorrelated terrain.
