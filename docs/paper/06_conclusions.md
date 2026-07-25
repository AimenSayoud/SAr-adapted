## 6. Conclusions

We assessed the applicability of C-band Sentinel-1 interferometry to vertical
displacement monitoring over a **floating peatland** — the configuration in
which the expected geomorphological signal is largest, and one that had not
previously been evaluated head-on.

**1. The limitation is physical, not algorithmic.** Six estimators resting on
distinct mathematical assumptions — up to maximum-likelihood phase linking —
fail identically. On the same 356-pair network the mat yields **5.4 %** of
pixels at temporal coherence ≥ 0.7, against **64.7 %** for land-cover-**matched**
vegetation on stable ground. The burden of proof shifts: there is no case for
seeking a seventh algorithm.

**2. The mat is a distinct radar target.** At matched cover and phenology its
coherence is significantly lower (Δ = −0.069, *p* = 2.2 × 10⁻⁴⁹, date-jackknife
sign-stable), its boundary is sharp, its scattering is volumetric, its
closure-phase dispersion is **3.2×** larger — and the drivers of its coherence
**reverse sign** relative to the adjacent grassland.

**3. The measurable seasonal signal is dielectric, and motion is bounded.**
Spatial aggregation recovers a seasonal amplitude of **3.29 mm** (*p* = 0.014)
invisible pixel by pixel. But the residual lake — which cannot breathe —
oscillates identically (2.63 mm, same phase); the mat-minus-lake difference
cancels (0.90 mm, *p* = 0.45); and the amplitude is one to two orders of
magnitude too small for flotation. **Seasonal vertical motion of the mat is
≤ 3.9 mm**, against 10–40 mm published for raised bogs.

**4. What the sensor tracks is surface wetness.** On deseasonalised anomalies
the aggregated phase co-varies with Sentinel-2 optical wetness (*r* = 0.42–0.45,
*p* ≤ 0.022) at **near-zero lag** (12 days), while temperature does not survive
deseasonalisation. Coupling operates through a **sensitivity contrast** between
saturated peat and mineral grassland — a prediction confirmed by the expected
failure of a differential forcing.

### Answer to the question posed

> **What does C-band Sentinel-1 InSAR actually measure over a floating
> peatland?**
>
> Not surface motion, which stays below 4 mm — but the **hydrological state** of
> that surface, through a **penetration-depth** effect. At C-band the floating
> peatland behaves as a **wet scattering volume with a quasi-random phase**,
> whose phase centre moves with moisture rather than with the substrate.

### Methodological contributions

- A **change of observable**: spatially aggregating a target that deforms as a
  unit drops the noise as 1/√N_eff and reveals a signal inaccessible pixel by
  pixel.
- A **weak-signal test protocol** — size-matched nulls, a null distribution
  rather than a single realisation, identical treatment of the null — which
  invalidated two intermediate conclusions and falsified two explicit
  predictions of our own.
- A **lightweight phase-linking implementation** applicable directly to standard
  interferometric products without an SLC processing chain.

### Outlook

The most direct route is **L-band**. At λ = 24 cm the radar penetrates the canopy
and reaches the mat surface. **NISAR** data are now global and free, and their
GUNW product is the direct analogue of our interferograms, so the chain developed
here applies without modification. The return of Sentinel-1 to a 6-day revisit
will reduce temporal decorrelation but **will not change the wavelength**, on
which volumetric decorrelation depends.

**In-situ validation** — surface laser and measured water table — remains
necessary to confirm the 3.9 mm bound and to separate, within the decorrelation
mechanism, dielectric variability from non-rigid micro-movement.

---

## Data and code availability

Processing code, notebooks and unit tests: [repository DOI]. Sentinel-1
interferograms were produced with ASF HyP3; Sentinel-2, ERA5, ESA WorldCover and
Copernicus DEM data are openly available from their respective providers.

## Author contributions

[CRediT statement]

## Declaration of competing interests

The authors declare no competing interests.

## Acknowledgements

[Funding, field support, data providers]
