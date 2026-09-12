## 6. Conclusions

We assessed the applicability of C-band Sentinel-1 interferometry to vertical
displacement monitoring over a floating peatland — the configuration in
which the expected geomorphological signal is largest, and one that had not
previously been evaluated head-on.

**1. The displacement-retrieval failure persists across tested standard burst-product approaches.** Six
estimators resting on distinct mathematical assumptions — up to eigenvalue-decomposition
phase linking — fail identically on pairwise multi-looked products. On the same 356-pair
network the mat yields 5.4 % of pixels at temporal coherence ≥ 0.7, against 64.7 % for
land-cover-matched vegetation on stable ground. Holding processing constant, this demonstrates
that per-pixel failure is a target property under standard burst processing; operational gains
require full-covariance single-look algorithms (e.g. SqueeSAR) or L-band observations rather
than an alternative pairwise burst inversion.

**2. The mat is a distinct radar target.** At matched cover and phenology its
coherence is significantly lower (mean Δ = −0.081, lower in 89 % of pairs,
date-jackknife leave-one-out range [−0.0842, −0.0774] sign invariant; 95 % CI [−0.109, −0.052]),
its boundary is sharp, its scattering is volumetric, its closure-phase dispersion is 3.2×
larger — and environmental predictors that actively modulate coherence in the mat have
no detectable effect in matched grassland.

**3. The seasonal signal is not uniquely mechanical, and motion remains unconstrained.**
Spatial aggregation recovers a seasonal point estimate of 3.29 mm LOS (95 % upper bound
7.32 mm LOS / 8.66 mm vertical; permutation p = 0.026 with Monte Carlo 95 % CI [0.021, 0.031]
against ≈4 600 reference-matched draws) invisible pixel by pixel. But the residual
open-water lake exhibits a consistent phase and amplitude trajectory (2.63 mm, within 9 days);
the mat-minus-lake difference cancels (0.90 mm, p = 0.45). Attributing the entire differential
signal to vertical motion gives an upper bound on differential apparent phase-centre displacement
between mat and grassland of ≤ 8.7 mm vertical. Because the phase centre is not rigidly coupled
to the peat, this does not bound mat motion: at a coupling fraction of 0.5 it admits 17 mm,
within the 10–40 mm published for raised bogs. What C-band establishes here is an inability
to distinguish absence of motion from peatland-scale breathing.

**4. What the sensor tracks is surface wetness.** On deseasonalised anomalies
the aggregated phase co-varies with Sentinel-2 optical wetness (lag-0 r = 0.42–0.45,
p ≤ 0.022) with an association occurring within Sentinel-1 12-day sampling, while temperature
does not survive deseasonalisation. Because hydrostatic coupling in a buoyant mat is itself
rapid, this 12-day association is consistent with a rapid dielectric response without excluding
rapid mechanical response. Coupling operates through a sensitivity contrast between saturated
peat and mineral grassland — a prediction confirmed by the expected failure of a differential
forcing.

### Answer to the question posed

> **What does C-band Sentinel-1 InSAR actually measure over a floating
> peatland?**
>
> Across standard burst interferograms, C-band InSAR measures a severely
> decorrelated volume whose phase cannot be resolved into per-pixel surface
> motion by tested pairwise or phase-linking algorithms. When spatially
> aggregated over the entire mat, the radar yields a reproducible seasonal
> phase (3.29 mm LOS) that co-varies with optical surface moisture anomalies
> at sub-revisit lag. However, because an identical seasonal trajectory appears
> over the adjacent open-water lake and cancels differentially, this observable
> cannot be uniquely attributed to vertical peat breathing. The available data
> are consistent with a moisture-dependent phase-centre migration in the top
> scattering volume, while mechanical peatland breathing remains unconstrained
> by satellite InSAR alone without in-situ datum anchoring.

### Methodological contributions

- A change of observable: spatially aggregating a target that deforms as a
  unit drops the noise as 1/√N_eff and reveals a signal inaccessible pixel by
  pixel.
- A weak-signal test protocol — size-matched nulls, a null distribution
  rather than a single realisation, identical treatment of the null — which
  invalidated two intermediate conclusions and falsified two explicit
  predictions of our own.
- A lightweight phase-linking implementation applicable directly to standard
  interferometric products without an SLC processing chain.

### Outlook

The most direct route is L-band. At λ = 24 cm the radar penetrates the canopy
and reaches the mat surface. NISAR provides global, open L-band observations.
The conceptual analysis is adaptable to NISAR GUNW products, although differences
in product geometry, spatial resolution, pairing strategy and processing require
validation before the chain transfers. The return of Sentinel-1 to a 6-day revisit
(ESA, 2024) will reduce temporal decorrelation but will not change the wavelength,
on which volumetric decorrelation depends.

**In-situ validation** — surface laser and measured water table — is now the
only route to constraining the phase-centre coupling on which any statement
about mat motion depends, and to separating, within the decorrelation mechanism,
dielectric variability from non-rigid micro-movement.

---

## Data and code availability

Processing code, configuration files, notebooks, and unit tests are available
in the repository at <https://github.com/AimenSayoud/SAr-adapted> (version `v1.0.0`,
archived on Zenodo: <https://doi.org/10.5281/zenodo.14999999>). Sentinel-1
interferograms were produced with ASF HyP3; Sentinel-2, ERA5, ESA WorldCover, and
Copernicus DEM data are openly available from their respective providers.

## Author contributions

**A. Sayoud:** Conceptualization, Methodology, Software, Validation, Formal analysis, Data Curation, Writing – original draft, Visualization. Co-authors to be credited for supervision, resources, and review & editing upon final manuscript submission.

## Declaration of competing interests

The authors declare no competing interests.

## Acknowledgements

The authors thank the Alaska Satellite Facility DAAC for processing Sentinel-1 interferograms via HyP3, the European Space Agency and Copernicus Programme for Sentinel-1 and Sentinel-2 imagery, and ECMWF for ERA5 atmospheric reanalysis data. We thank Poznań University of Life Sciences and the Rzecin Wetland research station teams for long-term site stewardship and environmental monitoring data.
