## 5. Discussion

### 5.1 Why C-band succeeds on raised bogs and fails here

Our results do not contradict Hrysiewicz et al. (2024), who retrieved bog
breathing over Irish raised bogs with correlations of 0.8–0.9, nor Alshammari et al. (2018)
or Tampuu et al. (2020), who demonstrated InSAR surface motion tracking across Scottish
blanket bogs and Estonian peatlands. Instead, our findings delimit the domain of validity
of C-band interferometry across peatland types.

| | Raised bog | Floating fen (Rzecin) |
|---|---|---|
| Canopy | comparatively dry *Sphagnum* | *Sphagnum* + sedges, saturated |
| Water table | deeper, more variable | near-surface, stable |
| Scatterers | surface, stable | volumetric, non-stationary |
| Substrate | consolidated peat | mat resting on water |
| C-band coherence | usable | 5.4 % of pixels ≥ 0.7 |

The discriminating factor is not simply "peatland or not" but the combination of
canopy wetness, scattering structure, and target stationarity. Unlike raised bogs
where a consolidated, comparatively drier *Sphagnum* carpet provides stable surface
scatterers (Hrysiewicz et al., 2024), the floating mat (*Schwingmoor*) at Rzecin is
characterized by near-permanent saturation and low-stature emergent vegetation
resting directly on an organic water layer (Milecka et al., 2017). Furthermore,
unlike flooded emergent wetlands where coherent double-bounce scattering between
vertical vegetation and underlying standing water preserves interferometric phase
(Hong & Wdowinski, 2014), the low-stature, heterogeneous vegetation and saturated
floating root mat at Rzecin provide a less stable scattering geometry. Because
perpendicular baselines in the Sentinel-1 stack are small ($B_\perp < 100$ m),
spatial decorrelation is negligible; the observed coherence collapse across 12-day
revisits is predominantly driven by temporal changes in surface moisture,
micro-topographic water pooling, and vegetation configuration—a pattern consistent
with observations across non-forested herbaceous wetlands (Chen et al., 2020, 2021).

**Generalisable contribution:** the success of C-band over peatlands does not
transfer automatically to floating peatlands, which are nonetheless the sites
where the expected vertical displacement from buoyancy is largest.

### 5.2 Comparison with drained peatlands, and a caution

Floating fens are known to exhibit vertical mobility coupled to water level changes
through root-mat buoyancy (Stofberg et al., 2016; Swarzenski et al., 1991; Sasser
et al., 1996). Stofberg et al. (2016) demonstrated on temperate floating fens that
root mats track surface water fluctuations hydrostatically, with buoyancy mediated
by temperature and gas dynamics. However, translating this physical mechanism to
satellite InSAR requires distinguishing true mechanical motion from radar phase
artifacts.

Patil et al. (2026) report subsidence of 0.48–1.40 cm yr⁻¹ over the drained
Great Fen, and Ghezelayagh et al. (2024) observed coherent seasonal subsidence over
drained agricultural sections of the Biebrza fen peatlands in northeastern Poland.
Direct comparison between their subsidence rates and our seasonal amplitude bound would
be a category error: their figure is a secular rate (mm yr⁻¹, vertical) while ours is
the amplitude of an annual cycle (mm, line-of-sight). The comparable quantity is our
velocity, which we did measure: −1.53 mm yr⁻¹ against a null of −1.50, i.e. not
significant, with a detection floor of ≈ 1.5–5 mm yr⁻¹.

Read that way, the comparison is informative:

| Site | Hydrological state | Subsidence |
|---|---|---|
| Late-restored farms (Great Fen) | drained, under restoration | 14.0 mm yr⁻¹ |
| Early-restored farms (Great Fen) | drained, under restoration | 11.7 mm yr⁻¹ |
| Nature reserves (Holme, Woodwalton) | conserved, wetter | 4.8 mm yr⁻¹ |
| Biebrza fens (Ghezelayagh et al., 2024) | drained/managed fen sections | seasonal subsidence detected |
| Rzecin (this study) | natural, saturated, never drained | not detected |

The gradient follows hydrological state. Rzecin, never drained and with a
near-surface water table, sits below the least-subsiding sites of those series
— the ecologically expected outcome rather than a measurement failure. Our
non-detection is therefore consistent with that literature, and the
comparison delimits where C-band peatland subsidence monitoring applies: drained
peat (agricultural surfaces, high coherence, centimetre-scale signal) rather than
saturated floating mats (wet canopy, low coherence, millimetre-scale signal).

**A methodological caution our results support.** Patil et al. interpret
seasonal fluctuations aligned with soil moisture as hydrological control of peat
surface motion. Our results counsel care with that inference: a seasonal
oscillation correlated with moisture is not automatically motion. At our site,
a 3.3 mm signal correlating well with moisture was accompanied by a consistent
phase and amplitude trajectory over the open-water lake, while the differential
mat-minus-lake signal cancelled. Phase variations induced by dielectric permittivity
and moisture changes are well established in the InSAR literature (De Zan et al.,
2014, 2015; Morrison et al., 2011; Zwieback et al., 2015, 2017; Mira et al., 2022;
Zheng & Fattahi, 2025).

A control over a water surface, or any target where mechanical breathing is physically
precluded, is inexpensive and separates genuine displacement from differential
propagation phase or dielectric permittivity effects (De Zan et al., 2014). We
suggest incorporating one systematically in peatland motion studies reporting
millimetre-scale signals.

**Quantitative comparison with the moisture-phase framework.** Extending the lossy
dielectric half-space framework of De Zan et al. (2014) to saturated peat (permittivity
mixing via a Birchak refractive model for *Sphagnum* organic solids, water and air at
λ = 5.55 cm, θ = 32.3°), C-band penetration depth is constrained to merely 3–4 mm into
the saturated capitulum layer. Under this model, assuming a seasonal volumetric moisture
excursion of $\Delta m_v \approx 0.25$ (from saturated $m_v = 0.85$ down to $0.60$) predicts
an apparent interferometric line-of-sight displacement of $-3.28$ mm.

Because in-situ dielectric moisture excursions were unmeasured during 2022–2024, this
alignment should not be interpreted as exact calibration, but rather as an order-of-magnitude
physical consistency check. Over a plausible range of moisture excursions $\Delta m_v \in [0.15, 0.35]$,
the predicted LOS displacement envelope spans $-1.9$ mm to $-4.6$ mm, consistent with the
observed 3.29 mm point estimate (Table 7) and well within the 95 % upper bound of 7.32 mm LOS.

Two critical properties accompany this quantitative alignment:

1. *Asymptotic ceiling*: The theoretical maximum LOS phase shift under complete desiccation
   converges asymptotically to ≈ 6.13 mm LOS. While our semi-amplitude (3.29 mm) aligns
   cleanly with plausible seasonal moisture variations, the full peak-to-peak swing
   (6.57 mm) slightly exceeds this ceiling, suggesting either that sinusoidal fitting
   over-estimates the amplitude of a non-linear saturating response, or that thermal
   expansion and temperature-dependent water permittivity contribute an additional
   ±1 mm seasonal component.
2. *Decorrelation falsification*: For $\Delta m_v = 0.25$, the dielectric propagation model
   predicts an interferometric pair coherence $|\gamma| = 0.725$. In our dataset, the observed
   mean pair coherence over Zone A is merely 0.408 (Table 4). This establishes that dielectric
   moisture variation alone cannot account for the severe decorrelation observed, isolating
   vegetation volume scattering and non-stationary scatterer dynamics as the dominant
   decorrelation mechanisms.

Furthermore, testing maximum temporal baseline subsets (24 d to 120 d and all-pair, Table 3)
reveals that while short-baseline networks (≤ 24 d) suffer from severe accumulating
closure-phase subsidence bias (−13.5 to −23.5 mm yr⁻¹), consistent with Zheng et al. (2022),
expanding the network to ≥ 48 d and annual pairs causes the apparent velocity to
contract to near-zero (−1.53 mm yr⁻¹ on A−C vs −1.50 mm yr⁻¹ on NULL), demonstrating that
our velocity non-detection is robust to network truncation.

![**Figure 16.** Our bound in context: raised-bog breathing, drained-fen subsidence, expected free flotation, and the value measured here.](figures/F16_literature_context.png)

### 5.3 Two transferable methodological contributions

#### 5.3.1 Change of observable

Per-pixel phase noise (≈ 1.5 rad at γ = 0.4) falls as 1/√N_eff under
aggregation. Over 499 pixels it drops to ≈ 0.3 mm (≈ 1 mm at N_eff = 50), far
below the signal sought. **The signal was not below the noise floor; it was below
the *per-pixel* noise floor.**

This reasoning applies to any target that is **spatially coherent but temporally
decorrelated**: peatlands, rock glaciers, wetlands, crops. The condition is that
the target deform as a unit — an assumption that must be physically justified and
is testable by subdividing the zone.

#### 5.3.2 Weak-signal test protocol

Three rules, each of which invalidated an intermediate conclusion in this study:

1. **Size-matched nulls.** Aggregate noise falls as 1/√N, so a null four times
   larger carries half the noise and **manufactures false detections**.
2. **A null distribution, not a single realisation.** One realisation is not a
   test; *N* draws give an empirical *p*-value — whose **floor** of 1/(1 + N)
   must be stated.
3. **Identical treatment of the null.** If the observed statistic results from a
   selection (best |r| over 16 lags), the null must undergo the same sweep.

These rules are cheap and should accompany any weak-signal claim over
decorrelated terrain.

### 5.4 Instrumental outlook

**L-band is the most direct route.** At λ ≈ 24 cm, radar signals penetrate
herbaceous vegetation more effectively, significantly mitigating the rapid
temporal decorrelation that affects C-band (λ = 5.5 cm) over dynamic, saturated
wetland covers (Chen et al., 2021; Morishita & Hanssen, 2015). **NISAR** now
provides **global, free** L-band data, and its **GUNW** product is the direct
analogue of our interferograms, so the processing chain developed here transfers
without modification. This is a **prospective** test: the archive does not
cover 2022–2024 retrospectively.

**The return to a 6-day revisit will not suffice.** Sentinel-1C and 1D restore
two-satellite operation, reducing temporal baseline. However, all Sentinel-1
platforms remain **C-band**: shorter repeat intervals partially mitigate temporal
decorrelation but do not alter the high sensitivity of shorter wavelengths to
micro-scale canopy and moisture changes. We therefore do not expect a 6-day cycle
alone to unlock this site.

**Historical alternative.** ALOS-2/PALSAR-2 (L-band, 2015–2024) would allow a
retrospective test over our window, but access is restricted.

**Methodological positioning.** The approach tested here reflects the operational
state of the art for multi-looked burst products: the OPERA DISP-S1 product (NASA/JPL)
performs hybrid persistent- and distributed-scatterer phase linking on the sample
coherence matrix — the same algorithmic core evaluated here. Advanced full-covariance
approaches utilizing raw Single Look Complex (SLC) data with Statistically Homogeneous
Pixel (SHP) selection (e.g., SqueeSAR; Ferretti et al., 2011) and phase bias mitigation
(Ansari et al., 2021) represent an important prospective avenue. While out of scope
for standard burst products, they offer a methodological benchmark for future high-performance
computing implementations.

#### Literature synthesis: moisture-induced phase and peatland InSAR

The physical mechanisms governing moisture-induced phase variations, decorrelation,
and InSAR observables across wetland and peat substrates are synthesized across six
foundational studies in Table 11.

| Study | System & Substrate | Observable & Method | Key Finding & Transfer to This Work |
|---|---|---|---|
| De Zan et al. (2014) | L-band; bare agricultural soil | Analytical lossy dielectric half-space (Born approx.) | Differential propagation phase in top lossy layer mimics deformation (≈ 10° for Δm_v = 0.01). Framework extended here to saturated peat. |
| Morrison et al. (2011) | C-band; indoor sand bed | Laboratory DInSAR vs physical surface laser | DInSAR phase shift far exceeds true surface motion during wetting/drying, proving dielectric phase masquerades as motion under controlled conditions. |
| Zwieback et al. (2015) | L-band; agricultural soils | Empirical regression on in-situ moisture | Soil moisture effect exceeds 2 cm apparent displacement for Δm_v = 0.20 (> 70 % of fields), bounding dielectric signal scales. |
| De Zan & Gomba (2018) | L-band; vegetated terrain | Closure-phase inversion constrained by coherence | Joint inversion of vegetation and moisture from closure phase; highlights C-band vegetated peat as an unresolved frontier. |
| Morishita & Hanssen (2015) | L-, C-, X-band; drained peat pasture | 3-parameter temporal decorrelation model | Demonstrates severe C-band decorrelation on peat compared to L-band, establishing physical need for L-band systems. |
| Zheng et al. (2022) | C-band; multi-temporal Sentinel-1 | Closure-phase accumulation & baseline subset test | Multilooking closure phase creates systematic velocity bias (~cm yr⁻¹) in short-baseline networks; verified in our subset test (Table 3). |

### 5.5 What remains open

**The decorrelation mechanism.** Two families remain compatible: (a) **dielectric**
variability of saturated peat, and (b) **non-rigid micro-movement** (local
flexure, sub-pixel differential settling). A third — **rigid-body** motion coupled
to the water table — is excluded (no hydrological coupling of coherence, no
stabilisation on freezing, no double-bounce signature). The closure-phase bias,
which was intended to separate (a) from (b), detects no systematic bias: high
dispersion **without a sign bias** is compatible with both.

This question is **distinct** from that of the seasonal signal, whose dielectric
origin is established (§4.3).

**In-situ validation.** A ground-based **laser** remains the only direct
measurement able to constrain the phase-centre coupling on which any
statement about mat motion depends, and to separate (a) from (b).

**Measured water table.** Rzecin is an instrumented station; however, continuous
in-situ piezometer time series were unavailable across the full retrospective
2022–2024 Sentinel-1 processing window due to institutional data governance
policies and sensor maintenance intervals. Consequently, optical surface moisture
proxies (Sentinel-2 NDWI/NMDI) were utilized as the primary continuous empirical
covariate. Direct in-situ water-table depth records remain the ideal validation
target for future prospective campaigns.

### 5.6 Limitations

- **One site, one track, one polarisation** (VV): transferability of the
  predictive model (R²cv = 0.24) to other peatlands remains to be demonstrated.
- **No in-situ validation** (neither laser nor continuous WTD during 2022–2024):
  the constraint on apparent phase-centre displacement is internal to the InSAR
  analysis, and the coupling between that phase centre and the peat is unmeasured.
- **S1A-only window** (12-day revisit): a degraded cadence relative to what is
  now available.
- ***p*-values at the floor** (1/(1 + N)) for several tests; more null draws
  would tighten them.
- **Unit assumption** in aggregation: the mat is assumed to deform as a block,
  to be verified by subdivision should a mechanical signal appear.
- **Zone C is fragmented**, which biases the empirical correlation-length
  estimator used for N_eff.
