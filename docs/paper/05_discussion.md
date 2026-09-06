## 5. Discussion

### 5.1 Why C-band succeeds on raised bogs and fails here

Our results do not contradict Hrysiewicz et al. (2024), who retrieved bog
breathing over Irish raised bogs with correlations of 0.8–0.9. They **delimit
the domain of validity** of the method.

| | Raised bog | Floating fen (Rzecin) |
|---|---|---|
| Canopy | comparatively **dry** *Sphagnum* | *Sphagnum* + sedges, **saturated** |
| Water table | deeper, more variable | **near-surface**, stable |
| Scatterers | **surface**, stable | **volumetric**, non-stationary |
| Substrate | consolidated peat | mat resting on water |
| C-band coherence | usable | 5.4 % of pixels ≥ 0.7 |

The discriminating factor is not simply "peatland or not" but the combination of
**canopy wetness, scattering structure, and target stationarity**. Unlike raised bogs
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

**Generalisable contribution:** the success of C-band over peatlands does **not**
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
Great Fen. Direct comparison with our bound would be a **category error**: their
figure is a **rate** (mm yr⁻¹, vertical) while ours is the **amplitude** of an
annual cycle (mm, line-of-sight). The comparable quantity is our **velocity**,
which we did measure: −1.53 mm yr⁻¹ against a null of −1.50, i.e. not
significant, with a detection floor of ≈ 1.5–5 mm yr⁻¹.

Read that way, the comparison is informative:

| Site | Hydrological state | Subsidence |
|---|---|---|
| Late-restored farms (Great Fen) | drained, under restoration | 14.0 mm yr⁻¹ |
| Early-restored farms (Great Fen) | drained, under restoration | 11.7 mm yr⁻¹ |
| Nature reserves (Holme, Woodwalton) | conserved, **wetter** | **4.8 mm yr⁻¹** |
| **Rzecin (this study)** | **natural, saturated, never drained** | **not detected** |

The gradient follows hydrological state. Rzecin, never drained and with a
near-surface water table, sits **below the least-subsiding site** of that series
— the ecologically expected outcome rather than a measurement failure. Our
non-detection is therefore **consistent with** that literature, and the
comparison delimits where C-band peatland subsidence monitoring applies: drained
peat (agricultural surfaces, high coherence, centimetre-scale signal) rather than
saturated floating mats (wet canopy, low coherence, millimetre-scale signal).

**A methodological caution our results support.** Patil et al. interpret
seasonal fluctuations aligned with soil moisture as hydrological control of peat
surface motion. Our results counsel care with that inference: **a seasonal
oscillation correlated with moisture is not automatically motion.** At our site a
3.3 mm signal correlating well with moisture proved **dielectric** — the residual
lake, which cannot breathe, oscillated at the same amplitude and phase. Phase
variations induced by dielectric permittivity and moisture changes are well
established in the InSAR literature (De Zan et al., 2014, 2015; Morrison et al.,
2011; Zwieback et al., 2015, 2017; Mira et al., 2022; Zheng & Fattahi, 2025). A
**control over a water surface**, or any target where motion is physically
excluded, is inexpensive and separates genuine displacement from a
penetration-depth or dielectric effect. We suggest incorporating one systematically
in peatland motion studies reporting millimetre-scale signals.

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

**Methodological positioning.** The approach tested here is that of the
operational state of the art: the **OPERA DISP-S1** product (NASA/JPL) performs
hybrid persistent- and distributed-scatterer phase linking on the sample
coherence matrix — the same algorithmic core. It is not applicable here (North
America only, and C-band, hence the same limitation), but it confirms that the
route explored is not marginal.

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

**Measured water table.** Rzecin is an instrumented site; an in-situ WTD series
would transform the H4 test, since its temporal structure differs from that of
temperature.

### 5.6 Limitations

- **One site, one track, one polarisation** (VV): transferability of the
  predictive model (R²cv = 0.24) to other peatlands remains to be demonstrated.
- **No in-situ validation** (neither laser nor WTD): the constraint on
  apparent phase-centre displacement is internal to the InSAR analysis, and the
  coupling between that phase centre and the peat is unmeasured.
- **S1A-only window** (12-day revisit): a degraded cadence relative to what is
  now available.
- ***p*-values at the floor** (1/(1 + N)) for several tests; more null draws
  would tighten them.
- **Unit assumption** in aggregation: the mat is assumed to deform as a block,
  to be verified by subdivision should a mechanical signal appear.
- **Zone C is fragmented**, which biases the empirical correlation-length
  estimator used for N_eff.
