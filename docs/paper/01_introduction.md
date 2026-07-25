## 1. Introduction

### 1.1 Motivation

Peatlands cover roughly 3 % of the land surface but store about one third of
global soil organic carbon. Their carbon balance is tightly coupled to water
table depth, which in turn expresses itself mechanically: peat swells when
saturated and subsides when it dries, a reversible oscillation known as *bog
breathing* or *Mooratmung*. Measuring this motion therefore provides an
indirect, spatially distributed observation of hydrological condition, and by
extension of carbon vulnerability.

Spaceborne interferometric SAR is, in principle, the ideal instrument: free,
continuous, spatially resolved and millimetre-sensitive. Sentinel-1 has provided
systematic coverage since 2014, making retrospective and operational monitoring
feasible without field infrastructure.

### 1.2 State of the art, and the tension it contains

**C-band works on raised bogs.** Hrysiewicz et al. (2024) retrieved bog
breathing over Irish raised bogs from Sentinel-1, reporting correlations of
0.8–0.9 against in-situ data, RMSE ≈ 7 mm yr⁻¹, and recovered amplitudes of
10–40 mm. This result refutes the widespread assumption that C-band is
unusable over peat.

**C-band also works on drained peat.** Patil et al. (2026) mapped subsidence of
0.48–1.40 cm yr⁻¹ across the Great Fen (UK), a drained lowland peatland under
restoration, and used the deformation field as a proxy for carbon flux.

**But peatlands are not one class of target.** Raised bogs carry a comparatively
dry *Sphagnum* cover with stable surface scatterers; drained agricultural peat
presents bare or cropped surfaces with high coherence. **Floating peatlands**
(*Schwingmoor*, transitional fens) are different in every respect that matters to
a radar: a dense low canopy, permanent saturation, and a mat resting directly on
water. Nothing guarantees that success on the former transfers to the latter.

**The case has not been tested head-on.** The InSAR wetland literature is
dominated by inundation mapping through double-bounce, coherence as a land-cover
descriptor, and drained or exploited peatlands. The floating mat — precisely
where the expected geomorphological signal is largest — remains poorly
documented.

### 1.3 The methodological trap we set out to avoid

A study that fails to retrieve displacement can always be attributed to a poor
processing choice. The InSAR literature offers a wide catalogue of algorithms
(SBAS, persistent scatterers, phase linking / distributed scatterers, networks
of varying baseline design), and it is tempting, when faced with failure, to try
one more indefinitely.

We adopt the opposite stance: **if estimators resting on different mathematical
assumptions fail in the same way, the burden of proof shifts.** The explanation
is then unlikely to be algorithmic, and the scientific task becomes
characterising the target rather than searching for a seventh algorithm.

### 1.4 Question

We therefore do not ask

> *Can Sentinel-1 InSAR measure vertical displacement of the Rzecin floating
> mat?*

— a closed question whose answer is a bare yes or no — but rather

> **What does C-band Sentinel-1 InSAR actually measure over a floating
> peatland, and why?**

### 1.5 Design

The analysis is organised around four falsifiable, competing hypotheses

![**Figure 1.** Study design. Four competing hypotheses, the test applied to each, and the verdict reached.](figures/F01_hypotheses.png)

(Fig. 1):

| | Hypothesis | Principal test | Section |
|---|---|---|---|
| **H1** | The failure originates in the **inversion algorithm** | six independent estimators on one network | §4.1 |
| **H2** | The **mat** is a distinct radar target | matched-cover comparison + multi-sensor validation | §4.2 |
| **H3** | The residual signal is surface **motion** | spatial aggregation + lake, null and magnitude controls | §4.3 |
| **H4** | The measured signal reflects **hydrological state** | correlation on deseasonalised anomalies | §4.4 |

### 1.6 Contributions

1. A **multi-method demonstration** that the limitation is physical rather than
   algorithmic, at a site where the expected signal is maximal.
2. A **quantitative upper bound** (≤ 3.9 mm seasonal vertical amplitude),
   directly comparable to published raised-bog values.
3. Evidence, from three independent arguments, that the detectable seasonal
   signal is **dielectric**.
4. Two transferable methodological contributions: a **change of observable**
   (spatial aggregation) and a **weak-signal test protocol** (size-matched nulls
   with empirical *p*-values).
