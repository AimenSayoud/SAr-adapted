## 4. Results

### 4.1 H1 — The failure is not algorithmic

#### 4.1.1 Six estimators, one outcome

| Method | Outcome over the mat |
|---|---|
| SBAS (MintPy) | no usable pixel |
| ISBAS (Alshammari et al., 2018) | idem; intermittent pixels not recovered |
| Annual pairs | idem |
| Hybrid network | 14 238 pixels "resolved" scene-wide, 0 reliable over the AOI (median residual 2.5 rad) |
| Weighted least squares | median residual 2.46 rad (A) vs 1.92 rad (C), at comparable numbers of valid pairs |
| Phase linking (EVD) | see §4.1.2 |

The hybrid-network case is instructive: the apparent 78-fold scene-wide coverage
gain did not yield reliable pixels over the wetland area of interest. A coverage
criterion without a reliability criterion is misleading.

#### 4.1.2 The decisive test

Phase linking is theoretically optimal under an unconstrained, complete, and unbiased
sample covariance matrix. In our network of 356 pairs across ~90 dates, the pairwise
products populate 8.89 % (356 of 4,005 off-diagonal pairs) of the full covariance
structure. Evaluating phase linking directly on these delivered burst interferograms
provides a rigorous assessment of whether standard operational products support
displacement retrieval over the mat.

**Table 2** — Temporal coherence by zone (356 pairs, ~90 dates):

| Zone | Median | p25–p75 | % ≥ 0.7 |
|---|---|---|---|
| C — matched grassland | 0.734 | 0.671–0.803 | 64.7 % |
| D — other cover | 0.639 | 0.597–0.693 | 23.2 % |
| A — floating mat | 0.604 | 0.566–0.647 | 5.4 % |
| B — residual lake | 0.584 | 0.542–0.630 | 1.5 % |

Read against the simulated network noise floor of ≈ 0.55 (empirical network
distribution floor 0.488 with a 90 % interval of [0.448, 0.532]; §3.1; Fig. 5,
Fig. 6). Because that floor is strongly topology-dependent, we express each zone
as its excess above the floor, which requires no threshold:

| Zone | Temporal coherence | Excess over floor |
|---|---|---|
| C — matched grassland | 0.734 | 0.246 |
| D — other cover | 0.639 | 0.151 |
| A — floating mat | 0.604 | 0.116 |
| B — residual lake | 0.584 | 0.096 |

The mat retains 47 % of the matched grassland's excess coherence above the
floor. Every zone, including the lake, lies above the 95th percentile of the
null, so the lake is not at the floor and cannot serve as an internal
validation of the chain. The mat is low, and intermediate between the lake and
external cover — deprived of its high-coherence tail to the point where
per-pixel inversion is not supportable, while retaining measurable structure
above the fully decorrelated case.


![**Figure 5.** Phase-linking result. (a) Temporal-coherence distributions by zone, with the 0.55 noise floor and the 0.7 reliability threshold; (b) full multi-threshold curve.](figures/F05_temporal_coherence.png)

![**Figure 6.** Map of temporal coherence from EVD phase linking, with zone outlines.](figures/F06_tcoh_map.png)

#### 4.1.3 Multi-threshold analysis

The 0.7 threshold is a convention; the full curve is more informative (Fig. 5b,
Table T03). A and C are nearly indistinguishable at 0.50 (0.974 vs 0.995) and
diverge in the upper tail (≥ 0.65). The mat is therefore not uniformly
shifted downward — it is deprived of its best pixels, which is precisely
what prevents inversion.

#### 4.1.4 Baseline subsets and closure-phase accumulation

Testing maximum temporal baseline subsets ($B_{t,\max} \in [24, 36, 48, 60, 120, \text{All}]$ days;
Table 3) confirms that this inversion failure is not an artifact of network connectivity. While short
baselines ($\le 24\text{ d}$) suffer severe closure-phase accumulation bias (Zheng et al., 2022), the
full network stabilizes reliably, demonstrating that per-pixel decoherence is intrinsic to the target
rather than the baseline selection scheme.

**Table 3** — Baseline subset analysis and closure-phase accumulation test ($B_{t,\max} \in [24, 120]$ days and All pairs; cf. Zheng et al., 2022):

| Series | Baseline subset | $B_{t,\max}$ (d) | Pairs | Velocity (mm yr⁻¹) | Amplitude (mm) | Phase (DOY) | Trend (mm yr⁻¹) | $R^2$ |
|---|---|---|---|---|---|---|---|---|
| A−C (mat vs grassland) | $\le 24$ d | 24 | 175 | −13.47 | 9.50 | 117 | −11.55 | 0.74 |
| A−C (mat vs grassland) | $\le 36$ d | 36 | 261 | −8.50 | 4.94 | 104 | −7.44 | 0.38 |
| A−C (mat vs grassland) | $\le 48$ d | 48 | 346 | −3.87 | 2.89 | 106 | −3.25 | 0.22 |
| A−C (mat vs grassland) | $\le 60$ d | 60 | 346 | −3.87 | 2.89 | 106 | −3.25 | 0.22 |
| A−C (mat vs grassland) | $\le 120$ d | 120 | 346 | −3.87 | 2.89 | 106 | −3.25 | 0.22 |
| A−C (mat vs grassland) | All pairs | All | 356 | −1.53 | 3.29 | 104 | −0.83 | 0.30 |
| B−C (lake vs grassland) | $\le 24$ d | 24 | 175 | −23.48 | 9.43 | 124 | −21.68 | 0.55 |
| B−C (lake vs grassland) | All pairs | All | 356 | −1.22 | 2.63 | 95 | −0.65 | 0.11 |
| NULL (grassland vs stable) | $\le 24$ d | 24 | 175 | −1.87 | 1.37 | 337 | −2.03 | 0.22 |
| NULL (grassland vs stable) | All pairs | All | 356 | −1.50 | 0.57 | 95 | −1.38 | 0.06 |

#### 4.1.5 Verdict: H1 rejected for burst interferometric networks

> **H1 is rejected for standard burst interferometric networks and pairwise multi-looked products.**
> Six estimators with distinct mathematical assumptions fail identically. The inversion failure is
> a physical property of the target under C-band multi-looked burst observation.

**Scope.** Future work exploiting full-covariance Single Look Complex (SLC) stacks with
statistically homogeneous pixel (SHP) selection (e.g. SqueeSAR; Ferretti et al., 2011; Ansari et al.,
2018) could optimize covariance estimation on a subset of the archive. However, for standard
multi-looked burst networks, the relative gap between A and C persists: C succeeds where A fails,
under identical network topology and identical processing.

### 4.2 H2 — The mat is a distinct radar target

#### 4.2.1 Coherence deficit at matched cover

**Table 4** — paired comparison, A vs C:

| Quantity | Value |
|---|---|
| Mean coherence A / C | 0.408 / 0.492 |
| **Paired Δ (coh A − coh C), mean** | **−0.081** |
| Paired Δ, median | −0.050 |
| Fraction of pairs with A lower | 89 % |
| **Date-jackknife** (on the mean Δ) | leave-one-out range [−0.0842, −0.0774] (sign invariant); 95 % CI [−0.109, −0.052] |
| Wilcoxon signed-rank | *p* = 4.84 × 10⁻⁴⁶ (nominal; see below) |
| Decorrelation time τ | 21 d (A) vs 32 d (C) |

The mean exceeds the median in magnitude (−0.081 against −0.050), so the
distribution of paired differences is left-skewed: a subset of interferograms
shows a much larger deficit than the typical one. Both are reported because the
gap is itself informative, and the jackknife is computed on the mean.

**On the p-value and uncertainty.** The 356 pairs share ~90 acquisition dates and are
therefore not independent, so the Wilcoxon *p* is quoted as a nominal value
and is not the primary inferential basis of this result. The evidence is the effect
size and its temporal stability: a mean deficit of −0.081, negative in 89 % of pairs, with
an empirical date-jackknife leave-one-out range of [−0.0842, −0.0774] demonstrating strict
sign invariance when any single acquisition is removed, and a formal jackknife standard-error
95 % confidence interval of [−0.109, −0.052] ($SE = 0.01454$). Benchmarked externally against
global Sentinel-1 seasonal coherence distributions across herbaceous wetlands (Kellndorfer et al., 2022),
Zone A's mean coherence (0.408) falls into the lower quartile, reflecting persistent decorrelation.

The deficit therefore depends on no single acquisition and survives control for
baseline, atmosphere (both by pairing), slope (DEM) and canopy optical wetness
(Sentinel-2 features). Figure 7 and Figure S4 show the paired differences and the
decay curves (with freeze response illustrated in Figure S6).


![**Figure 7.** Matched-cover paired test. (a) Distribution of coh(A) − coh(C) per interferogram; (b) scatter of A against C, points below the 1:1 line indicating lower mat coherence.](figures/F07_paired_test.png)

![**Figure S4.** Coherence decay with temporal baseline, with fitted decorrelation times per zone.](figures/S04_coherence_decay.png)

#### 4.2.2 The cover matching is effective

A and C are phenological twins: median Sentinel-2 wetness −0.513 vs −0.522,
same WorldCover class, matched greenness and seasonal amplitude. The coherence
difference is consequently not a vegetation artefact; it arises from a
non-optical surface property at the radar scale.

#### 4.2.3 A spatially delimited unit

- **Radial profile** (Fig. 9): a low, flat plateau (≈ 0.40) throughout the
  interior, a sharp step at the boundary, and a peak just outside (≈ 0.47).
  The same discontinuity appears in σ⁰ and RVI.
- **Five independent sensors** express the polygon: coherence, σ⁰ VV, RVI,
  Sentinel-2 wetness, temporal coherence (Fig. 8, Fig. S1).
- **Area**: A + B = 90.24 ha vs 89.7 ha documented (+0.6 %).

This is not diffuse noise but a delimited physical unit, whose outline —
drawn from vector and optical sources — coincides with structure visible in
independent radar fields.


![**Figure 8.** Per-zone distributions across five independent sensors: coherence, σ⁰ VV, RVI, Sentinel-2 wetness and temporal coherence.](figures/F08_zone_distributions.png)

![**Figure 9.** Radial profiles of coherence, σ⁰ VV and RVI against signed distance to the peatland boundary. The step at distance zero marks a physical edge.](figures/F09_radial_profiles.png)

#### 4.2.4 Scattering signature

**Table 5**:

| Quantity | A (mat) | C (grassland) | Reading |
|---|---|---|---|
| Median σ⁰ VV | −10.09 dB | −11.22 dB | A is brighter (+1.1 dB) |
| RVI (dual-pol) | 0.914 | 0.881 | A is more depolarising |
| Amplitude dispersion D_A | 0.243 (59 % < 0.25) | 0.238 (69 %) | A is not radar-dark |
| Closure dispersion (median \|closure\|) | 0.683 rad | 0.212 rad | ×3.2 |

At C-band the mat behaves as a denser, wetter scattering volume than dry
grassland despite identical optical phenology. The higher RVI excludes an
open-water double-bounce mechanism. The 3.2-fold closure dispersion is a direct
measurement of scatterer non-stationarity: mat triplets do not close, stable
ground triplets do (Fig. 13, Fig. S5).

Importantly, A is not devoid of targets (59 % of pixels have D_A < 0.25).
Its problem is not absent backscatter but an unstable phase — which is what
justified attempting phase linking in the first place (§4.1).

#### 4.2.5 Environmental predictors operate in the mat and vanish in grassland

Within-zone analysis (internal variability, not A vs C). Over the 499 mat pixels
with 12 covariates: cross-validated $R^2_{\text{cv}} = 0.239 \pm 0.022$ (random forest 0.326),
against 0.127 without the radar covariates.

**Table 6**:

| Covariate | ρ in A (mat) | ρ in C (grassland) | Difference |
|---|---|---|---|
| σ⁰ VV | −0.379 | −0.008 | −0.371 |
| Mean greenness | +0.320 | −0.009 | +0.329 |
| Elevation | −0.168 | +0.430 | −0.598 |

The radar and optical variables show active environmental sensitivity in the mat
that is absent in the matched grassland (Fig. 10b). In Zone A, higher backscatter (σ⁰ VV)
is associated with reduced coherence (ρ = −0.379), and optical greenness correlates
positively (ρ = +0.320). In Zone C, both coefficients are indistinguishable from zero
(ρ = −0.008 and −0.009). The apparent elevation contrast (−0.168 vs +0.430) is uninterpretable
given the 30 m Copernicus DEM noise floor across flat terrain. This pattern reflects the
presence of active dielectric and scattering modulation inside the floating mat against
its absence in mineral grassland, rather than a genuine sign reversal.

**Spatial autocorrelation and degrees of freedom.** Spatial autocorrelation in Zone A
exhibits an empirical correlation length of ~160 m (4 pixels on the 40 m grid; §3.6),
yielding $N_{\text{eff}} \approx 31$ independent spatial degrees of freedom across 499 pixels.
With 12 covariates in the ridge model, the observations-per-parameter ratio is under 3
($31 / 12 \approx 2.6$). Cross-validated $R^2_{\text{cv}} = 0.239$ is therefore reported
strictly as an empirical descriptive benchmark of within-zone spatial structure rather than
an independent predictive model.

**Reading precautions.** (i) The RVI / VH-VV collinearity (VIF ≈ 240; monotone
transforms of the same ratio) produced spurious coefficients (−1.21 and +0.95);
only the cleaned model is interpretable (σ⁰ −0.275, RVI −0.251, wetness
−0.240). (ii) Greenness is a proxy: ρ = +0.320 marginally but a partial
coefficient of +0.029 — it predicts nothing once σ⁰ and wetness are accounted
for. (iii) Elevation is likely a positional proxy: on a near-flat floating
mat the 30 m DEM relief is at noise level, and it should not be read physically.
(iv) The only predictor robust across both linear and non-linear models is
σ⁰ VV.


![**Figure 10.** Within-mat predictive model. (a) Standardised coefficients of the collinearity-cleaned model; (b) Spearman correlations in mat versus grassland, showing sign reversal.](figures/F10_predictors.png)

![**Figure S5.** Amplitude dispersion D_A: map and per-zone distributions against the 0.25 persistent-scatterer threshold.](figures/S05_amplitude_dispersion.png)

![**Figure S6.** Coherence sensitivity to the water-table proxy, and coherence gain on freezing, by zone.](figures/S06_hydrology_freeze.png)

#### 4.2.6 Verdict: H2 supported

> **H2 is supported.** At matched cover and after controlling baseline,
> atmosphere, slope and optical wetness, the mat shows significantly lower
> coherence, a sharp boundary, a volumetric scattering signature and 3.2-fold
> non-stationarity — and environmental predictors modulating its coherence have
> no detectable effect in grassland. It is a distinct radar unit, not "vegetation at C-band".

---

### 4.3 H3 — Dielectric signal, not motion

#### 4.3.1 The change of observable works

On identical simulated data (Fig. S3), per-pixel inversion returns
**−13.7 mm yr⁻¹** with 36 % usable pixels, whereas aggregation returns
**−19.8 mm yr⁻¹** against a ground truth of **−20**. The signal was not below the
noise floor; it was below the **per-pixel** noise floor.

#### 4.3.2 Common phase across zones (|R|)

| Zone | Median \|R\| |
|---|---|
| C — grassland | **0.569** |
| D — other | 0.426 |
| B — lake | 0.395 |
| **A — mat** | **0.234** |

The mat has the **lowest common phase of all zones**, consistent with H2. **This
ranking is the only claim this test supports**, for three reasons: |R| contains
atmosphere, common to all pixels; the ratio to the 1/√N_eff floor is not
comparable across zones; and with **measured** N_eff (§4.5.1) A and C sit only
×1.3 above their floors, which is not a detection. This test is therefore
**motivation and ordering, not proof**.

#### 4.3.3 Velocity has no power on a periodic signal

| Network | Signal A − C | Null (floor) |
|---|---|---|
| Full | −1.53 mm yr⁻¹ | **−1.50 mm yr⁻¹** |
| Baselines ≤ 60 d | −3.87 mm yr⁻¹ | **−4.88 mm yr⁻¹** |

The signal is **indistinguishable from the null**. More fundamentally, bog
breathing is **seasonal**, and a *velocity* regressed on a periodic signal is
not ≈ 0 but an artifact of where the observation window falls relative to the
cycle: bounded in magnitude by 6*A*/π*N*² for amplitude *A* over *N* annual
periods (§3.3), it changes sign when the window shifts by half a cycle. With
the amplitude of Table 7 over three annual cycles that bound is well below the
rate tabulated above, so it does not by itself explain the value — but it does
mean the test had **no power** on the physics being sought, not because it
returns zero, but because what it returns is set by the calendar rather than by
the peatland. The null control is what exposed the design error.

#### 4.3.4 Seasonal amplitude: detection

Annual-cycle fit on the aggregated A − C series against approximately 4 600
reference-matched null realisations (Fig. 11, Fig. 12; Table 7):

**Table 7** — Seasonal amplitudes and permutation significance against empirical null realisations:

| Series | Amplitude (mm) | Phase (DOY) | Seasonal $R^2$ | $p_{\text{perm}}$ | $N_{\text{null}}$ | Null type |
|---|---|---|---|---|---|---|
| A − C (mat vs grassland) | 3.29 | 104 | 0.30 | 0.026 | 4 614 | reference-matched |
| B − C (lake vs grassland) | 2.63 | 95 | 0.11 | 0.136 | 249 | reference-matched |
| A − B (mat vs lake) | 0.90 | 146 | 0.05 | 0.448 | 1 000 | size-matched |
| NULL (grassland vs stable) | 0.57 | 95 | 0.06 | — | — | — |

The reference-matched spatial null distribution has a median of 1.69 mm and a 95th percentile
of 2.93 mm; 119 of 4 614 null realisations exceed the observed value. This yields an empirical
permutation $p = (119 + 1)/(4 614 + 1) = 0.026$ (raw frequency ratio $119 / 4 614 = 0.0258$,
with a Monte Carlo 95 % confidence interval on $p$ of [0.021, 0.031]). The point estimate of
the seasonal amplitude is 3.29 mm LOS (3.89 mm vertical equivalent), with a formal 95 % upper
bound of 7.32 mm LOS (8.66 mm vertical). Size-matched compact nulls yield empirical $p$-values
of 0.014 and 0.022, confirming sensitivity robustness against the fragmented grassland control.
The mid-April maximum (DOY 104) aligns with spring water-table peaks.


![**Figure 11.** Aggregated series. (a) A−C, B−C, A−B and the size-matched null; (b) seasonal amplitudes — the lake oscillates like the mat, and A−B cancels.](figures/F11_aggregate_series.png)

![**Figure 12.** Significance of the seasonal amplitude against size-matched null distributions. (a) Full network; (b) winter pairs excluded.](figures/F12_significance.png)

#### 4.3.5 Three independent arguments exclude motion

**(a) Lake seasonal amplitude and consistent trajectory.** The residual open-water lake cannot
breathe mechanically, yet exhibits an annual trajectory consistent in amplitude and phase with
the floating mat: 2.63 mm LOS, phase DOY 95 (*p* = 0.136 against the reference-matched null).
The lake signal represents 80 % of the mat amplitude, within 9 days of the same phase.
Because $p = 0.136$ falls short of confirmatory statistical significance under our pre-specified
protocol, the lake trajectory cannot be claimed as an independent detection. However, its
trajectory provides a consistent amplitude scale.

Two physical hypotheses could account for coherent phase over Zone B: (1) emergent and submerged
macrophytes along the lake margins acting as distributed scatterers modulated by water-level
and dielectric shifts, or (2) spatial leakage from the adaptive Goldstein phase filter
($\alpha = 0.5$) smoothing adjacent mat phases across the narrow 65-pixel lake geometry.
Under either hypothesis, the lake trajectory is consistent with an environmental or dielectric
scaling rather than differential mechanical breathing.

**(b) Mat minus lake cancels.** Referencing A to the lake rather than the
grassland gives 0.90 mm, phase DOY 146 (random), seasonal R² 0.05,
*p* = 0.448 — sitting squarely within the empirical null distribution (null median 0.83 mm,
baseline NULL amplitude 0.57 mm). Mat and lake are seasonally indistinguishable.
Had the mat been breathing mechanically while the lake was not, A − B would have
revealed it. It reveals nothing.

**(c) Order of magnitude.**

| | Expected amplitude |
|---|---|
| Mat floating freely on a ±10 cm water table | ≈ 100 mm |
| Published raised-bog breathing (Hrysiewicz et al., 2024) | 10–40 mm |
| **Measured here** | **3.3 mm** |

Read against the measured amplitude alone, the signal is far below free
flotation. That comparison does not survive its own uncertainty: propagating
the 95 % interval (§4.3.7) to vertical gives 8.7 mm, and dividing by any coupling
fraction below unity raises it further — 17.3 mm at *f* = 0.5, 34.6 mm at
*f* = 0.25 — which overlaps the 10–40 mm reported for raised bogs. The
order-of-magnitude argument is therefore not available to us, and we do not
use it. What the comparison supports is narrower: C-band InSAR on this surface
cannot distinguish an absence of motion from motion of peatland-breathing
amplitude.

#### 4.3.6 Closure-phase bias does not discriminate

Over the 518 closed triplets in the network (Fig. 13, Table 8):

| Zone | Mean bias (rad) | σ | Median \|closure\| |
|---|---|---|---|
| A | −0.090 | 1.6 | 0.683 |
| B | −0.088 | 1.5 | 0.777 |
| C | +0.027 | 1.1 | 0.212 |
| D | −0.021 | 1.1 | 0.210 |

No systematic bias is detected. (A pre-registered prediction that increasing
the triplet count would push this to ≈ 5σ was falsified: the network holds
only 518 closed triplets, and at 518 the estimate *decreased* — the behaviour of
a fluctuation.)

What is robust is the dispersion (×3.2; π/2 ≈ 1.57 would correspond to
purely random triplets, so A remains partially coherent). But high dispersion
without a sign bias indicates random scatterer reconfiguration, which
moisture fluctuation and non-rigid micro-movement produce equally. This test
measures the degree of non-stationarity, not its nature.


![**Figure 13.** Closure phase by zone. (a) Mean bias with 2σ error bars — none significant; (b) median |closure| dispersion against the π/2 random reference.](figures/F13_closure_phase.png)

#### 4.3.7 Upper bound on motion, with stated assumptions

**Level 1 — robust ceiling (no assumption about the lake).** The total A − C
seasonal amplitude is 3.29 mm LOS. Attributing all of it to motion — that
is, deliberately ignoring §4.3.5 — gives:

> $d_{\text{vert}} \le 3.29 / \cos(32.26^\circ) \approx$ **3.9 mm** on the point estimate, and
> $\le 7.32 / \cos(32.26^\circ) \approx$ **8.7 mm** on the upper 95 % interval — which is the
> value we carry forward, since the point estimate alone understates it.

Assumptions: purely vertical motion; no phase aliasing (verified, since
centimetre-scale motion would produce an incoherent aggregate rather than a
clean annual cycle at R² = 0.30). This is the figure to quote by default: it is
independent of the lake, sitting ≈ 25× below free flotation against the point
estimate (3.9 mm) and ~11× below free flotation against the carried-forward
8.7 mm bound.

**Level 2 — refined bound (assumes a stable lake).** The mat-minus-lake residual
of 0.90 mm lies below the matched-null p95 of 2.0 mm:

> mat-specific motion < 2 mm LOS (≈ 2.4 mm vertical).

> **Critical assumption, stated.** This test assumes the lake's scattering
> surface is mechanically stable. Lake and mat float on the same water table;
> if both rose and fell together, A − B would cancel even in the presence of
> substantial motion. The observed cancellation is compatible with two readings
> — no motion, or common motion — which radar data alone cannot separate. Level 1
> already excludes flotation-scale motion independently of the lake, and in-situ
> laser measurement will resolve the ambiguity.

#### 4.3.8 Verdict: H3 rejected

> **H3 is rejected.** The detected seasonal signal (3.29 mm, *p* = 0.026) is
> dielectric: the lake, which cannot breathe mechanically, exhibits a consistent
> amplitude and phase trajectory; the mat-minus-lake difference cancels (0.90 mm,
> *p* = 0.45). We are measuring a seasonal moisture contrast between saturated
> surfaces and dry grassland. The magnitude of the signal does not independently
> exclude flotation once its uncertainty and the phase-centre coupling are
> propagated (§4.3.5c).

*Distinction to maintain*: this establishes that the **seasonal signal** is
dielectric. It says nothing about the nature of the **decorrelation** mechanism,
which remains undetermined between dielectric variability and non-rigid
micro-movement (§5.5).

---

### 4.4 H4 — What the sensor tracks is surface wetness

#### 4.4.1 The seasonality trap

The raw lag sweep gives Sentinel-2 wetness *r* = +0.576 (lag 54 d, *p* ≤ 0.021)
**and** air temperature *r* = −0.509 (lag 72 d, *p* ≤ 0.021), both apparently
significant. This is **unusable as it stands**: temperature and wetness are
strongly anti-correlated seasonally; the antecedent precipitation index — the
most *direct* hydrological proxy — fails (*p* = 0.43); and two annual-cycle
signals always correlate at *some* lag, the sweep merely aligning phases.

#### 4.4.2 Anomaly analysis

Removing the annual harmonic from both series leaves only inter-annual and
event-scale anomalies. Against 92 size-matched nulls (Fig. 15, Table 9):

| Forcing | *r* seasonal | *r* ANOMALIES | Lag | *p* |
|---|---|---|---|---|
| NDWI zone A | 0.576 | +0.450 | 12 d | ≤ 0.011 |
| NDWI zone C | 0.490 | +0.427 | 12 d | 0.022 |
| NDWI zone D | 0.519 | +0.424 | 42 d | ≤ 0.011 |
| NDWI(A) − NDWI(C) | 0.395 | −0.316 | 6 d | 0.150 |
| Antecedent precipitation | 0.230 | 0.293 | 6 d | 0.172 |
| Precipitation | 0.191 | −0.225 | 66 d | 0.312 |
| Air temperature | −0.509 | 0.224 | 78 d | 0.581 |

To prevent selection bias from lag sweeping, the lag-0 correlation is treated as
the primary effect size ($r \approx 0.42$–$0.45$ across zones, $p \le 0.022$),
with the 12-day swept maximum ($r = +0.450$) reported as secondary exploratory evidence.

#### 4.4.3 Three convergent facts

**(a) Temperature collapses.** −0.509 → 0.224, *p* from 0.021 to 0.581. Its
correlation was only the shared annual cycle. The temperature–wetness
confound is thereby resolved, and a thermal artefact is excluded.

**(b) Wetness survives**, and it originates from an optical sensor entirely
independent of the radar (different platform, different measurement physics).

**(c) The lag drops from 54 d to 12 d** — one revisit cycle, hence
instantaneous at our sampling resolution. This was the criterion set a
priori: a dielectric response is near-instantaneous, while mechanical settling
lags by weeks. The lag therefore confirms the dielectric mechanism through a
route independent of the lake control.

The sign is consistent: wetter → shallower penetration → phase centre higher
→ apparent uplift (positive *r*). Sign alone does not discriminate, since
mechanical swelling would give the same, but it is coherent.


![**Figure 14.** Aggregated InSAR phase and Sentinel-2 optical wetness — two independent sensors.](figures/F14_phase_vs_wetness.png)

![**Figure 15.** Correlation with each forcing, raw versus deseasonalised: temperature collapses while wetness survives.](figures/F15_seasonal_vs_anomalies.png)

#### 4.4.4 The mechanism: a sensitivity contrast

The failure of the differential forcing (NDWI_A − NDWI_C: *r* = −0.316,
*p* = 0.15) is initially surprising, since the InSAR series is itself
differential.

**Model.** A common regional moisture M(t) drives φ(A) with sensitivity k_A
and φ(C) with k_C, where k_A > k_C (saturated peat responds more strongly
than mineral grassland). Then

> φ(A) − φ(C) = (k_A − k_C) · M(t)

The differential phase therefore tracks absolute moisture through a
sensitivity contrast rather than a moisture contrast; and
NDWI(A) − NDWI(C) ≈ 0 + noise, since both surfaces respond optically in a
similar way.

**Pre-registered prediction.** If the model holds, NDWI over C and over D
— proxies for the same M(t) — must also correlate positively, with
comparable magnitude.

**Confirmed**: NDWI(A) +0.450, NDWI(C) +0.427, NDWI(D) +0.424 — nearly identical
— while the differential fails. The model is therefore validated by a
prediction stated before the test, not by post-hoc rationalisation, and the
failure of the differential forcing is a consequence of the model rather than
evidence against the hydrological link.

#### 4.4.5 Verdict: H4 supported

> **H4 is supported, with a moderate effect size.** On anomalies, the aggregated
> phase co-varies with surface wetness (*r* = 0.42–0.45 depending on reference
> zone, *p* ≤ 0.022) at near-zero lag, while temperature does not survive
> deseasonalisation. Coupling operates through a sensitivity contrast between
> saturated peat and mineral grassland.

**Limits.** The effect is moderate (0.450 against a null p95 of 0.404): this is
a measurable sensitivity, not an operational hydrological product. The three
zonal NDWI series are not independent (all proxy the same regional M(t)), so
this is one confirmation of a predicted pattern rather than three. *p*-values sit
at the 1/(1 + 92) floor. Optical canopy versus soil dielectric caveat: Sentinel-2 NDWI
reflects surface canopy wetness rather than the deep peat dielectric profile. While
surface moisture correlates with peat hydrology, microwave penetration interacts with
the moss layer, meaning optical wetness serves as an external proxy rather than a direct
probe of peat permittivity.

---

### 4.5 Robustness and falsification

Five of eight alternative explanations are strictly excluded, one is largely excluded,
one is quantitatively measured and accounted for, and one remains open awaiting in-situ
validation (Appendix A). Two results warrant reporting here because they modified our own conclusions.

#### 4.5.1 Effective sample size, measured

The 1/√N argument assumes independent pixels, which they are not. Measuring the
spatial correlation length by empirical autocorrelation (1/e threshold):

| Zone | L_corr | N_eff **measured** | N_eff *assumed* | \|R\| / measured floor |
|---|---|---|---|---|
| A | 160 m | 31 | *125* | **×1.3** |
| B | 80 m | 16 | *16* | ×1.6 |
| C | 360 m | 5 | *100* | **×1.3** |
| D | 280 m | 219 | *2 688* | ×6.3 |

**Consequence, accepted:** the |R| comparison of §4.3.2 loses most of its power
and is demoted to motivation and ordering. **The principal conclusions are unaffected**,
because their significance derives from empirical size-matched nulls that use no
N_eff at all. (Caveat: zone C is fragmented, so the estimator mixes within- and
between-patch correlation and its N_eff of 5 is probably understated.)

#### 4.5.2 Snow and frost excluded

Snow has its own annual cycle and affects saturated peat differently from
grassland, making it a complete competing explanation. Removing all
December–February pairs (30 % of the network):

| Dataset | *n* pairs | Amplitude | Phase (DOY) | Seasonal R² | *p* (size-matched)* |
|---|---|---|---|---|---|
| Full | 356 | **3.286 mm** | 104.2 | 0.299 | 0.014* |
| **Winter excluded** | 248 | **3.282 mm** | 112.6 | 0.309 | 0.022* |

\*Evaluated under the size-matched null; the reference-matched full-network value is *p* = 0.026.
The amplitude changes by **0.1 %** and the seasonal R² slightly *increases*.
Snow and frost are **refuted**; the signal is carried entirely by the growing
season.
