## 5. Discussion

### 5.1 Why C-band succeeds on raised bogs and fails here

Our results do not contradict Hrysiewicz et al. (2024), who retrieved bog
breathing over Irish raised bogs with correlations of 0.8–0.9, nor Alshammari et al. (2018)
or Tampuu et al. (2023), who demonstrated InSAR surface motion tracking across Scottish
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
In drained or degraded peatlands undergoing consolidation, oxidation, and compaction,
InSAR readily tracks secular subsidence on the order of centimetres per year
(e.g. Hoyt et al., 2020; Patil et al., 2026).
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

Given our aggregate velocity detection floor of $\approx 1.5\text{--}5\text{ mm yr}^{-1}$ (§4.3.3),
our data exclude multi-year linear subsidence exceeding $\approx 5\text{ mm yr}^{-1}$ over Rzecin,
but cannot distinguish slower subsidence from zero motion. The comparison is therefore informative
not as a demonstration of 'lower subsidence', but as an empirical boundary delimiting where C-band
peatland interferometry succeeds: drained peat (agricultural surfaces, high coherence, centimetre-scale signal)
rather than near-permanently saturated floating mats (wet canopy, low coherence, millimetre-scale signal).

**A methodological caution our results support.** Patil et al. interpret
seasonal fluctuations aligned with soil moisture as hydrological control of peat
surface motion. Our results counsel care with that inference: a seasonal
oscillation correlated with moisture is not automatically motion. At our site,
a 3.3 mm signal correlating with surface wetness can be reproduced entirely by
seasonal dielectric permittivity variations in the moss layer (forward modeling in §5.2),
demonstrating that moisture variations induce substantial apparent LOS motion without
mechanical displacement. Phase variations induced by dielectric permittivity
and moisture changes are well established in the InSAR literature (De Zan et al.,
2014, 2015; Morrison et al., 2011; Nolan & Fatland, 2003; Nolan et al., 2003;
Rabus et al., 2010; Ranjbar et al., 2021; Karamvasis & Karathanassi, 2023;
Zwieback et al., 2015, 2017; Mira et al., 2022; Zheng & Fattahi, 2025).

A control over a water surface, or any target where mechanical breathing is physically
precluded, is inexpensive and separates genuine displacement from differential
propagation phase or dielectric permittivity effects (De Zan et al., 2014). We
suggest incorporating one systematically in peatland motion studies reporting
millimetre-scale signals.

**Forward-model consistency analysis.** To evaluate whether the observed seasonal phase
oscillation can be physically accounted for by dielectric permittivity changes without
requiring mechanical peat displacement, we extend the lossy dielectric half-space framework of
De Zan et al. (2014) to saturated organic peat (cf. Nolan & Fatland, 2003; Rabus et al., 2010;
Ranjbar et al., 2021). Rather than presenting this as an exact calibration, we treat it as an
order-of-magnitude consistency analysis.

*Material parameters and mixing model.* Peat soils are structurally distinct from mineral soils:
they comprise approximately 90–95% pore space with very low solid bulk density. Standard empirical
polynomial models calibrated on mineral sand and clay fractions (e.g., Hallikainen et al., 1985)
carry no organic term and are inapplicable. We therefore implement a complex refractive mixing model
(Birchak et al., 1974; Ulaby & Long, 2014) with exponent $\alpha = 0.5$:
$$\epsilon_{\text{eff}}^{1/2} = m_v \epsilon_w^{1/2} + v_s \epsilon_s^{1/2} + (1 - v_s - m_v) \epsilon_{\text{air}}^{1/2}$$
where $\epsilon_{\text{air}} = 1.0$, dry organic solids have permittivity $\epsilon_s \approx 2.2$
with volume fraction $v_s = 0.07$ (consistent with peat bulk density $\rho_b \approx 0.10\text{ g cm}^{-3}$,
particle density $\rho_s \approx 1.45\text{ g cm}^{-3}$, and porosity $\phi \approx 93\%$), and the
complex permittivity of free water $\epsilon_w$ is evaluated via Debye relaxation (Klein & Swift, 1977;
Stogryn, 1971; Ulaby & Long, 2014) at $f = 5.405\text{ GHz}$ and mean in-situ growing season
temperature $T = 15^\circ\text{C}$ ($\epsilon_w = 72.9 - 24.9j$). This yields an effective loss tangent
$\tan \delta = \epsilon'' / \epsilon' \approx 0.25 - 0.35$.

*Geometry and layer structure.* For the Sentinel-1 radar geometry ($f = 5.405\text{ GHz}$,
$\lambda = 5.5466\text{ cm}$, measured local incidence angle $\theta = 32.26^\circ$ from HyP3 burst
metadata, VV polarization), the complex vertical wavenumber is $k_z(\epsilon) = \sqrt{k_0^2 \epsilon_{\text{eff}} - (k_0 \sin \theta)^2}$.
The resulting power penetration depth ($1/e$) is $\delta_p = 1 / (2 |\text{Im}(k_z)|) \approx 3.4 - 4.4\text{ mm}$
across saturated conditions ($m_v \in [0.70, 0.90]$; $3.6\text{ mm}$ at $m_v = 0.85$). This shallow skin depth
justifies treating the living moss carpet as an effective lossy dielectric half-space for coherent propagation:
the tightly packed *Sphagnum* capitula form a continuous, near-saturated matrix whose element spacing
($\sim 1\text{--}5\text{ mm}$) is much smaller than the radar wavelength ($\lambda_0 = 55.5\text{ mm}$),
fixing the coherent reflection boundary within the topmost 3–4 mm. Concurrently, the sparse, emergent sedge
canopy (*Carex limosa*, *Eriophorum*) extending above this carpet lacks sufficient bulk permittivity to form
a distinct dielectric boundary, but its non-stationary geometric reconfiguration and branch-scale scattering
induce significant volume scattering and wave depolarisation (RVI = 0.914). Coherent propagation phase and
interferometric coherence therefore probe distinct physical processes of the same target: phase tracks the
dielectric moisture shift of the capitulum skin layer, whereas canopy volume scattering and microstructural
rearrangement drive the severe temporal decorrelation.

*Moisture excursion in a saturated fen.* A volumetric moisture change of $\Delta m_v \approx 0.25$
(from near-saturated $m_v = 0.85$ in spring down to $0.60$ in late summer) in the uppermost 3–4 mm
skin layer is physically realistic and does not contradict the "hydrologically stable" water table
(0–30 cm depth) described in §2.1. *Sphagnum* mosses are non-vascular plants lacking stomata and roots;
water is drawn to the capitula entirely via capillary forces through the external stem leaf wicks.
During warm, sunny summer periods of high atmospheric evaporative demand, evaporation from the exposed
capitulum tips readily exceeds the rate of capillary replenishment from below, creating transient moisture
depletions of $\Delta m_v \approx 0.20 - 0.30$ in the uppermost few millimetres of the living capitula
(McCarter & Price, 2014; Strack et al., 2009) even while the underlying peat column remains permanently
waterlogged.

*Predicted apparent displacement envelope.* Leading with the uncertainty envelope across plausible
near-surface moisture excursions $\Delta m_v \in [0.15, 0.35]$, the model predicts an apparent
interferometric line-of-sight displacement spanning $-1.9\text{ mm}$ to $-4.6\text{ mm}$ (incorporating
temperature sensitivity across $2^\circ\text{C}$ to $25^\circ\text{C}$; $-2.1\text{ mm}$ to $-4.2\text{ mm}$
at $15^\circ\text{C}$). The observed 3.29 mm point estimate (Table 7) falls comfortably within this predicted
envelope. Rather than claiming a fine-tuned parameter alignment, the defensible physical conclusion is that
the observed amplitude is indistinguishable from the dielectric prediction, which is all the non-mechanical
hypothesis requires.

Two critical properties accompany this forward-model consistency analysis:

1. *Asymptotic ceiling and model limitations (Table T16)*: The theoretical maximum apparent LOS
   displacement under complete desiccation ($m_v \to 0$) converges asymptotically to $6.13\text{ mm}$ LOS.
   While our observed semi-amplitude ($3.29\text{ mm}$) falls safely below this ceiling, the full
   peak-to-peak swing under an unconstrained linear sinusoidal fit ($2 \times 3.286 = 6.57\text{ mm}$) slightly
   exceeds it. This exceedance represents an honest limitation of applying an unstratified half-space model
   with a simple sinusoid to a non-linear physical system: true capillary drying saturates. Fitting an empirical
   saturating seasonal model (Table T16) contracts the peak-to-peak excursion to $5.64\text{--}5.83\text{ mm}$
   (semi-amplitudes of $2.82\text{--}2.91\text{ mm}$), below the $6.13\text{ mm}$ ceiling. *(Note on $R^2$ values)*:
   Table T11 reports $R^2 = 0.299$ for the harmonic fit evaluated directly on the pairwise baseline series (356 pairs),
   whereas Table T16 reports $R^2 = 0.364$ for the harmonic fit on the inverted cumulative time series (90 dates).
2. *Separation of phase and coherence mechanisms*: The dielectric propagation model is sufficient to explain
   a phase perturbation of the observed order of magnitude (predicting apparent LOS displacements of $-1.9$
   to $-4.6$ mm), but a propagation model alone does not account for the observed coherence loss; additional
   scattering and target-instability processes are therefore required. The propagation model evaluates the
   coherent phase perturbation induced across a lossy dielectric boundary, whereas the observed temporal
   coherence ($\bar{\gamma} = 0.408$, Table 4) carries the entire scattering budget—including volume
   decorrelation, non-stationary vegetation growth, thermal noise, and temporal decorrelation. This clarifies
   that while dielectric fluctuations provide a physically consistent explanation for the seasonal phase
   oscillation, they operate alongside distinct physical processes (primarily canopy volume scattering and
   structural reorganizations) that drive the severe coherence loss.


Furthermore, testing maximum temporal baseline subsets (24 d to 120 d and all-pair; Table 3, Table T15)
reveals that while short-baseline networks (≤ 24 d) suffer from severe accumulating
closure-phase subsidence bias (−13.5 to −23.5 mm yr⁻¹), consistent with Zheng et al. (2022) fading signal bias,
expanding the network to ≥ 48 d and annual pairs causes the apparent velocity to
contract to near-zero (−1.53 mm yr⁻¹ on A−C vs −1.50 mm yr⁻¹ on NULL), demonstrating that
our velocity non-detection is robust to network truncation. Notably, while linear velocity is
highly sensitive to this fading bias, the seasonal harmonic amplitude stabilizes reliably across
all intermediate and full network configurations: the seasonal amplitude settles between 2.89 mm
(at $\le 48$ d) and 3.29 mm (all pairs), while seasonal phase locks tightly between DOY 104 and 106.
Harmonic inversion on aggregated multi-temporal networks is thus protected from the short-baseline
fading bias that severely distorts linear deformation rates.

Our estimated displacement bound sits in stark contrast with both raised-bog breathing and drained-fen subsidence scales documented in the literature (Figure S14).

### 5.3 Two transferable methodological contributions

#### 5.3.1 Change of observable

Per-pixel phase noise (≈ 1.5 rad at γ = 0.4) falls as 1/√N_eff under
aggregation. Over 499 pixels it drops to ≈ 0.3 mm (≈ 1 mm at N_eff = 50), far
below the signal sought. **The signal was not below the noise floor; it was below
the *per-pixel* noise floor.**

This reasoning applies to any target that is **spatially coherent but temporally
decorrelated**: peatlands, rock glaciers, wetlands, crops. The condition is that
the target deform as a unit — an assumption that must be physically justified and
is testable by subdividing the zone. Here, subdividing Zone A into concentric distance bands
(inner core, deep core, and outer margin; §4.3.4, Table T11) revealed no differential behaviour
at the resolvable scale: all sub-zones exhibit synchronized phase locking (DOY 100–108) and
consistent seasonal amplitudes (2.89–3.78 mm), which is consistent with spatially coherent seasonal
moisture variation across the mat without asserting proven kinematic rigidity.

#### 5.3.2 Weak-signal test protocol

A protocol structured around four analytical safeguards that invalidated intermediate conclusions or corrected numerical bounds during this study (Appendix A.12), distilled into three operational rules:

1. **Size-matched nulls.** Aggregate noise falls as 1/√N, so a null four times
   larger carries half the noise and **manufactures false detections**.
2. **A null distribution, not a single realisation.** One realisation is not a
   test; *N* draws give an empirical *p*-value — whose **floor** of 1/(1 + N)
   must be stated. Across our 10,750-pixel reservoir (Zone D), spatial autocorrelation
   limits effective independent draws to $N_{\text{eff}} \approx 219$ (Table 10), which
   remains fully adequate for establishing the empirical null distribution (§3.4).
3. **Identical treatment of the null.** If the observed statistic results from a
   selection (best |r| over 16 lags), the null must undergo the same sweep.

These rules are cheap, fully specified algorithmically (§3.4), and should accompany
any weak-signal claim over decorrelated terrain.

### 5.4 Instrumental outlook

**L-band is the most direct route.** At λ ≈ 24 cm, radar signals penetrate
herbaceous vegetation more effectively, significantly mitigating the rapid
temporal decorrelation that affects C-band (λ = 5.5 cm) over dynamic, saturated
wetland covers (Chen et al., 2021; Morishita & Hanssen, 2015; Ranjbar et al., 2021).
**NISAR** provides global, open L-band observations. The conceptual analysis is
adaptable to NISAR GUNW products, although differences in product geometry, spatial
resolution, pairing strategy and processing require validation before the chain
transfers. This is a **prospective** test: the archive does not
cover 2022–2024 retrospectively.

**The return to a 6-day revisit will not suffice.** Sentinel-1C and 1D restore
two-satellite operation, reducing temporal baseline (ESA, 2024). However, all Sentinel-1
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

#### Multi-geometry replication protocol and 2-LOS decomposition

A fundamental geometric constraint of single-geometry InSAR is the one-dimensional line-of-sight (LOS) projection. The displacement bound derived in §4.3.7 assumes purely vertical displacement ($d_{\text{LOS}} = d_{\text{vert}} \cos \theta$). However, apparent LOS phase shifts can theoretically arise from horizontal displacement components or non-isotropic scattering perturbations.

To resolve vertical deformation ($d_{\text{vert}}$) from east-west horizontal motion ($d_{\text{east}}$) without unverified geometric assumptions, an independent descending-track acquisition geometry is required. For ascending and descending LOS measurements ($d_{\text{LOS}}^{\text{asc}}, d_{\text{LOS}}^{\text{desc}}$) with local incidence angles $\theta_{\text{asc}}, \theta_{\text{desc}}$ and satellite track heading angles $\alpha_{\text{asc}}, \alpha_{\text{desc}}$:

$$\begin{pmatrix} d_{\text{LOS}}^{\text{asc}} \\ d_{\text{LOS}}^{\text{desc}} \end{pmatrix} = \begin{pmatrix} \cos \theta_{\text{asc}} & -\sin \theta_{\text{asc}} \cos \alpha_{\text{asc}} \\ \cos \theta_{\text{desc}} & \sin \theta_{\text{desc}} \cos \alpha_{\text{desc}} \end{pmatrix} \begin{pmatrix} d_{\text{vert}} \\ d_{\text{east}} \end{pmatrix}$$

Solving this system uncouples vertical peat breathing from potential lateral drift. Furthermore, because a genuine seasonal mechanical deformation must reproduce in phase across orbits and scale rigorously with incidence-angle projection ($\Delta \phi \propto \cos \theta$), whereas processing artifacts, local multipath, and orbit-specific atmospheric screens do not, executing the identical phase-linking and aggregation pipeline over the Sentinel-1 descending track represents the most direct and cost-effective empirical replication test available.

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
to the water table — is not supported by any of three independent indicators (no
hydrological coupling of coherence, no stabilisation on freezing, no dominant
double-bounce signature). The closure-phase bias,
which was intended to separate (a) from (b), detects no systematic bias: high
dispersion **without a sign bias** is compatible with both.

This question is **distinct** from that of the seasonal signal, whose dominant
dielectric/propagation interpretation is the best-supported reading of the
available controls (§4.3).

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

### 5.6 Scope of inference: instrumental limit versus capability

Our results answer **two separate questions**, which must not be conflated (Table 12; Figure 7):

| | Question | Answer | Primary evidence |
|---|---|---|---|
| **Q1** | Can Sentinel-1 measure vertical displacement of the floating mat? | **No**; apparent phase-centre displacement $\le 8.7$ mm, true mat motion unconstrained | §4.1, §4.2, §4.3 |
| **Q2** | Can Sentinel-1 inform on seasonal wetland hydrological state? | **Possibly yes**, moderate but measurable association | §4.3.4, §4.4 |

Q1 is an **instrumental-limit** result: across all six tested standard burst-product inversion families, per-pixel phase retrieval collapses, and the residual aggregate signal cannot be uniquely attributed to surface motion because forward dielectric modeling demonstrates that upper-moss moisture fluctuations alone reproduce the 3.3 mm seasonal amplitude. Q2 is a **capability** result: spatial aggregation successfully recovers an annual phase cycle ($3.29$ mm, $p = 0.026$) that covaries with independent optical surface wetness anomalies at near-zero lag ($r \approx 0.45$).

Conflating Q1 and Q2 would weaken both findings. Placing all claims into the four-level observable hierarchy (Figure 7a) clarifies that failure to invert Level 1 physical peat deformation from Level 3 InSAR observables does not preclude using Level 3 observables to track Level 2 moisture-induced phase-centre shifts.


![**Figure 7.** Conceptual framework for wetland InSAR. (a) Four-level observable hierarchy distinguishing physical quantities, radar interactions, raw InSAR observables, and inferred quantities; (b) mechanistic causal chain at Rzecin, demonstrating how hydrologically driven phase-centre variability produces apparent displacement without requiring true mechanical peat breathing.](figures/F07_conceptual_framework.png)

### 5.7 Limitations

- **One site, one track, one polarisation** (VV): transferability of the
  predictive model (R²cv = 0.24) to other peatlands remains to be demonstrated.
- **No in-situ validation** (neither laser nor continuous WTD during 2022–2024):
  the constraint on apparent phase-centre displacement is internal to the InSAR
  analysis, and the coupling between that phase centre and the peat is unmeasured.
- **S1A-only window** (12-day revisit): a degraded cadence relative to what is
  now available.
- ***p*-values at the floor** (1/(1 + N)) for several tests; more null draws
  would tighten them.
- **Unit assumption** in aggregation: the mat is assumed to deform as a block;
  our concentric core-versus-margin subdivision test confirmed synchronous phase locking
  (DOY 100–108) and consistent amplitudes across sub-zones (§4.3.4, Table T11), though
  sub-kilometer flexural modes below our spatial resolution cannot be ruled out.
- **Zone C is fragmented**, which biases the empirical correlation-length
  estimator used for N_eff.
- **Absence of alternative water controls in the burst footprint**: An exhaustive search across
  the Sentinel-1 burst frame confirmed that no alternative isolated water body $\ge 3 \times L_{\text{corr}} = 480\text{ m}$
  from peat margins exists within the scene. Open water bodies in this morainic terrain are either
  narrow bog pools subject to boundary-filter leakage or drainage ditches below radar resolution.
  Stable non-vegetated surfaces (such as rock outcrops or bare mineral soil parcels) are likewise
  absent in this forested wetland basin.
- **Grassland reference zone fragmentation**: A burst-wide survey confirmed that non-mat mineral
  grasslands in this region are physically fragmented into small agricultural clearings interspersed
  among pine plantations. The effective sample size ($N_{\text{eff}} \approx 5$) of Zone C is thus an
  inherent constraint of the post-glacial landscape rather than a sampling artifact, justifying reliance
  on the regional multi-zone null reservoir (Zone D, $N_{\text{eff}} \approx 219$).

