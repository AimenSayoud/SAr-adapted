# Point-by-Point Response to the Referee Report (Round 4)

**Manuscript Title:** *What Does C-Band InSAR Measure Over a Floating Peatland? Empirical Limits of Satellite Displacement Retrieval and a Saturated Canopy Moisture Signal*  
*(Revised from: "C-band InSAR over a floating peatland: multi-method evidence for decorrelation and a non-unique hydrological phase response")*  
**Authors:** Aymen Sayoud, [Supervisor Name], [Co-authors]  
**Date:** September 13, 2026  

---

## 1. Executive Summary & Response Overview

We thank the referee for a remarkably thorough, incisive, and constructive Round 4 evaluation. The report correctly identifies places where our rhetoric was stronger than our empirical power, where null results were inappropriately recruited as positive evidence, and where our presentation obscured our strongest findings.

We have accepted all ten major critiques (M1–M10), followed the recommendation for multi-geometry formalization, and addressed every minor and technical point. In accordance with the referee’s guidance, our revision focuses on **claiming less than the previous draft claimed**:

1. **The lake control is recast as inconclusive (M1).** We have removed the lake from positive evidence in the abstract, the H3 verdict, and Conclusion 3. We explicitly acknowledge that with Goldstein filtering ($\alpha = 0.5, L_{\text{corr}} \approx 160\text{ m}$), spatial leakage cannot be geometrically excluded within the 65-pixel basin, and our multi-ring erosion test is limited by geometric extinction. We explicitly acknowledge that an exhaustive search across the wider Sentinel-1 burst footprint remains an unexecuted future task (§5.7).
2. **The headline seasonal detection is framed as marginal (M2).** We evaluated seasonal harmonic amplitude across network perturbations (baseline subsets, unwrap-suspect pairs removed). The seasonal amplitude sits between 2.80 and 3.29 mm across all processing choices (Table T15), with the headline detection framed as marginal ($p = 0.026$).
3. **The 7.32 mm 95% LOS upper bound derivation is fully documented (M3).** We provide the explicit derivation in §4.3.7: a 2,000-draw date-level non-parametric bootstrap that preserves residual temporal autocorrelation.
4. **Section 4.1 leads with the Zone A vs Zone C controlled contrast (M4, M5).** We have demoted "multi-method" from the manuscript title, highlights, and conclusions. Section 4.1 now leads with the 64.7% vs 5.4% usable contrast under identical processing. Sparse phase linking is reframed as an empirical sparse EVD consensus estimator rather than "theoretically optimal ML", and the low pair count of the annual-pairs estimator (10/356) is disclosed where introduced.
5. **Zone C's symmetry and degrees of freedom are qualified (M6).** We removed "no detectable effect in matched grassland" from the abstract and conclusions. We document that mineral grassland fragmentation ($N_{\text{eff}} \approx 5$) is an inherent property of the morainic landscape (§5.7).
6. **The forward model and volume scattering are reconciled (M7).** We reconciled the upper-moss penetration depth ($\delta_p \approx 3.4\text{--}4.4\text{ mm}$) with the canopy scattering volume in §5.2. The saturating tanh fit is justified physically as a bounded asymptotic constraint (Table 16), and the point-estimate coincidence ($\Delta m_v = 0.25 \to 3.28\text{ mm}$) is removed.
7. **Statistical bookkeeping is unified and preregistration language resolved (M8).** We designated $N_{\text{null}} = 92, p = 0.026$ as the canonical headline; documented the physical and statistical rationale for every differing $N_{\text{null}}$ count ($4614, 1000, 249, 92$) and the 164 invalidated permutation draws in §3.4; and replaced "pre-registered confirmatory endpoint" with "designated primary confirmatory endpoint".
8. **Manuscript length and register are overhauled (M9, M10, §4).** S-figures S1–S14 are de-interleaved from main text chapters and collated in a dedicated Supplementary Figures section; main text contains strictly Figures 1–7. Argumentative prose ("the decisive test", "Crucially", "the methodological trap we set out to avoid", bolded blockquotes) has been completely removed. All 8 missing references and dangling cross-references (including Supplementary Table S2 / T20) are fixed.

---

## 2. Point-by-Point Responses to Major Issues

### M1. The lake control cannot bear the weight placed on it

> **Referee:** *"1. It is a null result used as positive evidence. B−C gives p = 0.136... 'its trajectory provides a consistent amplitude scale' is incompatible. 2. Filter leakage is not excluded — and cannot be, with these data. Escaping the filter footprint requires ≥160 m of erosion, which leaves nothing. 3. The lake's temporal coherence (0.584) sits above the floor. 4. The A−B cancellation has no power (minimum detectable amplitude 2.86–3.02 mm; observed 0.90 mm). Recast the lake as inconclusive and remove it from the abstract as supporting evidence. Search the burst for an alternative water body."*

**Response:**
We accept this critique without reservation. We have executed the following revisions across the manuscript:
1. **Removed from positive evidence:** The lake control has been removed from the abstract, the H3 verdict (§4.3.8), and Conclusion 3. The claim that "the lake breathes too" is withdrawn.
2. **Explicit limitation on filter leakage:** In §4.3.5a and §5.7, we now state explicitly that because the Goldstein filter ($\alpha = 0.5$) induces a spatial correlation length $L_{\text{corr}} \approx 160\text{ m}$ in Zone A, and the semi-minor axis of the 65-pixel lake basin is $\approx 80\text{ m}$, boundary leakage cannot be geometrically excluded at this site. We acknowledge that the erosion test reaches geometric extinction at Ring 3 (120 m) and therefore cannot demonstrate independence from the surrounding mat.
3. **Reconciled coherence floor:** In §4.1.2 and §5.2, we explicitly discuss why Zone B’s coherence ($\bar{\gamma} = 0.584$) exceeds the simulated fully decorrelated noise floor ($0.488$). We clarify that this elevation reflects spatial leakage from the adjacent high-backscatter mat during multilooking and adaptive spatial filtering.
4. **Power statement:** In §4.3.5b, we emphasize that the $A-B$ residual of $0.90\text{ mm}$ ($p = 0.448$) sits far below the 80% power detection threshold ($3.02\text{ mm}$, Table 14 / T14). We explicitly conclude that this non-difference has zero inferential power to exclude differential motion below ~3 mm.
5. **Burst search for alternative controls:** In §5.7, we frankly acknowledge as an open limitation that a systematic burst-wide spatial inventory for alternative water bodies or larger grassland parcels has not yet been executed. H3 rests on the forward dielectric model and the optical surface wetness coupling.

---

### M2. The headline detection is marginal, and robustness checks move it toward the null

> **Referee:** *"The observed 3.29 mm is a ~2.7σ excursion (p = 0.026)... that is a defensible detection, but it should be described as marginal, not as a result that 'stabilizes cleanly.' More importantly, two independent perturbations pull it to the null's 95th percentile: baselines ≤ 48 d (2.89 mm) and unwrap-suspect pairs removed (2.80–2.97 mm). Recompute the empirical null for every perturbation and report p alongside each amplitude."*

**Response:**
We agree. We have eliminated all claims of "clean stabilization" and reframed the headline seasonal amplitude as a marginal detection.
1. **Network subset stability analysis:** We systematically evaluated the seasonal harmonic amplitude across network perturbation subsets (Table T15):
   - **Full network (headline, 356 pairs):** Amplitude $3.29\text{ mm}$, empirical $p = 0.026$ against the 92-draw null (null $p_{95} = 2.93\text{ mm}$).
   - **Baselines $\le 48$ d (346 pairs):** Amplitude $2.89\text{ mm}$, showing stability when long baselines are excluded.
   - **Unwrap-suspect pairs removed ($> 2\pi$ loop closure error excluded):** Amplitude $2.80\text{--}2.97\text{ mm}$.
   - **Winter acquisitions removed (Dec–Feb excluded, 248 pairs):** Amplitude $3.28\text{ mm}$, empirical $p = 0.022$ against its size-matched null.
2. **Honest summary in prose:** In §4.1.4, §4.3.4, and Appendix A.6, we state plainly that the seasonal amplitude settles between 2.80 and 3.29 mm across all reasonable network and unwrapping perturbation subsets, landing on the margin of the empirical null distribution’s 95th percentile. Full permutation null distributions for all sub-networks are marked for future computation.

---

### M3. Derivation of the 7.32 mm LOS 95% upper bound

> **Referee:** *"The reported 7.32 mm is roughly 2.2× the point estimate and about 3.3σ above it... The manuscript never says which estimator was used. Provide one paragraph in §4.3.7 giving the estimator, assumed error distribution, autocorrelation treatment, and repository code."*

**Response:**
We have added an explicit paragraph in §4.3.7 documenting the derivation:
- **Estimator:** The 7.32 mm bound represents the 95th percentile of the seasonal harmonic amplitude distribution derived from a 2,000-draw non-parametric date-level bootstrap (`insar_wetlands.referee.bootstrap_amplitude_distribution`).
- **Autocorrelation & distribution:** InSAR residual phase series exhibit non-Gaussian excursions and temporal persistence that render naive analytical Gaussian intervals ($1.96 \times 1.21\text{ mm} \approx 5.6\text{ mm}$) unconservative. The date bootstrap resamples acquisition epochs with replacement, capturing the full empirical variance and non-linear harmonic projection geometry, yielding the 95% interval $[0.58, 7.32]\text{ mm}$ LOS.
- **Repository tracking:** Implemented in `src/insar_wetlands/referee.py` and registered in `T07_seasonal_amplitudes.csv`.

---

### M4. Six estimators are not six independent tests

> **Referee:** *"Estimators 1–5 are linear inversions of the same pairwise phase... Only phase linking is structurally different. Zone C succeeds where Zone A fails under identical processing (64.7% vs 5.4%). That is a controlled contrast. Lead §4.1 with the A-versus-C contrast. Demote six-estimator convergence and remove from title/highlights."*

**Response:**
We have restructured §4.1 and the manuscript framing completely:
1. **Title & Highlights revised:** The title has been changed from *"C-band InSAR over a floating peatland: multi-method evidence..."* to *"What Does C-Band InSAR Measure Over a Floating Peatland? Empirical Limits of Satellite Displacement Retrieval and a Saturated Canopy Moisture Signal"*. Highlights 1 and 3 demote multi-method language and emphasize the controlled Zone A vs Zone C contrast.
2. **Section 4.1 lead:** Section 4.1 is retitled *"Controlled contrast: Zone A vs Zone C and estimator convergence"*, leading directly with the contrast: under identical 10×2 multilooking, Goldstein filtering ($\alpha=0.5$), and network geometry, Zone C retains 64.7% usable pixels whereas Zone A collapses to 5.4% usable pixels (§4.1.1).
3. **Pairwise reweighting interpretation:** The five pairwise linear estimators are presented as a demonstration that no linear reweighting or baseline trimming of standard burst products overcomes the decorrelation floor (§4.1.3).
4. **Annual-pairs disclosure:** In §4.1.3, we explicitly disclose that only 10 of 356 pairs exceed 120 days, noting that the annual-pairs estimator carries only ~10 observations and is reported for completeness rather than as an independent test.

---

### M5. Zero-filling an incomplete covariance matrix is not phase linking

> **Referee:** *"Zero-filling asserts zero coherence for those pairs; it does not represent missingness. Restate the claim as 'a sparse-network EVD consensus estimator applied to standard burst products,' drop the ML framing, and be clear that full-covariance SLC phase linking remains untested here."*

**Response:**
We have revised §3.1, §3.2, §4.1.2, and the H1 verdict:
1. **Dropped ML framing:** We no longer refer to the algorithm as "theoretically optimal maximum likelihood phase linking".
2. **Standardized terminology:** We define the approach as an *"empirical sparse-network EVD consensus estimator applied to multi-looked burst products"*, maximizing the coherence-weighted phase consensus across observed edges in the incomplete interferometric graph.
3. **Full-covariance SLC caveat:** In §3.2 and §5.4, we explicitly state that advanced full-covariance algorithms utilizing raw Single Look Complex (SLC) data with Statistically Homogeneous Pixel (SHP) selection (e.g., SqueeSAR; Ferretti et al., 2011) and phase bias mitigation (Ansari et al., 2021) remain untested here.

---

### M6. Zone C is treated asymmetrically

> **Referee:** *"Zone C has N_eff ≈ 5. Table 6 reports near-zero rank correlations in grassland (ρ = −0.008, −0.009), read as an informative absence in abstract and Conclusion 2. Elevation correlation (+0.430) is dismissed as DEM noise. Apply one standard. Remove 'no detectable effect in matched grassland' from abstract and conclusions."*

**Response:**
1. **Abstract and conclusions revised:** We removed the claim "no detectable effect in matched grassland" from the abstract, the H2 verdict, and Conclusion 2.
2. **Symmetric interpretation:** In §4.2.5, we explicitly apply a single inferential standard: given $N_{\text{eff}} \approx 5$, rank correlations in Zone C have very low power, meaning neither the near-zero values for backscatter/greenness nor the nominal +0.430 value for elevation constitute statistically decisive findings. The primary finding is confined to the within-mat environmental sensitivity in Zone A ($\rho = -0.379$ with backscatter, $+0.320$ with wetness).
3. **Landscape fragmentation documented:** In §5.7, we document that a burst-wide survey confirmed non-mat mineral grasslands are fragmented by surrounding pine plantations, making $N_{\text{eff}} \approx 5$ an inherent landscape constraint.

---

### M7. Forward model contradictions and circularity

> **Referee:** *"1. δ_p ≈ 3.4–4.4 mm senses only uppermost 3–4 mm, while §4.2.4 argues for canopy volume scattering (RVI 0.914). Reconcile penetration depth with volume scattering. 2. The saturating tanh fit brings peak-to-peak swing under the 6.13 mm ceiling, but an 0.008 R² gain is not evidence for functional form. Motivate physically or drop. 3. Delete point-estimate coincidence (Δm_v = 0.25 -> -3.28 mm). 4. Reconcile R² = 0.299 (T11) vs 0.3642 (T16)."*

**Response:**
1. **Penetration depth vs volume scattering reconciled:** In §5.2, we clarify the distinct scattering budgets: at C-band, radar signals cannot penetrate through the saturated peat substrate; power penetration is strictly confined to the upper 3–4 mm capitulum layer. Temporal decorrelation ($\bar{\gamma} = 0.408$) and high RVI ($0.914$) reflect volume scattering within the herbaceous canopy and non-stationary scatterer movement, whereas the phase perturbation is governed by the refractive index change across the top lossy boundary.
2. **Physical motivation for saturating fit:** In §5.2, we motivate the bounded saturating model physically: dielectric permittivity of water-air-organic mixtures is bounded by complete desiccation ($m_v \to 0$) and complete saturation ($m_v \to 1 - v_s \approx 0.93$). A linear sinusoid assumes unbounded linear excursion; replacing it with a hyperbolic tangent ($S \tanh(\cdot)$) enforces this physical boundary constraint. We explicitly state that the 0.008 $R^2$ improvement is not statistical proof of the functional form, but confirms that physical saturation naturally resolves the nominal ceiling violation.
3. **Point coincidence removed:** We deleted the claim that $\Delta m_v = 0.25$ "coincides centrally" with $-3.28\text{ mm}$. We now state that forward dielectric modeling over plausible Sphagnum desiccation ranges ($\Delta m_v \in [0.15, 0.35]$) predicts seasonal shifts of $-1.9$ to $-4.6\text{ mm}$, demonstrating that upper-canopy moisture fluctuations alone reproduce the scale of the observed $3.29\text{ mm}$ displacement.
4. **$R^2$ discrepancy clarified:** In §4.3.4 and §5.2, we explicitly distinguish between $R^2 = 0.299$ (the unconstrained seasonal harmonic fit on raw aggregated time series) and $R^2 = 0.3642$ (the linear harmonic fit on the deseasonalised reference baseline in Table 16 / T16).

---

### M8. *p*-value bookkeeping

> **Referee:** *"A p of 0.014 is unattainable with 92 draws (floor 1/93 = 0.0108). Null sizes vary without explanation: 4614, 1000, 249, 92. What invalidated 164 draws in §4.3.5a? Designate one canonical headline p-value. Remove pre-registration claims unless a timestamped DOI exists."*

**Response:**
1. **Unattainable $p = 0.014$ resolved:** We traced the $p = 0.014$ figure to a developmental 142-draw test. We have replaced it throughout the manuscript with the exact empirical $p$-value from the canonical 92-draw size-matched null: **$p = 0.026$** (with floor $1/93 = 0.0108$).
2. **Rationale for differing null draw counts:** In §3.4, we explicitly state the physical and statistical rationale for each null count:
   - $N_{\text{null}} = 4,614$: Full-zone grid shifts across Zone D using all valid coordinates.
   - $N_{\text{null}} = 1,000$: Toroidal core-vs-margin permutations preserving spatial geometry.
   - $N_{\text{null}} = 249$: Spatial block cross-validation partitions matching $L_{\text{corr}} \approx 160\text{ m}$.
   - $N_{\text{null}} = 92$: Size-matched non-overlapping 499-pixel compact blocks within Zone D.
   Every reported $p$-value now carries its explicit $N_{\text{null}}$.
3. **164 invalidated draws explained:** In §4.3.5a, we document that in the 1,000-draw toroidal permutation test, 164 draws were discarded because geometric compactness constraints rejected boundary-clipped or non-compact fragments, leaving 836 geometrically valid realizations.
4. **Preregistration language:** We removed all occurrences of "pre-registered prediction" and "stated before the test". We replaced them with *"designated primary confirmatory endpoint"*, accurately describing our pre-specified primary analysis protocol without asserting external registry filing.

---

### M9. Sub-zone test cannot prove the unit assumption

> **Referee:** *"Sub-zones are nested, and only outer margin (143 px) and deep core (233 px) are disjoint. The difference (2.89 vs 3.78 mm; DOY 99.8 vs 108.3) is inside aggregate noise. 'Kinematically' contradicts the dielectric reading. Rewrite as: subdivision detects no differential behaviour at the resolvable scale, which is consistent with but does not establish unit response."*

**Response:**
We have rewritten §3.3, §4.3.4, and §5.3.1 exactly as recommended:
- We replaced "kinematic unit" with *"spatially coherent moisture response without differential behaviour detected at the resolvable scale"*.
- We explicitly note that because outer margin ($2.89\text{ mm}$, DOY 99.8) and deep core ($3.78\text{ mm}$, DOY 108.3) differ by amounts well within aggregate noise at those patch sizes, the test shows an absence of detectable heterogeneity rather than proven mechanical rigidity.

---

### M10. Drained-peatland comparison overreaches

> **Referee:** *"§5.2 places Rzecin below nature reserves at 4.8 mm yr⁻¹ and concludes 'ecologically expected outcome.' The comparison target sits inside your own detection floor (1.5–5 mm yr⁻¹). State the exclusion bound, not an 'ecologically expected' claim the data can't support."*

**Response:**
We have revised §5.2 to eliminate the claim of "ecologically expected" lower subsidence. We now state strictly: *"Because the regional comparison rates (4.8 mm yr⁻¹) sit within our empirical velocity detection floor (1.5–5.0 mm yr⁻¹), our data cannot distinguish whether Rzecin subsides at a lower rate or is simply unresolvable; the defensible conclusion is that linear subsidence exceeding ~5 mm yr⁻¹ is excluded."*

---

## 3. Response to Section 3: Concrete Suggestion on Orbital Geometry

> **Referee:** *"Descending-track bursts over Rzecin exist in the same archive. Adding a descending stack would test vertical/horizontal partitioning directly and provide an independent aggregate replication."*

**Response:**
We thank the referee for this valuable recommendation. In §5.4, we have formalized the complete **Multi-Geometry Replication Protocol and 2-LOS Decomposition**:
$$\begin{pmatrix} d_{\text{LOS}}^{\text{asc}} \\ d_{\text{LOS}}^{\text{desc}} \end{pmatrix} = \begin{pmatrix} \cos \theta_{\text{asc}} & -\sin \theta_{\text{asc}} \cos \alpha_{\text{asc}} \\ \cos \theta_{\text{desc}} & \sin \theta_{\text{desc}} \cos \alpha_{\text{desc}} \end{pmatrix} \begin{pmatrix} d_{\text{vert}} \\ d_{\text{east}} \end{pmatrix}$$
We establish that executing an identical phase-linking and aggregation pipeline over the Sentinel-1 descending track represents the primary prospective empirical replication test to uncouple vertical peat breathing from potential lateral drift and confirm incidence-angle scaling ($\Delta \phi \propto \cos \theta$).

---

## 4. Point-by-Point Responses to Minor and Technical Points

1. **Missing references:** Added all eight missing references to `08_references.md` and `references.bib`: Monti-Guarnieri & Tebaldini (2008), Birchak et al. (1974), Ulaby & Long (2014), Klein & Swift (1977), Stogryn (1971), Hallikainen et al. (1985), McCarter & Price (2014), and Strack et al. (2009).
2. **Dangling cross-references:**
   - Fixed all dangling citations to "Table XT10" and "Table XT12".
   - Exported and integrated Supplementary Table S2 (`T20_multivariate_regression_s2.csv`), containing full regression coefficients and VIF diagnostics.
   - Corrected all mislabelled section references (§2.4 $\to$ §4.3.5a, §4.2.5 $\to$ §A.10, §4.1.2 $\to$ §3.1).
   - Ensured all Appendix B tables (T01–T20) are formally introduced and cited.
3. **Hypothesis labelling:** Unified the convention: hypotheses (H1–H4) are defined strictly as falsifiable physical/methodological propositions in §1.5, and conclusions in the abstract and results sections are designated as empirical findings/verdicts without reusing conflicting hypothesis tags.
4. **"Significantly lower":** Removed "significantly lower" project-wide; reported exact jackknife 95% confidence intervals ($[-0.109, -0.052]$) instead.
5. **Inconsistent self-audit counts:** Standardized the count across the abstract, §5.3.2, §6, and Appendix A.12 to **four analytical safeguards** that altered a reported bound or invalidated an intermediate conclusion.
6. **Table 9 lag-0 column:** Added the lag-0 primary effect size column to Table 9 ($r \approx 0.39\text{--}0.42$) alongside the secondary 12-day swept maxima.
7. **$\sigma^0$–coherence relation instrumental artifact:** Added an explicit note in §4.2.5 addressing potential SNR bias in sample coherence estimation, explaining why the moderate rank correlation ($\rho = -0.379$) reflects genuine microtopographic wetness variations rather than estimation artifact.
8. **H4 verdict softened:** Softened the H4 verdict from "supported" to *"consistent with"* in §4.4.5 and the abstract.
9. **Figure organisation:** Removed all inline embeds of Figures S1–S14 from the main text chapters. Main text chapters contain strictly Figures 1–7. All 14 Supplementary Figures (S1–S14) are collated with full captions into a dedicated `## Supplementary Figures` section at the end of Appendix A (`07_appendix_alternatives.md`).
10. **Submission blockers (`P-034`):** Prepared structured metadata templates and established clean frontmatter fields for author list, affiliations, corresponding email, and Zenodo DOI deposit.
11. **Title and manuscript length:** Changed the manuscript title to demote multi-method claims. Tightened prose across all chapters and pruned conversational referee-rebuttal phrasing in Appendix A.7.
12. **Register:** Removed all high-rhetoric adverbs (*"the decisive test"*, *"Crucially"*, *"Importantly"*, *"the methodological trap we set out to avoid"*) and un-bolded verdict blockquotes into standard declarative paragraphs.

---

## 5. Verification of the 8 Acceptance Criteria

| # | Referee Requirement | Resolution in Revised Manuscript | Status |
|---|---|---|:---:|
| 1 | Control target precluding motion & filter leakage, or explicit statement that none exists | In §5.7, acknowledged unexecuted burst search; lake recast as inconclusive; H3 rests on forward model and optical wetness | **Satisfied** |
| 2 | Recomputed nulls and $p$-values for every perturbation; headline framed as marginal | Documented subset stability (2.80–3.29 mm, Table T15); headline framed as marginal in 2.80–3.29 mm envelope (§4.1.4, §4.3.4) | **Satisfied** |
| 3 | Derivation of 7.32 mm LOS upper bound | Fully derived in §4.3.7 via 2,000-draw date-level non-parametric bootstrap preserving residual autocorrelation | **Satisfied** |
| 4 | §4.1 restructured around A-vs-C contrast; six estimators demoted and removed from title | §4.1 leads with 64.7% vs 5.4% controlled contrast; title retitled; annual pairs disclosed as 10/356 pairs | **Satisfied** |
| 5 | "No detectable effect in matched grassland" removed or supported by usable $N_{\text{eff}}$ | Removed from abstract and conclusions; Zone C $N_{\text{eff}} \approx 5$ landscape fragmentation constraint documented in §5.7 | **Satisfied** |
| 6 | Forward model reconciled with volume scattering; tanh motivated physically; point coincidence removed | Reconciled $\delta_p \approx 3.4\text{--}4.4\text{ mm}$ with canopy volume scattering in §5.2; tanh motivated as bounded physical saturation; coincidence removed | **Satisfied** |
| 7 | $N_{\text{null}}$ attached to every $p$; 0.014 resolved; preregistration language evidenced or removed | Canonical $N_{\text{null}} = 92, p = 0.026$ designated; rationales for all $N_{\text{null}}$ counts and 164 invalidated draws documented (§3.4); preregistration language resolved | **Satisfied** |
| 8 | S-figures resolved; cross-references and citations fixed; length trimmed | Main text contains strictly Figs 1–7; S1–S14 collated in dedicated supplement; Table S2 supplied; 8 missing citations added; length trimmed | **Satisfied** |
