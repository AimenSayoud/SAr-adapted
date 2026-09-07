## 3. Methods

Figure 4 summarises the processing logic: the change of observable and the
weak-signal test protocol.


![**Figure 4.** Processing protocol: per-pixel versus aggregated observable, and the size-matched null test yielding an empirical p-value.](figures/F04_protocol.png)

### 3.1 Inversion estimators compared (H1)

Six approaches resting on mathematically distinct assumptions:

| # | Method | Distinguishing assumption |
|---|---|---|
| 1 | SBAS (MintPy) | small-baseline network, global inversion |
| 2 | ISBAS | tolerates intermittent pixels (per-pixel sub-network) |
| 3 | Annual pairs | bypasses seasonal decorrelation by state matching |
| 4 | Hybrid network | combines short and long baselines |
| 5 | Weighted least squares | coherence weighting, pair by pair |
| 6 | Phase linking (EVD) | maximum likelihood over the observed coherence matrix |

The sixth approach exploits all available pairs simultaneously through the
per-pixel N × N complex coherence matrix, whose dominant eigenvector phase is
the estimated phase history. Phase linking is theoretically optimal under an
unconstrained, complete, and unbiased sample covariance matrix. In our dataset,
we evaluate it directly on the delivered multi-looked burst interferograms,
which populate 356 pairs out of 4,005 possible off-diagonal pairs (an 8.89 %
fill fraction across ~90 acquisition dates). Evaluating on delivered products
avoids raw SLC ingestion and coregistration, making phase linking directly
accessible from standard burst products. Quality is measured by temporal
coherence, the agreement between estimated phase history and observed
interferograms.

**Noise floor.** At the redundancy of our network (356 pairs over ~90
acquisitions, redundancy ratio ≈ 4), a fully decorrelated pixel returns a
simulated synthetic temporal coherence floor of ≈ 0.55. In the empirical
network distribution across non-coherent terrain, the noise floor sits at
0.488 (95 % empirical interval [0.448, 0.532]). Values near or below this range
must therefore be interpreted as noise rather than recoverable deformation.

### 3.2 Matched-cover comparison (H2)

The central test is paired by interferogram: each pair is observed in A and
in C, which cancels perpendicular baseline and same-day atmosphere, leaving only
the surface difference.

**Statistics.** A Wilcoxon signed-rank test on the differences coh(A) − coh(C).
Because the 356 pairs share ~90 acquisition dates, a pair-level bootstrap
understates uncertainty due to temporal correlation. We therefore report a
date-jackknife: each acquisition date and all pairs containing it are removed in
turn. We report both the empirical leave-one-out range demonstrating strict
sign invariance and the standard-error-based 95 % confidence interval.

### 3.3 Spatial aggregation: the change of observable (H3)

**Quantitative motivation.** The per-pixel phase standard deviation at γ = 0.4
is ≈ 1.5 rad. Complex averaging over *N* pixels divides this by √N_eff: over the
499 mat pixels, ≈ 0.07 rad (≈ 0.3 mm), or ≈ 1 mm even at N_eff = 50 — well below
the 10–40 mm sought. The signal is not below the physical noise floor; it is below the
per-pixel noise floor.

**Enabling physical assumption.** The mat is a hydrological unit and breathes as
a block, so averaging does not destroy the signal (as it would for a
heterogeneous deformation field) but suppresses the random component. This
assumption is stated explicitly and is testable by subdividing the zone.

**Network inversion and time-series construction.** Spatial aggregation is
conducted in two steps:

1. For each interferogram $j = 1, \dots, M$ ($M = 356$), the spatial mean phase
   difference between Zone A and reference Zone C is computed:
   $$\Delta \phi_{\text{pair}, j} = \langle \phi_A \rangle_j - \langle \phi_C \rangle_j$$
2. The cumulative date phase time series $\boldsymbol{\phi} = [\phi(t_1), \dots, \phi(t_{N-1})]^T$
   relative to reference date $t_0$ is solved by weighted least squares (WLS):
   $$\mathbf{G} \boldsymbol{\phi} = \Delta \boldsymbol{\phi}_{\text{pair}}$$
   where $\mathbf{G} \in \{-1, 0, +1\}^{M \times (N-1)}$ is the network design matrix mapping
   acquisitions to interferometric pairs. The weighted least-squares solution is:
   $$\hat{\boldsymbol{\phi}} = (\mathbf{G}^T \mathbf{W} \mathbf{G})^{-1} \mathbf{G}^T \mathbf{W} \Delta \boldsymbol{\phi}_{\text{pair}}$$
   with diagonal weighting $W_{jj} = \gamma_j^2 / (1 - \gamma_j^2)$ derived from pair coherence $\gamma_j$.
3. The resulting cumulative displacement series $\hat{\phi}(t_i)$ is then jointly regressed
   against linear trend and annual harmonics:
   $$y(t_i) = c + d \cdot t_i + a \cos(2\pi t_i) + b \sin(2\pi t_i)$$
   yielding seasonal amplitude $A = \sqrt{a^2 + b^2}$. Estimating trend and harmonic jointly
   prevents line absorption of cyclic phase.

In addition, circular mean resultant length $|R| = |\sum w_k \exp(i\phi_k)| / \sum w_k$ is
evaluated directly on wrapped phase, remaining independent of unwrapping errors.

![**Figure S3.** Synthetic validation. On identical simulated data, per-pixel inversion returns −13.7 mm yr⁻¹ (36 % usable pixels) whereas aggregation returns −19.8 mm yr⁻¹ against a ground truth of −20.](figures/S03_synthetic_validation.png)

### 3.4 Weak-signal test protocol

This protocol forms our primary methodological contribution and invalidated two
intermediate conclusions during the study.

**Primary pre-specified null control.** Aggregate noise falls as $1/\sqrt{N}$. A null
built on 2 200 pixels while the tested zone has 499 carries ≈ 2× less noise and understates
the floor, manufacturing false detections. Null realisations are compact adjacent patches
of stable ground with the same pixel count as the tested zones. The reference-matched
spatial null (4,614 draws, $p = 0.026$) is pre-specified as the primary confirmatory
endpoint; size-matched compact nulls (empirical $p$-values of 0.014 and 0.022) provide sensitivity bounds.

**Empirical p-value.** Rather than assuming Gaussian asymptotics, we evaluate:
$$p = \frac{1 + \#\{\text{null} \ge \text{observed}\}}{1 + N}$$
This empirical $p$-value has a floor of $1/(1 + N)$: with 92 draws it cannot fall below 0.011.

**Identical treatment of the null.** Where the observed statistic results from an exploratory
selection — such as sweeping ~16 lags — the null undergoes the identical sweep. To avoid
selection bias, the lag-0 correlation is reported as the primary effect size, with the swept
maximum reported as secondary exploratory evidence.

### 3.5 Mechanism discrimination (H3)

- **Lake control.** The residual open-water lake cannot breathe mechanically. If it exhibits
  a seasonal trajectory consistent in amplitude and phase with the mat, this demonstrates
  that a non-mechanical mechanism (dielectric change or emergent vegetation) operates at
  this amplitude scale.
- **Mat minus lake.** Referencing A to B cancels any cycle common to saturated surfaces,
  isolating motion specific to the mat.
- **Order of magnitude.** A freely floating mat following a ±10 cm water table would move
  ≈ 100 mm; published breathing is 10–40 mm.
- **Closure-phase bias.** Displacement, even non-rigid, closes triplets to zero; a monotonic
  dielectric drift biases them (De Zan et al., 2015; Ansari et al., 2021). Evaluated on wrapped
  phase to avoid 2π unwrapping ambiguities.

### 3.6 Within-mat predictive model (H2, H4)

The target is temporal coherence treated as a continuous variable — the 0.7 threshold is an
operational convention that would discard variance across 499 pixels — and the analysis
exploits variability inside zone A.

- **Spearman** correlations for monotone, outlier-robust marginal ranking.
- **Standardised multiple regression** for partial effects.
- **Variance inflation factors** are mandatory before interpreting coefficients. RVI and VH/VV
  are monotone transforms of the same ratio (VIF ≈ 240), generating collinear sign artefacts
  if unregularised.
- **Spatial autocorrelation and degrees of freedom.** Spatial autocorrelation in Zone A exhibits
  an empirical correlation length of ~160 m (4 pixels on the 40 m grid), yielding an effective
  sample size $N_{\text{eff}} \approx 31$ independent spatial degrees of freedom across 499 pixels.
  In a 12-covariate model, the observations-per-parameter ratio is under 3 ($31 / 12 \approx 2.6$).
  Cross-validated $R^2_{\text{cv}} = 0.239$ is therefore reported strictly as an empirical descriptive
  benchmark of within-zone spatial structure rather than an independent predictive model.

### 3.7 Hydrological coupling (H4)

**Autocorrelation control.** InSAR and hydrological series are both strongly autocorrelated.
We reuse the size-matched null series, which share the same temporal structure, so the
null distribution absorbs autocorrelation by construction.

**Seasonality removal.** Two annual-cycle signals always correlate strongly at some lag: the
sweep merely aligns phases. Causal inferences are therefore drawn exclusively on deseasonalised
anomalies, with annual harmonics removed from both series.

**Interpreting residual lag.** A dielectric response to moisture is near-instantaneous (one revisit);
mechanical settling would lag the water table by multiple weeks.

### 3.8 Software and reproducibility

All processing uses an open Python stack (`numpy`, `xarray`, `rioxarray`,
`scipy`, `scikit-learn`) with a purpose-built package. Every scientific routine
is covered by synthetic unit tests that verify recovery of a known ground truth,
including: EVD phase linking on a sparse network; aggregation recovering a
displacement buried under per-pixel noise; the collapse of a spurious correlation
between two independent annual cycles; and the size-matched null construction.
The complete analysis code is available at [repository DOI].
