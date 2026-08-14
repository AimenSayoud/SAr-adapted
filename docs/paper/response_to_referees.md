# Response to referees

**Manuscript:** *What does C-band InSAR measure over a floating peatland?
Multi-method evidence for a dielectric-dominated signal and an upper bound on
mat motion*

**Authors:** [A. Sayoud], [supervisor], [co-authors]

---

## Letter to the editor

We thank both referees for reports of unusual quality. Referee 1 worked
systematically through our statistics and protocol; Referee 2 attacked the
physics and the logical structure of our claims. Between them they identified several
places where our conclusions outran our evidence, and we have revised
accordingly.

We ran thirteen new tests in response to specific objections. Four of them
required us to withdraw or narrow a published claim, and we state those first
rather than leaving them to be found:

1. **The confidence interval on our headline amplitude is much wider than we had
   acknowledged.** A date-level bootstrap (2 000 resamples) gives 3.29 mm with a
   95 % interval of **[0.58, 7.32] mm**. Our published bound of "≤ 3.9 mm" was
   read off a point estimate. It is not supportable as written, and the refined
   "< 2.4 mm" second-level bound is withdrawn entirely.

2. **Referee 2's objection to our null construction was correct.** We built the
   null from two patches of stable ground, so the real reference zone never
   entered it. Rebuilding the null while keeping the actual reference zone gives
   a distribution **1.80× wider** in median, and the *p*-value moves from 0.014
   to **0.038**. The reference zone does carry most of the variance, exactly as
   the referee suspected. We now report the larger, correct *p*.

3. **The propagated interval kills our order-of-magnitude argument against
   flotation, and we had not noticed.** Carrying the upper interval through the
   line-of-sight-to-vertical conversion gives 8.66 mm, and through the coupling
   fraction of item 4 below it reaches 17 mm at *f* = 0.5 and 35 mm at
   *f* = 0.25 — overlapping the 10–40 mm range published for raised-bog
   breathing. The statement that our signal is "one to two orders of magnitude
   too small for flotation" is therefore **withdrawn**, along with the figure
   that plotted it on a logarithmic axis as though the separation were
   established. This is not a weakening of the bound; it is the loss of a
   separate argument, and we are grateful it was caught before submission.

4. **Our reading of the lake as "at the noise floor" does not hold.** Simulating
   the floor properly on our real network topology puts it at 0.488
   (90 % interval [0.448, 0.532]), and every zone — including the lake at
   0.584 — sits above that interval. The "internal validation" we claimed from
   the lake is withdrawn, and the dichotomy "mat ≈ noise, grassland ≈ signal" is
   replaced by the weaker claim the data support.

Referee 2 also identified a **logical inconsistency at the centre of the paper**,
which we accept in full, and which is the origin of item 3 above. If the phase centre sits in a saturated canopy volume
decoupled from the substrate — which is our own conclusion — then our phase does
not observe the peat, and a bound derived from it is a bound on *apparent
phase-centre displacement*. We have rewritten every statement of the bound
accordingly, changed the subtitle, and added the coupling assumption explicitly.

Two objections resolved in our favour, and both are new results:

- **The block-deformation assumption, which we had promised to test and had
  not, holds.** Repeating the amplitude retrieval at 499 → 250 → 125 → 60 pixels
  leaves it flat (×0.93) where pure noise would inflate it by ×2.88, while the
  scatter across draws rises as 1/√N exactly as expected.
- **The 27 mat pixels that survive our reliability threshold sit
  preferentially at the margin.** We present this as an observation with two
  live explanations rather than as a result: our own rank test is invalid on
  clustered pixels, and a mixed-pixel explanation is not excluded. See M6.

Finally, we have corrected the numerical inconsistencies both referees found.
They were transcription errors — numbers typed into the prose from an early run
and never revisited — and we have added a build-time check that compares every
reported figure against the exported data tables, in both directions, so the
manuscript cannot be built with a stale number again. We are grateful this was
caught; it was our own reproducibility discipline applied to the wrong layer.

A point-by-point response follows. We mark each item **[Done]**,
**[Conceded]**, **[Deferred]** (with the reason), or **[Disagree]** (with the
argument).

---

## Response to Referee 1

### Major issues

**M1. Text and Appendix B report different numbers. [Done]**

Corrected, and the cause addressed structurally. All values now come from a
single pinned run:

| Quantity | Was | Now | Statistic named |
|---|---|---|---|
| Wilcoxon *p*, A vs C | 2.2 × 10⁻⁴⁹ | **4.84 × 10⁻⁴⁶** | — |
| Paired Δ coherence | −0.069 | **−0.081** (mean); −0.050 (median) | both now reported |
| Date-jackknife range | [−0.0705, −0.0652] | **[−0.0842, −0.0774]** | on the mean |
| Fraction with A lower | 88 % | **89 %** | — |
| Zone D usable fraction | 23.1 % | **23.2 %** | — |

You correctly noted that Table 9 did not say *which* statistic −0.069 was. It
was neither: it came from a superseded run. We now report the mean and the
median together, because the gap between them (−0.081 against −0.050) is itself
informative — the distribution of paired differences is left-skewed, so a subset
of interferograms carries a much larger deficit than the typical one.

**One correction to the report, offered respectfully.** For zones C and D the
main text was right and Appendix B was wrong, not the reverse. Investigating the
discrepancy found a defect in our own code: the appendix computed the usable
fraction as `nanmean(values >= 0.7)`, and because `NaN >= 0.7` evaluates to
`False` rather than `NaN`, masked pixels were silently counted as failures. Zones
with complete coverage were unaffected — which is why A (0.054) and B (0.015)
agreed while fragmented C and edge-clipped D did not. The correct values are the
published **64.7 %** and **23.2 %**. Both paths now share one implementation.

The build now fails if any registered quantity in the prose disagrees with the
exported CSV, and separately if a superseded value is left behind anywhere —
because these figures repeat across abstract, results, conclusions and appendix,
and correcting four of five sites would otherwise pass silently.

**M2. Appendix A.4 declares atmosphere "excluded by construction". [Conceded in part; Done in part]**

We accept that a surface-type-dependent near-surface humidity anomaly is a
distinct mechanism from the common atmospheric screen, that it does not cancel in
our double difference, and that our nulls cannot absorb it because both null
patches share a surface type. It is now a separate entry in the A.2 alternatives
table rather than being folded into "atmosphere".

We have taken your discriminator. Our phase-to-displacement convention is
`d = −λ/(4π)·φ` (MintPy sign convention), so positive displacement is motion
**toward** the sensor. This is now stated explicitly in §3.3 with its mapping onto
the HyP3 `INSAR_ISCE_BURST` product, rather than left implicit. An increase in
atmospheric wet delay lengthens the path and produces apparent motion **away**
from the sensor; shallower penetration moves the phase centre **up**, toward it.
The observed correlation is positive — wetter surface, apparent uplift —
consistent with penetration depth and inconsistent with wet delay. We make this
argument explicitly in both §4.4.3 and A.4 now.

**We must flag one limitation honestly.** Our unit tests construct synthetic
phase *using* this same constant and then recover it, which is a self-consistent
round trip and cannot detect an error in the convention itself. Nothing in our
chain independently verifies that the delivered product's sign matches MintPy's.
The sign chain has several links (product → re-wrapping → eigen-decomposition →
double difference → scaling), any of which could invert it. Since your exclusion
turns entirely on this, we state the convention explicitly.

**We were, however, too pessimistic about what it would take to verify it.**
Three external checks need no field access, and we are doing the first before
submission: comparing the sign of our chain against **EGMS** line-of-sight
velocities over any non-zero-signal urban or industrial target in the same scene
— an afternoon's work, which also supplies an external reference. The sign of the
DEM-error term is fixed analytically by the perpendicular baseline recorded in
the metadata, and the sign of topographically correlated atmospheric delay is
likewise known; either gives independent confirmation. A chain of inference whose
first link the authors themselves declare unverified is not one we should ask you
to accept, and it will be verified rather than declared.

**M3. The synthetic validation validates the wrong observable. [Done]**

You are right that our only end-to-end validation recovered a velocity, the
quantity we elsewhere argue is uninformative. We have replaced it. Injecting a
3.29 mm annual cycle at our measured coherence (0.408) and effective sample size
(31):

| | Median recovered | IQR |
|---|---|---|
| Per pixel | 3.44 mm | 1.38 |
| Aggregated | 3.30 mm | 0.25 |

We report this as a **dispersion** result rather than a shift in the median,
because fitting an annual cycle to ~90 dates averages noise down even per pixel.
The accurate claim is that per-pixel estimation cannot *resolve* the signal — its
interquartile range is 5.5× wider — not that it cannot recover it. The previous
figure is withdrawn.

**M4. The pair-to-date inversion is not described. [Done]**

A fair criticism: we defined the observables per interferogram and then plotted a
date-indexed series without documenting the step between. §3.3 now specifies the
design matrix, the coherence weighting, the reference-date choice, and reports
post-fit residual statistics.

We also accept the consequence you drew. Our claim of immunity to unwrapping
errors applies to |R| and the closure-phase statistic, which are computed on
wrapped phase; it does **not** survive a network inversion to date-referenced
values, which requires phase consistency the wrapped input does not provide. The
claim is now scoped to the two statistics it actually covers.

**M5. The matched control has ~5 effective samples. [Partly done — the outstanding test is blocking]**

We took your option 3. Repeating the retrieval against eight independent compact
control patches gives amplitudes from 2.52 to 10.16 mm — a factor of four — with
seasonal phase spanning 140 days. Our observed value sits inside that range at
roughly the 43rd percentile.

We must report a limitation in our own test: those patches were drawn from all
external cover, not from land-cover-matched grassland. Part of the spread is
therefore the *control* varying rather than the mat, since different covers have
different dielectric seasonality. We are re-running the test restricted to the
WorldCover grassland class before drawing a conclusion, and will report it in the
next revision.

We should state the result without softening it. Four of the seven controls
give a **larger** amplitude than ours; our value sits at the 43rd percentile, an
implied *p* of about 0.6 against the variability of control choice. That is not
"our value lies within the range" — it is "our value is indistinguishable from
what an arbitrary control produces".

Our reading that the spread across *unmatched* cover is precisely why matching
was necessary may well be right, but you would be entitled to call it a
rationalisation rather than a test, because as things stand it is one. **We
therefore treat the grassland-restricted re-run not as a next-revision
improvement but as the test that decides whether H3 survives**, and we will not
submit before knowing its outcome:

- matched controls tight around ~1 mm while A − C gives 3.29 mm → matching does
  the work, H3 stands, and we gain the strongest available argument for our
  design;
- matched controls spread over 2–8 mm → the amplitude is an artefact of control
  selection and H3 falls.

**A consequence we had missed.** The seasonal phase across those controls spans
140 days. Our DOY 104 maximum is therefore not corroborative of anything, and we
have removed its use as physical confirmation ("consistent with spring swelling
at high water table") from §4.3.4 and the conclusions.

**M6. "A ≈ noise" is contradicted by our own multi-threshold table. [Conceded]**

Accepted, and independently confirmed by our noise-floor work (S5 below). Your
tail argument is correct: at threshold 0.65 the mat retains 23.6 % against 18.5 %
for the acknowledged-noise zone and 42.9 % for external cover — the mat is low
and *intermediate*, not noise. The dichotomy is removed. The weaker claim, which
is all our argument requires, now stands in its place: the mat is deprived of its
high-coherence tail to the point where per-pixel inversion is not supportable,
while retaining measurable structure above the fully decorrelated case.

**The 27 pixels above 0.7. [Partly done — an observation, not yet a result]**

You were right that a referee would ask, and right that it belonged in the
paper. The observation is that they lie preferentially at the margin: median
signed distance to the boundary 40 m, against 113 m for the mat as a whole.

**We must withdraw the statistic we first attached to it, and the reason is one
our own protocol should have caught.** We reported Mann–Whitney *p* = 0.004 on
27 pixels against 499, at a measured spatial correlation length of 160 m — four
pixels at our sampling. Contiguous pixels are not independent observations, so
the effective sample size is a handful and that *p* is badly inflated. This is
exactly the error our §5.3.2 rules were written to prevent; we applied them to
the seasonal null and not to our own new test. The test is being redone by
**toroidal permutation**, which shifts the surviving cluster rigidly inside the
zone so that every null realisation preserves its size, shape and internal
autocorrelation.

**And a second explanation is not excluded, which your own numbers point to.**
The median distance is 40 m — exactly one pixel. With multi-looking and
filtering applied upstream, the effective resolution cell is wider than the
sampling, so first-ring pixels are contaminated by stable surrounding terrain by
construction. "The anchored edge holds phase" and "edge pixels contain stable
grassland" predict the same observation, and we cannot separate them at this
resolution. It is presented as an observation with two live explanations, one of
them instrumental, and retained as a testable prediction for the coring
transects rather than as a result.

**The test we should have run instead, and now will.** We eroded the 65-pixel
lake, where it could not work. We did not erode the **mat**, which has 499
pixels and can easily afford it: one to two pixels of erosion leaves ~420 then
~350, and tests in a single pass whether the coherence deficit *increases* on
the pure interior, whether the seasonal amplitude survives, and whether closure
dispersion rises — while answering the margin question by interior-versus-ring
difference rather than by an invalid rank test. Our radial profile already
suggests it will go our way (interior plateau 0.40, exterior peak 0.47). This is
the highest-value test remaining to us and it needs no field access.

**M7. Rzecin is instrumented and we use an optical proxy. [Deferred, with a correction]**

**The correction first, because you identified a genuine hole.** The "water-table
proxy" in the appendix was never defined in the manuscript, and the column header
"coherence per WTD" implied a measured water-table depth we do not have. It is a
**30-day rolling cumulative-precipitation anomaly from ERA5**. It is now defined
where it is first used and relabelled throughout.

**And the consequential error you found is real.** §5.5 excluded rigid-body
motion partly on "no hydrological coupling of coherence", while our own appendix
shows coupling in every zone (*r* = −0.26 over the mat). The defensible statement
is that there is **no differential** coupling — it is present everywhere and
nearly identical across zones, so it does not single out the mat. Corrected.

On the in-situ series: we accept that "preferable" is not an adequate answer.
See *Work requiring field access* below for what we are doing and on what
timescale.

**M8. The subdivision test is promised twice and never performed. [Done — result favourable, one gap]**

This was the most productive item in either report. Results:

| Pixels aggregated | Median amplitude | Scatter across draws |
|---|---|---|
| 499 (whole zone) | 3.29 mm | — |
| 250 | 3.53 mm | 0.21 |
| 125 | 3.55 mm | 0.44 |
| 60 | 3.52 mm | 0.76 |

The amplitude is **flat** — ×0.93 from 499 to 60 pixels, where a noise-dominated
observable would inflate by ×2.88 — while the scatter across draws rises as
1/√N exactly as predicted. This is your first outcome: the signal is spatially
coherent across the mat, and the aggregation assumption is supported rather than
merely asserted.

We accept your point about the botanical survey. The 32 recognised plant
communities are in tension with a single "hydrological unit", and we now name
that tension explicitly.

**One limitation of this test, which we would rather state than have found.**
Every sub-patch of the mat keeps the *same* reference zone, the same regional
atmospheric screen and the same temporal window. An annual cycle carried by the
**reference** — whose effective sample size is 5 — would therefore be equally
flat under subdivision of the target. The test as run is blind to precisely the
alternative that still matters, and read that way it is the same fact as M5
rather than an independent confirmation.

The missing test is symmetric and costs an hour: **subdivide the reference at
fixed target**. Flat under both subdivisions demonstrates spatial coherence of
both terms; collapsing or exploding under subdivision of the reference
identifies the source. It is implemented and will be reported before
submission, and we do not present M8 as settled until it is done.

**M9. Site characterisation. [Deferred — requires field access]**

We accept all four sub-points and cannot answer (a) from the archive we hold. See
*Work requiring field access*. On (b), (c) and (d) we have revised the text now:
the site's drainage and extraction history is stated rather than implied, our
placement in the subsidence gradient is argued from the *observed condition
during 2022–2024* rather than from site history, the trophic classification is
presented as actively debated rather than settled, and trophic status and mat
buoyancy are kept as independent axes — the *Schwingmoor* character is asserted
only on buoyancy evidence, never inferred from fen status.

**M10. The lake control needs its own alternative-mechanism paragraph. [Conceded]**

Added, and we accept the position you propose: the match constrains mat motion
*conditional* on the lake being mechanically static and dielectrically analogous,
and the level-1 bound does not depend on it.

Your observation about the lake's RVI deserves a direct answer, because it is
sharp. The lake has the **highest** RVI of all four zones (0.993) while being
specular in σ⁰ (−15.41 dB). We use high RVI over the mat to exclude open-water
double-bounce; applied consistently, the same reasoning says the lake is not
behaving as open water either. We now state that the dominant scatterer in that
zone is most likely the vegetated fringe, which weakens "the lake cannot breathe"
as a clean physical control and is disclosed as such.

Referee 2 raised the related instrumental question, and our test there was
inconclusive — see §4.4 in the second response.

### Statistical issues

**S1. The Wilcoxon *p* assumes independent pairs. [Done]**

Removed from the abstract and the conclusions. The result now leads with the
effect size and its stability: a mean deficit of −0.081, negative in 89 % of
pairs, with a date-jackknife range of [−0.0842, −0.0774] that never changes sign
under removal of any single acquisition. The nominal *p* is retained once, in the
results table, explicitly labelled as uncorrected for pair dependence and
explicitly *not* the inferential basis. Your point about a 16 % relative
difference carrying an extreme *p*-value is taken and stated.

**S2. Several *p*-values sit at the resolution floor. [Partly done — one gap remains]**

The published null was raised to 2 000 draws (1 851 completed). The seasonal
amplitude gives *p* = **0.0216** against a floor of 0.00054 — a measurement
rather than a floor, and roughly 1.5× our published 0.014.

**But we must flag that this leaves the wrong null well-sampled and the right
one under-sampled.** The *p* that now matters is the one from the
reference-corrected null (R2 §4.5), and that rests on **184 draws with 7
exceedances**. Its exact binomial 95 % interval is **[0.015, 0.077]**, which
**straddles 0.05**. We cannot rest H3 on a *p* whose own interval contains the
conventional threshold. The corrected null is being raised to ≥ 5 000 draws —
machine time and nothing else — before submission.

We also owe you a note on the 7.5 % of draws that failed to place. A placement
failure rate that high means the null is spatially constrained, and if accepted
locations differ systematically (scene edge, particular covers or slopes) the
null is biased. We will report the map of accepted locations and the sensitivity
of *p* to the placement constraint.

**S3. No multiplicity control across the forcing family. [Deferred — with your third option adopted]**

We have taken your cheapest option and pre-registered the primary test. The
mat-referenced optical wetness at near-zero lag was our a priori mechanism, it is
now the single primary test, and the remaining six forcings are demoted to
exploratory and reported as such. The machinery for a family-level null (the best
forcing at its best lag, with the null undergoing the same maximisation) is
implemented but not yet run; we will report it if the editor prefers that route.

**S4. No confidence interval on the headline 3.29 mm. [Done — and it is the most consequential change in this revision]**

Date-level bootstrap, 2 000 resamples, consistent with the jackknife treatment
elsewhere:

> **3.29 mm, 95 % CI [0.58, 7.32] mm** (mean 3.45, sd 1.73)

We draw the consequences rather than reporting the interval and moving on:

- The lower bound sits **below** our null median (0.85 mm), so some resamples are
  indistinguishable from stable ground. This does not contradict *p* = 0.02 — the
  null test asks whether the amplitude exceeds what stable ground yields, the
  interval asks how well it is determined, and "detectably non-zero but poorly
  determined" is a coherent position — but it must be stated.
- **The "≤ 3.9 mm" bound is withdrawn as written.** Propagating the upper
  interval gives a bound of ≈ **7 mm** on apparent phase-centre displacement.
- **The refined "< 2.4 mm" second-level bound is withdrawn entirely.** It was
  derived from a point estimate and does not survive its own uncertainty.
- **Propagated to vertical and through the coupling fraction, it withdraws our
  order-of-magnitude argument as well** — see R2 §4.1.

**A cross-validation we failed to claim.** The bootstrap gives 3.29 mm with a
standard deviation of 1.73, i.e. **1.90σ from zero**, a one-sided *p* of 0.028.
Our corrected empirical null gives 0.038. These are two methodologically
independent routes — one resampling acquisitions, one permuting terrain — and
they agree to within a factor of 1.4. We had presented them merely as not
contradicting each other; they do better than that, and this is direct evidence
that the protocol is calibrated.

**A caveat on the scheme, which probably makes our interval too wide.** An
i.i.d. bootstrap over dates destroys the temporal sampling structure, and we are
fitting an *annual cycle*, whose conditioning depends directly on how evenly the
year is covered. Resamples that over- or under-represent parts of the year
degrade the harmonic fit and inflate the interval. [0.58, 7.32] is therefore
likely conservative. We will report three intervals and discuss the differences:
the analytic interval from the harmonic fit's covariance, a residual (wild)
bootstrap preserving temporal coverage, and the i.i.d. date bootstrap as the
conservative bound.

**S5. The 0.55 noise floor is asserted, not documented. [Done — and it changed a conclusion]**

Now documented with design, realisations, distribution and topology sensitivity.
Simulating the *actual estimator* — the eigen-decomposition scored against the
interferograms it was fitted to — on our real network topology, over 400
realisations of a fully decorrelated pixel:

> **median 0.488, 90 % interval [0.448, 0.532]**, at redundancy 4.0

Topology sensitivity, which you asked for explicitly:

| Pairs | Redundancy | Median floor | p95 |
|---|---|---|---|
| 178 | 2.0 | 0.681 | 0.748 |
| 267 | 3.0 | 0.569 | 0.603 |
| **356** | **4.0** | **0.504** | **0.540** |
| 534 | 6.0 | 0.413 | 0.443 |
| 712 | 8.0 | 0.366 | 0.387 |

The floor moves by nearly a factor of two across plausible networks, so we now
quote it with its network and its spread rather than as a bare number.

**Two consequences we accept.** First, our published 0.55 sits between
redundancy 3 and 4 in this table. We accept that "we cannot reconcile our own
published number" is not a tenable position in a revision, since the floor
governs the 0.7 threshold, hence the usable-pixel fractions, hence the central
argument. We are identifying which assumption differs and will report a single
reproducible value.

**We also accept that the 0.7 threshold is now post-hoc.** It was chosen to sit
comfortably above 0.55; with the floor at 0.488 [0.448, 0.532] that
justification is gone. Rather than defend a threshold we are removing the
privilege given to any: the full usable-fraction-against-threshold curve becomes
the primary result, with the floor and its interval shown as a shaded band. This
is more robust and more informative than a single number, and it answers M6 in
the same stroke.

**And we adopt your threshold-free reformulation**, which we think is the right
replacement for the dichotomy we withdrew:

| Zone | Temporal coherence | Excess over floor (0.488) |
|---|---|---|
| C — matched grassland | 0.734 | 0.246 |
| D — other cover | 0.639 | 0.151 |
| **A — mat** | **0.604** | **0.116** |
| B — lake | 0.584 | 0.096 |

> The mat retains **47 %** of the matched grassland's excess coherence above the
> network floor.

That is quantitative, defensible, uses no arbitrary threshold, and survives the
floor moving. Second, and more importantly, **every zone lies above the
95th percentile of this null** — including the lake at 0.584. Our claim that the
lake "sits at the floor" and thereby validated the chain internally is
**withdrawn**. This reaches your M6 from a second direction.

**S6. Reconcile the theoretical noise budget with the empirical null. [Partly done]**

The comparison is now made. Our published null has median 0.85 mm and p95
1.93 mm; the null rebuilt against the real reference zone has median 1.53 mm and
p95 2.94 mm. The gap between them **is** the answer to your question: the
empirical null is dominated by the reference term, not by the target, and the
1/√N argument alone does not capture it. We present this as validation of the
protocol's necessity rather than of the 1/√N heuristic.

**S7. The closed-triplet count needs explaining. [Done]**

You are right that leaving this open was disproportionate to its difficulty: the
number of closed triangles is fully determined by the pair list. We now compute
it two independent ways — trace(A³)/6 on the adjacency matrix and direct
enumeration — and require them to agree, alongside the count a random graph of
the same density would give. The gap between the two is the explanation: HyP3
pairs follow constrained perpendicular and temporal baselines, so the network is
far from a random graph and closes correspondingly fewer triangles at equal
density. Both figures are reported.

This also matters for how we describe the closure result, since we originally
called a prediction falsified *because* the network holds only 518 triplets —
see R2 §4.8, where that description changes.

### Moderate issues

**Mo1. Qualitative-only reporting of four estimators. [Deferred]** Accepted. Four
of six are reported as "no usable pixel" with no numbers, which does not support
an argument resting on six comparable failures. Producing usable-pixel counts and
residual statistics on the same network for all six requires re-running three
chains and is scheduled for the next revision.

**Mo2. The descending track is unused. [Deferred]** Accepted as the cheapest
independent check within reach, and we agree its absence will be raised. It is
scheduled; the products are free and the chain applies unchanged.

**Mo3. The freeze test is thin and over-interpreted. [Conceded]** "Cold" is now
defined (both acquisitions below 273.15 K, ERA5 2 m temperature), the subsample
size is stated (31 of 356 pairs), and the claim it supports in §5.5 is reduced
accordingly rather than defended.

**Mo4. The incidence-angle correction is over-framed. [Done]** You are right that
≈ 32° is the expected value for the near-range sub-swath and that presenting it
as a discovery invites the observation that using the correct sub-swath geometry
is standard practice. Reworded as a caution — read incidence from product
metadata rather than assuming a mid-swath value — with the 9 % propagation into
any displacement bound retained, since that part is worth flagging.

**Mo5. Imprecise counts. [Done]** The acquisition count is now stated once (90)
and used consistently; "89 date increments" is expressed as such where the
distinction matters. Winter pairs are 108 of 356, given as 30.3 %.

**Mo6. Zone C fragmentation disclosed but not resolved. [Conceded — see M5 and §4.5]**
Your observation about the asymmetry was the sharper half of this item, and
Referee 2 reached it independently. Testing it directly confirms it: the null
built on compact patches is 1.80× narrower than one keeping the real fragmented
reference. Control and null did not share construction geometry, and that
mattered.

**Mo7. Calibrate the certainty language. [Done]** All four instances revised:
"every observation points to a single mechanism" → consistent with; "the signal
**is** dielectric" → replaced by the three-control formulation you suggested,
which is both accurate and stronger; "excluded" for rigid-body motion → scoped to
*no differential* coupling; "excluded by construction" for atmosphere → replaced,
per M2. We have followed your instruction not to hedge generally — the emoji
verdict markers and the internal test labels are gone, replaced by precise scope
statements rather than vague qualifiers.

### Editorial

**E1. Length. [Deferred — planned, not improvised]** We accept the diagnosis and
the consolidation targets, including the exact duplications you identified
(Figure 1 against Table 1, Figure 5b against Table B3, Figure 11b against
Table 15, and the identical Tables 19 and A2). Target: a main text of 8 figures
and 6 tables, with Appendix B moved to supplementary material alongside the CSVs,
and the short 2–5 row tables rewritten as sentences.

**E2. Abstract length. [Deferred]** Will be cut to under 350 words. The
four-hypothesis structure is preserved, as you recommend; the numbers within each
hypothesis are what will go.

**E3. Highlights. [Deferred]** Accepted on both counts. They will be brought
under 85 characters and reordered so that what the sensor *does* measure leads,
with the six-estimator result second — which also matches the retitling below.

**E4. Question-form title. [Done — proposed]** A declarative alternative is
proposed at the end of this letter, incorporating Referee 2's requalification of
the bound.

**E6. Reference completeness. [Deferred]** The Hrysiewicz volume record will be
verified and the three `[to verify]` markers completed, Milecka in particular
since it now supports three separate claims.

**E7. Data and code availability. [Deferred — before submission]** The repository
will be archived with a DOI before submission rather than after. We accept that a
placeholder is inconsistent with claiming reproducibility as a contribution.

**E5. Figure 2b legend is in French. [Done]** Corrected, and the cause removed
rather than the symptom. The shared plotting code hard-coded French zone labels
and the English export never overrode them, so the leak affected **three**
figures and one table, not one — and could not have been fixed in the document,
because the text is rendered into the image. Labels are now set explicitly by the
calling notebook before any figure is drawn.

---

## Response to Referee 2

### 4.1 The bound on mat motion is logically inconsistent with our own conclusion. [Conceded in full]

This is the most important item in either report and we accept it without
reservation. If the phase centre lies in a saturated canopy volume decoupled from
the substrate — our own conclusion — then the phase does not observe the peat, and
our bound is a bound on **apparent phase-centre displacement**.

We have quantified the gap rather than merely acknowledging it. Writing *f* for
the fraction of the phase that follows mat motion, the observable is *f* × true
motion, so our measured bound implies:

| Coupling *f* | Implied bound on vertical mat motion |
|---|---|
| 1.00 (phase centre rigid with the surface) | 8.7 mm |
| 0.75 | 11.5 mm |
| 0.50 | 17.3 mm |
| 0.25 | 34.6 mm |
| 0.10 | 86.6 mm |

The table is anchored on the **upper 95 % interval converted to vertical**
(7.32 mm ÷ cos 32.26° = 8.66 mm), not on the point estimate. Anchoring it on
3.9 mm, as our first draft of this letter did, would have applied an uncertainty
correction to a quantity from which we had just removed the uncertainty.

**And *f* is not constrained by our data.**

**The consequence we initially failed to propagate.** Published raised-bog
breathing spans 10–40 mm. Our table reaches 10 mm at *f* = 0.87 and covers the
whole published range by *f* = 0.25. **Our measurement is therefore compatible
with peatland breathing of the published magnitude**, and the order-of-magnitude
argument is withdrawn wherever it appears — abstract, §4.3.5(c), the summary
table, the literature-context figure, and the conclusions. What survives is a
cleaner statement and, we think, a more useful one:

> C-band InSAR on this surface cannot distinguish an absence of motion from
> motion of peatland-breathing amplitude. Every statement of the bound is
reworded, the coupling assumption is added explicitly to the assumption list, and
the subtitle is changed.

On the attenuation bias, we investigated and found the situation is more
complicated than either of us stated — in a way that partly favours the paper.
Two effects act in opposite directions:

- **Estimator bias from noise.** Amplitude is √(a² + b²), a positive function of
  two noisy coefficients, so noise biases it **upward**, not toward zero. We
  measured it: ×1.070 at per-pixel noise, **×1.003 after aggregation**. The
  aggregated 3.29 mm is therefore *not* inflated by noise.
- **Signal loss from decoupling**, which is your point and which is real, but is
  a multiplicative loss on the signal rather than an additive noise effect. It is
  the coupling fraction *f* above.

We now report both, separately, because reporting either alone misstates the
uncertainty.

We think the reframing strengthens the paper. "C-band cannot constrain mat motion
over this surface" is a cleaner instrumental-limit result than a number, and it is
what our own appendix already anticipated in noting that a laser showing 20 mm
against our 3.3 mm would quantify C-band's insensitivity directly.

### 4.2 The "six estimators including the MLE" claim does not hold. [Conceded, with one partial defence]

We accept the substance and have taken your option (b), downgrading the claim to
"no processing chain available on delivered burst products succeeds" — which we
agree remains publishable and is what we actually demonstrated.

We measured the matrix fill you identified: **8.9 %** (356 of 4 005 possible
pairs across 90 dates, redundancy 4.0). An eigen-decomposition on a matrix that
sparse, assembled from multi-looked and filtered products rather than from an
adaptively estimated sample covariance matrix, cannot carry "if this fails,
nothing will". The abstract and conclusions now carry the reservation that
previously appeared only in the scope paragraph.

**One partial defence, offered for the record rather than to resist the point.**
Your third sub-argument — that all six estimators inherit the same upstream
unwrapping, so their common failure may have a common upstream cause — holds for
five of them. It does not hold for the phase-linking path, which **re-wraps**
before the decomposition; unwrapping errors are multiples of 2π and cancel under
re-wrapping. The multi-looking and filtering objection stands for all six, and the
sparsity objection stands. We state the narrower claim.

### 4.3 The "floating" premise is not established. [Conceded — requires field access]

Accepted, including that our current defensive paragraph mitigates without
resolving. We cannot quantify the buoyant fraction from the archive we hold. See
*Work requiring field access*.

One new piece of evidence bears on this from an unexpected direction: the only
mat pixels that hold interferometric phase are significantly concentrated at the
margin (Referee 1, M6). That is consistent with an anchored perimeter and a
decoupled interior, and it is a testable prediction for the coring transects.

### 4.4 Possible contamination of the lake control. [Conceded — our test was inconclusive]

We ran the erosion test you proposed and must report that it does not settle the
question, for a reason that is itself the answer:

| Erosion | Pixels | B − C | A − B |
|---|---|---|---|
| 0 px | 65 | 2.63 mm | 0.90 mm |
| 1 px | 26 | 2.22 mm | 1.27 mm |
| 2 px | 4 | — | too few pixels |

Eroding by one pixel removes 60 % of the zone. That alone raises the noise floor
by ×1.58, so the observed 16 % change in B − C is not resolvable, and two pixels
of erosion leaves four pixels. **A 65-pixel zone cannot be eroded and still
tested.** We therefore concede the instrumental contamination objection as a
declared limitation rather than claiming a passed test, and note that a
higher-resolution water mask is the route to settling it.

**But the physical argument does not depend on *this* lake, and we had not seen
that.** What the control requires is a surface that cannot breathe mechanically,
in the same scene, under the same atmospheric screen. There are larger water
bodies elsewhere in the burst, not enclosed by the mat, which would give a clean
"cannot breathe" control that is large enough to erode and free of the vegetated
fringe. This is pure reprocessing and it is the way to save our most elegant
control rather than concede it. It is scheduled.

On the RVI: before concluding "vegetated fringe" we will check the absolute VH
level against the Sentinel-1 IW noise-equivalent sigma-zero (≈ −22 dB). Over calm
water both polarisations approach the instrument noise floor and the VH/VV ratio
— hence RVI — becomes noise-dominated and tends high with no vegetation at all.
We note that σ⁰ VV = −15.41 dB at 32° incidence is too bright for calm open water
(−20 to −25 dB would be expected), which does support the vegetated-fringe or
roughened-water reading; but we would rather settle it with the noise-floor check
than concede on an intuition, and the conclusion is firmer either way.

### 4.5 Null construction: matched in size, but is it matched in structure? [Conceded — you were right]

Our null was built as (null₁ − null₂): a single compact patch of stable ground
cut into two adjacent halves, size-matched to the real zones. **The real
reference zone never entered it.** Your inference was correct.

Rebuilding the null so that only the target is randomised while the actual
reference is retained:

| Null construction | *n* | Median | p95 | *p*-value |
|---|---|---|---|---|
| Published (null₁ − null₂) | 188 | 0.85 mm | 1.93 mm | 0.021 |
| **Keeping the real reference** | 184 | **1.53 mm** | **2.94 mm** | **0.038** |

The correct null is **1.80× wider** in median. The fragmented reference zone does
carry most of the variance, and *p* moves from 0.014 as published to **0.038**.
We now report the larger value and the construction it comes from. The signal
survives, with a materially reduced margin.

On the corollary: we agree that an *N*_eff of 5 for the reference is itself
troubling. See Referee 1's M5 for the multi-control work, its result, and the
limitation in our own version of that test.

### 4.6 Numerical inconsistencies between text and Appendix B. [Done]

Corrected as set out under Referee 1's M1, including the code defect that made
the appendix wrong for two rows rather than the text. The 2.2 × 10⁻⁴⁹ figure is
removed from the abstract in favour of the jackknife interval, as you recommend.

Two rows in your table we should address specifically. The null-median
discrepancy (0.86, 0.83, 0.568) was three different quantities reported with the
same label — the amplitude null median, the same after winter exclusion, and the
exported table's own value; they are now distinguished. The *R*²cv difference
(0.239 ± 0.022 against 0.244 ± 0.018) was a figure regenerated after a
collinearity fix while the text was not; now both come from one run.

### 4.7 Missing literature, and it is central. [Deferred — and we agree it is the highest-return item]

We accept this fully. Building the dielectric argument without the literature
that models moisture-induced phase directly is a real gap, and we agree that a
quantitative comparison of our 3.29 mm against the predictions of that framework
is the strongest and cheapest strengthening available to us. It is the first item
in our revision queue. The peatland-InSAR review and the British ISBAS work will
also be incorporated.

We flag one dependency: the comparison is only meaningful once we have the
predicted magnitude for our wavelength and moisture range, which is a reading
task we have not completed. We prefer to defer it one revision cycle and do it
properly than to cite the framework without engaging its numbers.

### 4.8 Unresolved tension: no closure-phase bias. [Conceded — and it changes how we describe the result]

You are right that this is a tension with our retained conclusion and that we
treated it only as a falsified prediction. We have done the power calculation you
asked for:

| Zone | Triplets | SE (rad) | Detectable at 2σ | Equivalent |
|---|---|---|---|---|
| **A (mat)** | 518 | 0.0568 | **0.114 rad** | 0.50 mm |
| B (lake) | 518 | 0.0584 | 0.117 rad | 0.52 mm |
| C (grassland) | 518 | 0.0249 | 0.050 rad | 0.22 mm |

The observed bias over the mat is 0.090 rad — **1.6σ**, below our detection
threshold.

**This changes the description materially, and not in our favour rhetorically.**
We had reported "predicted a bias at ≈ 5σ, found 1.6σ" as a *falsified
prediction*. The observed bias is in the expected direction at 1.6σ — weak
positive evidence, not a refutation.

**We must be careful not to over-correct, and the honest position is that the
label depends on a number we do not yet have.** If the moisture-phase framework
predicts a bias of order 0.5 rad for our wavelength and moisture range, then our
non-detection at a 0.114 rad limit **refutes** the dielectric hypothesis — the
test is not underpowered, it is discriminating and negative. If the prediction is
of order 0.05 rad, the test is genuinely underpowered and "inconclusive" is
right. **We therefore cannot describe this test before completing §4.7**, which
is why that item is not deferrable: it is the calibration that decides whether
our most direct test supports our main conclusion or contradicts it.

One figure we had and did not use: the mat's closure standard error (0.0568 rad)
is **2.28×** the grassland's (0.0249). That is the same information as our
published ×3.2 dispersion ratio, reached independently.

### 4.9 The expected order of magnitude is unsourced. [Conceded — requires the site archive]

Accepted. The ±10 cm water-table amplitude was assumed rather than sourced, and
you are right that if the true amplitude over 2022–2024 is 2–3 cm the expectation
drops to ~25 mm and the argument loses much of its force. This is bound up with
the in-situ series request below.

### 4.10 No in-situ validation on an instrumented site. [Deferred — see below]

We accept the framing: on a paper with a hydrological conclusion, the absence of
any ground data from an instrumented site reads as insufficient effort rather
than an inherent limit. See *Work requiring field access*.

### 4.11 Reproducibility and finish. [Partly done]

- Repository archiving with a DOI: before submission, not after. Accepted.
- The unresolved `[X]` and the three `[to verify]` markers: being completed.
- The aggregated-series construction: now documented (Referee 1, M4).
- The 0.55 floor simulation: now documented with distribution and topology
  sensitivity (Referee 1, S5).
- Optical/radar temporal alignment and the delivered product's looks and
  filtering: now stated.
- Null draws raised to 2 000 (Referee 1, S2).

### 4.12 Structure and balance. [Partly done]

We accept the redundancy diagnosis and will move §4.5 into the appendix and
§4.3.7 to the head of the Discussion, as you suggest.

On the appendix documenting errors corrected during the analysis: we have taken
your advice on framing while keeping the content, since Referee 1 asked
explicitly that it be preserved. It is now presented as robustness checks that
modified our conclusions, without the vocabulary of fault. We note that this
revision adds three more entries to it, all of which weakened our own claims.

### 4.13 Style and scientific English. [Partly done]

- Bold: reduced substantially. The manuscript carried 483 emphasised segments,
  188 in the Results alone; we accept that this reads as typographic insistence.
- Rhetorical section titles: neutralised, as you identify them.
- "A ≈ noise, C ≈ signal": removed and replaced with quantities and their
  uncertainties (Referee 1, M6).
- Terminology: "floating peatland", "floating fen" and "floating mat" harmonised
  on one defined term; *Schwingmoor* defined at first use in the abstract.
- "idem" in the estimator table replaced; em-dash density reduced.

---

## Summary of new analyses in this revision

| Test | Objection | Outcome | Effect on the paper |
|---|---|---|---|
| Coupling-scaled bound | R2 §4.1 | bound is on phase centre; *f* unconstrained | claim **withdrawn**, replaced by a scoped statement |
| Bootstrap CI | R1 S4 | 3.29 mm [0.58, 7.32] | bound **withdrawn**; order-of-magnitude argument withdrawn |
| Null keeping the real reference | R2 §4.5 | 1.80× wider; *p* → 0.038 | *p* **restated**; more draws required |
| Noise-floor simulation | R1 S5 | 0.488 [0.448, 0.532] | "lake at floor" **withdrawn**; threshold replaced by a curve |
| Matrix fill | R2 §4.2 | 8.9 %, redundancy 4.0 | MLE claim **narrowed** to burst products |
| Eight independent controls | R1 M5 | 2.5–10.2 mm, phase spread 140 d | H3 **pending** the matched re-run; DOY corroboration removed |
| 2 000 null draws | R1 S2 | *p* = 0.0216 | *p* **restated**, off the floor |
| Estimator noise bias | R2 §4.1 | ×1.003 after aggregation | amplitude **confirmed** not noise-inflated |
| Subdivision, 499→60 px | R1 M8 | flat ×0.93 against ×2.88 for noise | aggregation assumption **supported**, pending the symmetric test |
| Seasonal synthetic validation | R1 M3 | aggregation IQR 5.5× tighter | validation **replaced** on the right observable |
| Closure detection limit | R2 §4.8 | limit 0.114 rad; observed 1.6σ | description **pending** the predicted magnitude |
| Lake-mask erosion | R2 §4.4 | zone too small to erode | limitation **declared**; larger control scheduled |
| The 27 surviving pixels | R1 M6 | preferentially marginal | statistic **withdrawn**; re-test scheduled |

Six claims withdrawn or narrowed, three restated, two supported, two pending.

---

## Work requiring field access or archive requests

Both referees identified, correctly, that the items which most limit this paper
cannot be fixed by reprocessing. We separate them here so the editor can judge
what is achievable within a revision cycle.

**Requested and pending — needed for the next revision**

1. **In-situ water-table series, 2022–2024** *(R1 M7, R2 §4.9, §4.10)*. Requested
   from the site's long-term monitoring programme. This single dataset addresses
   three separate objections: it replaces the optical proxy in our hydrological
   test, it sources the water-table amplitude underpinning our order-of-magnitude
   argument, and it provides a forcing whose temporal structure differs from
   temperature. If it proves unavailable for our window we will state the specific
   reason — gap, spatial representativeness, or access — rather than describing it
   as merely preferable.

2. **Stratigraphic and probing records** *(R1 M9a, R2 §4.3)*. Requested from the
   site archive, to quantify what fraction of the 79.84 ha study zone is confirmed
   buoyant mat. Our new margin result gives this a specific prediction to test: an
   anchored perimeter with a decoupled interior. Should a significant part prove
   to be grounded peat, our bound applies to a mixture, which makes it
   conservative — but we will say so rather than leave it implied.

**Scheduled fieldwork — beyond this revision cycle**

3. **Surface laser** *(R2 §4.10, and the coupling fraction of §4.1)*. This is now
   the single most valuable measurement available to the project, and its role has
   changed because of this revision. It is no longer a cross-validation of a
   displacement we claim to measure; it is the only way to constrain the coupling
   fraction *f* on which our entire bound depends. As our appendix already argues,
   a laser showing 20 mm of real motion against our 3.3 mm signal would quantify
   C-band's insensitivity to this surface directly — a stronger and more useful
   result than any successful agreement.

4. **UAV microtopography**, to test the hummock-hollow model and the internal
   heterogeneity that our subdivision test probes only at the radar scale.

**Reprocessing, scheduled, no field access needed**

5. **Descending-track reproduction** *(R1 Mo2)*. Products are free and the chain
   applies unchanged.
6. **Quantitative outcomes for all six estimators** *(R1 Mo1)*.
7. **Multi-control test restricted to matched grassland** *(R1 M5)*.
8. **Phase linking on an SLC sample covariance matrix** *(R2 §4.2 option a)*,
   should the editor consider option (b) insufficient.

**Literature, in progress**

9. **The moisture-induced-phase framework** *(R2 §4.7)*, with a quantitative
   comparison against our 3.29 mm, and the closure-bias magnitude it predicts
   feeding directly into §4.8.

---

## Changes to title and framing

Both reports converge on the point that our subtitle claims more than we
demonstrate. Referee 1 notes that some journals discourage interrogative titles;
Referee 2 asks for the bound to be requalified. We propose:

> **C-band Sentinel-1 cannot constrain vertical motion of a saturated floating
> peatland: coherence deficit, phase-centre decoupling, and a weak-signal test
> protocol**

Our first draft of this response proposed "*tracks surface wetness rather than
mat motion*". We have withdrawn it. It would have promoted H4 — a moderate
correlation on anomalies, an optical proxy, a sign convention we had declared
unverified, and no measured water table — to the title, making our least secure
conclusion the most prominent claim. The second half, "rather than mat motion",
also sits badly with the coupling table: we cannot exclude mat motion, only
report that our phase does not observe it.

The formulation above states a non-capacity assertively, which is a result, and
names the three things we actually demonstrate. "Floating peatland" is retained
for indexing. "An upper bound on mat motion" is removed from the subtitle
entirely, and the bound appears in the text only as a constraint on apparent
phase-centre displacement under a stated coupling assumption.

We are grateful to both referees. The paper makes weaker claims than it did, and
we think it is a better paper for it.

We should be accurate about what has survived so far, since one claim in our
first draft of this letter did not. The spatial coherence of the aggregated
signal under subdivision is a genuine new result, pending the symmetric test on
the reference. The marginal concentration of the surviving pixels is **not** yet
a result — it has not been correctly tested, and we would have been claiming
exactly the pattern we have just spent this revision correcting.

What has been most useful is not any single test but the fact that the protocol
turned on its authors. It was built to keep us honest about weak signals; applied
to our own conclusions, it withdrew four of them. That is the strongest evidence
we can offer that it works, and it is the part of this paper we would most like
to see used elsewhere.
