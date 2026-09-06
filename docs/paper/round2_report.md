# Round-2 response report — Phase L

**Rzecin floating peatland, Sentinel-1 C-band InSAR**
Status report on the nine tests built to answer the second-round critique of
the response letter and its twelve tests.

---

## 0. What this document is

The second-round critique made one structural charge and several specific ones.
The structural charge was that our twelve tests had been graded by their
authors, and that the design on which the headline result rests — matching a
target zone against a land-cover-matched reference — had been *asserted* rather
than *tested*. Phase L was built to answer that charge first and everything else
second: **L1 is a gate**, and the remaining tests were only to be read if the
gate passed.

It passed. What follows reports all nine outcomes, including two where the
result is weaker than the first run suggested, one that is still open, and
three places where **we conceded too much in round 1** and the data says so.

Every number below is produced by code in `src/insar_wetlands/referee.py`,
exercised by `tests/test_referee.py`, and reproduced by
`notebooks/05_robustness/phaseL_gate.ipynb`.

---

## 1. The gate — L1 *(answers R1 M5, corrected)*

**The objection.** The seasonal amplitude is a difference between the mat and a
reference zone. If the reference is not genuinely comparable, the difference
measures the mismatch, not the mat. Round 1 asserted the match; it never tested
whether a *differently chosen* matched reference would give the same answer.

**The test.** Build 25 independent control patches drawn from the
land-cover-matched pool — same WorldCover class, same flatness and water
screens, outside the mat, disjoint from the published reference C — and
recompute the seasonal amplitude against each.

**Result.**

| Quantity | Phase K (unmatched controls) | Phase L (matched controls) |
|---|---|---|
| Median amplitude | — | **3.115 mm** |
| Spread | range 2.52 – 10.16 mm, sd 2.92 | **IQR 0.482** |
| Controls clearing the null p95 (2.025 mm) | — | **25 / 25 (100 %)** |
| Empirical *p* | — | **0.0033** |

**Verdict: GATE PASSED.** The amplitude is stable across independently chosen
matched references — the spread collapses from sd 2.92 mm to an IQR of 0.48 mm
once the control pool is matched on land cover. The design is now tested rather
than asserted, which was the critique's central demand.

**A diagnostic worth keeping.** The matched null (0.855 mm) agrees closely with
the published D-based null (0.850 mm), but a null built against the *real* zone
C gives 1.533 mm. The extra variance is therefore a property of **C's
fragmentation**, not of the null construction. This is a real limitation of the
published reference and should be stated as such.

---

## 2. L2 — Erode the mat *(the untested high-value test)*

**The objection.** Zone A's border pixels contain surrounding terrain through
the resolution cell and sidelobes. If the coherence deficit and the seasonal
amplitude are contamination artefacts, both should move when the border is
stripped.

**Result.** Under successive one-pixel erosions:

| Quantity | 0 px | 1 px | 2 px | Reading |
|---|---|---|---|---|
| Coherence deficit | −0.0809 | −0.0832 | −0.0789 | **stable**, non-monotonic, 5 % swing |
| Closure dispersion (× reference) | 3.22 | 3.35 | 3.58 | rises **monotonically** |
| Seasonal amplitude | ×1.00 | — | **×1.15** | against ×1.46 predicted by noise alone |

**Verdict.** The coherence deficit is **not** a border artefact: it does not
weaken under erosion, and it swings non-monotonically by 5 %, which is
consistent with resampling noise rather than with a contamination gradient. The
closure dispersion *rises* monotonically, the opposite of what contamination
predicts. The amplitude rises by ×1.15 where pure noise predicts ×1.46 — so the
signal is partly, but not wholly, explained by the shrinking sample.

---

## 3. L3 — Subdivide the reference *(makes M8 conclusive)*

**The objection.** Round 1 subdivided only the *target*. A flat amplitude under
target subdivision is consistent with a coherent signal, but also with a
reference-side artefact that subdivision never touched.

**Result.** Amplitude ratio under subdivision, both directions:

| Subdivision | Ratio | Noise prediction |
|---|---|---|
| Target | 0.93 | rise |
| **Reference** | **1.02** | rise |

**Verdict: M8 is now conclusive.** The amplitude is flat under subdivision of
*either* side. A noise-floor artefact must rise as either patch shrinks; it does
not. This closes the objection completely and is the cleanest positive result in
Phase L.

---

## 4. L4 — Closed triplets *(R1 S7)*

**The objection.** The closure-phase test found no significant bias. If the
network contained few closed triplets, that null result would be uninformative —
an absence of power, not an absence of bias.

**Result.** **518** closed triplets. Verified two independent ways —
trace(A³)/6 and direct enumeration — which agree exactly. That is **6.28×
more** than a random network of the same node and edge count.

**Verdict.** The network is unusually well closed, not poorly closed. The
closure null result is **not** explained by a shortage of triplets. Note this
inverts the direction the objection assumed.

> **This does not yet settle §4.8.** Whether the closure result is *discriminating
> and negative* or *underpowered for the effect size in question* depends on the
> moisture-phase literature (De Zan 2014; Zwieback 2015, 2017), which we have not
> yet worked through. See §10.

---

## 5. L5 — Amplitude without the network inversion *(closes the M4 hole)*

**The objection.** The published amplitude comes from the network inversion, so
our "wrapped phase is immune to unwrapping errors" defence does not cover it.

**Result.** A seasonal amplitude fitted on **wrapped** phase, with no network
inversion at all:

| Estimator | Amplitude |
|---|---|
| Wrapped, no inversion | 2.32 mm |
| Network inversion (published) | 3.29 mm |

The difference is **1.94 σ** of the wrapped estimator, whose fit quality is
*r*² = 0.0053.

**Verdict: INCONCLUSIVE, and reported as such.** The wrapped estimator is far
too noisy to either confirm or refute the network value — an *r*² of 0.005
means it explains essentially none of the variance. 1.94 σ is not a
contradiction, but neither is it a corroboration. The M4 hole is **narrowed,
not closed.**

---

## 6. L6 — Corrected null at 5 000 draws *(gated on L1)*

**Result.** 4 614 valid draws, empirical *p* = **0.0260**, exact 95 % CI
**[0.0214, 0.0308]**.

**Verdict.** The interval no longer straddles 0.05. In round 1 the draw count
was low enough that the confidence interval on the *p*-value itself crossed the
threshold, which we flagged. It no longer does.

---

## 7. L7 — The 27 pixels *(R1 M6)* — **OPEN**

This test went through four corrections. Each is recorded because the sequence
is itself the most instructive part of Phase L.

**7.1 — The original (round 1).** Mann–Whitney on 27 pixels, *p* = 0.004.
Invalid: the pixels were selected by a threshold on temporal coherence and
scored against a spatially autocorrelated distance field.

**7.2 — Toroidal permutation, two-sided: *p* = 0.262.** We initially read this
as "the concentration does not survive". **That reading was wrong.** Signed
distance is negative inside the zone, so *closer to the margin* means a
**higher** value — and the observed −40.0 m sits **above** the null p95 of
−56.6 m. An observation beyond the null's 95th percentile cannot have *p* =
0.262 for a directional question. The statistic was two-sided on |deviation
from the null median|, and the null is strongly left-skewed
(|p05 − median| = 166.7 m against |observed − median| = 49.4 m), so the long
inward tail — draws deviating in the **opposite** direction to the hypothesis —
inflated the count.

**7.3 — One-sided: *p* pinned at the floor.** `p_greater` = 1/4001 exactly:
**zero of 4 000** draws reached −40 m. The permutation therefore *bounds* p
from above and cannot estimate it. Correct statement: **p < 2.5 × 10⁻⁴
(0/4000)** — never "p = 0.00025".

**7.4 — The shape confound, and what it revealed.** The null had fallen back to
compact blobs because no rigid shift of the observed set fits inside zone A.
Shape descriptors show why that control is invalid:

| Descriptor | Observed | Compact null |
|---|---|---|
| Pixels | 27 | 27 |
| **Connected components** | **15** | 1 |
| Radius of gyration | 11.81 | 2.15 (**×5.50**) |
| Fill fraction | 0.029 | 0.643 |
| Bounding box | 25 × 37 px (1000 × 1480 m) | — |

A set of 15 scattered fragments can reach parts of a distance field that a
single 27-pixel blob cannot occupy at *any* position. The *p*-value therefore
measured **shape as much as position**, which is exactly how it came to be
pinned at the floor.

**7.5 — The correction to our own concession.** The M6 objection was that 27
*clustered* pixels at a 160 m correlation length are not 27 independent
observations. **That premise is wrong.** The 27 pixels are 15 separate fragments
— mostly isolated pixels and pairs — spread over a 1000 × 1480 m box. Most
fragments are separated by *more* than the correlation length, so the effective
sample size is nearer **15 than 1**: the Mann–Whitney was inflated by roughly a
factor of two, not by a factor of 27. It was still not a valid test, but
conceding a 27-fold inflation was an overcorrection on our part.

**7.6 — Current status: OPEN, pending re-run.** A component-matched null is now
implemented — it reproduces the observed fragment count and fragment sizes at
random disjoint positions, preserving what fixes the effective sample size while
randomising position. **It has not yet been run on the real data**; the last
Colab execution predates it. Until it runs, the honest statement is:

> The Mann–Whitney was invalid, though less so than we conceded. The marginal
> concentration is in the predicted direction and is large (−40 m against a null
> median of −89 m), but no null matched on both size and shape has yet been
> constructed, so it stands as an **observation**, not as a significant
> permutation result.

**A limitation no statistic can remove.** The observed median distance is 40 m —
exactly one pixel. An anchored margin that holds phase, and border pixels
contaminated by stable surrounding terrain, predict the same thing at this
resolution. L7 cannot separate them whatever *p*-value it returns.

---

## 8. L8 — Three intervals, not one *(R1 S4, refined)*

**The objection.** An i.i.d. bootstrap over dates destroys the temporal sampling
structure and will understate the interval.

**Result.**

| Method | 95 % interval (LOS, mm) | Width |
|---|---|---|
| Analytic | [2.22, 4.35] | 2.13 |
| Bootstrap | [0.58, 7.32] | 6.74 (**×3.2**) |

Propagating the analytic upper limit: 4.35 mm LOS → **5.14 mm vertical**
(incidence 32.3°, factor 1.183) → **10.3 mm** at an assumed coupling fraction
*f* = 0.5.

**Verdict.** The two intervals differ by a factor of 3.2, which is itself the
finding: the bound is far more method-dependent than round 1 acknowledged.

> **Decision required before submission.** Which interval the manuscript quotes
> must be stated explicitly and justified. Quoting the analytic interval while
> the bootstrap is 3.2× wider needs a reason in the text; the conservative
> choice for an upper *bound* is the wider one. This is an authors' call, not a
> computational result, and it is currently unstated.

---

## 9. Claims withdrawn from the manuscript

All instances removed and verified absent from `docs/paper/`:

- "one to two orders of magnitude too small for flotation"
- "≤ 3.9 mm" and "< 2.4 mm" bounds
- "A ≈ noise, C ≈ signal"
- "B sits at the floor"
- "below 4 mm"

These rested on the pre-Phase-L bound and on the N_eff assumptions since
measured. The surviving quantitative claim is the ≤ 8.7 mm ceiling on apparent
**phase-centre** displacement, with the explicit statement that it does not
transfer to the peat surface.

---

## 10. Outstanding — and honestly classified

### 10.1 Not deferrable, not compute-bound

| Item | Why it blocks | Effort |
|---|---|---|
| **De Zan (2014); Zwieback (2015, 2017)** moisture-phase literature | Decides whether the closure result is *discriminating and negative* or *underpowered*. §4.8 cannot be written without it. L4 establishes the network is well closed but not what effect size it could detect. | Reading, days |
| **EGMS sign check** | H4's first inferential link is unverified by our own admission. | An afternoon |
| **L7 re-run** with the component-matched null | The only Phase L result still open. | Minutes |

### 10.2 Requires field access or archive

- In-situ water-table depth series, 2022–2024
- Stratigraphic / coring records constraining the floating fraction
- Laser and UAV campaigns (see Appendix A.11 for what each would test)
- Descending-track reproduction (Mo2)
- Quantitative outcomes for all six estimators (Mo1)

### 10.3 Editorial, before submission

- Cut to ≈ 8 figures and 6 tables
- Abstract under 350 words; highlights under 85 characters each
- Complete the reference list; mint the Zenodo DOI
- Replace author and affiliation placeholders in `response_to_referees.md`
  (`[A. Sayoud]`, `[supervisor]`, `[co-authors]`)

---

## 11. Summary table

| Test | Objection | Outcome | Status |
|---|---|---|---|
| **L1** | Matched design asserted, not tested | 25/25 controls clear the null; IQR 0.482; *p* = 0.0033 | **PASSED** |
| **L2** | Border contamination | Deficit stable under erosion; closure dispersion rises | **Answered** |
| **L3** | Only the target was subdivided | Flat both ways (0.93 / 1.02) | **Conclusive** |
| **L4** | Closure null may lack power | 518 triplets, 6.28× a random network | **Answered**, §4.8 still open |
| **L5** | Amplitude depends on the inversion | 2.32 vs 3.29 mm, 1.94 σ, *r*² = 0.005 | **Inconclusive** |
| **L6** | *p*-value CI straddles 0.05 | *p* = 0.0260, CI [0.0214, 0.0308] | **Resolved** |
| **L7** | 27 clustered pixels not independent | Premise wrong (15 fragments); fair null pending | **OPEN** |
| **L8** | Bootstrap understates the interval | Analytic vs bootstrap differ ×3.2 | **Answered**, choice unstated |

**Three results changed direction once tested properly**: L4 (the network is
better closed than a random one, not worse), L7.5 (the pixels are not clustered,
so we over-conceded), and L7.2 (the concentration is in the predicted direction,
not absent). **Two results are weaker than the first pass suggested**: L5 is
inconclusive, and L7 has no valid *p*-value yet.

---

## 12. Methodological note for the referees

Four of the corrections in Phase L were to **our own** tests, found by our own
diagnostics, after the tests had already produced a publishable-looking number:

1. A two-sided statistic answering a directional question on a skewed null
   (L7.2) — it reversed the verdict.
2. A *p*-value pinned at the permutation floor, reported as if estimated
   (L7.3).
3. A null that did not match the observation on the property doing the work
   (L7.4).
4. A first version of the fairness check that used whole-set spread, which for
   a fragmented set is a function of *position* — the hypothesis under test — and
   would have suppressed a true result. Caught by a synthetic case before it
   touched the data.

We report these rather than quietly fixing them, for the same reason Appendix
A.13 records the errors corrected during the original analysis: a protocol that
never surprises its authors is not being tested.
