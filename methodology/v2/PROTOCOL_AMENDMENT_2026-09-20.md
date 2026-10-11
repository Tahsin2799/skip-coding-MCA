# Protocol amendment, 20 September 2026

Amends `protocol.md`, frozen the same morning at commit `206e1d0` and unchanged
on disk. This document does not rewrite it. It records what the review of
`PAPER_IMPROVEMENT_PLAN_2026-09-20.md` changed, why, and under which of three
labels: a **correction** of something implemented or described wrongly, a
**scientific extension** that adds a question, or a **reporting decision** about
how results are read. Sections of `protocol.md` not listed here stand.

**This amendment was written after the 17, 18 and 19 September pilots and after
a read-only review of them. It is not a preregistration**, and the paper will
describe it as an amendment informed by pilot runs. No amended-design result
existed when it was written: the v03 code produced nothing but a 2-replicate
smoke run used to check that the arms execute, and the specification below was
fixed before any v03 estimate was read.

Schema string for amended runs: `v03-2026-09-20`. The unchanged `v2-pilot-0.1`
is retired, so no v03 output can be mistaken for a v2 one.

## 1. What the study is about

Unchanged in substance from `protocol.md` section 1, with one correction and one
sharpening.

**Correction.** `protocol.md` calls the four advice items part of "this
project's own published index". The index is not published: the manuscript built
on it was rejected by PLOS ONE (PONE-D-26-44638, 10 September 2026) and only a
source archive is deposited on Zenodo. The paper says "an index we built and
submitted" and nothing stronger. `PUBLISHED_INDEX_EVIDENCE.csv` holds the
applied evidence, and it is thin: one exact instance, ours, and broader examples
of healthcare-dependent predictors elsewhere. The manuscript states the limited
applied evidence rather than implying the practice is widespread.

**The primary decision.** Community prioritization by underlying structural
vulnerability. Individual outreach at a fixed budget and detection of
unrecognized disease are secondary decisions with their own endpoints, reported
separately. Recovering `V` and predicting `Y` are different objectives, and the
paper keeps them apart in every sentence.

## 2. Constructs

| Symbol | What it is |
|---|---|
| `V` | simulated structural vulnerability, the construct the index is meant to rank |
| `Y` | simulated current disease |
| `C` | contact with care, the access variable the design manipulates |
| `R` | recognition, that is diagnosis, possible only given contact |
| `G` | eligibility for the advice questions |
| `A*` | factual receipt of advice |
| `A` | the observed item, which is `A*` when `G = 1` and missing otherwise |

**Unrecognized disease is `Y = 1` and `R = 0`**, not `R = 0`, which also holds
everyone without disease. The v2 protocol never defined it; the v03 metric suite
measures capture of it directly.

**What an access intervention holds fixed.** `V` and `Y` are imposed before
access and do not move with it. Invariance claims are restricted to them and to
arms that read no gated item. Nothing here claims that an index is invariant in
any wider sense.

**Correction to a statement in `protocol.md` section 3.** It says "true
vulnerability does not depend on age in this generator, so a difference across
bands is the index performing differently and not the bands differing in need."
The first clause is right only in the narrow sense that age has no separate term
in the equation for `V`. The donor covariates do vary with age, and older
respondents have less schooling, so mean need genuinely differs across bands.
Age-band results are therefore read as changes across `q1_age_slope` and against
the oracle's own band profile, never as levels.

**The oracle is an `X`-information benchmark.** It knows the rung scores and the
true loadings but not the residual part of `V`, so it is not a universal upper
bound, and not a bound at all for a method with information it lacks. It is
never described as the best achievable.

**What the generator does not model.** Any health benefit of outreach or
treatment. Yield and coverage are decision-relevant outcomes; they are not
estimates of prevented events, and no sentence in the paper may imply otherwise.

## 3. Mechanisms: what changed in the equations

### 3.1 Rung spacing is canonical (correction)

`rung_matrix()` scored a category by its position among the levels **present in
the pool**, after `droplevels()`. The charcoal rung is empty in all four donors,
so kerosene-or-coal scored .50 where the canonical four-rung ladder puts it at
1/3, and a pool missing any rung silently rescaled every other rung in it. Rung
scores now come from the canonical ladders in `harmonize_donors.R`, so `V` means
the same thing in every population. This changes `V` slightly everywhere, which
is why every intercept was re-solved.

### 3.2 `V` is standardized on the pool (correction)

The signal was standardized inside each draw, so the 1,000-person development
half and the 20,000-person target stood on two scales that differed by the
sampling error of about fifty clusters. It is now standardized with constants
computed once from the whole pool for that specification.

### 3.3 The access settings share random numbers (correction, restoring v1)

Contact, recognition and advice were drawn with `rbinom` separately in each
setting, so a person's movement across the access grid mixed the access change
with a fresh coin. Each person now carries fixed uniforms `U_C`, `U_R` and
`U_Ak`, drawn once with the population, and

    C_s = 1{U_C < P(C | s)},  R_s = C_s * 1{U_R < p_recog},  A*_k = C_s * 1{U_Ak < p_k}

so contact and recognition are nested across c20, c50 and c80 and any migration
is attributable to access. v1 did this; v2 lost it. A tie-break uniform and a
cluster-level tie-break uniform are drawn at the same time and see nothing else.

### 3.4 The questionnaire gate is explicit (scientific extension)

`protocol.md` has one gate, `G = C`: the items are asked of everyone who reached
care. The motivating survey has the other one. BDHS 2017-18 Q317B and Q318 ask
the advice questions only of respondents who report a previous diagnosis, and
Kibria et al (2024) analyse exactly that universe, 1,710 respondents who had
been told they had high blood pressure. Both gates are now in the design:

    contact gate:    G = C
    diagnosis gate:  G = R

and `gate` is a factor of the core grid, not a stress. Factual receipt and
observation are separate: `A*` still requires contact, and the gate decides only
whether the item is recorded. Under the diagnosis gate every observed respondent
is recognized, so the recognition coefficient cannot be estimated from observed
advice, and an imputed response for an unrecognized person is a **hypothetical
completion under a different eligibility history**, not recovered advice. That
sentence is in the manuscript's Methods.

Complete-case selection and imputation training now key on the observed item set
rather than on contact. The v2 code dropped on `contacted`, which would have
implemented the diagnosis gate wrongly.

Ordinary item nonresponse inside the eligible group is defined and deliberately
absent from the core: there is no specific hypothesis about it here, and adding
it would confound the structural skip with an ordinary missingness mechanism. It
is available as a labelled stress if a reviewer asks.

### 3.5 Advice parameters are named (housekeeping)

The advice odds ratio, the advice loading on `V` and the item mode are spec
fields rather than constants inside the generator, because three of the new
controls move them.

## 4. Parameter-to-source audit

For every mechanism constant: the quantity, the conditional population it refers
to, the source, and whether the value is **estimated**, **calibrated**,
**illustrative** or a **design choice**. The corrections are the point of the
table.

| Parameter | Quantity | Source and what it actually measures | Label |
|---|---|---|---|
| `OUTCOME_PREVALENCE` = .25 | `P(Y = 1)` in adults | Geldsetzer 2019 pools 17.5 percent across 44 LMICs; Bangladeshi adult women about 28.9 percent | calibrated to a plausible level, not an estimate for any population |
| `age_log_or_per_decade` = log 1.74 | odds of `Y` per decade | Abba et al 2021, adjusted multilevel model on 12 DHS countries | estimated, and population-invariant by construction |
| `sex_slope` = .1 | male-female log odds | none | design choice, carries no weight |
| `sigma2` (ICC .05) | between-cluster variance of `Y` | Ghosh et al 2023 and this project's own fits, .036 to .053 | calibrated |
| contact grid .20/.50/.80 | `P(C = 1)` | Geldsetzer 2019: 74 percent of hypertensive adults had ever had blood pressure measured, with wide regional spread | **manipulated**, a designed range bracketing observed coverage, not an estimate |
| `q1` = .50 | `P(R = 1 | C = 1, Y = 1)` | Geldsetzer 2019 awareness .392 is **population** awareness among hypertensive adults, not conditional on contact. Conditioning on contact raises it | design anchor bracketing the published figure. **Corrected**: the ledger previously read the two as the same quantity. The stress battery moves it to .20 and .80 |
| `q0` = .05 | `P(R = 1 | C = 1, Y = 0)` | self-report specificity .945 to .964 implies 3.6 to 5.5 percent false labelling | **Corrected**: a self-report false-positive rate is measured against a measured-blood-pressure standard and depends on treatment status and the outcome definition, so it is not exactly this conditional. Design choice inside a plausible range; `q0_high` = .15 is the pessimistic stress |
| `advice_or` = 16 | odds of advice, recognized against unrecognized, among the contacted | **No survey observes this contrast**: BDHS asks only the diagnosed. Kibria et al compare advice after a qualified doctor's diagnosis with advice after another provider's, which is a different comparison. **Corrected**: the previous justification in `CITATION_LEDGER.md` mis-stated it | design choice. What is anchored is the implied level given recognition, .46 to .74, against Kibria's 54 to 77 percent salt advice among the diagnosed |
| `advice_probabilities` | baseline advice among contacted unrecognized people | none directly | design choice; `items_two`, `items_eight` vary the count |
| `advice_v_loading` = .25 | direct path from `V` to advice | none | design choice; the `advice_strong` control raises it to 1.0 |
| `access_kappa` = .736 | between-cluster SD of contact, logit scale | Nepal 2022 `wbp15`, 6,099 adults in 476 clusters | **approximate calibration**: a noisy cluster-level estimate of a lifetime "ever measured" question, not a coefficient of the latent conditional model the generator uses |
| `access_rho` = -.25 | correlation of cluster access supply with cluster deprivation | Nepal 2022, -.23 Pearson and -.24 Spearman, on a four-item deprivation count of our own | approximate calibration; attenuated by noise in both cluster means, so treated as a floor. `rho_zero` and `rho_strong` bracket it |
| `loadings`, `r2`, `BMI_RISK` | construction of `V` | v1 design values; the double burden of malnutrition justifies the non-monotone BMI shape (Popkin 2020) | design choices. `r2` is a grid axis, `loadings_flat` and `bmi_monotone` are controls |

**Time windows.** Contact is "reached care in the relevant window", recognition
is "ever told", and advice is "currently receiving advice". The sources measure
three different windows, and the generator collapses them into one period. That
is a design simplification, stated in the limitations rather than smoothed over.

**Survey weights.** The generator resamples clusters with equal probability from
a donor pool and treats every synthetic person as one person. The estimand is
therefore a property of the **finite empirical pool**, not of any national
population, and no result is described as representative of a country. A
design-weighted resampling sensitivity is deferred: the DHS biomarker and
anthropometry subsamples carry their own weights whose documentation
`../BDHS_AUDIT.md` records as unresolved, and applying donor weights after
weighted cluster resampling would double-count them. Recorded as decision D2.8.

## 5. Estimands

`protocol.md` section 3 named `cluster_rho_all`, `misallocation_all` and
`capture_all`. The amendment keeps the whole-target rule and changes what the
headline quantities are, because the policy the paper talks about is an
allocation at a budget rather than a threshold.

**Primary, community prioritization:** `cluster_rho_all`, and `cl10_misalloc`,
the share of the clusters selected under an explicit 10 percent cluster budget
whose true mean vulnerability is below the median. The quantile-based
`misallocation_all` of the frozen protocol is retained as a secondary column;
its top decile is defined by a quantile and so has a size that moves with ties,
which is why the budget version leads.

**Secondary, individual outreach:** `b10_hn_capture` at the prespecified 10
percent individual budget, with `b10_hn_shortfall_per1000` as the interpretable
unit and 5 and 20 percent as the budget curve.

**Secondary, unrecognized disease:** `b10_unrec_capture` and
`unrec_capture_tier`.

**Secondary, fixed threshold:** `capture_all`, read beside `tier_flag_rate`.

**Comparator:** the apparent-performance family, always reported beside
`coverage`.

Every definition, denominator and failure rule is in `METRIC_SPECIFICATION.md`,
which is hashed into the run manifest.

## 6. Methods compared

Ten arms. `protocol.md` had seven; the changes are three additions and one
split, each with a reason.

| Arm | Information it uses | Target-time adaptation | Status |
|---|---|---|---|
| `zero` | six inputs, items with the skip at the reference level | none | the untested convention |
| `category` | six inputs, items with a "not asked" level | none | missing-indicator analogue |
| `complete_case` | six inputs and items, eligible people only | none | drops the ineligible |
| `drop_items` | six inputs | none | discards the gated block |
| `impute_frozen` | six inputs, items completed by **development** item models | none | **new**: transport without refitting |
| `impute_adaptive` | six inputs, items completed by models refit on the **target's** eligible respondents | refits item models | the v2 `impute` arm, relabelled |
| `passive_mca` | six inputs, items with the skip category excluded from the analysis | none | **new**: subset MCA (Greenacre and Pardo 2006), the remedy the MCA literature offers |
| `simple_score` | six inputs, as a count of externally defined deprivations | none | **new**: a deployable rule with no fitting |
| `raw_items` | six inputs and items as separate covariates, supervised | none | reference |
| `oracle` | generator rung scores and true loadings | none | reference, `X`-information only |

Points the plan required to be explicit.

- **The imputation arms use age and sex** in the item models although the MCA
  input list does not. That is deliberate and is described as a deployment
  package: age and sex are available wherever the index would be built. It is
  recorded as an information advantage over the arms that do not use them.
- **The imputation is improper.** Responses are drawn from fixed fitted
  probabilities; model parameters are not drawn. No Rubin variance is claimed,
  and `m > 1` is not called proper multiple imputation.
- **The deployable rule for an imputation arm is one averaged score per
  person**, not an average of per-imputation performance. Each imputation gives
  its own MCA axis, oriented and standardized on its own development
  coordinates; the person's score is the mean of those standardized scores, and
  tiers, budgets and every metric are computed once from it. v2 averaged the
  metrics instead, which is the performance of no decision anyone could take.
- **Items are imputed one at a time**, so their dependence given the covariates
  is not preserved. The paper says so rather than claiming imputation repairs
  structural skips in general.
- `impute_m20` asks whether five draws are enough for the averaged score to have
  settled. It cannot resolve identification.
- **`simple_score` is a count of six binary deprivations**, each defined
  outside this study: global MPI cut-offs for water, sanitation, cooking fuel
  and schooling, an occupation-class convention for work, and both WHO BMI tails
  for nutrition. It uses none of the generator's rung spacing or loadings. Its
  BMI rule is tested against a truth it was not written for by the
  `bmi_monotone` control.
- **`raw_items` folds are whole clusters.** Every target cluster is unseen at
  deployment, so a person-level fold would flatter the out-of-fold development
  score.
- **`oracle` carries development scaling to the target.** It previously
  standardized target scores on the target itself, which preserves ranks but
  moves people across fixed development cuts. That is a different policy, and
  the fixed-threshold policy needs the development constants.
- **Regularized iterative MCA (`missMDA`) is deferred**, not rejected. The
  package is not installed, and Josse and Husson (2016) recommend the passive
  treatment precisely when the missing entry is a category in itself, which a
  structural skip is. Installing a package is not a justification for its
  assumptions. Recorded as decision D4.7.

## 7. Design

52 configurations, one directory, every factor in the id.

| Block | Grid | Configurations |
|---|---|---:|
| Core | {Q3, bd_allage} x r2 {.4, .6, .8} x slope {0, .5} x gate {contact, diagnosis} | 24 |
| Stress | eleven one-at-a-time variants at the reference cell under the contact gate, plus `q1_low` and `q1_high` mirrored under the diagnosis gate, where `q1` sets how many people are eligible at all | 13 |
| Controls | `recog_flat`, `proxy_block`, `advice_strong`, `bmi_monotone` | 4 |
| Quadrant | Q1, Q2 and Q4 at the reference cell | 3 |
| Donor | four donors at 18-49 x both gates | 8 |

The Q3 reference cell is the core cell itself and is **not** run a second time
for the quadrant block; the quadrant comparison reuses it. Stress conclusions
are deliberately restricted to the contact gate except for the two recognition
stresses, and the paper says so.

The four controls are there to make the conclusion falsifiable:

- `recog_flat` removes the `Y -> R` link by recognizing cases and non-cases at
  the same rate, so the items carry no disease signal beyond `V`. The direct
  `V -> advice` and access paths remain, and the paper says that too.
- `proxy_block` replaces the advice items with four always-observed risk proxies
  of the same count, functions of `V` only. It separates "more items" from
  "items behind a care gate". This is v1's `stress_baseline_proxy_4` idea
  rebuilt for the plasmode, not the v1 generator copied over.
- `advice_strong` makes the gated block carry much more information about need
  than the six inputs do. If dropping the items still wins there, the
  recommendation is stronger; if it loses, that is the boundary condition and it
  is reported as one.
- `bmi_monotone` replaces the non-monotone BMI truth with a monotone one.

Deferred, with the claim narrowed accordingly: a multidimensional-truth control
(so no claim is made that one-dimensional measurement is robust), a control
separating substantive item count from repeated gate indicators (so duplication
is described as an explanation, not a demonstrated mechanism), leave-one-donor-out
transport, Nigeria as a fifth donor, an ordinal-method tournament, and any
machine-learning comparison.

**Randomness.** One substream per purpose: the population, each arm family's
development half, and each arm-by-setting scoring step. Adding an arm cannot
move the population or another arm's draws. The stream key leaves the gate out,
so the two gates of a cell see the same people and the same uniforms and the
gate contrast is paired. Pilot and production carry different stream labels, so
no pilot replicate can be appended to production.

**Provenance.** Every run directory holds a manifest that hashes the scientific
sources, this amendment, `METRIC_SPECIFICATION.md`, the grid with every
configuration's resolved specification, the arm list, the donor cache, the
frozen quadrant membership, the package versions and the RNG policy. A resume is
refused unless the manifest signature and the configuration hash both match, and
the refusal names the component that moved. Readers take their file list from
the manifest, not from a glob. One writer per directory, through a lock that
records its owner and is reclaimed only when that owner is provably gone, so an
interrupted run resumes without a manual unlock while two live writers are still
refused. `RUN_MANIFEST_SPEC.md` section 6 states the rule.

## 8. The claim rule

Unchanged in spirit from `protocol.md` section 6 and tightened in two places.

A contrast between two arms is paired within replicate and computed **inside one
scenario cell**, a cell being a configuration by access setting. For paired
differences `d_r`, MCSE is `sd(d_r) / sqrt(R)`.

- **Robust superiority** requires the same arm to lead by more than three MCSE
  in **every** required cell. Required cells are the 24 core configurations by
  three settings.
- A **sign change** across cells is reported as scenario dependence, naming the
  factor it depends on.
- **Insufficient precision** is reported as unresolved. Unresolved is not
  equivalence.
- **Equivalence** may be claimed only against a margin fixed in advance, and
  only when the three-MCSE interval lies inside the margin in every required
  cell. Proposed margins are .02 for `cluster_rho_all` and `cl10_misalloc` and
  .01 for the capture metrics, which is ten people per thousand high-need
  people. **They are proposals pending the author's confirmation** (D3.5); until
  confirmed, the summary reports the range of differences and claims no
  equivalence.
- Exact ties, all-missing contrasts and incomplete pairs are handled explicitly.
  A zero difference with zero MCSE is a tie, not a division by zero, and the
  number of retained and dropped pairs is reported.
- Validity is metric-specific: a separated tier regression removes the odds
  ratio, not the targeting metrics computed on the same rows.

v2's summary compared only adjacent arms in a hardcoded order, at c50 only, and
oriented the misallocation family by a column name that no longer existed. All
three are fixed, and the fixes have tests whose expected answers are known in
advance.

## 9. Replicate count

`protocol.md` fixed 1,000 from the 17 September pilot. That count does not carry
over: the metrics, the arms and the gates have changed. The amended design runs
a pilot first, and the production replicate count is set from the observed
paired spread as `ceiling((sd_pilot_difference / target_mcse)^2)` for each
required contrast, with target MCSEs recorded before production in
`REVISION_PILOT_REPORT.md`. Pilot replicates are excluded from the production
stream by construction.

## 10. What this amendment freezes

Sections 2 to 8 above, together with `METRIC_SPECIFICATION.md`. Not frozen: the
replicate count, which the pilot sets; the equivalence margins, which need the
author's confirmation; and whether the age bands are a headline or a secondary
analysis, which stays a reporting decision.
