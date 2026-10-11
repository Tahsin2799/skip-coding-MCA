# v2 protocol: how the structural skip in an advice-based index should be coded

20 September 2026. This is the analysis protocol for the final simulation. It is
written to be frozen before the production run, and every mechanism is stated as
an equation here before it appears as code. Where a number is sourced, the
source is in `CITATION_LEDGER.md`; where it is a design choice, it says so.

No EQUATOR checklist covers simulation studies, so the reporting follows ADEMP
(Morris, White and Crowther, *Stat Med* 2019;38:2074-2102,
doi:10.1002/sim.8086), and the paper will say that rather than claim a checklist
it does not have.

## 1. Aim

An index built by multiple correspondence analysis on survey items is a common
way to rank people or places by need. Four of the ten inputs in this project's
own published index are clinical advice items, which are asked only of people
who reached care. Everyone else is a structural skip, and v1 and
`R/analysis_main.R:79` both send that skip to the reference level without ever
testing the choice.

The question is not whether the index discriminates. It is whether the coding of
the skip decides **who the index finds**, and whether the performance statistics
a validation study would report can detect the difference. The answer the pilot
gives is that they cannot, and that is the paper.

## 2. Data-generating mechanism

Covariates are resampled from real DHS records and never modelled, which is what
makes this plasmode rather than parametric. Everything downstream of the
covariates is imposed from a known conditional. Access is simulated and never
resampled: resampling a treatment-like variable alongside the covariates breaks
positivity (Shaw et al 2025, arXiv:2504.11740).

### 2.1 The graph

```
        (household and person, resampled)
                 X = water, sanitation, fuel, bmi, education, work
                 age, sex, cluster c
                          |
                          v
   r2 ------------->      V  (true vulnerability)  <---- e, individual noise
                        / | \
                       /  |  \
                      /   |   \
      Vbar_c --> A_c /    |    \                  u_c ~ N(0, sigma2)
    (cluster access)      |     \                      |
            \             |      \                     |
             \            v       ------------------>  v
              ------->   C  (contact) <--- eps        Y  (outcome)
                          |                            |
                          v                            |
                          R  (recognized) <------------
                          |
                          v
                 A1..A4  (advice items, NA unless C = 1)
                          |
                          v
                     index arm ----> tier ----> every metric
```

Two edges carry the whole study. `Y -> R` is why the advice items look
predictive, and `A_c -> C` with `Vbar_c -> A_c` is why the missingness is
cluster-structured rather than individual noise. The outcome is imposed before
access, so `Y` does not depend on the setting and an arm that uses no advice
item cannot move across the access grid. That is the negative control, and it is
asserted in the test suite rather than assumed.

### 2.2 Equations

Rung scores `r(X) in [0,1]^6`, worst = 1, education reversed. BMI is
deliberately non-monotone, `(underweight, normal, overweight, obese) =
(.70, .00, .50, 1.00)`, so an index that orders the four classes monotonically
must invert one tail.

True vulnerability, with `z()` for standardize and `lambda = (.8, .7, .6, .6, .5, .5)`:

    V_i = sqrt(r2) * z( r(X_i)' lambda ) + sqrt(1 - r2) * e_i,    e_i ~ N(0, 1)

`r2` is the share of V the six inputs carry, and it bounds what any index of
them could recover. It is a design choice, which is why it is a grid axis rather
than a constant.

Cluster terms, one draw per cluster:

    u_c ~ N(0, sigma2),        sigma2 = (.05 / .95) * pi^2 / 3      (ICC .05, sourced)
    A_c = rho * z(Vbar_c) + sqrt(1 - rho^2) * eta_c,   eta_c ~ N(0, 1),  rho = -.25

`A_c` is cluster access supply and `rho` is the inverse care law: where mean
deprivation is higher, supply is lower.

Outcome, imposed before access:

    logit P(Y_i = 1) = alpha_p + .8 * V_i + u_c + log(1.74) * (age_i - mean(age)) / 10
                       + .1 * 1{male}

`alpha_p` is solved per population so that `E[Y] = .25`; it has to absorb that
population's covariate and age distribution, and nothing above it does. The age
term is per decade, straight from the published odds ratio, so it is a property
of age and not of the population's age spread.

Contact, one equation per setting `s in {c20, c50, c80}`:

    logit P(C_i = 1) = gamma_s + .4 * eps_i + .3 * V_i + .736 * A_c

`gamma_s` is solved so that `E[C] = .20, .50, .80`. It is solved against the
full model including `A_c`: the 17 September pilot solved it without that term
and realized .218 and .780.

Recognition, given disease and given contact:

    q1_i = expit( logit(q1) + s_age * z(age_i) ),    q1 = .50
    P(R_i = 1 | C_i = 1) = q1_i if Y_i = 1, else q0 = .05
    R_i = C_i * Bernoulli(.)

`s_age` is `q1_age_slope`. Its direction is sourced (Foti and Hardy, *Am J
Hypertens* 2025;38:529-531: adults 18-39 with regular care access are diagnosed
more slowly than those 40 and over). Its magnitude is a design choice, so it is
a grid axis with 0 as its own control, and 0 reproduces the constant-q1
generator exactly.

Advice items, `k = 1..4`, `p = (.10, .10, .05, .15)`:

    logit P(A_ik = 1) = logit(p_k) + log(16) * R_i + .25 * V_i
    A_ik observed only if C_i = 1, otherwise NA

The NA is the structural skip. It is the object of the study.

### 2.3 Populations

| Population | What it is | Records | Clusters |
|---|---|---:|---:|
| `Q1`-`Q4` | India strata, crossed on deprivation level against coherence of the deprivation structure, membership frozen in `quadrant_cells.csv` | 27,644 / 24,772 / 82,643 / 28,486 | 1,336 to 4,396 |
| `bd_allage` | Bangladesh 2017-18, 18 to 95, both sexes, the only donor whose six inputs reach past 54 | 12,921 | 672 |
| `donor_*` | India, Nepal, Albania, Bangladesh at 18-49, the replication check | 167,715 / 6,510 / 13,542 / 9,508 | 9,085 / 476 / 713 / 672 |

`donor_nepal` is a women-only population: 2,415 Nepali men aged 18 to 49 have a
measured BMI and 4,273 match to the men's recode, but the intersection is
exactly zero, because the anthropometry and men's interview subsamples are
disjoint. That is recorded in `DONOR_INVENTORY_2026-09-16.md` and it is a
property of the survey, not of the harmonization. It means `sex` is constant
there, so the sex term drops out of every model in that population, and the
Nepal replication is a replication in women.

Development is 1,000 people, target is 20,000, both drawn as whole clusters with
replacement so the real within-cluster correlation of water, sanitation and fuel
survives. A cluster drawn twice becomes two clusters. Pools smaller than the
target redraw clusters: `bd_allage` about 1.5 times and `donor_nepal` about 3.1.
That is recorded rather than corrected, because holding the target size constant
is what makes the donor arm a replication check.

## 3. Estimands

The index is a ranking device, so the estimand is about the ranking, not about
the fit.

**Primary, targeting.** Over the whole target population:

- `cluster_rho_all`, the Spearman correlation between cluster mean index score
  and cluster mean true vulnerability.
- `misallocation_all`, the share of the clusters in the index's top decile whose
  true mean vulnerability is below the median.

**Secondary, early detection.** `capture_all`, and `capture_<band>` in
18-29, 30-44, 45-59 and 60+: of the people genuinely in the top fifth of true
vulnerability, the share the index places in its top tier. True vulnerability
does not depend on age in this generator, so a difference across bands is the
index performing differently and not the bands differing in need.

**Apparent performance**, which is what a validation study would report and is
the comparator rather than the target: `or_tier5_target`, `auc_target`,
`cal_slope`, `cal_intercept`, and the Mundlak within and between slopes.

### 3.1 The estimand is the population the index is applied to

Every arm is scored on the whole target, including the rows it cannot score. An
arm that cannot score a person does not flag that person, and a cluster it
cannot score cannot be prioritized.

This is not a technicality. `complete_case` drops everyone who made no contact,
so at c20 it is applied to a fifth of the population, and that fifth is enriched
for vulnerability: mean V .159 against .000. Measured only on the rows it can
score it reads `cluster_rho` .595 and `capture` .431; measured on the population
it is meant to serve it reads .374 and .100. Those are means over 50 replicates
at Q3, `r2` .6 and slope 0, with an MCSE of .007 on the second pair. The first
number is not wrong, it
answers a different question, and reporting only the first makes the comparison
between arms a comparison between different populations.

Sperrin et al (*J Clin Epidemiol* 2020;125:183-187,
doi:10.1016/j.jclinepi.2020.03.028) set the rule this follows: for prediction,
what matters is performance in the population where the model will be applied,
and missing data at deployment is part of the model rather than a nuisance to be
removed first. Sisk et al (*Stat Methods Med Res* 2023;32:1461-1477,
doi:10.1177/09622802231165001) treat "missing data permitted at deployment"
against "prohibited" as a design scenario and evaluate both; they also find that
missing indicators, which is this study's `category` arm, help in many settings
but are harmful under outcome-dependent missingness, which is exactly the
mechanism here. Hoogland et al (*Stat Med* 2020;39:3591-3607,
doi:10.1002/sim.8682) give the validation-time version of the same argument.

Both versions are stored. The `_all` columns are primary; the unsuffixed ones
are what a scoreable-rows-only validation would report, and the gap between them
is itself a result. For an arm that scores everyone they are identical by
construction, which is asserted in the test suite.

## 4. Methods compared

Five coding arms, all building the same MCA index on the same six inputs plus
the four advice items, differing only in how the skip is coded:

| Arm | The skip becomes |
|---|---|
| `zero` | the reference level. What v1 and `R/analysis_main.R:79` do, untested |
| `category` | its own level, "not asked" |
| `complete_case` | nothing: the person is dropped |
| `drop_items` | nothing: the items are dropped |
| `impute` | a draw from a model on the six inputs, m = 5 |

Two reference arms bound the comparison. `oracle` scores with the generator's
own rung scores and true loadings and is the best any index of these six inputs
could do. `raw_items` enters every input and item as its own covariate with no
composite, taking its development score out of fold since unlike MCA it has seen
the outcome; it answers "why build a composite at all".

`drop_items` and `oracle` use no advice item, so neither can depend on access.
They are the negative controls, and if either moves across the access grid the
run is not trustworthy.

## 5. Design

One directory per run. A configuration id leaves out a factor that is held at a
single value, so two runs can otherwise collide on one file and share one seed
stream; the runner now refuses to resume a file written under a different
specification, and warns when the output directory holds a configuration outside
the current grid.

| Run | Grid | Configs | Directory |
|---|---|---:|---|
| Core | {Q3, bd_allage} x r2 {.4, .6, .8} x q1_age_slope {0, .5} | 12 | `production_v02` |
| Stress | eleven one-at-a-time variants at Q3, r2 .6, slope 0 | 11 | `production_v02` |
| Quadrant robustness | Q1-Q4 at r2 .6, slope 0 | 4 | `quadrant_v02` |
| Donor replication | four donors at 18-49, r2 .6, slope 0 | 4 | `donor_v02` |

Access setting and coding arm are within-replicate: all three settings and all
seven arms share one generated population, so they multiply fits rather than
configurations. 1,000 replicates everywhere. At the measured 26.5 minutes per
1,000-replicate configuration on six workers that is about 14 hours in total.

The quadrant axis is a robustness check and not a factor: measured on the pilot
it explained .1 to 2.1 percent of variance against the coding arm's 82 to 95,
with no interaction. The donor axis answers the JCE scope rule, which rejects
work "limited to one country, jurisdiction, or clinical area with limited
generalizability"; the quadrants are built from India's own strata, so without
the donor run the design is a one-country design.

The eleven stress variants and what each one is for are in
`pilot_v2_config.R`. Four of them change a term inside the contact model or the
between-cluster spread of V and carry their own solved contact intercepts; the
test suite asserts that every variant still realizes the nominal .20 / .50 / .80
grid and the .25 prevalence, because a variant that quietly lands at .18 is a
confounded comparison rather than a stress test.

## 6. Performance measures, Monte Carlo error, and what may be claimed

Each metric is averaged over replicates and reported with its Monte Carlo
standard error.

**The development half of every replicate is common to the three access
settings.** Development is generated once, at c50, and does not depend on
access, so the coding frames, the MCA fit, the tier cuts, the tier model and the
Mundlak fit are built once per replicate and handed to all three settings. Until
19 September each setting rebuilt them, which redrew the imputation draws and
the `raw_items` fold split, and `or_tier5`, a development quantity, moved by up
to 10.9 odds ratio units across settings in 15,997 of 16,000 replicates although
the development people were identical. Morris et al give the reason this
matters: applying every method to the same simulated data is what leaves the
contrast between methods free of that Monte Carlo error.

**Comparisons between arms are paired within replicate.** The arms see the same
people, so most of the Monte Carlo error cancels in the difference. Unpaired,
the MCSE on `capture_18_29` at 1,000 replicates is .0058 and half the adjacent
pairs look tied; paired it is about .001 and they are not.

**An ordering is claimed only when it survives everywhere it is measured.** The
rule, fixed before the run: a pair is ordered if the paired contrast exceeds
three MCSE in every cell of the core grid, and is reported as a tie otherwise. A
pair whose sign changes across cells is reported as a dependence on that factor,
not as an ordering. This is the discipline Nießl et al (*WIREs Data Min Knowl
Discov* 2022;12:e1441, doi:10.1002/widm.1441) argue for: a benchmark ranking
that is not checked against the multiplicity of design options will be
over-optimistic about itself.

Measured on the 18 September pilot, this rule already removes three orderings
that an earlier run had claimed:

- `impute` against `drop_items` on targeting. The difference runs from -.009
  (t = -9.0) to +.010 (t = +14.4) depending on population, `r2` and
  `q1_age_slope`. They are equivalent, and that is the finding: repairing the
  items and discarding them do the same work.
- `impute` against `zero` on the odds ratio and on AUC. Flips sign between
  `r2` .4 and .6 at t of 3.4 and 7.7.
- `zero` against `complete_case`. In Q3 `zero` is worse by .065 to .097 at
  t = -39 to -53; in `bd_allage` it is better at r2 .4 and .6 and worse at .8.

What the pilot does support, and what the production run is expected to confirm:
`category` collapses rather than degrades; the advice-free and advice-repaired
arms lead on targeting; `complete_case` leads on apparent performance and not on
targeting; and the two families are close to uncorrelated.

Because 31 configurations cannot be tabled inside 3,000 words, the report is a
meta-model of each performance measure on the design factors, as Skrondal
recommends (*Multivar Behav Res* 2000;35:137-167,
doi:10.1207/s15327906mbr3502_1), plus a nested loop plot for the surviving grid
(Rücker and Schwarzer, *BMC Med Res Methodol* 2014;14:129,
doi:10.1186/1471-2288-14-129).

### 6.1 Replicate count

1,000, set from the worst MCSE on the 17 September pilot and confirmed against
the metrics added since. Worst at 1,000 was .0167 on `log_or` and .0136 on
`log_or_target`. The age-band metrics are noisier because they are computed on
subsets: unpaired MCSE .0058 to .0068 on `capture_<band>`. Paired, which is how
they are read, every adjacent arm pair separates at t = 10 to 106. No metric
needs more replicates.

### 6.2 Failures

A failure is never random, so the count is reported rather than dropped
silently. Two known sources are handled:

- A constant input. `glm` aborts on a one-level factor, and Bangladesh 18-95 is
  92 percent one water rung, so a 1,000-row draw can hold a single level. The
  18 September run lost three replicates this way, on `impute` and `raw_items`
  only and on the least variable draws. The constant term is now dropped from
  the fit and the index still sees the column.
- A constant adjuster. `sex` is not a variable in every population: Nepal at
  18-49 is women only, so `y ~ tier + age + sex` aborts for every arm including
  the two reference arms. A probe of the donor arm lost all 84 of its Nepal rows
  to this before it was fixed. Constant adjusters are dropped from the tier
  model, the Mundlak fit and the target refit alike.
- Separation. Flagged on v1's rule, and those rows are dropped from the summary
  rather than winsorized. Ratios are summarized on the log scale, because one
  separated fit at 1.4e8 otherwise destroys the mean.

## 7. What is frozen

Frozen by this document: the equations in section 2, the populations in 2.3, the
estimands in section 3 including the common-population rule in 3.1, the arms in
section 4, the grid in section 5, and the claim rule in section 6.

Not frozen, and deliberately so: the choice of `or_tier5` or `or_tier5_target`
for any given sentence in the paper, which must simply be stated each time; and
whether the `bd_allage` age bands are the headline or a secondary analysis,
which is a reporting decision best made with the bands in hand.

Held open: Nigeria 2024 as a fifth donor, the only candidate that would populate
the charcoal rung, which is empty in all four current donors.

After this document is frozen, `harmonize_donors.R`, `plasmode_generate.R`,
`index_arms.R`, `quadrant_populations.R`, `pilot_v2_config.R` and
`run_pilot_v2.R` may not be edited while a run is in progress. The runner
sources the first five and a resumed run would source the edited files.

## 8. Deviations from the pilot, and why

| Change | Reason |
|---|---|
| Development built once per replicate | `or_tier5` moved across settings for `impute` and `raw_items`; common random numbers |
| Targeting read on the whole target | `complete_case` was scored on a contact-selected fifth of it at c20 |
| `coverage` measured against the target the arm was handed | it was measured inside the arm's own subset and read 1.000 whatever was dropped |
| Constant inputs dropped from fits | three replicates lost, not at random |
| Constant adjusters dropped from the metric models | Nepal is women only and every arm aborted there |
| Resume refuses a spec mismatch | two runs can share one configuration id and one seed stream |
| Donor populations added | an all-India design fails the JCE scope rule |
| Stress battery added | the pilot ran the core grid only |
