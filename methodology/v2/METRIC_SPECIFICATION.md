# Metric specification, v03

20 September 2026. Every column a v03 run writes: what it is, whose it is, what
it divides by, and what it does when it cannot be computed. Written before the
implementation it describes and hashed into every run manifest.

Two rules run through all of it.

**The denominator is the population the index is applied to.** A person the arm
cannot score is a person it does not flag, and a cluster it cannot score cannot
be prioritized. Metrics carrying `_all`, and every budget metric, divide by the
whole target. Metrics carrying `_scored` are what a validation study restricted
to scoreable respondents would report; they are diagnostics, and the gap between
the two families is itself a result.

**Thresholds come from the whole target, once.** The high-need set is the top
fifth of true vulnerability over all target people, computed before any arm is
scored, and never recomputed inside an arm's scoreable subset.

## 1. Units and populations

| Symbol | Meaning |
|---|---|
| `N` | target people, 20,000 unless the configuration says otherwise |
| `n_scored` | target people the arm gives a usable score |
| `n_hn` | target people in the top fifth of `V`, the high-need set |
| `n_unrec` | target people with `Y = 1` and `R = 0`, unrecognized disease |
| `n_clusters` | resampled target clusters |

`n_unrec` depends on the access setting, because recognition does. It is stored
per row so that any rate built on it can be recomputed.

## 2. Coverage

| Column | Definition |
|---|---|
| `coverage` | `n_scored / N`, measured against the target the arm was handed |
| `cluster_score_coverage` | share of target clusters with at least one scored resident |
| `cluster_resident_coverage` | mean over clusters of the share of residents scored |

Coverage is not a diagnostic that should read 1 for every arm. `complete_case`
scores only the eligible, and the MCA arms drop the few target rows carrying a
category development never saw.

## 3. Fixed-threshold policy: development tiers

The development score is cut at its own quintiles and oriented so the top tier
carries the most outcome. Target scores are clamped into the development range
and cut at the same breaks, so the flagged share in the target is whatever the
development cuts produce. A score with heavy ties supports fewer than five
tiers; `n_tiers` says how many, and the top tier is then wider than a fifth.

| Column | Definition | Direction |
|---|---|---|
| `n_tiers` | tiers the development score supports, 1 to 5 | context |
| `tier_flag_rate` | share of the whole target in the top tier | context |
| `capture_all` | flagged high-need people / `n_hn` | higher better |
| `unrec_capture_tier` | flagged unrecognized cases / `n_unrec` | higher better |
| `capture_scored` | scored-only version: the top fifth of `V` recomputed inside the scored rows | diagnostic |

`capture_all` and `tier_flag_rate` must be read together: an arm that flags a
quarter of the population captures more of the high-need than one that flags a
tenth, without being better at ranking.

## 4. Fixed-budget policies

Budgets are 5, 10 and 20 percent of target people. **The prespecified primary
individual budget is 10 percent.** Slots are `floor(b * N)`, an integer rule
recorded here rather than left to rounding.

Selection ranks scoreable people by score, descending, and breaks ties on a
uniform drawn once per person in the generator, before any arm exists. The
tie-break sees no truth, no outcome, no recognition, no age and no arm. Under a
pure complete-case policy only scoreable people can be selected: if fewer can be
scored than the budget allows, the remaining slots stay **unfilled** and are
reported. Filling them from a fallback score is a different policy and would
have to be a separate arm.

For each budget `b` in `b05`, `b10`, `b20`:

| Column | Definition | Direction |
|---|---|---|
| `{b}_slots` | `floor(b * N)` | context |
| `{b}_used` | people actually selected | context |
| `{b}_hn` | selected people who are high need (count) | context |
| `{b}_unrec` | selected people with unrecognized disease (count) | context |
| `{b}_hn_capture` | `{b}_hn / n_hn` | higher better |
| `{b}_hn_ppv` | `{b}_hn / {b}_used`, NA when nothing is selected | higher better |
| `{b}_unrec_capture` | `{b}_unrec / n_unrec`, NA when `n_unrec` is 0 | higher better |
| `{b}_hn_shortfall_per1000` | `1000 * (min(slots, n_hn) - {b}_hn) / N` | lower better |

The shortfall is the interpretable unit: high-need people missed per 1,000
target people, against the ideal allocation at the same budget, which is zero by
construction. It needs no monetary assumption.

## 5. Cluster targeting

A cluster's estimated priority is the mean score of its **scored** residents.
That is a deployable estimate, not a mean over every resident, and the two
coverage columns say how thin it is. A cluster with no scored resident cannot be
chosen: it is ranked tied below every scored cluster in the correlation, and it
is never selected under the budget policy. With no scored cluster anywhere, the
cluster columns are NA rather than zero.

Cluster truth is the mean of `V` over **all** residents.

| Column | Definition | Direction |
|---|---|---|
| `cluster_rho_all` | Spearman correlation of estimated priority with cluster truth, over all target clusters | higher better |
| `cl10_slots` | `floor(0.10 * n_clusters)` | context |
| `cl10_used` | clusters selected | context |
| `cl10_misalloc` | selected clusters whose true mean `V` is below the median cluster / `cl10_used` | lower better |
| `cl10_precision` | selected clusters that are in the true top decile / `cl10_used` | higher better |
| `cl10_hn_reached` | high-need people living in selected clusters / `n_hn` | higher better |
| `cl10_residents` | residents of selected clusters | context |
| `misallocation_all` | legacy quantile version: clusters at or above the 90th percentile of estimated priority whose true mean is below the median | lower better |
| `cluster_rho_scored`, `misallocation_scored` | the same two, computed on scored rows only | diagnostic |

Equal cluster counts do not mean equal residents served or equal cost;
`cl10_residents` is reported so that is visible.

## 6. Disease model and apparent performance

These are properties of a **specified predictor among its stated evaluation
population**, not of the index alone, and they are the comparator in this study
rather than the target.

| Column | Definition |
|---|---|
| `or_top` | odds ratio of the top development tier against the bottom, from `y ~ tier + age + sex` fitted in development; NA when the score supports one tier |
| `or_top_target` | the same association refit in the target's scored rows |
| `beta_within`, `beta_between` | Mundlak within and between slopes of the standardized development score |
| `auc_model_scored` | AUC of the development-fitted `tier + age + sex` model among scored target rows |
| `auc_index_scored` | AUC of the index score alone among the same rows |
| `auc_agesex_scored` | AUC of a development-fitted age-and-sex model among the same rows, the baseline the other two are read against |
| `cal_slope`, `cal_intercept` | calibration of the development model's linear predictor among scored target rows |
| `n_dev_used` | development rows the index was built on |

`or_top` is a development quantity: development is always generated at c50, it
is built once per replicate and shared by the three access settings, so it
cannot carry a transport claim. `or_top_target` is the one that moves with
access. Whichever the manuscript reports, it says which it means.

No prediction is invented for an unscored row to make an AUC denominator equal.
Where an arm scores fewer people, its apparent performance is read on the people
it scores, and `coverage` stands beside it.

## 7. Age bands

Bands are 18-29, 30-44, 45-59 and 60 and over, fixed rather than
population-specific, so an 18-49 pool simply returns NA above 49. Bands rather
than a linear age term because discrimination and calibration are not linear in
age (Kartoun et al, *Sci Rep* 2022;12:12542).

| Column | Definition |
|---|---|
| `n_hn_<band>` | high-need people in the band, over the whole target: always reported |
| `cap_tier_<band>` | flagged high-need people in the band / `n_hn_<band>` |
| `cap_b10_<band>` | high-need people in the band selected at the 10 percent budget / `n_hn_<band>` |
| `capture_scored_<band>`, `auc_scored_<band>` | scored-only diagnostics |

The high-need threshold is the whole-target one in every band, so a band is not
re-ranked against itself. A band with fewer than 100 high-need people returns NA
for the rate and still reports the count. **Unscored people stay in the band
denominator**: that is how unequal coverage across ages becomes visible, and it
is the defect the v2 age metric had.

True vulnerability has no direct age term in the generator, but the donor
covariates do vary with age, so bands differ in mean need. Read the change
across `q1_age_slope`, and against the oracle, not the level.

## 8. Failure, validity and missing values

| Column | Meaning |
|---|---|
| `failed` | the arm produced no row for this replicate and setting |
| `stage` | `development` or `target`, where it failed |
| `message` | the condition message, kept rather than summarized |
| `separation_dev`, `separation_target` | the tier regression separated |
| `level_fixups` | unseen levels sent to the modal seen level during imputation |

Validity is metric-specific. A separated tier regression invalidates `or_top`
and the Mundlak slopes computed from it; it does not invalidate the cluster
correlation or the budget metrics computed on the same rows, and the summary
drops rows per metric through the `needs` column of its metadata table rather
than dropping the replicate.

An empty subgroup returns NA with its count, never a zero success rate. An
all-NA band returns NA rather than NaN. A rate with a zero denominator is NA.

## 9. Realized rates, stored per replicate and setting

`contact_rate`, `recognized_rate`, `eligible_rate`, `unrecognized_case_rate`,
`outcome_rate`. Eligibility is not contact under the diagnosis gate and must not
be checked against the contact grid.

## 10. Primary and secondary

**Primary, community prioritization:** `cluster_rho_all` and `cl10_misalloc`.

**Secondary, individual outreach at a budget:** `b10_hn_capture`, with
`b10_hn_shortfall_per1000` as the interpretable unit, and `b05`/`b20` as the
budget curve.

**Secondary, detection of unrecognized disease:** `b10_unrec_capture` and
`unrec_capture_tier`.

**Secondary, fixed-threshold policy:** `capture_all` with `tier_flag_rate`.

**Comparator, apparent performance:** section 6, always beside `coverage`.

**Diagnostics:** every `_scored` column, the age bands, coverage and the rates.
