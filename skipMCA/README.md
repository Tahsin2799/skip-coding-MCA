# skipMCA

Questionnaire skips in multiple correspondence analysis (MCA) indices.

A block of survey items asked only of people who pass one filter question
(advice asked only of the diagnosed, child items asked only of recent mothers)
shares a pattern: who was asked. MCA builds its first axis from the strongest
shared pattern, so the coding of the skip decides how much of the filter
enters an index, and whether the filtered items dominate it.

```r
library(skipMCA)
dev <- simulate_filtered_survey(contact = 0.5, seed = 1)
items <- paste0("advice", 1:4)
inputs <- paste0("x", 1:6)

skip_audit(dev, items, inputs, budget = 0.1)          # before building the index
fit <- skip_mca(dev, items, inputs, coding = "fuzzy_average", orient = dev$need)
target <- simulate_filtered_survey(contact = 0.2, seed = 2)
targeting_metrics(predict(fit, target)[, 1], target$need, target$cluster)
```

| Function | What it does |
|---|---|
| `skip_audit()` | Each coding's squared correlation ratio with eligibility (closed form and empirical), the inertia ratio, the R² of the index with the filtered block and with eligibility, budget check |
| `skip_mca()`, `predict()` | MCA index under missing fuzzy average, subset (specific) MCA, missing single (the skip as its own category), the skip as "no", or the items dropped; new respondents placed by the transition formula |
| `eligibility_eta2()` | Theorem 1, the closed-form squared correlation ratio of a coded item with eligibility |
| `targeting_metrics()` | Cluster rank correlation, inclusion error, top-decile precision and the high-need share selected at a budget |
| `simulate_filtered_survey()` | Example data with a contact or a diagnosis filter |

The theory and its simulation test are in the accompanying paper, "Coding
questionnaire skips in composite health indices: a plasmode simulation study". Missing fuzzy average is an existing option (van
der Heijden and Escofier 2003; Josse et al. 2012); the package's contribution is
the audit.

Status: 0.1.0, not on CRAN. The maintainer email in `DESCRIPTION` is a
placeholder.
