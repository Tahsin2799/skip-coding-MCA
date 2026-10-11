# Coding questionnaire skips in composite health indices: code

Code for *Coding questionnaire skips in composite health indices: a plasmode
simulation study* (submitted to the Journal of Clinical Epidemiology). The paper derives how much of a questionnaire filter
each coding of a skipped item carries into an index built by multiple
correspondence analysis (MCA), and tests the theory in a plasmode simulation
that resamples real covariates from four Demographic and Health Surveys,
imposes known disease, contact with care, diagnosis and advice items, and
compares ten prespecified ways of building the index.

This release contains the code, the files the code reads, the summary results
from which every table and figure of the paper is drawn, and the R package
`skipMCA`, which implements the filter audit and missing fuzzy average for
index builders. It contains no survey data.

Author: A. F. M. Tahsin Shahriar (ORCID
[0009-0003-0565-0803](https://orcid.org/0009-0003-0565-0803)), Institute of
Applied Statistics and Data Science, University of Dhaka. Citation metadata are
in `CITATION.cff`.

## What is here

| Path | What it does |
|---|---|
| `methodology/v2/harmonize_donors.R` | reads the DHS person recodes and harmonizes the six index inputs by value label |
| `methodology/v2/build_donor_cache.R` | writes the harmonized covariate cache every run reads |
| `methodology/v2/quadrant_populations.R`, `build_quadrant_cells.R`, `quadrant_cells.csv` | the four structural strata of deprivation level against coherence, pooled from Indian, Albanian and Bangladeshi region-by-residence cells, and their frozen membership |
| `methodology/v2/plasmode_generate.R` | the data-generating mechanism |
| `methodology/v2/index_arms.R` | the ten methods and every metric |
| `methodology/v2/pilot_v2_config.R` | the design grid, the solved intercepts and the variations |
| `methodology/v2/solve_intercepts.R` | solves the outcome and contact intercepts |
| `methodology/v2/run_pilot_v2.R`, `run_manifest.R` | the runner, its manifest and its resumable lock |
| `methodology/v2/verify_from_artifacts.R` | recomputes metrics from saved person-level data in base R |
| `methodology/v2/run_checks.R`, `summarize_convergence.R`, `claim_ledger.R`, `size_production.R` | the readers: checks, summary, the prespecified decision rule, and replicate sizing |
| `methodology/v2/protocol.md`, `PROTOCOL_AMENDMENT_2026-09-20.md`, `METRIC_SPECIFICATION.md` | the frozen protocol, its amendment and the metric definitions; the runner hashes them into every manifest |
| `methodology/v2/tests/` | contract tests |
| `JCE/theory/filter_capture.R` | the theory's quantities: squared correlation ratios with eligibility, the inertia ratio, block and eligibility R², and missing fuzzy average (post hoc) |
| `JCE/theory/theory_production_link.R`, `theory_link_summary.R` | regenerates every production replicate from its seeds, computes the theory diagnostics and scores missing fuzzy average, then pairs it with the ten methods under the frozen rule (post hoc) |
| `JCE/theory/theory_checks.R`, `josse2012_check.R`, `theory_figure.R` | probe checks of the propositions, the rebuilt simulation of Josse et al. (2012), and Figure 1 |
| `JCE/extensions/` | the three additional simulations of section 4.4 (advice loading, development size, development contact), with regularized iterative MCA, Firth-penalized item models and principal component indices; `README.md` is their runbook (post hoc) |
| `JCE/analysis/review_analysis.R`, `test_review_analysis.R` | the documented reporting correction: recomputes the validation-targeting correlations with the specified separation exclusions, after checking every input against the production manifest, and its regression test |
| `JCE/analysis/editor_review_checks.R` | the mechanism and item-model checks quoted in supplement A3 and A4 |
| `JCE/figures/figures.R`, `JCE/supplement/supplement_tables.R`, `supplement_methods.md` | the manuscript's figures, supplementary tables and technical appendix text |
| `skipMCA/` | the R package |
| `outputs/methodology/production_v03_displays.rds`, `production_v03_claims.csv` | summary results of the production run: the displays the tables and figures use, and the 132 classified contrasts |
| `JCE/analysis/review_analysis.rds`, `correlation_correction.csv`, `correlation_exclusions.csv` | results of the reporting correction |
| `JCE/theory/theory_link_summary_r1000.rds` | the theory diagnostics on the 51,000 regenerated production replicates with filtered items, and missing fuzzy average under the frozen rule (post hoc) |
| `JCE/extensions/extensions_summary.rds` | summaries of the three additional simulations (post hoc) |

## Requirements

R 4.5.2 with FactoMineR 2.13, dplyr 1.2.0, tidyr 1.3.2, tibble 3.3.1, digest
0.6.39, cli 3.6.5 and haven 2.5.5, the versions the production run recorded in
its manifest. matrixStats 1.5.0 is needed only to rebuild the quadrant
membership, MASS for the rebuilt Josse et al. simulation, and ggplot2 4.0.2 and
patchwork 1.3.2 for the figures. The additional simulations also need logistf
1.26.1 and testthat 3.3.2. Regularized iterative MCA needs missMDA 1.23, which
requires FactoMineR 2.17; FactoMineR 2.17 moves the production MCA methods in
the fifth or sixth decimal place, so keep 2.13 in the main library and install
missMDA with FactoMineR 2.17 in a separate library that only that pass puts
first on its path (`--lib`, below). The runs use six PSOCK workers; do not
switch to forked workers (`mclapply`), which crash inside `La.svd` when MCA
runs in them. Limit BLAS threads before R starts, as the production commands
below do.

## Data

The DHS files are available from the DHS Program (https://dhsprogram.com) after
registration and are not redistributed here. Four person recode (PR) files are
needed, plus the women's and men's recodes (IR, MR) for occupation in India,
Nepal and Albania:

| Survey | Files, as the folders are named on download |
|---|---|
| India NFHS-5 2019-21 | `IAPR7EDT/IAPR7EFL.DTA`, `IAIR7EDT/IAIR7EFL.DTA`, `IAMR7EDT/IAMR7EFL.DTA` |
| Nepal 2022 | `NPPR82DT/NPPR82FL.DTA`, `NPIR82DT/NPIR82FL.DTA`, `NPMR82DT/NPMR82FL.DTA` |
| Albania 2017-18 | `ALPR71DT/ALPR71FL.DTA`, `ALIR71DT/ALIR71FL.DTA`, `ALMR71DT/ALMR71FL.DTA` |
| Bangladesh 2017-18 | `BDPR7RDT/BDPR7RFL.DTA` |

Put the folders in one directory and build the covariate cache:

```sh
Rscript methodology/v2/build_donor_cache.R --dhs /path/to/that/directory
```

It prints the cache's sha256, which should be
`c9ac7c45a02be2563164ed6187353f8f16557f855a732ac253e9ad9ba4605e11`, the value
the production run hashed. `harmonize_donors.R` still names the author's data
directory on its line 16; it is left as it was run, because it is one of the
hashed sources, and `--dhs` overrides it.

## Running it

All commands run from the repository root.

**1. Contract tests**, which need the cache:

```sh
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 Rscript methodology/v2/tests/run_all.R
```

**2. Production**, 52 configurations at 1,000 replicates, about 24 hours on
six workers of an Apple M2:

```sh
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 \
  Rscript methodology/v2/run_pilot_v2.R --replicates 1000 --workers 6 \
    --stream production --blocks core,stress,controls,quadrant,donor \
    --artifact-reps 2 --out outputs/methodology/production_v03
```

The run is resumable: issue the same command again after an interruption and it
continues from the last finished replicate. It refuses to resume if any hashed
source, the grid, the packages or the donor cache has changed, and names what
changed. Verify it before reading any estimate, then build the claim table and
the displays:

```sh
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 Rscript methodology/v2/verify_from_artifacts.R outputs/methodology/production_v03
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 Rscript methodology/v2/run_checks.R outputs/methodology/production_v03
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 Rscript methodology/v2/summarize_convergence.R outputs/methodology/production_v03 --margins
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 Rscript methodology/v2/claim_ledger.R outputs/methodology/production_v03 --margins
```

**3. The theory on the production replicates** (post hoc), about 3 hours. The
first command checks that regeneration reproduces production; the second
regenerates every replicate from its seeds:

```sh
Rscript JCE/theory/theory_production_link.R --verify
Rscript JCE/theory/theory_production_link.R --reps 1000 \
  --blocks core,stress,controls,quadrant,donor --workers 6 \
  --out outputs/methodology/theory_link_r1000
Rscript JCE/theory/theory_link_summary.R outputs/methodology/theory_link_r1000
Rscript JCE/theory/theory_checks.R
Rscript JCE/theory/josse2012_check.R
Rscript JCE/analysis/editor_review_checks.R
```

**4. The additional simulations** (post hoc), about 25 hours in two passes; the
first two commands test the code and check that the reference configuration
reproduces production exactly. `JCE/extensions/README.md` gives the details.

```sh
Rscript JCE/extensions/test_extensions.R
Rscript JCE/extensions/run_extensions.R --verify
Rscript JCE/extensions/run_extensions.R --reps 1000 \
  --studies sweep,devsize,devcontact --workers 6 --out outputs/methodology/extensions_v1
Rscript JCE/extensions/run_extensions.R --reps 1000 \
  --studies sweep,devsize,devcontact --arms rimca --lib /path/to/missMDA-library --workers 6 \
  --out outputs/methodology/extensions_v1
Rscript JCE/extensions/extensions_summary.R outputs/methodology/extensions_v1
```

**Tables and figures from the included results.** The summary results are
included, so the supplementary tables and the figures can be rebuilt without
the survey data or any simulation run:

```sh
Rscript JCE/supplement/supplement_tables.R
Rscript JCE/figures/figures.R
```

**5. Tables, figures and appendix from your own runs.** The reporting correction reads the
production run only; it regenerates no replicate and changes neither the claim
table nor the rule that produced it. The figures and supplement read its
output, so it runs first:

```sh
Rscript JCE/analysis/test_review_analysis.R
Rscript JCE/analysis/review_analysis.R
Rscript JCE/supplement/supplement_tables.R
Rscript JCE/figures/figures.R
```

**The package.** `skipMCA` needs base R only; FactoMineR, knitr, rmarkdown and
testthat are suggested for its checks and vignette:

```sh
R CMD build skipMCA && R CMD check --no-manual skipMCA_0.1.0.tar.gz
```

## What the runs produced

Production: 52 configurations, 1,000 replicates each, 1,560,000 method-level
rows and no failed fit, under manifest signature `572bb27bcdda`. 31,200 metrics
recomputed from saved person-level data agreed with the run exactly. The
prespecified decision rule, with equivalence margins of .02, .02 and .01,
ordered 83 of the 132 contrasts consistently, found 38 orderings that depend on
the scenario and one equivalence, and left 10 undecided.

Theory (post hoc): regeneration reproduced production exactly in all 52
configurations; Theorem 1 matched the empirical squared correlation ratio in
all 204,000 coding-by-replicate rows to within 6e-15.

Additional simulations (post hoc): 16 configurations at 1,000 replicates in two
passes with no failed chunk. The reference configuration reproduced production
exactly over all 1,000 replicates under both filters (60,000 rows, 81 metrics).
18 of 822,000 method rows failed: both Firth-penalized imputers in 3 replicates
where the development sample held about 30 diagnosed respondents.

Random numbers are seeded per purpose from the stream label and the
configuration, so a rerun with the same code, packages and cache draws the
same random numbers; floating-point results can still differ slightly across
platforms.

## License

MIT, see `LICENSE`.
