# Index construction under the skip-coding arms, the deployable comparators, and
# the metrics each one is scored on. Amended for v03; see
# PROTOCOL_AMENDMENT_2026-09-20.md and METRIC_SPECIFICATION.md.
#
# The advice items are asked only of people inside the questionnaire gate.
# Everyone else is a structural skip. v1 and R/analysis_main.R:79 both send that
# skip to the reference level without testing the choice; these arms test it.
#
#   source("methodology/v2/index_arms.R")

# Depends on rung_matrix() from plasmode_generate.R for the oracle arm, and on
# INDEX_INPUTS from harmonize_donors.R; the runner and the tests source both.
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(FactoMineR)
  library(cli)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# Coding arms build the same MCA index on the six inputs plus the advice items
# and differ only in how the skip is coded. `impute_frozen` completes target
# items with the development item models; `impute_adaptive` refits them on the
# target's own eligible respondents, which is what the v2 `impute` arm did.
# `passive_mca` gives the skip its own category and then excludes that category
# from the analysis (subset MCA, Greenacre and Pardo 2006; FactoMineR `excl`).
CODING_ARMS <- c("zero", "category", "complete_case", "drop_items",
                 "impute_frozen", "impute_adaptive", "passive_mca")

# A transparent count of externally defined deprivations. No advice item, no
# fitted weight, nothing learned from the generator.
PRACTICAL_ARMS <- c(CODING_ARMS, "simple_score")

# Reference arms bound the comparison and are not deployable recommendations.
# `raw_items` is supervised: every input and item as its own covariate. `oracle`
# knows the generator's rung scores and loadings but not the residual part of V,
# so it is an X-information benchmark, not a ceiling for every metric.
REFERENCE_ARMS <- c("raw_items", "oracle")
ALL_ARMS <- c(PRACTICAL_ARMS, REFERENCE_ARMS)

# Arms that never read an advice item. Under the advice mechanisms none of them
# can move across the access settings; that is asserted, not assumed.
ADVICE_FREE_ARMS <- c("drop_items", "simple_score", "oracle")

# Both imputation arms share one development half: the development item models,
# the m completed development frames and the m MCA fits.
prep_family <- function(arm) {
  if (arm %in% c("impute_frozen", "impute_adaptive")) "impute" else arm
}

# ---- Skip coding --------------------------------------------------------

advice_names <- function(d) grep("^advice[0-9]+$", names(d), value = TRUE)

item_label <- function(x, nm) factor(x, levels = 0:1, labels = paste0(nm, "_", 0:1))

# Rows whose items were actually asked. The gate decides this, not contact:
# under the diagnosis gate a contacted but unrecognized person is a skip too.
observed_rows <- function(d, items = advice_names(d)) {
  if (!length(items)) return(rep(TRUE, nrow(d)))
  stats::complete.cases(d[, items, drop = FALSE])
}

code_items <- function(d, how) {
  for (v in advice_names(d)) {
    d[[v]] <- switch(how,
      zero = item_label(tidyr::replace_na(d[[v]], 0L), v),
      category = factor(if_else(is.na(d[[v]]), 2L, d[[v]]), levels = 0:2,
                        labels = paste0(v, "_", c("0", "1", "notasked"))),
      observed = item_label(d[[v]], v))
  }
  d
}

notasked_levels <- function(d) paste0(advice_names(d), "_notasked")

# The imputation models are fit on eligible respondents only, so a category that
# is rare in the population can be absent from that subset while still appearing
# among the people being imputed. `predict.glm` then aborts with "factor has new
# levels"; the 2026-09-17 pilot lost 231 of 4,000 replicates that way, and not at
# random. For the prediction only, an unseen level is sent to the modal seen
# level, and the count is returned so the substitution is reported.
align_to_fit <- function(newdata, fitted_on, vars) {
  n_fix <- 0L
  for (v in vars) {
    if (!is.factor(newdata[[v]])) next
    seen <- levels(droplevels(fitted_on[[v]]))
    bad <- !as.character(newdata[[v]]) %in% seen
    if (!any(bad)) next
    modal <- names(which.max(table(droplevels(fitted_on[[v]]))))
    newdata[[v]] <- factor(replace(as.character(newdata[[v]]), bad, modal),
                           levels = seen)
    n_fix <- n_fix + sum(bad)
  }
  list(newdata = newdata, n_fix = n_fix)
}

# Terms that can enter a formula. A constant column carries no information and
# `glm` aborts on a one-level factor, which cost three replicates on
# 18 September, all on the least variable draws. The constant term is dropped
# from the fit; the index still sees the original column.
varying <- function(d, vars) {
  keep <- vapply(vars, function(v) {
    x <- d[[v]]
    if (is.factor(x)) nlevels(droplevels(x)) > 1L else length(unique(x[!is.na(x)])) > 1L
  }, logical(1))
  if (!any(keep)) "1" else vars[keep]
}

# One logistic model per item on the eligible rows. It predicts the response a
# person outside the gate would give if asked, which is a hypothetical
# completion, not recovered advice: under the diagnosis gate every eligible
# respondent is recognized, so the model is extrapolating to people who were
# never diagnosed. Items are modelled separately, so their dependence given the
# covariates is not preserved.
item_models <- function(d) {
  items <- advice_names(d)
  obs <- d[observed_rows(d, items), , drop = FALSE]
  rhs <- varying(obs, c(INDEX_INPUTS, "age", "sex"))
  fits <- lapply(items, \(v) stats::glm(stats::reformulate(rhs, v), data = obs,
                                         family = stats::binomial))
  names(fits) <- items
  list(fits = fits, obs = obs, rhs = rhs, items = items)
}

# m completed frames. Responses are drawn from fixed fitted probabilities; model
# parameters are not drawn, so this is improper imputation and no Rubin
# variance is claimed from it. Observed responses are never replaced.
complete_items <- function(d, mod, m) {
  n_fix <- 0L
  probs <- list()
  miss <- list()
  for (v in mod$items) {
    miss[[v]] <- which(is.na(d[[v]]))
    al <- align_to_fit(d[miss[[v]], , drop = FALSE], mod$obs, intersect(mod$rhs, INDEX_INPUTS))
    n_fix <- n_fix + al$n_fix
    probs[[v]] <- if (length(miss[[v]])) {
      stats::predict(mod$fits[[v]], newdata = al$newdata, type = "response")
    } else numeric(0)
  }
  frames <- lapply(seq_len(m), function(i) {
    di <- d
    for (v in mod$items) {
      di[[v]][miss[[v]]] <- stats::rbinom(length(miss[[v]]), 1L, probs[[v]])
      di[[v]] <- item_label(di[[v]], v)
    }
    di
  })
  list(frames = frames, n_fix = n_fix)
}

# ---- Tiers --------------------------------------------------------------

# Quintile-cut the development score and orient so the top tier carries the most
# outcome. Do not remove the flip: dropping an input silently flipped the dim-1
# sign in v1 and produced a Tier-5 OR of 0.06. A score with heavy ties can
# support fewer than five tiers; the top tier is then larger than a fifth, which
# `tier_flag_rate` reports.
dev_tiers <- function(s_dev, dev) {
  cuts <- function(s) {
    brk <- unique(stats::quantile(s, probs = seq(0, 1, 0.2)))
    tier <- if (length(brk) > 1L) cut(s, brk, include.lowest = TRUE, labels = FALSE) else rep(1L, length(s))
    list(brk = brk, tier = tier)
  }
  k <- cuts(s_dev)
  flip <- mean(dev$y[k$tier == max(k$tier)]) < mean(dev$y[k$tier == 1L])
  if (flip) {
    s_dev <- -s_dev
    k <- cuts(s_dev)
  }
  list(score_dev = s_dev, breaks = k$brk, tier_dev = k$tier, flip = flip,
       n_tiers = max(k$tier))
}

# Target scores outside the development range are clamped into it, not left to
# `cut` as NA. A target person more deprived than anyone in development belongs
# in the top tier, and an NA tier propagates silently.
project_tiers <- function(dt, s_tgt) {
  if (dt$flip) s_tgt <- -s_tgt
  brk <- dt$breaks
  tier <- if (length(brk) > 1L) {
    cut(pmin(pmax(s_tgt, min(brk)), max(brk)), brk, include.lowest = TRUE, labels = FALSE)
  } else rep(1L, length(s_tgt))
  list(score = s_tgt, tier = tier)
}

# ---- MCA ----------------------------------------------------------------

# Subset (specific) MCA projection, verified in COMPARATOR_FEASIBILITY.md.
# FactoMineR's predict.MCA ignores `excl`: on the package's own poison data it
# misses the fitted coordinates by .244 and returns undefined cosines. The
# transition formula below is what CA() applies to supplementary rows when
# categories are excluded, and it reproduces the fitted coordinates of the
# training rows to 5e-15.
passive_projector <- function(fit, dev_f, excl) {
  Z <- tab.disjonctif(dev_f)
  fj <- colSums(Z) / (nrow(Z) * ncol(dev_f))
  keep <- as.numeric(!colnames(Z) %in% excl)
  v <- fit$svd$V[, 1]
  list(cols = colnames(Z), q = ncol(dev_f), w = v * keep, offset = sum(fj * v * keep))
}

passive_coord <- function(pj, frame_f) {
  Z <- tab.disjonctif(frame_f)[, pj$cols, drop = FALSE]
  as.numeric((Z / pj$q) %*% pj$w) - pj$offset
}

mca_dev <- function(dev, vars, excl = NULL) {
  dev_f <- droplevels(dev[, vars, drop = FALSE])
  ex <- intersect(excl, unlist(lapply(dev_f, levels), use.names = FALSE))
  fit <- MCA(dev_f, ncp = 1, graph = FALSE, excl = if (length(ex)) ex else NULL)
  list(fit = fit, vars = vars, levels = lapply(dev_f, levels), excl = ex,
       coord_dev = unname(fit$ind$coord[, 1]),
       proj = if (length(ex)) passive_projector(fit, dev_f, ex))
}

# A target row carrying a level development never saw cannot be placed on the
# development axis. It is left unscored and counted against coverage.
mca_coords <- function(mp, tgt) {
  keep <- Reduce(`&`, lapply(mp$vars, \(v) as.character(tgt[[v]]) %in% mp$levels[[v]]))
  coord <- rep(NA_real_, nrow(tgt))
  if (any(keep)) {
    f <- as.data.frame(lapply(mp$vars, \(v) factor(as.character(tgt[[v]][keep]),
                                                   levels = mp$levels[[v]])))
    names(f) <- mp$vars
    coord[keep] <- if (length(mp$excl)) passive_coord(mp$proj, f) else
      unname(stats::predict(mp$fit, newdata = f)$coord[, 1])
  }
  coord
}

# ---- Simple count score -------------------------------------------------

# Equal-weight count of six binary deprivations, each defined outside this
# study. Water, sanitation and fuel follow the global MPI cut-offs (OPHI 2018
# revision): unimproved or surface water; unimproved, open or shared sanitation;
# polluting cooking fuel, where kerosene counts as polluting because the
# harmonized rung merges it with coal. Education is the individual version of
# the MPI schooling cut-off: no schooling or incomplete primary. Work marks manual
# and agricultural occupations, the occupation-class convention of
# socioeconomic position measures. BMI flags both WHO tails, underweight (< 18.5)
# and obese (>= 30); overweight is not flagged. None of these uses the
# generator's rung spacing or loadings.
SIMPLE_DEPRIVED <- c(
  "water_unimproved", "water_surface",
  "sanitation_improved_shared", "sanitation_unimproved", "sanitation_open",
  "fuel_kero_coal", "fuel_charcoal", "fuel_biomass",
  "education_none", "education_incomplete_primary",
  "work_manual", "work_agriculture",
  "bmi_underweight", "bmi_obese")

simple_count <- function(d) {
  rowSums(vapply(INDEX_INPUTS, \(v) as.character(d[[v]]) %in% SIMPLE_DEPRIVED,
                 logical(nrow(d))))
}

# ---- Development half, one per arm family --------------------------------

# Everything that depends on development alone is built once per replicate and
# handed to every access setting: coding, fits, tier cuts, the tier model and
# the Mundlak fit. Morris, White and Crowther (Stat Med 2019;38:2074-2102): the
# contrast between methods is free of that Monte Carlo error only if every
# method is applied to the same simulated data.
#
# Each family returns `dev` (the rows the index was built on), `dt` (tiers) and
# `project(tgt, arm)`, which returns a score on the pre-flip development scale
# with NA for every target row the arm cannot score.
prepare_family <- function(dev_raw, family, spec) {
  out <- switch(family,
    zero = , category = , complete_case = , drop_items = , passive_mca = {
      how <- switch(family, zero = "zero", category = "category",
                    complete_case = "observed", drop_items = "observed",
                    passive_mca = "category")
      dev <- if (family == "complete_case") dev_raw[observed_rows(dev_raw), , drop = FALSE] else dev_raw
      dev <- code_items(dev, how)
      items <- if (family == "drop_items") character(0) else advice_names(dev_raw)
      mp <- mca_dev(dev, c(INDEX_INPUTS, items),
                    excl = if (family == "passive_mca") notasked_levels(dev_raw))
      list(dev = dev, s_dev = mp$coord_dev, n_fix = 0L, project = function(tgt, arm) {
        t <- if (family == "complete_case") tgt[observed_rows(tgt), , drop = FALSE] else tgt
        s <- rep(NA_real_, nrow(tgt))
        s[t$row_id] <- mca_coords(mp, code_items(t, how))
        list(score = s, n_fix = 0L)
      })
    },
    impute = {
      mod <- item_models(dev_raw)
      comp <- complete_items(dev_raw, mod, spec$m_impute)
      mps <- lapply(comp$frames, mca_dev, vars = c(INDEX_INPUTS, advice_names(dev_raw)))
      # Each imputation has its own axis. Orient and standardize each on its own
      # development coordinates, then average: the deployable decision is one
      # averaged score per person, not an average of per-imputation metrics.
      std <- lapply(mps, function(mp) {
        sgn <- if (dev_tiers(mp$coord_dev, dev_raw)$flip) -1 else 1
        list(sgn = sgn, mu = mean(mp$coord_dev), sd = stats::sd(mp$coord_dev))
      })
      z <- function(coords) Reduce(`+`, Map(\(x, st) st$sgn * (x - st$mu) / st$sd, coords, std)) / length(coords)
      list(dev = dev_raw, s_dev = z(lapply(mps, `[[`, "coord_dev")), n_fix = comp$n_fix,
           project = function(tgt, arm) {
             tm <- if (arm == "impute_adaptive") item_models(tgt) else mod
             tc <- complete_items(tgt, tm, spec$m_impute)
             list(score = z(Map(mca_coords, mps, tc$frames)), n_fix = tc$n_fix)
           })
    },
    simple_score = list(dev = dev_raw, s_dev = simple_count(dev_raw), n_fix = 0L,
                        project = function(tgt, arm) list(score = simple_count(tgt), n_fix = 0L)),
    raw_items = raw_items_dev(code_items(dev_raw, "zero"),
                              c(INDEX_INPUTS, advice_names(dev_raw))),
    oracle = {
      raw <- function(d) as.numeric(rung_matrix(d, spec$bmi_risk) %*% spec$loadings)
      r_dev <- raw(dev_raw)
      # Development centre and scale are stored and reused. Standardizing the
      # target on itself preserves its ranks but moves people across the fixed
      # development tier cuts, which is a different policy.
      mu <- mean(r_dev)
      sdv <- stats::sd(r_dev)
      list(dev = dev_raw, s_dev = (r_dev - mu) / sdv, n_fix = 0L,
           project = function(tgt, arm) list(score = (raw(tgt) - mu) / sdv, n_fix = 0L))
    },
    cli_abort("Unknown arm family {.val {family}}."))
  out$family <- family
  out$dt <- dev_tiers(out$s_dev, out$dev)
  out$dm <- dev_metrics(out$dev, out$dt)
  out
}

# Every input and advice item as its own covariate, no composite. The score is a
# fitted linear predictor that has seen the development outcome, so the
# development score is taken out of fold. Folds are whole clusters because every
# target cluster is unseen at deployment. Age and sex are excluded, as they are
# from every MCA score; the tier model adds them for every arm alike.
raw_items_dev <- function(dev, vars, n_folds = 5L) {
  for (v in vars) dev[[v]] <- droplevels(dev[[v]])
  rhs <- varying(dev, vars)
  cl <- unique(dev$cluster)
  folds <- sample(rep_len(seq_len(n_folds), length(cl)))[match(dev$cluster, cl)]
  s_dev <- numeric(nrow(dev))
  for (k in seq_len(n_folds)) {
    train <- dev[folds != k, , drop = FALSE]
    rk <- varying(train, rhs[rhs != "1"])
    fit <- stats::glm(stats::reformulate(rk, "y"), data = train, family = stats::binomial)
    ok <- align_to_fit(dev[folds == k, , drop = FALSE], train, intersect(rk, vars))
    s_dev[folds == k] <- stats::predict(fit, newdata = ok$newdata, type = "link")
  }
  full <- stats::glm(stats::reformulate(rhs, "y"), data = dev, family = stats::binomial)
  lv <- lapply(dev[vars], levels)
  list(dev = dev, s_dev = s_dev, n_fix = 0L, project = function(tgt, arm) {
    t <- code_items(tgt, "zero")
    keep <- Reduce(`&`, lapply(vars, \(v) as.character(t[[v]]) %in% lv[[v]]))
    s <- rep(NA_real_, nrow(t))
    if (any(keep)) {
      tk <- t[keep, , drop = FALSE]
      for (v in vars) tk[[v]] <- factor(as.character(tk[[v]]), levels = lv[[v]])
      s[keep] <- stats::predict(full, newdata = tk, type = "link")
    }
    list(score = s, n_fix = 0L)
  })
}

# ---- Development metrics -------------------------------------------------

# The tier model, an age-and-sex baseline and the Mundlak fit, all on the rows
# the index was built on. Constant terms are dropped: Nepal at 18-49 is women
# only, so `sex` does not vary there. The Mundlak frame is built from scratch:
# `dev` carries the generator's latents, and a formula term sharing a name with
# one of them resolves to the generator's column (the 2026-09-17 pilot did this).
dev_metrics <- function(dev, dt) {
  dev$tier <- factor(dt$tier_dev, levels = seq_len(dt$n_tiers))
  adj <- varying(dev, c("age", "sex"))
  rhs <- varying(dev, c("tier", "age", "sex"))
  fit <- stats::glm(stats::reformulate(rhs, "y"), data = dev, family = stats::binomial)
  base <- stats::glm(stats::reformulate(adj, "y"), data = dev, family = stats::binomial)
  top <- paste0("tier", dt$n_tiers)
  or_top <- if (top %in% names(stats::coef(fit))) exp(unname(stats::coef(fit)[top])) else NA_real_

  mf <- data.frame(y = dev$y, age = dev$age, sex = dev$sex,
                   score = as.numeric(scale(dt$score_dev)))
  mf$zbar <- ave(mf$score, dev$cluster, FUN = mean)
  mf$w <- mf$score - mf$zbar
  mund <- stats::glm(stats::reformulate(c(varying(mf, c("w", "zbar")), adj), "y"),
                     data = mf, family = stats::binomial)
  cf <- stats::coef(mund)
  list(fit = fit, base = base,
       separation = any(fit$fitted.values < 1e-8 | fit$fitted.values > 1 - 1e-8),
       row = tibble(or_top = or_top, n_tiers = dt$n_tiers, n_dev_used = nrow(dev),
                    beta_within = if ("w" %in% names(cf)) unname(cf["w"]) else NA_real_,
                    beta_between = if ("zbar" %in% names(cf)) unname(cf["zbar"]) else NA_real_))
}

# ---- Target metrics ------------------------------------------------------

AGE_BREAKS <- c(17, 29, 44, 59, 200)
AGE_LABELS <- c("18_29", "30_44", "45_59", "60p")
MIN_BAND <- 100L
BUDGETS <- c(b05 = 0.05, b10 = 0.10, b20 = 0.20)
CLUSTER_BUDGET <- 0.10

fast_auc <- function(y, p) {
  ok <- !is.na(p)
  y <- y[ok]
  p <- p[ok]
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  (sum(rank(p)[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

ratio <- function(num, den) if (is.na(den) || den == 0) NA_real_ else num / den

# Truth about the target that does not depend on any arm. The high-need
# threshold is the top fifth of V over the WHOLE target and is computed once;
# no arm's scoreable subset ever redefines it. Cluster truth uses every resident.
target_truth <- function(tgt) {
  cl_ids <- sort(unique(tgt$cluster))
  ci <- match(tgt$cluster, cl_ids)
  cl_v <- as.numeric(tapply(tgt$V, ci, mean))
  n_top <- floor(CLUSTER_BUDGET * length(cl_ids))
  list(N = nrow(tgt),
       hn = tgt$V >= stats::quantile(tgt$V, 0.8),
       unrec = tgt$y == 1L & tgt$recognized == 0L,
       band = cut(tgt$age, AGE_BREAKS, labels = AGE_LABELS),
       ci = ci, n_cl = length(cl_ids), cl_v = cl_v,
       cl_size = tabulate(ci, length(cl_ids)),
       cl_tiebreak = as.numeric(tapply(tgt$cluster_tiebreak, ci, `[`, 1L)),
       cl_median = stats::median(cl_v),
       cl_top = order(-cl_v)[seq_len(n_top)])
}

# Fixed-budget selection among scoreable units. Ties are broken by a uniform
# drawn in the generator, blind to V, Y, recognition, age and arm. If fewer units
# can be scored than the budget allows, the remaining slots stay unfilled: a
# fallback score would be a different policy.
select_top <- function(score, tiebreak, slots) {
  ok <- which(!is.na(score))
  if (!length(ok) || slots < 1L) return(integer(0))
  ok[order(-score[ok], tiebreak[ok])][seq_len(min(slots, length(ok)))]
}

# Cluster-level targeting on the whole target. A cluster's priority is the mean
# score of its scored residents, which is a deployable estimate, not a mean over
# every resident; the coverage columns say how thin it is. A cluster with no
# scored resident cannot be chosen: it ranks tied below every scored cluster in
# the Spearman correlation and is never selected. With no scored cluster at all
# the targeting columns are NA.
cluster_metrics <- function(s_all, tr) {
  n_scored_cl <- tabulate(tr$ci[!is.na(s_all)], tr$n_cl)
  cs <- rep(NA_real_, tr$n_cl)
  has <- n_scored_cl > 0L
  if (any(has)) {
    m <- tapply(s_all[!is.na(s_all)], tr$ci[!is.na(s_all)], mean)
    cs[as.integer(names(m))] <- as.numeric(m)
  }
  out <- list(cluster_score_coverage = mean(has),
              cluster_resident_coverage = mean(n_scored_cl / tr$cl_size))
  if (!any(has)) {
    return(c(out, list(cluster_rho_all = NA_real_, misallocation_all = NA_real_,
                       cl10_slots = floor(CLUSTER_BUDGET * tr$n_cl), cl10_used = 0L,
                       cl10_misalloc = NA_real_, cl10_precision = NA_real_,
                       cl10_hn_reached = 0, cl10_residents = 0L)))
  }
  s_rank <- ifelse(has, cs, min(cs, na.rm = TRUE) - 1)
  rho <- if (stats::sd(s_rank) > 0) stats::cor(s_rank, tr$cl_v, method = "spearman") else NA_real_
  top_q <- s_rank >= stats::quantile(s_rank, 0.9)
  slots <- floor(CLUSTER_BUDGET * tr$n_cl)
  sel <- select_top(cs, tr$cl_tiebreak, slots)
  in_sel <- tr$ci %in% sel
  c(out, list(
    cluster_rho_all = rho,
    misallocation_all = mean(tr$cl_v[top_q] < tr$cl_median),
    cl10_slots = slots, cl10_used = length(sel),
    cl10_misalloc = if (length(sel)) mean(tr$cl_v[sel] < tr$cl_median) else NA_real_,
    cl10_precision = if (length(sel)) mean(sel %in% tr$cl_top) else NA_real_,
    cl10_hn_reached = ratio(sum(tr$hn & in_sel), sum(tr$hn)),
    cl10_residents = sum(tr$cl_size[sel])))
}

# What a validation restricted to scoreable rows would report: the top fifth of
# V is re-drawn inside the scored rows, and clusters are ranked on scored rows
# only. Kept as diagnostics, never as the comparison between arms.
scored_diagnostics <- function(tgt, s, flagged) {
  vul <- tgt$V >= stats::quantile(tgt$V, 0.8)
  cl <- tibble(cluster = tgt$cluster, s = s, v = tgt$V) |>
    group_by(cluster) |>
    summarise(s = mean(s), v = mean(v), .groups = "drop")
  top <- cl$s >= stats::quantile(cl$s, 0.9)
  list(capture_scored = mean(flagged[vul]),
       cluster_rho_scored = if (nrow(cl) > 2L && stats::sd(cl$s) > 0) {
         stats::cor(cl$s, cl$v, method = "spearman")
       } else NA_real_,
       misallocation_scored = mean(cl$v[top] < stats::median(cl$v)))
}

# Every metric for one arm at one access setting. See METRIC_SPECIFICATION.md
# for the definition, denominator and failure behaviour of each column.
target_metrics <- function(prep, sc, tgt, tr) {
  dm <- prep$dm
  n_t <- prep$dt$n_tiers
  scored <- which(!is.na(sc$score))
  pt <- project_tiers(prep$dt, sc$score[scored])
  s_all <- rep(NA_real_, tr$N)
  s_all[scored] <- pt$score
  tier_all <- rep(NA_integer_, tr$N)
  tier_all[scored] <- pt$tier
  flagged <- !is.na(tier_all) & tier_all == n_t
  n_hn <- sum(tr$hn)
  n_unrec <- sum(tr$unrec)

  out <- list(N = tr$N, n_scored = length(scored), coverage = length(scored) / tr$N,
              n_hn = n_hn, n_unrec = n_unrec, n_clusters = tr$n_cl,
              level_fixups = prep$n_fix + sc$n_fix,
              tier_flag_rate = mean(flagged),
              capture_all = ratio(sum(flagged & tr$hn), n_hn),
              unrec_capture_tier = ratio(sum(flagged & tr$unrec), n_unrec))

  for (b in names(BUDGETS)) {
    slots <- floor(BUDGETS[[b]] * tr$N)
    sel <- select_top(s_all, tgt$tiebreak, slots)
    hn_sel <- sum(tr$hn[sel])
    ur_sel <- sum(tr$unrec[sel])
    out[paste0(b, c("_slots", "_used", "_hn", "_unrec", "_hn_capture", "_hn_ppv",
                    "_unrec_capture", "_hn_shortfall_per1000"))] <- list(
      slots, length(sel), hn_sel, ur_sel, ratio(hn_sel, n_hn), ratio(hn_sel, length(sel)),
      ratio(ur_sel, n_unrec), 1000 * (min(slots, n_hn) - hn_sel) / tr$N)
    if (b == "b10") sel10 <- sel
  }

  in10 <- seq_len(tr$N) %in% sel10
  for (b in AGE_LABELS) {
    i <- !is.na(tr$band) & tr$band == b
    den <- sum(tr$hn & i)
    ok <- den >= MIN_BAND
    out[[paste0("n_hn_", b)]] <- den
    out[[paste0("cap_tier_", b)]] <- if (ok) sum(flagged & tr$hn & i) / den else NA_real_
    out[[paste0("cap_b10_", b)]] <- if (ok) sum(in10 & tr$hn & i) / den else NA_real_
  }

  out <- c(out, cluster_metrics(s_all, tr))

  # The disease model and the scored-only diagnostics need at least a handful of
  # scored rows with both outcomes present.
  ts <- tgt[scored, , drop = FALSE]
  ts$tier <- factor(pt$tier, levels = seq_len(n_t))
  enough <- length(scored) >= MIN_BAND && length(unique(ts$y)) == 2L
  if (enough) {
    lp <- stats::predict(dm$fit, newdata = ts, type = "link")
    cal <- stats::coef(stats::glm(ts$y ~ lp, family = stats::binomial))
    fit_t <- stats::glm(stats::reformulate(varying(ts, c("tier", "age", "sex")), "y"),
                        data = ts, family = stats::binomial)
    top <- paste0("tier", n_t)
    out <- c(out, list(
      or_top_target = if (top %in% names(stats::coef(fit_t))) exp(unname(stats::coef(fit_t)[top])) else NA_real_,
      separation_target = any(fit_t$fitted.values < 1e-8 | fit_t$fitted.values > 1 - 1e-8),
      auc_model_scored = fast_auc(ts$y, lp),
      auc_index_scored = fast_auc(ts$y, pt$score),
      auc_agesex_scored = fast_auc(ts$y, stats::predict(dm$base, newdata = ts, type = "link")),
      cal_slope = unname(cal[2]), cal_intercept = unname(cal[1])),
      scored_diagnostics(ts, pt$score, flagged[scored]))
    band_s <- tr$band[scored]
    vul_s <- ts$V >= stats::quantile(ts$V, 0.8)
    for (b in AGE_LABELS) {
      i <- which(!is.na(band_s) & band_s == b)
      okb <- length(i) >= MIN_BAND && length(unique(ts$y[i])) == 2L
      out[[paste0("auc_scored_", b)]] <- if (okb) fast_auc(ts$y[i], lp[i]) else NA_real_
      j <- i[vul_s[i]]
      out[[paste0("capture_scored_", b)]] <- if (length(j) >= MIN_BAND) mean(flagged[scored][j]) else NA_real_
    }
  } else {
    out <- c(out, list(or_top_target = NA_real_, separation_target = NA,
                       auc_model_scored = NA_real_, auc_index_scored = NA_real_,
                       auc_agesex_scored = NA_real_, cal_slope = NA_real_,
                       cal_intercept = NA_real_, capture_scored = NA_real_,
                       cluster_rho_scored = NA_real_, misallocation_scored = NA_real_))
    for (b in AGE_LABELS) out[paste0(c("auc_scored_", "capture_scored_"), b)] <- NA_real_
  }
  bind_cols(dm$row, tibble(separation_dev = dm$separation), as_tibble(out))
}

# ---- One arm ------------------------------------------------------------

# One arm at one access setting. The runner calls prepare_family once per
# replicate and this once per arm and setting.
score_prepared <- function(prep, arm, tgt, spec) {
  tgt$row_id <- seq_len(nrow(tgt))
  target_metrics(prep, prep$project(tgt, arm), tgt, target_truth(tgt))
}

# Convenience for tests and one-off scoring.
score_arm <- function(dev_raw, tgt_raw, arm, spec) {
  score_prepared(prepare_family(dev_raw, prep_family(arm), spec), arm, tgt_raw, spec)
}
