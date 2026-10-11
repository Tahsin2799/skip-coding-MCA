# Plasmode generator for the v2 simulation, amended for v03.
#
# Covariates are resampled from real DHS records and never modelled. True
# vulnerability, the outcome, cluster access supply, contact, recognition and the
# advice items are imposed from known conditionals. Access is simulated and never
# resampled: resampling a treatment-like variable alongside the covariates breaks
# positivity.
#
#   source("methodology/v2/plasmode_generate.R")
#
# Needs the canonical ladders from harmonize_donors.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(cli)
})

# ---- True vulnerability -------------------------------------------------

# Both tails of BMI mark risk, so this input is deliberately non-monotone. An
# index that orders the four classes monotonically must invert one tail. The
# `bmi_monotone` control swaps in a monotone truth through `spec$bmi_risk`.
BMI_RISK <- c(underweight = 0.70, normal = 0.00, overweight = 0.50, obese = 1.00)

# Rung scores on [0, 1], worst = 1, looked up on the canonical ladder rather than
# on whatever levels a pool happens to hold. Until v03 the spacing came from
# nlevels() after droplevels(), so a rung empty in a pool changed the score of
# every other rung: charcoal is empty in all four donors, which put kero_coal at
# .50 instead of 1/3. A common construct across populations needs common spacing.
RUNG_LADDERS <- list(water = WATER_LEVELS, sanitation = SANITATION_LEVELS,
                     fuel = FUEL_LEVELS, education = rev(EDUCATION_LEVELS),
                     work = WORK_LEVELS)

rung_matrix <- function(d, bmi_risk = BMI_RISK) {
  rung <- function(v) {
    lad <- RUNG_LADDERS[[v]]
    pos <- match(sub(paste0("^", v, "_"), "", as.character(d[[v]])), lad)
    (pos - 1) / (length(lad) - 1)
  }
  cbind(water = rung("water"), sanitation = rung("sanitation"), fuel = rung("fuel"),
        education = rung("education"), work = rung("work"),
        bmi = unname(bmi_risk[sub("^bmi_", "", as.character(d$bmi))]))
}

# The signal is standardized with constants from the whole pool, not from each
# draw, so V means the same thing in the development and target halves of a
# replicate. Computed once per population and specification.
signal_scale <- function(pool, spec) {
  s <- rung_matrix(pool, spec$bmi_risk) %*% spec$loadings
  c(center = mean(s), scale = stats::sd(s))
}

# V is a weighted sum of the rung scores plus individual noise. The index sees
# the categories only, never the rung scores and never the weights, so it has to
# recover the ordering instead of being handed it. `r2` is the share of V the six
# inputs carry; the remainder is vulnerability they cannot see.
true_vulnerability <- function(d, spec, scl) {
  m <- rung_matrix(d, spec$bmi_risk)
  bad <- which(!stats::complete.cases(m))
  if (length(bad)) {
    cli_abort(c(
      "Rung matrix has {length(bad)} incomplete row{?s} out of {nrow(m)}.",
      i = "First offending row is {bad[1]}; check that every input carries a canonical level."
    ))
  }
  signal <- (as.numeric(m %*% spec$loadings) - scl[["center"]]) / scl[["scale"]]
  sqrt(spec$r2) * signal + sqrt(1 - spec$r2) * stats::rnorm(nrow(m))
}

# ---- Population draw ----------------------------------------------------

# Clusters are drawn with replacement so the real within-cluster correlation of
# water, sanitation and fuel survives. A cluster drawn twice becomes two distinct
# clusters, otherwise the duplicate rows would look like one oversized cluster.
draw_population <- function(pool, n_rows, tag) {
  ids <- unique(pool$cluster)
  sizes <- as.integer(table(pool$cluster)[ids])
  take <- sample.int(length(ids), size = ceiling(2 * n_rows / mean(sizes)), replace = TRUE)
  run <- cumsum(sizes[take])
  take <- take[seq_len(max(which(run <= n_rows), 1L))]
  idx <- split(seq_len(nrow(pool)), pool$cluster)
  rows <- unlist(idx[ids[take]], use.names = FALSE)
  out <- pool[rows, , drop = FALSE]
  out$cluster <- rep(paste0(tag, "_", seq_along(take)), times = sizes[take])
  out <- out[seq_len(min(nrow(out), n_rows)), , drop = FALSE]
  out$row_id <- seq_len(nrow(out))
  out
}

# ---- Mechanisms ---------------------------------------------------------

# Cluster access supply. `rho` is its correlation with cluster mean deprivation:
# negative is the inverse care law, anchored at about -.25 on Nepal 2022 wbp15.
cluster_access <- function(cluster, v, rho) {
  bar <- as.numeric(scale(tapply(v, cluster, mean)))
  raw <- stats::rnorm(length(bar))
  a <- rho * bar + sqrt(1 - rho^2) * raw
  unname(a[match(cluster, sort(unique(cluster)))])
}

impose_outcome <- function(d, spec) {
  u <- stats::rnorm(length(unique(d$cluster)), sd = sqrt(spec$sigma2))
  d$u <- u[match(d$cluster, sort(unique(d$cluster)))]
  # Named `eps_contact`, not `z`: a column called `z` in this frame is silently
  # picked up by any downstream glm whose formula names `z`, which is how the
  # pilot's Mundlak fit ended up regressing the outcome on this noise term.
  d$eps_contact <- stats::rnorm(nrow(d))
  # Age enters per decade from the published odds ratio, so the coefficient is a
  # property of age and not of the population's age spread.
  d$risk <- stats::plogis(spec$outcome_intercept[[spec$population]] +
                            spec$v_loading * d$V + d$u +
                            spec$age_log_or_per_decade * (d$age - mean(d$age)) / 10 +
                            spec$sex_slope * (d$sex == "male"))
  d$y <- stats::rbinom(nrow(d), 1L, d$risk)
  d
}

# Everything that decides who is contacted, recognized and advised is a fixed
# uniform drawn once per person, then compared with a setting-specific
# probability. So the access settings are common random numbers: anyone
# contacted at c20 is contacted at c50 and c80, and a person's status changes
# across settings only because access changed. v1 did this; v2 redrew each
# setting with rbinom, which made individual migration across settings partly
# noise. The tie-break uniforms are for budget selection and see nothing else.
draw_uniforms <- function(d, n_items) {
  d$u_contact <- stats::runif(nrow(d))
  d$u_recog <- stats::runif(nrow(d))
  for (k in seq_len(n_items)) d[[paste0("u_item", k)]] <- stats::runif(nrow(d))
  d$tiebreak <- stats::runif(nrow(d))
  cl <- sort(unique(d$cluster))
  d$cluster_tiebreak <- stats::runif(length(cl))[match(d$cluster, cl)]
  d
}

# Contact, recognition, then advice conditional on recognition.
#
# Factual advice receipt needs contact. Whether the item is OBSERVED depends on
# the questionnaire gate, which is separate: `contact` asks everyone who reached
# care (G = C), `diagnosis` asks only people told they have the condition
# (G = R), which is what BDHS 2017-18 Q317B-Q318 and Kibria et al 2024 describe.
# Everyone outside the gate is a structural skip, and how that skip is coded is
# the thing the study is about. Under the diagnosis gate an observed respondent
# is always recognized, so the recognition term cannot be estimated from observed
# advice.
#
# `proxy` items are the always-observed risk-proxy control: functions of V only,
# asked of everyone, so access cannot move them.
impose_access <- function(d, setting, spec) {
  p_contact <- stats::plogis(spec$contact_intercepts[[setting]] + 0.4 * d$eps_contact +
                               0.3 * d$V + spec$access_kappa * d$A)
  d$contacted <- as.integer(d$u_contact < p_contact)
  # Recognition given disease and given contact rises with age (Foti and Hardy,
  # Am J Hypertens 2025;38:529-531). The direction is sourced; the slope is a
  # design choice, and `q1_age_slope = 0` gives a constant q1.
  q1_i <- stats::plogis(stats::qlogis(spec$q1) +
                          spec$q1_age_slope * as.numeric(scale(d$age)))
  p_recog <- if_else(d$y == 1L, q1_i, spec$q0)
  d$recognized <- d$contacted * as.integer(d$u_recog < p_recog)
  d$eligible <- switch(spec$gate,
                       contact = d$contacted,
                       diagnosis = d$recognized,
                       cli_abort("Unknown gate {.val {spec$gate}}."))
  for (k in seq_along(spec$advice_probabilities)) {
    base <- stats::qlogis(spec$advice_probabilities[k])
    u <- d[[paste0("u_item", k)]]
    d[[paste0("advice", k)]] <- if (spec$item_mode == "proxy") {
      as.integer(u < stats::plogis(base + spec$advice_v_loading * d$V))
    } else {
      got <- d$contacted * as.integer(u < stats::plogis(base + log(spec$advice_or) * d$recognized +
                                                         spec$advice_v_loading * d$V))
      if_else(d$eligible == 1L, got, NA_integer_)
    }
  }
  if (spec$item_mode == "proxy") d$eligible <- 1L
  d$setting <- setting
  d
}

# ---- One replicate ------------------------------------------------------

# `scl` is signal_scale(pool, spec); the runner computes it once per population.
generate_replicate <- function(pool, spec, n_dev, n_target, scl = signal_scale(pool, spec)) {
  dev <- draw_population(pool, n_dev, "dev")
  tgt <- draw_population(pool, n_target, "tgt")
  for (nm in c("dev", "tgt")) {
    d <- get(nm)
    d$V <- true_vulnerability(d, spec, scl)
    d$A <- cluster_access(d$cluster, d$V, spec$access_rho)
    d <- impose_outcome(d, spec)
    assign(nm, draw_uniforms(d, length(spec$advice_probabilities)))
  }
  list(development = impose_access(dev, "c50", spec),
       targets = lapply(spec$settings, \(s) impose_access(tgt, s, spec)) |>
         setNames(spec$settings))
}
