# Configuration for the v2 plasmode simulation, amended for v03.
#
#   source("methodology/v2/pilot_v2_config.R")
#
# Needs BMI_RISK from plasmode_generate.R for the default truth.

suppressPackageStartupMessages({library(tibble); library(tidyr); library(cli)})

SCHEMA_VERSION <- "v03-2026-09-20"

# Outcome prevalence target, and the intercepts that hit it. The intercept is
# population-specific because it has to absorb that population's covariate and
# age distribution; the slopes above it are not. Solved by solve_intercepts.R,
# and tests/test_generator.R asserts on independent draws that the realized rates
# still match. Rerun the solver after any change to the pools, the rung ladders,
# the construction of V, the age term, `v_loading`, `sigma2` or `sex_slope`.
OUTCOME_PREVALENCE <- 0.25

# Solved 20 September on canonical rung spacing and pool-standardized V, 20
# replicates of 20,000 per population.
OUTCOME_INTERCEPTS <- c(Q1 = -1.4018739424817,
                        Q2 = -1.3874000585407,
                        Q3 = -1.3979337882166,
                        Q4 = -1.3897578168946,
                        bd_allage = -1.5156911575790,
                        donor_india = -1.3925435620646,
                        donor_nepal = -1.3550533785203,
                        donor_albania = -1.3936592700210,
                        donor_bangladesh = -1.3962307474406)

# One set for every population: the solved contact intercepts differ by at most
# .04 on the logit scale, which the realized-rate test absorbs. bd_allage needs
# its own outcome intercept because its age spread is 16.2 years against the
# quadrants' 9.1 and age enters per decade.
CONTACT_INTERCEPTS <- c(c20 = -1.5831718277591, c50 = 0.0030644116036, c80 = 1.5892429868911)

# The populations the runner knows how to build. `Q1`-`Q4` are India strata
# crossed on deprivation level and coherence; `bd_allage` is Bangladesh 18-95;
# the donors are each country at 18-49, the window in which all six inputs exist
# for every one of them.
QUADRANTS  <- c("Q1", "Q2", "Q3", "Q4")
DONOR_POPULATIONS <- c(donor_india = "india", donor_nepal = "nepal",
                       donor_albania = "albania", donor_bangladesh = "bangladesh")
ALL_POPULATIONS <- c(QUADRANTS, "bd_allage", names(DONOR_POPULATIONS))

GATES <- c("contact", "diagnosis")

pilot_spec <- function(q1 = 0.50, access_rho = -0.25, m_impute = 5L, r2 = 0.6,
                       q1_age_slope = 0, population = "Q3", gate = "contact") {
  list(
    version            = SCHEMA_VERSION,
    population         = population,
    gate               = gate,
    outcome_intercept  = OUTCOME_INTERCEPTS,
    v_loading          = 0.8,
    # Log odds of the outcome per decade of age. Abba et al, BMC Public Health
    # 2021;21:2218: OR 1.74 per 10 years across 12 DHS countries.
    age_log_or_per_decade = log(1.74),
    sex_slope          = 0.1,                        # design
    sigma2             = (0.05 / 0.95) * pi^2 / 3,   # cluster ICC .05, sourced
    contact_intercepts = CONTACT_INTERCEPTS,
    settings           = c("c20", "c50", "c80"),
    # Recognition among contacted non-cases. A design value inside the false
    # positive range implied by self-report specificity of .945 to .964; those
    # studies measure a different conditional, see PROTOCOL_AMENDMENT section 4.
    q0                 = 0.05,
    # Recognition among contacted cases. A design anchor bracketing pooled LMIC
    # awareness of .392 (Geldsetzer 2019), which is not conditional on contact.
    q1                 = q1,
    # Log-odds of recognition per SD of age. Direction sourced (Foti and Hardy
    # 2025); magnitude is design.
    q1_age_slope       = q1_age_slope,
    advice_probabilities = c(0.10, 0.10, 0.05, 0.15),
    # Advice odds ratio for recognized against unrecognized contacted people. A
    # design choice: BDHS never asks the unrecognized, so no survey observes this
    # contrast. The implied probabilities given recognition, .46 to .74, sit
    # beside Kibria et al's 54 to 77 percent salt advice among the diagnosed.
    advice_or          = 16,
    advice_v_loading   = 0.25,                       # design
    item_mode          = "advice",
    bmi_risk           = BMI_RISK,
    # Approximate calibration against Nepal 2022 wbp15, not estimates of the
    # latent access model's coefficients.
    access_kappa       = 0.736,
    access_rho         = access_rho,
    loadings           = c(0.8, 0.7, 0.6, 0.6, 0.5, 0.5),   # design, v1 values
    r2                 = r2,                         # design; a grid axis
    m_impute           = m_impute
  )
}

# ---- Variants -----------------------------------------------------------

# One-at-a-time variations around the reference cell (Q3, r2 .6, slope 0,
# contact gate). Each asks whether a conclusion survives a mechanism moving.
# Variants that change a term inside p_contact or the construction of V carry
# their own solved intercepts; solve_intercepts.R --stress solves them.
STRESS <- list(
  q1_low     = list(q1 = 0.20),
  q1_high    = list(q1 = 0.80),
  q0_high    = list(q0 = 0.15),
  rho_zero   = list(access_rho = 0),
  rho_strong = list(access_rho = -0.5),
  kappa_zero = list(access_kappa = 0),
  items_two   = list(advice_probabilities = c(0.10, 0.10)),
  items_eight = list(advice_probabilities = rep(c(0.10, 0.10, 0.05, 0.15), 2L)),
  loadings_flat = list(loadings = rep(0.65, 6)),
  ndev_small  = list(),
  impute_m20  = list(m_impute = 20L)
)

# Mechanism controls, built to challenge the conclusion rather than confirm it.
CONTROLS <- list(
  # Recognition no longer depends on disease: cases and non-cases are recognized
  # at the same rate, .1625, which keeps the expected recognized share among the
  # contacted at the reference value when prevalence is .25. The advice items
  # then carry no disease signal beyond V; the direct V and access paths remain.
  recog_flat    = list(q1 = 0.1625, q0 = 0.1625),
  # Four always-observed risk proxies instead of advice: functions of V only,
  # prevalence .5, asked of everyone. v1's stress_baseline_proxy_4, rebuilt for
  # the plasmode. Separates "more items" from "items behind a care gate".
  proxy_block   = list(item_mode = "proxy", advice_probabilities = rep(0.5, 4),
                       advice_v_loading = 1.2),
  # The advice block carries much more information about need beyond the six
  # inputs. The boundary of dropping the items: here it should lose.
  advice_strong = list(advice_v_loading = 1.0),
  # A monotone BMI truth, so the simple score's two-tail rule is tested against
  # a relationship it was not written for.
  bmi_monotone  = list(bmi_risk = c(underweight = 0, normal = 1 / 3,
                                    overweight = 2 / 3, obese = 1))
)

VARIANT_NDEV <- c(ndev_small = 500L)

`%||%` <- function(x, y) if (is.null(x)) y else x

# Variants that change a term inside p_contact or the construction of V need
# their own contact intercepts, solved at Q3 on 20 replicates of 20,000 by
# solve_intercepts.R --stress. That run also reports what each variant does to
# the outcome intercept: every shift is under .005 on the logit scale, inside
# the .02 prevalence tolerance the generator test enforces, so the population
# outcome intercept stands. Re-key it if a future variant moves it past .01.
VARIANT_INTERCEPTS <- list(
  rho_zero      = c(c20 = -1.5940221939971, c50 = 0.0042973493238, c80 = 1.6014012231968),
  rho_strong    = c(c20 = -1.5978230961461, c50 = -0.0233018728466, c80 = 1.5523706500620),
  kappa_zero    = c(c20 = -1.4580612563225, c50 = 0.0002043784744, c80 = 1.4584134522060),
  loadings_flat = c(c20 = -1.6008579919599, c50 = -0.0143568055093, c80 = 1.5725634087540),
  bmi_monotone  = c(c20 = -1.6009325417957, c50 = -0.0146981182331, c80 = 1.5705678407246)
)

variant_overrides <- function(variant) {
  ov <- c(STRESS, CONTROLS)[[variant]] %||% list()
  ci <- VARIANT_INTERCEPTS[[variant]]
  if (!is.null(ci)) ov$contact_intercepts <- ci
  ov
}

# ---- The v03 grid -------------------------------------------------------

# Every factor is in the configuration id, including the ones held at a single
# value, so no two specifications can share a file or a seed stream. The stream
# id is the configuration id without the gate: the two gates of one cell see the
# same resampled people and the same uniforms, so the gate contrast is paired.
#
# Blocks: `core` crosses population, r2, recognition slope and gate; `stress`
# and `controls` sit at the reference cell under the contact gate, with the two
# recognition stresses mirrored under the diagnosis gate because q1 sets the
# size of the eligible group there; `quadrant` adds Q1, Q2 and Q4 at the
# reference cell (Q3 is the core cell itself, not a second copy); `donor` runs
# each country at 18-49 under both gates.
v03_grid <- function(blocks = c("core", "stress", "controls", "quadrant", "donor"),
                     n_dev = 1000L, n_target = 20000L) {
  row <- function(block, population, r2 = 0.6, q1_age_slope = 0, gate = "contact",
                  variant = "base") {
    tibble(block = block, population = population, r2 = r2,
           q1_age_slope = q1_age_slope, gate = gate, variant = variant)
  }
  parts <- list()
  if ("core" %in% blocks) {
    g <- tidyr::expand_grid(population = c("Q3", "bd_allage"), r2 = c(0.4, 0.6, 0.8),
                            q1_age_slope = c(0, 0.5), gate = GATES)
    parts$core <- row("core", g$population, g$r2, g$q1_age_slope, g$gate)
  }
  if ("stress" %in% blocks) {
    parts$stress <- dplyr::bind_rows(
      row("stress", "Q3", variant = names(STRESS)),
      row("stress", "Q3", gate = "diagnosis", variant = c("q1_low", "q1_high")))
  }
  if ("controls" %in% blocks) parts$controls <- row("controls", "Q3", variant = names(CONTROLS))
  if ("quadrant" %in% blocks) parts$quadrant <- row("quadrant", c("Q1", "Q2", "Q4"))
  if ("donor" %in% blocks) {
    g <- tidyr::expand_grid(population = names(DONOR_POPULATIONS), gate = GATES)
    parts$donor <- row("donor", g$population, gate = g$gate)
  }
  g <- dplyr::bind_rows(parts)
  g$n_dev <- unname(ifelse(g$variant %in% names(VARIANT_NDEV), VARIANT_NDEV[g$variant], n_dev))
  g$n_target <- n_target
  stem <- sprintf("%s_nd%d_nt%d_r%02.0f_q%02.0f", g$population, g$n_dev, g$n_target,
                  g$r2 * 100, g$q1_age_slope * 100)
  g$stream_id <- paste0(stem, "_", g$variant)
  g$config_id <- paste0(stem, "_g", substr(g$gate, 1, 3), "_", g$variant)
  g$spec <- lapply(seq_len(nrow(g)), function(i) {
    sp <- pilot_spec(r2 = g$r2[i], q1_age_slope = g$q1_age_slope[i],
                     population = g$population[i], gate = g$gate[i])
    modifyList(sp, variant_overrides(g$variant[i]))
  })
  if (anyDuplicated(g$config_id)) {
    cli::cli_abort(c("Configuration ids are not unique.",
                     i = "Duplicated: {.val {unique(g$config_id[duplicated(g$config_id)])}}"))
  }
  g
}
