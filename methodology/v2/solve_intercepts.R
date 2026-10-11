# Solve the outcome and contact intercepts against the realized generator.
#
# Each intercept is solved against the full model, cluster access term and all,
# holding a set of drawn populations fixed so the root find sees a smooth curve.
# The solver seeds are its own ("solve_..."); tests/test_generator.R checks the
# realized rates on independent draws, so the solving sample never doubles as
# the validation sample.
#
#   Rscript methodology/v2/solve_intercepts.R --reps 20
#   Rscript methodology/v2/solve_intercepts.R --reps 20 --stress

suppressPackageStartupMessages({
  library(dplyr)
  library(digest)
})

root <- c(".", "..", "../..") |>
  Filter(f = \(p) file.exists(file.path(p, "methodology/v2/plasmode_generate.R"))) |>
  (\(x) x[1])()

source(file.path(root, "methodology/v2/harmonize_donors.R"))
source(file.path(root, "methodology/v2/plasmode_generate.R"))
source(file.path(root, "methodology/v2/index_arms.R"))
source(file.path(root, "methodology/v2/quadrant_populations.R"))
source(file.path(root, "methodology/v2/pilot_v2_config.R"))

outcome_at <- function(alpha, parts, spec) {
  mean(vapply(parts, \(p) mean(stats::plogis(
    alpha + spec$v_loading * p$V + p$u +
      spec$age_log_or_per_decade * (p$age - mean(p$age)) / 10 +
      spec$sex_slope * p$male)), numeric(1)))
}

contact_at <- function(alpha, parts, spec) {
  mean(vapply(parts, \(p) mean(stats::plogis(
    alpha + 0.4 * p$eps_contact + 0.3 * p$V + spec$access_kappa * p$A)), numeric(1)))
}

# One population's draws, with every term of both linear predictors except the
# intercept being solved.
draw_parts <- function(pool, spec, n_rep, n_target, tag) {
  scl <- signal_scale(pool, spec)
  lapply(seq_len(n_rep), function(i) {
    set.seed(digest2int(paste0("solve_", tag, "_", i)))
    d <- draw_population(pool, n_target, "tgt")
    d$V <- true_vulnerability(d, spec, scl)
    d$A <- cluster_access(d$cluster, d$V, spec$access_rho)
    u <- stats::rnorm(length(unique(d$cluster)), sd = sqrt(spec$sigma2))
    list(V = d$V, A = d$A, age = d$age, male = as.numeric(d$sex == "male"),
         u = u[match(d$cluster, sort(unique(d$cluster)))],
         eps_contact = stats::rnorm(nrow(d)))
  })
}

solve_for <- function(f, target, parts, spec) {
  stats::uniroot(\(a) f(a, parts, spec) - target, interval = c(-8, 8), tol = 1e-8)$root
}

build_pools <- function() {
  harmonized <- lapply(readRDS(file.path(root, "outputs/methodology/donors_harmonized_geo.rds")),
                       prefix_levels)
  keep <- function(d) droplevels(as.data.frame(d[, c(INDEX_INPUTS, "age", "sex", "cluster")]))
  pools <- lapply(build_quadrants(harmonized,
    quadrant_cells(file.path(root, "methodology/v2/quadrant_cells.csv"))), keep)
  pools[["bd_allage"]] <- keep(build_age_pool(harmonized, "bangladesh", c(18L, Inf)))
  for (nm in names(DONOR_POPULATIONS)) {
    pools[[nm]] <- keep(build_age_pool(harmonized, DONOR_POPULATIONS[[nm]], c(18L, 49L)))
  }
  pools
}

if (sys.nframe() == 0L) {
  argv <- commandArgs(trailingOnly = TRUE)
  i <- match("--reps", argv)
  n_rep <- if (is.na(i)) 10L else as.integer(argv[i + 1L])
  targets <- c(c20 = 0.20, c50 = 0.50, c80 = 0.80)
  pools <- build_pools()

  # Outcome and contact intercepts depend on V, age, sex and the cluster terms,
  # none of which the gate or the r2 grid touches except through V, which is
  # standardized; r2 = .6 is the reference and the tests check .4 and .8.
  spec <- pilot_spec()
  rows <- list()
  out_int <- numeric(0)
  con_int <- list()
  for (nm in names(pools)) {
    parts <- draw_parts(pools[[nm]], spec, n_rep, 20000L, nm)
    out_int[nm] <- solve_for(outcome_at, OUTCOME_PREVALENCE, parts, spec)
    con_int[[nm]] <- vapply(targets, solve_for, numeric(1), f = contact_at, parts = parts, spec = spec)
    rows[[nm]] <- data.frame(population = nm, outcome_intercept = out_int[[nm]],
                             c20 = con_int[[nm]][["c20"]], c50 = con_int[[nm]][["c50"]],
                             c80 = con_int[[nm]][["c80"]])
  }
  cat(sprintf("\nsolved on %d replicates of 20,000 per population\n\n", n_rep))
  print(do.call(rbind, rows), row.names = FALSE, digits = 6)
  cm <- do.call(rbind, con_int)
  cat(sprintf("\ncontact intercepts vary across populations by at most %.4f on the logit scale\n",
              max(apply(cm, 2, \(x) diff(range(x))))))

  cat("\n# paste into pilot_v2_config.R\n")
  cat("OUTCOME_INTERCEPTS <- c(",
      paste(sprintf("%s = %.13f", names(out_int), out_int), collapse = ",\n                        "),
      ")\n", sep = "")
  cat("CONTACT_INTERCEPTS <- c(",
      paste(sprintf("%s = %.13f", colnames(cm), colMeans(cm)), collapse = ", "), ")\n", sep = "")

  # Variants that change the contact model or the construction of V, solved at
  # the reference population. The outcome intercept is reported for each one; a
  # shift above .01 on the logit scale gets its own value in the config.
  if ("--stress" %in% argv) {
    cat("\nvariants solved at Q3\n\n")
    for (nm in names(VARIANT_INTERCEPTS)) {
      ov <- c(STRESS, CONTROLS)[[nm]]
      sp <- modifyList(pilot_spec(population = "Q3"), ov)
      parts <- draw_parts(pools[["Q3"]], sp, n_rep, 20000L, paste0("Q3_", nm))
      a_out <- solve_for(outcome_at, OUTCOME_PREVALENCE, parts, sp)
      a_con <- vapply(targets, solve_for, numeric(1), f = contact_at, parts = parts, spec = sp)
      cat(sprintf("  %s = list(contact_intercepts = c(c20 = %.13f, c50 = %.13f, c80 = %.13f),\n    outcome_Q3 = %.13f),  # outcome shift %.4f\n",
                  nm, a_con[["c20"]], a_con[["c50"]], a_con[["c80"]], a_out,
                  a_out - out_int[["Q3"]]))
    }
  }
}
