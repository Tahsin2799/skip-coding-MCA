# Turn a verified run directory into the tables the manuscript may quote.
#
#   Rscript methodology/v2/claim_ledger.R <run directory> [--margins]
#
# Writes, beside the run directory:
#   <dir>_claims.csv     every prespecified contrast, its effect, its
#                        uncertainty, the cells it holds in, and its verdict
#   <dir>_displays.rds   the data behind each main figure, one tibble per
#                        display, so no figure is ever drawn from a typed number
#
# Nothing here decides what is interesting. It applies the frozen claim rule to
# the frozen list of contrasts and reports what survives.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

root <- c(".", "..", "../..") |>
  Filter(f = \(p) file.exists(file.path(p, "methodology/v2/plasmode_generate.R"))) |>
  (\(x) x[1])()
source(file.path(root, "methodology/v2/summarize_convergence.R"))

# The contrasts the paper is allowed to make. Every practical arm against every
# other practical arm, plus each one against the two reference arms. Fixed here
# rather than chosen after seeing the run.
PRACTICAL <- c("zero", "category", "complete_case", "drop_items",
               "impute_frozen", "impute_adaptive", "passive_mca", "simple_score")
REFERENCE <- c("raw_items", "oracle")

prespecified_pairs <- function(arms) {
  p <- intersect(PRACTICAL, arms)
  r <- intersect(REFERENCE, arms)
  c(utils::combn(p, 2L, simplify = FALSE),
    unlist(lapply(r, \(x) lapply(p, \(y) c(y, x))), recursive = FALSE))
}

claim_table <- function(d, metrics = METRIC_META$metric[METRIC_META$primary],
                        margins = NULL) {
  core <- d |> filter(block == "core")
  cells <- unique(core$cell)
  pairs <- prespecified_pairs(levels(droplevels(core$arm)))
  bind_rows(lapply(metrics, function(m) {
    ct <- cellwise_contrasts(core, m, pairs = pairs, cells = cells)
    classify_pairs(ct, margin = if (!is.null(margins)) margins[[m]] else NULL) |>
      mutate(required_cells = length(cells))
  }))
}

# One tibble per display. Figures are drawn from these and from nothing else.
displays <- function(d) {
  core <- d |> filter(block == "core", !failed)
  list(
    # Figure 2: what a validation restricted to scoreable rows would report,
    # against the same arm read on the population it serves, with coverage.
    coverage_gap = core |>
      group_by(gate, setting, arm) |>
      summarise(coverage = mean(coverage),
                rho_scored = mean(cluster_rho_scored, na.rm = TRUE),
                rho_all = mean(cluster_rho_all, na.rm = TRUE),
                capture_scored = mean(capture_scored, na.rm = TRUE),
                capture_all = mean(capture_all, na.rm = TRUE),
                flag_rate = mean(tier_flag_rate, na.rm = TRUE),
                mcse_rho_all = stats::sd(cluster_rho_all, na.rm = TRUE) / sqrt(n()),
                .groups = "drop"),
    # Figure 3: who a programme reaches at a stated budget.
    budget = core |>
      select(gate, setting, arm, rep, config_id,
             b05_hn_capture, b10_hn_capture, b20_hn_capture,
             b10_unrec_capture, b10_hn_shortfall_per1000) |>
      pivot_longer(starts_with("b"), names_to = "metric", values_to = "value") |>
      group_by(gate, setting, arm, metric) |>
      summarise(mean = mean(value, na.rm = TRUE),
                mcse = stats::sd(value, na.rm = TRUE) / sqrt(sum(!is.na(value))),
                .groups = "drop"),
    # Figure 4: every core cell, so a reversal cannot hide inside an average.
    cellwise = core |>
      group_by(gate, population, r2, q1_age_slope, setting, arm) |>
      summarise(cluster_rho_all = mean(cluster_rho_all, na.rm = TRUE),
                cl10_misalloc = mean(cl10_misalloc, na.rm = TRUE),
                b10_hn_capture = mean(b10_hn_capture, na.rm = TRUE),
                .groups = "drop"),
    # Optional figure 5: where the apparent and the true families disagree.
    families = family_correlation(core),
    # Supplement: the stress battery and the mechanism controls against the
    # reference cell they were built around.
    variants = d |>
      filter(!failed, block %in% c("stress", "controls", "quadrant", "donor")) |>
      group_by(block, population, variant, gate, setting, arm) |>
      summarise(cluster_rho_all = mean(cluster_rho_all, na.rm = TRUE),
                b10_hn_capture = mean(b10_hn_capture, na.rm = TRUE),
                coverage = mean(coverage), .groups = "drop"),
    # Supplement: age bands, on whole-target denominators.
    age_bands = d |>
      filter(!failed, !is.na(n_hn_60p), n_hn_60p > 0) |>
      group_by(population, gate, q1_age_slope, setting, arm) |>
      summarise(across(starts_with("cap_tier_"), \(x) mean(x, na.rm = TRUE)),
                across(starts_with("cap_b10_"), \(x) mean(x, na.rm = TRUE)),
                .groups = "drop"),
    # The proxy control asks its items of everyone, so it is left out here as
    # in run_checks.R; pooled, the contact gate's eligibility reads above its
    # contact rate.
    rates = d |>
      filter(variant != "proxy_block") |>
      distinct(config_id, rep, setting, gate, contact_rate, recognized_rate,
               eligible_rate, outcome_rate, unrecognized_case_rate) |>
      group_by(gate, setting) |>
      summarise(across(ends_with("rate"), mean), .groups = "drop"))
}

if (sys.nframe() == 0L) {
  argv <- commandArgs(trailingOnly = TRUE)
  dir <- argv[1]
  margins <- if ("--margins" %in% argv) MARGINS else NULL
  d <- read_run(dir)
  cl <- claim_table(d, margins = margins)
  ds <- displays(d)
  utils::write.csv(cl, paste0(dir, "_claims.csv"), row.names = FALSE)
  saveRDS(ds, paste0(dir, "_displays.rds"), version = 3)

  cat("\nclaims written: ", paste0(dir, "_claims.csv"), "\n", sep = "")
  cat("displays written:", paste0(dir, "_displays.rds"), "\n")
  cat("\nverdicts by metric\n")
  print(as.data.frame(cl |> count(metric, claim)), row.names = FALSE)
  cat("\nwhat survives the claim rule\n")
  print(as.data.frame(cl |> filter(claim %in% c("robust_superiority", "equivalent")) |>
    arrange(metric, desc(min_diff)) |>
    select(metric, better, worse, claim, min_diff, max_diff, worst_mcse)),
    digits = 3, row.names = FALSE)
}
