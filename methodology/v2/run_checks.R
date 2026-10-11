# Scientific readiness checks on a v03 run directory, beyond the contract suite.
#
#   Rscript methodology/v2/run_checks.R <run directory>
#
# A green test suite says the code does what it was written to do. This says the
# run realized the mechanisms it was supposed to realize, that the contrasts the
# design depends on are actually paired, and how many replicates the production
# run needs to reach the target Monte Carlo error.

suppressPackageStartupMessages({library(dplyr); library(tidyr)})

root <- c(".", "..", "../..") |>
  Filter(f = \(p) file.exists(file.path(p, "methodology/v2/plasmode_generate.R"))) |>
  (\(x) x[1])()
source(file.path(root, "methodology/v2/summarize_convergence.R"))

# Target Monte Carlo standard errors for the paired differences, in the units of
# each metric. Fixed before production: a margin of .01 on capture needs an MCSE
# well below it for an equivalence statement to be decidable, and a fifth of the
# margin is the rule used here.
TARGET_MCSE <- c(cluster_rho_all = 0.005, cl10_misalloc = 0.005,
                 b10_hn_capture = 0.002, capture_all = 0.002)

if (sys.nframe() == 0L) {
  dir <- commandArgs(trailingOnly = TRUE)[1]
  d <- read_run(dir)
  m <- attr(d, "manifest")

  cat("\n== completeness\n")
  reps <- d |> distinct(config_id, rep) |> count(config_id, name = "reps")
  cat("configurations on disk:", nrow(reps), "of", nrow(m$grid), "\n")
  cat("replicates per configuration:", paste(range(reps$reps), collapse = " to "), "\n")
  cat("rows:", nrow(d), "| failed:", sum(d$failed), "\n")
  if (sum(d$failed)) {
    print(as.data.frame(d |> filter(failed) |> count(config_id, arm, stage, name = "n")),
          row.names = FALSE)
  }

  # The gate is meant to change what is recorded, not who is in the population.
  # Only cells that ran under both gates can be compared.
  cat("\n== the gate contrast is paired\n")
  pair <- d |>
    filter(block %in% c("core", "donor", "stress"), !failed, arm == "oracle", setting == "c50") |>
    mutate(stem = sub("_g(con|dia)_", "_", config_id)) |>
    select(stem, gate, rep, N, n_hn, outcome_rate, capture_all)
  both <- pair |> count(stem, gate) |> count(stem) |> filter(n == 2L)
  pair <- pair |> filter(stem %in% both$stem) |>
    pivot_wider(names_from = gate, values_from = c(N, n_hn, outcome_rate, capture_all))
  if (!nrow(pair)) {
    cat("no cell ran under both gates; nothing to compare\n")
  } else {
    cat("cells compared:", nrow(pair), "\n")
    cat("identical target size:   ", isTRUE(all.equal(pair$N_contact, pair$N_diagnosis)), "\n")
    cat("identical high-need set: ", isTRUE(all.equal(pair$n_hn_contact, pair$n_hn_diagnosis)), "\n")
    cat("identical prevalence:    ", isTRUE(all.equal(pair$outcome_rate_contact, pair$outcome_rate_diagnosis)), "\n")
    cat("identical advice-free capture:", isTRUE(all.equal(pair$capture_all_contact, pair$capture_all_diagnosis)), "\n")
  }

  # The proxy control is left out: its items are asked of everyone by
  # construction, so pooling it would make the contact gate's eligibility read
  # above its contact rate.
  cat("\n== realized rates, by gate and setting (excluding the proxy control)\n")
  print(as.data.frame(d |>
    filter(variant != "proxy_block") |>
    distinct(config_id, rep, gate, setting, contact_rate, recognized_rate,
             eligible_rate, outcome_rate, unrecognized_case_rate) |>
    group_by(gate, setting) |>
    summarise(across(ends_with("rate"), mean), .groups = "drop")), digits = 4, row.names = FALSE)

  cat("\n== coverage, tiers, level substitutions and development size\n")
  print(as.data.frame(d |> filter(!failed) |> group_by(gate, arm) |>
    summarise(coverage = mean(coverage), n_tiers = mean(n_tiers),
              flag_rate = mean(tier_flag_rate), fixups = mean(level_fixups),
              n_dev = mean(n_dev_used), .groups = "drop") |>
    arrange(gate, arm)), digits = 3, row.names = FALSE)

  cat("\n== separation and unusable metrics\n")
  print(as.data.frame(d |> filter(!failed) |> group_by(gate, arm) |>
    summarise(sep_dev = mean(separation_dev), sep_target = mean(separation_target, na.rm = TRUE),
              or_missing = mean(is.na(or_top)), auc_missing = mean(is.na(auc_model_scored)),
              rho_missing = mean(is.na(cluster_rho_all)), .groups = "drop")),
    digits = 3, row.names = FALSE)

  cat("\n== negative controls: worst spread across the access grid in any replicate\n")
  print(as.data.frame(d |> filter(!failed, arm %in% ADVICE_FREE_ARMS_V03) |>
    group_by(config_id, rep, arm) |>
    summarise(rho = diff(range(cluster_rho_all)), capture = diff(range(capture_all)),
              budget = diff(range(b10_hn_capture)), .groups = "drop") |>
    group_by(arm) |>
    summarise(across(c(rho, capture, budget), max), .groups = "drop")), row.names = FALSE)

  cat("\n== paired precision on the core grid\n")
  for (metric in names(TARGET_MCSE)) {
    ct <- cellwise_contrasts(d |> filter(block == "core"), metric)
    if (!nrow(ct)) next
    need <- required_reps(ct, TARGET_MCSE[[metric]])
    cat(sprintf("  %-16s worst paired sd %.4f | worst MCSE here %.4f | target %.3f | needs %d replicates (%s vs %s)\n",
                metric, max(ct$mcse * sqrt(ct$n), na.rm = TRUE), max(ct$mcse, na.rm = TRUE),
                TARGET_MCSE[[metric]], need$need[1], need$better[1], need$worse[1]))
  }

  cat("\n== runtime\n")
  per <- d |> distinct(config_id, rep, elapsed)
  cat("processor seconds per replicate: mean", round(mean(per$elapsed), 2),
      "| max", round(max(per$elapsed), 2), "\n")
  cat("processor hours in this run:", round(sum(per$elapsed) / 3600, 2), "\n")
}
