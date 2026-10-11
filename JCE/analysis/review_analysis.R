# Correct the manuscript correlations without changing the frozen simulation.
# Rscript JCE/analysis/review_analysis.R
source("methodology/v2/summarize_convergence.R")

APPARENT_REVIEW <- c("auc_index_scored", "auc_model_scored", "log_or_top",
                     "log_or_top_target", "beta_within", "beta_between")

corrected_family_correlation <- function(d, apparent = APPARENT_REVIEW,
                                         targeting = NULL) {
  for (m in intersect(apparent, names(d))) {
    flag <- METRIC_META$needs[match(m, METRIC_META$metric)]
    if (!is.na(flag)) d[[m]][d[[flag]] %in% TRUE] <- NA_real_
  }
  family_correlation(d, apparent, targeting)
}

if (sys.nframe() == 0L) {
  run <- "outputs/methodology/production_v03"
  mf <- manifest_files(run)
  source_paths <- file.path("methodology/v2", names(mf$manifest$source_hashes))
  source_hashes <- setNames(vapply(source_paths, \(p) digest::digest(file = p, algo = "sha256"),
                                  character(1)), names(mf$manifest$source_hashes))
  stopifnot(identical(source_hashes, mf$manifest$source_hashes))
  stopifnot(identical(digest::digest(file = "methodology/v2/quadrant_cells.csv", algo = "sha256"),
                      mf$manifest$quadrant_cells))
  ids <- mf$manifest$grid$config_id[mf$manifest$grid$block == "core"]
  paths <- mf$paths[tools::file_path_sans_ext(basename(mf$paths)) %in% ids]
  targeting <- METRIC_META$metric[METRIC_META$family %in% c("targeting", "budget")]
  cols <- unique(c("config_id", "rep", "arm", "setting", "gate", "failed",
                   "separation_dev", "separation_target", "or_top", "or_top_target",
                   "N", "n_dev_used", "coverage", "b10_slots", "b10_used",
                   "cl10_slots", "cl10_used", "cl10_residents",
                   setdiff(APPARENT_REVIEW, c("log_or_top", "log_or_top_target")), targeting))
  raw_core <- bind_rows(lapply(paths, function(p) {
    x <- readRDS(p)
    stopifnot(identical(attr(x, "provenance")$signature, mf$manifest$signature))
    x[, cols]
  }))
  core <- raw_core |>
    mutate(cell = paste(config_id, setting, sep = "|"),
           log_or_top = log(or_top), log_or_top_target = log(or_top_target))
  stopifnot(nrow(core) == 720000L, n_distinct(core$cell) == 72L, !any(core$failed))

  original <- family_correlation(core, APPARENT_REVIEW, targeting)
  corrected <- corrected_family_correlation(core, APPARENT_REVIEW, targeting)
  cached <- readRDS(paste0(run, "_displays.rds"))$families |>
    filter(a %in% APPARENT_REVIEW) |> arrange(a, t)
  stopifnot(isTRUE(all.equal(original |> arrange(a, t), cached, tolerance = 1e-12)))

  exclusions <- bind_rows(lapply(APPARENT_REVIEW, function(m) {
    flag <- METRIC_META$needs[match(m, METRIC_META$metric)]
    tibble(metric = m, validity_flag = flag,
           flagged_rows = sum(core[[flag]] %in% TRUE),
           finite_flagged_rows = sum(core[[flag]] %in% TRUE & is.finite(core[[m]])))
  }))
  valid <- core
  for (m in APPARENT_REVIEW) {
    flag <- METRIC_META$needs[match(m, METRIC_META$metric)]
    valid[[m]][valid[[flag]] %in% TRUE] <- NA_real_
  }
  principal <- c("cluster_rho_all", "cl10_misalloc", "b10_hn_capture")
  method_means <- valid |> group_by(cell, arm) |>
    summarise(across(all_of(c(APPARENT_REVIEW, principal)), \(x) mean(x, na.rm = TRUE)),
              failed = FALSE, .groups = "drop")
  sensitivity <- family_correlation(method_means, APPARENT_REVIEW, principal)
  separation <- core |> filter(arm == "complete_case") |>
    distinct(config_id, rep, gate, separation_dev, n_dev_used) |>
    group_by(gate) |> summarise(replicates = n(), separated = sum(separation_dev),
                                fraction = mean(separation_dev), n_dev = mean(n_dev_used),
                                .groups = "drop")
  policies <- core |> group_by(gate, setting, arm) |>
    summarise(across(c(cl10_misalloc, cl10_precision, cl10_hn_reached), mean),
              individual_budget_used = mean(b10_used / b10_slots),
              cluster_budget_used = mean(cl10_used / cl10_slots),
              residents_in_selected_clusters = mean(cl10_residents / N),
              coverage = mean(coverage), .groups = "drop")
  sample_sizes <- core |> distinct(config_id, rep, setting, N) |>
    summarise(min_target = min(N), mean_target = mean(N), max_target = max(N))

  source("methodology/v2/harmonize_donors.R")
  source("methodology/v2/quadrant_populations.R")
  source("methodology/v2/pilot_v2_config.R")
  cache_path <- "outputs/methodology/donors_harmonized_geo.rds"
  stopifnot(identical(digest::digest(file = cache_path, algo = "sha256"), mf$manifest$donor_cache))
  harmonized <- lapply(readRDS(cache_path), prefix_levels)
  pools <- build_quadrants(harmonized)
  pools$bd_allage <- build_age_pool(harmonized, "bangladesh", c(18L, Inf))
  for (nm in names(DONOR_POPULATIONS)) {
    pools[[nm]] <- build_age_pool(harmonized, DONOR_POPULATIONS[[nm]], c(18L, 49L))
  }
  donor_counts <- bind_rows(lapply(names(pools), function(nm) {
    d <- pools[[nm]]
    tibble(population = nm, people = nrow(d), clusters = n_distinct(d$cluster),
           min_age = min(d$age), max_age = max(d$age), women = sum(d$sex == "female"),
           men = sum(d$sex == "male"))
  }))
  pool_composition <- bind_rows(lapply(names(pools), function(nm) {
    pools[[nm]] |> count(donor, name = "people") |> mutate(population = nm)
  }))
  donor_flow <- bind_rows(lapply(names(harmonized), function(nm) {
    d <- harmonized[[nm]] |> filter(age >= 18, age <= 49)
    tibble(donor = nm, age_eligible = nrow(d),
           six_inputs_complete = sum(complete.cases(d[, INDEX_INPUTS])))
  }))
  comparison <- original |> select(a, t, original = median_rho) |>
    left_join(corrected |> select(a, t, corrected = median_rho), by = c("a", "t")) |>
    mutate(change = corrected - original)
  result <- list(signature = mf$manifest$signature, source_hashes = source_hashes,
                 families = corrected, original_families = original, comparison = comparison,
                 exclusions = exclusions, method_mean_sensitivity = sensitivity,
                 separation = separation, policies = policies, sample_sizes = sample_sizes,
                 donor_counts = donor_counts, donor_flow = donor_flow,
                 pool_composition = pool_composition)
  saveRDS(result, "JCE/analysis/review_analysis.rds", version = 3)
  write.csv(comparison, "JCE/analysis/correlation_correction.csv", row.names = FALSE)
  write.csv(exclusions, "JCE/analysis/correlation_exclusions.csv", row.names = FALSE)
  cat("Core rows:", nrow(core), "  corrected pairs:", nrow(corrected),
      "  negative medians:", sum(corrected$median_rho < 0), "\n")
  cat("Corrected median range:", range(corrected$median_rho), "\n")
  cat("Largest median change:", max(abs(comparison$change)), "\n")
  cat("Frozen source hashes match production; wrote JCE review summaries.\n")
}
