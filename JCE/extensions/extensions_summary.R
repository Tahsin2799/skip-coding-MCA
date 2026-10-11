# Summaries of the post hoc extensions (run_extensions.R output).
#
#   Rscript JCE/extensions/extensions_summary.R outputs/methodology/extensions_v1
#
# Writes JCE/extensions/extensions_summary.rds, or the path given second. Every contrast is paired within
# replicate and contact level, with its Monte Carlo standard error; these are
# descriptive sensitivity analyses, not the frozen decision rule.

suppressPackageStartupMessages({library(dplyr); library(tidyr)})

dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(dir)) dir <- "outputs/methodology/extensions_v1"
out_file <- commandArgs(trailingOnly = TRUE)[2]
if (is.na(out_file)) out_file <- "JCE/extensions/extensions_summary.rds"
files <- list.files(dir, pattern = "[.]rds$", full.names = TRUE)
stopifnot(length(files) > 0)
L <- lapply(files, readRDS)
stopifnot(all(vapply(L, \(x) length(x$errors) == 0L, logical(1))))
cat("FactoMineR by pass:", paste(unique(vapply(L, `[[`, "", "factominer")), collapse = ", "), "\n")

parse_id <- function(id) {
  m <- regmatches(id, regexec("^ext_([a-z]+)_l([0-9]{3})_nd([0-9]+)_dev(c[0-9]{2})_g([a-z]{3})$", id))[[1]]
  tibble(config_id = id, study = m[2], loading = as.numeric(m[3]) / 100, n_dev = as.integer(m[4]),
         dev_setting = m[5], gate = c(con = "contact", dia = "diagnosis")[m[6]])
}
metrics <- bind_rows(lapply(L, `[[`, "metrics"))
info <- bind_rows(lapply(unique(metrics$config_id), parse_id))
metrics <- metrics |> left_join(info, by = "config_id")
diag_mca <- bind_rows(lapply(L, `[[`, "diag_mca")) |> left_join(info, by = "config_id")
diag_pca <- bind_rows(lapply(L, `[[`, "diag_pca")) |> left_join(info, by = "config_id")
size <- bind_rows(lapply(L, `[[`, "size")) |> left_join(info, by = "config_id")
reps <- metrics |> distinct(config_id, arm, rep) |> count(config_id, arm, name = "replicates")
cat("Replicates per configuration and arm:", paste(range(reps$replicates), collapse = " to "), "\n")
cat("Failed rows:", sum(metrics$failed), "of", nrow(metrics), "\n")

num <- \(x, k = 3) round(x, k)
show <- \(x) print(as.data.frame(x), row.names = FALSE)

# Paired difference of every arm against dropping the items, per configuration
# and contact level.
paired <- function(d, metric = "cluster_rho_all", ref = "drop_items") {
  base <- d |> filter(arm == ref, !failed) |> select(config_id, rep, setting, ref_value = all_of(metric))
  d |> filter(!failed) |> select(config_id, rep, setting, arm, value = all_of(metric)) |>
    inner_join(base, by = c("config_id", "rep", "setting")) |>
    mutate(diff = value - ref_value)
}
summ <- function(x) {
  # The paired spread and tail are taken before the mean replaces `diff`.
  x |> summarise(mcse = sd(diff) / sqrt(n()), median_diff = median(diff), behind_10 = mean(diff < -0.10),
                 n = n(), mean = mean(value), diff = mean(diff), .groups = "drop")
}
pr <- paired(metrics) |> left_join(info, by = "config_id")

# Sweep (review B4): where does each method overtake dropping as the advice
# items carry more information about need?
sweep <- pr |> filter(study == "sweep") |> group_by(gate, loading, arm) |> summ()
cat("\n== Sweep: cluster rank correlation and paired difference from dropping, pooled over contact ==\n")
show(sweep |> filter(arm %in% c("drop_items", "impute_adaptive", "impute_frozen", "impute_flic_frozen",
                                "impute_flic_adaptive", "fuzzy_average", "rimca", "simple_score", "oracle")) |>
       transmute(gate, loading, arm, mean = num(mean), diff = num(diff), mcse = num(mcse, 4)) |>
       arrange(gate, arm, loading))
crossover <- sweep |> filter(arm != "drop_items") |> group_by(gate, arm) |> arrange(loading, .by_group = TRUE) |>
  summarise(crossover = {
    i <- which(diff[-1] > 0 & diff[-length(diff)] <= 0)
    if (length(i)) loading[i[1]] + (0 - diff[i[1]]) * (loading[i[1] + 1] - loading[i[1]]) / (diff[i[1] + 1] - diff[i[1]])
    else NA_real_
  }, ahead_everywhere = all(diff > 0), .groups = "drop")
cat("\nLoading at which each method's mean difference from dropping turns positive (linear interpolation):\n")
show(crossover |> mutate(crossover = num(crossover)))
# coding_diagnostics() (JCE/theory/filter_capture.R) names the inertia ratio
# `kappa` and the R-squared with eligibility `filter_share`.
kappa_sweep <- diag_mca |> filter(study == "sweep") |> group_by(gate, loading, coding) |>
  summarise(inertia_ratio = num(mean(kappa, na.rm = TRUE)), r2_eligibility = num(mean(filter_share)),
            .groups = "drop")

# Development size (A6): the imputers under the diagnosis filter.
devsize_cfg <- info |> filter(gate == "diagnosis", loading == 0.25, dev_setting == "c50",
                              study %in% c("sweep", "devsize"))
devsize <- pr |> filter(config_id %in% devsize_cfg$config_id,
                        arm %in% c("impute_adaptive", "impute_frozen", "impute_flic_adaptive",
                                   "impute_flic_frozen", "fuzzy_average", "rimca")) |>
  group_by(n_dev, arm) |> summ()
cat("\n== Development size, diagnosis filter: difference from dropping the items ==\n")
show(devsize |> mutate(across(c(mean, diff, median_diff, behind_10), num), mcse = num(mcse, 4)))
devsize_models <- size |> filter(config_id %in% devsize_cfg$config_id) |> group_by(n_dev) |>
  summarise(dev_eligible = num(mean(dev_eligible), 1), events_per_parameter = num(mean(events_per_parameter), 2),
            .groups = "drop")
show(devsize_models)

# Development contact (A1): is the loss of the codings that read the skip
# largest where development and target contact match, or where half the target
# is eligible?
devc_cfg <- info |> filter(loading == 0.25, n_dev == 1000L, study %in% c("sweep", "devcontact"))
devcontact <- metrics |> filter(config_id %in% devc_cfg$config_id, !failed,
                                arm %in% c("zero", "passive_mca", "category", "drop_items", "fuzzy_average")) |>
  group_by(gate, dev_setting, setting, arm) |>
  summarise(rho = mean(cluster_rho_all), mcse = sd(cluster_rho_all) / sqrt(n()), .groups = "drop")
cat("\n== Development contact by target contact: cluster rank correlation ==\n")
show(devcontact |> select(-mcse) |> mutate(rho = num(rho)) |>
       pivot_wider(names_from = setting, values_from = rho) |> arrange(gate, arm, dev_setting))
worst <- devcontact |> filter(arm %in% c("zero", "passive_mca", "category")) |>
  group_by(gate, dev_setting, arm) |> slice_min(rho, n = 1) |> select(gate, dev_setting, arm, worst_target = setting)
cat("\nTarget contact level at which each skip coding ranks worst, by development contact:\n")
show(worst)
# The development share asked sets each coding's squared correlation ratio, so
# the inertia ratio and the eligibility R-squared are read by development contact.
devcontact_theory <- diag_mca |> filter(config_id %in% devc_cfg$config_id) |>
  group_by(gate, dev_setting, coding) |>
  summarise(p = num(mean(p_dev)), inertia_ratio = num(mean(kappa, na.rm = TRUE), 2),
            r2_block = num(mean(block_share)), r2_eligibility = num(mean(filter_share)), .groups = "drop")
show(devcontact_theory)

# Principal component indices: the closed forms and the capture ratio.
cat("\n== Principal component indices: closed-form check and dominance, by filter ==\n")
# The closed forms are checked in every development sample; the means are taken
# over the sweep, where development is the reference (1,000 people, 50 percent
# contact), since development contact and size move p.
pca_formula <- diag_pca |> group_by(gate, coding) |>
  summarise(max_formula_diff = signif(max(eta2_max_diff), 2), samples = n(), .groups = "drop")
pca_theory_loading <- diag_pca |> filter(study == "sweep") |> group_by(gate, loading, coding) |>
  summarise(p = num(mean(p_dev)), eta2 = num(mean(eta2_mean)), inertia_ratio = num(mean(inertia_ratio)),
            r2_block = num(mean(r2_block)), r2_eligibility = num(mean(r2_eligibility)), .groups = "drop")
pca_theory <- diag_pca |> filter(study == "sweep") |> group_by(gate, coding) |>
  summarise(p = num(mean(p_dev)), eta2 = num(mean(eta2_mean)), inertia_ratio = num(mean(inertia_ratio)),
            r2_block = num(mean(r2_block)), r2_eligibility = num(mean(r2_eligibility)), .groups = "drop") |>
  left_join(pca_formula, by = c("gate", "coding"))
show(pca_theory)
show(pca_theory_loading)
pca_vs_mca <- metrics |> filter(!failed, study == "sweep") |>
  filter(arm %in% c("zero", "pca_no", "fuzzy_average", "pca_fuzzy", "category", "pca_category",
                    "drop_items", "pca_drop")) |>
  group_by(gate, loading, arm) |> summarise(rho = num(mean(cluster_rho_all)), .groups = "drop") |>
  pivot_wider(names_from = arm, values_from = rho)
cat("\nPCA against MCA counterparts, cluster rank correlation:\n")
show(pca_vs_mca)

saveRDS(list(info = info, sweep = sweep, crossover = crossover, kappa_sweep = kappa_sweep,
             devsize = devsize, devsize_models = devsize_models, devcontact = devcontact, worst = worst,
             devcontact_theory = devcontact_theory,
             pca_theory = pca_theory, pca_theory_loading = pca_theory_loading, pca_vs_mca = pca_vs_mca,
             replicates = reps, failed = sum(metrics$failed)),
        out_file)
