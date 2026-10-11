# Checks of the mechanisms and item-model sizes quoted in supplement A3 and A4.
#
# Part 1 regenerates fresh replicates with the production generator and arms to
# measure what each MCA axis tracks, the realized access-need association and
# the size of the imputation models. Part 2 reads the saved core production
# rows and the claim table. Nothing in methodology/v2 or outputs/ is written.
#
#   Rscript JCE/analysis/editor_review_checks.R          # 25 fresh replicates
#   NREP=5 Rscript JCE/analysis/editor_review_checks.R   # quicker

suppressPackageStartupMessages({library(dplyr); library(tidyr)})
source("methodology/v2/harmonize_donors.R")
source("methodology/v2/plasmode_generate.R")
source("methodology/v2/index_arms.R")
source("methodology/v2/quadrant_populations.R")
source("methodology/v2/pilot_v2_config.R")

# ---- Part 1: mechanisms in fresh replicates ------------------------------

harmonized <- lapply(readRDS("outputs/methodology/donors_harmonized_geo.rds"), prefix_levels)
keep_cols <- function(d) droplevels(as.data.frame(d[, c(INDEX_INPUTS, "age", "sex", "cluster")]))
pools <- lapply(build_quadrants(harmonized, quadrant_cells("methodology/v2/quadrant_cells.csv")), keep_cols)
pools$bd_allage <- keep_cols(build_age_pool(harmonized, "bangladesh", c(18L, Inf)))

cluster_cor <- function(t, x) {
  stats::cor(tapply(t[[x]], t$cluster, mean), tapply(t$V, t$cluster, mean), method = "spearman")
}

axis_profile <- function(dev, how) {
  d <- code_items(dev, if (how == "subset") "category" else how)
  mp <- mca_dev(d, c(INDEX_INPUTS, advice_names(dev)),
                excl = if (how == "subset") notasked_levels(dev))
  g <- observed_rows(dev)
  c(eig1 = unname(mp$fit$eig[1, 1]),
    cor_G = abs(stats::cor(mp$coord_dev, g)),
    cor_V = abs(stats::cor(mp$coord_dev, dev$V)))
}

probe_one <- function(pop, gate, r) {
  spec <- pilot_spec(population = pop, gate = gate, r2 = 0.6)
  set.seed(1000 * r + nchar(pop))
  sim <- generate_replicate(pools[[pop]], spec, 1000L, 20000L, signal_scale(pools[[pop]], spec))
  dev <- sim$development
  g <- observed_rows(dev)
  obs <- dev[g, ]
  rhs <- varying(obs, c(INDEX_INPUTS, "age", "sex"))
  n_par <- length(stats::coef(suppressWarnings(
    stats::glm(stats::reformulate(rhs, "advice1"), data = obs, family = stats::binomial))))
  min_class <- min(vapply(advice_names(dev), \(v) min(table(factor(obs[[v]], levels = 0:1))), numeric(1)))
  sub <- mca_dev(code_items(dev, "category"), c(INDEX_INPUTS, advice_names(dev)),
                 excl = notasked_levels(dev))
  w <- setNames(sub$proj$w, sub$proj$cols)
  same_sign <- mean(vapply(advice_names(dev), \(v) sign(w[[paste0(v, "_0")]]) == sign(w[[paste0(v, "_1")]]),
                           logical(1)))
  tibble(pop, gate, rep = r,
         contact_need_c20 = cluster_cor(sim$targets$c20, "contacted"),
         contact_need_c50 = cluster_cor(sim$targets$c50, "contacted"),
         contact_need_c80 = cluster_cor(sim$targets$c80, "contacted"),
         eligible_need_c20 = cluster_cor(sim$targets$c20, "eligible"),
         eligible_need_c80 = cluster_cor(sim$targets$c80, "eligible"),
         recognized_c80 = mean(sim$targets$c80$recognized),
         !!!setNames(as.list(axis_profile(dev, "zero")), paste0(c("eig1", "corG", "corV"), "_no")),
         !!!setNames(as.list(axis_profile(dev, "category")), paste0(c("eig1", "corG", "corV"), "_category")),
         !!!setNames(as.list(axis_profile(dev, "subset")), paste0(c("eig1", "corG", "corV"), "_subset")),
         eig1_dropped = unname(mca_dev(code_items(dev, "observed"), INDEX_INPUTS)$fit$eig[1, 1]),
         subset_yes_no_same_sign = same_sign,
         dev_eligible = sum(g), item_model_parameters = n_par, smaller_class = min_class,
         events_per_parameter = min_class / (n_par - 1))
}

n_rep <- as.integer(Sys.getenv("NREP", "25"))
probe <- bind_rows(lapply(seq_len(n_rep), \(r) bind_rows(
  lapply(c("Q3", "bd_allage"), \(p) bind_rows(lapply(GATES, \(g) probe_one(p, g, r)))))))
probe_summary <- probe |> select(-rep) |> group_by(pop, gate) |>
  summarise(across(everything(), \(x) round(mean(x), 3)), .groups = "drop")

# ---- Part 2: saved production rows ---------------------------------------

files <- list.files("outputs/methodology/production_v03", pattern = "_base[.]rds$", full.names = TRUE)
files <- files[grepl("^(Q3|bd_allage)_nd1000", basename(files))]
rows <- bind_rows(lapply(files, \(p) {
  x <- readRDS(p)
  x[!x$failed, c("arm", "setting", "gate", "rep", "config_id", "cluster_rho_all",
                 "b10_hn_capture", "cl10_misalloc")]
}))
stopifnot(length(files) == 24L, nrow(rows) == 720000L)

# Arms that read no filtered item are identical in every cell of a replicate.
invariant_spread <- rows |> filter(arm %in% ADVICE_FREE_ARMS) |>
  group_by(arm, rep, cfg = sub("_g(con|dia)_", "_", config_id)) |>
  summarise(spread = diff(range(b10_hn_capture)) + diff(range(cluster_rho_all)), .groups = "drop") |>
  group_by(arm) |> summarise(max_spread = max(spread))

# The matched setting, target contact equal to development contact, against the
# shifted ones.
matched_worst <- rows |> group_by(gate, config_id, arm, setting) |>
  summarise(r = mean(cluster_rho_all, na.rm = TRUE), .groups = "drop") |>
  group_by(gate, config_id, arm) |>
  summarise(c50_worst = setting[which.min(r)] == "c50", .groups = "drop") |>
  group_by(gate, arm) |> summarise(configs = n(), c50_worst = sum(c50_worst), .groups = "drop") |>
  filter(arm %in% c("zero", "category", "passive_mca"))

# Imputers against dropping under the diagnosis filter, replicate by replicate.
paired <- rows |> filter(gate == "diagnosis") |>
  select(config_id, rep, setting, arm, cluster_rho_all) |>
  pivot_wider(names_from = arm, values_from = cluster_rho_all)
imputer_tail <- bind_rows(lapply(c("impute_frozen", "impute_adaptive"), \(a) {
  dif <- paired[[a]] - paired$drop_items
  tibble(arm = a, mean_diff = mean(dif, na.rm = TRUE), median_diff = median(dif, na.rm = TRUE),
         share_beating_drop = mean(dif > 0, na.rm = TRUE),
         share_losing_by_10 = mean(dif < -0.10, na.rm = TRUE),
         median_c20 = median(dif[paired$setting == "c20"], na.rm = TRUE),
         median_c80 = median(dif[paired$setting == "c80"], na.rm = TRUE))
}))

inclusion_floor <- rows |> filter(gate == "contact", setting == "c20") |>
  group_by(arm) |> summarise(mean_inclusion_error = mean(cl10_misalloc, na.rm = TRUE))

claims <- read.csv("outputs/methodology/production_v03_claims.csv")
margin <- c(cluster_rho_all = 0.02, cl10_misalloc = 0.02, b10_hn_capture = 0.01)
robust <- claims[grepl("^robust", claims$claim), ]
robust$sub_margin_cell <- pmin(abs(robust$min_diff), abs(robust$max_diff)) < margin[robust$metric]
robust_small <- robust |> group_by(metric) |>
  summarise(robust = n(), with_sub_margin_cell = sum(sub_margin_cell), .groups = "drop")

out <- list(probe = probe, probe_summary = probe_summary, invariant_spread = invariant_spread,
            matched_worst = matched_worst, imputer_tail = imputer_tail,
            inclusion_floor = inclusion_floor, robust_small = robust_small)
saveRDS(out, "JCE/analysis/editor_review_checks.rds")
for (nm in setdiff(names(out), "probe")) {
  cat("\n==", nm, "==\n")
  print(as.data.frame(out[[nm]]), row.names = FALSE)
}
