# Summaries of the production-linked theory run (theory_production_link.R).
#
# Joins the development diagnostics of each regenerated production replicate to
# the production outcomes of the same replicate, and pairs the missing fuzzy average
# coding with the ten production methods under the frozen claim rule and
# margins, used read-only. Everything about the new coding is post hoc.
#
#   Rscript JCE/theory/theory_link_summary.R outputs/methodology/theory_link_r100

suppressPackageStartupMessages({library(dplyr); library(tidyr)})
source("methodology/v2/summarize_convergence.R")

link_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(link_dir)) link_dir <- "outputs/methodology/theory_link_r100"
prod_dir <- "outputs/methodology/production_v03"
signature <- readRDS(file.path(prod_dir, "manifest.rds"))$signature

L <- lapply(list.files(link_dir, pattern = "[.]rds$", full.names = TRUE), readRDS)
stopifnot(length(L) > 0, all(vapply(L, \(x) identical(x$manifest_signature, signature), logical(1))),
          all(vapply(L, \(x) length(x$errors) == 0L, logical(1))))
diag <- bind_rows(lapply(L, `[[`, "diag"))
ep <- bind_rows(lapply(L, `[[`, "metrics"))
# Runs before 8 October named the coding "eligible_profile".
ep$arm[ep$arm == "eligible_profile"] <- "fuzzy_average"
diag$coding[diag$coding == "eligible_profile"] <- "fuzzy_average"
n_rep <- max(ep$rep)

prod <- bind_rows(lapply(unique(ep$config_id), \(cid) {
  x <- readRDS(file.path(prod_dir, paste0(cid, ".rds")))
  x[x$rep <= n_rep, ]
}))
info <- distinct(prod, config_id, block, variant, gate, population)
d <- bind_rows(prod, ep |> left_join(info, by = "config_id")) |>
  mutate(arm = factor(arm, levels = c(ARM_ORDER, "fuzzy_average")),
         cell = paste(config_id, setting, sep = "|"))
num <- \(x, k = 3) round(mean(x, na.rm = TRUE), k)
show <- \(x) print(as.data.frame(x), row.names = FALSE)
CODING_ARM <- c(category = "category", subset = "passive_mca", zero = "zero")

cat("Production replicates per configuration:", n_rep, " configurations:", length(L), "\n")

cat("\n== Theorem 1: filter loading equals the chi-square formula ==\n")
cat("max |phi - formula| over", nrow(diag), "rows:", signif(max(abs(diag$phi - diag$phi_formula), na.rm = TRUE), 2), "\n")
ord <- diag |> filter(coding %in% c("zero", "subset")) |> select(config_id, rep, coding, phi) |>
  pivot_wider(names_from = coding, values_from = phi)
cat("Replicates where skip-as-no carries less filter than subset MCA:",
    sum(ord$zero < ord$subset), "of", nrow(ord), "\n")

# Targeting of each skip coding on the same replicates, joined to its diagnostics
tau <- d |> filter(!failed) |> select(config_id, rep, setting, arm, cluster_rho_all) |>
  pivot_wider(names_from = arm, values_from = cluster_rho_all)
cap <- diag |> filter(coding %in% names(CODING_ARM)) |>
  left_join(tau |> filter(setting == "c50"), by = c("config_id", "rep")) |>
  mutate(tau_coding = case_when(coding == "category" ~ category, coding == "subset" ~ passive_mca,
                                TRUE ~ zero),
         loss = drop_items - tau_coding) |>
  left_join(info, by = "config_id")
by_cfg <- cap |> group_by(block, variant, gate, population, config_id, coding) |>
  summarise(kappa = num(kappa), block_share = num(block_share), filter_share = num(filter_share),
            loss = num(loss), tau_sd = round(sd(tau_coding), 3), gap = num(eig1 - eig2), .groups = "drop")

cat("\n== Proposition 2: capture threshold (configuration means, 50 percent contact) ==\n")
show(by_cfg |> mutate(side = ifelse(kappa > 1, "kappa > 1", "kappa <= 1")) |> group_by(coding, side) |>
       summarise(configs = n(), block_share = num(block_share), loss = num(loss), .groups = "drop"))
cat("Spearman, kappa against loss over configuration-coding pairs:",
    round(stats::cor(by_cfg$kappa, by_cfg$loss, method = "spearman"), 3),
    "\nSpearman, block share against loss:",
    round(stats::cor(by_cfg$block_share, by_cfg$loss, method = "spearman"), 3), "\n")
cat("Replicate level, kappa against loss: Spearman",
    round(stats::cor(cap$kappa, cap$loss, method = "spearman", use = "complete"), 3), "\n")

cat("\n== Proposition 4: replicate SD of targeting against distance from the threshold ==\n")
show(by_cfg |> mutate(band = cut(abs(kappa - 1), c(0, .15, .3, .6, Inf), include.lowest = TRUE)) |>
       group_by(band) |> summarise(pairs = n(), tau_sd = num(tau_sd), gap = num(gap), .groups = "drop"))
cat("Spearman, |kappa - 1| against SD:",
    round(stats::cor(abs(by_cfg$kappa - 1), by_cfg$tau_sd, method = "spearman"), 3), "\n")

cat("\n== Proposition 6: budget saturation under the diagnosis filter ==\n")
show(d |> filter(gate == "diagnosis", arm %in% c("zero", "category", "passive_mca", "drop_items")) |>
       mutate(recognized_minus_budget = round(recognized_rate - 0.10, 2) > 0) |>
       group_by(setting, recognized_above_budget = recognized_minus_budget, arm) |>
       summarise(recognized_rate = num(recognized_rate), unrec_capture = num(b10_unrec_capture, 5),
                 .groups = "drop"))

cat("\n== Remedy: missing fuzzy average against every production method (post hoc) ==\n")
show(d |> filter(!failed, block == "core") |> group_by(gate, arm) |>
       summarise(cluster_rho = num(cluster_rho_all), inclusion_error = num(cl10_misalloc, 4),
                 hn_capture = num(b10_hn_capture), unrec_capture = num(b10_unrec_capture), .groups = "drop") |>
       arrange(gate, desc(cluster_rho)))
show(d |> filter(!failed, block != "core", arm %in% c("fuzzy_average", "drop_items", "impute_adaptive",
                                                       "impute_frozen", "simple_score")) |>
       group_by(variant, gate, arm) |> summarise(cluster_rho = num(cluster_rho_all), .groups = "drop") |>
       pivot_wider(names_from = arm, values_from = cluster_rho))

core <- d |> filter(block == "core")
others <- setdiff(levels(droplevels(core$arm)), "fuzzy_average")
claims <- bind_rows(lapply(c("cluster_rho_all", "cl10_misalloc", "b10_hn_capture"), \(m)
  classify_pairs(cellwise_contrasts(core, m, pairs = lapply(others, \(a) c("fuzzy_average", a)),
                                    cells = unique(core$cell)), margin = MARGINS[[m]])))
cat("\nFrozen rule, 72 core cells,", n_rep, "replicates each (preliminary below 1,000):\n")
show(claims |> select(metric, worse, claim, min_diff, max_diff, worst_mcse))

# Aggregates for the manuscript displays. Production arms are recomputed from the
# same rows so they can be checked against the production displays file.
METRIC_COLS <- c("cluster_rho_all", "cl10_misalloc", "cl10_precision", "b10_hn_capture",
                 "b10_unrec_capture")
cellwise <- core |> filter(!failed) |>
  group_by(gate, population, config_id, setting, arm) |>
  summarise(across(all_of(METRIC_COLS), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
by_setting <- cellwise |> group_by(gate, setting, arm) |>
  summarise(across(all_of(METRIC_COLS), mean), .groups = "drop")

# The theory against the production replicates, by coding and filter, over the
# core configurations. Loss is against dropping the items at 50 percent contact,
# the development level; for missing fuzzy average kappa is the corollary's bound.
P_INPUTS <- 6L
stopifnot(isTRUE(all.equal(cap$kappa, cap$k * cap$lamA / (P_INPUTS * cap$lamX))))
fa_loss <- tau |> filter(setting == "c50") |>
  transmute(config_id, rep, loss = drop_items - fuzzy_average)
theory_table <- bind_rows(
  cap |> filter(block == "core") |> group_by(gate, coding) |>
    summarise(p_dev = mean(p_dev), phi = mean(phi), phi_formula = mean(phi_formula),
              kappa = mean(kappa), block_share = mean(block_share),
              filter_share = mean(filter_share), loss_sd = sd(loss), loss = mean(loss),
              .groups = "drop"),
  diag |> filter(coding == "fuzzy_average") |> left_join(info, by = "config_id") |>
    filter(block == "core") |> left_join(fa_loss, by = c("config_id", "rep")) |>
    group_by(gate, coding) |>
    summarise(p_dev = mean(p_dev), phi = 0, phi_formula = 0,
              kappa_bound = mean(k * p_dev / (P_INPUTS * lamX)),
              filter_share = mean(filter_share), loss_sd = sd(loss), loss = mean(loss),
              .groups = "drop"))
cat("\n== Theory against the core production replicates ==\n")
show(theory_table |> mutate(across(where(is.numeric), \(x) round(x, 3))))

# Missing fuzzy average in the 45 population, filter and contact combinations of
# the nine-pool comparison (reference scenario for the two core populations).
pools_fa <- d |>
  filter(!failed, arm == "fuzzy_average",
         block %in% c("quadrant", "donor") | (block == "core" & grepl("_r60_q00_", config_id))) |>
  group_by(population, gate, setting) |>
  summarise(cluster_rho_all = mean(cluster_rho_all), .groups = "drop")

# Supplement A6.3, A6.6 and A6.9 (added 19 October 2026): the inertia ratio and
# instability within each coding, the exact eigenvalue of the filtered block,
# budget cells above the share diagnosed, and missing fuzzy average with eight
# advice items.
cat("\n== Within codings, exact block eigenvalues, budget cells ==\n")
show(by_cfg |> group_by(coding) |>
       summarise(pairs = n(), spearman_kappa_loss = round(stats::cor(kappa, loss, method = "spearman"), 2),
                 sd_within_.3 = num(tau_sd[abs(kappa - 1) <= .3]), sd_beyond_.3 = num(tau_sd[abs(kappa - 1) > .3]),
                 min_kappa = round(min(kappa), 2), .groups = "drop"))
show(cap |> group_by(coding, gate) |>
       summarise(samples = n(), lamA_is_1 = mean(abs(lamA - 1) < 1e-9),
                 lamA_is_1_minus_p = num(abs(lamA - (1 - p_dev)) < 1e-9),
                 lamA_at_least_1_minus_p = mean(lamA >= 1 - p_dev - 1e-9), .groups = "drop"))
show(d |> filter(!failed, gate == "diagnosis", arm %in% c("zero", "category", "passive_mca", "drop_items")) |>
       group_by(config_id, setting, arm) |>
       summarise(recognized_rate = mean(recognized_rate), unrec = mean(b10_unrec_capture), .groups = "drop") |>
       filter(recognized_rate > 0.10) |>
       group_by(arm) |>
       summarise(cells_above_budget = n(), settings = paste(unique(setting), collapse = " "),
                 mean = signif(mean(unrec), 2), max = signif(max(unrec), 2), .groups = "drop"))
show(d |> filter(!failed, variant == "items_eight", arm %in% c("fuzzy_average", "drop_items")) |>
       group_by(arm, setting) |> summarise(hn_capture = num(b10_hn_capture), .groups = "drop"))

saveRDS(list(n_rep = n_rep, by_cfg = by_cfg, cap = cap, claims = claims, cellwise = cellwise,
             by_setting = by_setting, theory_table = theory_table, pools_fa = pools_fa),
        file.path("JCE/theory", sprintf("theory_link_summary_r%d.rds", n_rep)))
