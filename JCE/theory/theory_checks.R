# Numerical checks of the propositions of supplement A6 of the paper.
#
# Probe replicates from the unchanged production generator and arms, never
# production estimates. Writes only JCE/theory/theory_checks.rds.
#
#   Rscript JCE/theory/theory_checks.R            # 12 replicates per cell
#   NREP=3 Rscript JCE/theory/theory_checks.R     # smoke

suppressPackageStartupMessages({library(dplyr); library(tidyr); library(FactoMineR)})
source("methodology/v2/harmonize_donors.R")
source("methodology/v2/plasmode_generate.R")
source("methodology/v2/index_arms.R")
source("methodology/v2/quadrant_populations.R")
source("methodology/v2/pilot_v2_config.R")
source("JCE/theory/filter_capture.R")

harmonized <- lapply(readRDS("outputs/methodology/donors_harmonized_geo.rds"), prefix_levels)
keep_cols <- function(d) droplevels(as.data.frame(d[, c(INDEX_INPUTS, "age", "sex", "cluster")]))
pools <- lapply(build_quadrants(harmonized, quadrant_cells("methodology/v2/quadrant_cells.csv")), keep_cols)
pools$bd_allage <- keep_cols(build_age_pool(harmonized, "bangladesh", c(18L, Inf)))

cluster_tau <- function(score, t) {
  ok <- !is.na(score)
  m <- tapply(score[ok], t$cluster[ok], mean)
  stats::cor(m, tapply(t$V, t$cluster, mean)[names(m)], method = "spearman")
}

orient <- function(s_dev, dev) if (dev_tiers(s_dev, dev)$flip) -1 else 1

# One probe replicate: development diagnostics, the two-direction prediction of
# targeting, and targeting of the missing fuzzy average against dropping and
# the adaptive imputer on the same data.
check_one <- function(pop, gate, r, variant = "base") {
  spec <- modifyList(pilot_spec(population = pop, gate = gate, r2 = 0.6), variant_overrides(variant))
  set.seed(7000 * r + nchar(pop) + nchar(variant))
  sim <- generate_replicate(pools[[pop]], spec, 1000L, 20000L, signal_scale(pools[[pop]], spec))
  dev <- sim$development
  items <- advice_names(dev)
  diag <- coding_diagnostics(dev) |> mutate(pop, gate, variant, rep = r)

  mpX <- mca_dev(dev, INDEX_INPUTS)
  muX <- mean(mpX$coord_dev); sdX <- sqrt(mean((mpX$coord_dev - muX)^2))
  hX <- (mpX$coord_dev - muX) / sdX
  sgX <- orient(mpX$coord_dev, dev)
  out <- list()
  for (nm in names(CODINGS)) {
    cd <- CODINGS[[nm]]
    dc <- code_items(dev, cd$how)
    ex <- if (cd$excl) notasked_levels(dev) else character(0)
    mpA <- mca_dev(dc, items, excl = if (cd$excl) ex)
    muA <- mean(mpA$coord_dev); sdA <- sqrt(mean((mpA$coord_dev - muA)^2))
    hA <- (mpA$coord_dev - muA) / sdA
    mp <- mca_dev(dc, c(INDEX_INPUTS, items), excl = if (cd$excl) ex)
    vars <- c(INDEX_INPUTS, items)
    # Two-direction reduction: the generalized eigenproblem of the MCA
    # criterion restricted to span{hX, hA}.
    M <- matrix(c(b_form(dc, vars, hX, hX, ex), b_form(dc, vars, hX, hA, ex),
                  b_form(dc, vars, hA, hX, ex), b_form(dc, vars, hA, hA, ex)), 2)
    rXA <- mean(hX * hA)
    e <- eigen(solve(matrix(c(1, rXA, rXA, 1), 2), M))
    w <- Re(e$vectors[, which.max(Re(e$values))])
    s_pred <- w[1] * hX + w[2] * hA
    s1 <- unit(mp$coord_dev)
    sg_pred <- orient(s_pred, dev)
    sg_act <- orient(mp$coord_dev, dev)
    for (s in spec$settings) {
      t <- sim$targets[[s]]
      tc <- code_items(t, cd$how)
      xT <- (mca_coords(mpX, t) - muX) / sdX
      aT <- (mca_coords(mpA, tc) - muA) / sdA
      out[[length(out) + 1]] <- tibble(
        coding = nm, setting = s,
        span_r2 = summary(stats::lm(s1 ~ hX + hA))$r.squared,
        pred_axis_cor = abs(stats::cor(unit(s_pred), s1)),
        tau_actual = cluster_tau(sg_act * mca_coords(mp, tc), t),
        tau_pred = cluster_tau(sg_pred * (w[1] * xT + w[2] * aT), t),
        tau_filter = cluster_tau(t$eligible, t),
        var_cluster_G = stats::var(tapply(t$eligible, t$cluster, mean)))
    }
  }
  prep_fa <- fuzzy_average_prep(dev)
  imp <- prepare_family(dev, "impute", spec)
  for (s in spec$settings) {
    t <- sim$targets[[s]]
    fa <- prep_fa$project(t, "fuzzy_average")$score
    ia <- imp$project(t, "impute_adaptive")$score
    out[[length(out) + 1]] <- tibble(
      coding = "fuzzy_average", setting = s,
      fit_check = max(abs(abs(prep_fa$fit$coord_dev) - abs(prep_fa$fit$check))),
      tau_actual = cluster_tau(if (prep_fa$dt$flip) -fa else fa, t),
      tau_impute_adaptive = cluster_tau(if (imp$dt$flip) -ia else ia, t),
      tau_drop = cluster_tau(sgX * (mca_coords(mpX, t) - muX) / sdX, t))
  }
  list(diag = diag, target = bind_rows(out) |> mutate(pop, gate, variant, rep = r, k = length(items)))
}

n_rep <- as.integer(Sys.getenv("NREP", "12"))
runs <- bind_rows(
  expand.grid(pop = c("Q3", "bd_allage"), gate = GATES, variant = "base", rep = seq_len(n_rep),
              stringsAsFactors = FALSE),
  expand.grid(pop = "Q3", gate = "contact",
              variant = c("advice_strong", "rho_zero", "rho_strong", "items_two", "items_eight", "recog_flat"),
              rep = seq_len(n_rep), stringsAsFactors = FALSE),
  expand.grid(pop = "Q3", gate = "diagnosis", variant = c("q1_low", "q1_high"), rep = seq_len(n_rep),
              stringsAsFactors = FALSE))
res <- lapply(seq_len(nrow(runs)), \(i) check_one(runs$pop[i], runs$gate[i], runs$rep[i], runs$variant[i]))
diag <- bind_rows(lapply(res, `[[`, "diag"))
target <- bind_rows(lapply(res, `[[`, "target"))
saveRDS(list(diag = diag, target = target), "JCE/theory/theory_checks.rds")

num <- \(x) round(mean(x, na.rm = TRUE), 3)
show <- \(x) print(as.data.frame(x), row.names = FALSE)

cat("\n== P1. Filter loading: empirical against the chi-square formula ==\n")
show(diag |> filter(variant == "base", coding != "fuzzy_average") |> group_by(gate, pop, coding) |>
       summarise(p_dev = num(p_dev), phi = num(phi), phi_formula = num(phi_formula),
                 max_abs_diff = signif(max(abs(phi - phi_formula)), 2), .groups = "drop"))

cat("\n== P2. Capture ratio against the share of the axis taken by the block ==\n")
show(diag |> filter(coding != "fuzzy_average") |> group_by(variant, gate, pop, coding) |>
       summarise(kappa = num(kappa), block_share = num(block_share), filter_share = num(filter_share),
                 need_share = num(need_share), gap = num(eig1 - eig2), .groups = "drop"))

cat("\n== P3. Targeting predicted by the two-direction reduction ==\n")
tp <- target |> filter(coding != "fuzzy_average")
show(tp |> filter(variant == "base") |> group_by(gate, coding, setting) |>
       summarise(tau_actual = num(tau_actual), tau_pred = num(tau_pred), tau_filter = num(tau_filter),
                 var_cluster_G = signif(mean(var_cluster_G), 3), span_r2 = num(span_r2), .groups = "drop"))
kap <- diag |> select(variant, gate, pop, rep, coding, kappa)
tpk <- tp |> left_join(kap, by = c("variant", "gate", "pop", "rep", "coding"))
far <- abs(tpk$kappa - 1) > 0.25
cat("\nAll replicates: r =", round(stats::cor(tpk$tau_pred, tpk$tau_actual, use = "complete"), 3),
    " mean |error| =", round(mean(abs(tpk$tau_pred - tpk$tau_actual), na.rm = TRUE), 3),
    "\n|kappa - 1| > .25: r =", round(stats::cor(tpk$tau_pred[far], tpk$tau_actual[far], use = "complete"), 3),
    " mean |error| =", round(mean(abs(tpk$tau_pred - tpk$tau_actual)[far], na.rm = TRUE), 3),
    "\n|kappa - 1| <= .25: mean |error| =",
    round(mean(abs(tpk$tau_pred - tpk$tau_actual)[!far], na.rm = TRUE), 3), "\n")

cat("\n== P4. Instability near the threshold: replicate SD of targeting against |kappa - 1| ==\n")
show(tpk |> filter(setting == "c50") |> group_by(variant, gate, pop, coding) |>
       summarise(kappa = num(kappa), tau_mean = num(tau_actual), tau_sd = round(sd(tau_actual), 3),
                 .groups = "drop") |> arrange(abs(kappa - 1)))

cat("\n== R1. Missing fuzzy average against dropping and the adaptive imputer ==\n")
show(target |> filter(coding == "fuzzy_average") |> group_by(variant, gate, pop) |>
       summarise(fit_check = signif(max(fit_check), 2), fuzzy_average = num(tau_actual),
                 drop = num(tau_drop), impute_adaptive = num(tau_impute_adaptive),
                 fa_minus_drop = num(tau_actual - tau_drop),
                 fa_minus_impute = num(tau_actual - tau_impute_adaptive),
                 fa_sd = round(sd(tau_actual), 3), impute_sd = round(sd(tau_impute_adaptive), 3),
                 .groups = "drop"))
show(diag |> filter(coding == "fuzzy_average") |> group_by(variant, gate) |>
       summarise(filter_share = num(filter_share), .groups = "drop"))

cat("\n== Corollary: fuzzy-average block eigenvalue against the share asked ==\n")
show(diag |> filter(coding == "fuzzy_average") |> group_by(variant, gate, pop) |>
       summarise(p_dev = num(p_dev), lamA = num(lamA), max_ratio = round(max(lamA / p_dev), 3),
                 kappa = num(kappa), .groups = "drop"))
