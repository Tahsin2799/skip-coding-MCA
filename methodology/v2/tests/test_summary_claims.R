# Comparison, direction and equivalence contracts.
# Run: Rscript methodology/v2/tests/test_summary_claims.R
#
# These are the rules that decide what the paper may say, so they are tested on
# fixtures whose answer is known in advance rather than on a run.

suppressPackageStartupMessages({library(dplyr); library(tibble)})
root <- c(".", "..", "../..") |>
  Filter(f = \(p) file.exists(file.path(p, "methodology/v2/plasmode_generate.R"))) |>
  (\(x) x[1])()
source(file.path(root, "methodology/v2/summarize_convergence.R"))

# A frame shaped like a run: two arms, R replicates, and whatever per-cell
# differences the caller asks for.
# `shared` is the replicate effect both arms see, which pairing removes;
# `pair_sd` is the part that does not cancel and so drives the paired MCSE.
make_run <- function(effects, R = 200, pair_sd = 0.01, shared = 0.05, extra = list()) {
  set.seed(7)
  rows <- lapply(names(effects), function(nm) {
    base <- rnorm(R, 0.5, shared)
    eff <- effects[[nm]]
    bind_rows(
      tibble(cell = nm, rep = seq_len(R), arm = "A", value = base + eff + rnorm(R, 0, pair_sd)),
      tibble(cell = nm, rep = seq_len(R), arm = "B", value = base + rnorm(R, 0, pair_sd)))
  })
  d <- bind_rows(rows) |>
    mutate(config_id = sub("\\|.*", "", cell), setting = sub(".*\\|", "", cell),
           failed = FALSE, block = "core", cluster_rho_all = value,
           separation_dev = FALSE, separation_target = FALSE)
  for (nm in names(extra)) d[[nm]] <- extra[[nm]]
  d$arm <- factor(d$arm, levels = c("A", "B"))
  d
}

# ---- A method that wins in one cell and loses in another -----------------
# Its pooled mean is positive, which is exactly the trap: a hardcoded ranking
# checked at one setting would claim the ordering.
mixed <- make_run(list(`c1|c50` = 0.05, `c2|c50` = -0.02, `c3|c50` = 0.04))
ct <- cellwise_contrasts(mixed, "cluster_rho_all")
cl <- classify_pairs(ct)
stopifnot("a sign change is scenario dependence, not an ordering" =
            cl$claim == "scenario_dependent")
stopifnot("the pooled mean would have been positive" = mean(ct$diff) > 0)
stopifnot("both directions are counted" = cl$pos == 2 && cl$neg == 1)
stopifnot("the range of differences is reported" =
            cl$min_diff < 0 && cl$max_diff > 0)

# ---- A method that wins everywhere ---------------------------------------
won <- classify_pairs(cellwise_contrasts(
  make_run(list(`c1|c50` = 0.05, `c2|c50` = 0.04, `c3|c50` = 0.06)), "cluster_rho_all"))
stopifnot("winning in every cell by more than three MCSE is robust superiority" =
            won$claim == "robust_superiority" && won$pos == 3)

# ---- Near-zero effects with wide uncertainty -----------------------------
noisy <- make_run(list(`c1|c50` = 0.001, `c2|c50` = -0.001), R = 30, pair_sd = 0.4)
nc <- classify_pairs(cellwise_contrasts(noisy, "cluster_rho_all"), margin = 0.01)
stopifnot("imprecision is unresolved, never equivalence" = nc$claim == "unresolved")

# ---- A genuinely negligible difference, precisely measured ---------------
tiny <- make_run(list(`c1|c50` = 0.0005, `c2|c50` = -0.0005), R = 4000, pair_sd = 0.02)
eq <- classify_pairs(cellwise_contrasts(tiny, "cluster_rho_all"), margin = 0.01)
stopifnot("equivalence needs the whole interval inside the margin, in every cell" =
            eq$claim == "equivalent" && eq$inside == 2)
stopifnot("and without a margin the same data are only unresolved" =
            classify_pairs(cellwise_contrasts(tiny, "cluster_rho_all"))$claim == "unresolved")

# ---- An exact tie is a tie, not a division by zero -----------------------
zero <- make_run(list(`c1|c50` = 0), pair_sd = 0)
zero$cluster_rho_all <- 0.5
zt <- cellwise_contrasts(zero, "cluster_rho_all")
stopifnot("zero difference with zero spread reads as t = 0" =
            zt$diff == 0 && zt$mcse == 0 && zt$t == 0)
stopifnot("and it is not claimed as an ordering" =
            classify_pairs(zt)$claim == "unresolved")

# ---- The sign convention comes from the metric metadata ------------------
# Lower misallocation is better, so an arm with the smaller value must be the
# one the contrast favours. The v2 summary checked the old column name and so
# silently oriented this family the wrong way.
mis <- make_run(list(`c1|c50` = -0.05, `c2|c50` = -0.04))
mis$cl10_misalloc <- mis$value
stopifnot("a lower-is-better metric is oriented by its metadata" =
            classify_pairs(cellwise_contrasts(mis, "cl10_misalloc"))$claim == "robust_superiority")
stopifnot("every metric the summary reports carries a direction" =
            all(METRIC_META$direction %in% c(-1, 1)))
stopifnot("misallocation is one of them, under its current name" =
            "misallocation_all" %in% METRIC_META$metric &&
            METRIC_META$direction[METRIC_META$metric == "misallocation_all"] == -1)

# ---- Validity flags are metric-specific ----------------------------------
# A separated tier regression invalidates the odds ratio computed from it. It
# does not invalidate the cluster correlation computed on the same rows.
sep <- make_run(list(`c1|c50` = 0.05))
sep$log_or_top <- sep$value
sep$separation_dev <- sep$rep <= 50
stopifnot("separated rows are dropped from the odds ratio" =
            max(cellwise_contrasts(sep, "log_or_top")$n) == 150)
stopifnot("but not from the targeting metric" =
            max(cellwise_contrasts(sep, "cluster_rho_all")$n) == 200)

# ---- Incomplete pairs are reported, not silently dropped -----------------
gap <- make_run(list(`c1|c50` = 0.05))
gap$cluster_rho_all[gap$arm == "A" & gap$rep <= 20] <- NA_real_
gc <- cellwise_contrasts(gap, "cluster_rho_all")
stopifnot("only complete pairs are used" = gc$n == 180,
          "and the dropped ones are counted" = gc$n_dropped == 20)

# ---- Oriented association between the two families -----------------------
# Rising AUC with falling misallocation is agreement, so the oriented
# correlation must be positive.
set.seed(11)
fam <- tibble(cell = "c1|c50", config_id = "c1", setting = "c50", rep = 1:200,
              arm = factor("A"), failed = FALSE,
              auc_model_scored = seq(0.6, 0.8, length.out = 200) + rnorm(200, 0, 0.005))
fam$misallocation_all <- 0.5 - (fam$auc_model_scored - 0.6) + rnorm(200, 0, 0.005)
fam$separation_target <- FALSE
fc <- family_correlation(fam, apparent = "auc_model_scored", targeting = "misallocation_all")
stopifnot("agreement between the families reads positive" = fc$median_rho > 0.9)

# ---- Replicate planning --------------------------------------------------
need <- required_reps(cellwise_contrasts(
  make_run(list(`c1|c50` = 0.01), R = 200, pair_sd = 0.05), "cluster_rho_all"), 0.002)
stopifnot("the required replicate count scales as (sd / target)^2" =
            need$need > 500 && need$need < 5000)

cat("test_summary_claims.R: all contracts hold\n")
