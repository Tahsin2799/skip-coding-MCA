# Recompute a run's metrics from the saved row-level artifacts, without calling
# the helpers the runner used.
#
#   Rscript methodology/v2/verify_from_artifacts.R <run directory>
#
# The runner saves, for the first `--artifact-reps` replicates of a
# configuration, every arm's target score and the development cuts it is read
# against. This file recomputes coverage, the tier and budget metrics and the
# cluster metrics in plain base R and compares them with what the run stored. A
# check that reuses target_metrics() would only prove the function is
# deterministic.

suppressPackageStartupMessages({library(dplyr)})

root <- c(".", "..", "../..") |>
  Filter(f = \(p) file.exists(file.path(p, "methodology/v2/plasmode_generate.R"))) |>
  (\(x) x[1])()
source(file.path(root, "methodology/v2/run_manifest.R"))

recompute_one <- function(art, arm, setting) {
  a <- art$scores[[paste(arm, setting, sep = "_")]]
  if (is.null(a)) return(NULL)
  truth <- art$truth
  acc <- art$access[[setting]]
  n <- nrow(truth)

  s <- a$score
  if (a$flip) s <- -s
  brk <- a$breaks
  cl <- pmin(pmax(s, min(brk)), max(brk))
  tier <- rep(NA_integer_, n)
  ok <- !is.na(s)
  tier[ok] <- as.integer(cut(cl[ok], brk, include.lowest = TRUE, labels = FALSE))
  flagged <- !is.na(tier) & tier == a$n_tiers

  hn <- truth$V >= quantile(truth$V, 0.8)
  unrec <- truth$y == 1L & acc$recognized == 0L

  # Budget selection: highest scores first, ties broken on the generator's
  # uniform, unscoreable people never selected.
  slots <- floor(0.10 * n)
  cand <- which(ok)
  cand <- cand[order(-s[cand], truth$tiebreak[cand])]
  sel <- head(cand, min(slots, length(cand)))

  # Cluster priority: mean score of the scored residents; a cluster with none
  # ranks below every cluster that has one.
  cid <- as.integer(factor(truth$cluster))
  ncl <- max(cid)
  cs <- rep(NA_real_, ncl)
  for (j in seq_len(ncl)) {
    v <- s[cid == j & ok]
    if (length(v)) cs[j] <- mean(v)
  }
  cl_v <- tapply(truth$V, cid, mean)
  cl_tb <- tapply(truth$cluster_tiebreak, cid, function(x) x[1])
  rank_in <- ifelse(is.na(cs), min(cs, na.rm = TRUE) - 1, cs)
  cslots <- floor(0.10 * ncl)
  ccand <- which(!is.na(cs))
  ccand <- ccand[order(-cs[ccand], cl_tb[ccand])]
  csel <- head(ccand, min(cslots, length(ccand)))

  c(coverage = sum(ok) / n,
    capture_all = sum(flagged & hn) / sum(hn),
    tier_flag_rate = mean(flagged),
    unrec_capture_tier = if (sum(unrec)) sum(flagged & unrec) / sum(unrec) else NA_real_,
    b10_used = length(sel),
    b10_hn_capture = sum(hn[sel]) / sum(hn),
    b10_unrec_capture = if (sum(unrec)) sum(unrec[sel]) / sum(unrec) else NA_real_,
    cluster_rho_all = suppressWarnings(cor(rank_in, as.numeric(cl_v), method = "spearman")),
    cl10_misalloc = if (length(csel)) mean(cl_v[csel] < median(cl_v)) else NA_real_,
    cl10_hn_reached = sum(hn & cid %in% csel) / sum(hn))
}

if (sys.nframe() == 0L) {
  dir <- commandArgs(trailingOnly = TRUE)[1]
  arts <- list.files(file.path(dir, "artifacts"), pattern = "\\.rds$", full.names = TRUE)
  if (!length(arts)) stop("no artifacts in ", dir, call. = FALSE)
  mf <- manifest_files(dir)
  results <- bind_rows(lapply(mf$paths, readRDS))

  rows <- list()
  for (p in arts) {
    art <- readRDS(p)
    for (arm in unique(results$arm)) {
      for (setting in unique(results$setting)) {
        mine <- recompute_one(art, arm, setting)
        if (is.null(mine)) next
        theirs <- results |>
          filter(config_id == art$config_id, rep == art$rep, arm == !!arm, setting == !!setting)
        if (!nrow(theirs)) next
        for (k in names(mine)) {
          rows[[length(rows) + 1L]] <- data.frame(
            config_id = art$config_id, rep = art$rep, arm = arm, setting = setting,
            metric = k, mine = mine[[k]], stored = as.numeric(theirs[[k]]))
        }
      }
    }
  }
  cmp <- bind_rows(rows) |>
    mutate(diff = abs(mine - stored),
           agree = is.na(mine) & is.na(stored) | (!is.na(diff) & diff < 1e-9))

  cat("\nindependent recomputation from artifacts\n")
  cat("  files:      ", length(arts), "\n")
  cat("  comparisons:", nrow(cmp), "\n")
  cat("  disagreeing:", sum(!cmp$agree), "\n")
  print(as.data.frame(cmp |> group_by(metric) |>
    summarise(n = n(), worst_abs_diff = max(diff, na.rm = TRUE),
              disagree = sum(!agree), .groups = "drop")), digits = 3, row.names = FALSE)
  if (any(!cmp$agree)) {
    print(as.data.frame(head(cmp[!cmp$agree, ], 20)), digits = 10, row.names = FALSE)
    quit(status = 1L)
  }
}
