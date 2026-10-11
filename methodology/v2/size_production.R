# How many replicates the production run needs, read off a pilot directory.
#
#   Rscript methodology/v2/size_production.R <pilot directory> [N ...]
#
# `run_checks.R` reports the replicates needed to reach a target Monte Carlo
# error on the worst of all arm pairs. That is the wrong quantity for sizing:
# the paper claims only the prespecified contrasts, and the frozen rule decides
# a contrast by the sign of its per-cell t, not by an absolute MCSE. So this
# holds each cell's observed mean difference and spread fixed, shrinks the MCSE
# as sqrt(n_pilot / N), re-applies `classify_pairs`, and counts what is decided.
#
# A projection, not a measurement: it assumes the per-cell spreads carry to the
# production stream, which they should, since only the replicate count changes.

suppressPackageStartupMessages({library(dplyr); library(tidyr)})

root <- c(".", "..", "../..") |>
  Filter(f = \(p) file.exists(file.path(p, "methodology/v2/plasmode_generate.R"))) |>
  (\(x) x[1])()
source(file.path(root, "methodology/v2/summarize_convergence.R"))
source(file.path(root, "methodology/v2/claim_ledger.R"))

# Wall seconds per replicate on six workers, measured over the whole v03 grid.
SEC_PER_REP <- 1.67

project_claims <- function(core, N) {
  cells <- unique(core$cell)
  pairs <- prespecified_pairs(levels(droplevels(core$arm)))
  bind_rows(lapply(METRIC_META$metric[METRIC_META$primary], function(m) {
    ct <- cellwise_contrasts(core, m, pairs = pairs, cells = cells) |>
      mutate(mcse = mcse * sqrt(n / N),
             # An exact tie stays a tie under any N; a difference with no spread
             # stays decisive. Same convention as cellwise_contrasts().
             t = ifelse(is.na(diff) | is.na(mcse), NA_real_,
                 ifelse(mcse == 0, ifelse(diff == 0, 0, sign(diff) * Inf), diff / mcse)))
    classify_pairs(ct, margin = MARGINS[[m]]) |> mutate(N = N)
  }))
}

if (sys.nframe() == 0L) {
  argv <- commandArgs(trailingOnly = TRUE)
  dir <- argv[1]
  Ns <- if (length(argv) > 1L) as.integer(argv[-1]) else c(100L, 250L, 500L, 1000L, 2000L, 4000L)
  d <- read_run(dir)
  core <- d |> filter(block == "core")
  n_config <- nrow(attr(d, "manifest")$grid)
  res <- bind_rows(lapply(Ns, \(N) project_claims(core, N)))

  cat("\n== projected verdicts on the prespecified contrasts\n")
  tab <- res |> count(N, claim) |> pivot_wider(names_from = claim, values_from = n, values_fill = 0L)
  tab$decided <- rowSums(tab[, setdiff(names(tab), c("N", "unresolved"))])
  tab$hours <- round(n_config * tab$N * SEC_PER_REP / 3600, 1)
  print(as.data.frame(tab), row.names = FALSE)

  cat("\n== still unresolved at the largest count tried\n")
  print(as.data.frame(res |> filter(N == max(Ns), claim == "unresolved") |>
    select(metric, better, worse, pos, neg, cells, min_diff, max_diff)),
    digits = 3, row.names = FALSE)

  cat("\n== verdicts that change as N grows\n")
  moved <- res |> group_by(metric, better, worse) |>
    filter(n_distinct(claim) > 1L) |> ungroup()
  print(as.data.frame(moved |> select(metric, better, worse, N, claim) |>
    pivot_wider(names_from = N, values_from = claim)), row.names = FALSE)
}
