# Read a v03 run directory, verify it, and apply the claim rule.
#
#   Rscript methodology/v2/summarize_convergence.R <run directory> [--margins]
#
# Nothing here pools archived v2 runs: their schema and their estimands differ,
# and the reader refuses a directory without a v03 manifest.

options(width = 150)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

root <- c(".", "..", "../..") |>
  Filter(f = \(p) file.exists(file.path(p, "methodology/v2/plasmode_generate.R"))) |>
  (\(x) x[1])()
source(file.path(root, "methodology/v2/run_manifest.R"))

# Direction is metadata, not something the reader remembers: +1 means larger is
# better. `needs` names a validity flag that must be FALSE for the value to
# count, which keeps a separated tier regression from discarding the targeting
# metrics computed on the same rows.
METRIC_META <- tibble::tribble(
  ~metric,                 ~direction, ~family,      ~primary, ~needs,
  "cluster_rho_all",        1,        "targeting",   TRUE,     NA,
  "cl10_misalloc",         -1,        "targeting",   TRUE,     NA,
  "cl10_precision",         1,        "targeting",   FALSE,    NA,
  "cl10_hn_reached",        1,        "targeting",   FALSE,    NA,
  "misallocation_all",     -1,        "targeting",   FALSE,    NA,
  "b10_hn_capture",         1,        "budget",      TRUE,     NA,
  "b10_hn_shortfall_per1000", -1,     "budget",      FALSE,    NA,
  "b05_hn_capture",         1,        "budget",      FALSE,    NA,
  "b20_hn_capture",         1,        "budget",      FALSE,    NA,
  "b10_unrec_capture",      1,        "detection",   FALSE,    NA,
  "capture_all",            1,        "detection",   FALSE,    NA,
  "unrec_capture_tier",     1,        "detection",   FALSE,    NA,
  "coverage",               1,        "coverage",    FALSE,    NA,
  "auc_model_scored",       1,        "apparent",    FALSE,    "separation_target",
  "auc_index_scored",       1,        "apparent",    FALSE,    "separation_target",
  "log_or_top",             1,        "apparent",    FALSE,    "separation_dev",
  "log_or_top_target",      1,        "apparent",    FALSE,    "separation_target",
  "beta_within",            1,        "apparent",    FALSE,    "separation_dev",
  "beta_between",           1,        "apparent",    FALSE,    "separation_dev",
  "cluster_rho_scored",     1,        "apparent",    FALSE,    NA,
  "capture_scored",         1,        "apparent",    FALSE,    NA,
  "cal_slope",              1,        "calibration", FALSE,    "separation_target",
  "cal_intercept",          1,        "calibration", FALSE,    "separation_target"
)

# Equivalence margins in the units of the metric. They are claims about what is
# substantively negligible for the intended decision, so they are set before any
# amended result is read and are used only when `--margins` is passed. See
# REVISION_DECISIONS.md D3.5: unconfirmed margins mean no equivalence claim and
# a reported range of differences instead.
MARGINS <- c(cluster_rho_all = 0.02, cl10_misalloc = 0.02, b10_hn_capture = 0.01,
             capture_all = 0.01, b10_unrec_capture = 0.01)

ARM_ORDER <- c("oracle", "raw_items", "drop_items", "simple_score", "passive_mca",
               "impute_frozen", "impute_adaptive", "complete_case", "zero", "category")

read_run <- function(dir) {
  mf <- manifest_files(dir)
  if (!length(mf$paths)) stop("no manifest-listed result files in ", dir, call. = FALSE)
  d <- bind_rows(lapply(mf$paths, function(p) {
    x <- readRDS(p)
    prov <- attr(x, "provenance")
    if (!identical(prov$signature, mf$manifest$signature)) {
      stop(basename(p), " does not carry this directory's manifest signature", call. = FALSE)
    }
    x
  }))
  d$arm <- factor(d$arm, levels = intersect(ARM_ORDER, unique(d$arm)))
  d$log_or_top <- log(d$or_top)
  d$log_or_top_target <- log(d$or_top_target)
  d$cell <- paste(d$config_id, d$setting, sep = "|")
  attr(d, "manifest") <- mf$manifest
  attr(d, "missing_configs") <- mf$missing
  d
}

# Completeness, duplicates, leftovers and realized rates, before any estimate is
# read. A verification that runs after the numbers have been looked at is not a
# verification.
verify_run <- function(dir, d = read_run(dir), n_rep = NULL) {
  m <- attr(d, "manifest")
  ok <- d |> filter(!failed)
  reps <- d |> distinct(config_id, rep) |> count(config_id, name = "reps")
  n_rep <- n_rep %||% max(reps$reps)
  n_arms <- length(m$arms)
  list(
    manifest_signature = m$signature,
    configurations = c(in_manifest = nrow(m$grid), on_disk = length(unique(d$config_id)),
                       missing = length(attr(d, "missing_configs"))),
    short_configurations = reps$config_id[reps$reps != n_rep],
    duplicate_rows = sum(duplicated(d[, c("config_id", "rep", "arm", "setting")])),
    incomplete_replicates = d |> count(config_id, rep, name = "rows") |>
      filter(rows != n_arms * 3L) |> nrow(),
    failed_rows = sum(d$failed),
    failures = d |> filter(failed) |> count(config_id, arm, setting, stage, name = "n"),
    separated = c(dev = sum(ok$separation_dev, na.rm = TRUE),
                  target = sum(ok$separation_target, na.rm = TRUE)),
    level_fixups = sum(ok$level_fixups, na.rm = TRUE),
    rates = d |> distinct(config_id, rep, setting, contact_rate, outcome_rate,
                          recognized_rate, eligible_rate) |>
      group_by(setting) |>
      summarise(across(c(contact_rate, outcome_rate, recognized_rate, eligible_rate), mean),
                .groups = "drop"),
    stale_locks = list.files(dir, pattern = "^run\\.lock$", include.dirs = TRUE),
    partial_files = list.files(dir, pattern = "^pending-"),
    sources_match = identical(
      m$source_hashes,
      setNames(vapply(file.path(root, "methodology/v2", names(m$source_hashes)),
                      \(p) digest::digest(file = p, algo = "sha256"), character(1)),
               names(m$source_hashes))),
    negative_controls = access_invariance(ok)
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# The advice-free arms read no gated item and the outcome is imposed before
# access, so nothing they produce may move across the access grid. Any spread is
# a contaminated contrast, not a small effect.
access_invariance <- function(d) {
  cols <- intersect(c("cluster_rho_all", "capture_all", "b10_hn_capture", "auc_index_scored"),
                    names(d))
  d |>
    filter(arm %in% ADVICE_FREE_ARMS_V03) |>
    group_by(config_id, rep, arm) |>
    summarise(across(all_of(cols), \(x) diff(range(x, na.rm = TRUE))), .groups = "drop") |>
    group_by(arm) |>
    summarise(across(all_of(cols), \(x) max(x, na.rm = TRUE)), .groups = "drop")
}

ADVICE_FREE_ARMS_V03 <- c("drop_items", "simple_score", "oracle")

# ---- The claim rule -----------------------------------------------------

# Paired within replicate and computed inside one scenario cell. The arms see the
# same people, so most Monte Carlo error cancels in the difference; pooling cells
# first would hide a reversal between them.
cellwise_contrasts <- function(d, metric, pairs = NULL, cells = NULL) {
  meta <- METRIC_META[METRIC_META$metric == metric, ]
  if (!nrow(meta)) stop("no metadata for metric ", metric, call. = FALSE)
  x <- d |> filter(!failed)
  if (!is.na(meta$needs)) {
    flag <- x[[meta$needs]]
    x <- x[is.na(flag) | !flag, ]
  }
  if (!is.null(cells)) x <- x[x$cell %in% cells, ]
  arms <- levels(droplevels(x$arm))
  if (is.null(pairs)) pairs <- utils::combn(arms, 2L, simplify = FALSE)
  w <- x |>
    select(cell, config_id, setting, rep, arm, value = all_of(metric)) |>
    pivot_wider(names_from = arm, values_from = value)
  out <- list()
  for (p in pairs) {
    if (!all(p %in% names(w))) next
    for (g in unique(w$cell)) {
      wg <- w[w$cell == g, ]
      dif <- meta$direction * (wg[[p[1]]] - wg[[p[2]]])
      keep <- is.finite(dif)
      n <- sum(keep)
      dif <- dif[keep]
      se <- if (n > 1L) stats::sd(dif) / sqrt(n) else NA_real_
      out[[length(out) + 1L]] <- tibble(
        metric = metric, better = p[1], worse = p[2], cell = g,
        config_id = wg$config_id[1], setting = wg$setting[1],
        diff = if (n) mean(dif) else NA_real_, mcse = se, n = n, n_dropped = nrow(wg) - n,
        # An exact tie has a zero difference and a zero MCSE; that is a tie, not
        # a division by zero. A non-zero difference with no spread is decisive.
        t = if (!n || is.na(se)) NA_real_ else if (se == 0) {
          if (mean(dif) == 0) 0 else sign(mean(dif)) * Inf
        } else mean(dif) / se)
    }
  }
  bind_rows(out)
}

# One verdict per pair over the cells it was measured in. Superiority has to
# hold in every required cell; a sign change is scenario dependence and is
# reported as such; equivalence needs a prespecified margin and the whole
# three-MCSE interval inside it everywhere.
classify_pairs <- function(ct, margin = NULL) {
  ct |>
    group_by(metric, better, worse) |>
    summarise(
      cells = n(), usable = sum(!is.na(t)),
      pos = sum(!is.na(t) & t > 3), neg = sum(!is.na(t) & t < -3),
      min_diff = min(diff, na.rm = TRUE), max_diff = max(diff, na.rm = TRUE),
      worst_mcse = suppressWarnings(max(mcse, na.rm = TRUE)),
      inside = if (is.null(margin)) NA_integer_ else
        sum(!is.na(mcse) & abs(diff) + 3 * mcse < margin),
      .groups = "drop") |>
    mutate(claim = case_when(
      usable < cells ~ "incomplete",
      pos == cells ~ "robust_superiority",
      neg == cells ~ "robust_inferiority",
      pos > 0 & neg > 0 ~ "scenario_dependent",
      !is.na(inside) & inside == cells ~ "equivalent",
      TRUE ~ "unresolved"))
}

# Replicates needed for a target MCSE on each paired difference, from the
# observed per-cell spread: ceiling((sd / target)^2).
required_reps <- function(ct, target_mcse) {
  ct |>
    mutate(sd_pair = mcse * sqrt(n),
           need = ceiling((sd_pair / target_mcse)^2)) |>
    group_by(metric, better, worse) |>
    summarise(worst_cell = cell[which.max(need)], need = max(need), .groups = "drop") |>
    arrange(desc(need))
}

# Spearman between what a validation study reports and what the index is for,
# oriented so that positive means the two agree. Computed inside a scenario cell
# first: a correlation pooled over cells is a correlation over a mixture of
# designs, which is a different quantity.
family_correlation <- function(d, apparent = NULL, targeting = NULL) {
  apparent <- apparent %||% METRIC_META$metric[METRIC_META$family == "apparent"]
  targeting <- targeting %||% METRIC_META$metric[METRIC_META$family %in% c("targeting", "budget")]
  dirs <- setNames(METRIC_META$direction, METRIC_META$metric)
  x <- d |> filter(!failed)
  g <- expand.grid(a = intersect(apparent, names(x)), t = intersect(targeting, names(x)),
                   stringsAsFactors = FALSE)
  per_cell <- lapply(split(x, x$cell), function(xc) {
    g$rho <- mapply(\(a, t) {
      if (stats::sd(xc[[a]], na.rm = TRUE) == 0 || stats::sd(xc[[t]], na.rm = TRUE) == 0) return(NA_real_)
      dirs[[a]] * dirs[[t]] * stats::cor(xc[[a]], xc[[t]], method = "spearman",
                                         use = "pairwise.complete.obs")
    }, g$a, g$t)
    g$cell <- xc$cell[1]
    g
  })
  bind_rows(per_cell) |>
    group_by(a, t) |>
    summarise(cells = sum(!is.na(rho)), median_rho = stats::median(rho, na.rm = TRUE),
              min_rho = min(rho, na.rm = TRUE), max_rho = max(rho, na.rm = TRUE),
              .groups = "drop")
}

if (sys.nframe() == 0L) {
  argv <- commandArgs(trailingOnly = TRUE)
  dir <- if (length(argv)) argv[1] else "outputs/methodology/production_v03"
  use_margins <- "--margins" %in% argv
  d <- read_run(dir)
  v <- verify_run(dir, d)
  core <- d |> filter(block == "core")
  core_cells <- unique(core$cell)

  cat("\nrun directory: ", dir, "\n", sep = "")
  cat("manifest:      ", substr(v$manifest_signature, 1, 12), " (", attr(d, "manifest")$stream, ")\n", sep = "")
  cat("configurations:", v$configurations[["on_disk"]], "of", v$configurations[["in_manifest"]],
      "| short:", length(v$short_configurations), "| duplicate rows:", v$duplicate_rows,
      "| incomplete replicates:", v$incomplete_replicates, "\n")
  cat("rows failed:   ", v$failed_rows, " | separated dev/target: ", v$separated[["dev"]],
      "/", v$separated[["target"]], " | level substitutions: ", v$level_fixups, "\n", sep = "")
  cat("sources match manifest:", v$sources_match, "| partial files:", length(v$partial_files), "\n")
  if (nrow(v$failures)) {
    cat("\nfailures, which are never random\n")
    print(as.data.frame(v$failures), row.names = FALSE)
  }

  cat("\nrealized rates by setting\n")
  print(as.data.frame(v$rates), digits = 4, row.names = FALSE)

  cat("\nadvice-free arms across the access grid (must be exactly zero)\n")
  print(as.data.frame(v$negative_controls), digits = 3, row.names = FALSE)

  cat("\nshare of the target each arm can score, by gate and setting\n")
  print(as.data.frame(d |> filter(!failed) |> group_by(arm, gate, setting) |>
    summarise(coverage = mean(coverage), .groups = "drop") |>
    pivot_wider(names_from = c(gate, setting), values_from = coverage)),
    digits = 3, row.names = FALSE)

  for (m in METRIC_META$metric[METRIC_META$primary]) {
    cat("\n", m, ": paired contrasts over the core grid\n", sep = "")
    ct <- cellwise_contrasts(core, m, cells = core_cells)
    cl <- classify_pairs(ct, margin = if (use_margins) MARGINS[[m]] else NULL)
    print(as.data.frame(cl |> arrange(claim, desc(min_diff))), digits = 3, row.names = FALSE)
  }

  cat("\nwhole target against the rows each arm can score, c20 and c80\n")
  print(as.data.frame(d |> filter(!failed, setting %in% c("c20", "c80")) |>
    group_by(gate, arm, setting) |>
    summarise(coverage = mean(coverage), rho_scored = mean(cluster_rho_scored, na.rm = TRUE),
              rho_all = mean(cluster_rho_all, na.rm = TRUE),
              capture_scored = mean(capture_scored, na.rm = TRUE),
              capture_all = mean(capture_all), .groups = "drop") |>
    arrange(gate, setting, desc(rho_all))), digits = 3, row.names = FALSE)

  cat("\napparent performance against targeting, Spearman within cell\n")
  print(as.data.frame(family_correlation(core)), digits = 3, row.names = FALSE)
  cat("\n  Read the range, not one number: these are correlations between what a\n")
  cat("  validation study reports and what the index is for, one per scenario cell.\n")
}
