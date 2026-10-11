# Quadrant populations for the v2 plasmode simulation.
#
# The structural axis is not the donor. It is a 2x2 of population deprivation
# level against the coherence of the deprivation structure, because those two are
# what decide whether a one-dimensional categorical index can work at all, and
# because they are close to orthogonal across real subpopulations (r = -0.12).
#
# Populations are built from DHS region by urban/rural strata, pooled only within
# a quadrant. Pooling strata that DIFFER in deprivation level inflates the first
# MCA axis: across the twelve largest Indian states the state-level adjusted
# inertia exceeds both its urban and its rural part, twelve times out of twelve.
# So membership is frozen in quadrant_cells.csv rather than recomputed.

suppressPackageStartupMessages({library(dplyr); library(FactoMineR); library(cli)})

INPUTS <- c("water", "sanitation", "fuel", "bmi", "education", "work")

deprivation_score <- function(d) {
  rowMeans(cbind((as.integer(d$water) - 1) / 3,
                 (as.integer(d$sanitation) - 1) / 3,
                 (as.integer(d$fuel) - 1) / 3,
                 (6 - as.integer(d$education)) / 5))
}

# Bergsma (2013) bias-corrected Cramer's V. The uncorrected statistic is biased
# upward in small samples, and the bias is large at the sizes these cells run to:
# one fixed Indian cell read .096 at n = 700 and .063 at its full n = 16,281.
# Since the quadrant split sits near .104, that made the coherence axis partly a
# cell-size axis. Measured on the 52 frozen cells, cor(n, coherence) was -.467
# uncorrected and -.136 corrected, and 18 of 63 profiled cells change which side
# of the split they fall on. Do not revert this to the raw statistic.
cramers_v <- function(a, b) {
  t <- table(a, b)
  if (nrow(t) < 2 || ncol(t) < 2) return(NA_real_)
  n <- sum(t); r <- nrow(t); k <- ncol(t)
  if (n < 2) return(NA_real_)
  phi2 <- as.numeric(suppressWarnings(chisq.test(t, correct = FALSE)$statistic)) / n
  phi2 <- max(0, phi2 - (k - 1) * (r - 1) / (n - 1))
  rt <- r - (r - 1)^2 / (n - 1)
  kt <- k - (k - 1)^2 / (n - 1)
  sqrt(phi2 / min(kt - 1, rt - 1))
}

coherence <- function(d) {
  d <- droplevels(d)
  mean(vapply(combn(INPUTS, 2, simplify = FALSE),
              function(p) cramers_v(d[[p[1]]], d[[p[2]]]), numeric(1)), na.rm = TRUE)
}

# Greenacre-adjusted share of the first MCA axis.
adjusted_inertia_dim1 <- function(d) {
  d <- droplevels(d[, INPUTS])
  Q <- ncol(d); J <- sum(vapply(d, nlevels, integer(1)))
  lam <- MCA(d, ncp = 1, graph = FALSE)$eig[, 1]
  adj <- (Q / (Q - 1))^2 * (lam[lam > 1 / Q] - 1 / Q)^2
  if (!length(adj)) return(NA_real_)
  unname(100 * adj[1] / ((Q / (Q - 1)) * (sum(lam^2) - (J - Q) / Q^2)))
}

quadrant_cells <- function(path = "methodology/v2/quadrant_cells.csv") {
  read.csv(path, stringsAsFactors = FALSE)
}

# The age range is an argument, not a constant. It defaults to 18-49 because
# that is the widest window in which all six inputs exist for the donors the
# quadrants are built from: `bmi` and `work` come from the individual
# questionnaires, whose universe is women 15-49 and men 15-54, so both are
# exactly zero above 54 and Indian records aged 50-54 are all male. See
# AGE_WINDOW_2026-09-18.md. Widening it changes `age_slope`, which is
# window-specific.
cell_pool <- function(harmonized, age_range = c(18L, 49L)) {
  bind_rows(lapply(names(harmonized), function(k)
    harmonized[[k]] |>
      filter(if_all(all_of(INPUTS), ~ !is.na(.x)),
             age >= age_range[1], age <= age_range[2]) |>
      mutate(cell = paste0(k, " | ", region, " | ",
                           ifelse(urban == 1, "urban", "rural")))))
}

build_quadrants <- function(harmonized, cells = quadrant_cells(),
                            age_range = c(18L, 49L)) {
  pool <- cell_pool(harmonized, age_range)
  split(pool, cells$quadrant[match(pool$cell, cells$cell)])
}

# The age arm. Bangladesh 2017-18 is the only donor that spans adult life with
# both sexes, because occupation comes from `sb308` inside the PR file so the
# ever-married universe never binds, and the biomarker subsample covers 18 and
# over rather than the individual questionnaire's ceiling. 12,921 complete
# records across 672 clusters, 18 to 95, sex balance holding in every band.
build_age_pool <- function(harmonized, donor = "bangladesh",
                           age_range = c(18L, Inf)) {
  d <- cell_pool(harmonized[donor], age_range)
  if (!nrow(d)) {
    cli::cli_abort("{.val {donor}} has no records complete on all six inputs in {age_range[1]}-{age_range[2]}.")
  }
  d
}
