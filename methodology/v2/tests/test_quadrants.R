# Population construction contracts. Run: Rscript methodology/v2/tests/test_quadrants.R

suppressPackageStartupMessages({library(dplyr)})
root <- c(".", "..", "../..") |>
  Filter(f = \(p) file.exists(file.path(p, "methodology/v2/plasmode_generate.R"))) |>
  (\(x) x[1])()
for (f in c("harmonize_donors", "quadrant_populations")) {
  source(file.path(root, sprintf("methodology/v2/%s.R", f)))
}

h <- readRDS(file.path(root, "outputs/methodology/donors_harmonized_geo.rds"))
cells <- quadrant_cells(file.path(root, "methodology/v2/quadrant_cells.csv"))
pools <- build_quadrants(lapply(h, prefix_levels), cells)

# ---- Pools are what the design says they are ----------------------------
stopifnot("four quadrants" = length(pools) == 4L,
          "named Q1 to Q4" = identical(sort(names(pools)), paste0("Q", 1:4)))
for (q in names(pools)) {
  stopifnot("every pool clears the 20,000 target size" = nrow(pools[[q]]) >= 20000L,
            "every input is complete" =
              !any(vapply(INPUTS, \(v) anyNA(pools[[q]][[v]]), logical(1))))
}
stopifnot("the membership file has the columns the builder writes" =
            all(c("cell", "quadrant", "n", "deprivation", "coherence") %in% names(cells)))

# ---- Stored profiles must match what the functions compute ---------------
# The file used to be written ad hoc, so its stored coherence came from the
# uncorrected statistic while `coherence()` was later expected to reproduce it.
# Nothing checked, and 18 of 63 cells were on the wrong side of the split.
# Only build_quadrant_cells.R may write this file.
pool_all <- cell_pool(lapply(h, prefix_levels))
parts <- split(pool_all, pool_all$cell)
check <- sample(cells$cell, min(8L, nrow(cells)))
for (cn in check) {
  d <- parts[[cn]]
  stopifnot("the stored cell size matches the pool" =
              nrow(d) == cells$n[cells$cell == cn])
  stopifnot("the stored coherence matches coherence()" =
              abs(coherence(d) - cells$coherence[cells$cell == cn]) < 5e-4)
  stopifnot("the stored deprivation matches deprivation_score()" =
              abs(mean(deprivation_score(d)) - cells$deprivation[cells$cell == cn]) < 5e-4)
}

# ---- The two axes must be properties of the population, not of its size --
# `coherence` is a mean Cramer's V, which is biased upward in small samples.
# Uncorrected, the same population read .096 at n = 700 and .063 at n = 16,281,
# and the quadrant split sits at ~.104, so the axis was partly a size axis:
# 18 of 24 cells under n = 2,000 landed in a high-coherence quadrant.
set.seed(21)
big <- pools[["Q3"]]
big <- big[big$cell == names(which.max(table(big$cell))), ]
sizes <- c(700L, 2000L, 10000L, nrow(big))
coh <- vapply(sizes, \(m) mean(replicate(3,
  coherence(big[sample.int(nrow(big), min(m, nrow(big))), ]))), numeric(1))
dep <- vapply(sizes, \(m) mean(replicate(3,
  mean(deprivation_score(big[sample.int(nrow(big), min(m, nrow(big))), ])))), numeric(1))

stopifnot("deprivation does not depend on sample size" =
            max(abs(dep - dep[length(dep)])) < 0.01)
stopifnot("coherence does not depend on sample size" =
            max(abs(coh - coh[length(coh)])) < 0.01)

# ---- Quadrant boxes are disjoint on the two axes ------------------------
b <- cells |> group_by(quadrant) |>
  summarise(dep_lo = min(deprivation), dep_hi = max(deprivation),
            coh_lo = min(coherence), coh_hi = max(coherence), .groups = "drop")
stopifnot("the deprivation split separates {Q1,Q2} from {Q3,Q4}" =
            max(b$dep_hi[b$quadrant %in% c("Q1", "Q2")]) <
              min(b$dep_lo[b$quadrant %in% c("Q3", "Q4")]))
stopifnot("the coherence split separates {Q1,Q3} from {Q2,Q4}" =
            max(b$coh_hi[b$quadrant %in% c("Q1", "Q3")]) <
              min(b$coh_lo[b$quadrant %in% c("Q2", "Q4")]))

cat("test_quadrants.R: all contracts hold\n")
