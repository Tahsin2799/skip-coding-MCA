# Subset (specific) MCA projection contracts.
# Run: Rscript methodology/v2/tests/test_passive_mca.R
#
# FactoMineR's predict.MCA ignores `excl`. The evidence and the alternative are
# in COMPARATOR_FEASIBILITY.md; this file holds the line.

suppressPackageStartupMessages({library(dplyr); library(digest); library(FactoMineR)})
root <- c(".", "..", "../..") |>
  Filter(f = \(p) file.exists(file.path(p, "methodology/v2/plasmode_generate.R"))) |>
  (\(x) x[1])()
for (f in c("harmonize_donors", "plasmode_generate", "index_arms",
            "quadrant_populations", "pilot_v2_config")) {
  source(file.path(root, sprintf("methodology/v2/%s.R", f)))
}

# ---- The package-level discrepancy, on the package's own data ------------
data(poison, package = "FactoMineR")
X <- poison[, 3:8]
fit <- MCA(X, ncp = 1, excl = c(1, 3), graph = FALSE)
fitted <- unname(fit$ind$coord[, 1])
hand <- {
  Z <- tab.disjonctif(X)
  fj <- colSums(Z) / (nrow(Z) * ncol(X))
  v <- fit$svd$V[, 1]
  keep <- rep(1, length(v))
  keep[c(1, 3)] <- 0
  as.numeric((Z / ncol(X)) %*% (v * keep)) - sum(fj * v * keep)
}
predicted <- suppressWarnings(predict(fit, newdata = X)$coord[, 1])
stopifnot("the transition formula reproduces the fitted coordinates" =
            max(abs(fitted - hand)) < 1e-10)
stopifnot("predict.MCA does not, so it must never be used with excl" =
            max(abs(fitted - predicted)) > 0.1)
sup <- MCA(rbind(X, X), ncp = 1, excl = c(1, 3),
           ind.sup = (nrow(X) + 1):(2 * nrow(X)), graph = FALSE)
stopifnot("the documented supplementary route agrees with the formula" =
            max(abs(unname(sup$ind.sup$coord[, 1]) - hand)) < 1e-10)
stopifnot("and supplementary rows leave the active fit untouched" =
            isTRUE(all.equal(fit$eig, sup$eig)))

# ---- The same checks on the study's own frames ---------------------------
pools <- readRDS(file.path(root, "outputs/methodology/donors_harmonized_geo.rds")) |>
  lapply(prefix_levels) |>
  build_quadrants(quadrant_cells(file.path(root, "methodology/v2/quadrant_cells.csv"))) |>
  lapply(\(d) droplevels(as.data.frame(d[, c(INDEX_INPUTS, "age", "sex", "cluster")])))
spec <- pilot_spec(population = "Q3")
set.seed(digest2int("test_passive_1"))
sim <- generate_replicate(pools[["Q3"]], spec, 1000L, 20000L)

prep <- prepare_family(sim$development, "passive_mca", spec)
dev_cat <- code_items(sim$development, "category")
vars <- c(INDEX_INPUTS, advice_names(sim$development))
mp <- mca_dev(dev_cat, vars, excl = notasked_levels(sim$development))
stopifnot("the not-asked categories are the excluded ones" =
            setequal(mp$excl, notasked_levels(sim$development)))

dev_f <- droplevels(dev_cat[, vars, drop = FALSE])
stopifnot("projecting the training rows returns their fitted coordinates" =
            max(abs(passive_coord(mp$proj, dev_f) - mp$coord_dev)) < 1e-9)

tgt <- sim$targets$c20
coord <- mca_coords(mp, code_items(tgt, "category"))
# A target row carrying a level the development draw never saw is left unscored
# and counted against coverage; here that is a rare fuel rung.
ok <- !is.na(coord)
stopifnot("most rows are scoreable and the rest are unscored, not wrong" =
            mean(ok) > 0.95 && any(!ok))

# Permutation, duplication and a changed case mix must not move anyone's score:
# the projection reads one row at a time against frozen training quantities.
gap <- function(a, b) max(abs(a - b), na.rm = TRUE)
perm <- sample.int(nrow(tgt))
stopifnot("projection is invariant to target row order" =
            gap(mca_coords(mp, code_items(tgt[perm, ], "category")), coord[perm]) < 1e-12)
dup <- rbind(tgt[1:100, ], tgt[1:100, ])
cd <- mca_coords(mp, code_items(dup, "category"))
stopifnot("duplicated rows get identical coordinates" = gap(cd[1:100], cd[101:200]) < 1e-12)
contacted_only <- tgt[tgt$contacted == 1L, ]
stopifnot("a different target case mix does not move the survivors" =
            gap(mca_coords(mp, code_items(contacted_only, "category")),
                coord[tgt$contacted == 1L]) < 1e-12)
stopifnot("and does not change who can be scored" =
            identical(is.na(mca_coords(mp, code_items(contacted_only, "category"))),
                      is.na(coord[tgt$contacted == 1L])))

# The excluded categories carry no weight, so a skipped item moves nobody, while
# an observed response still does. That is what "passive" means here.
stopifnot("excluded categories have exactly zero weight" =
            all(mp$proj$w[mp$proj$cols %in% mp$excl] == 0),
          "and the categories that remain do not" =
            any(mp$proj$w[!mp$proj$cols %in% mp$excl] != 0))
skipped <- tgt[1:50, ]
for (v in advice_names(skipped)) skipped[[v]] <- NA_integer_
answered <- skipped
for (v in advice_names(answered)) answered[[v]] <- 1L
stopifnot("an observed response moves the coordinate" =
            any(abs(mca_coords(mp, code_items(answered, "category")) -
                      mca_coords(mp, code_items(skipped, "category"))) > 1e-6))

# ---- Levels the training data never saw ----------------------------------
odd <- sim$targets$c20
before <- score_prepared(prep, "passive_mca", odd, spec)
lv <- levels(odd$water)
odd$water <- factor(as.character(odd$water), levels = c(lv, "water_unseen"))
odd$water[which(!is.na(coord))[1:5]] <- "water_unseen"
after <- score_prepared(prep, "passive_mca", odd, spec)
stopifnot("unseen levels cost coverage rather than aborting" =
            after$n_scored == before$n_scored - 5L && after$coverage < before$coverage)

# ---- A near-degenerate input ---------------------------------------------
deg <- sim$development
deg$water <- factor(rep(levels(droplevels(deg$water))[1], nrow(deg)))
p2 <- try(prepare_family(deg, "passive_mca", spec), silent = TRUE)
stopifnot("a constant input does not abort the fit" = !inherits(p2, "try-error"))
stopifnot("and the arm still scores the target" =
            is.finite(score_prepared(p2, "passive_mca", sim$targets$c50, spec)$cluster_rho_all))

# ---- Cost -----------------------------------------------------------------
t0 <- Sys.time()
invisible(mca_coords(mp, code_items(sim$targets$c50, "category")))
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
stopifnot("projecting 20,000 target rows is cheap" = elapsed < 10)

cat(sprintf("test_passive_mca.R: all contracts hold (20,000-row projection %.2fs)\n", elapsed))
