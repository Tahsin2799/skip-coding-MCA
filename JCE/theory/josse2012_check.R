# Does Theorem 1 explain a published result? Josse, Chavent, Liquet and Husson
# (J Classif 2012;29:91-116, section 5.2, Table 3) found that missing passive
# modified margin, which equals subset MCA, collapsed under their "non-random"
# missingness pattern (rows 7 and 8: modified RV .43 and .71, against .77 and
# .91 for missing fuzzy average) without saying why. In that pattern the first
# 60 of 100 individuals miss variables 1 to 3 together and the last 60 miss
# variables 9 and 10 together: two blocks sharing one skip each, the setting of
# Theorem 1.
#
# This rebuilds their design and records, beside their outcome, how much of the
# method's first two axes the two missingness indicators explain.
#
#   Rscript JCE/theory/josse2012_check.R      # 500 replications per cell

suppressPackageStartupMessages({library(FactoMineR); library(MASS)})

n_rep <- as.integer(Sys.getenv("NREP", "500"))

# Their protocol: 100 x 10 multivariate normal, two-block diagonal correlation
# (variables 1-6 and 7-10), within-block correlation .4 (low) or .8 (strong),
# each variable cut into three equal-count categories.
simulate <- function(rho) {
  S <- matrix(0, 10, 10)
  S[1:6, 1:6] <- rho
  S[7:10, 7:10] <- rho
  diag(S) <- 1
  x <- mvrnorm(100, rep(0, 10), S)
  d <- as.data.frame(lapply(seq_len(10), \(j) factor(paste0("v", j, "_",
    cut(x[, j], quantile(x[, j], c(0, 1/3, 2/3, 1)), include.lowest = TRUE, labels = FALSE)))))
  names(d) <- paste0("v", 1:10)
  d
}

# Their non-random MCAR pattern at 30 percent.
make_missing <- function(d) {
  d[1:60, 1:3] <- NA
  d[41:100, 9:10] <- NA
  d
}

# Smilde et al. (2009) modified RV between two configurations.
rv_mod <- function(X, Y) {
  a <- tcrossprod(X); diag(a) <- 0
  b <- tcrossprod(Y); diag(b) <- 0
  sum(a * b) / sqrt(sum(a * a) * sum(b * b))
}

single_coded <- function(d) {
  as.data.frame(lapply(names(d), \(v) {
    x <- as.character(d[[v]])
    x[is.na(x)] <- paste0(v, "_NA")
    factor(x)
  }), col.names = names(d))
}

one_rep <- function(rho) {
  full <- simulate(rho)
  miss <- make_missing(full)
  truth <- MCA(full, ncp = 2, graph = FALSE)$ind$coord[, 1:2]
  sc <- single_coded(miss)
  na_levels <- paste0(names(miss), "_NA")
  na_levels <- intersect(na_levels, unlist(lapply(sc, levels)))
  coords <- list(
    single = MCA(sc, ncp = 2, graph = FALSE)$ind$coord[, 1:2],
    passive_modified_margin = MCA(sc, ncp = 2, graph = FALSE, excl = na_levels)$ind$coord[, 1:2],
    fuzzy_average = as.matrix(CA(tab.disjonctif.prop(miss), ncp = 2, graph = FALSE)$row$coord[, 1:2]))
  G <- cbind(g1 = as.numeric(!is.na(miss$v1)), g2 = as.numeric(!is.na(miss$v9)))
  do.call(rbind, lapply(names(coords), \(m) data.frame(
    link = ifelse(rho == 0.4, "low", "strong"), method = m,
    rv_mod = rv_mod(truth, coords[[m]]),
    # Mean share of the first two axes explained by the two missingness indicators
    filter_share = mean(vapply(1:2, \(a) summary(lm(coords[[m]][, a] ~ G))$r.squared, numeric(1))))))
}

set.seed(20121)
res <- do.call(rbind, lapply(rep(c(0.4, 0.8), each = n_rep), one_rep))
out <- aggregate(cbind(rv_mod, filter_share) ~ link + method, data = res, FUN = median)
out[, 3:4] <- round(out[, 3:4], 2)
print(out, row.names = FALSE)
cat("\nPublished medians (Table 3, rows 7-8): low  Passive .43, Average .77, NA .30;",
    "strong Passive .71, Average .91, NA .45\n")
# Theorem 1 for each item of the blocks: own category 1, subset 1 - p with
# p = .4 asked, fuzzy average 0.
saveRDS(res, "JCE/theory/josse2012_check.rds")
