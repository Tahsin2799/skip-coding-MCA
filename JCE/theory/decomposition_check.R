# Theorem 2 of supplement A6: the joint analysis as a rank-one update.
#
#   Q B = P B_X + c k B_E + k theta u u'
#
# with (c, theta) = (1, 1) under missing single and (1, 1 - p) under subset MCA,
# B_E the analysis of the filtered items among those asked (zero for the
# skipped) and u the unit eligibility score. Checked on probe replicates built as in
# theory_checks.R (unchanged production generator, never production estimates)
# with indicator matrices built here, against FactoMineR through
# coding_diagnostics(). Also the multiple factor analysis weighting under
# missing single. Writes only JCE/theory/decomposition_check.rds.
#
#   Rscript JCE/theory/decomposition_check.R            # 50 replicates per cell
#   NREP=2 Rscript JCE/theory/decomposition_check.R     # smoke

suppressPackageStartupMessages({library(dplyr); library(FactoMineR)})
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

# Columns (z_j - f_j 1) / sqrt(w n_j), so that their cross-product is the matrix
# B of Lemma A6.0 for a block of w items. Empty columns are dropped.
tcols <- function(Z, w) {
  Z <- Z[, colSums(Z) > 0, drop = FALSE]
  nj <- colSums(Z)
  sweep(sweep(Z, 2, nj / nrow(Z)), 2, sqrt(w * nj), "/")
}
top <- function(Tm) { s <- svd(Tm, nv = 0); list(l = s$d[1]^2, h = s$u[, 1]) }

pieces <- function(dev) {
  items <- advice_names(dev)
  G <- as.numeric(observed_rows(dev))
  E <- G == 1
  n <- nrow(dev); p <- mean(G); k <- length(items); P <- length(INDEX_INPUTS)
  ZX <- do.call(cbind, lapply(INDEX_INPUTS, \(v) {
    x <- droplevels(dev[[v]])
    outer(as.character(x), levels(x), `==`) * 1
  }))
  ans <- lapply(items, \(v) cbind(as.numeric(dev[[v]] %in% 0L), as.numeric(dev[[v]] %in% 1L)))
  TE <- matrix(0, n, 0)
  te <- tcols(do.call(cbind, lapply(ans, \(a) a[E, , drop = FALSE])), k)
  TE <- matrix(0, n, ncol(te)); TE[E, ] <- te
  list(n = n, p = p, k = k, P = P, Q = P + k, G = G, ZX = ZX, TE = TE,
       u = (G - p) / sqrt(n * p * (1 - p)),
       Z = list(category = do.call(cbind, lapply(ans, \(a) cbind(a, 1 - G))),
                subset = do.call(cbind, ans)))
}

# Largest root of the secular equation of B_0 + rho u u', and the squared
# correlation of its eigenvector with u (Golub 1973).
secular <- function(T0, u, rho) {
  s <- svd(T0, nv = 0)
  keep <- s$d^2 > 1e-12
  mu <- c(s$d[keep]^2, 0)
  c2 <- as.numeric(crossprod(s$u[, keep, drop = FALSE], u))^2
  c2 <- c(c2, max(0, 1 - sum(c2)))
  f <- function(l) rho * sum(c2 / (l - mu)) - 1
  l1 <- stats::uniroot(f, c(max(mu) + 1e-12, max(mu) + rho + 1e-9), tol = 1e-15)$root
  list(l1 = l1, r2 = 1 / (rho^2 * sum(c2 / (l1 - mu)^2)), mu1 = max(mu))
}

check_dev <- function(dev) {
  x <- pieces(dev)
  fm <- coding_diagnostics(dev)
  TX <- tcols(x$ZX, x$P)
  lamX <- top(TX)$l
  lamE <- if (ncol(x$TE)) top(x$TE)$l else 0
  T0 <- cbind(sqrt(x$P / x$Q) * TX, sqrt(x$k / x$Q) * x$TE)
  B0 <- tcrossprod(T0)
  hX <- top(TX)$h
  rows <- lapply(c("category", "subset"), function(nm) {
    theta <- if (nm == "category") 1 else 1 - x$p
    rho <- x$k * theta / x$Q
    Tj <- tcols(cbind(x$ZX, x$Z[[nm]]), x$Q)
    tj <- top(Tj)
    sc <- secular(T0, x$u, rho)
    f <- fm[fm$coding == nm, ]
    tibble::tibble(
      coding = nm, p = x$p, k = x$k,
      decomposition_gap = max(abs(tcrossprod(Tj) - B0 - rho * tcrossprod(x$u))),
      eig1 = tj$l, eig1_secular = sc$l1, eig1_factominer = f$eig1,
      r2_elig = stats::cor(tj$h, x$G)^2, r2_elig_secular = sc$r2, r2_elig_factominer = f$filter_share,
      P_lamX = x$P * lamX, k_lamE = x$k * lamE, k_theta = x$k * theta, Q_mu1 = x$Q * sc$mu1,
      kappa = f$kappa)
  })
  # Multiple factor analysis weighting, missing single: each block divided by
  # its own first eigenvalue (Escofier and Pages 1994).
  TA <- tcols(x$Z$category, x$k)
  lamA <- top(TA)$l
  hm <- top(cbind(TX / sqrt(lamX), TA / sqrt(lamA)))
  mfa <- tibble::tibble(lamA_single = lamA, mfa_eig1 = hm$l, mfa_r2_elig = stats::cor(hm$h, x$G)^2,
                        mfa_r2_inputs = stats::cor(hm$h, hX)^2)
  list(rows = dplyr::bind_rows(rows), mfa = mfa)
}

# FactoMineR's own MFA on one development sample, as a check of the weighting.
mfa_factominer <- function(dev) {
  items <- advice_names(dev)
  dc <- code_items(dev, "category")
  d <- droplevels(as.data.frame(lapply(dc[, c(INDEX_INPUTS, items)], factor)))
  FactoMineR::MFA(d, group = c(length(INDEX_INPUTS), length(items)), type = c("n", "n"),
                  ncp = 2, graph = FALSE)$eig[1, 1]
}

n_rep <- as.integer(Sys.getenv("NREP", "50"))
runs <- bind_rows(
  expand.grid(pop = c("Q3", "bd_allage"), gate = GATES, variant = "base", rep = seq_len(n_rep),
              stringsAsFactors = FALSE),
  expand.grid(pop = "Q3", gate = "contact", variant = c("items_two", "items_eight"),
              rep = seq_len(n_rep), stringsAsFactors = FALSE))
res <- lapply(seq_len(nrow(runs)), function(i) {
  r <- runs[i, ]
  spec <- modifyList(pilot_spec(population = r$pop, gate = r$gate, r2 = 0.6), variant_overrides(r$variant))
  set.seed(7000 * r$rep + nchar(r$pop) + nchar(r$variant))
  dev <- generate_replicate(pools[[r$pop]], spec, 1000L, 20000L, signal_scale(pools[[r$pop]], spec))$development
  out <- check_dev(dev)
  if (i == 1) out$mfa$mfa_eig1_factominer <- mfa_factominer(dev)
  list(rows = cbind(r, out$rows, row.names = NULL), mfa = cbind(r, out$mfa, row.names = NULL))
})
rows <- bind_rows(lapply(res, `[[`, "rows"))
mfa <- bind_rows(lapply(res, `[[`, "mfa"))
saveRDS(list(rows = rows, mfa = mfa), "JCE/theory/decomposition_check.rds")

num <- \(x, k = 3) round(mean(x, na.rm = TRUE), k)
show <- \(x) print(as.data.frame(x), row.names = FALSE)

cat("\n== Theorem 2 on", nrow(rows) / 2, "development samples, two codings ==\n")
cat("largest |Q B - (P B_X + k B_E + k theta u u')| / Q:", signif(max(rows$decomposition_gap), 2), "\n")
cat("largest |eig1 - secular root|:", signif(max(abs(rows$eig1 - rows$eig1_secular)), 2),
    "  against FactoMineR:", signif(max(abs(rows$eig1 - rows$eig1_factominer)), 2), "\n")
cat("largest |R2 with eligibility - secular value|:", signif(max(abs(rows$r2_elig - rows$r2_elig_secular)), 2),
    "  against FactoMineR:", signif(max(abs(rows$r2_elig - rows$r2_elig_factominer)), 2), "\n")

cat("\n== The three quantities competing for the first axis ==\n")
show(rows |> group_by(variant, gate, pop, coding) |>
       summarise(samples = n(), filter_first = num(k_theta > Q_mu1), kappa_above_1 = num(kappa > 1),
                 elig_over_half = num(r2_elig > .5), r2_elig_sd = round(sd(r2_elig), 3),
                 p = num(p), P_lamX = num(P_lamX), k_lamE = num(k_lamE), k_theta = num(k_theta),
                 Q_mu1 = num(Q_mu1), kappa = num(kappa), r2_elig = num(r2_elig), .groups = "drop"))
agree <- \(a, b) round(mean(a == b), 3)
cat("\nEligibility explains more than half the first axis in", num(rows$r2_elig > .5), "of sample-codings.\n")
cat("Agreement with that: k theta > Q mu1(B0)", agree(rows$k_theta > rows$Q_mu1, rows$r2_elig > .5),
    "; kappa > 1", agree(rows$kappa > 1, rows$r2_elig > .5), "\n")
sc <- rows[rows$coding == "subset" & rows$gate == "contact" & rows$variant == "base", ]
cat("Subset MCA, contact filter, four items (", nrow(sc), "samples): k theta > Q mu1(B0)",
    agree(sc$k_theta > sc$Q_mu1, sc$r2_elig > .5), "; kappa > 1", agree(sc$kappa > 1, sc$r2_elig > .5), "\n")

cat("\n== Multiple factor analysis weighting under missing single ==\n")
cat("block eigenvalue under missing single:", signif(min(mfa$lamA_single), 8), "to", signif(max(mfa$lamA_single), 8), "\n")
cat("first eigenvalue, here against FactoMineR::MFA (first sample):", signif(mfa$mfa_eig1[1], 8), "against",
    signif(mfa$mfa_eig1_factominer[1], 8), "\n")
show(mfa |> group_by(variant, gate, pop) |>
       summarise(r2_elig = num(mfa_r2_elig), r2_elig_min = round(min(mfa_r2_elig), 3),
                 r2_elig_max = round(max(mfa_r2_elig), 3), r2_inputs = num(mfa_r2_inputs), .groups = "drop"))
cat("All samples: mean", num(mfa$mfa_r2_elig), " range", round(min(mfa$mfa_r2_elig), 3), "to",
    round(max(mfa$mfa_r2_elig), 3), "\n")
