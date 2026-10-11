# Filter-capture quantities for an MCA index whose filtered items share one
# eligibility indicator G. Supplement A6 of the paper gives the definitions and proofs.
#
#   source("JCE/theory/filter_capture.R")
#
# Needs index_arms.R (code_items, mca_dev, mca_coords, dev_tiers, dev_metrics,
# observed_rows, advice_names, notasked_levels, INDEX_INPUTS) and FactoMineR.

# ---- Homogeneity algebra ----------------------------------------------------

# Centered score with unit population variance.
unit <- function(x) {
  x <- x - mean(x)
  x / sqrt(mean(x^2))
}

# h'P_q g / n for one categorical variable: sum over categories of
# n_c m_c(h) m_c(g) / n, for centered h and g. Subset analysis drops the excluded
# category from the sum but keeps the overall centering.
cross_eta <- function(f, h, g, drop = character(0)) {
  f <- as.character(f)
  nc <- tapply(h, f, length)
  keep <- !names(nc) %in% drop
  sum((nc * tapply(h, f, mean) * tapply(g, f, mean))[keep]) / length(h)
}

# The MCA criterion h'Bg / n, the mean of cross_eta over variables. With h = g of
# unit variance it is the mean squared correlation ratio, which the first axis
# maximizes.
b_form <- function(d, vars, h, g, excl = character(0)) {
  mean(vapply(vars, \(v) cross_eta(d[[v]], h, g, excl), numeric(1)))
}

# Filter loading of an item whose skipped respondents all receive the answer
# profile v (a distribution over the analysed categories): p(1 - p) times the
# chi-square distance of v from the eligible profile, weighted by overall column
# frequency. Own category, subset MCA and skip-as-"no" are special cases.
filter_loading <- function(p, f_eligible, v, analysed = rep(TRUE, length(v))) {
  f <- p * f_eligible + (1 - p) * v
  keep <- analysed & f > 0
  p * (1 - p) * sum((v[keep] - f_eligible[keep])^2 / f[keep])
}

# ---- Missing fuzzy average -------------------------------------------------

# Missing fuzzy average (van der Heijden and Escofier 2003, section 12.3.6;
# Josse et al. 2012, section 3.1): a missing cell of the indicator matrix gets
# the category proportions of the non-missing objects. For a skip those are the
# answers of the people asked, which by Theorem 1 leaves no filter direction.
# MCA of a fuzzy indicator matrix is correspondence analysis of that matrix;
# any row, in development or target, is placed with the transition formula on
# the development column standard coordinates.
asked_profile <- function(d, items) {
  setNames(lapply(items, \(v) {
    x <- d[[v]][!is.na(d[[v]])]
    if (!length(x)) return(c(0.5, 0.5))
    c(mean(x == 0L), mean(x == 1L))
  }), items)
}

fuzzy_items <- function(d, items, prof) {
  do.call(cbind, lapply(items, \(v) {
    x <- d[[v]]
    m <- cbind(as.numeric(x %in% 0L), as.numeric(x %in% 1L))
    sk <- is.na(x)
    m[sk, ] <- matrix(prof[[v]], sum(sk), 2, byrow = TRUE)
    colnames(m) <- paste0(v, c("_0", "_1"))
    m
  }))
}

fuzzy_z <- function(d, lv, items, prof) {
  zx <- do.call(cbind, lapply(INDEX_INPUTS, \(v) {
    m <- outer(as.character(d[[v]]), lv[[v]], `==`) * 1
    colnames(m) <- lv[[v]]
    m
  }))
  cbind(zx, fuzzy_items(d, items, prof))
}

# First eigenvalue of the filtered block analysed alone under missing fuzzy
# average. The corollary to Theorem 1 bounds it by the share asked.
fuzzy_block_eig <- function(dev, items) {
  za <- fuzzy_items(dev, items, asked_profile(dev, items))
  za <- za[, colSums(za) > 0, drop = FALSE]
  FactoMineR::CA(as.data.frame(za), ncp = 1, graph = FALSE)$eig[1, 1]
}

fuzzy_fit <- function(dev, items) {
  lv <- lapply(setNames(INDEX_INPUTS, INDEX_INPUTS), \(v) levels(droplevels(dev[[v]])))
  z <- fuzzy_z(dev, lv, items, asked_profile(dev, items))
  z <- z[, colSums(z) > 0, drop = FALSE]
  fit <- FactoMineR::CA(as.data.frame(z), ncp = 1, graph = FALSE)
  gamma <- as.matrix(fit$col$coord)[, 1] / sqrt(fit$eig[1, 1])
  q <- length(INDEX_INPUTS) + length(items)
  list(lv = lv, cols = colnames(z), gamma = gamma, q = q, eig = fit$eig[1, 1],
       coord_dev = as.numeric((z / q) %*% gamma), check = as.matrix(fit$row$coord)[, 1])
}

# Target profiles come from the target's own eligible respondents, so the
# filtered columns stay orthogonal to eligibility wherever the index is used.
# No model is fitted. Rows with an input category development never saw are
# left unscored, as for every other MCA arm.
fuzzy_score <- function(ff, d, items, prof = asked_profile(d, items)) {
  ok <- Reduce(`&`, lapply(INDEX_INPUTS, \(v) as.character(d[[v]]) %in% ff$lv[[v]]))
  z <- fuzzy_z(d, ff$lv, items, prof)[, ff$cols, drop = FALSE]
  s <- as.numeric((z / ff$q) %*% ff$gamma)
  s[!ok] <- NA_real_
  s
}

# A prepared family in the shape score_prepared() expects, so the new coding is
# scored by the production metric code unchanged.
fuzzy_average_prep <- function(dev) {
  items <- advice_names(dev)
  ff <- fuzzy_fit(dev, items)
  out <- list(dev = dev, s_dev = ff$coord_dev, n_fix = 0L, family = "fuzzy_average",
              fit = ff, project = function(tgt, arm) {
                s <- rep(NA_real_, nrow(tgt))
                s[seq_len(nrow(tgt))] <- fuzzy_score(ff, tgt, items)
                list(score = s, n_fix = 0L)
              })
  out$dt <- dev_tiers(out$s_dev, out$dev)
  out$dm <- dev_metrics(out$dev, out$dt)
  out
}

# ---- Development diagnostics ---------------------------------------------------

CODINGS <- list(category = list(how = "category", excl = FALSE),
                subset   = list(how = "category", excl = TRUE),
                zero     = list(how = "zero",     excl = FALSE))

# Everything here uses development data only, so it is available before an
# index is deployed and needs no truth.
coding_diagnostics <- function(dev) {
  items <- advice_names(dev)
  k <- length(items)
  P <- length(INDEX_INPUTS)
  G <- as.numeric(observed_rows(dev))
  p <- mean(G)
  sG <- unit(G)
  mpX <- mca_dev(dev, INDEX_INPUTS)
  hX <- unit(mpX$coord_dev)
  lamX <- b_form(dev, INDEX_INPUTS, hX, hX)
  rows <- lapply(names(CODINGS), function(nm) {
    cd <- CODINGS[[nm]]
    dc <- code_items(dev, cd$how)
    ex <- if (cd$excl) notasked_levels(dev) else character(0)
    phi <- vapply(items, \(v) cross_eta(dc[[v]], sG, sG, ex), numeric(1))
    pi_k <- vapply(items, \(v) mean(dev[[v]][G == 1] == 1L), numeric(1))
    phi_formula <- vapply(pi_k, \(pk) switch(nm,
      category = filter_loading(p, c(1 - pk, pk, 0), c(0, 0, 1)),
      subset   = filter_loading(p, c(1 - pk, pk, 0), c(0, 0, 1), c(TRUE, TRUE, FALSE)),
      zero     = filter_loading(p, c(1 - pk, pk), c(1, 0))), numeric(1))
    mpA <- mca_dev(dc, items, excl = if (cd$excl) ex)
    hA <- unit(mpA$coord_dev)
    lamA <- b_form(dc, items, hA, hA, ex)
    mp <- mca_dev(dc, c(INDEX_INPUTS, items), excl = if (cd$excl) ex)
    s1 <- unit(mp$coord_dev)
    tibble::tibble(coding = nm, k = k, p_dev = p, phi = mean(phi), phi_formula = mean(phi_formula),
                   lamX = lamX, lamA = lamA, kappa = (k * lamA) / (P * lamX),
                   eig1 = mp$fit$eig[1, 1], eig2 = mp$fit$eig[2, 1],
                   filter_share = stats::cor(s1, G)^2, block_share = stats::cor(s1, hA)^2,
                   need_share = stats::cor(s1, hX)^2)
  })
  ff <- fuzzy_fit(dev, items)
  lamA_fa <- fuzzy_block_eig(dev, items)
  dplyr::bind_rows(rows, tibble::tibble(coding = "fuzzy_average", k = k, p_dev = p, phi = 0,
                                        phi_formula = 0, lamX = lamX, lamA = lamA_fa,
                                        kappa = (k * lamA_fa) / (P * lamX), eig1 = ff$eig,
                                        filter_share = stats::cor(ff$coord_dev, G)^2,
                                        need_share = stats::cor(ff$coord_dev, hX)^2))
}
