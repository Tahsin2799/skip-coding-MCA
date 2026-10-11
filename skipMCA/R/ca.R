# Correspondence analysis of an indicator matrix, crisp or fuzzy, with an
# optional subset of analysed columns that keep their full-data masses
# (Greenacre and Pardo 2006). Every row sums to q, the number of items.

ca_core <- function(z, q, keep = rep(TRUE, ncol(z)), ncp = 1) {
  n <- nrow(z)
  cm <- colSums(z) / (n * q)
  used <- keep & cm > 0
  # Standardized residuals: S = D_r^(-1/2) (P - r c') D_c^(-1/2) with r = 1/n.
  a <- z[, used, drop = FALSE] / q
  s <- sweep(sweep(a, 2, cm[used]), 2, sqrt(cm[used]), "/") / sqrt(n)
  ncp <- min(ncp, ncol(s), nrow(s))
  sv <- svd(s, nu = 0, nv = ncp)
  gamma <- matrix(0, ncol(z), ncp, dimnames = list(colnames(z), paste0("dim", seq_len(ncp))))
  gamma[used, ] <- sv$v[, seq_len(ncp), drop = FALSE] / sqrt(cm[used])
  # The sign of an axis is arbitrary; fix it so that the column with the
  # largest absolute standard coordinate is positive.
  flip <- apply(gamma, 2, function(g) sign(g[which.max(abs(g))]))
  gamma <- sweep(gamma, 2, flip, "*")
  list(eig = sv$d^2, gamma = gamma, cm = cm, used = used, q = q)
}

# Row principal coordinates by the transition formula,
# f_i = sum over analysed columns of (z_ij / q - c_j) gamma_j. With every column
# analysed the centring term vanishes; under subset analysis it does not.
ca_project <- function(core, z) {
  u <- core$used
  a <- z[, u, drop = FALSE] / core$q
  sweep(a, 2, core$cm[u]) %*% core$gamma[u, , drop = FALSE]
}

# Squared correlation ratio of one coded item with a score h, centred and
# scaled to unit population variance, over the analysed columns.
eta2 <- function(z, h, keep = rep(TRUE, ncol(z))) {
  h <- h - mean(h)
  h <- h / sqrt(mean(h^2))
  nj <- colSums(z)
  ok <- keep & nj > 0
  sum((colSums(z[, ok, drop = FALSE] * h))^2 / nj[ok]) / length(h)
}
