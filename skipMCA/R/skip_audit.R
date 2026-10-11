#' Squared correlation ratio of a coded item with eligibility (Theorem 1)
#'
#' The squared correlation ratio between eligibility and a filtered item (in
#' homogeneity analysis, its discrimination measure on the eligibility score),
#' given the share asked, the answer distribution of those asked, and the mean
#' profile the coding gives the skipped:
#' \deqn{\varphi = p(1-p)\sum_{j \in S}\frac{(\bar v_j - f_{j|G})^2}{f_j},
#'   \qquad f_j = p f_{j|G} + (1-p)\bar v_j.}
#'
#' @param p Share of respondents asked.
#' @param asked Answer distribution of those asked over the item's columns;
#'   zero for a column only the skipped use, such as a skip category.
#' @param skipped Mean profile given to the skipped over the same columns.
#' @param analysed Logical, which columns the analysis keeps.
#' @return The squared correlation ratio, between 0 and 1.
#' @references Michailidis G, de Leeuw J (1998). The Gifi system of descriptive
#'   multivariate analysis. *Statistical Science* 13, 307-336.
#' @examples
#' # A binary item answered "yes" by 30 percent of the 40 percent asked.
#' p <- 0.4; pi <- 0.3
#' eligibility_eta2(p, c(1 - pi, pi, 0), c(0, 0, 1))                  # missing single: 1
#' eligibility_eta2(p, c(1 - pi, pi, 0), c(0, 0, 1), c(TRUE, TRUE, FALSE)) # subset: 1 - p
#' eligibility_eta2(p, c(1 - pi, pi), c(1, 0))                        # "no"
#' eligibility_eta2(p, c(1 - pi, pi), c(1 - pi, pi))                  # fuzzy average: 0
#' @export
eligibility_eta2 <- function(p, asked, skipped, analysed = rep(TRUE, length(asked))) {
  if (length(asked) != length(skipped) || length(analysed) != length(asked)) {
    stop("`asked`, `skipped` and `analysed` must have the same length.", call. = FALSE)
  }
  f <- p * asked + (1 - p) * skipped
  keep <- analysed & f > 0
  p * (1 - p) * sum((skipped[keep] - asked[keep])^2 / f[keep])
}

#' Audit questionnaire skips before building an MCA index
#'
#' Computes, from development data alone and for each coding of the skip, the
#' quantities that say whether a block of filtered items will hand the index
#' to the filter:
#' \describe{
#'   \item{squared correlation ratio with eligibility}{how much of the filter
#'     each filtered item carries, empirically and by [eligibility_eta2()];
#'     they agree exactly.}
#'   \item{inertia ratio}{\eqn{\kappa = k\lambda_A/(P\lambda_X)}: each block's
#'     maximum axial inertia (the first eigenvalue of its separate analysis, as
#'     in multiple factor analysis) times its number of items, filtered against
#'     substantive. The filtered block dominates the first axis as
#'     \eqn{\kappa} passes 1, and near 1 the axis is unstable between samples.}
#'   \item{R-squared with the block and with eligibility}{the squared
#'     correlation of the index's first axis with the filtered block's own first
#'     axis and with eligibility.}
#'   \item{budget check}{if the index ranks the eligible first and the budget
#'     fraction is at most the eligible share, no one outside the filter can
#'     be selected.}
#' }
#' Whether a dominant filter harms is a judgment the audit cannot make: it is benign when
#' passing the filter is part of what the index should measure, and harmful
#' when the filter measures access, as contact with care does in a need index.
#'
#' @inheritParams skip_mca
#' @param codings Codings to compare.
#' @param budget Optional budget fraction for the budget check.
#' @return An object of class `skip_audit`: `codings`, one row per coding;
#'   `items`, the squared correlation ratio of every item with eligibility
#'   under every coding; and the
#'   share asked, the item counts and \eqn{\lambda_X}.
#' @examples
#' dev <- simulate_filtered_survey(n_clusters = 100, seed = 1)
#' items <- grep("^advice", names(dev), value = TRUE)
#' skip_audit(dev, items, paste0("x", 1:6), budget = 0.1)
#' @export
skip_audit <- function(data, items, substantive,
                       codings = c("category", "subset", "no", "fuzzy_average"),
                       eligible = NULL, no_level = NULL, budget = NULL) {
  codings <- unique(match.arg(codings, setdiff(CODINGS, "drop"), several.ok = TRUE))
  data <- as.data.frame(data)
  check_substantive(data, substantive)
  g <- filter_status(data, items, eligible)
  spec <- make_spec(data, items, substantive, g, no_level)
  p <- spec$p
  s_g <- as.numeric(g)
  k <- length(items)
  big_p <- length(substantive)

  x_blk <- index_block(data, spec, g, "drop", include_items = FALSE)
  x_core <- ca_core(x_blk$z, x_blk$q, x_blk$keep)
  lambda_x <- x_core$eig[1]

  item_rows <- list()
  coding_rows <- list()
  for (cd in codings) {
    phi <- phi_formula <- numeric(k)
    for (i in seq_len(k)) {
      v <- items[i]
      b <- item_block(data[[v]], g, cd, spec$levels[[v]], spec$profile[[v]], spec$no_level[[v]], v)
      phi[i] <- eta2(b$z, s_g, b$keep)
      asked <- spec$profile[[v]]
      nl <- length(asked)
      skipped <- switch(cd,
        fuzzy_average = asked,
        no = as.numeric(spec$levels[[v]] == spec$no_level[[v]]),
        c(rep(0, nl), 1))
      if (cd %in% c("category", "subset")) asked <- c(asked, 0)
      phi_formula[i] <- eligibility_eta2(p, asked, skipped, b$keep)
    }
    item_rows[[cd]] <- data.frame(item = items, coding = cd, eta2 = phi, formula = phi_formula)
    a_blk <- index_block(data, spec, g, cd, include_substantive = FALSE)
    a_core <- ca_core(a_blk$z, a_blk$q, a_blk$keep)
    h_a <- ca_project(a_core, a_blk$z)[, 1]
    full <- index_block(data, spec, g, cd)
    f_core <- ca_core(full$z, full$q, full$keep)
    axis <- ca_project(f_core, full$z)[, 1]
    coding_rows[[cd]] <- data.frame(
      coding = cd, eta2 = mean(phi), lambda_items = a_core$eig[1],
      inertia_ratio = k * a_core$eig[1] / (big_p * lambda_x),
      r2_block = stats::cor(axis, h_a)^2, r2_eligibility = stats::cor(axis, s_g)^2)
  }
  structure(list(codings = do.call(rbind, coding_rows), items = do.call(rbind, item_rows),
                 p = p, k = k, P = big_p, lambda_substantive = lambda_x, budget = budget),
            class = "skip_audit")
}

#' @export
print.skip_audit <- function(x, digits = 3, ...) {
  cat(sprintf("Filter audit: %d filtered and %d substantive items; %.1f%% of respondents asked\n",
              x$k, x$P, 100 * x$p))
  cat(sprintf("First eigenvalue of the substantive items alone: %.4f\n\n", x$lambda_substantive))
  tab <- x$codings
  rownames(tab) <- NULL
  num <- vapply(tab, is.numeric, logical(1))
  tab[num] <- lapply(tab[num], round, digits)
  print(tab, row.names = FALSE)
  dominated <- x$codings$coding[x$codings$inertia_ratio >= 1]
  near <- x$codings$coding[abs(x$codings$inertia_ratio - 1) < 0.3]
  cat("\n")
  if (length(dominated)) {
    cat("Inertia ratio at or above 1, the filtered items are likely to dominate the first axis:",
        paste(dominated, collapse = ", "), "\n")
  }
  if (length(near)) {
    cat("Inertia ratio within .3 of 1, the axis may swing between samples:",
        paste(near, collapse = ", "), "\n")
  }
  if (!is.null(x$budget)) {
    if (x$budget <= x$p) {
      cat(sprintf(paste0("Budget %.1f%% is at most the eligible share %.1f%%: an index that ranks ",
                         "the eligible first selects no one outside the filter.\n"),
                  100 * x$budget, 100 * x$p))
    } else {
      cat(sprintf(paste0("Budget %.1f%% exceeds the eligible share %.1f%%: at most %.0f%% of slots ",
                         "can go outside the filter if the eligible are ranked first.\n"),
                  100 * x$budget, 100 * x$p, 100 * (x$budget - x$p) / x$budget))
    }
  }
  cat("Is passing the filter part of what the index should measure? If not, prefer",
      "missing fuzzy average or drop the items.\n")
  invisible(x)
}
