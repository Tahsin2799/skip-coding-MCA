#' Fit an MCA index with a chosen coding of questionnaire skips
#'
#' Builds a multiple correspondence analysis (MCA) index from substantive
#' categorical items and a block of filtered items that only respondents who
#' passed one filter question were asked. Skipped answers are `NA`.
#'
#' The codings are:
#' \describe{
#'   \item{`"fuzzy_average"`}{missing fuzzy average (van der Heijden and
#'     Escofier 2003): a skipped respondent receives the answer distribution of
#'     those asked. Its squared correlation ratio with eligibility is zero, and
#'     it is the only fixed profile with that property.}
#'   \item{`"subset"`}{subset (specific) MCA: the skip is a category whose
#'     column is left out of the analysis while keeping the full-data masses.
#'     It gives the same solution as missing passive modified margin and keeps
#'     a squared correlation ratio with eligibility of \eqn{1 - p} per item.}
#'   \item{`"category"`}{missing single, the skip as a category of its own;
#'     squared correlation ratio with eligibility 1.}
#'   \item{`"no"`}{the skip coded as the reference answer `no_level`; squared
#'     correlation ratio with eligibility \eqn{\pi(1-p)/(1-p\pi)} for a binary
#'     item answered "yes" with probability \eqn{\pi}.}
#'   \item{`"drop"`}{the filtered items are left out.}
#' }
#'
#' Rows are placed by the transition formula on the development column
#' standard coordinates, which reproduces the fitted row coordinates exactly,
#' so development and new respondents are scored the same way.
#'
#' @param data A data frame, the development sample.
#' @param items Names of the filtered items. `NA` marks a skip; every eligible
#'   respondent must have answered every item.
#' @param substantive Names of the other categorical items, without missing
#'   values.
#' @param coding One of `"fuzzy_average"`, `"subset"`, `"category"`, `"no"` or
#'   `"drop"`.
#' @param eligible Optional logical vector, or the name of a logical column,
#'   marking who passed the filter. By default, respondents who answered the
#'   first filtered item.
#' @param no_level For `coding = "no"`, the answer given to the skipped, one
#'   value or one per item. Defaults to the first observed level of each item.
#' @param ncp Number of dimensions to keep.
#' @param orient Optional numeric vector, one value per row of `data`. If
#'   given, each axis is signed so that it correlates positively with it, for
#'   example a disease indicator in development. The sign of an MCA axis is
#'   otherwise arbitrary.
#'
#' @return An object of class `skip_mca` with the development scores
#'   (`scores`), the eigenvalues (`eig`), the share asked (`p`) and what
#'   [predict.skip_mca()] needs.
#' @seealso [skip_audit()] to decide whether the filtered items will dominate
#'   the index under a coding; [predict.skip_mca()]; [targeting_metrics()].
#' @references van der Heijden PGM, Escofier B (2003). Multiple correspondence
#'   analysis with missing data. In Escofier B (ed.), *Analyse des
#'   correspondances: recherches au coeur de l'analyse des donnees*, 152-170.
#'   Presses Universitaires de Rennes.
#'
#'   Greenacre M, Pardo R (2006). Subset correspondence analysis. *Sociological
#'   Methods & Research* 35, 193-218.
#'
#'   Josse J, Chavent M, Liquet B, Husson F (2012). Handling missing values
#'   with regularized iterative multiple correspondence analysis. *Journal of
#'   Classification* 29, 91-116.
#' @examples
#' dev <- simulate_filtered_survey(n_clusters = 100, seed = 1)
#' items <- grep("^advice", names(dev), value = TRUE)
#' inputs <- paste0("x", 1:6)
#' fit <- skip_mca(dev, items, inputs, coding = "fuzzy_average", orient = dev$need)
#' fit
#' @export
skip_mca <- function(data, items, substantive, coding = "fuzzy_average", eligible = NULL,
                     no_level = NULL, ncp = 1, orient = NULL) {
  coding <- match.arg(coding, CODINGS)
  data <- as.data.frame(data)
  check_substantive(data, substantive)
  g <- filter_status(data, items, eligible)
  spec <- make_spec(data, items, substantive, g, no_level)
  blk <- index_block(data, spec, g, coding)
  core <- ca_core(blk$z, blk$q, blk$keep, ncp)
  scores <- ca_project(core, blk$z)
  if (!is.null(orient)) {
    if (length(orient) != nrow(data)) stop("`orient` needs one value per row.", call. = FALSE)
    flip <- apply(scores, 2, function(s) {
      r <- stats::cor(s, orient, use = "complete.obs")
      if (is.na(r) || r >= 0) 1 else -1
    })
    core$gamma <- sweep(core$gamma, 2, flip, "*")
    scores <- sweep(scores, 2, flip, "*")
  }
  structure(list(coding = coding, items = items, substantive = substantive, spec = spec,
                 core = core, eig = core$eig, scores = scores, p = spec$p, eligible = g,
                 n = nrow(data)),
            class = "skip_mca")
}

#' Score new respondents on a fitted index
#'
#' Places respondents on the development axes by the transition formula.
#' Under missing fuzzy average, the skipped in `newdata` receive by default the
#' answer distribution of the respondents asked in `newdata` itself, so the
#' filtered columns carry no filter wherever the index is used; no model is
#' fitted.
#'
#' @param object A `skip_mca` fit.
#' @param newdata A data frame with the same items. If missing, the
#'   development scores are returned.
#' @param eligible As in [skip_mca()], for `newdata`.
#' @param profile For missing fuzzy average, `"target"` (the answers of those
#'   asked in `newdata`) or `"development"`.
#' @param ... Unused.
#' @return A matrix with one column per dimension. Rows with a category the
#'   development sample never saw are `NA`.
#' @examples
#' dev <- simulate_filtered_survey(n_clusters = 100, contact = 0.5, seed = 1)
#' target <- simulate_filtered_survey(n_clusters = 100, contact = 0.2, seed = 2)
#' items <- grep("^advice", names(dev), value = TRUE)
#' fit <- skip_mca(dev, items, paste0("x", 1:6), orient = dev$need)
#' head(predict(fit, target))
#' @export
predict.skip_mca <- function(object, newdata, eligible = NULL,
                             profile = c("target", "development"), ...) {
  if (missing(newdata)) return(object$scores)
  profile <- match.arg(profile)
  newdata <- as.data.frame(newdata)
  spec <- object$spec
  check_substantive(newdata, spec$substantive)
  g <- filter_status(newdata, spec$items, eligible)
  prof <- spec$profile
  if (object$coding == "fuzzy_average" && profile == "target") {
    if (any(g)) {
      prof <- asked_profile(newdata, spec$items, spec$levels, g)
    } else {
      warning("No one in `newdata` passed the filter; using the development profile.",
              call. = FALSE)
    }
  }
  blk <- index_block(newdata, spec, g, object$coding, profile = prof)
  out <- ca_project(object$core, blk$z)
  out[blk$unseen, ] <- NA_real_
  out
}

#' @export
print.skip_mca <- function(x, ...) {
  cat("MCA index, skip coding:", x$coding, "\n")
  cat(sprintf("%d substantive and %d filtered items; %d respondents, %.1f%% asked\n",
              length(x$substantive), length(x$items), x$n, 100 * x$p))
  cat("First eigenvalues:", paste(formatC(utils::head(x$eig, 3), format = "f", digits = 4),
                                 collapse = ", "), "\n")
  invisible(x)
}
