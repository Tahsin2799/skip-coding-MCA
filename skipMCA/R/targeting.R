#' Allocation outcomes of an index against a known need criterion
#'
#' Scores an index the way a programme would use it, over everyone it serves.
#' Unscored people (`NA`) are never selected. Ties are broken by row order.
#'
#' Individual budget: the top `budget` fraction of all people by score are
#' selected; reported are the share of the high-need group (the top
#' `high_need` fraction by `truth`) selected, which is one minus the exclusion
#' error, and the share of slots used.
#'
#' Cluster budget, when `cluster` is given: a cluster's priority is the mean
#' score of its scored residents and its truth the mean of `truth` over all
#' residents. Clusters with no scored resident rank last. Reported are the
#' Spearman rank correlation between priority and truth, the inclusion error
#' (share of selected clusters below the median truth) and the top-decile
#' precision (share of selected clusters in the true top tenth).
#'
#' @param score Index scores, higher meaning higher priority.
#' @param truth The need criterion, measured independently of the index.
#' @param cluster Optional cluster identifiers.
#' @param budget Budget fraction, for people and for clusters.
#' @param high_need Fraction of people, by `truth`, who form the high-need group.
#' @return A named numeric vector.
#' @references Brown C, Ravallion M, van de Walle D (2018). A poor means test?
#'   Econometric targeting in Africa. *Journal of Development Economics* 134,
#'   109-124.
#' @examples
#' set.seed(1)
#' truth <- rnorm(1000)
#' cluster <- rep(1:50, each = 20)
#' targeting_metrics(truth + rnorm(1000), truth, cluster)
#' @export
targeting_metrics <- function(score, truth, cluster = NULL, budget = 0.1, high_need = 0.2) {
  n <- length(score)
  if (length(truth) != n || anyNA(truth)) stop("`truth` needs one non-missing value per person.",
                                                call. = FALSE)
  slots <- floor(budget * n)
  scored <- which(!is.na(score))
  ranked <- scored[order(-score[scored], scored)]
  picked <- utils::head(ranked, slots)
  hn <- truth >= stats::quantile(truth, 1 - high_need, type = 1, names = FALSE)
  out <- c(high_need_selected = sum(hn[picked]) / sum(hn),
           slots_used = if (slots > 0) length(picked) / slots else NA_real_,
           share_scored = length(scored) / n)
  if (is.null(cluster)) return(out)

  if (length(cluster) != n) stop("`cluster` needs one value per person.", call. = FALSE)
  ids <- unique(cluster)
  priority <- vapply(ids, function(j) {
    s <- score[cluster == j]
    if (all(is.na(s))) NA_real_ else mean(s, na.rm = TRUE)
  }, numeric(1))
  need <- vapply(ids, function(j) mean(truth[cluster == j]), numeric(1))
  ranking <- priority
  ranking[is.na(ranking)] <- min(priority, na.rm = TRUE) - 1
  n_sel <- floor(budget * length(ids))
  sel_pool <- which(!is.na(priority))
  sel <- utils::head(sel_pool[order(-priority[sel_pool], sel_pool)], n_sel)
  top_decile <- need >= stats::quantile(need, 0.9, type = 1, names = FALSE)
  c(out,
    cluster_rank_correlation = stats::cor(ranking, need, method = "spearman"),
    inclusion_error = if (length(sel)) mean(need[sel] < stats::median(need)) else NA_real_,
    top_decile_precision = if (length(sel)) mean(top_decile[sel]) else NA_real_)
}
