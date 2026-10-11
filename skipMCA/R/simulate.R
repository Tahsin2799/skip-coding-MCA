#' Simulate a survey with items asked only after a filter
#'
#' A small generator for examples and tests, after the plasmode design of the
#' paper but with synthetic covariates. Six three-level deprivation items
#' depend weakly on a cluster and an individual deprivation, so that, as in
#' the survey populations of the paper, their first MCA eigenvalue is near .27;
#' latent need combines them with noise. Contact with care depends on cluster
#' access, which is unrelated to need, and slightly on need; disease depends on
#' need; recognition requires contact; advice requires contact and is likelier
#' once recognized. Advice items are asked only of the contacted (`filter =
#' "contact"`) or only of the recognized (`filter = "diagnosis"`), and are `NA`
#' for everyone else.
#'
#' @param n_clusters Number of clusters.
#' @param cluster_size People per cluster.
#' @param contact Approximate share contacted.
#' @param filter `"contact"` or `"diagnosis"`.
#' @param n_items Number of advice items.
#' @param seed Optional seed.
#' @return A data frame with `cluster`, `need` (the latent need, the truth for
#'   targeting), `disease`, `contacted`, `recognized`, `eligible`, the
#'   substantive items `x1` to `x6` and the advice items.
#' @examples
#' d <- simulate_filtered_survey(n_clusters = 20, seed = 1)
#' table(asked = !is.na(d$advice1), d$contacted)
#' @export
simulate_filtered_survey <- function(n_clusters = 200, cluster_size = 25, contact = 0.5,
                                     filter = c("contact", "diagnosis"), n_items = 4,
                                     seed = NULL) {
  filter <- match.arg(filter)
  if (!is.null(seed)) set.seed(seed)
  n <- n_clusters * cluster_size
  cluster <- rep(seq_len(n_clusters), each = cluster_size)
  depriv <- stats::rnorm(n_clusters)[cluster] + stats::rnorm(n)
  x <- vapply(1:6, function(i) {
    latent <- 0.3 * depriv + stats::rnorm(n)
    as.integer(cut(latent, stats::quantile(latent, c(0, 1 / 3, 2 / 3, 1)), include.lowest = TRUE))
  }, integer(n))
  need <- as.numeric(scale(rowMeans(x))) * sqrt(0.6) + stats::rnorm(n, sd = sqrt(0.4))
  disease <- stats::rbinom(n, 1, stats::plogis(-1.4 + 0.8 * need))
  access <- stats::rnorm(n_clusters)[cluster]
  contacted <- stats::runif(n) < stats::plogis(stats::qlogis(contact) + 0.75 * access + 0.3 * need)
  recognized <- contacted & stats::runif(n) < ifelse(disease == 1, 0.5, 0.05)
  eligible <- if (filter == "contact") contacted else recognized
  base <- rep(c(0.10, 0.10, 0.05, 0.15), length.out = n_items)
  advice <- vapply(seq_len(n_items), function(k) {
    advised <- contacted &
      stats::runif(n) < stats::plogis(stats::qlogis(base[k]) + log(16) * recognized + 0.25 * need)
    ifelse(eligible, ifelse(advised, "yes", "no"), NA_character_)
  }, character(n))
  out <- data.frame(cluster = cluster, need = need, disease = disease, contacted = contacted,
                    recognized = recognized, eligible = eligible)
  for (i in 1:6) out[[paste0("x", i)]] <- factor(x[, i], levels = 1:3, labels = c("low", "mid", "high"))
  for (k in seq_len(n_items)) out[[paste0("advice", k)]] <- factor(advice[, k], levels = c("no", "yes"))
  out
}
