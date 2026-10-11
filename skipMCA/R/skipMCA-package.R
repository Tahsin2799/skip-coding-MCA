#' skipMCA: questionnaire skips in multiple correspondence analysis indices
#'
#' A block of items asked only of respondents who pass one filter question
#' shares a pattern, who was asked, whatever the answers. Multiple
#' correspondence analysis builds its first axis from the strongest shared
#' pattern, so the coding of the skip decides how much of the filter enters an
#' index and whether the filter takes the index over.
#'
#' * [skip_audit()]: each coding's squared correlation ratio with eligibility,
#'   the inertia ratio, the R-squared with eligibility and a budget check, from
#'   development data alone.
#' * [skip_mca()] and [predict.skip_mca()]: fit an index under missing fuzzy
#'   average, subset MCA, missing single (the skip as its own category) or the
#'   skip as "no", and score
#'   new respondents.
#' * [eligibility_eta2()]: the closed form of Theorem 1.
#' * [targeting_metrics()]: allocation outcomes against a need criterion.
#' * [simulate_filtered_survey()]: example data.
#'
#' @keywords internal
"_PACKAGE"
