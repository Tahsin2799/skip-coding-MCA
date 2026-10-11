# Separation must remove only the affected regression metric.
# Rscript JCE/analysis/test_review_analysis.R
source("methodology/v2/summarize_convergence.R")
correlate <- family_correlation
if (file.exists("JCE/analysis/review_analysis.R")) {
  source("JCE/analysis/review_analysis.R")
  correlate <- corrected_family_correlation
}
fixture <- tibble(cell = "one", failed = FALSE,
                  separation_dev = c(FALSE, FALSE, FALSE, TRUE),
                  separation_target = FALSE,
                  log_or_top = c(1, 2, 3, -100),
                  auc_index_scored = c(1, 2, 3, 0),
                  cluster_rho_all = c(1, 2, 3, 4))
ans <- correlate(fixture, apparent = c("log_or_top", "auc_index_scored"),
                 targeting = "cluster_rho_all")
stopifnot("separated regression must not reverse a valid association" =
            abs(ans$median_rho[ans$a == "log_or_top"] - 1) < 1e-12,
          "development separation must not remove target AUC" =
            abs(ans$median_rho[ans$a == "auc_index_scored"] + .2) < 1e-12)
fixture$separation_dev <- FALSE
fixture$separation_target <- c(FALSE, FALSE, FALSE, TRUE)
fixture$log_or_top <- c(1, 2, 3, 0)
ans <- correlate(fixture, apparent = c("log_or_top", "auc_index_scored"),
                 targeting = "cluster_rho_all")
stopifnot("target separation must mask target AUC" =
            abs(ans$median_rho[ans$a == "auc_index_scored"] - 1) < 1e-12,
          "target separation must preserve valid development metrics" =
            abs(ans$median_rho[ans$a == "log_or_top"] + .2) < 1e-12)
cat("Separation-filter regression test passed.\n")
