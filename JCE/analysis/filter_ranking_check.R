# Cluster ranking under a dominant filter, quoted in section 4.1 and supplement
# A6.11. Reads the person-level artifacts production saved for two replicates of
# every configuration (24 per filter in the core design) and compares each
# skip coding's cluster ranking with the ranking given by the cluster
# eligibility rate itself. Nothing is written.
#
#   Rscript JCE/analysis/filter_ranking_check.R

dir <- "outputs/methodology/production_v03/artifacts"
files <- list.files(dir, full.names = TRUE,
  pattern = "^(Q3|bd_allage)_nd1000_nt20000_r[468]0_q[05]0_g(con|dia)_base_rep[0-9]+\\.rds$")
stopifnot(length(files) == 48)

spearman <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))

# Cluster priority as in verify_from_artifacts.R: mean score of the scored
# residents, clusters with none ranked below every other.
cluster_priority <- function(a, cid) {
  s <- if (a$flip) -a$score else a$score
  cs <- tapply(s, cid, function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE))
  ifelse(is.na(cs), min(cs, na.rm = TRUE) - 1, cs)
}

# Share of the 10% individual budget given to diagnosed people, with the
# verifier's ordering and tie-break.
diagnosed_share <- function(a, truth, acc) {
  s <- if (a$flip) -a$score else a$score
  cand <- which(!is.na(s))
  cand <- cand[order(-s[cand], truth$tiebreak[cand])]
  sel <- head(cand, floor(0.10 * nrow(truth)))
  mean(acc$recognized[sel] == 1L)
}

codings <- c(single = "category", subset = "passive_mca", no = "zero")
rows <- list()
for (path in files) {
  art <- readRDS(path)
  cid <- as.integer(factor(art$truth$cluster))
  need <- tapply(art$truth$V, cid, mean)
  for (st in c("c20", "c50", "c80")) {
    acc <- art$access[[st]]
    rate <- tapply(acc$eligible, cid, mean)
    none <- rate == 0
    pr <- lapply(codings, function(m) cluster_priority(art$scores[[paste0(m, "_", st)]], cid))
    rows[[length(rows) + 1]] <- data.frame(
      filter = if (grepl("_gdia_", path)) "diagnosis" else "contact",
      contact = st,
      share_eligible = mean(acc$eligible),
      clusters_none_eligible = mean(none),
      rho_eligibility_rate = spearman(rate, need),
      rho_single = spearman(pr$single, need),
      rho_subset = spearman(pr$subset, need),
      rho_no = spearman(pr$no, need),
      rho_single_none_eligible = if (sum(none) > 10) spearman(pr$single[none], need[none]) else NA_real_,
      diagnosed_selected_single = diagnosed_share(art$scores[[paste0("category_", st)]], art$truth, acc),
      diagnosed_selected_subset = diagnosed_share(art$scores[[paste0("passive_mca_", st)]], art$truth, acc),
      diagnosed_selected_no = diagnosed_share(art$scores[[paste0("zero_", st)]], art$truth, acc),
      diagnosed_selected_dropped = diagnosed_share(art$scores[[paste0("drop_items_", st)]], art$truth, acc))
  }
}
d <- do.call(rbind, rows)

means <- aggregate(. ~ filter + contact, data = d, FUN = mean, na.action = na.pass)
ranges <- aggregate(cbind(rho_eligibility_rate, rho_single) ~ filter + contact, data = d, FUN = range)
print(format(means, digits = 3), row.names = FALSE)
print(ranges, row.names = FALSE, digits = 3)
