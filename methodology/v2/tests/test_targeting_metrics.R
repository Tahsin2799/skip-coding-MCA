# Denominator and budget contracts, on fixtures small enough to check by hand.
# Run: Rscript methodology/v2/tests/test_targeting_metrics.R

suppressPackageStartupMessages({library(dplyr); library(tibble)})
root <- c(".", "..", "../..") |>
  Filter(f = \(p) file.exists(file.path(p, "methodology/v2/plasmode_generate.R"))) |>
  (\(x) x[1])()
for (f in c("harmonize_donors", "plasmode_generate", "index_arms")) {
  source(file.path(root, sprintf("methodology/v2/%s.R", f)))
}

# Twenty people, V = 1:20, two per cluster, ten clusters. The top fifth of V is
# rows 17 to 20, which is clusters 9 and 10. Rows 19 and 20 are the only
# unrecognized cases.
fixture <- function(score) {
  tgt <- tibble(
    row_id = 1:20, V = as.numeric(1:20), age = rep(c(25, 35, 50, 70), 5),
    sex = factor(rep(c("male", "female"), 10), levels = c("male", "female")),
    y = rep(c(0L, 1L), 10), recognized = 0L,
    cluster = paste0("c", rep(1:10, each = 2)),
    tiebreak = seq(0.01, 0.20, by = 0.01),
    cluster_tiebreak = rep(seq(0.1, 1, by = 0.1), each = 2))
  tgt$y <- c(rep(0L, 18), 1L, 1L)
  prep <- list(
    n_fix = 0L,
    dt = list(breaks = c(0, 10, 20), flip = FALSE, n_tiers = 2L),
    dm = list(separation = FALSE,
              row = tibble(or_top = NA_real_, n_tiers = 2L, n_dev_used = 0L,
                           beta_within = NA_real_, beta_between = NA_real_)))
  list(tgt = tgt, prep = prep, out = target_metrics(prep, list(score = score, n_fix = 0L),
                                                    tgt, target_truth(tgt)))
}

# Everyone in clusters 1 to 9 is scored; cluster 10 cannot be scored at all.
score_9 <- c(as.numeric(1:18), NA, NA)
f <- fixture(score_9)
o <- f$out

stopifnot("coverage counts every target person, not the arm's own frame" =
            o$coverage == 0.9 && o$n_scored == 18L && o$N == 20L)
stopifnot("the high-need denominator is the whole target" = o$n_hn == 4L)
stopifnot("the unrecognized-case denominator is counted too" = o$n_unrec == 2L)

# Tier 2 is score above 10, so rows 11 to 18 are flagged: eight people, of whom
# rows 17 and 18 are high need.
stopifnot("the tier policy flags the right share" = o$tier_flag_rate == 8 / 20)
stopifnot("whole-target capture divides by everyone in the true top fifth" =
            o$capture_all == 2 / 4)
stopifnot("unscored high-need people count against the arm" = o$capture_all < 1)

# Ten percent of twenty people is two slots; the two highest scores are rows 18
# and 17, both high need, which is the best any arm could do at this budget.
stopifnot("the budget is floor(b * N)" = o$b10_slots == 2L && o$b20_slots == 4L && o$b05_slots == 1L)
stopifnot("both slots are used and both find a high-need person" =
            o$b10_used == 2L && o$b10_hn == 2L && o$b10_hn_capture == 0.5 && o$b10_hn_ppv == 1)
stopifnot("the shortfall against the ideal allocation is zero here" =
            o$b10_hn_shortfall_per1000 == 0)
stopifnot("neither unrecognized case is reachable, both being in cluster 10" =
            o$b10_unrec == 0L && o$b10_unrec_capture == 0)

# Cluster targeting: nine of ten clusters can be ranked, and the tenth, which is
# the truly most deprived, can never be selected.
stopifnot("cluster score coverage is reported" =
            o$cluster_score_coverage == 0.9 && o$cluster_resident_coverage == 0.9)
stopifnot("one cluster slot at ten percent of ten clusters" = o$cl10_slots == 1L)
stopifnot("the selected cluster is the best one the arm can see, not the truly worst" =
            o$cl10_used == 1L && o$cl10_precision == 0)
stopifnot("that cluster is above the median, so it is not a misallocation" =
            o$cl10_misalloc == 0)
stopifnot("half the high-need people live in the selected cluster" =
            o$cl10_hn_reached == 0.5 && o$cl10_residents == 2L)
stopifnot("an unscored cluster ranks below every scored one, so the correlation is not 1" =
            o$cluster_rho_all < 1)

# Age bands keep unscored people in the denominator and report the count even
# when it is too small to yield a rate.
stopifnot("band denominators are always reported" =
            o$n_hn_18_29 + o$n_hn_30_44 + o$n_hn_45_59 + o$n_hn_60p == 4L)
stopifnot("a band below the minimum is unavailable, not zero" =
            is.na(o$cap_tier_18_29) && is.na(o$cap_b10_60p))

# ---- Unequal coverage across bands must show in whole-target capture ------
# Everyone is scored identically, except that the arm cannot score the two
# oldest high-need people. Scored-only capture cannot see the difference;
# whole-target capture must.
drop_old <- score_9
drop_old[17:18] <- NA
g <- fixture(drop_old)
stopifnot("losing high-need people lowers whole-target capture" =
            g$out$capture_all == 0 && g$out$coverage == 0.8)

# ---- A wholly unscored target -------------------------------------------
none <- fixture(rep(NA_real_, 20))$out
stopifnot("no scores means no coverage and no flags" =
            none$coverage == 0 && none$tier_flag_rate == 0 && none$capture_all == 0)
stopifnot("cluster targeting is unavailable rather than zero" =
            is.na(none$cluster_rho_all) && is.na(none$cl10_misalloc) && none$cl10_used == 0L)
stopifnot("the budget reports unused capacity" =
            none$b10_slots == 2L && none$b10_used == 0L &&
            none$b10_hn_shortfall_per1000 == 1000 * 2 / 20)

# ---- An empty subgroup is unavailable, not a zero -------------------------
recognized_all <- fixture(score_9)
t2 <- recognized_all$tgt
t2$recognized <- 1L
o2 <- target_metrics(recognized_all$prep, list(score = score_9, n_fix = 0L), t2, target_truth(t2))
stopifnot("with no unrecognized case the capture is NA, not 0" =
            o2$n_unrec == 0L && is.na(o2$unrec_capture_tier) && is.na(o2$b10_unrec_capture))

# ---- Ties are broken blind to the truth -----------------------------------
# Every score is identical, so only the tie-break decides. Changing the truth
# must not change who is selected.
tied <- fixture(rep(5, 20))
t3 <- tied$tgt
t3$V <- rev(t3$V)
t3$y <- rev(t3$y)
sel_a <- select_top(rep(5, 20), tied$tgt$tiebreak, 3L)
sel_b <- select_top(rep(5, 20), t3$tiebreak, 3L)
stopifnot("the same scores and tie-breaks select the same people" = identical(sel_a, sel_b))
stopifnot("exactly the budget is filled" = length(sel_a) == 3L)
stopifnot("selection uses the tie-break, not the row order" =
            identical(sel_a, order(tied$tgt$tiebreak)[1:3]))

# ---- Fewer scoreable people than slots ------------------------------------
thin <- c(rep(NA_real_, 17), 1, 2, 3)
o4 <- fixture(thin)$out
stopifnot("a pure complete-case policy leaves unfilled capacity" =
            o4$b20_slots == 4L && o4$b20_used == 3L)
stopifnot("and the unfilled slots are visible in the shortfall" =
            o4$b20_hn_shortfall_per1000 > 0)

cat("test_targeting_metrics.R: all contracts hold\n")
