# Index and arm contracts. Run: Rscript methodology/v2/tests/test_index_arms.R
#
# These check implementation identities, not scientific conclusions. A test that
# requires a favoured arm to win would pass by construction and tell us nothing.

suppressPackageStartupMessages({library(dplyr); library(tibble); library(digest)})
root <- c(".", "..", "../..") |>
  Filter(f = \(p) file.exists(file.path(p, "methodology/v2/plasmode_generate.R"))) |>
  (\(x) x[1])()
for (f in c("harmonize_donors", "plasmode_generate", "index_arms",
            "quadrant_populations", "pilot_v2_config")) {
  source(file.path(root, sprintf("methodology/v2/%s.R", f)))
}

pools <- readRDS(file.path(root, "outputs/methodology/donors_harmonized_geo.rds")) |>
  lapply(prefix_levels) |>
  build_quadrants(quadrant_cells(file.path(root, "methodology/v2/quadrant_cells.csv"))) |>
  lapply(\(d) droplevels(as.data.frame(d[, c(INDEX_INPUTS, "age", "sex", "cluster")])))
spec <- pilot_spec(population = "Q1")
set.seed(digest2int("test_arms_1"))
sim <- generate_replicate(pools[[1]], spec, 1000L, 20000L)

prep_of <- function(arm, s = spec) prepare_family(sim$development, prep_family(arm), s)

# ---- The arm registry ---------------------------------------------------
stopifnot("ten arms in total" = length(ALL_ARMS) == 10L,
          "the two imputation arms share one development half" =
            prep_family("impute_frozen") == prep_family("impute_adaptive"),
          "every arm has a family" =
            all(vapply(ALL_ARMS, \(a) is.character(prep_family(a)), logical(1))))

# ---- The metric frame is insulated from the data frame -------------------
# THE regression test. dev_metrics builds `w`, `zbar` and `score` for the
# Mundlak fit. If it ever again adds them to the development frame instead of a
# private one, a generator column of the same name wins and the slope measures
# that column.
z <- prep_of("zero")
clean <- dev_metrics(z$dev, z$dt)
decoy <- z$dev
for (nm in c("z", "w", "zbar", "score")) decoy[[nm]] <- rnorm(nrow(decoy))
stopifnot("metrics ignore decoy columns named like model terms" =
            isTRUE(all.equal(as.list(clean$row), as.list(dev_metrics(decoy, z$dt)$row))))
stopifnot("the within slope is not zero for a working index" =
            abs(clean$row$beta_within) > 0.1)

# ---- Separation is detected ---------------------------------------------
fake <- z$dev
fake$y <- as.integer(z$dt$tier_dev >= 4L)
stopifnot("a perfectly separating development fit is flagged" =
            isTRUE(dev_metrics(fake, z$dt)$separation))
stopifnot("an ordinary fit is not flagged" = isFALSE(clean$separation))

# ---- Imputation survives a level absent from the eligible subset ---------
# 231 of 4,000 replicates were lost to this, and not at random: the draws that
# failed were the structurally heterogeneous ones.
d <- sim$development
tw <- table(d$water)
rare <- names(sort(tw[tw > 0]))[1]
d$eligible[d$water == rare] <- 0L
for (v in advice_names(d)) d[[v]][d$eligible == 0L] <- NA_integer_
mod <- try(item_models(d), silent = TRUE)
stopifnot("the item models still fit" = !inherits(mod, "try-error"))
out <- try(complete_items(d, mod, 2L), silent = TRUE)
stopifnot("completion does not abort on an unseen level" = !inherits(out, "try-error"))
stopifnot("every row is retained" = nrow(out$frames[[1]]) == nrow(d))
stopifnot("no advice item is left missing" =
            !any(vapply(advice_names(d), \(v) anyNA(out$frames[[1]][[v]]), logical(1))))
stopifnot("the substitution count is reported" = out$n_fix > 0L)
stopifnot("observed responses are never replaced" =
            all(out$frames[[1]][[advice_names(d)[1]]][d$eligible == 1L] ==
                  item_label(d[[advice_names(d)[1]]], advice_names(d)[1])[d$eligible == 1L]))

# ---- A constant input does not cost the replicate ------------------------
# Bangladesh 18-95 is 92 percent one water rung, so a 1,000-row draw can hold a
# single level and `glm` aborts on a one-level factor.
deg <- sim$development
deg$water <- factor(rep(levels(droplevels(deg$water))[1], nrow(deg)))
stopifnot("a constant column is not a model term" =
            !("water" %in% varying(deg, INDEX_INPUTS)),
          "and a frame with no varying term still yields a formula" =
            identical(varying(deg[, "water", drop = FALSE], "water"), "1"))
for (arm in c("impute_frozen", "raw_items")) {
  p <- try(prepare_family(deg, prep_family(arm), spec), silent = TRUE)
  stopifnot("a degenerate input does not abort the arm" = !inherits(p, "try-error"))
  r <- try(score_prepared(p, arm, sim$targets$c50, spec), silent = TRUE)
  stopifnot("and it still scores the target" =
              !inherits(r, "try-error") && is.finite(r$cluster_rho_all))
}

# ---- A single-sex population is a population, not a failure --------------
# Nepal 2022 has 6,510 women and no men in the 18-49 window once the six inputs
# are complete, and `y ~ tier + age + sex` aborts on that for every arm.
one_sex <- lapply(list(dev = sim$development, tgt = sim$targets$c50), function(d) {
  d <- d[d$sex == levels(droplevels(d$sex))[1], , drop = FALSE]
  d$sex <- droplevels(d$sex)
  d$row_id <- seq_len(nrow(d))
  d
})
stopifnot("the fixture really is single sex" = nlevels(one_sex$dev$sex) == 1L)
for (arm in ALL_ARMS) {
  r <- try(score_arm(one_sex$dev, one_sex$tgt, arm, spec), silent = TRUE)
  stopifnot("every arm scores a single-sex population" = !inherits(r, "try-error"),
            "and returns a usable targeting estimate" = is.finite(r$cluster_rho_all))
}

# ---- Development is common to the three access settings ------------------
for (arm in c("impute_frozen", "impute_adaptive", "raw_items", "zero")) {
  prep <- prep_of(arm)
  rows <- lapply(spec$settings, \(s) score_prepared(prep, arm, sim$targets[[s]], spec))
  for (col in c("or_top", "beta_within", "beta_between", "n_dev_used")) {
    v <- vapply(rows, \(r) r[[col]], numeric(1))
    stopifnot("a development quantity is identical across access settings" =
                isTRUE(all.equal(v[1], v[2])) && isTRUE(all.equal(v[1], v[3])))
  }
}

# ---- The advice-free arms cannot depend on access ------------------------
# They read no gated item and the outcome is imposed before access, so these are
# the negative controls. If one moves, the access contrast is contaminated.
for (arm in ADVICE_FREE_ARMS) {
  r <- lapply(spec$settings, \(s) score_arm(sim$development, sim$targets[[s]], arm, spec))
  same <- function(col) length(unique(vapply(r, \(x) x[[col]], numeric(1)))) == 1L
  stopifnot("an advice-free arm is invariant to the access setting" =
              all(vapply(c("cluster_rho_all", "capture_all", "b10_hn_capture",
                           "coverage", "cl10_misalloc"), same, logical(1))))
}

# ---- Coverage and the common estimand ------------------------------------
cc20 <- score_arm(sim$development, sim$targets$c20, "complete_case", spec)
full <- score_arm(sim$development, sim$targets$c50, "drop_items", spec)
stopifnot("complete_case coverage tracks the share the gate admits" =
            cc20$coverage > 0.1 && cc20$coverage < 0.35)
stopifnot("a full-coverage arm reads exactly 1" = full$coverage == 1)
stopifnot("scored and whole-target metrics agree when every row is scored" =
            isTRUE(all.equal(full$capture_scored, full$capture_all)) &&
            isTRUE(all.equal(full$cluster_rho_scored, full$cluster_rho_all)) &&
            isTRUE(all.equal(full$misallocation_scored, full$misallocation_all)))
stopifnot("an arm that cannot score a person cannot flag them" =
            cc20$capture_all <= cc20$capture_scored)
stopifnot("denominators are reported alongside the rates" =
            full$n_hn > 0 && full$N == nrow(sim$targets$c50) &&
            full$n_scored == full$N)

# ---- Tier orientation and the tier scale ---------------------------------
# Dropping an input silently flipped the dim-1 sign in v1 and produced a
# Tier-5 OR of 0.06. The auto-orient guard must survive.
for (arm in c("zero", "drop_items", "complete_case", "simple_score", "passive_mca")) {
  p <- prep_of(arm)
  pr <- tapply(p$dev$y, p$dt$tier_dev, mean)
  stopifnot("the top tier carries more outcome than the bottom" =
              pr[[length(pr)]] >= pr[[1]])
  ix <- p$project(sim$targets$c50, arm)
  tier <- project_tiers(p$dt, ix$score[!is.na(ix$score)])$tier
  stopifnot("every scoreable target row gets a tier" = !anyNA(tier))
}
stopifnot("a score with heavy ties reports how many tiers it supports" =
            prep_of("simple_score")$dt$n_tiers <= 5L)

# ---- The simple score is deployable and reads no advice item -------------
ss <- prep_of("simple_score")
stopifnot("the simple score is a count of six deprivations" =
            all(ss$s_dev %in% 0:6),
          "it scores every target row" =
            !anyNA(ss$project(sim$targets$c20, "simple_score")$score),
          "and it uses no advice item" =
            !any(grepl("advice", SIMPLE_DEPRIVED)))

# ---- The oracle carries development scaling to the target ----------------
# Standardizing the target on itself preserves ranks but moves people across the
# fixed development cuts, which is a different policy.
orc <- prep_of("oracle")
s_tgt <- orc$project(sim$targets$c50, "oracle")$score
raw_dev <- as.numeric(rung_matrix(sim$development, spec$bmi_risk) %*% spec$loadings)
raw_tgt <- as.numeric(rung_matrix(sim$targets$c50, spec$bmi_risk) %*% spec$loadings)
stopifnot("the target is standardized with development constants" =
            isTRUE(all.equal(s_tgt, (raw_tgt - mean(raw_dev)) / stats::sd(raw_dev))))
stopifnot("so the target mean is not forced to zero" = abs(mean(s_tgt)) > 1e-10)
stopifnot("and the development score is its own standardization" =
            isTRUE(all.equal(abs(orc$dt$score_dev), abs((raw_dev - mean(raw_dev)) / stats::sd(raw_dev)))))

# ---- Level fixups are counted, not silent --------------------------------
imp <- score_arm(sim$development, sim$targets$c50, "impute_adaptive", spec)
stopifnot("the substitution count is carried into the row" =
            is.numeric(imp$level_fixups) && imp$level_fixups >= 0)

cat("test_index_arms.R: all contracts hold\n")
