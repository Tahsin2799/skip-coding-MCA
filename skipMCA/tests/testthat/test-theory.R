dev <- simulate_filtered_survey(n_clusters = 80, cluster_size = 20, contact = 0.5, seed = 11)
items <- paste0("advice", 1:4)
inputs <- paste0("x", 1:6)
p <- mean(dev$eligible)

test_that("the empirical squared correlation ratio equals Theorem 1 for every coding", {
  a <- skip_audit(dev, items, inputs)
  expect_equal(a$items$eta2, a$items$formula, tolerance = 1e-10)
  by_coding <- split(a$items, a$items$coding)
  expect_equal(by_coding$category$eta2, rep(1, 4), tolerance = 1e-10)
  expect_equal(by_coding$subset$eta2, rep(1 - p, 4), tolerance = 1e-10)
  expect_lt(max(abs(by_coding$fuzzy_average$eta2)), 1e-12)
  yes <- vapply(items, function(v) mean(dev[[v]][dev$eligible] == "yes"), numeric(1))
  expect_equal(by_coding$no$eta2, unname(yes * (1 - p) / (1 - p * yes)), tolerance = 1e-10)
})

test_that("eligibility_eta2 reproduces the corollary values", {
  pi <- 0.3
  expect_equal(eligibility_eta2(0.4, c(1 - pi, pi, 0), c(0, 0, 1)), 1)
  expect_equal(eligibility_eta2(0.4, c(1 - pi, pi, 0), c(0, 0, 1), c(TRUE, TRUE, FALSE)), 0.6)
  expect_equal(eligibility_eta2(0.4, c(1 - pi, pi), c(1, 0)), pi * 0.6 / (1 - 0.4 * pi))
  expect_equal(eligibility_eta2(0.4, c(1 - pi, pi), c(1 - pi, pi)), 0)
  expect_error(eligibility_eta2(0.4, c(0.5, 0.5), c(1, 0, 0)))
})

test_that("under missing fuzzy average the filtered block's eigenvalue is at most the share asked", {
  a <- skip_audit(dev, items, inputs, codings = "fuzzy_average")
  expect_lte(a$codings$lambda_items, p + 1e-12)
  expect_lt(a$codings$r2_eligibility, 0.05)
})

test_that("missing single gives a first eigenvalue of at least k/Q", {
  expect_gte(skip_mca(dev, items, inputs, coding = "category")$eig[1], 4 / 10 - 1e-12)
})

test_that("prediction reproduces the development scores and ignores row order", {
  for (cd in c("fuzzy_average", "subset", "category", "no", "drop")) {
    fit <- skip_mca(dev, items, inputs, coding = cd, orient = dev$need)
    expect_equal(unname(predict(fit, dev)), unname(fit$scores), tolerance = 1e-10)
    perm <- sample(nrow(dev))
    expect_equal(unname(predict(fit, dev[perm, ])), unname(fit$scores[perm, , drop = FALSE]),
                 tolerance = 1e-10)
  }
})

test_that("target rows with a category development never saw are left unscored", {
  tgt <- simulate_filtered_survey(n_clusters = 20, cluster_size = 20, contact = 0.2, seed = 12)
  tgt$x1 <- factor(as.character(tgt$x1), levels = c("low", "mid", "high", "new"))
  tgt$x1[1] <- "new"
  s <- predict(skip_mca(dev, items, inputs, orient = dev$need), tgt)
  expect_true(is.na(s[1, 1]))
  expect_false(anyNA(s[-1, 1]))
})

test_that("inconsistent skip patterns are refused", {
  bad <- dev
  bad$advice2[which(dev$eligible)[1]] <- NA
  expect_error(skip_mca(bad, items, inputs), "eligible respondents have missing answers")
  bad <- dev
  bad$advice2[which(!dev$eligible)[1]] <- "no"
  expect_error(skip_mca(bad, items, inputs), "share one filter")
})

test_that("a perfect score targets perfectly", {
  m <- targeting_metrics(dev$need, dev$need, dev$cluster)
  expect_equal(m[["cluster_rank_correlation"]], 1)
  expect_equal(m[["inclusion_error"]], 0)
  expect_equal(m[["top_decile_precision"]], 1)
  expect_equal(m[["high_need_selected"]], 0.5, tolerance = 0.01)
})

test_that("the codings agree with FactoMineR", {
  skip_if_not_installed("FactoMineR")
  fz <- skip_mca(dev, items, inputs, coding = "fuzzy_average")
  zfa <- skipMCA:::index_block(dev, fz$spec, fz$eligible, "fuzzy_average")$z
  ca <- FactoMineR::CA(as.data.frame(zfa), ncp = 1, graph = FALSE)
  expect_equal(fz$eig[1], ca$eig[1, 1], tolerance = 1e-8)
  expect_equal(abs(fz$scores[, 1]), abs(as.matrix(ca$row$coord)[, 1]), tolerance = 1e-8,
               ignore_attr = TRUE)
  d <- dev[, c(inputs, items)]
  for (v in items) d[[v]] <- factor(ifelse(is.na(d[[v]]), "skipped", as.character(d[[v]])))
  own <- skip_mca(dev, items, inputs, coding = "category")
  mca <- FactoMineR::MCA(d, ncp = 1, graph = FALSE)
  expect_equal(own$eig[1], mca$eig[1, 1], tolerance = 1e-8)
  expect_equal(abs(own$scores[, 1]), abs(as.matrix(mca$ind$coord)[, 1]), tolerance = 1e-8,
               ignore_attr = TRUE)
  sub <- skip_mca(dev, items, inputs, coding = "subset")
  skipcols <- which(unlist(lapply(d, levels)) == "skipped")
  smca <- FactoMineR::MCA(d, ncp = 1, excl = skipcols, graph = FALSE)
  expect_equal(sub$eig[1], smca$eig[1, 1], tolerance = 1e-8)
  expect_equal(abs(sub$scores[, 1]), abs(as.matrix(smca$ind$coord)[, 1]), tolerance = 1e-8,
               ignore_attr = TRUE)
})
