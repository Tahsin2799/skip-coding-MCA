# Contract tests for the post hoc extensions. Run from the repository root:
#   Rscript JCE/extensions/test_extensions.R
# The runner's --verify mode separately checks that the reference configuration
# reproduces production exactly.

suppressPackageStartupMessages({library(dplyr); library(tibble); library(digest); library(FactoMineR); library(testthat)})
for (f in c("harmonize_donors", "plasmode_generate", "index_arms", "quadrant_populations", "pilot_v2_config")) {
  source(file.path("methodology/v2", paste0(f, ".R")))
}
source("JCE/theory/filter_capture.R")
source("JCE/extensions/extensions.R")
harmonized <- lapply(readRDS("outputs/methodology/donors_harmonized_geo.rds"), prefix_levels)
keep_cols <- function(d) droplevels(as.data.frame(d[, c(INDEX_INPUTS, "age", "sex", "cluster")]))
pools <- lapply(build_quadrants(harmonized, quadrant_cells("methodology/v2/quadrant_cells.csv")), keep_cols)["Q3"]
scales <- list(Q3 = signal_scale(pools$Q3, pilot_spec(population = "Q3")))
grid <- extension_grid()

sim_for <- function(cfg, rep = 1L) {
  set.seed(digest2int(paste0(paste(cfg$stream, cfg$stream_id, rep, sep = "|"), "|gen")))
  generate_extension(pools$Q3, cfg$spec[[1]], cfg$n_dev, cfg$n_target, cfg$dev_setting, scales$Q3)
}

test_that("the generator with development at 50 percent contact is the production generator", {
  cfg <- grid[grid$study == "sweep" & grid$loading == 0.25 & grid$gate == "diagnosis", ]
  set.seed(digest2int(paste0("production|", REFERENCE_STREAM, "|3|gen")))
  a <- generate_replicate(pools$Q3, cfg$spec[[1]], 1000L, 20000L, scales$Q3)
  expect_identical(sim_for(cfg, 3L), a)
})

test_that("configurations sharing the production stream see the same people", {
  a <- sim_for(grid[grid$study == "sweep" & grid$loading == 0 & grid$gate == "contact", ])
  b <- sim_for(grid[grid$study == "sweep" & grid$loading == 1 & grid$gate == "diagnosis", ])
  c <- sim_for(grid[grid$study == "devcontact" & grid$dev_setting == "c80" & grid$gate == "contact", ])
  expect_identical(a$targets$c20$V, b$targets$c20$V)
  expect_identical(a$development$V, c$development$V)
  expect_gt(mean(c$development$contacted), mean(a$development$contacted))
})

test_that("FLIC item models are logistf::flic fits and match the event rate", {
  cfg <- grid[grid$study == "sweep" & grid$loading == 0.25 & grid$gate == "diagnosis", ]
  dev <- sim_for(cfg)$development
  mod <- item_models_flic(dev)
  for (v in mod$items) {
    f <- mod$fits[[v]]
    expect_equal(mean(stats::predict(f, newdata = mod$obs, type = "response")), mean(mod$obs[[v]]),
                 tolerance = 1e-6)
    expect_true(all(is.finite(stats::coef(f)[!is.na(stats::coef(f))])))
  }
})

test_that("principal component filter loadings match their closed forms", {
  for (gate in GATES) {
    cfg <- grid[grid$study == "sweep" & grid$loading == 0.25 & grid$gate == gate, ]
    d <- pca_diagnostics(sim_for(cfg)$development)
    expect_lt(max(d$eta2_max_diff), 1e-10)
    expect_lt(d$r2_eligibility[d$coding == "fuzzy"], 0.05)
  }
})

test_that("every extension arm scores a replicate", {
  cfg <- grid[grid$study == "sweep" & grid$loading == 0.5 & grid$gate == "contact", ]
  r <- run_ext_replicate(cfg, 1L, arms = setdiff(EXT_ARMS, c("raw_items", "impute_frozen", "impute_adaptive")))
  expect_false(any(r$metrics$failed))
  expect_true(all(is.finite(r$metrics$cluster_rho_all)))
})
cat("extension tests passed\n")
