# Theory diagnostics and the missing fuzzy average on the exact production
# replicates, regenerated from their production seeds.
#
# Post hoc. The coding was derived from the theory after production, so its
# results are reported apart from the 132 frozen verdicts, under the same rule
# and margins, and labelled as post hoc. Nothing in methodology/v2 or in
# outputs/methodology/production_v03 is written.
#
#   Rscript JCE/theory/theory_production_link.R --verify
#   Rscript JCE/theory/theory_production_link.R --reps 200 --blocks core --workers 6 \
#     --out outputs/methodology/theory_link_v1
#
# One file per configuration in --out; a configuration already written is
# skipped, so an interrupted run resumes.

suppressPackageStartupMessages({library(dplyr); library(tibble); library(digest); library(parallel)})

argv <- commandArgs(trailingOnly = TRUE)
get1 <- function(flag, default) {
  i <- match(flag, argv)
  if (is.na(i)) default else argv[i + 1L]
}
verify_only <- "--verify" %in% argv
n_rep <- as.integer(get1("--reps", "200"))
workers <- as.integer(get1("--workers", "6"))
blocks <- strsplit(get1("--blocks", "core"), ",")[[1]]
out_dir <- get1("--out", "outputs/methodology/theory_link_v1")
root <- normalizePath(".")

setup <- function(root) {
  setwd(root)
  suppressPackageStartupMessages({library(dplyr); library(tibble); library(digest); library(FactoMineR)})
  source("methodology/v2/harmonize_donors.R")
  source("methodology/v2/plasmode_generate.R")
  source("methodology/v2/index_arms.R")
  source("methodology/v2/quadrant_populations.R")
  source("methodology/v2/pilot_v2_config.R")
  source("JCE/theory/filter_capture.R")
  # Pools exactly as run_pilot_v2.R builds them.
  harmonized <- lapply(readRDS("outputs/methodology/donors_harmonized_geo.rds"), prefix_levels)
  keep_cols <- function(d) droplevels(as.data.frame(d[, c(INDEX_INPUTS, "age", "sex", "cluster")]))
  pools <- lapply(build_quadrants(harmonized, quadrant_cells("methodology/v2/quadrant_cells.csv")), keep_cols)
  pools[["bd_allage"]] <- keep_cols(build_age_pool(harmonized, "bangladesh", c(18L, Inf)))
  for (nm in names(DONOR_POPULATIONS)) {
    pools[[nm]] <- keep_cols(build_age_pool(harmonized, DONOR_POPULATIONS[[nm]], c(18L, 49L)))
  }
  assign("pools", pools, envir = globalenv())
  invisible(TRUE)
}
setup(root)

grid <- v03_grid()
grid <- grid[grid$block %in% blocks, ]
manifest <- readRDS("outputs/methodology/production_v03/manifest.rds")

regenerate <- function(cfg, rep) {
  spec <- cfg$spec[[1]]
  set.seed(digest2int(paste0(paste("production", cfg$stream_id, rep, sep = "|"), "|gen")))
  generate_replicate(pools[[cfg$population]], spec, cfg$n_dev, cfg$n_target,
                     signal_scale(pools[[cfg$population]], spec))
}

link_one <- function(cfg, rep) {
  spec <- cfg$spec[[1]]
  sim <- regenerate(cfg, rep)
  diag <- if (spec$item_mode == "advice") coding_diagnostics(sim$development) |>
    mutate(config_id = cfg$config_id, rep = rep)
  prep <- fuzzy_average_prep(sim$development)
  met <- bind_rows(lapply(spec$settings, \(s)
    score_prepared(prep, "fuzzy_average", sim$targets[[s]], spec) |> mutate(setting = s))) |>
    mutate(config_id = cfg$config_id, rep = rep, arm = "fuzzy_average", failed = FALSE)
  list(diag = diag, metrics = met)
}

# Regeneration check: the true V of the saved artifact rows and one production
# metric must be reproduced exactly before any result is trusted.
verify <- function(cfg, rep = 1L) {
  sim <- regenerate(cfg, rep)
  art <- readRDS(sprintf("outputs/methodology/production_v03/artifacts/%s_rep%04d.rds", cfg$config_id, rep))
  t0 <- sim$targets[[cfg$spec[[1]]$settings[1]]]
  prod <- readRDS(sprintf("outputs/methodology/production_v03/%s.rds", cfg$config_id))
  want <- prod$cluster_rho_all[prod$arm == "drop_items" & prod$rep == rep & prod$setting == "c20"]
  got <- score_prepared(prepare_family(sim$development, "drop_items", cfg$spec[[1]]), "drop_items",
                        sim$targets$c20, cfg$spec[[1]])$cluster_rho_all
  tibble(config_id = cfg$config_id, rep = rep,
         max_abs_V_diff = max(abs(t0$V[art$truth$row_id] - art$truth$V)),
         drop_items_rho_diff = abs(got - want))
}

if (verify_only) {
  chk <- bind_rows(lapply(seq_len(nrow(grid)), \(i) verify(grid[i, ])))
  print(as.data.frame(chk), row.names = FALSE)
  if (any(chk$max_abs_V_diff > 0 | chk$drop_items_rho_diff > 1e-12)) stop("Regeneration does not reproduce production.")
  cat("Regeneration reproduces production in all", nrow(chk), "configurations checked.\n")
  quit(save = "no")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cl <- makeCluster(workers, type = "PSOCK")
invisible(clusterCall(cl, setup, root))
clusterExport(cl, c("regenerate", "link_one", "grid"))
for (i in seq_len(nrow(grid))) {
  cfg <- grid[i, ]
  f <- file.path(out_dir, paste0(cfg$config_id, ".rds"))
  if (file.exists(f)) next
  t0 <- Sys.time()
  res <- parLapply(cl, seq_len(n_rep), \(r, cid) tryCatch(link_one(grid[grid$config_id == cid, ], r),
                                                         error = \(e) list(error = conditionMessage(e), rep = r)),
                   cid = cfg$config_id)
  bad <- vapply(res, \(x) !is.null(x$error), logical(1))
  saveRDS(list(config_id = cfg$config_id, reps = n_rep, manifest_signature = manifest$signature,
               diag = bind_rows(lapply(res[!bad], `[[`, "diag")),
               metrics = bind_rows(lapply(res[!bad], `[[`, "metrics")),
               errors = res[bad]), f)
  cat(sprintf("%s: %d replicates, %d errors, %.1f min\n", cfg$config_id, n_rep, sum(bad),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
stopCluster(cl)
