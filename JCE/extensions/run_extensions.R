# Runner for the post hoc extensions (extensions.R). Writes only to --out.
#
#   Rscript JCE/extensions/run_extensions.R --verify
#   Rscript JCE/extensions/run_extensions.R --reps 1000 --studies sweep,devsize,devcontact \
#     --workers 6 --out outputs/methodology/extensions_v1
#   Rscript JCE/extensions/run_extensions.R --reps 1000 --arms rimca \
#     --lib ~/R/missMDA-lib --workers 6 --out outputs/methodology/extensions_v1
#
# The second command scores regularized iterative MCA alone, in processes that
# put the library holding missMDA and FactoMineR 2.17 first on the path; every
# other arm runs under the FactoMineR 2.13 that production used.
#
# --verify regenerates replicates 1 and 2 of the reference configuration
# (loading .25, development at 50 percent contact) under both filters and
# requires every production arm to reproduce its saved production metrics
# exactly. Results are saved in chunks of 50 replicates per configuration; a
# chunk already written is skipped, so an interrupted run resumes.

suppressPackageStartupMessages({library(dplyr); library(tibble); library(digest); library(parallel)})

argv <- commandArgs(trailingOnly = TRUE)
get1 <- function(flag, default) {
  i <- match(flag, argv)
  if (is.na(i)) default else argv[i + 1L]
}
verify_only <- "--verify" %in% argv
n_rep <- as.integer(get1("--reps", "1000"))
workers <- as.integer(get1("--workers", "6"))
studies <- strsplit(get1("--studies", "sweep,devsize,devcontact"), ",")[[1]]
only <- get1("--configs", NA_character_)
arm_set <- get1("--arms", "main")
lib <- get1("--lib", NA_character_)
out_dir <- get1("--out", "outputs/methodology/extensions_v1")
chunk <- 50L
root <- normalizePath(".")

setup <- function(root, lib = NA_character_) {
  setwd(root)
  if (!is.na(lib)) .libPaths(c(path.expand(lib), .libPaths()))
  suppressPackageStartupMessages({library(dplyr); library(tibble); library(digest); library(FactoMineR)})
  source("methodology/v2/harmonize_donors.R")
  source("methodology/v2/plasmode_generate.R")
  source("methodology/v2/index_arms.R")
  source("methodology/v2/quadrant_populations.R")
  source("methodology/v2/pilot_v2_config.R")
  source("JCE/theory/filter_capture.R")
  source("JCE/extensions/extensions.R")
  harmonized <- lapply(readRDS("outputs/methodology/donors_harmonized_geo.rds"), prefix_levels)
  keep_cols <- function(d) droplevels(as.data.frame(d[, c(INDEX_INPUTS, "age", "sex", "cluster")]))
  pools <- lapply(build_quadrants(harmonized, quadrant_cells("methodology/v2/quadrant_cells.csv")), keep_cols)
  assign("pools", pools["Q3"], envir = globalenv())
  assign("scales", list(Q3 = signal_scale(pools$Q3, pilot_spec(population = "Q3"))), envir = globalenv())
  invisible(TRUE)
}
setup(root, lib)
arms <- if (arm_set == "main") NULL else strsplit(arm_set, ",")[[1]]
tag <- if (arm_set == "main") "" else paste0("_", gsub(",", "-", arm_set))

grid <- extension_grid(studies)
if (!is.na(only)) grid <- grid[grid$config_id %in% strsplit(only, ",")[[1]], ]
manifest <- readRDS("outputs/methodology/production_v03/manifest.rds")

if (verify_only) {
  ref <- extension_grid("sweep")
  ref <- ref[ref$loading == 0.25, ]
  chk <- bind_rows(lapply(seq_len(nrow(ref)), function(i) {
    cfg <- ref[i, ]
    prod <- readRDS(sprintf("outputs/methodology/production_v03/Q3_nd1000_nt20000_r60_q00_g%s_base.rds",
                            substr(cfg$gate, 1, 3)))
    bind_rows(lapply(1:2, function(r) {
      got <- run_ext_replicate(cfg, r, arms = ALL_ARMS)$metrics
      want <- prod[prod$rep == r, ]
      j <- inner_join(got, want, by = c("arm", "setting"), suffix = c("", ".prod"))
      tibble(gate = cfg$gate, rep = r, rows = nrow(j),
             rho = max(abs(j$cluster_rho_all - j$cluster_rho_all.prod), na.rm = TRUE),
             misalloc = max(abs(j$cl10_misalloc - j$cl10_misalloc.prod), na.rm = TRUE),
             capture = max(abs(j$b10_hn_capture - j$b10_hn_capture.prod), na.rm = TRUE))
    }))
  }))
  print(as.data.frame(chk), row.names = FALSE)
  if (any(chk$rows != 30L) || any(unlist(chk[, c("rho", "misalloc", "capture")]) > 1e-12)) {
    stop("The extension runner does not reproduce production.")
  }
  cat("The reference configuration reproduces production exactly for every production arm.\n")
  quit(save = "no")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cl <- makeCluster(workers, type = "PSOCK")
invisible(clusterCall(cl, setup, root, lib))
clusterExport(cl, c("grid", "arms"))
for (i in seq_len(nrow(grid))) {
  cfg <- grid[i, ]
  for (start in seq(1L, n_rep, by = chunk)) {
    reps <- start:min(n_rep, start + chunk - 1L)
    f <- file.path(out_dir, sprintf("%s%s_reps%04d-%04d.rds", cfg$config_id, tag, min(reps), max(reps)))
    if (file.exists(f)) next
    t0 <- Sys.time()
    res <- parLapply(cl, reps, \(r, cid) tryCatch({
                       cfg_r <- grid[grid$config_id == cid, ]
                       run_ext_replicate(cfg_r, r, if (is.null(arms)) ext_arms(cfg_r) else arms)
                     },
                                                  error = \(e) list(error = conditionMessage(e), rep = r)),
                     cid = cfg$config_id)
    bad <- vapply(res, \(x) !is.null(x$error), logical(1))
    saveRDS(list(config_id = cfg$config_id, reps = reps, arms = if (is.null(arms)) ext_arms(cfg) else arms,
                 manifest_signature = manifest$signature,
                 factominer = as.character(packageVersion("FactoMineR")),
                 metrics = bind_rows(lapply(res[!bad], `[[`, "metrics")),
                 diag_mca = bind_rows(lapply(res[!bad], `[[`, "diag_mca")),
                 diag_pca = bind_rows(lapply(res[!bad], `[[`, "diag_pca")),
                 size = bind_rows(lapply(res[!bad], `[[`, "size")),
                 errors = res[bad]), f)
    cat(sprintf("%s%s reps %d-%d: %d errors, %.1f min\n", cfg$config_id, tag, min(reps), max(reps), sum(bad),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
stopCluster(cl)
