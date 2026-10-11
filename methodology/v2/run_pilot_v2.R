# Runner for the v03 plasmode simulation. Resumable, one file per configuration,
# one manifest per directory.
#
# Never mclapply: forked workers segfault inside La.svd when MCA runs in them.
# PSOCK only. Export VECLIB_MAXIMUM_THREADS=1 before R starts, not inside it.
#
#   Rscript methodology/v2/run_pilot_v2.R --replicates 100 --workers 6 \
#     --stream pilot --blocks core,stress,controls,quadrant,donor \
#     --out outputs/methodology/revision_pilot_v03
#
# Every factor is in the configuration id, so one directory can hold the whole
# grid. The manifest lists it; readers take their file list from the manifest
# rather than globbing, and a resume is refused when any scientific input moved.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(parallel)
  library(digest)
  library(cli)
})

root <- c(".", "..", "../..") |>
  Filter(f = \(p) file.exists(file.path(p, "methodology/v2/plasmode_generate.R"))) |>
  (\(x) x[1])()

source(file.path(root, "methodology/v2/harmonize_donors.R"))
source(file.path(root, "methodology/v2/plasmode_generate.R"))
source(file.path(root, "methodology/v2/index_arms.R"))
source(file.path(root, "methodology/v2/quadrant_populations.R"))
source(file.path(root, "methodology/v2/pilot_v2_config.R"))
source(file.path(root, "methodology/v2/run_manifest.R"))

# ---- Arguments ----------------------------------------------------------

read_args <- function(argv) {
  get1 <- function(flag, default) {
    i <- match(flag, argv)
    if (is.na(i)) return(default)
    if (i == length(argv)) cli_abort("{.arg {flag}} was given with no value.")
    argv[i + 1L]
  }
  n_rep <- suppressWarnings(as.integer(get1("--replicates", "100")))
  workers <- suppressWarnings(as.integer(get1("--workers", "4")))
  if (is.na(n_rep) || n_rep < 1L) {
    cli_abort("{.arg --replicates} must be a positive integer, not {.val {get1('--replicates', NA)}}.")
  }
  if (is.na(workers) || workers < 1L || workers > 8L) {
    cli_abort("{.arg --workers} must be between 1 and 8, not {.val {get1('--workers', NA)}}.")
  }
  stream <- get1("--stream", "pilot")
  if (!stream %in% c("pilot", "production", "smoke")) {
    cli_abort(c("{.arg --stream} must be {.val pilot}, {.val production} or {.val smoke}.",
                i = "The stream is part of every seed, so the three never share replicates."))
  }
  blocks <- strsplit(get1("--blocks", "core,stress,controls,quadrant,donor"), ",")[[1]]
  unknown <- setdiff(blocks, c("core", "stress", "controls", "quadrant", "donor"))
  if (length(unknown)) cli_abort("No such block: {.val {unknown}}.")
  configs <- get1("--configs", NA_character_)
  arts <- suppressWarnings(as.integer(get1("--artifact-reps", "0")))
  list(n_rep = n_rep, workers = workers, stream = stream, blocks = blocks,
       configs = if (is.na(configs)) NULL else strsplit(configs, ",")[[1]],
       artifact_reps = if (is.na(arts)) 0L else arts,
       # A scratch default on purpose: every run meant to be kept names its own
       # directory, so a probe cannot land on top of one that is.
       out = get1("--out", file.path(root, "outputs/methodology/probe_v03")))
}

args <- read_args(commandArgs(trailingOnly = TRUE))

# The manifest always declares the whole grid, whatever this launch runs: a
# directory holds one grid, several commands fill it, and the verification can
# tell a configuration that has not run yet from one that was never declared.
grid <- v03_grid()
run_grid <- grid[grid$block %in% args$blocks, ]
if (!is.null(args$configs)) {
  miss <- setdiff(args$configs, grid$config_id)
  if (length(miss)) cli_abort("No configuration named {.val {miss}}.")
  run_grid <- run_grid[run_grid$config_id %in% args$configs, ]
}
if (!nrow(run_grid)) cli_abort("The selected blocks and configurations are empty.")

unsolved <- vapply(grid$spec, \(s) anyNA(unlist(s[c("outcome_intercept", "contact_intercepts")])),
                   logical(1))
if (any(unsolved)) {
  cli_abort(c("{sum(unsolved)} configuration{?s} carr{?ies/y} an unsolved intercept.",
              i = "Run {.code Rscript methodology/v2/solve_intercepts.R --reps 20 --stress} and paste its output into pilot_v2_config.R."))
}

# ---- Pools --------------------------------------------------------------

harmonized <- lapply(readRDS(file.path(root, "outputs/methodology/donors_harmonized_geo.rds")),
                     prefix_levels)
keep_cols <- function(d) droplevels(as.data.frame(
  d[, c(INDEX_INPUTS, "age", "sex", "cluster")]))

needed <- unique(run_grid$population)
pools <- list()
if (any(needed %in% QUADRANTS)) {
  pools <- lapply(build_quadrants(harmonized,
    quadrant_cells(file.path(root, "methodology/v2/quadrant_cells.csv"))), keep_cols)
}
if ("bd_allage" %in% needed) {
  pools[["bd_allage"]] <- keep_cols(build_age_pool(harmonized, "bangladesh", c(18L, Inf)))
}
for (nm in intersect(needed, names(DONOR_POPULATIONS))) {
  pools[[nm]] <- keep_cols(build_age_pool(harmonized, DONOR_POPULATIONS[[nm]], c(18L, 49L)))
}

# The signal scale is a property of the pool and the specification, not of a
# draw, so V means the same thing in the development and target halves.
scales <- lapply(seq_len(nrow(run_grid)), \(i) signal_scale(pools[[run_grid$population[i]]],
                                                            run_grid$spec[[i]]))
names(scales) <- run_grid$config_id

# ---- One replicate ------------------------------------------------------

# Independent substreams: the population, each arm family's development half and
# each arm-setting scoring step draw from their own seeds. Adding an arm cannot
# move the population or another arm's draws, and the stream id leaves the gate
# out, so the two gates of a cell see the same people and the same uniforms.
run_replicate <- function(rep_id, cfg, spec, scl, stream, artifact_dir = NULL,
                          artifact_reps = 0L) {
  if (rep_id > artifact_reps) artifact_dir <- NULL
  key <- paste(stream, cfg$stream_id, rep_id, sep = "|")
  started <- Sys.time()
  set.seed(digest2int(paste0(key, "|gen")))
  sim <- generate_replicate(pools[[cfg$population]], spec, cfg$n_dev, cfg$n_target, scl)

  preps <- list()
  for (fam in unique(vapply(ALL_ARMS, prep_family, character(1)))) {
    set.seed(digest2int(paste0(key, "|prep|", fam)))
    preps[[fam]] <- tryCatch(prepare_family(sim$development, fam, spec), error = \(e) e)
  }
  rows <- vector("list", length(spec$settings) * length(ALL_ARMS))
  i <- 0L
  arts <- list()
  for (arm in ALL_ARMS) {
    prep <- preps[[prep_family(arm)]]
    for (s in spec$settings) {
      i <- i + 1L
      set.seed(digest2int(paste0(key, "|score|", arm, "|", s)))
      rows[[i]] <- if (inherits(prep, "error")) {
        tibble(setting = s, arm = arm, failed = TRUE, stage = "development",
               message = conditionMessage(prep))
      } else tryCatch(
        score_prepared(prep, arm, sim$targets[[s]], spec) |>
          mutate(setting = s, arm = arm, failed = FALSE, stage = NA_character_,
                 message = NA_character_),
        error = \(e) tibble(setting = s, arm = arm, failed = TRUE, stage = "target",
                            message = conditionMessage(e)))
      # Enough to recompute every metric from outside: the score each arm gives
      # each target person, and the development cuts those scores are read
      # against. The check is only a check if it does not call the same helper.
      if (!is.null(artifact_dir) && !inherits(prep, "error")) {
        set.seed(digest2int(paste0(key, "|score|", arm, "|", s)))
        arts[[paste(arm, s, sep = "_")]] <- list(
          score = prep$project(sim$targets[[s]], arm)$score,
          breaks = prep$dt$breaks, flip = prep$dt$flip, n_tiers = prep$dt$n_tiers)
      }
    }
  }
  rates <- tibble(
    setting = spec$settings,
    contact_rate = vapply(spec$settings, \(s) mean(sim$targets[[s]]$contacted), numeric(1)),
    recognized_rate = vapply(spec$settings, \(s) mean(sim$targets[[s]]$recognized), numeric(1)),
    eligible_rate = vapply(spec$settings, \(s) mean(sim$targets[[s]]$eligible), numeric(1)),
    unrecognized_case_rate = vapply(spec$settings,
      \(s) mean(sim$targets[[s]]$y == 1L & sim$targets[[s]]$recognized == 0L), numeric(1)),
    outcome_rate = vapply(spec$settings, \(s) mean(sim$targets[[s]]$y), numeric(1)))

  if (!is.null(artifact_dir)) {
    t0 <- sim$targets[[spec$settings[1]]]
    saveRDS(list(config_id = cfg$config_id, rep = rep_id, scores = arts,
                 truth = t0[, c("row_id", "cluster", "age", "sex", "V", "y", "tiebreak",
                                "cluster_tiebreak")],
                 access = lapply(sim$targets, \(t) t[, c("contacted", "recognized", "eligible")]),
                 dev_eligible = sum(observed_rows(sim$development))),
            file.path(artifact_dir, sprintf("%s_rep%04d.rds", cfg$config_id, rep_id)))
  }

  bind_rows(rows) |>
    left_join(rates, by = "setting") |>
    mutate(config_id = cfg$config_id, rep = rep_id, population = cfg$population,
           r2 = spec$r2, q1_age_slope = spec$q1_age_slope, gate = spec$gate,
           variant = cfg$variant, block = cfg$block, n_dev = cfg$n_dev,
           elapsed = as.numeric(difftime(Sys.time(), started, units = "secs")))
}

# ---- Drive --------------------------------------------------------------

# Inside a function, not at top level: `on.exit` registered at top level of an
# Rscript never fires, so the lock outlived the run that took it and the worker
# cluster was never stopped.
main <- function() {

manifest <- run_manifest(grid, ALL_ARMS, args$stream, root)
initialize_run(args$out, manifest)
lock <- claim_lock(args$out)
on.exit(release_lock(lock), add = TRUE)

artifact_dir <- if (args$artifact_reps > 0L) {
  d <- file.path(args$out, "artifacts")
  dir.create(d, showWarnings = FALSE)
  d
}

cl <- makeCluster(args$workers, type = "PSOCK")
on.exit(stopCluster(cl), add = TRUE)
clusterExport(cl, c("root", "pools", "run_replicate"), envir = globalenv())
invisible(clusterEvalQ(cl, {
  suppressPackageStartupMessages({library(dplyr); library(tibble); library(digest)})
  source(file.path(root, "methodology/v2/harmonize_donors.R"))
  source(file.path(root, "methodology/v2/plasmode_generate.R"))
  source(file.path(root, "methodology/v2/index_arms.R"))
}))

started <- Sys.time()
done <- 0L
chunk <- max(args$workers * 5L, 25L)

for (k in seq_len(nrow(run_grid))) {
  cfg <- run_grid[k, ]
  spec <- run_grid$spec[[k]]
  path <- file.path(args$out, paste0(cfg$config_id, ".rds"))
  have <- read_verified(path, manifest$signature, config_hash(cfg, spec, ALL_ARMS))
  todo <- setdiff(seq_len(args$n_rep), unique(have$rep))
  if (!length(todo)) next

  out <- have
  for (part in split(todo, ceiling(seq_along(todo) / chunk))) {
    res <- parLapply(cl, part, run_replicate, cfg = cfg, spec = spec,
                     scl = scales[[cfg$config_id]], stream = args$stream,
                     artifact_dir = artifact_dir, artifact_reps = args$artifact_reps)
    out <- bind_rows(out, bind_rows(res))
    attr(out, "provenance") <- list(
      signature = manifest$signature, config_hash = config_hash(cfg, spec, ALL_ARMS),
      config_id = cfg$config_id, stream = args$stream, spec = spec,
      r_version = R.version.string, platform = R.version$platform,
      host = Sys.info()[["nodename"]], finished = Sys.time(), workers = args$workers)
    save_atomic(out, path)
    done <- done + length(part)
  }
  cli_alert_success("{cfg$config_id}: {length(todo)} replicate{?s} written, {length(unique(out$rep))} held")
}

elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
per_rep <- if (done > 0L) elapsed / done else NA_real_
cat(sprintf(
  "\n%s %s\n  configurations %d of %d in the manifest\n  replicates     %d new in %.1f min\n  sec/replicate  %.2f wall on %d worker(s)\n  projected      %.1f h for %d replicates x %d configurations\n",
  manifest$schema, args$stream, nrow(run_grid), nrow(grid), done, elapsed / 60,
  per_rep, args$workers, per_rep * args$n_rep * nrow(grid) / 3600, args$n_rep, nrow(grid)))

}

main()
