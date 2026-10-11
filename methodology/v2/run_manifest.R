# Provenance for a v03 run directory: what was run, from which sources, under
# which package versions, and the rules for resuming one safely.
#
# The pattern is v1's (methodology/run_pilot.R), kept deliberately: hash every
# scientific source, sign the manifest, refuse a resume that does not match, and
# let readers take their file list from the manifest rather than from a glob.
#
#   source("methodology/v2/run_manifest.R")

suppressPackageStartupMessages({library(digest); library(cli)})

MANIFEST_SCHEMA <- "v03-plasmode"

HASHED_SOURCES <- c("harmonize_donors.R", "plasmode_generate.R", "index_arms.R",
                    "quadrant_populations.R", "pilot_v2_config.R", "run_pilot_v2.R",
                    "run_manifest.R", "solve_intercepts.R", "summarize_convergence.R",
                    "protocol.md", "PROTOCOL_AMENDMENT_2026-09-20.md",
                    "METRIC_SPECIFICATION.md")

MANIFEST_PACKAGES <- c("FactoMineR", "dplyr", "tidyr", "tibble", "digest", "cli", "haven")

# Everything that shapes one configuration's results, including the factors held
# at a single value across a run. A spec-only hash would miss the sample sizes
# and the arm list.
config_hash <- function(cfg, spec, arms) {
  digest(list(config_id = cfg$config_id, stream_id = cfg$stream_id,
              population = cfg$population, n_dev = cfg$n_dev, n_target = cfg$n_target,
              gate = spec$gate, variant = cfg$variant, spec = spec, arms = arms),
         algo = "sha256")
}

run_manifest <- function(grid, arms, stream, root) {
  src <- file.path(root, "methodology/v2", HASHED_SOURCES)
  present <- file.exists(src)
  body <- list(
    schema = MANIFEST_SCHEMA,
    stream = stream,
    grid = grid[, setdiff(names(grid), "spec")],
    config_hashes = vapply(seq_len(nrow(grid)),
                           \(i) config_hash(grid[i, ], grid$spec[[i]], arms), character(1)),
    arms = arms,
    metric_spec = "METRIC_SPECIFICATION.md v03-2026-09-20",
    source_hashes = setNames(vapply(src[present], \(p) digest(file = p, algo = "sha256"),
                                    character(1)), HASHED_SOURCES[present]),
    missing_sources = basename(src[!present]),
    donor_cache = digest(file = file.path(root, "outputs/methodology/donors_harmonized_geo.rds"),
                         algo = "sha256"),
    quadrant_cells = digest(file = file.path(root, "methodology/v2/quadrant_cells.csv"),
                            algo = "sha256"),
    packages = vapply(MANIFEST_PACKAGES, \(p) as.character(utils::packageVersion(p)), character(1)),
    r_version = R.version.string,
    platform = R.version$platform,
    rng = paste("Mersenne-Twister; set.seed(digest2int(stream|stream_id|rep|purpose));",
                "purposes: gen, prep|family, score|arm|setting"))
  c(body, list(signature = digest(body, algo = "sha256"),
               created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)))
}

# Names the component that moved rather than only refusing. A source edit, a
# package upgrade, a changed donor cache and a changed grid are different
# problems with different answers.
manifest_diff <- function(previous, current) {
  fields <- c("schema", "stream", "arms", "metric_spec", "source_hashes", "donor_cache",
              "quadrant_cells", "packages", "r_version", "platform", "rng", "config_hashes")
  moved <- fields[!vapply(fields, \(f) identical(previous[[f]], current[[f]]), logical(1))]
  if ("source_hashes" %in% moved) {
    a <- previous$source_hashes
    b <- current$source_hashes
    files <- union(names(a), names(b))
    moved <- c(moved, paste0("  source: ", basename(files[!vapply(files,
      \(f) identical(a[[f]], b[[f]]), logical(1))])))
  }
  if ("config_hashes" %in% moved && nrow(previous$grid) == nrow(current$grid)) {
    i <- which(previous$config_hashes != current$config_hashes)
    moved <- c(moved, paste0("  configuration: ", current$grid$config_id[i]))
  }
  moved
}

initialize_run <- function(output_dir, manifest) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (length(manifest$missing_sources)) {
    cli_abort("Missing hashed source{?s}: {.file {manifest$missing_sources}}.")
  }
  path <- file.path(output_dir, "manifest.rds")
  if (file.exists(path)) {
    previous <- readRDS(path)
    if (!identical(previous$signature, manifest$signature)) {
      cli_abort(c("{.file {output_dir}} was written under a different manifest.",
                  i = "Moved: {.val {manifest_diff(previous, manifest)}}",
                  i = "Use a new output directory, or restore the recorded inputs."))
    }
  } else saveRDS(manifest, path, version = 3)
  invisible(path)
}

# Who holds a lock, or NULL if the record is unusable. Written by claim_lock().
lock_owner <- function(lock) {
  f <- file.path(lock, "owner")
  if (!file.exists(f)) return(NULL)
  ln <- readLines(f, warn = FALSE)
  pid <- suppressWarnings(as.integer(sub("^pid ", "", grep("^pid ", ln, value = TRUE)[1])))
  host <- sub("^host ", "", grep("^host ", ln, value = TRUE)[1])
  if (is.na(pid) || is.na(host)) return(NULL)
  list(pid = pid, host = host)
}

# Alive *and* an R process. A pid alone is not enough: pids are recycled, and a
# recycled one held by something unrelated would keep a dead run's lock forever.
owner_alive <- function(pid) {
  out <- suppressWarnings(system2("ps", c("-p", pid, "-o", "command="),
                                  stdout = TRUE, stderr = FALSE))
  length(out) > 0L && any(grepl("R\\.framework|Rscript|run_pilot_v2", out))
}

# A directory lock, not a file lock: one writer per run directory. A lock whose
# owner is provably gone is reclaimed, because a run that is killed or
# interrupted cannot release it (a signal does not run `on.exit`), and requiring
# a manual delete before every resume is how a resumable run stops being one.
# A lock whose owner is alive, or that cannot be checked, is still refused: two
# writers racing on one configuration file is the failure this exists to stop.
claim_lock <- function(output_dir) {
  lock <- file.path(output_dir, "run.lock")
  if (!dir.create(lock, showWarnings = FALSE)) {
    own <- lock_owner(lock)
    if (is.null(own)) {
      cli_abort(c("{.file {lock}} exists and records no usable owner.",
                  i = "Remove it by hand once no run is in progress: {.code rm -r {lock}}"))
    }
    if (!identical(own$host, Sys.info()[["nodename"]])) {
      cli_abort(c("{.file {lock}} is held by pid {own$pid} on {.val {own$host}}.",
                  i = "This host cannot tell whether that run is alive.",
                  i = "Remove it by hand once no run is in progress: {.code rm -r {lock}}"))
    }
    if (owner_alive(own$pid)) {
      cli_abort(c("{.file {lock}} is held by pid {own$pid}, which is running.",
                  i = "One writer per run directory. Stop it, or wait for it."))
    }
    # Rename before deleting, so that two processes reclaiming at once have
    # exactly one winner: the loser's rename fails and it aborts below.
    stale <- paste0(lock, ".stale-", Sys.getpid())
    if (!file.rename(lock, stale) || !dir.create(lock, showWarnings = FALSE)) {
      cli_abort(c("{.file {lock}} could not be reclaimed.",
                  i = "Another run may have taken it first."))
    }
    unlink(stale, recursive = TRUE)
    cli_alert_info("Reclaimed {.file {lock}}: pid {own$pid} is gone.")
  }
  writeLines(c(paste("pid", Sys.getpid()), paste("host", Sys.info()[["nodename"]]),
               format(Sys.time(), tz = "UTC", usetz = TRUE)),
             file.path(lock, "owner"))
  lock
}

release_lock <- function(lock) unlink(lock, recursive = TRUE)

# Results are read back only through their own signature and configuration hash,
# so a file written under a different specification is refused rather than
# appended to.
read_verified <- function(path, signature, hash) {
  if (!file.exists(path)) return(NULL)
  x <- readRDS(path)
  p <- attr(x, "provenance")
  if (!identical(p$signature, signature) || !identical(p$config_hash, hash)) {
    cli_abort(c("{.file {basename(path)}} was written under a different specification.",
                i = "Give this run its own {.arg --out} directory, or delete the file if it is stale."))
  }
  x
}

save_atomic <- function(x, path) {
  tmp <- tempfile("pending-", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp, version = 3)
  if (!file.rename(tmp, path)) cli_abort("Atomic rename failed for {.file {path}}.")
  invisible(path)
}

# The file list a reader may use. A stray rds in the directory is not a result.
manifest_files <- function(output_dir) {
  m <- readRDS(file.path(output_dir, "manifest.rds"))
  paths <- file.path(output_dir, paste0(m$grid$config_id, ".rds"))
  list(manifest = m, paths = paths[file.exists(paths)],
       missing = m$grid$config_id[!file.exists(paths)])
}
