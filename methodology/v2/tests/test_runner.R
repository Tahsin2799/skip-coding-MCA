# Runner and provenance contracts. Run: Rscript methodology/v2/tests/test_runner.R
#
# These cost a couple of minutes because they start the runner for real. They
# are worth it: each one guards a way of losing or silently corrupting a run
# that takes most of a day.

suppressPackageStartupMessages({library(dplyr); library(digest)})
root <- c(".", "..", "../..") |>
  Filter(f = \(p) file.exists(file.path(p, "methodology/v2/plasmode_generate.R"))) |>
  (\(x) x[1])()
source(file.path(root, "methodology/v2/harmonize_donors.R"))
source(file.path(root, "methodology/v2/plasmode_generate.R"))
source(file.path(root, "methodology/v2/pilot_v2_config.R"))
source(file.path(root, "methodology/v2/run_manifest.R"))

runner <- file.path(root, "methodology/v2/run_pilot_v2.R")
out <- file.path(tempdir(), paste0("runner_test_", as.integer(Sys.time())))
dir.create(out, showWarnings = FALSE)
on.exit(unlink(out, recursive = TRUE), add = TRUE)

run <- function(...) {
  res <- suppressWarnings(system2("Rscript", c(shQuote(runner), ...),
                                  stdout = TRUE, stderr = TRUE,
                                  env = c("VECLIB_MAXIMUM_THREADS=1", "OMP_NUM_THREADS=1")))
  list(status = if (is.null(attr(res, "status"))) 0L else attr(res, "status"),
       text = paste(res, collapse = "\n"))
}

CFG <- "Q3_nd1000_nt20000_r60_q00_gcon_base"
base <- c("--replicates", "2", "--workers", "2", "--blocks", "core",
          "--configs", CFG, "--stream", "smoke", "--out", shQuote(out))

# ---- The grid is identified, not implied --------------------------------
# A configuration id used to leave out any factor held at one value, so two
# specifications could share a file and a seed stream.
g <- v03_grid()
stopifnot("every configuration id is unique" = !anyDuplicated(g$config_id),
          "the id carries the gate" = all(grepl("_gcon_|_gdia_", g$config_id)),
          "and the population, both sample sizes, r2 and the slope" =
            all(grepl("^[A-Za-z0-9_]+_nd[0-9]+_nt[0-9]+_r[0-9]{2}_q[0-9]{2}_", g$config_id)))
stopifnot("the two gates of a cell share one population stream" =
            sum(duplicated(g$stream_id)) == length(unique(g$stream_id[duplicated(g$stream_id)])))
stopifnot("no configuration carries an unsolved intercept" =
            !any(vapply(g$spec, \(s) anyNA(unlist(s[c("outcome_intercept", "contact_intercepts")])),
                        logical(1))))

# ---- An unknown block or configuration is refused ------------------------
bad <- run(c("--blocks", "nonsense", "--out", shQuote(out)))
stopifnot("an unknown block aborts" = bad$status != 0L && grepl("nonsense", bad$text))
bad2 <- run(c("--blocks", "core", "--configs", "no_such_config", "--out", shQuote(out)))
stopifnot("an unknown configuration aborts" = bad2$status != 0L)
bad3 <- run(c("--replicates", "2", "--workers", "2", "--blocks", "core",
              "--configs", CFG, "--stream", "experimental", "--out", shQuote(out)))
stopifnot("an unknown stream aborts" = bad3$status != 0L && grepl("stream", bad3$text))

# ---- A run writes a manifest and what it says it writes ------------------
first <- run(base)
stopifnot("the run completes" = first$status == 0L)
f <- file.path(out, paste0(CFG, ".rds"))
stopifnot("it writes one file per configuration" = file.exists(f),
          "and a manifest" = file.exists(file.path(out, "manifest.rds")))
m <- readRDS(file.path(out, "manifest.rds"))
stopifnot("the manifest carries the whole declared grid, not just what ran" =
            nrow(m$grid) == nrow(v03_grid()),
          "the hashed sources" = length(m$source_hashes) >= 10L,
          "the donor cache identity" = nchar(m$donor_cache) == 64L,
          "the package versions" = "FactoMineR" %in% names(m$packages),
          "the RNG policy" = grepl("substream|purposes", m$rng),
          "and the stream" = identical(m$stream, "smoke"))

d <- readRDS(f)
stopifnot("two replicates" = length(unique(d$rep)) == 2L,
          "ten arms at three settings" = nrow(d) == 2L * 10L * 3L,
          "nothing failed" = sum(d$failed) == 0L,
          "the result carries the manifest signature" =
            identical(attr(d, "provenance")$signature, m$signature),
          "and its own configuration hash" =
            nchar(attr(d, "provenance")$config_hash) == 64L)
for (col in c("coverage", "capture_all", "cluster_rho_all", "cl10_misalloc",
              "b10_hn_capture", "b10_slots", "n_hn", "n_unrec", "level_fixups",
              "contact_rate", "eligible_rate", "gate", "block", "variant")) {
  stopifnot("every reported column is present" = col %in% names(d))
}

# ---- The lock is a lock --------------------------------------------------
lock <- claim_lock(out)
held <- run(base)
stopifnot("a second writer is refused while the directory is locked" =
            held$status != 0L && grepl("lock", held$text))
stopifnot("and the lock is not removed for it" = dir.exists(lock))
release_lock(lock)
stopifnot("releasing removes it" = !dir.exists(lock))

# ---- A lock whose owner is gone does not end the run ---------------------
# A killed or interrupted run cannot release its lock, because a signal does not
# run `on.exit`. If that lock were permanent, every resume would need a manual
# delete first and the run would not be resumable in any useful sense.
a_gone_pid <- function() {
  for (i in 1:200) {
    p <- sample.int(90000L, 1L) + 9000L
    none <- !length(suppressWarnings(system2("ps", c("-p", p, "-o", "pid="),
                                             stdout = TRUE, stderr = FALSE)))
    if (none) return(p)
  }
  stop("could not find a pid with no process on it")
}
write_lock <- function(pid, host) {
  l <- file.path(out, "run.lock")
  dir.create(l, showWarnings = FALSE)
  writeLines(c(paste("pid", pid), paste("host", host), "when"), file.path(l, "owner"))
  l
}

stopifnot("this process reads as alive" = owner_alive(Sys.getpid()))

l <- write_lock(a_gone_pid(), Sys.info()[["nodename"]])
reclaimed <- run(base)
stopifnot("a lock left by a dead owner is reclaimed and the run proceeds" =
            reclaimed$status == 0L && grepl("Reclaim", reclaimed$text))
stopifnot("no stale copy is left behind" =
            !length(list.files(out, pattern = "run\\.lock\\.stale")))
release_lock(file.path(out, "run.lock"))

l <- write_lock(Sys.getpid(), Sys.info()[["nodename"]])
live <- run(base)
stopifnot("a lock held by a live R process is still refused" =
            live$status != 0L && grepl("running", live$text))
stopifnot("and it is left alone" = dir.exists(l))
release_lock(l)

l <- write_lock(a_gone_pid(), "some-other-machine")
foreign <- run(base)
stopifnot("a lock from another host is refused, alive or not" =
            foreign$status != 0L && grepl("host", foreign$text))
release_lock(l)

l <- file.path(out, "run.lock")
dir.create(l, showWarnings = FALSE)
ownerless <- run(base)
stopifnot("a lock with no owner record is refused" =
            ownerless$status != 0L && grepl("owner", ownerless$text))
release_lock(l)

# ---- The same specification resumes, a different one does not ------------
again <- run(base)
stopifnot("re-running a finished configuration is a no-op" =
            again$status == 0L && grepl("replicates     0 new", again$text))
more <- run(c("--replicates", "3", "--workers", "2", "--blocks", "core",
              "--configs", CFG, "--stream", "smoke", "--out", shQuote(out)))
stopifnot("a longer run appends rather than restarts" =
            more$status == 0L && length(unique(readRDS(f)$rep)) == 3L)

# A different stream is a different seed stream, so it must not append.
other <- run(c("--replicates", "2", "--workers", "2", "--blocks", "core",
               "--configs", CFG, "--stream", "pilot", "--out", shQuote(out)))
stopifnot("a different stream is refused in the same directory" =
            other$status != 0L && grepl("manifest", other$text))
stopifnot("and the completed file is untouched" = length(unique(readRDS(f)$rep)) == 3L)

# ---- A moved scientific source is refused, and named ---------------------
m2 <- m
m2$source_hashes[["index_arms.R"]] <- paste(rep("0", 64), collapse = "")
m2$signature <- digest(m2[setdiff(names(m2), c("signature", "created_utc"))], algo = "sha256")
moved <- manifest_diff(m2, m)
stopifnot("the diff names the component that moved" =
            any(grepl("source_hashes", moved)) && any(grepl("index_arms", moved)))
tmp_dir <- file.path(tempdir(), "manifest_clash")
dir.create(tmp_dir, showWarnings = FALSE)
on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
saveRDS(m2, file.path(tmp_dir, "manifest.rds"), version = 3)
clash <- try(initialize_run(tmp_dir, m), silent = TRUE)
stopifnot("and initialisation refuses it" = inherits(clash, "try-error"))

# ---- Reading is driven by the manifest, not by a glob --------------------
invisible(file.copy(f, file.path(out, "stray_results.rds")))
mf <- manifest_files(out)
stopifnot("a stray rds is not a result" =
            !any(grepl("stray", mf$paths)),
          "and the configurations still to run are listed" =
            length(mf$missing) == nrow(m$grid) - 1L)
unlink(file.path(out, "stray_results.rds"))

# ---- Atomic writes -------------------------------------------------------
# A half-written rds must never appear under a configuration's name.
blocked <- file.path(out, "blocked.rds")
dir.create(blocked)
failed_write <- try(save_atomic(data.frame(a = 1), blocked), silent = TRUE)
stopifnot("a failed rename aborts rather than leaving a partial file" =
            inherits(failed_write, "try-error"),
          "and leaves no pending file behind" =
            !length(list.files(out, pattern = "^pending-")))
unlink(blocked, recursive = TRUE)

# ---- Result files are verified on the way in -----------------------------
stopifnot("a wrong signature is refused" =
            inherits(try(read_verified(f, "not-the-signature", attr(d, "provenance")$config_hash),
                         silent = TRUE), "try-error"))
stopifnot("a wrong configuration hash is refused" =
            inherits(try(read_verified(f, m$signature, "not-the-hash"), silent = TRUE),
                     "try-error"))
stopifnot("and the right pair is accepted" =
            !is.null(read_verified(f, m$signature, attr(d, "provenance")$config_hash)))

cat("test_runner.R: all contracts hold\n")
