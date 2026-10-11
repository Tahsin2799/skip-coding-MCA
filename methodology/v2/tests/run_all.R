# Every v2 contract. Run: Rscript methodology/v2/tests/run_all.R
#
# Run this before any production job. Each file is a separate process, so one
# failure does not hide the others.

root <- c(".", "..", "../..") |>
  Filter(f = \(p) file.exists(file.path(p, "methodology/v2/plasmode_generate.R"))) |>
  (\(x) x[1])()
files <- sort(list.files(file.path(root, "methodology/v2/tests"),
                         pattern = "^test_.*\\.R$", full.names = TRUE))
started <- Sys.time()
failed <- character(0)
for (f in files) {
  out <- system2("Rscript", shQuote(f), stdout = TRUE, stderr = TRUE)
  ok <- !is.null(attr(out, "status")) && attr(out, "status") != 0L
  if (ok) {
    failed <- c(failed, basename(f))
    cat(sprintf("FAIL  %s\n", basename(f)))
    cat(paste0("      ", utils::tail(out, 6), collapse = "\n"), "\n")
  } else {
    cat(sprintf("ok    %s\n", basename(f)))
  }
}
cat(sprintf("\n%d file%s, %d failed, %.0fs\n", length(files),
            if (length(files) == 1) "" else "s", length(failed),
            as.numeric(difftime(Sys.time(), started, units = "secs"))))
if (length(failed)) quit(status = 1L)
