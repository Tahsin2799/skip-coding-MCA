# Rebuild outputs/methodology/donors_harmonized_geo.rds, the covariate cache that
# every v03 run reads, from the four DHS person recode files.
#
#   Rscript methodology/v2/build_donor_cache.R [--dhs <directory>] [--out <file>]
#
# <directory> holds the DHS Stata files in the folders donor_registry() names,
# IAPR7EDT/IAPR7EFL.DTA and so on; without --dhs it is DONOR_ROOT from
# harmonize_donors.R. That file is hashed into every run manifest, so the path is
# overridden here rather than edited there.
#
# The cache keeps every person aged 18 and over, stores factor levels without
# the variable prefix that harmonize_donor() adds (the runner re-applies it on
# load), and carries region and urban as its last two columns. Rebuilt this way
# on 25 September 2026, it reproduced the file hashed into the production
# manifest byte for byte: sha256 c9ac7c45a02be2563164ed6187353f8f16557f855a732ac253e9ad9ba4605e11.

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
})

root <- c(".", "..", "../..") |>
  Filter(f = \(p) file.exists(file.path(p, "methodology/v2/harmonize_donors.R"))) |>
  (\(x) x[1])()
source(file.path(root, "methodology/v2/harmonize_donors.R"))

argv <- commandArgs(trailingOnly = TRUE)
arg <- function(flag, default) {
  i <- match(flag, argv)
  if (is.na(i)) default else argv[i + 1]
}
DONOR_ROOT <- arg("--dhs", DONOR_ROOT)
out <- arg("--out", file.path(root, "outputs/methodology/donors_harmonized_geo.rds"))

strip_prefix <- function(d) {
  for (v in INDEX_INPUTS) levels(d[[v]]) <- sub(paste0("^", v, "_"), "", levels(d[[v]]))
  d
}

donors <- donor_registry()$donor
cache <- lapply(setNames(donors, donors), function(n) {
  harmonize_donor(n, age_min = 18, age_max = Inf) |>
    strip_prefix() |>
    relocate(region, urban, .after = last_col())
})

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(cache, out)
cat(out, "\nsha256", digest::digest(file = out, algo = "sha256"), "\n")
