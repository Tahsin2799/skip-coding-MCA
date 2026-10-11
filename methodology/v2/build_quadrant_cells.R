# Rebuild methodology/v2/quadrant_cells.csv from the harmonized donors.
#
# The membership file used to be produced ad hoc, which meant the coherence
# values stored in it could not be checked against the function that computes
# them. They could not: the file was written with the uncorrected Cramer's V,
# so 18 of 63 cells sat on the wrong side of the split. This script makes the
# file reproducible and is the only thing that should write it.
#
#   Rscript methodology/v2/build_quadrant_cells.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(cli)
})

root <- c(".", "..", "../..") |>
  Filter(f = \(p) file.exists(file.path(p, "methodology/v2/quadrant_populations.R"))) |>
  (\(x) x[1])()
source(file.path(root, "methodology/v2/harmonize_donors.R"))
source(file.path(root, "methodology/v2/quadrant_populations.R"))

MIN_CELL <- 700L          # below this the MCA on six inputs is not stable
TARGET   <- 20000L        # the design's target population size

harmonized <- lapply(readRDS(file.path(root, "outputs/methodology/donors_harmonized_geo.rds")),
                     prefix_levels)
pool <- cell_pool(harmonized)
parts <- split(pool, pool$cell)
parts <- parts[vapply(parts, nrow, integer(1)) >= MIN_CELL]
cli_alert_info("profiling {length(parts)} cells at n >= {MIN_CELL}")

prof <- bind_rows(lapply(names(parts), \(cn) {
  d <- parts[[cn]]
  tibble(cell = cn, n = nrow(d),
         deprivation = round(mean(deprivation_score(d)), 4),
         coherence   = round(coherence(d), 4),
         adj_inertia = round(adjusted_inertia_dim1(d), 1))
})) |> filter(is.finite(coherence), is.finite(deprivation))

# Split at the record-weighted medians so each quadrant is populated rather
# than at round numbers chosen by eye.
dep_cut <- as.numeric(matrixStats::weightedMedian(prof$deprivation, prof$n))
coh_cut <- as.numeric(matrixStats::weightedMedian(prof$coherence, prof$n))

prof <- prof |>
  mutate(quadrant = case_when(
    deprivation <  dep_cut & coherence <  coh_cut ~ "Q1",
    deprivation <  dep_cut & coherence >= coh_cut ~ "Q2",
    deprivation >= dep_cut & coherence <  coh_cut ~ "Q3",
    TRUE                                          ~ "Q4"))

sizes <- prof |> count(quadrant, wt = n, name = "records")
cli_alert_info("split at deprivation {round(dep_cut,4)}, coherence {round(coh_cut,4)}")
print(as.data.frame(prof |> group_by(quadrant) |>
  summarise(cells = n(), records = sum(n),
            dep = sprintf("%.3f-%.3f", min(deprivation), max(deprivation)),
            coh = sprintf("%.3f-%.3f", min(coherence), max(coherence)),
            .groups = "drop")), row.names = FALSE)

short <- sizes$quadrant[sizes$records < TARGET]
if (length(short)) {
  cli_abort(c("Quadrant{?s} {.val {short}} fall{?s/} short of the {TARGET} target size.",
              i = "Lower MIN_CELL or move a split before freezing."))
}

out <- file.path(root, "methodology/v2/quadrant_cells.csv")
write.csv(prof |> select(cell, quadrant, n, deprivation, coherence, adj_inertia),
          out, row.names = FALSE)
cli_alert_success("wrote {.path {out}} with {nrow(prof)} cells")
