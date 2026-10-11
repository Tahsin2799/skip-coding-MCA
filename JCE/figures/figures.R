# Figures 1 to 3 and Supplementary Figure S1, from production summaries, the
# production-linked theory summary and the documented correlation correction.
#
#   Rscript JCE/figures/figures.R
#
# Writes JCE/figures/fig{1,2,3,S1}.pdf (vector, for submission) and .png (300
# dpi). Figure 1 is drawn by JCE/theory/theory_figure.R. Missing fuzzy average in
# Figure 2 comes from the post hoc production-linked run.
# Colours: the contact levels are an ordinal one-hue ramp and the two series of
# Figure 2 are categorical slots 1 and 2 of the dataviz reference palette, both
# validated with its validate_palette.js on 25 September 2026.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

ds <- readRDS("outputs/methodology/production_v03_displays.rds")
review <- readRDS("JCE/analysis/review_analysis.rds")
theory <- readRDS("JCE/theory/theory_link_summary_r1000.rds")
source("JCE/theory/theory_figure.R")
out <- "JCE/figures"
dir.create(out, showWarnings = FALSE)

METHODS <- c(
  oracle = "Oracle (reference)",
  simple_score = "Deprivation count",
  impute_adaptive = "Imputation, target models",
  impute_frozen = "Imputation, development models",
  fuzzy_average = "Missing fuzzy average (post hoc)",
  drop_items = "Items dropped",
  raw_items = "Supervised model (reference)",
  zero = "Skip coded as \"no\"",
  passive_mca = "Subset MCA",
  complete_case = "Complete-case analysis",
  category = "Missing single (skip as its own category)"
)
FILTERS <- c(contact = "Contact filter", diagnosis = "Diagnosis filter")
POPULATIONS <- c(Q3 = "Pooled stratum Q3, 18 to 49", bd_allage = "Bangladesh, 18 to 95")
CONTACT <- c(c20 = "20%", c50 = "50%", c80 = "80%")
CONTACT_COLOURS <- c("20%" = "#86b6ef", "50%" = "#2a78d6", "80%" = "#104281")
INK <- "#0b0b0b"
GRID <- "#e6e5e1"

label_methods <- function(d) {
  d |> mutate(method = factor(METHODS[as.character(arm)], levels = rev(METHODS)))
}

theme_fig <- function() {
  theme_minimal(base_size = 9) +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(colour = GRID, linewidth = 0.3),
          axis.text = element_text(colour = INK),
          strip.text = element_text(colour = INK, face = "bold", hjust = 0),
          legend.position = "top",
          legend.justification = "left",
          panel.spacing.x = unit(1.4, "lines"),
          plot.background = element_rect(fill = "white", colour = NA))
}

# cairo_pdf embeds the fonts, which the journal requires of vector figures;
# the default pdf device only references Helvetica.
save_fig <- function(p, name, width, height) {
  ggsave(file.path(out, paste0(name, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(out, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
}

# Figure 1: the capture ratio against axis share and targeting loss.
save_fig(capture_plot(theory), "fig1", 7, 3.4)

# Figure 2: every core scenario, so a reversal cannot hide inside an average.
f1 <- bind_rows(ds$cellwise |> mutate(arm = as.character(arm)),
                theory$cellwise |> filter(arm == "fuzzy_average") |>
                  mutate(arm = as.character(arm))) |>
  label_methods() |>
  mutate(filter = factor(FILTERS[gate], levels = FILTERS),
         population = factor(POPULATIONS[population], levels = POPULATIONS),
         contact = factor(CONTACT[as.character(setting)], levels = CONTACT))

p1 <- ggplot(f1, aes(cluster_rho_all, method, fill = contact)) +
  geom_point(shape = 21, size = 1.8, stroke = 0.25, colour = "white",
             position = position_dodge(width = 0.6)) +
  facet_grid(filter ~ population) +
  scale_fill_manual(values = CONTACT_COLOURS, name = "Contact with care") +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  scale_x_continuous(limits = c(-0.05, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "Rank correlation with true cluster vulnerability, whole target population",
       y = NULL) +
  theme_fig()
save_fig(p1, "fig2", 7.2, 6.8)

# Figure 3: what a validation restricted to scoreable respondents would report.
f2 <- ds$coverage_gap |>
  label_methods() |>
  mutate(method = droplevels(method)) |>
  mutate(filter = factor(FILTERS[gate], levels = FILTERS),
         contact = factor(paste(CONTACT[as.character(setting)], "contact"),
                          levels = paste(CONTACT, "contact"))) |>
  select(filter, contact, method, `Whole target population` = rho_all,
         `People the method can score` = rho_scored)

# The scored-only value is a hollow ring drawn first and the whole-target value a
# filled dot drawn on top, so where the two agree both stay visible.
READ_ON <- c("Whole target population", "People the method can score")
f2_long <- f2 |>
  pivot_longer(all_of(READ_ON), names_to = "read_on", values_to = "rho") |>
  mutate(read_on = factor(read_on, levels = rev(READ_ON))) |>
  arrange(read_on)

p2 <- ggplot(f2, aes(y = method)) +
  geom_segment(aes(x = `Whole target population`, xend = `People the method can score`),
               colour = "#b9b8b3", linewidth = 0.5) +
  geom_point(data = f2_long, aes(x = rho, colour = read_on, shape = read_on, size = read_on),
             stroke = 0.9) +
  facet_grid(filter ~ contact) +
  scale_colour_manual(values = c(`Whole target population` = "#2a78d6",
                                 `People the method can score` = "#eb6834"),
                      breaks = READ_ON, name = "Rank correlation computed on") +
  scale_shape_manual(values = c(`Whole target population` = 16,
                                `People the method can score` = 1),
                     breaks = READ_ON, name = "Rank correlation computed on") +
  scale_size_manual(values = c(`Whole target population` = 1.5,
                               `People the method can score` = 2.8),
                    breaks = READ_ON, name = "Rank correlation computed on") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Rank correlation with true cluster vulnerability", y = NULL) +
  theme_fig()
save_fig(p2, "fig3", 7.2, 5.6)

# Figure S1: do the statistics a validation study reports track targeting?
APPARENT <- c(
  auc_index_scored = "AUC of the index",
  auc_model_scored = "AUC of the model",
  log_or_top = "Top-tier odds ratio\nDevelopment",
  log_or_top_target = "Top-tier odds ratio\nTarget",
  beta_within = "Within-cluster slope",
  beta_between = "Between-cluster slope"
)
TARGETING <- c(
  cluster_rho_all = "Cluster rank correlation",
  cl10_misalloc = "Inclusion error, 10% of clusters",
  misallocation_all = "Inclusion error, top tenth of clusters",
  cl10_precision = "Selected clusters in the true top tenth",
  cl10_hn_reached = "Most vulnerable fifth in selected clusters",
  b05_hn_capture = "Most vulnerable fifth selected, 5% budget",
  b10_hn_capture = "Most vulnerable fifth selected, 10% budget",
  b10_hn_shortfall_per1000 = "Shortfall per 1,000, 10% budget",
  b20_hn_capture = "Most vulnerable fifth selected, 20% budget"
)
f3 <- review$families |>
  filter(a %in% names(APPARENT)) |>
  mutate(apparent = factor(APPARENT[a], levels = APPARENT),
         targeting = factor(TARGETING[t], levels = rev(TARGETING)))
stopifnot(!anyNA(f3$apparent), !anyNA(f3$targeting), nrow(f3) == 54)

p3 <- ggplot(f3, aes(y = targeting)) +
  geom_vline(xintercept = 0, colour = "#8a8984", linewidth = 0.4) +
  geom_segment(aes(x = min_rho, xend = max_rho), colour = "#b9b8b3", linewidth = 0.5) +
  geom_point(aes(x = median_rho), shape = 21, size = 1.9, stroke = 0.25,
             colour = "white", fill = INK) +
  facet_wrap(~apparent, ncol = 3) +
  scale_x_continuous(limits = c(-0.7, 0.7), breaks = seq(-0.6, 0.6, 0.3)) +
  labs(x = "Correlation across methods and replicates within a scenario",
       y = NULL) +
  theme_fig()
save_fig(p3, "figS1", 7.2, 5.6)

# Figure 4: the post hoc extensions at the Q3 reference configuration
# (JCE/extensions). Series take categorical slots 1 to 4 of the reference
# palette (validated for adjacent pairs, light surface, 11 October 2026); slots
# 3 and 4 sit below 3:1 contrast, so each series also has its own point shape
# and Tables S11 and S12 give the values.
ext <- readRDS("JCE/extensions/extensions_summary.rds")
EXT_SERIES <- c(impute_adaptive = "Imputation, target models",
                impute_frozen = "Imputation, development models",
                rimca = "Regularized iterative MCA",
                fuzzy_average = "Missing fuzzy average")
EXT_COLOURS <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100"), EXT_SERIES)
EXT_SHAPES <- setNames(c(16, 17, 15, 18), EXT_SERIES)
ITEM_FIT <- c("Maximum likelihood", "FLIC")
ext_series <- function(d) {
  d |> filter(arm %in% c(names(EXT_SERIES), "impute_flic_adaptive", "impute_flic_frozen")) |>
    mutate(series = factor(EXT_SERIES[sub("flic_", "", arm)], levels = EXT_SERIES),
           fit = factor(ifelse(grepl("flic", arm), ITEM_FIT[2], ITEM_FIT[1]), levels = ITEM_FIT))
}
ext_scales <- list(
  geom_hline(yintercept = 0, colour = "#8a8984", linewidth = 0.4),
  geom_line(linewidth = 0.6),
  geom_point(size = 1.8),
  scale_colour_manual(values = EXT_COLOURS, name = NULL),
  scale_shape_manual(values = EXT_SHAPES, name = NULL),
  scale_linetype_manual(values = c("solid", "22"), name = "Item models"),
  guides(colour = guide_legend(nrow = 2, order = 1), shape = guide_legend(nrow = 2, order = 1),
         linetype = guide_legend(nrow = 2, order = 2, override.aes = list(colour = INK, shape = NA)))
)

f4a <- ext$sweep |> ext_series() |> mutate(filter = factor(FILTERS[gate], levels = FILTERS))
p4a <- ggplot(f4a, aes(loading, diff, colour = series, shape = series, linetype = fit,
                       group = interaction(series, fit))) +
  geom_vline(xintercept = 0.25, colour = "#b9b8b3", linewidth = 0.4, linetype = "22") +
  ext_scales +
  facet_wrap(~filter) +
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = c("0", ".25", ".5", ".75", "1")) +
  labs(x = "Advice loading on vulnerability (production .25)",
       y = "Difference from dropping the items") +
  theme_fig()

f4b <- ext$devsize |> ext_series() |>
  left_join(ext$devsize_models, by = "n_dev") |>
  mutate(size = factor(sprintf("%s\n(%.1f)", format(n_dev, big.mark = ","), events_per_parameter)))
p4b <- ggplot(f4b, aes(size, diff, colour = series, shape = series, linetype = fit,
                       group = interaction(series, fit))) +
  ext_scales +
  facet_wrap(~"Diagnosis filter") +
  labs(x = "Development sample\n(events per item-model parameter)", y = NULL) +
  theme_fig()

CODINGS <- c(category = "Missing single", passive_mca = "Subset MCA", zero = "Skip coded as \"no\"")
DEV_CONTACT <- c(c20 = "20%", c50 = "50%", c80 = "80%")
f4c <- ext$devcontact |> filter(arm %in% names(CODINGS)) |>
  mutate(coding = factor(CODINGS[arm], levels = CODINGS),
         filter = factor(FILTERS[gate], levels = FILTERS),
         target = as.numeric(sub("c", "", setting)),
         development = factor(DEV_CONTACT[dev_setting], levels = DEV_CONTACT))
dropped <- ext$devcontact |> filter(arm == "drop_items") |> summarise(rho = mean(rho)) |> pull(rho)
p4c <- ggplot(f4c, aes(target, rho, colour = development)) +
  geom_hline(yintercept = dropped, colour = "#8a8984", linewidth = 0.4, linetype = "22") +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.8) +
  facet_grid(filter ~ coding) +
  scale_colour_manual(values = setNames(unname(CONTACT_COLOURS), DEV_CONTACT),
                      name = "Contact with care in development") +
  scale_x_continuous(breaks = c(20, 50, 80), labels = c("20%", "50%", "80%"), expand = expansion(0.08)) +
  scale_y_continuous(limits = c(-0.05, 0.9), breaks = seq(0, 0.8, 0.2)) +
  labs(x = "Contact with care in the target", y = "Cluster rank correlation") +
  theme_fig()

p4 <- patchwork::wrap_plots(
  patchwork::wrap_plots(p4a, p4b, widths = c(2, 1), guides = "collect") &
    theme(legend.position = "top", legend.justification = "left", legend.title.position = "top"),
  p4c, ncol = 1, heights = c(1, 1.15)) +
  patchwork::plot_annotation(tag_levels = "A")
save_fig(p4, "fig4", 7.2, 7.6)

# Figure S2: the imposed data-generating structure, laid out by Graphviz from
# figS2_dag.dot (supplement A3).
dag <- file.path(out, "figS2_dag.dot")
stopifnot(system2("dot", c("-Tpdf", dag, "-o", file.path(out, "figS2.pdf"))) == 0,
          system2("dot", c("-Tpng", "-Gdpi=300", dag, "-o", file.path(out, "figS2.png"))) == 0)
