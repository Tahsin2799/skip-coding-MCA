# Capture ratio against the share of the index axis taken by the filtered items,
# and against targeting loss, for every production configuration and skip coding.
# Figure 1 of the manuscript; JCE/figures/figures.R sources this file.
#
#   Rscript JCE/theory/theory_figure.R JCE/theory/theory_link_summary_r1000.rds
#
# Run directly, it writes JCE/theory/fig_capture.{pdf,png}. Colours are
# categorical slots 1 to 3 of the dataviz reference palette, validated on
# 7 October 2026 (aqua is below 3:1 against the surface, so codings also differ
# by shape and carry a legend).

suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})

capture_plot <- function(s) {
  coding_labels <- c(category = "Missing single", subset = "Subset MCA",
                     zero = "Skip coded as \"no\"")
  colours <- setNames(c("#2a78d6", "#eb6834", "#1baf7a"), coding_labels)
  shapes <- setNames(c(16, 17, 15), coding_labels)
  ink <- "#0b0b0b"
  grid <- "#e6e5e1"

  d <- s$by_cfg |>
    mutate(coding = factor(coding_labels[coding], levels = coding_labels)) |>
    pivot_longer(c(block_share, loss), names_to = "panel", values_to = "y") |>
    mutate(panel = factor(panel, levels = c("block_share", "loss"),
                          labels = c("A. R² of the index axis with the filtered block's axis",
                                     "B. Targeting loss against dropping the items")))

  ggplot(d, aes(kappa, y, colour = coding, shape = coding)) +
    geom_vline(xintercept = 1, colour = ink, linewidth = 0.3, linetype = "dashed") +
    geom_point(size = 2.2, alpha = 0.85) +
    facet_wrap(~panel, ncol = 2, scales = "free_y") +
    scale_colour_manual(values = colours, name = NULL) +
    scale_shape_manual(values = shapes, name = NULL) +
    scale_x_log10(breaks = c(0.75, 1, 1.5, 2, 3, 5)) +
    labs(x = expression("Inertia ratio " * kappa == k * lambda[A] / (P * lambda[X]) * ", development data only"),
         y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = grid, linewidth = 0.3),
          axis.text = element_text(colour = ink),
          strip.text = element_text(colour = ink, face = "bold", hjust = 0),
          legend.position = "top", legend.justification = "left",
          plot.background = element_rect(fill = "white", colour = NA))
}

if (sys.nframe() == 0L) {
  src <- commandArgs(trailingOnly = TRUE)[1]
  if (is.na(src)) src <- "JCE/theory/theory_link_summary_r1000.rds"
  p <- capture_plot(readRDS(src))
  ggsave("JCE/theory/fig_capture.pdf", p, width = 7, height = 3.4)
  ggsave("JCE/theory/fig_capture.png", p, width = 7, height = 3.4, dpi = 300)
}
