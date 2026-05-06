# ============================================================
# shared_PCA_MEANS_ONLY_SICB_colored.R
# Mean positions only in shared PCA space:
#   1) PC1 vs PC2
#   2) PC2 vs PC3
#   3) PC3 vs PC4
#
# Shared dataset: 1950 habitats + CT time series
# Styling/coloring matched to watershed_landscape_SICB.R
#
# Saves PDFs to:
#   Figures/shared_PCA_means_only_SICB/
# ============================================================

suppressPackageStartupMessages({
  library(geomorph)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

# ============================================================
# Canonical upstream
# ============================================================

source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")
source("R/methods/01_specimen_sampling_study_design/04_subset_CT_timeseries_plus_1950habitats.R")

# Output folder
fig_dir <- file.path("Figures", "shared_PCA_means_only_SICB")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Palette — copied from watershed_landscape_SICB.R logic
# ============================================================

group_levels <- levels(gdf_Fig6$group)

group_palette <- setNames(rep(NA_character_, length(group_levels)), group_levels)
group_palette["CT_1950"] <- "steelblue"
group_palette["CT_1956"] <- "dodgerblue4"
group_palette["CT_1970"] <- "navy"

other <- setdiff(group_levels, c("CT_1950", "CT_1956", "CT_1970"))

if (length(other) > 0) {
  if (
    exists("hab_palette", inherits = TRUE) &&
    is.vector(hab_palette) &&
    !is.null(names(hab_palette))
  ) {
    hab_names <- sub("_1950$", "", other)
    cols <- hab_palette[hab_names]

    if (any(is.na(cols))) {
      missing <- hab_names[is.na(cols)]
      fallback <- scales::hue_pal()(length(missing))
      names(fallback) <- missing
      cols[is.na(cols)] <- fallback[missing]
    }

    names(cols) <- other
    group_palette[other] <- cols
  } else {
    group_palette[other] <- scales::hue_pal()(length(other))
  }
}

# Optional safety check
if (any(is.na(group_palette))) {
  stop("Some group colors are NA. Check group levels and hab_palette names.")
}

# ============================================================
# PCA on residual coordinates in shared space
# ============================================================

pca_Fig7 <- geomorph::gm.prcomp(coords_resid_Fig6)

pct_Fig7 <- 100 * (pca_Fig7$sdev^2 / sum(pca_Fig7$sdev^2))

pc1_lab <- sprintf("PC1 (%.1f%%)", pct_Fig7[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", pct_Fig7[2])
pc3_lab <- sprintf("PC3 (%.1f%%)", pct_Fig7[3])
pc4_lab <- sprintf("PC4 (%.1f%%)", pct_Fig7[4])

pca_Fig7_df <- data.frame(
  specimen = rownames(pca_Fig7$x),
  PC1      = pca_Fig7$x[, 1],
  PC2      = pca_Fig7$x[, 2],
  PC3      = pca_Fig7$x[, 3],
  PC4      = pca_Fig7$x[, 4],
  group    = gdf_Fig6$group,
  habitat  = gdf_Fig6$habitat,
  year     = gdf_Fig6$year,
  SL_mm    = gdf_Fig6$SL_mm,
  logSL    = gdf_Fig6$logSL,
  stringsAsFactors = FALSE
)

stopifnot(identical(pca_Fig7_df$specimen, gdf_Fig6$specimen))

# ============================================================
# Mean positions only
# ============================================================

pca_means_df <- pca_Fig7_df %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    PC3 = mean(PC3, na.rm = TRUE),
    PC4 = mean(PC4, na.rm = TRUE),
    n   = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(group = factor(group, levels = group_levels))

# ============================================================
# Plot helper: means only, journal-size 3 x 2.36 inches
# ============================================================

make_mean_pca_plot <- function(df, xvar, yvar, xlab, ylab) {

  ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_point(
      aes(fill = group),
      shape = 21,
      color = "black",
      size = 1.65,
      stroke = 0.25,
      alpha = 0.95
    ) +
    scale_fill_manual(values = group_palette, guide = "none", drop = FALSE) +
    labs(
      x = xlab,
      y = ylab
    ) +
    coord_equal() +
    theme_classic(base_family = "Helvetica", base_size = 6) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      legend.position = "none",
      axis.title = element_text(size = 6),
      axis.text  = element_text(size = 5),
      axis.line = element_line(linewidth = 0.22),
      axis.ticks = element_line(linewidth = 0.22),
      axis.ticks.length = unit(1.1, "mm"),
      plot.margin = margin(2, 2, 2, 2, unit = "pt")
    )
}

# ============================================================
# Make mean-only plots
# ============================================================

p_pc12_means <- make_mean_pca_plot(
  df = pca_means_df,
  xvar = "PC1",
  yvar = "PC2",
  xlab = pc1_lab,
  ylab = pc2_lab
)

p_pc23_means <- make_mean_pca_plot(
  df = pca_means_df,
  xvar = "PC2",
  yvar = "PC3",
  xlab = pc2_lab,
  ylab = pc3_lab
)

p_pc34_means <- make_mean_pca_plot(
  df = pca_means_df,
  xvar = "PC3",
  yvar = "PC4",
  xlab = pc3_lab,
  ylab = pc4_lab
)

# ============================================================
# Save PDFs: 3 x 2.36 inches
# ============================================================

ggsave(
  filename = file.path(fig_dir, "shared_PCA_means_PC1_PC2_SICB_3in.pdf"),
  plot = p_pc12_means,
  width = 3,
  height = 2.36,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(fig_dir, "shared_PCA_means_PC2_PC3_SICB_3in.pdf"),
  plot = p_pc23_means,
  width = 3,
  height = 2.36,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(fig_dir, "shared_PCA_means_PC3_PC4_SICB_3in.pdf"),
  plot = p_pc34_means,
  width = 3,
  height = 2.36,
  units = "in",
  device = cairo_pdf
)

# Optional PNGs for checking quickly

ggsave(
  filename = file.path(fig_dir, "shared_PCA_means_PC1_PC2_SICB_3in.png"),
  plot = p_pc12_means,
  width = 3,
  height = 2.36,
  units = "in",
  dpi = 1200,
  bg = "white"
)

ggsave(
  filename = file.path(fig_dir, "shared_PCA_means_PC2_PC3_SICB_3in.png"),
  plot = p_pc23_means,
  width = 3,
  height = 2.36,
  units = "in",
  dpi = 1200,
  bg = "white"
)

ggsave(
  filename = file.path(fig_dir, "shared_PCA_means_PC3_PC4_SICB_3in.png"),
  plot = p_pc34_means,
  width = 3,
  height = 2.36,
  units = "in",
  dpi = 1200,
  bg = "white"
)

message("Saved mean-only 3 x 2.36 inch SICB shared PCA plots to: ",
        normalizePath(fig_dir, winslash = "/"))
