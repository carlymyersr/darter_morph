# ============================================================
# Scripts/watershed_landscape_SICB.R
# PCA of residual coordinates:
#   1) PC1 vs PC2
#   2) PC2 vs PC3
#   3) PC3 vs PC4
#
# Journal-size version:
#   - Saves PDFs at 3 x 2.36 inches
#   - No title
#   - No legend
#   - Reduced point size, hull linewidth, axis text, axis lines, and ticks
#
# Saves PDFs to:
#   Figures/watershed_landscape_SICB/
# ============================================================

suppressPackageStartupMessages({
  library(geomorph)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

# Canonical upstream
source("R/00_setup_morpho.R")
source("R/01_build_metadata.R")
source("R/04_subset_CT_timeseries_plus_1950habitats.R")

# Output folder
fig_dir <- file.path("Figures", "watershed_landscape_SICB")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helper: convex hulls
# ============================================================

make_hulls <- function(df, group_col, x, y, min_n = 3) {
  group_vals <- unique(df[[group_col]])
  
  hull_list <- lapply(group_vals, function(g) {
    sub <- df[df[[group_col]] == g, , drop = FALSE]
    if (nrow(sub) < min_n) return(NULL)
    
    h <- chull(sub[[x]], sub[[y]])
    sub[h, , drop = FALSE]
  })
  
  hull_df <- do.call(rbind, hull_list)
  rownames(hull_df) <- NULL
  hull_df
}

# ============================================================
# Palette
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

# ============================================================
# PCA on residual coordinates
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
# Plot helper: journal-size 3 x 2.36 inches
# ============================================================

make_pca_plot <- function(df, xvar, yvar, xlab, ylab) {
  
  hull_df <- make_hulls(
    df = df,
    group_col = "group",
    x = xvar,
    y = yvar,
    min_n = 3
  )
  
  ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_polygon(
      data = hull_df,
      aes(group = group, fill = group),
      alpha = 0.12,
      color = NA
    ) +
    geom_polygon(
      data = hull_df,
      aes(group = group, color = group),
      fill = NA,
      linewidth = 0.22
    ) +
    geom_point(
      aes(color = group),
      size = 0.70,
      alpha = 0.85
    ) +
    scale_color_manual(values = group_palette, guide = "none") +
    scale_fill_manual(values = group_palette, guide = "none") +
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
# Make plots
# ============================================================

p_pc12 <- make_pca_plot(
  df = pca_Fig7_df,
  xvar = "PC1",
  yvar = "PC2",
  xlab = pc1_lab,
  ylab = pc2_lab
)

p_pc23 <- make_pca_plot(
  df = pca_Fig7_df,
  xvar = "PC2",
  yvar = "PC3",
  xlab = pc2_lab,
  ylab = pc3_lab
)

p_pc34 <- make_pca_plot(
  df = pca_Fig7_df,
  xvar = "PC3",
  yvar = "PC4",
  xlab = pc3_lab,
  ylab = pc4_lab
)

# ============================================================
# Save PDFs: journal-size 3 x 2.36 inches
# ============================================================

ggsave(
  filename = file.path(fig_dir, "watershed_landscape_PC1_PC2_SICB_3in.pdf"),
  plot = p_pc12,
  width = 3,
  height = 2.36,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(fig_dir, "watershed_landscape_PC2_PC3_SICB_3in.pdf"),
  plot = p_pc23,
  width = 3,
  height = 2.36,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(fig_dir, "watershed_landscape_PC3_PC4_SICB_3in.pdf"),
  plot = p_pc34,
  width = 3,
  height = 2.36,
  units = "in",
  device = cairo_pdf
)

message("Saved 3 x 2.36 inch SICB watershed landscape PCA plots to: ",
        normalizePath(fig_dir, winslash = "/"))