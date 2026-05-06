# ============================================================
# Swift_Quabbin_1950_PCA_clean.R
# 
# Clean PCA plots:
#   - PC1 vs PC2
#   - PC2 vs PC3
# 
# PCA computed on ALL 1950 habitats
# Plots show Quabbin + Swift only
# ============================================================

suppressPackageStartupMessages({
  library(geomorph)
  library(ggplot2)
  library(dplyr)
})

# ---------------------------
# Load canonical data
# ---------------------------
source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")
source("R/methods/01_specimen_sampling_study_design/02_subset_1950.R")

# ---------------------------
# Output directory
# ---------------------------
OUTDIR <- file.path("Figures", "Swift and Quabbin")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# PCA on ALL 1950 residual data
# ============================================================

pca_1950 <- gm.prcomp(coords_resid_1950)
pct_1950 <- 100 * (pca_1950$sdev^2 / sum(pca_1950$sdev^2))

pc_df <- data.frame(
  specimen = rownames(pca_1950$x),
  PC1 = pca_1950$x[, 1],
  PC2 = pca_1950$x[, 2],
  PC3 = pca_1950$x[, 3],
  habitat = gdf_1950$habitat,
  collection_id = gdf_1950$collection_id,
  stringsAsFactors = FALSE
)

# ============================================================
# Subset: Quabbin + Swift
# ============================================================

plot_df <- pc_df %>%
  filter(habitat %in% c("Quabbin", "Swift River")) %>%
  mutate(
    point_group = case_when(
      habitat == "Swift River" ~ "Swift River",
      habitat == "Quabbin" & collection_id == "F2379" ~ "bottom",
      habitat == "Quabbin" & collection_id == "F2374" ~ "middle",
      habitat == "Quabbin" & collection_id == "F2385" ~ "top",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(point_group))

# ============================================================
# Palettes
# ============================================================

point_palette <- c(
  "bottom"      = "#EC7D3C",
  "middle"      = "#F0E442",
  "top"         = "#A6761D",
  "Swift River" = "tomato"
)

fill_palette <- c(
  "Quabbin"     = "darkgoldenrod2",
  "Swift River" = "tomato"
)

line_palette <- c(
  "Quabbin"     = "darkgoldenrod2",
  "Swift River" = "tomato3"
)

# ============================================================
# Hull function
# ============================================================

make_closed_hulls <- function(df, xvar, yvar) {
  df %>%
    group_by(habitat) %>%
    slice(chull(.data[[xvar]], .data[[yvar]])) %>%
    slice(c(1:n(), 1)) %>%
    ungroup()
}

# ============================================================
# Base plot function (NO LEGEND)
# ============================================================

make_plot <- function(plot_df, hull_df, xvar, yvar, xlab, ylab) {
  
  ggplot() +
    
    geom_polygon(
      data = hull_df,
      aes(x = .data[[xvar]], y = .data[[yvar]], group = habitat, fill = habitat),
      alpha = 0.12,
      color = NA
    ) +
    
    geom_path(
      data = hull_df,
      aes(x = .data[[xvar]], y = .data[[yvar]], group = habitat, color = habitat),
      linewidth = 0.6,
      show.legend = FALSE
    ) +
    
    geom_point(
      data = plot_df,
      aes(x = .data[[xvar]], y = .data[[yvar]], color = point_group),
      size = 1.2,
      alpha = 0.9
    ) +
    
    scale_fill_manual(values = fill_palette) +
    scale_color_manual(values = c(line_palette, point_palette)) +
    
    labs(x = xlab, y = ylab) +
    
    coord_equal() +
    
    theme_classic(base_size = 6) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 6),
      axis.text  = element_text(size = 5),
      axis.line  = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      plot.margin = margin(2, 2, 2, 2, unit = "pt")
    )
}

# ============================================================
# PC1 vs PC2
# ============================================================

hull_12 <- make_closed_hulls(plot_df, "PC1", "PC2")

p12 <- make_plot(
  plot_df,
  hull_12,
  "PC1", "PC2",
  sprintf("PC1 (%.1f%%)", pct_1950[1]),
  sprintf("PC2 (%.1f%%)", pct_1950[2])
)

# ============================================================
# PC2 vs PC3
# ============================================================

hull_23 <- make_closed_hulls(plot_df, "PC2", "PC3")

p23 <- make_plot(
  plot_df,
  hull_23,
  "PC2", "PC3",
  sprintf("PC2 (%.1f%%)", pct_1950[2]),
  sprintf("PC3 (%.1f%%)", pct_1950[3])
)

# ============================================================
# Save
# ============================================================

ggsave(
  file.path(OUTDIR, "Fig_1950_QS_PC12.pdf"),
  p12,
  width = 3,
  height = 2.4,
  units = "in"
)

ggsave(
  file.path(OUTDIR, "Fig_1950_QS_PC12.png"),
  p12,
  width = 3,
  height = 2.4,
  units = "in",
  dpi = 1200
)

ggsave(
  file.path(OUTDIR, "Fig_1950_QS_PC23.pdf"),
  p23,
  width = 3,
  height = 2.4,
  units = "in"
)

ggsave(
  file.path(OUTDIR, "Fig_1950_QS_PC23.png"),
  p23,
  width = 3,
  height = 2.4,
  units = "in",
  dpi = 1200
)

message("Saved to: ", normalizePath(OUTDIR))