# ============================================================
# scripts/1950_raw__centroid_PCA_CsizeColor_habHull.R
# FIGURE 2 (1950): RAW morphospace structured by centroid size + habitat
#   - PCA on coords_1950 (raw; not size-corrected residuals)
#   - Hulls filled by habitat
#   - Points colored by centroid size
# Saves: figures/Fig2_1950_raw_PCA_CsizeColor_habHull.pdf + .png
# ============================================================

# ---------------------------
# Load canonical project objects
# ---------------------------
source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")
source("R/methods/01_specimen_sampling_study_design/02_subset_1950.R")

suppressPackageStartupMessages({
  library(geomorph)
  library(ggplot2)
  library(dplyr)
  library(viridis)   # for scale_color_viridis_c
  library(scales)    # for squish
})

# ---------------------------
# Fallback hull helper (only if not defined in 00_setup_morpho.R)
# ---------------------------
if (!exists("make_hulls")) {
  make_hulls <- function(df, group_col = "habitat", x = "PC1", y = "PC2", min_n = 3) {
    split_df <- split(df, df[[group_col]])
    hulls <- lapply(names(split_df), function(g) {
      d <- split_df[[g]]
      d <- d[complete.cases(d[[x]], d[[y]]), , drop = FALSE]
      if (nrow(d) < min_n) return(NULL)
      idx <- chull(d[[x]], d[[y]])
      out <- d[idx, , drop = FALSE]
      out[[group_col]] <- g
      out
    })
    do.call(rbind, hulls)
  }
}

# ============================================================
# 1950 RAW PCA
# ============================================================

# --- Compute PCA ---
pca_raw_1950 <- gm.prcomp(coords_1950)
pct_raw_1950 <- 100 * (pca_raw_1950$sdev^2 / sum(pca_raw_1950$sdev^2))

pc1_lab <- sprintf("PC1 (%.1f%%)", pct_raw_1950[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", pct_raw_1950[2])

# --- Canonical plotting df ---
pca1950_raw_df <- data.frame(
  specimen = rownames(pca_raw_1950$x),
  PC1      = pca_raw_1950$x[, 1],
  PC2      = pca_raw_1950$x[, 2],
  habitat  = factor(gdf_1950$habitat),
  Csize    = gdf_1950$Csize,
  logCsize = gdf_1950$logCsize,
  stringsAsFactors = FALSE
)

stopifnot(identical(pca1950_raw_df$specimen, gdf_1950$specimen))

# --- Hulls by habitat ---
hulls_1950 <- make_hulls(
  df = pca1950_raw_df,
  group_col = "habitat",
  x = "PC1",
  y = "PC2"
)

# ============================================================
# FIGURE 2: centroid-size-colored points + habitat hulls
# ============================================================

p_fig2 <- ggplot(pca1950_raw_df, aes(PC1, PC2)) +
  
  geom_polygon(
    data = hulls_1950,
    aes(fill = habitat, group = habitat),
    alpha = 0.18,
    color = NA
  ) +
  
  geom_polygon(
    data = hulls_1950,
    aes(group = habitat),
    fill = NA,
    linewidth = 0.6,
    color = "black"
  ) +
  
  geom_point(
    aes(color = Csize),
    size = 2.2,
    alpha = 0.90
  ) +
  
  scale_color_viridis_c(
    name   = "Centroid size",
    option = "C",
    limits = range(pca1950_raw_df$Csize, na.rm = TRUE),
    oob    = scales::squish
  ) +
  
  scale_fill_manual(
    name   = "Habitat",
    values = hab_palette,
    drop   = FALSE
  ) +
  
  labs(
    title = "Raw morphospace structured by centroid size and habitat (1950)",
    x = pc1_lab,
    y = pc2_lab
  ) +
  
  coord_equal() +
  theme_classic(base_family = "Helvetica", base_size = 9) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 9),
    axis.text  = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8)
  )

print(p_fig2)

# ============================================================
# Save (PDF + PNG) into figures/
# ============================================================

fig_dir <- "figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

fig_base <- file.path(fig_dir, "Fig2_1950_raw_PCA_CsizeColor_habHull")

ggsave(
  filename = paste0(fig_base, ".pdf"),
  plot     = p_fig2,
  device   = cairo_pdf,
  width    = 7,
  height   = 5,
  units    = "in"
)

ggsave(
  filename = paste0(fig_base, ".png"),
  plot     = p_fig2,
  width    = 7,
  height   = 5,
  units    = "in",
  dpi      = 600
)

if (exists("VERBOSE") && isTRUE(VERBOSE)) {
  cat("\nSaved Figure 2 to:\n")
  cat("  ", paste0(fig_base, ".pdf"), "\n", sep = "")
  cat("  ", paste0(fig_base, ".png"), "\n", sep = "")
}

FIG2_OBJECTS <- c("pca_raw_1950", "pct_raw_1950", "pca1950_raw_df", "hulls_1950", "p_fig2")
if (exists("VERBOSE") && isTRUE(VERBOSE)) {
  cat("\nFigure 2 objects created:\n")
  print(FIG2_OBJECTS)
}