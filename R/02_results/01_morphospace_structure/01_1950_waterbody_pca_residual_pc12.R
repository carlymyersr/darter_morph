# ============================================================
# scripts/1950_PCA_residual_Fig3.R
# FIGURE 3: Habitat-associated variation after size correction
# PCA on size-corrected shapes (1950; residuals from procD.lm)
# Size variable: gdf_1950$size_for_allometry (set in 01_build_metadata.R)
# Residuals computed in: R/02_subset_1950.R
# ============================================================

# ---------------------------
# Load canonical project objects
# ---------------------------
source("R/00_setup_morpho.R")
source("R/01_build_metadata.R")
source("R/02_subset_1950.R")

suppressPackageStartupMessages({
  library(geomorph)
  library(ggplot2)
  library(dplyr)
})

# ---------------------------
# Output directory
# ---------------------------
FIG_DIR <- "figures"
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

# ---------------------------
# Preconditions (now that residuals are created in 02_subset_1950)
# ---------------------------
if (!exists("coords_resid_1950")) stop("coords_resid_1950 not found. Did R/02_subset_1950.R run fully?")
if (!exists("fit_allo_1950"))     stop("fit_allo_1950 not found. Did R/02_subset_1950.R run fully?")
if (!exists("make_hulls"))        stop("make_hulls() not found. Keep it in R/00_setup_morpho.R.")
stopifnot(identical(dimnames(coords_resid_1950), dimnames(coords_1950)))

# Nice label for plot subtitle
size_label <- if ("size_label" %in% names(gdf_1950) && !all(is.na(gdf_1950$size_label))) {
  unique(na.omit(gdf_1950$size_label))[1]
} else if ("size_label" %in% names(gdf) && !all(is.na(gdf$size_label))) {
  unique(na.omit(gdf$size_label))[1]
} else {
  "size_for_allometry"
}

# ============================================================
# 1950 RESIDUAL PCA
# ============================================================

pca_resid_1950 <- gm.prcomp(coords_resid_1950)
pct_resid_1950 <- 100 * (pca_resid_1950$sdev^2 / sum(pca_resid_1950$sdev^2))

pc1_resid_lab <- sprintf("PC1 (%.1f%%)", pct_resid_1950[1])
pc2_resid_lab <- sprintf("PC2 (%.1f%%)", pct_resid_1950[2])

# --- Canonical plotting df ---
pca1950_resid_df <- data.frame(
  specimen = rownames(pca_resid_1950$x),
  PC1      = pca_resid_1950$x[, 1],
  PC2      = pca_resid_1950$x[, 2],
  habitat  = factor(gdf_1950$habitat),
  stringsAsFactors = FALSE
)
stopifnot(identical(pca1950_resid_df$specimen, gdf_1950$specimen))

# --- Hulls by habitat ---
hull1950_resid_df <- make_hulls(pca1950_resid_df, group_col = "habitat", x = "PC1", y = "PC2")

# ============================================================
# FIGURE 3
# ============================================================

p_fig3 <- ggplot(pca1950_resid_df, aes(PC1, PC2)) +
  geom_polygon(
    data = hull1950_resid_df,
    aes(group = habitat, fill = habitat),
    alpha = 0.15,
    color = NA
  ) +
  geom_polygon(
    data = hull1950_resid_df,
    aes(group = habitat, color = habitat),
    fill = NA,
    linewidth = 0.6
  ) +
  geom_point(aes(color = habitat), size = 2, alpha = 0.85) +
  scale_color_manual(values = hab_palette) +
  scale_fill_manual(values = hab_palette) +
  labs(
    title = "Habitat-associated variation after size correction (1950; residual shape space)",
    subtitle = paste0("Allometry: shape ~ ", size_label),
    x = pc1_resid_lab,
    y = pc2_resid_lab,
    color = "Habitat",
    fill  = "Habitat"
  ) +
  coord_equal() +
  theme_classic(base_family = "Helvetica", base_size = 9) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.text  = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8)
  )

print(p_fig3)

# ============================================================
# Save outputs
# ============================================================

pdf_path <- file.path(FIG_DIR, "Fig3_PCA_sizeCorrected_habitat_1950.pdf")
png_path <- file.path(FIG_DIR, "Fig3_PCA_sizeCorrected_habitat_1950.png")

ggsave(pdf_path, p_fig3, width = 7, height = 5, units = "in")
ggsave(png_path, p_fig3, width = 7, height = 5, units = "in", dpi = 600)

cat("Saved:\n  ", normalizePath(pdf_path), "\n  ", normalizePath(png_path), "\n")