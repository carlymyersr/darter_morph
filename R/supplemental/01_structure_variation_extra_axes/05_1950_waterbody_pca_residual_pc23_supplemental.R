# ============================================================
# scripts/1950_Locality_Resid_PC3.R
# PCA on size-corrected shapes (1950; residuals from procD.lm)
# Plot PC2 vs PC3 by habitat/locality grouping
# Standalone script modeled after 1950_PCA_residual_Fig3.R
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
})

# ---------------------------
# Helper: convex hulls by group
# ---------------------------
make_hulls <- function(df, group_col, x = "PC2", y = "PC3") {
  group_vals <- unique(df[[group_col]])

  hull_list <- lapply(group_vals, function(g) {
    sub <- df[df[[group_col]] == g, , drop = FALSE]

    # Need at least 3 non-collinear points for a hull
    if (nrow(sub) < 3) return(NULL)

    h <- chull(sub[[x]], sub[[y]])
    sub[h, , drop = FALSE]
  })

  hull_df <- do.call(rbind, hull_list)
  rownames(hull_df) <- NULL
  hull_df
}

# ---------------------------
# Output directory
# ---------------------------
FIG_DIR <- "figures"
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

# ---------------------------
# Preconditions
# ---------------------------
if (!exists("coords_resid_1950")) stop("coords_resid_1950 not found. Did R/methods/01_specimen_sampling_study_design/02_subset_1950.R run fully?")
if (!exists("fit_allo_1950"))     stop("fit_allo_1950 not found. Did R/methods/01_specimen_sampling_study_design/02_subset_1950.R run fully?")
stopifnot(identical(dimnames(coords_resid_1950), dimnames(coords_1950)))
stopifnot(identical(dimnames(coords_resid_1950)[[3]], gdf_1950$specimen))

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

pc2_resid_lab <- sprintf("PC2 (%.1f%%)", pct_resid_1950[2])
pc3_resid_lab <- sprintf("PC3 (%.1f%%)", pct_resid_1950[3])

# --- Canonical plotting df ---
pca1950_resid_df <- data.frame(
  specimen = rownames(pca_resid_1950$x),
  PC2      = pca_resid_1950$x[, 2],
  PC3      = pca_resid_1950$x[, 3],
  habitat  = factor(gdf_1950$habitat),
  stringsAsFactors = FALSE
)
stopifnot(identical(pca1950_resid_df$specimen, gdf_1950$specimen))

# --- Hulls by habitat ---
hull1950_resid_df <- make_hulls(pca1950_resid_df, group_col = "habitat", x = "PC2", y = "PC3")

# ============================================================
# Plot
# ============================================================

p_locality_pc23 <- ggplot(pca1950_resid_df, aes(PC2, PC3)) +
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
    title = "Habitat/locality variation after size correction (1950; PC2 vs PC3)",
    subtitle = paste0("Allometry: shape ~ ", size_label),
    x = pc2_resid_lab,
    y = pc3_resid_lab,
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

print(p_locality_pc23)

# ============================================================
# Save outputs
# ============================================================

pdf_path <- file.path(FIG_DIR, "Fig_PCA_sizeCorrected_locality_PC2_PC3_1950.pdf")
png_path <- file.path(FIG_DIR, "Fig_PCA_sizeCorrected_locality_PC2_PC3_1950.png")

ggsave(pdf_path, p_locality_pc23, width = 7, height = 5, units = "in")
ggsave(png_path, p_locality_pc23, width = 7, height = 5, units = "in", dpi = 600)

cat("Saved:\n  ", normalizePath(pdf_path), "\n  ", normalizePath(png_path), "\n")
