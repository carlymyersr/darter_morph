# ============================================================
# scripts/1950_Locality_Resid_PC4.R
# PCA on size-corrected shapes (1950; residuals from procD.lm)
# Plot PC3 vs PC4 by habitat/locality grouping
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
make_hulls <- function(df, group_col, x = "PC3", y = "PC4") {
  group_vals <- unique(df[[group_col]])
  
  hull_list <- lapply(group_vals, function(g) {
    sub <- df[df[[group_col]] == g, , drop = FALSE]
    
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
if (!exists("coords_resid_1950")) stop("coords_resid_1950 not found.")
if (!exists("fit_allo_1950"))     stop("fit_allo_1950 not found.")
stopifnot(identical(dimnames(coords_resid_1950), dimnames(coords_1950)))
stopifnot(identical(dimnames(coords_resid_1950)[[3]], gdf_1950$specimen))

# Label
size_label <- if ("size_label" %in% names(gdf_1950) && !all(is.na(gdf_1950$size_label))) {
  unique(na.omit(gdf_1950$size_label))[1]
} else {
  "size_for_allometry"
}

# ============================================================
# PCA
# ============================================================

pca_resid_1950 <- gm.prcomp(coords_resid_1950)
pct <- 100 * (pca_resid_1950$sdev^2 / sum(pca_resid_1950$sdev^2))

pc3_lab <- sprintf("PC3 (%.1f%%)", pct[3])
pc4_lab <- sprintf("PC4 (%.1f%%)", pct[4])

# ---------------------------
# Data frame
# ---------------------------
pca_df <- data.frame(
  specimen = rownames(pca_resid_1950$x),
  PC3      = pca_resid_1950$x[, 3],
  PC4      = pca_resid_1950$x[, 4],
  habitat  = factor(gdf_1950$habitat),
  stringsAsFactors = FALSE
)

stopifnot(identical(pca_df$specimen, gdf_1950$specimen))

# ---------------------------
# Hulls
# ---------------------------
hull_df <- make_hulls(pca_df, group_col = "habitat", x = "PC3", y = "PC4")

# ============================================================
# Plot
# ============================================================

p_plot <- ggplot(pca_df, aes(PC3, PC4)) +
  geom_polygon(
    data = hull_df,
    aes(group = habitat, fill = habitat),
    alpha = 0.15,
    color = NA
  ) +
  geom_polygon(
    data = hull_df,
    aes(group = habitat, color = habitat),
    fill = NA,
    linewidth = 0.6
  ) +
  geom_point(aes(color = habitat), size = 2, alpha = 0.85) +
  scale_color_manual(values = hab_palette) +
  scale_fill_manual(values = hab_palette) +
  labs(
    title = "Habitat/locality variation after size correction (1950; PC3 vs PC4)",
    subtitle = paste0("Allometry: shape ~ ", size_label),
    x = pc3_lab,
    y = pc4_lab,
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

print(p_plot)

# ============================================================
# Save
# ============================================================

pdf_path <- file.path(FIG_DIR, "Fig_PCA_sizeCorrected_locality_PC3_PC4_1950.pdf")
png_path <- file.path(FIG_DIR, "Fig_PCA_sizeCorrected_locality_PC3_PC4_1950.png")

ggsave(pdf_path, p_plot, width = 7, height = 5, units = "in")
ggsave(png_path, p_plot, width = 7, height = 5, units = "in", dpi = 600)

cat("Saved:\n", normalizePath(pdf_path), "\n", normalizePath(png_path), "\n")