# ============================================================
# scripts/1950_bigsmall_figure.R
# PCA on size-corrected shapes (1950; residuals)
# Big water vs small water grouping
# Formatting matched to 1950_PCA_residual_Fig3.R
# Includes:
#   - PC1 vs PC2
#   - PC2 vs PC3
#   - PC3 vs PC4
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

# Make hulls for figure
make_hulls <- function(df, group_col, x = "PC1", y = "PC2") {
  group_vals <- unique(df[[group_col]])
  
  hull_list <- lapply(group_vals, function(g) {
    sub <- df[df[[group_col]] == g, , drop = FALSE]
    
    # Need at least 3 non-collinear points for a hull
    if (nrow(sub) < 3) return(NULL)
    
    h <- chull(sub[[x]], sub[[y]])
    sub[h, , drop = FALSE]
  })
  
  hull_df <- do.call(rbind, hull_list)
  if (!is.null(hull_df)) rownames(hull_df) <- NULL
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
if (!exists("coords_resid_1950")) stop("coords_resid_1950 not found. Did R/02_subset_1950.R run fully?")
if (!exists("fit_allo_1950"))     stop("fit_allo_1950 not found. Did R/02_subset_1950.R run fully?")
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

# ---------------------------
# Define big water / small water grouping
# ---------------------------
gdf_1950 <- gdf_1950 %>%
  mutate(
    water_size = case_when(
      habitat %in% c("Connecticut River", "Quabbin") ~ "Big water",
      habitat %in% c("Sawmill River", "Swift River", "Fort River") ~ "Small water",
      TRUE ~ NA_character_
    ),
    water_size = factor(water_size, levels = c("Big water", "Small water"))
  )

# Drop unexpected NA rows
keep <- !is.na(gdf_1950$water_size)
gdf_bw <- gdf_1950[keep, , drop = FALSE]
coords_bw <- coords_resid_1950[, , keep, drop = FALSE]

stopifnot(identical(dimnames(coords_bw)[[3]], gdf_bw$specimen))

# Inspect grouping
print(table(gdf_bw$habitat, gdf_bw$water_size))
print(table(gdf_bw$water_size))

# Optional custom colors for the binary grouping
bw_palette <- c(
  "Big water"   = "#1B9E77",  # teal-green (colorblind safe)
  "Small water" = "#F4A261"   # warm orange
)

# ============================================================
# PCA on residuals for big/small subset
# ============================================================

pca_bw <- gm.prcomp(coords_bw)
pct_bw <- 100 * (pca_bw$sdev^2 / sum(pca_bw$sdev^2))

pc1_lab <- sprintf("PC1 (%.1f%%)", pct_bw[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", pct_bw[2])
pc3_lab <- sprintf("PC3 (%.1f%%)", pct_bw[3])
pc4_lab <- sprintf("PC4 (%.1f%%)", pct_bw[4])

# Plotting df
pca_bw_df <- data.frame(
  specimen   = rownames(pca_bw$x),
  PC1        = pca_bw$x[, 1],
  PC2        = pca_bw$x[, 2],
  PC3        = pca_bw$x[, 3],
  PC4        = pca_bw$x[, 4],
  water_size = gdf_bw$water_size,
  habitat    = gdf_bw$habitat,
  stringsAsFactors = FALSE
)

stopifnot(identical(pca_bw_df$specimen, gdf_bw$specimen))

# Hulls by big/small group
hull_bw_pc12 <- make_hulls(pca_bw_df, group_col = "water_size", x = "PC1", y = "PC2")
hull_bw_pc23 <- make_hulls(pca_bw_df, group_col = "water_size", x = "PC2", y = "PC3")
hull_bw_pc34 <- make_hulls(pca_bw_df, group_col = "water_size", x = "PC3", y = "PC4")

# ============================================================
# Plot 1: PC1 vs PC2
# ============================================================

p_bigsmall_pc12 <- ggplot(pca_bw_df, aes(PC1, PC2)) +
  geom_polygon(
    data = hull_bw_pc12,
    aes(group = water_size, fill = water_size),
    alpha = 0.15,
    color = NA
  ) +
  geom_polygon(
    data = hull_bw_pc12,
    aes(group = water_size, color = water_size),
    fill = NA,
    linewidth = 0.6
  ) +
  geom_point(aes(color = water_size), size = 2, alpha = 0.85) +
  scale_color_manual(values = bw_palette) +
  scale_fill_manual(values = bw_palette) +
  labs(
    title = "Big-water vs small-water variation after size correction (1950; residual shape space)",
    subtitle = paste0("PC1 vs PC2 | Allometry: shape ~ ", size_label),
    x = pc1_lab,
    y = pc2_lab,
    color = "Water type",
    fill  = "Water type"
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

print(p_bigsmall_pc12)

# ============================================================
# Plot 2: PC2 vs PC3
# ============================================================

p_bigsmall_pc23 <- ggplot(pca_bw_df, aes(PC2, PC3)) +
  geom_polygon(
    data = hull_bw_pc23,
    aes(group = water_size, fill = water_size),
    alpha = 0.15,
    color = NA
  ) +
  geom_polygon(
    data = hull_bw_pc23,
    aes(group = water_size, color = water_size),
    fill = NA,
    linewidth = 0.6
  ) +
  geom_point(aes(color = water_size), size = 2, alpha = 0.85) +
  scale_color_manual(values = bw_palette) +
  scale_fill_manual(values = bw_palette) +
  labs(
    title = "Big-water vs small-water variation after size correction (1950; residual shape space)",
    subtitle = paste0("PC2 vs PC3 | Allometry: shape ~ ", size_label),
    x = pc2_lab,
    y = pc3_lab,
    color = "Water type",
    fill  = "Water type"
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

print(p_bigsmall_pc23)

# ============================================================
# Plot 3: PC3 vs PC4
# ============================================================

p_bigsmall_pc34 <- ggplot(pca_bw_df, aes(PC3, PC4)) +
  geom_polygon(
    data = hull_bw_pc34,
    aes(group = water_size, fill = water_size),
    alpha = 0.15,
    color = NA
  ) +
  geom_polygon(
    data = hull_bw_pc34,
    aes(group = water_size, color = water_size),
    fill = NA,
    linewidth = 0.6
  ) +
  geom_point(aes(color = water_size), size = 2, alpha = 0.85) +
  scale_color_manual(values = bw_palette) +
  scale_fill_manual(values = bw_palette) +
  labs(
    title = "Big-water vs small-water variation after size correction (1950; residual shape space)",
    subtitle = paste0("PC3 vs PC4 | Allometry: shape ~ ", size_label),
    x = pc3_lab,
    y = pc4_lab,
    color = "Water type",
    fill  = "Water type"
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

print(p_bigsmall_pc34)

# ============================================================
# Save outputs
# ============================================================

pdf_path_pc12 <- file.path(FIG_DIR, "Fig_PCA_sizeCorrected_bigSmallWater_1950_PC12.pdf")
png_path_pc12 <- file.path(FIG_DIR, "Fig_PCA_sizeCorrected_bigSmallWater_1950_PC12.png")

pdf_path_pc23 <- file.path(FIG_DIR, "Fig_PCA_sizeCorrected_bigSmallWater_1950_PC23.pdf")
png_path_pc23 <- file.path(FIG_DIR, "Fig_PCA_sizeCorrected_bigSmallWater_1950_PC23.png")

pdf_path_pc34 <- file.path(FIG_DIR, "Fig_PCA_sizeCorrected_bigSmallWater_1950_PC34.pdf")
png_path_pc34 <- file.path(FIG_DIR, "Fig_PCA_sizeCorrected_bigSmallWater_1950_PC34.png")

ggsave(pdf_path_pc12, p_bigsmall_pc12, width = 7, height = 5, units = "in")
ggsave(png_path_pc12, p_bigsmall_pc12, width = 7, height = 5, units = "in", dpi = 600)

ggsave(pdf_path_pc23, p_bigsmall_pc23, width = 7, height = 5, units = "in")
ggsave(png_path_pc23, p_bigsmall_pc23, width = 7, height = 5, units = "in", dpi = 600)

ggsave(pdf_path_pc34, p_bigsmall_pc34, width = 7, height = 5, units = "in")
ggsave(png_path_pc34, p_bigsmall_pc34, width = 7, height = 5, units = "in", dpi = 600)

cat("Saved:\n")
cat("  ", normalizePath(pdf_path_pc12), "\n")
cat("  ", normalizePath(png_path_pc12), "\n")
cat("  ", normalizePath(pdf_path_pc23), "\n")
cat("  ", normalizePath(png_path_pc23), "\n")
cat("  ", normalizePath(pdf_path_pc34), "\n")
cat("  ", normalizePath(png_path_pc34), "\n")