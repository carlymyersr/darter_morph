# ============================================================
# scripts/1950_variation_source_alt_mainstem_vs_smalltrib.R
# PCA on size-corrected shapes (1950; residuals)
# Alternative hypothesis of variation source:
#   1) Mainstem = Connecticut River (1950)
#   2) Quabbin system = Quabbin + Swift River
#   3) Tributaries = Sawmill River + Fort River
# Includes:
#   - procD.lm on shape ~ variation_source_alt
#   - PC1 vs PC2 plot
#   - PC2 vs PC3 plot
#   - PC3 vs PC4 plot
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
make_hulls <- function(df, group_col, x, y) {
  group_vals <- unique(df[[group_col]])
  
  hull_list <- lapply(group_vals, function(g) {
    sub <- df[df[[group_col]] == g, , drop = FALSE]
    
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

OUT_DIR <- file.path("Outputs", "1950_variation_source_alt_quabbinSystem")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

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

# ---------------------------
# Define alternative variation-source grouping
# ---------------------------
gdf_1950 <- gdf_1950 %>%
  mutate(
    variation_source_alt = case_when(
      habitat == "Connecticut River" ~ "Mainstem",
      habitat %in% c("Quabbin", "Swift River") ~ "Quabbin system",
      habitat %in% c("Sawmill River", "Fort River") ~ "Tributaries",
      TRUE ~ NA_character_
    ),
    variation_source_alt = factor(
      variation_source_alt,
      levels = c("Mainstem", "Quabbin system", "Tributaries")
    )
  )

# Drop unexpected NA rows
keep <- !is.na(gdf_1950$variation_source_alt)
gdf_vs <- gdf_1950[keep, , drop = FALSE]
coords_vs <- coords_resid_1950[, , keep, drop = FALSE]

stopifnot(identical(dimnames(coords_vs)[[3]], gdf_vs$specimen))

cat("\nHabitat x variation_source_alt table:\n")
print(table(gdf_vs$habitat, gdf_vs$variation_source_alt))
cat("\nVariation source alt counts:\n")
print(table(gdf_vs$variation_source_alt))

# ---------------------------
# Custom colors
# ---------------------------
vs_palette <- c(
  "Mainstem" = "#1B9E77",
  "Quabbin system" = "#E6AB02",
  "Tributaries" = "#7570B3"
)

# ============================================================
# procD: test alternative variation-source grouping
# ============================================================

procD_formula <- coords_vs ~ variation_source_alt

fit_procD_alt <- procD.lm(
  f1 = procD_formula,
  data = gdf_vs,
  iter = 999,
  RRPP = TRUE
)

cat("\n============================================================\n")
cat("procD.lm results: shape ~ variation_source_alt\n")
cat("============================================================\n")
print(summary(fit_procD_alt))

capture.output(
  {
    cat("============================================================\n")
    cat("procD.lm results: shape ~ variation_source_alt\n")
    cat("============================================================\n")
    print(summary(fit_procD_alt))
  },
  file = file.path(OUT_DIR, "procD_variation_source_alt_summary.txt")
)

saveRDS(fit_procD_alt, file = file.path(OUT_DIR, "fit_procD_variation_source_alt.rds"))

# ============================================================
# PCA on residuals for alternative grouping subset
# ============================================================

pca_vs <- gm.prcomp(coords_vs)
pct_vs <- 100 * (pca_vs$sdev^2 / sum(pca_vs$sdev^2))

pc1_lab <- sprintf("PC1 (%.1f%%)", pct_vs[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", pct_vs[2])
pc3_lab <- sprintf("PC3 (%.1f%%)", pct_vs[3])
pc4_lab <- sprintf("PC4 (%.1f%%)", pct_vs[4])

pca_vs_df <- data.frame(
  specimen             = rownames(pca_vs$x),
  PC1                  = pca_vs$x[, 1],
  PC2                  = pca_vs$x[, 2],
  PC3                  = pca_vs$x[, 3],
  PC4                  = pca_vs$x[, 4],
  variation_source_alt = gdf_vs$variation_source_alt,
  habitat              = gdf_vs$habitat,
  stringsAsFactors     = FALSE
)

stopifnot(identical(pca_vs_df$specimen, gdf_vs$specimen))

# Hulls
hull_pc12 <- make_hulls(pca_vs_df, group_col = "variation_source_alt", x = "PC1", y = "PC2")
hull_pc23 <- make_hulls(pca_vs_df, group_col = "variation_source_alt", x = "PC2", y = "PC3")
hull_pc34 <- make_hulls(pca_vs_df, group_col = "variation_source_alt", x = "PC3", y = "PC4")

# ============================================================
# Plot 1: PC1 vs PC2
# ============================================================

p_pc12 <- ggplot(pca_vs_df, aes(PC1, PC2)) +
  geom_polygon(
    data = hull_pc12,
    aes(group = variation_source_alt, fill = variation_source_alt),
    alpha = 0.15,
    color = NA
  ) +
  geom_polygon(
    data = hull_pc12,
    aes(group = variation_source_alt, color = variation_source_alt),
    fill = NA,
    linewidth = 0.6
  ) +
  geom_point(aes(color = variation_source_alt), size = 2, alpha = 0.85) +
  scale_color_manual(values = vs_palette) +
  scale_fill_manual(values = vs_palette) +
  labs(
    title = "Alternative variation-source hypothesis in 1950 after size correction",
    subtitle = paste0("PC1 vs PC2 | Allometry: shape ~ ", size_label),
    x = pc1_lab,
    y = pc2_lab,
    color = "Variation source",
    fill  = "Variation source"
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

print(p_pc12)

# ============================================================
# Plot 2: PC2 vs PC3
# ============================================================

p_pc23 <- ggplot(pca_vs_df, aes(PC2, PC3)) +
  geom_polygon(
    data = hull_pc23,
    aes(group = variation_source_alt, fill = variation_source_alt),
    alpha = 0.15,
    color = NA
  ) +
  geom_polygon(
    data = hull_pc23,
    aes(group = variation_source_alt, color = variation_source_alt),
    fill = NA,
    linewidth = 0.6
  ) +
  geom_point(aes(color = variation_source_alt), size = 2, alpha = 0.85) +
  scale_color_manual(values = vs_palette) +
  scale_fill_manual(values = vs_palette) +
  labs(
    title = "Alternative variation-source hypothesis in 1950 after size correction",
    subtitle = paste0("PC2 vs PC3 | Allometry: shape ~ ", size_label),
    x = pc2_lab,
    y = pc3_lab,
    color = "Variation source",
    fill  = "Variation source"
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

print(p_pc23)

# ============================================================
# Plot 3: PC3 vs PC4
# ============================================================

p_pc34 <- ggplot(pca_vs_df, aes(PC3, PC4)) +
  geom_polygon(
    data = hull_pc34,
    aes(group = variation_source_alt, fill = variation_source_alt),
    alpha = 0.15,
    color = NA
  ) +
  geom_polygon(
    data = hull_pc34,
    aes(group = variation_source_alt, color = variation_source_alt),
    fill = NA,
    linewidth = 0.6
  ) +
  geom_point(aes(color = variation_source_alt), size = 2, alpha = 0.85) +
  scale_color_manual(values = vs_palette) +
  scale_fill_manual(values = vs_palette) +
  labs(
    title = "Alternative variation-source hypothesis in 1950 after size correction",
    subtitle = paste0("PC3 vs PC4 | Allometry: shape ~ ", size_label),
    x = pc3_lab,
    y = pc4_lab,
    color = "Variation source",
    fill  = "Variation source"
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

print(p_pc34)

# ============================================================
# Save outputs
# ============================================================

pdf_path_pc12 <- file.path(FIG_DIR, "Fig_PCA_variationSourceAlt_1950_PC12.pdf")
png_path_pc12 <- file.path(FIG_DIR, "Fig_PCA_variationSourceAlt_1950_PC12.png")

pdf_path_pc23 <- file.path(FIG_DIR, "Fig_PCA_variationSourceAlt_1950_PC23.pdf")
png_path_pc23 <- file.path(FIG_DIR, "Fig_PCA_variationSourceAlt_1950_PC23.png")

pdf_path_pc34 <- file.path(FIG_DIR, "Fig_PCA_variationSourceAlt_1950_PC34.pdf")
png_path_pc34 <- file.path(FIG_DIR, "Fig_PCA_variationSourceAlt_1950_PC34.png")

ggsave(pdf_path_pc12, p_pc12, width = 7, height = 5, units = "in")
ggsave(png_path_pc12, p_pc12, width = 7, height = 5, units = "in", dpi = 600)

ggsave(pdf_path_pc23, p_pc23, width = 7, height = 5, units = "in")
ggsave(png_path_pc23, p_pc23, width = 7, height = 5, units = "in", dpi = 600)

ggsave(pdf_path_pc34, p_pc34, width = 7, height = 5, units = "in")
ggsave(png_path_pc34, p_pc34, width = 7, height = 5, units = "in", dpi = 600)

cat("\nSaved figures:\n")
cat("  ", normalizePath(pdf_path_pc12), "\n")
cat("  ", normalizePath(png_path_pc12), "\n")
cat("  ", normalizePath(pdf_path_pc23), "\n")
cat("  ", normalizePath(png_path_pc23), "\n")
cat("  ", normalizePath(pdf_path_pc34), "\n")
cat("  ", normalizePath(png_path_pc34), "\n")
cat("\nSaved procD outputs:\n")
cat("  ", normalizePath(file.path(OUT_DIR, "procD_variation_source_alt_summary.txt")), "\n")
cat("  ", normalizePath(file.path(OUT_DIR, "fit_procD_variation_source_alt.rds")), "\n")