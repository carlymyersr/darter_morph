# ============================================================
# scripts/mainsteam_alltime_vs_tribandres.R
# PCA on size-corrected shapes (residuals)
# Variation-source grouping:
#   1) Mainstem    = Connecticut River (1950, 1956, 1970)
#   2) Reservoir System   = Quabbin + Swift River
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
source("R/00_setup_morpho.R")
source("R/01_build_metadata.R")

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
# Output directories
# ---------------------------
FIG_DIR <- "figures"
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

OUT_DIR <- file.path("Outputs", "mainsteam_alltime_vs_tribandres")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ---------------------------
# Preconditions
# ---------------------------
if (!exists("coords_gpa")) stop("coords_gpa not found. Did R/00_setup_morpho.R run fully?")
if (!exists("gdf"))       stop("gdf not found. Did R/01_build_metadata.R run fully?")
if (!exists("subset_coords_to_gdf")) {
  stop("subset_coords_to_gdf() not found. It should be defined in R/00_setup_morpho.R")
}
if (!exists("allometry_residuals")) {
  stop("allometry_residuals() not found. It should be defined in R/01_build_metadata.R")
}

# ---------------------------
# Build combined subset
# ---------------------------
# Mainstem includes Connecticut River across all three sampled years.
# Reservoir groups are retained as 1950 habitats, matching the original reference structure.

gdf_vs <- gdf %>%
  filter(
    !is.na(habitat),
    (
      habitat == "Connecticut River" & year %in% c(1950, 1956, 1970)
    ) |
      habitat %in% c("Quabbin", "Swift River", "Sawmill River", "Fort River")
  ) %>%
  mutate(
    variation_source_alt = case_when(
      habitat == "Connecticut River" ~ "Mainstem",
      habitat %in% c("Quabbin", "Swift River") ~ "Reservoir System",
      habitat %in% c("Sawmill River", "Fort River") ~ "Tributaries",
      TRUE ~ NA_character_
    ),
    variation_source_alt = factor(
      variation_source_alt,
      levels = c("Mainstem", "Reservoir System", "Tributaries")
    )
  ) %>%
  filter(!is.na(variation_source_alt)) %>%
  droplevels()

stopifnot(nrow(gdf_vs) >= 3)

cat("\nHabitat x year table:\n")
print(table(gdf_vs$habitat, gdf_vs$year))

cat("\nHabitat x variation_source_alt table:\n")
print(table(gdf_vs$habitat, gdf_vs$variation_source_alt))

cat("\nVariation source alt counts:\n")
print(table(gdf_vs$variation_source_alt))

# ---------------------------
# Subset coordinates to match metadata order
# ---------------------------
coords_vs_raw <- subset_coords_to_gdf(coords_gpa, gdf_vs)
stopifnot(identical(dimnames(coords_vs_raw)[[3]], gdf_vs$specimen))

# ---------------------------
# Allometry correction on this combined subset
# ---------------------------
size_vec <- setNames(gdf_vs$size_for_allometry, gdf_vs$specimen)
allo_vs <- allometry_residuals(coords_vs_raw, size_vec)
coords_vs <- allo_vs$residuals

stopifnot(identical(dimnames(coords_vs)[[3]], gdf_vs$specimen))

size_label <- if ("size_label" %in% names(gdf_vs) && !all(is.na(gdf_vs$size_label))) {
  unique(na.omit(gdf_vs$size_label))[1]
} else if ("size_label" %in% names(gdf) && !all(is.na(gdf$size_label))) {
  unique(na.omit(gdf$size_label))[1]
} else {
  "size_for_allometry"
}

# ---------------------------
# Custom colors
# ---------------------------
vs_palette <- c(
  "Mainstem" = "#1B9E77",
  "Reservoir System" = "#E6AB02",
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
saveRDS(allo_vs,        file = file.path(OUT_DIR, "allometry_fit_combined_subset.rds"))

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
  year                 = gdf_vs$year,
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
    title = "Variation-source grouping with Connecticut River across time after size correction",
    subtitle = paste0("PC1 vs PC2 | Mainstem = CT 1950 + 1956 + 1970 | Allometry: shape ~ ", size_label),
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
    title = "Variation-source grouping with Connecticut River across time after size correction",
    subtitle = paste0("PC2 vs PC3 | Mainstem = CT 1950 + 1956 + 1970 | Allometry: shape ~ ", size_label),
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
    title = "Variation-source grouping with Connecticut River across time after size correction",
    subtitle = paste0("PC3 vs PC4 | Mainstem = CT 1950 + 1956 + 1970 | Allometry: shape ~ ", size_label),
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
# Mahalanobis distance test:
# Are individuals outside the Mainstem-defined landscape?
# Reference = Mainstem in PC1-PC4 space
# Threshold = empirical 95th percentile of observed Mainstem distances
# ============================================================

library(dplyr)
library(ggplot2)

MAHAL_PCS <- c("PC1", "PC2", "PC3", "PC4")

# Mainstem reference scores
mainstem_ref <- pca_vs_df %>%
  filter(variation_source_alt == "Mainstem") %>%
  select(all_of(MAHAL_PCS))

# Safety checks
if (nrow(mainstem_ref) <= length(MAHAL_PCS)) {
  stop("Not enough Mainstem specimens for stable covariance estimation with PC1-PC4.")
}

# Mainstem mean and covariance matrix
mainstem_center <- colMeans(mainstem_ref)
mainstem_cov <- cov(mainstem_ref)

# Optional but helpful: check whether covariance matrix is invertible
if (det(mainstem_cov) < .Machine$double.eps) {
  stop("Mainstem covariance matrix is nearly singular. Try fewer PCs, e.g. PC1-PC3.")
}

# Calculate Mahalanobis distance for every individual
pca_vs_df <- pca_vs_df %>%
  mutate(
    mahal_d2_mainstem_PC1_PC4 = mahalanobis(
      x = select(., all_of(MAHAL_PCS)),
      center = mainstem_center,
      cov = mainstem_cov
    ),
    mahal_d_mainstem_PC1_PC4 = sqrt(mahal_d2_mainstem_PC1_PC4)
  )

# Empirical 95th percentile threshold based only on observed Mainstem distances
mainstem_95_threshold_d2 <- quantile(
  pca_vs_df$mahal_d2_mainstem_PC1_PC4[pca_vs_df$variation_source_alt == "Mainstem"],
  probs = 0.95,
  na.rm = TRUE
)

mainstem_95_threshold_d <- sqrt(mainstem_95_threshold_d2)

# Flag individuals outside the observed Mainstem reference distribution
pca_vs_df <- pca_vs_df %>%
  mutate(
    outside_mainstem_95_PC1_PC4 =
      mahal_d2_mainstem_PC1_PC4 > mainstem_95_threshold_d2
  )

# Summary by group
mahal_summary_by_group <- pca_vs_df %>%
  group_by(variation_source_alt) %>%
  summarise(
    n = n(),
    n_outside_mainstem_95 = sum(outside_mainstem_95_PC1_PC4, na.rm = TRUE),
    prop_outside_mainstem_95 = n_outside_mainstem_95 / n,
    mean_mahal_d = mean(mahal_d_mainstem_PC1_PC4, na.rm = TRUE),
    median_mahal_d = median(mahal_d_mainstem_PC1_PC4, na.rm = TRUE),
    max_mahal_d = max(mahal_d_mainstem_PC1_PC4, na.rm = TRUE),
    .groups = "drop"
  )

print(mahal_summary_by_group)

# Individual-level output
mahal_individuals <- pca_vs_df %>%
  select(
    specimen,
    habitat,
    year,
    variation_source_alt,
    all_of(MAHAL_PCS),
    mahal_d2_mainstem_PC1_PC4,
    mahal_d_mainstem_PC1_PC4,
    outside_mainstem_95_PC1_PC4
  ) %>%
  arrange(desc(mahal_d_mainstem_PC1_PC4))

print(mahal_individuals)

# Save outputs
write.csv(
  mahal_individuals,
  file = file.path(OUT_DIR, "mahalanobis_individuals_PC1_PC4_mainstem_reference.csv"),
  row.names = FALSE
)

write.csv(
  mahal_summary_by_group,
  file = file.path(OUT_DIR, "mahalanobis_summary_by_group_PC1_PC4_mainstem_reference.csv"),
  row.names = FALSE
)

capture.output(
  {
    cat("============================================================\n")
    cat("Mahalanobis distance test: Mainstem reference distribution\n")
    cat("PCs used:", paste(MAHAL_PCS, collapse = ", "), "\n")
    cat("Reference group: Mainstem = Connecticut River 1950, 1956, 1970\n")
    cat("Threshold: empirical 95th percentile of observed Mainstem distances\n")
    cat("95th percentile squared distance:", mainstem_95_threshold_d2, "\n")
    cat("95th percentile distance:", mainstem_95_threshold_d, "\n\n")
    
    cat("Summary by group:\n")
    print(mahal_summary_by_group)
    
    cat("\nIndividuals ranked by Mahalanobis distance:\n")
    print(mahal_individuals)
  },
  file = file.path(OUT_DIR, "mahalanobis_PC1_PC4_mainstem_reference_summary.txt")
)

cat("\nMahalanobis outputs saved to:\n")
cat("  ", normalizePath(file.path(OUT_DIR, "mahalanobis_individuals_PC1_PC4_mainstem_reference.csv")), "\n")
cat("  ", normalizePath(file.path(OUT_DIR, "mahalanobis_summary_by_group_PC1_PC4_mainstem_reference.csv")), "\n")
cat("  ", normalizePath(file.path(OUT_DIR, "mahalanobis_PC1_PC4_mainstem_reference_summary.txt")), "\n")

# ============================================================
# Save outputs
# ============================================================

pdf_path_pc12 <- file.path(FIG_DIR, "Fig_mainsteam_alltime_vs_tribandres_PC12.pdf")
png_path_pc12 <- file.path(FIG_DIR, "Fig_mainsteam_alltime_vs_tribandres_PC12.png")

pdf_path_pc23 <- file.path(FIG_DIR, "Fig_mainsteam_alltime_vs_tribandres_PC23.pdf")
png_path_pc23 <- file.path(FIG_DIR, "Fig_mainsteam_alltime_vs_tribandres_PC23.png")

pdf_path_pc34 <- file.path(FIG_DIR, "Fig_mainsteam_alltime_vs_tribandres_PC34.pdf")
png_path_pc34 <- file.path(FIG_DIR, "Fig_mainsteam_alltime_vs_tribandres_PC34.png")

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
cat("\nSaved procD/allometry outputs:\n")
cat("  ", normalizePath(file.path(OUT_DIR, "procD_variation_source_alt_summary.txt")), "\n")
cat("  ", normalizePath(file.path(OUT_DIR, "fit_procD_variation_source_alt.rds")), "\n")
cat("  ", normalizePath(file.path(OUT_DIR, "allometry_fit_combined_subset.rds")), "\n")
