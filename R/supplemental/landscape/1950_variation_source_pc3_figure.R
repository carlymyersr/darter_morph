# ============================================================
# scripts/1950_variation_source_pc3_figure.R
# PCA on size-corrected shapes (1950; residuals)
# Variation source grouping:
#   1) Connecticut River
#   2) Quabbin
#   3) Small tributaries (Sawmill, Swift, Fort)
# Standalone script modeled after 1950_variation_source_figure.R
# PC2 vs PC3 version
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
# Define variation source grouping
# ---------------------------
gdf_1950 <- gdf_1950 %>%
  mutate(
    variation_source = case_when(
      habitat == "Connecticut River" ~ "Mainstem",
      habitat == "Quabbin"           ~ "Reservoir",
      habitat %in% c("Sawmill River", "Swift River", "Fort River") ~ "Small tributaries",
      TRUE ~ NA_character_
    ),
    variation_source = factor(
      variation_source,
      levels = c("Mainstem", "Reservoir", "Small tributaries")
    )
  )

# Drop unexpected NA rows
keep <- !is.na(gdf_1950$variation_source)
gdf_vs <- gdf_1950[keep, , drop = FALSE]
coords_vs <- coords_resid_1950[, , keep, drop = FALSE]

stopifnot(identical(dimnames(coords_vs)[[3]], gdf_vs$specimen))

# Inspect grouping
cat("\nHabitat x variation_source table:\n")
print(table(gdf_vs$habitat, gdf_vs$variation_source))
cat("\nVariation source counts:\n")
print(table(gdf_vs$variation_source))

# ---------------------------
# Custom colors
# ---------------------------
vs_palette <- c(
  "Mainstem"           = "#1B9E77",  # teal-green
  "Reservoir"          = "#F4A261",  # lighter orange
  "Small tributaries"  = "#7570B3"   # purple
)

# ============================================================
# PCA on residuals for variation-source subset
# ============================================================

pca_vs <- gm.prcomp(coords_vs)
pct_vs <- 100 * (pca_vs$sdev^2 / sum(pca_vs$sdev^2))

if (length(pct_vs) < 3) stop("PCA returned fewer than 3 components; cannot plot PC2 vs PC3.")

pc2_lab <- sprintf("PC2 (%.1f%%)", pct_vs[2])
pc3_lab <- sprintf("PC3 (%.1f%%)", pct_vs[3])

# Plotting df
pca_vs_df <- data.frame(
  specimen         = rownames(pca_vs$x),
  PC1              = pca_vs$x[, 1],
  PC2              = pca_vs$x[, 2],
  PC3              = pca_vs$x[, 3],
  variation_source = gdf_vs$variation_source,
  habitat          = gdf_vs$habitat,
  stringsAsFactors = FALSE
)

stopifnot(identical(pca_vs_df$specimen, gdf_vs$specimen))

# Hulls by variation-source group
hull_vs_df <- make_hulls(pca_vs_df, group_col = "variation_source", x = "PC2", y = "PC3")

# ============================================================
# Plot
# ============================================================

p_variation_source_pc23 <- ggplot(pca_vs_df, aes(PC2, PC3)) +
  geom_polygon(
    data = hull_vs_df,
    aes(group = variation_source, fill = variation_source),
    alpha = 0.15,
    color = NA
  ) +
  geom_polygon(
    data = hull_vs_df,
    aes(group = variation_source, color = variation_source),
    fill = NA,
    linewidth = 0.6
  ) +
  geom_point(aes(color = variation_source), size = 2, alpha = 0.85) +
  scale_color_manual(values = vs_palette) +
  scale_fill_manual(values = vs_palette) +
  labs(
    title = "Variation sources in 1950 after size correction",
    subtitle = paste0("Allometry: shape ~ ", size_label),
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

print(p_variation_source_pc23)

# ============================================================
# Save outputs
# ============================================================

pdf_path <- file.path(FIG_DIR, "Fig_PCA_variationSource_PC2_PC3_1950.pdf")
png_path <- file.path(FIG_DIR, "Fig_PCA_variationSource_PC2_PC3_1950.png")

ggsave(pdf_path, p_variation_source_pc23, width = 7, height = 5, units = "in")
ggsave(png_path, p_variation_source_pc23, width = 7, height = 5, units = "in", dpi = 600)

cat("\nSaved:\n  ", normalizePath(pdf_path), "\n  ", normalizePath(png_path), "\n")
