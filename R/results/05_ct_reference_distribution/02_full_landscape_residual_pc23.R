# ============================================================
# Scripts/Figure7_RESID_CTtimeseries_plus_1950habitats_PC23.R
# PCA of residual coordinates (PC2 vs PC3)
# CT time series + 1950 habitats
# Saves to: Figures/
# ============================================================

suppressPackageStartupMessages({
  library(geomorph)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

# Canonical upstream
source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")
source("R/methods/01_specimen_sampling_study_design/04_subset_CT_timeseries_plus_1950habitats.R")  # provides gdf_Fig6, coords_resid_Fig6, allo_Fig6

fig_dir <- "Figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helper: convex hulls
# ============================================================

make_hulls <- function(df, group_col, x = "PC2", y = "PC3", min_n = 3) {
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
# Palette / labels
# ============================================================

group_levels <- levels(gdf_Fig6$group)
group_palette <- setNames(rep(NA_character_, length(group_levels)), group_levels)
group_palette["CT_1950"] <- "steelblue"
group_palette["CT_1956"] <- "dodgerblue4"
group_palette["CT_1970"] <- "navy"

other <- setdiff(group_levels, c("CT_1950","CT_1956","CT_1970"))
if (length(other) > 0) {
  if (exists("hab_palette", inherits = TRUE) &&
      is.vector(hab_palette) &&
      !is.null(names(hab_palette))) {
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

group_labels <- setNames(group_levels, group_levels)
group_labels["CT_1950"] <- "Connecticut River (1950)"
group_labels["CT_1956"] <- "Connecticut River (1956)"
group_labels["CT_1970"] <- "Connecticut River (1970)"
for (g in other) group_labels[g] <- gsub("_1950$", " (1950)", g)

# ============================================================
# PCA on residual coordinates
# ============================================================

pca_Fig7 <- geomorph::gm.prcomp(coords_resid_Fig6)
pct_Fig7 <- 100 * (pca_Fig7$sdev^2 / sum(pca_Fig7$sdev^2))

pc2_lab <- sprintf("PC2 (%.1f%%)", pct_Fig7[2])
pc3_lab <- sprintf("PC3 (%.1f%%)", pct_Fig7[3])

pca_Fig7_df <- data.frame(
  specimen = rownames(pca_Fig7$x),
  PC2      = pca_Fig7$x[, 2],
  PC3      = pca_Fig7$x[, 3],
  group    = gdf_Fig6$group,
  habitat  = gdf_Fig6$habitat,
  year     = gdf_Fig6$year,
  SL_mm    = gdf_Fig6$SL_mm,
  logSL    = gdf_Fig6$logSL,
  stringsAsFactors = FALSE
)

stopifnot(identical(pca_Fig7_df$specimen, gdf_Fig6$specimen))

hull_Fig7_df <- make_hulls(
  df = pca_Fig7_df,
  group_col = "group",
  x = "PC2",
  y = "PC3",
  min_n = 3
)

# ============================================================
# Plot: PC2 vs PC3
# ============================================================

p_fig7_pc23 <- ggplot(pca_Fig7_df, aes(PC2, PC3)) +
  geom_polygon(
    data = hull_Fig7_df,
    aes(group = group, fill = group),
    alpha = 0.15,
    color = NA
  ) +
  geom_polygon(
    data = hull_Fig7_df,
    aes(group = group, color = group),
    fill = NA,
    linewidth = 0.7
  ) +
  geom_point(aes(color = group), size = 2.1, alpha = 0.85) +
  scale_color_manual(values = group_palette, labels = group_labels, name = "Group") +
  scale_fill_manual(values = group_palette, labels = group_labels, name = "Group") +
  labs(
    title = "PCA of size-corrected lateral shape (PC2 vs PC3): CT time series + 1950 habitats",
    subtitle = "Residuals from procD.lm(shape ~ logCsize) computed on the combined subset",
    x = pc2_lab,
    y = pc3_lab
  ) +
  coord_equal() +
  theme_classic(base_family = "Helvetica", base_size = 9) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text  = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8)
  )

print(p_fig7_pc23)

# ============================================================
# Save outputs
# ============================================================

pdf_path <- file.path(fig_dir, "Fig7B_PCA_CTtimeseries_plus_1950habitats_groupHulls_RESID_logCsize_PC2_PC3.pdf")
png_path <- file.path(fig_dir, "Fig7B_PCA_CTtimeseries_plus_1950habitats_groupHulls_RESID_logCsize_PC2_PC3.png")

ggsave(pdf_path, p_fig7_pc23, width = 7, height = 5, units = "in")
ggsave(png_path, p_fig7_pc23, width = 7, height = 5, units = "in", dpi = 600)

message("Saved PC2 vs PC3 figure to: ", normalizePath(fig_dir, winslash = "/"))