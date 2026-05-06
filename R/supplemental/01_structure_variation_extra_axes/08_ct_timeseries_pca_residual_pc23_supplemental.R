# ============================================================
# Scripts/Figure5_PCA_CT_timeseries_RESID_logCsize_PC2_PC3.R
# FIGURE 5 variant (RESID): PCA CT time-series (1950, 1956, 1970)
#   - Residuals from procD.lm(shape ~ logCsize) computed on CT subset
#   - PCA on residual coords
#   - Plot PC2 vs PC3
#   - Points colored by GROUP + convex hulls by GROUP
# Saves to: Figures/
# ============================================================

suppressPackageStartupMessages({
  library(geomorph)
  library(ggplot2)
  library(dplyr)
})

# Canonical upstream
source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")
source("R/methods/01_specimen_sampling_study_design/03_subset_CT_timeseries.R")   # creates coords_resid_CT3 + gdf_CT3

# ---- 0) PCA on RESID subset ----
pca_CT3_resid <- geomorph::gm.prcomp(coords_resid_CT3)
pct_CT3_resid <- 100 * (pca_CT3_resid$sdev^2 / sum(pca_CT3_resid$sdev^2))

pc2_lab <- sprintf("PC2 (%.1f%%)", pct_CT3_resid[2])
pc3_lab <- sprintf("PC3 (%.1f%%)", pct_CT3_resid[3])

# ---- 1) Plotting df ----
pca_CT3_resid_df <- data.frame(
  specimen = rownames(pca_CT3_resid$x),
  PC2      = pca_CT3_resid$x[, 2],
  PC3      = pca_CT3_resid$x[, 3],
  group    = gdf_CT3$group,
  habitat  = gdf_CT3$habitat,
  year     = gdf_CT3$year,
  SL_mm    = gdf_CT3$SL_mm,
  logSL    = gdf_CT3$logSL,
  stringsAsFactors = FALSE
)
stopifnot(identical(pca_CT3_resid_df$specimen, gdf_CT3$specimen))

# ---- 2) Hulls by GROUP ----
if (!exists("make_hulls", inherits = TRUE)) {
  stop("make_hulls() not found. It should be defined in R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
}

hull_CT3_resid_df <- make_hulls(
  df = pca_CT3_resid_df,
  group_col = "group",
  x = "PC2",
  y = "PC3",
  min_n = 3
)

# ---- 3) Palette + labels ----
group_palette <- c(
  "CT_1950" = "steelblue",
  "CT_1956" = "dodgerblue4",
  "CT_1970" = "navy",
  "CT_1979" = "orange"
)

group_labels <- c(
  "CT_1950" = "Connecticut River (1950)",
  "CT_1956" = "Connecticut River (1956)",
  "CT_1970" = "Connecticut River (1970)",
  "CT_1979" = "Connecticut River (1979)"
)

missing_cols <- setdiff(levels(pca_CT3_resid_df$group), names(group_palette))
if (length(missing_cols) > 0) {
  stop("group_palette missing colors for: ", paste(missing_cols, collapse = ", "))
}

# ---- 4) Plot ----
p_fig5_pc23 <- ggplot(pca_CT3_resid_df, aes(PC2, PC3)) +
  geom_polygon(
    data = hull_CT3_resid_df,
    aes(group = group, fill = group),
    alpha = 0.15,
    color = NA
  ) +
  geom_polygon(
    data = hull_CT3_resid_df,
    aes(group = group, color = group),
    fill = NA,
    linewidth = 0.7
  ) +
  geom_point(aes(color = group), size = 2.1, alpha = 0.85) +
  scale_color_manual(values = group_palette, labels = group_labels, name = "Group") +
  scale_fill_manual(values = group_palette, labels = group_labels, name = "Group") +
  labs(
    title = "PCA of size-corrected lateral shape (residuals ~ logCsize): Connecticut River (1950 / 1956 / 1970)",
    subtitle = "PC2 vs. PC3; residuals from procD.lm(shape ~ logCsize) computed within CT subset",
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

print(p_fig5_pc23)

# ---- 5) Save to Figures/ ----
fig_dir <- "Figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pdf_path <- file.path(fig_dir, "Fig5_PCA_CT1950_CT1956_CT1970_groupHulls_RESID_logCsize_PC2_PC3.pdf")
png_path <- file.path(fig_dir, "Fig5_PCA_CT1950_CT1956_CT1970_groupHulls_RESID_logCsize_PC2_PC3.png")

ggsave(pdf_path, p_fig5_pc23, width = 7, height = 5, units = "in")
ggsave(png_path, p_fig5_pc23, width = 7, height = 5, units = "in", dpi = 600)

message("Saved:\n  ", normalizePath(pdf_path, winslash = "/"),
        "\n  ", normalizePath(png_path, winslash = "/"))