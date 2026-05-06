# ============================================================
# Scripts/Figure5_PCA_CT_timeseries_RESID_logCsize.R
# FIGURE 5 (RESID): PCA CT time-series (1950, 1956, 1970)
#   - Residuals from procD.lm(shape ~ logCsize) computed on CT subset
#   - PCA on residual coords
#   - Points colored by GROUP + convex hulls by GROUP
#   - PC1 inverted for plotting
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

pc1_lab <- sprintf("-PC1 (%.1f%%)", pct_CT3_resid[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", pct_CT3_resid[2])

# ---- 1) Plotting df ----
pca_CT3_resid_df <- data.frame(
  specimen = rownames(pca_CT3_resid$x),
  PC1      = pca_CT3_resid$x[, 1],
  PC2      = pca_CT3_resid$x[, 2],
  group    = gdf_CT3$group,
  habitat  = gdf_CT3$habitat,
  year     = gdf_CT3$year,
  SL_mm    = gdf_CT3$SL_mm,
  logSL    = gdf_CT3$logSL,
  stringsAsFactors = FALSE
)
stopifnot(identical(pca_CT3_resid_df$specimen, gdf_CT3$specimen))

# ---- 1b) Invert PC1 BEFORE hull computation ----
pca_CT3_resid_df$PC1 <- -pca_CT3_resid_df$PC1

# ---- 2) Hulls by GROUP ----
if (!exists("make_hulls", inherits = TRUE)) {
  stop("make_hulls() not found. It should be defined in R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
}

hull_CT3_resid_df <- make_hulls(
  df = pca_CT3_resid_df,
  group_col = "group",
  x = "PC1",
  y = "PC2",
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
p_fig5 <- ggplot(pca_CT3_resid_df, aes(PC1, PC2)) +
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
    subtitle = "Residuals from procD.lm(shape ~ logCsize) computed within CT subset",
    x = pc1_lab,
    y = pc2_lab
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

print(p_fig5)

# ---- 5) Save to Figures/ ----
fig_dir <- "Figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pdf_path <- file.path(fig_dir, "Fig5_PCA_CT1950_CT1956_CT1970_groupHulls_RESID_logCsize.pdf")
png_path <- file.path(fig_dir, "Fig5_PCA_CT1950_CT1956_CT1970_groupHulls_RESID_logCsize.png")

ggsave(pdf_path, p_fig5, width = 7, height = 5, units = "in")
ggsave(png_path, p_fig5, width = 7, height = 5, units = "in", dpi = 600)

message("Saved:\n  ", normalizePath(pdf_path, winslash = "/"),
        "\n  ", normalizePath(png_path, winslash = "/"))