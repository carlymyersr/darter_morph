# ============================================================
# Scripts/Figure4_PCA_CT_timeseries_RAW.R
# FIGURE 4 (RAW): PCA CT time-series (1950, 1956, 1970)
#   - Points colored by GROUP
#   - Convex hulls by GROUP
# Saves to: Figures/
# ============================================================

suppressPackageStartupMessages({
  library(geomorph)
  library(ggplot2)
  library(dplyr)
  library(viridis)
})

# Canonical upstream
source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")
source("R/methods/01_specimen_sampling_study_design/03_subset_CT_timeseries.R")

# ---- 0) PCA on combined CT subset (RAW) ----
pca_CT3 <- geomorph::gm.prcomp(coords_CT3)
pct_CT3 <- 100 * (pca_CT3$sdev^2 / sum(pca_CT3$sdev^2))

pc1_lab <- sprintf("PC1 (%.1f%%)", pct_CT3[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", pct_CT3[2])

# ---- 1) Plotting df ----
pca_CT3_df <- data.frame(
  specimen = rownames(pca_CT3$x),
  PC1      = pca_CT3$x[, 1],
  PC2      = pca_CT3$x[, 2],
  group    = gdf_CT3$group,
  habitat  = gdf_CT3$habitat,
  year     = gdf_CT3$year,
  Csize    = gdf_CT3$Csize,
  logCsize = gdf_CT3$logCsize,
  SL_mm    = gdf_CT3$SL_mm,
  logSL    = gdf_CT3$logSL,
  stringsAsFactors = FALSE
)
stopifnot(identical(pca_CT3_df$specimen, gdf_CT3$specimen))

# ---- 2) Hulls by GROUP ----
if (!exists("make_hulls", inherits = TRUE)) {
  stop("make_hulls() not found. It should be defined in R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
}

hull_CT3_df <- make_hulls(
  df = pca_CT3_df,
  group_col = "group",
  x = "PC1",
  y = "PC2"
)

# ---- 3) Palette + labels ----
group_palette <- c(
  "CT_1950" = "steelblue",
  "CT_1956" = "dodgerblue4",
  "CT_1970" = "navy"
)

group_labels <- c(
  "CT_1950" = "Connecticut River (1950)",
  "CT_1956" = "Connecticut River (1956)",
  "CT_1970" = "Connecticut River (1970)"
)

missing_cols <- setdiff(levels(pca_CT3_df$group), names(group_palette))
if (length(missing_cols) > 0) stop("group_palette missing colors for: ", paste(missing_cols, collapse = ", "))

# ---- 4) Plot ----
p_fig4 <- ggplot(pca_CT3_df, aes(PC1, PC2)) +
  geom_polygon(
    data = hull_CT3_df,
    aes(group = group, fill = group),
    alpha = 0.15,
    color = NA
  ) +
  geom_polygon(
    data = hull_CT3_df,
    aes(group = group, color = group),
    fill = NA,
    linewidth = 0.7
  ) +
  geom_point(aes(color = group), size = 2.1, alpha = 0.85) +
  scale_color_manual(values = group_palette, labels = group_labels, name = "Group") +
  scale_fill_manual(values = group_palette,  labels = group_labels, name = "Group") +
  labs(
    title = "PCA of raw GPA lateral shape: Connecticut River (1950 / 1956 / 1970)",
    x = pc1_lab,
    y = pc2_lab
  ) +
  coord_equal() +
  theme_classic(base_family = "Helvetica", base_size = 9) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 11, face = "bold"),
    axis.title = element_text(size = 9),
    axis.text  = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8)
  )

print(p_fig4)

# ---- 4b) Plot: raw PCA with points colored by centroid size ----
p_fig4_csize <- ggplot(pca_CT3_df, aes(PC1, PC2)) +
  geom_polygon(
    data = hull_CT3_df,
    aes(group = group, fill = group),
    alpha = 0.15,
    color = NA
  ) +
  geom_polygon(
    data = hull_CT3_df,
    aes(group = group),
    fill = NA,
    linewidth = 0.7,
    color = "black"
  ) +
  geom_point(
    aes(color = Csize),
    size = 2.1,
    alpha = 0.90
  ) +
  scale_color_viridis_c(
    name = "Centroid size",
    option = "C"
  ) +
  scale_fill_manual(values = group_palette, labels = group_labels, name = "Group") +
  labs(
    title = "PCA of raw GPA lateral shape: Connecticut River (1950 / 1956 / 1970)",
    subtitle = "Points colored by centroid size",
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

print(p_fig4_csize)

# ---- 5b) Save centroid-size-colored version ----
pdf_path_csize <- file.path(fig_dir, "Fig4_PCA_CT1950_CT1956_CT1970_groupHulls_RAW_pointsByCsize.pdf")
png_path_csize <- file.path(fig_dir, "Fig4_PCA_CT1950_CT1956_CT1970_groupHulls_RAW_pointsByCsize.png")

ggsave(pdf_path_csize, p_fig4_csize, width = 7, height = 5, units = "in")
ggsave(png_path_csize, p_fig4_csize, width = 7, height = 5, units = "in", dpi = 600)


# ---- 6) Centroid size summary + stats across CT years ----

ct_csize_df <- pca_CT3_df %>%
  mutate(
    year = factor(year, levels = c(1950, 1956, 1970))
  ) %>%
  filter(!is.na(year), is.finite(Csize))

# summary table
ct_csize_summary <- ct_csize_df %>%
  group_by(year) %>%
  summarise(
    n = n(),
    mean_Csize = mean(Csize, na.rm = TRUE),
    sd_Csize = sd(Csize, na.rm = TRUE),
    se_Csize = sd_Csize / sqrt(n),
    min_Csize = min(Csize, na.rm = TRUE),
    max_Csize = max(Csize, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  ct_csize_summary,
  file.path(fig_dir, "CT_centroid_size_summary.csv"),
  row.names = FALSE
)

# ANOVA
fit_csize_aov <- aov(Csize ~ year, data = ct_csize_df)
csize_aov_tab <- as.data.frame(anova(fit_csize_aov))
csize_aov_tab$term <- rownames(csize_aov_tab)
rownames(csize_aov_tab) <- NULL

write.csv(
  csize_aov_tab,
  file.path(fig_dir, "CT_centroid_size_ANOVA.csv"),
  row.names = FALSE
)

# Tukey pairwise
csize_tukey <- TukeyHSD(fit_csize_aov, "year")
csize_tukey_df <- as.data.frame(csize_tukey$year)
csize_tukey_df$comparison <- rownames(csize_tukey_df)
rownames(csize_tukey_df) <- NULL

write.csv(
  csize_tukey_df,
  file.path(fig_dir, "CT_centroid_size_Tukey.csv"),
  row.names = FALSE
)

# console output
cat("\nCentroid size ANOVA:\n")
print(csize_aov_tab)

cat("\nCentroid size Tukey pairwise comparisons:\n")
print(csize_tukey_df)



# ---- 7) Boxplot of centroid size by year ----
p_csize_box <- ggplot(ct_csize_df, aes(x = year, y = Csize)) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = "white",
    color = "black"
  ) +
  geom_jitter(
    width = 0.08,
    height = 0,
    alpha = 0.75,
    size = 2
  ) +
  labs(
    title = "Centroid size across Connecticut River time points",
    subtitle = paste0(
      "ANOVA p = ",
      formatC(csize_aov_tab$`Pr(>F)`[csize_aov_tab$term == "year"], format = "e", digits = 2)
    ),
    x = "Year",
    y = "Centroid size"
  ) +
  theme_classic(base_family = "Helvetica", base_size = 9) +
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 8),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8)
  )

print(p_csize_box)

pdf_path_csize_box <- file.path(fig_dir, "Fig4_CT_centroid_size_boxplot.pdf")
png_path_csize_box <- file.path(fig_dir, "Fig4_CT_centroid_size_boxplot.png")



# ---- 5) Save to Figures/ ----
fig_dir <- "Figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(pdf_path_csize_box, p_csize_box, width = 5.5, height = 4.5, units = "in")
ggsave(png_path_csize_box, p_csize_box, width = 5.5, height = 4.5, units = "in", dpi = 600)

message(
  "Saved centroid-size boxplot + stats:\n  ",
  normalizePath(pdf_path_csize_box, winslash = "/"),
  "\n  ",
  normalizePath(png_path_csize_box, winslash = "/"),
  "\n  ",
  normalizePath(file.path(fig_dir, "CT_centroid_size_summary.csv"), winslash = "/"),
  "\n  ",
  normalizePath(file.path(fig_dir, "CT_centroid_size_ANOVA.csv"), winslash = "/"),
  "\n  ",
  normalizePath(file.path(fig_dir, "CT_centroid_size_Tukey.csv"), winslash = "/")
)

message(
  "Saved centroid-size-colored figure:\n  ",
  normalizePath(pdf_path_csize, winslash = "/"),
  "\n  ",
  normalizePath(png_path_csize, winslash = "/")
)

pdf_path <- file.path(fig_dir, "Fig4_PCA_CT1950_CT1956_CT1970_groupHulls_RAW.pdf")
png_path <- file.path(fig_dir, "Fig4_PCA_CT1950_CT1956_CT1970_groupHulls_RAW.png")

ggsave(pdf_path, p_fig4, width = 7, height = 5, units = "in")
ggsave(png_path, p_fig4, width = 7, height = 5, units = "in", dpi = 600)

message("Saved:\n  ", normalizePath(pdf_path, winslash = "/"),
        "\n  ", normalizePath(png_path, winslash = "/"))