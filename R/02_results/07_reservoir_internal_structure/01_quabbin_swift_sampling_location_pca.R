# ============================================================
# Scripts/1950_Quabbin_vs_Swift_PC12_hullsHabitat_pointsByCollection.R
#
# PCA plots for Quabbin vs Swift River, plus one additional figure
# including Connecticut River time-series individuals (1950, 1956, 1970)
#
# Saves:
#   1) PC1 vs PC2 (Quabbin + Swift only)
#   2) PC2 vs PC3 (Quabbin + Swift only)
#   3) PC1 vs PC2 with enlarged outlier
#   4) PC1 vs PC2 including CT time-series individuals
#      saved as quabbin_swift_individuals_plustime.pdf/.png
# ============================================================

suppressPackageStartupMessages({
  library(geomorph)
  library(ggplot2)
  library(dplyr)
})

source("R/00_setup_morpho.R")
source("R/01_build_metadata.R")
source("R/02_subset_1950.R")
source("R/03_subset_CT_timeseries.R")   # << included as requested

# ---- 1) PCA on ALL 1950 residual data ----
pca_1950 <- gm.prcomp(coords_resid_1950)
pct_1950 <- 100 * (pca_1950$sdev^2 / sum(pca_1950$sdev^2))

pc_df <- data.frame(
  specimen      = rownames(pca_1950$x),
  PC1           = pca_1950$x[, 1],
  PC2           = pca_1950$x[, 2],
  PC3           = pca_1950$x[, 3],
  habitat       = gdf_1950$habitat,
  collection_id = gdf_1950$collection_id,
  year          = gdf_1950$year,
  stringsAsFactors = FALSE
)

stopifnot(identical(pc_df$specimen, gdf_1950$specimen))

# ---- 2) subset to Quabbin + Swift ----
plot_df <- pc_df %>%
  filter(habitat %in% c("Quabbin", "Swift River")) %>%
  mutate(
    point_group = case_when(
      habitat == "Swift River" ~ "Swift River",
      habitat == "Quabbin" & collection_id == "F2379" ~ "bottom",
      habitat == "Quabbin" & collection_id == "F2374" ~ "middle",
      habitat == "Quabbin" & collection_id == "F2385" ~ "top",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(point_group)) %>%
  mutate(
    is_outlier = specimen == "F2384_2_side"
  )

plot_df$point_group <- factor(
  plot_df$point_group,
  levels = c("bottom", "middle", "top", "Swift River")
)

# ---- 3) palettes ----
point_palette <- c(
  "bottom"      = "#E69F00",
  "middle"      = "#F0E442",
  "top"         = "#A6761D",
  "Swift River" = "tomato"
)

fill_palette <- c(
  "Quabbin"     = "darkgoldenrod2",
  "Swift River" = "tomato"
)

line_palette <- c(
  "Quabbin"     = "darkgoldenrod4",
  "Swift River" = "tomato3"
)

ct_palette <- c(
  "CT_1950" = "#1B9E77",
  "CT_1956" = "#7570B3",
  "CT_1970" = "#66A61E"
)

# ---- 4) helper to make closed hulls ----
make_closed_hulls <- function(df, xvar, yvar) {
  hull_df <- df %>%
    group_by(habitat) %>%
    slice(chull(.data[[xvar]], .data[[yvar]])) %>%
    ungroup()
  
  hull_df <- hull_df %>%
    group_by(habitat) %>%
    slice(c(1:n(), 1)) %>%
    ungroup()
  
  hull_df
}

# ---- 5) helper to make standard plot ----
make_pca_plot <- function(plot_df, hull_df, xvar, yvar, xlab, ylab, title, subtitle) {
  ggplot() +
    geom_polygon(
      data = hull_df,
      aes(x = .data[[xvar]], y = .data[[yvar]], group = habitat, fill = habitat),
      alpha = 0.12,
      color = NA
    ) +
    geom_path(
      data = hull_df,
      aes(x = .data[[xvar]], y = .data[[yvar]], group = habitat, color = habitat),
      linewidth = 0.9,
      show.legend = FALSE
    ) +
    geom_point(
      data = plot_df,
      aes(x = .data[[xvar]], y = .data[[yvar]], color = point_group),
      size = 2.6,
      alpha = 0.95
    ) +
    scale_fill_manual(values = fill_palette, name = "Habitat hull") +
    scale_color_manual(
      values = c(line_palette, point_palette),
      breaks = c("bottom", "middle", "top", "Swift River"),
      labels = c(
        "bottom" = "Quabbin - bottom",
        "middle" = "Quabbin - middle",
        "top" = "Quabbin - top",
        "Swift River" = "Swift River"
      ),
      name = "Points"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab
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
}

# ---- 6) helper to make PC1 vs PC2 plot with enlarged outlier ----
make_pca_plot_outlier_big <- function(plot_df, hull_df, xvar, yvar, xlab, ylab, title, subtitle) {
  ggplot() +
    geom_polygon(
      data = hull_df,
      aes(x = .data[[xvar]], y = .data[[yvar]], group = habitat, fill = habitat),
      alpha = 0.12,
      color = NA
    ) +
    geom_path(
      data = hull_df,
      aes(x = .data[[xvar]], y = .data[[yvar]], group = habitat, color = habitat),
      linewidth = 0.9,
      show.legend = FALSE
    ) +
    geom_point(
      data = plot_df,
      aes(
        x = .data[[xvar]],
        y = .data[[yvar]],
        color = point_group,
        size = is_outlier
      ),
      alpha = 0.95
    ) +
    scale_size_manual(
      values = c(`FALSE` = 2.6, `TRUE` = 5.0),
      guide = "none"
    ) +
    scale_fill_manual(values = fill_palette, name = "Habitat hull") +
    scale_color_manual(
      values = c(line_palette, point_palette),
      breaks = c("bottom", "middle", "top", "Swift River"),
      labels = c(
        "bottom" = "Quabbin - bottom",
        "middle" = "Quabbin - middle",
        "top" = "Quabbin - top",
        "Swift River" = "Swift River"
      ),
      name = "Points"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = xlab,
      y = ylab
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
}

# ---- 7) PC1 vs PC2 ----
hull_df_12 <- make_closed_hulls(plot_df, "PC1", "PC2")

p12 <- make_pca_plot(
  plot_df = plot_df,
  hull_df = hull_df_12,
  xvar = "PC1",
  yvar = "PC2",
  xlab = sprintf("PC1 (%.1f%%)", pct_1950[1]),
  ylab = sprintf("PC2 (%.1f%%)", pct_1950[2]),
  title = "PCA of size-corrected lateral shape: Quabbin vs Swift River (1950)",
  subtitle = "PC1 vs PC2; PCA computed on all 1950 residual data; hulls by habitat"
)

print(p12)

# ---- 8) PC2 vs PC3 ----
hull_df_23 <- make_closed_hulls(plot_df, "PC2", "PC3")

p23 <- make_pca_plot(
  plot_df = plot_df,
  hull_df = hull_df_23,
  xvar = "PC2",
  yvar = "PC3",
  xlab = sprintf("PC2 (%.1f%%)", pct_1950[2]),
  ylab = sprintf("PC3 (%.1f%%)", pct_1950[3]),
  title = "PCA of size-corrected lateral shape: Quabbin vs Swift River (1950)",
  subtitle = "PC2 vs PC3; PCA computed on all 1950 residual data; hulls by habitat"
)

print(p23)

# ---- 9) separate PC1 vs PC2 with enlarged outlier ----
p12_outlier_big <- make_pca_plot_outlier_big(
  plot_df = plot_df,
  hull_df = hull_df_12,
  xvar = "PC1",
  yvar = "PC2",
  xlab = sprintf("PC1 (%.1f%%)", pct_1950[1]),
  ylab = sprintf("PC2 (%.1f%%)", pct_1950[2]),
  title = "PCA of size-corrected lateral shape: Quabbin vs Swift River (1950)",
  subtitle = "PC1 vs PC2; outlier F2384_2_side shown with larger point"
)

print(p12_outlier_big)

# ============================================================
# ---- 10) NEW FIGURE: Quabbin + Swift individuals PLUS CT time series ----
# Rebuild a combined subset and run a PCA on that combined residual dataset
# so the time points are in the same ordination as Quabbin/Swift.
# ============================================================

gdf_qs_ct <- gdf %>%
  filter(
    (!is.na(habitat) & habitat %in% c("Quabbin", "Swift River")) |
      (!is.na(habitat) & habitat == "Connecticut River" & year %in% c(1950, 1956, 1970))
  ) %>%
  mutate(
    point_group = case_when(
      habitat == "Swift River" ~ "Swift River",
      habitat == "Quabbin" & collection_id == "F2379" ~ "bottom",
      habitat == "Quabbin" & collection_id == "F2374" ~ "middle",
      habitat == "Quabbin" & collection_id == "F2385" ~ "top",
      habitat == "Connecticut River" & year == 1950 ~ "CT_1950",
      habitat == "Connecticut River" & year == 1956 ~ "CT_1956",
      habitat == "Connecticut River" & year == 1970 ~ "CT_1970",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(point_group)) %>%
  droplevels()

coords_qs_ct <- subset_coords_to_gdf(coords_gpa, gdf_qs_ct)
stopifnot(identical(dimnames(coords_qs_ct)[[3]], gdf_qs_ct$specimen))

size_vec_qs_ct <- setNames(gdf_qs_ct$size_for_allometry, gdf_qs_ct$specimen)
allo_qs_ct <- allometry_residuals(coords_qs_ct, size_vec_qs_ct)
coords_resid_qs_ct <- allo_qs_ct$residuals

pca_qs_ct <- gm.prcomp(coords_resid_qs_ct)
pct_qs_ct <- 100 * (pca_qs_ct$sdev^2 / sum(pca_qs_ct$sdev^2))

pc_df_qs_ct <- data.frame(
  specimen      = rownames(pca_qs_ct$x),
  PC1           = pca_qs_ct$x[, 1],
  PC2           = pca_qs_ct$x[, 2],
  PC3           = pca_qs_ct$x[, 3],
  habitat       = gdf_qs_ct$habitat,
  collection_id = gdf_qs_ct$collection_id,
  year          = gdf_qs_ct$year,
  point_group   = gdf_qs_ct$point_group,
  stringsAsFactors = FALSE
)

stopifnot(identical(pc_df_qs_ct$specimen, gdf_qs_ct$specimen))

# hulls still only for Quabbin + Swift
hull_df_qs_ct <- pc_df_qs_ct %>%
  filter(habitat %in% c("Quabbin", "Swift River")) %>%
  group_by(habitat) %>%
  slice(chull(PC1, PC2)) %>%
  ungroup() %>%
  group_by(habitat) %>%
  slice(c(1:n(), 1)) %>%
  ungroup()

# separate point data for plotting
qs_points <- pc_df_qs_ct %>%
  filter(point_group %in% c("bottom", "middle", "top", "Swift River")) %>%
  mutate(
    point_group = factor(point_group, levels = c("bottom", "middle", "top", "Swift River"))
  )

ct_points <- pc_df_qs_ct %>%
  filter(point_group %in% c("CT_1950", "CT_1956", "CT_1970")) %>%
  mutate(
    point_group = factor(point_group, levels = c("CT_1950", "CT_1956", "CT_1970"))
  )

p12_plustime <- ggplot() +
  geom_polygon(
    data = hull_df_qs_ct,
    aes(x = PC1, y = PC2, group = habitat, fill = habitat),
    alpha = 0.12,
    color = NA
  ) +
  geom_path(
    data = hull_df_qs_ct,
    aes(x = PC1, y = PC2, group = habitat, color = habitat),
    linewidth = 0.9,
    show.legend = FALSE
  ) +
  geom_point(
    data = qs_points,
    aes(x = PC1, y = PC2, color = point_group),
    size = 2.6,
    alpha = 0.95
  ) +
  geom_point(
    data = ct_points,
    aes(x = PC1, y = PC2, color = point_group, shape = point_group),
    size = 2.8,
    alpha = 0.95
  ) +
  scale_fill_manual(values = fill_palette, name = "Habitat hull") +
  scale_color_manual(
    values = c(point_palette, ct_palette),
    breaks = c("bottom", "middle", "top", "Swift River", "CT_1950", "CT_1956", "CT_1970"),
    labels = c(
      "bottom" = "Quabbin - bottom",
      "middle" = "Quabbin - middle",
      "top" = "Quabbin - top",
      "Swift River" = "Swift River",
      "CT_1950" = "Connecticut River 1950",
      "CT_1956" = "Connecticut River 1956",
      "CT_1970" = "Connecticut River 1970"
    ),
    name = "Points"
  ) +
  scale_shape_manual(
    values = c(
      "CT_1950" = 15,
      "CT_1956" = 17,
      "CT_1970" = 18
    ),
    guide = "none"
  ) +
  labs(
    title = "PCA of size-corrected lateral shape: Quabbin, Swift River, and CT time series",
    subtitle = "PC1 vs PC2; PCA computed on combined residual data for Quabbin + Swift + Connecticut River (1950, 1956, 1970)",
    x = sprintf("PC1 (%.1f%%)", pct_qs_ct[1]),
    y = sprintf("PC2 (%.1f%%)", pct_qs_ct[2])
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

print(p12_plustime)



# ============================================================
# PCA: Quabbin + Swift ONLY
# Goal: visualize internal structure within Quabbin
#       with hulls by habitat
# ============================================================


# ---- 1) subset to Quabbin + Swift ----
gdf_qs <- gdf %>%
  filter(habitat %in% c("Quabbin", "Swift River")) %>%
  mutate(
    point_group = case_when(
      habitat == "Swift River" ~ "Swift River",
      habitat == "Quabbin" & collection_id == "F2379" ~ "bottom",
      habitat == "Quabbin" & collection_id == "F2374" ~ "middle",
      habitat == "Quabbin" & collection_id == "F2385" ~ "top",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(point_group)) %>%
  droplevels()

# ---- 2) subset coords ----
coords_qs <- subset_coords_to_gdf(coords_gpa, gdf_qs)
stopifnot(identical(dimnames(coords_qs)[[3]], gdf_qs$specimen))

# ---- 3) allometry correction ----
size_vec <- setNames(gdf_qs$size_for_allometry, gdf_qs$specimen)
allo_qs <- allometry_residuals(coords_qs, size_vec)
coords_resid_qs <- allo_qs$residuals

# ---- 4) PCA ----
pca_qs <- gm.prcomp(coords_resid_qs)
pct_qs <- 100 * (pca_qs$sdev^2 / sum(pca_qs$sdev^2))

pc_df <- data.frame(
  specimen      = rownames(pca_qs$x),
  PC1           = pca_qs$x[, 1],
  PC2           = pca_qs$x[, 2],
  PC3           = pca_qs$x[, 3],
  habitat       = gdf_qs$habitat,
  point_group   = gdf_qs$point_group,
  collection_id = gdf_qs$collection_id,
  stringsAsFactors = FALSE
)

stopifnot(identical(pc_df$specimen, gdf_qs$specimen))

# ---- 5) palettes ----
point_palette <- c(
  "bottom"      = "#E69F00",
  "middle"      = "#F0E442",
  "top"         = "#A6761D",
  "Swift River" = "tomato"
)

fill_palette <- c(
  "Quabbin"     = "darkgoldenrod2",
  "Swift River" = "tomato"
)

line_palette <- c(
  "Quabbin"     = "darkgoldenrod4",
  "Swift River" = "tomato3"
)

# ---- 6) helper to make closed hulls ----
make_closed_hulls <- function(df, xvar, yvar) {
  hull_df <- df %>%
    group_by(habitat) %>%
    slice(chull(.data[[xvar]], .data[[yvar]])) %>%
    ungroup()
  
  hull_df <- hull_df %>%
    group_by(habitat) %>%
    slice(c(1:n(), 1)) %>%
    ungroup()
  
  hull_df
}

# ---- 7) build hulls for PC1 vs PC2 ----
hull_df_12 <- make_closed_hulls(pc_df, "PC1", "PC2")

# ---- 8) plot PC1 vs PC2 ----
p12 <- ggplot() +
  geom_polygon(
    data = hull_df_12,
    aes(x = PC1, y = PC2, group = habitat, fill = habitat),
    alpha = 0.12,
    color = NA
  ) +
  geom_path(
    data = hull_df_12,
    aes(x = PC1, y = PC2, group = habitat, color = habitat),
    linewidth = 0.9,
    show.legend = FALSE
  ) +
  geom_point(
    data = pc_df,
    aes(x = PC1, y = PC2, color = point_group),
    size = 2.8,
    alpha = 0.95
  ) +
  scale_fill_manual(values = fill_palette, name = "Habitat hull") +
  scale_color_manual(
    values = c(line_palette, point_palette),
    breaks = c("bottom", "middle", "top", "Swift River"),
    labels = c(
      "bottom" = "Quabbin - bottom",
      "middle" = "Quabbin - middle",
      "top" = "Quabbin - top",
      "Swift River" = "Swift River"
    ),
    name = "Points"
  ) +
  labs(
    title = "PCA of size-corrected lateral shape: Quabbin + Swift River",
    subtitle = "PC1 vs PC2; PCA computed only on Quabbin and Swift individuals; hulls by habitat",
    x = sprintf("PC1 (%.1f%%)", pct_qs[1]),
    y = sprintf("PC2 (%.1f%%)", pct_qs[2])
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

print(p12)

# ---- 9) save ----
fig_dir <- "Figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pdf_path_12 <- file.path(fig_dir, "quabbin_swift_local_PCA_PC12_hulls.pdf")
png_path_12 <- file.path(fig_dir, "quabbin_swift_local_PCA_PC12_hulls.png")

ggsave(pdf_path_12, p12, width = 7, height = 5, units = "in")
ggsave(png_path_12, p12, width = 7, height = 5, units = "in", dpi = 600)

message(
  "Saved:\n  ",
  normalizePath(pdf_path_12, winslash = "/"),
  "\n  ",
  normalizePath(png_path_12, winslash = "/")
)


# ---- 11) save ----
fig_dir <- "Figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pdf_path_12 <- file.path(fig_dir, "Fig_1950_Quabbin_vs_Swift_onFull1950PC12_simple.pdf")
png_path_12 <- file.path(fig_dir, "Fig_1950_Quabbin_vs_Swift_onFull1950PC12_simple.png")

pdf_path_23 <- file.path(fig_dir, "Fig_1950_Quabbin_vs_Swift_onFull1950PC23_simple.pdf")
png_path_23 <- file.path(fig_dir, "Fig_1950_Quabbin_vs_Swift_onFull1950PC23_simple.png")

pdf_path_12_outlier <- file.path(fig_dir, "Fig_1950_Quabbin_vs_Swift_onFull1950PC12_outlierBig.pdf")
png_path_12_outlier <- file.path(fig_dir, "Fig_1950_Quabbin_vs_Swift_onFull1950PC12_outlierBig.png")

# requested new filename
pdf_path_plustime <- file.path(fig_dir, "quabbin_swift_individuals_plustime.pdf")
png_path_plustime <- file.path(fig_dir, "quabbin_swift_individuals_plustime.png")

ggsave(pdf_path_12, p12, width = 7, height = 5, units = "in")
ggsave(png_path_12, p12, width = 7, height = 5, units = "in", dpi = 600)

ggsave(pdf_path_23, p23, width = 7, height = 5, units = "in")
ggsave(png_path_23, p23, width = 7, height = 5, units = "in", dpi = 600)

ggsave(pdf_path_12_outlier, p12_outlier_big, width = 7, height = 5, units = "in")
ggsave(png_path_12_outlier, p12_outlier_big, width = 7, height = 5, units = "in", dpi = 600)

ggsave(pdf_path_plustime, p12_plustime, width = 7, height = 5, units = "in")
ggsave(png_path_plustime, p12_plustime, width = 7, height = 5, units = "in", dpi = 600)

message(
  "Saved:\n  ",
  normalizePath(pdf_path_12, winslash = "/"),
  "\n  ",
  normalizePath(png_path_12, winslash = "/"),
  "\n  ",
  normalizePath(pdf_path_23, winslash = "/"),
  "\n  ",
  normalizePath(png_path_23, winslash = "/"),
  "\n  ",
  normalizePath(pdf_path_12_outlier, winslash = "/"),
  "\n  ",
  normalizePath(png_path_12_outlier, winslash = "/"),
  "\n  ",
  normalizePath(pdf_path_plustime, winslash = "/"),
  "\n  ",
  normalizePath(png_path_plustime, winslash = "/")
)