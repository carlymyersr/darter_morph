# ============================================================
# Scripts/Figure6_PCA_CTtimeseries_plus_1950habitats_RAW.R
# FIGURE 6 (RAW): PCA combining:
#   - CT time series (1950/1956/1970)
#   - All 1950 habitats
#   - Points colored by GROUP, hulls by GROUP
# Saves to: Figures/
# ============================================================

suppressPackageStartupMessages({
  library(geomorph)
  library(ggplot2)
  library(dplyr)
})

# Canonical upstream
source("R/00_setup_morpho.R")
source("R/01_build_metadata.R")
source("R/04_subset_CT_timeseries_plus_1950habitats.R")

# ---- 0) PCA on combined subset (RAW) ----
pca_Fig6 <- geomorph::gm.prcomp(coords_Fig6)
pct_Fig6 <- 100 * (pca_Fig6$sdev^2 / sum(pca_Fig6$sdev^2))

pc1_lab <- sprintf("PC1 (%.1f%%)", pct_Fig6[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", pct_Fig6[2])

# ---- 1) Plotting df ----
pca_Fig6_df <- data.frame(
  specimen = rownames(pca_Fig6$x),
  PC1      = pca_Fig6$x[, 1],
  PC2      = pca_Fig6$x[, 2],
  group    = gdf_Fig6$group,
  habitat  = gdf_Fig6$habitat,
  year     = gdf_Fig6$year,
  SL_mm    = gdf_Fig6$SL_mm,
  logSL    = gdf_Fig6$logSL,
  stringsAsFactors = FALSE
)
stopifnot(identical(pca_Fig6_df$specimen, gdf_Fig6$specimen))

# ---- 2) Hulls by GROUP ----
if (!exists("make_hulls", inherits = TRUE)) {
  stop("make_hulls() not found. It should be defined in R/00_setup_morpho.R")
}

hull_Fig6_df <- make_hulls(
  df = pca_Fig6_df,
  group_col = "group",
  x = "PC1",
  y = "PC2",
  min_n = 3
)

# ---- 3) Group palette + labels ----
# CT years: blues; other 1950 habitats: pull from hab_palette if available, else default colors
group_levels <- levels(pca_Fig6_df$group)

group_palette <- setNames(rep(NA_character_, length(group_levels)), group_levels)

# CT time-series colors (fixed)
group_palette["CT_1950"] <- "steelblue"
group_palette["CT_1956"] <- "dodgerblue4"
group_palette["CT_1970"] <- "navy"

# Other habitats 1950
other_1950 <- setdiff(group_levels, c("CT_1950","CT_1956","CT_1970"))
if (length(other_1950) > 0) {
  if (exists("hab_palette", inherits = TRUE) && is.vector(hab_palette) && !is.null(names(hab_palette))) {
    # map "Quabbin_1950" -> "Quabbin"
    hab_names <- sub("_1950$", "", other_1950)
    cols <- hab_palette[hab_names]
    # If any are missing in hab_palette, fill them deterministically
    if (any(is.na(cols))) {
      missing <- hab_names[is.na(cols)]
      fallback <- scales::hue_pal()(length(missing))
      names(fallback) <- missing
      cols[is.na(cols)] <- fallback[missing]
    }
    names(cols) <- other_1950
    group_palette[other_1950] <- cols
  } else {
    # fallback if hab_palette doesn't exist
    group_palette[other_1950] <- scales::hue_pal()(length(other_1950))
  }
}

# labels
group_labels <- setNames(group_levels, group_levels)
group_labels["CT_1950"] <- "Connecticut River (1950)"
group_labels["CT_1956"] <- "Connecticut River (1956)"
group_labels["CT_1970"] <- "Connecticut River (1970)"
for (g in other_1950) {
  group_labels[g] <- gsub("_1950$", " (1950)", g)
}

missing_cols <- names(group_palette)[is.na(group_palette)]
if (length(missing_cols) > 0) stop("group_palette still missing colors for: ", paste(missing_cols, collapse = ", "))

# ---- 4) Plot ----
p_fig6 <- ggplot(pca_Fig6_df, aes(PC1, PC2)) +
  geom_polygon(
    data = hull_Fig6_df,
    aes(group = group, fill = group),
    alpha = 0.15,
    color = NA
  ) +
  geom_polygon(
    data = hull_Fig6_df,
    aes(group = group, color = group),
    fill = NA,
    linewidth = 0.7
  ) +
  geom_point(aes(color = group), size = 2.1, alpha = 0.85) +
  scale_color_manual(values = group_palette, labels = group_labels, name = "Group") +
  scale_fill_manual(values = group_palette,  labels = group_labels, name = "Group") +
  labs(
    title = "PCA of raw GPA lateral shape: CT time series + 1950 habitats",
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

print(p_fig6)

# ---- 5) Save to Figures/ ----
fig_dir <- "Figures"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pdf_path <- file.path(fig_dir, "Fig6_PCA_CTtimeseries_plus_1950habitats_groupHulls_RAW.pdf")
png_path <- file.path(fig_dir, "Fig6_PCA_CTtimeseries_plus_1950habitats_groupHulls_RAW.png")

ggsave(pdf_path, p_fig6, width = 7, height = 5, units = "in")
ggsave(png_path, p_fig6, width = 7, height = 5, units = "in", dpi = 600)

message("Saved:\n  ", normalizePath(pdf_path, winslash = "/"),
        "\n  ", normalizePath(png_path, winslash = "/"))