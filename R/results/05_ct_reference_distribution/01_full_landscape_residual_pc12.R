# ============================================================
# Scripts/Figure7_RESID_CTtimeseries_plus_1950habitats.R
# FIGURE 7:
#   A) Residual diagnostic plot: regression score vs logCsize
#   B) PCA of residual coordinates (residuals ~ logCsize) with hulls by group
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
# A) Residual diagnostic plot (regression score vs logCsize)
#   - Uses geomorph::procD.lm regression scores from allo_Fig6$fit
# ============================================================

fit <- allo_Fig6$fit

# Regression scores (1D summary of shape variation along fitted axis)
# geomorph stores these in different places across versions; we’ll extract robustly.
get_reg_scores <- function(fit_obj) {
  # Try common slots:
  cand <- list(
    tryCatch(fit_obj$RS, error = function(e) NULL),
    tryCatch(fit_obj$RegScore, error = function(e) NULL),
    tryCatch(fit_obj$reg.score, error = function(e) NULL),
    tryCatch(fit_obj$regscore, error = function(e) NULL)
  )
  cand <- Filter(function(x) !is.null(x), cand)
  if (length(cand) == 0) return(NULL)
  as.numeric(cand[[1]])
}

rs <- get_reg_scores(fit)

# If regression scores unavailable, fall back to PC1 of RAW as a size check
if (is.null(rs) || length(rs) != nrow(gdf_Fig6)) {
  warning("Could not extract procD.lm regression scores; using PC1 of RAW coords as fallback for size diagnostic.")
  pca_tmp <- gm.prcomp(coords_Fig6)
  rs <- as.numeric(pca_tmp$x[, 1])
}

diag_df <- data.frame(
  specimen = gdf_Fig6$specimen,
  group    = gdf_Fig6$group,
  year     = gdf_Fig6$year,
  habitat  = gdf_Fig6$habitat,
  logCsize = gdf_Fig6$logCsize,
  score    = rs,
  stringsAsFactors = FALSE
)

# Palette/labels re-used
group_levels <- levels(diag_df$group)
group_palette <- setNames(rep(NA_character_, length(group_levels)), group_levels)
group_palette["CT_1950"] <- "steelblue"
group_palette["CT_1956"] <- "dodgerblue4"
group_palette["CT_1970"] <- "navy"

other <- setdiff(group_levels, c("CT_1950","CT_1956","CT_1970"))
if (length(other) > 0) {
  if (exists("hab_palette", inherits = TRUE) && is.vector(hab_palette) && !is.null(names(hab_palette))) {
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

p_resid_diag <- ggplot(diag_df, aes(x = logCsize, y = score)) +
  geom_point(aes(color = group), size = 2.1, alpha = 0.85) +
  scale_color_manual(values = group_palette, labels = group_labels, name = "Group") +
  labs(
    title = "Residual diagnostic: regression score vs logCsize",
    subtitle = "Combined dataset (CT time series + 1950 habitats); score from procD.lm(shape ~ logCsize)",
    x = "logCsize (centroid size from GPA)",
    y = "Regression score (shape ~ size)"
  ) +
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

print(p_resid_diag)

ggsave(file.path(fig_dir, "Fig7A_residual_diagnostic_score_vs_logCsize.pdf"),
       p_resid_diag, width = 7, height = 5, units = "in")
ggsave(file.path(fig_dir, "Fig7A_residual_diagnostic_score_vs_logCsize.png"),
       p_resid_diag, width = 7, height = 5, units = "in", dpi = 600)

# Also save the allometry ANOVA text (handy for writeup)
anova_txt <- file.path(fig_dir, "Fig7A_allometry_ANOVA_shape_vs_logCsize.txt")
writeLines(capture.output(anova(fit)), con = anova_txt)
message("Wrote: ", normalizePath(anova_txt, winslash = "/"))

# ============================================================
# B) PCA on residual coordinates (Fig 7)
# ============================================================

pca_Fig7 <- geomorph::gm.prcomp(coords_resid_Fig6)
pct_Fig7 <- 100 * (pca_Fig7$sdev^2 / sum(pca_Fig7$sdev^2))

pc1_lab <- sprintf("PC1 (%.1f%%)", pct_Fig7[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", pct_Fig7[2])

pca_Fig7_df <- data.frame(
  specimen = rownames(pca_Fig7$x),
  PC1      = pca_Fig7$x[, 1],
  PC2      = pca_Fig7$x[, 2],
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
  x = "PC1",
  y = "PC2",
  min_n = 3
)

p_fig7 <- ggplot(pca_Fig7_df, aes(PC1, PC2)) +
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
  scale_fill_manual(values = group_palette,  labels = group_labels, name = "Group") +
  labs(
    title = "PCA of size-corrected lateral shape (residuals ~ logCsize): CT time series + 1950 habitats",
    subtitle = "Residuals from procD.lm(shape ~ logCsize) computed on the combined subset",
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

print(p_fig7)

ggsave(file.path(fig_dir, "Fig7B_PCA_CTtimeseries_plus_1950habitats_groupHulls_RESID_logCsize.pdf"),
       p_fig7, width = 7, height = 5, units = "in")
ggsave(file.path(fig_dir, "Fig7B_PCA_CTtimeseries_plus_1950habitats_groupHulls_RESID_logCsize.png"),
       p_fig7, width = 7, height = 5, units = "in", dpi = 600)

message("Saved Fig7A + Fig7B to: ", normalizePath(fig_dir, winslash = "/"))