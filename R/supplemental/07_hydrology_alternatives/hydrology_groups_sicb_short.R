# ============================================================
# scripts/mainstem_alltime_vs_tribandres_SICB.R
# PCA on size-corrected shapes (residuals)
# Variation-source grouping:
#   1) Mainstem = Connecticut River (1950, 1956, 1970)
#   2) Reservoir System = Quabbin + Swift River
#   3) Tributaries = Sawmill River + Fort River
# ============================================================

source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")

suppressPackageStartupMessages({
  library(geomorph)
  library(ggplot2)
  library(dplyr)
})

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

FIG_DIR <- file.path("Figures", "variation_source_SICB")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_DIR <- file.path("Outputs", "mainstem_alltime_vs_tribandres")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("coords_gpa")) stop("coords_gpa not found.")
if (!exists("gdf")) stop("gdf not found.")
if (!exists("subset_coords_to_gdf")) stop("subset_coords_to_gdf() not found.")
if (!exists("allometry_residuals")) stop("allometry_residuals() not found.")

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

coords_vs_raw <- subset_coords_to_gdf(coords_gpa, gdf_vs)
stopifnot(identical(dimnames(coords_vs_raw)[[3]], gdf_vs$specimen))

size_vec <- setNames(gdf_vs$size_for_allometry, gdf_vs$specimen)
allo_vs <- allometry_residuals(coords_vs_raw, size_vec)
coords_vs <- allo_vs$residuals

stopifnot(identical(dimnames(coords_vs)[[3]], gdf_vs$specimen))

# New palette — intentionally different from watershed landscape figure
vs_palette <- c(
  "Mainstem" = "#4b6043",         
  "Reservoir System" = "#028A0F",  
  "Tributaries" = "#98BF64"        
)

fit_procD_alt <- procD.lm(
  f1 = coords_vs ~ variation_source_alt,
  data = gdf_vs,
  iter = 999,
  RRPP = TRUE
)

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
saveRDS(allo_vs, file = file.path(OUT_DIR, "allometry_fit_combined_subset.rds"))

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

# ============================================================
# Journal-size PCA plot helper: tuned for 3 x 2.36 inch PDFs
# ============================================================

make_pca_plot <- function(df, xvar, yvar, xlab, ylab) {
  
  hull_df <- make_hulls(
    df = df,
    group_col = "variation_source_alt",
    x = xvar,
    y = yvar
  )
  
  ggplot(df, aes(x = .data[[xvar]], y = .data[[yvar]])) +
    geom_polygon(
      data = hull_df,
      aes(group = variation_source_alt, fill = variation_source_alt),
      alpha = 0.13,
      color = NA
    ) +
    geom_polygon(
      data = hull_df,
      aes(group = variation_source_alt, color = variation_source_alt),
      fill = NA,
      linewidth = 0.22
    ) +
    geom_point(
      aes(color = variation_source_alt),
      size = 0.75,
      alpha = 0.85
    ) +
    scale_color_manual(values = vs_palette, guide = "none") +
    scale_fill_manual(values = vs_palette, guide = "none") +
    labs(
      x = xlab,
      y = ylab
    ) +
    coord_equal() +
    theme_classic(base_family = "Helvetica", base_size = 6) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      legend.position = "none",
      
      axis.title = element_text(size = 6),
      axis.text  = element_text(size = 5),
      
      axis.line = element_line(linewidth = 0.22),
      axis.ticks = element_line(linewidth = 0.22),
      axis.ticks.length = unit(1.2, "mm"),
      
      plot.margin = margin(2, 2, 2, 2)
    )
}

p_pc12 <- make_pca_plot(pca_vs_df, "PC1", "PC2", pc1_lab, pc2_lab)
p_pc23 <- make_pca_plot(pca_vs_df, "PC2", "PC3", pc2_lab, pc3_lab)
p_pc34 <- make_pca_plot(pca_vs_df, "PC3", "PC4", pc3_lab, pc4_lab)

ggsave(
  file.path(FIG_DIR, "Fig_mainstem_alltime_vs_tribandres_PC12_SICB_3in.pdf"),
  p_pc12,
  width = 3,
  height = 2.36,
  units = "in",
  device = cairo_pdf
)

ggsave(
  file.path(FIG_DIR, "Fig_mainstem_alltime_vs_tribandres_PC23_SICB_3in.pdf"),
  p_pc23,
  width = 3,
  height = 2.36,
  units = "in",
  device = cairo_pdf
)

ggsave(
  file.path(FIG_DIR, "Fig_mainstem_alltime_vs_tribandres_PC34_SICB_3in.pdf"),
  p_pc34,
  width = 3,
  height = 2.36,
  units = "in",
  device = cairo_pdf
)

cat("\nSaved SICB figures to:\n")
cat("  ", normalizePath(FIG_DIR, winslash = "/"), "\n")