# ============================================================
# CT3_figures_SICB.R
# Clean script for:
#   - PCA plot (PC1 vs PC2)
#   - TPS plots (scaled for figure panels)
# ============================================================

suppressPackageStartupMessages({
  library(geomorph)
  library(ggplot2)
  library(dplyr)
  library(showtext)
  
  font_add("Arial", "/System/Library/Fonts/Supplemental/Arial.ttf")
  showtext_auto()
})

# ---------------------------
# Load canonical data
# ---------------------------
source("R/00_setup_morpho.R")
source("R/01_build_metadata.R")
source("R/03_subset_CT_timeseries.R")

# ---------------------------
# Output directory
# ---------------------------
OUTDIR <- file.path("Figures", "CT3_figures")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# PCA (RESID) — FINAL CLEAN VERSION
# ============================================================

pca <- gm.prcomp(coords_resid_CT3)
pct <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

df <- data.frame(
  specimen = rownames(pca$x),
  PC1 = -pca$x[,1],   # invert
  PC2 =  pca$x[,2],
  group = gdf_CT3$group,
  stringsAsFactors = FALSE
)

# ---- Hulls (only groups with >=3 points) ----
df_hull <- df %>%
  dplyr::group_by(group) %>%
  dplyr::filter(n() >= 3) %>%
  dplyr::ungroup()

hulls <- make_hulls(
  df = df_hull,
  group_col = "group",
  x = "PC1",
  y = "PC2"
)

# ---- Color palette ----
group_palette <- c(
  "CT_1950" = "steelblue",
  "CT_1956" = "dodgerblue4",
  "CT_1970" = "navy",
  "CT_1979" = "orange"
)

# ---- PCA plot ----
p <- ggplot(df, aes(PC1, PC2)) +
  
  geom_polygon(
    data = hulls,
    aes(group = group, fill = group),
    alpha = 0.15,
    color = NA
  ) +
  
  geom_polygon(
    data = hulls,
    aes(group = group, color = group),
    fill = NA,
    linewidth = 0.25
  ) +
  
  geom_point(
    aes(color = group),
    size = 0.65,
    alpha = 0.85
  ) +
  
  scale_color_manual(values = group_palette, name = NULL) +
  scale_fill_manual(values = group_palette, name = NULL) +
  
  coord_equal() +
  
  theme_classic(base_size = 6, base_family = "Arial") +
  theme(
    axis.title = element_text(size = 6),
    axis.text  = element_text(size = 5),
    
    axis.line = element_line(linewidth = 0.25),
    axis.ticks = element_line(linewidth = 0.25),
    axis.ticks.length = unit(1.5, "pt"),
    
    legend.position = "right",
    legend.text = element_text(size = 5),
    legend.key.size = unit(0.18, "in"),
    legend.spacing.y = unit(0.02, "in"),
    
    plot.margin = margin(2, 2, 2, 2, unit = "pt")
  ) +
  
  guides(
    color = guide_legend(
      override.aes = list(size = 1.2, alpha = 0.85)
    ),
    fill = guide_legend(
      override.aes = list(alpha = 0.15)
    )
  ) +
  
  labs(
    x = paste0("-PC1 (", round(pct[1], 1), "%)"),
    y = paste0("PC2 (", round(pct[2], 1), "%)")
  )

# ---- Save ----
ggsave(
  file.path(OUTDIR, "PCA_CT3_PC1_PC2_hulls.png"),
  p,
  width = 3,
  height = 2.4,
  units = "in",
  dpi = 1200,
  bg = "white"
)

ggsave(
  file.path(OUTDIR, "PCA_CT3_PC1_PC2_hulls.pdf"),
  p,
  width = 3,
  height = 2.4,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)
# ============================================================
# TPS SETTINGS
# ============================================================

TPS_WIDTH  <- 1.5
TPS_HEIGHT <- 1.5
TPS_RES    <- 600

TPS_GRID_COL <- adjustcolor("black", alpha.f = 0.60)
TPS_PT_COL   <- adjustcolor("black", alpha.f = 0.90)

TPS_GRID_LWD <- 0.30
TPS_PT_CEX   <- 0.6
TPS_MAG      <- 2
TPS_N_COL_CELL <- 10

# ---------------------------
# Reference shape: RAW CT_1950
# ---------------------------
idx_ref <- which(gdf_CT3$group == "CT_1950")

ref_raw <- mshape(coords_CT3[, , idx_ref, drop = FALSE])
ref_resid <- mshape(coords_resid_CT3[, , idx_ref, drop = FALSE])

ref_raw <- matrix(
  ref_raw,
  ncol = 2,
  dimnames = list(dimnames(coords_CT3)[[1]], c("x", "y"))
)

ref_resid <- matrix(
  ref_resid,
  ncol = 2,
  dimnames = list(dimnames(coords_resid_CT3)[[1]], c("x", "y"))
)

# ---------------------------
# Clean filename helper
# prevents .png.png / .png.pdf
# ---------------------------
clean_filebase <- function(x) {
  x <- sub("\\.png$", "", x)
  x <- sub("\\.pdf$", "", x)
  x
}

# ---------------------------
# TPS plotting helper
# ---------------------------

draw_tps_one_device <- function(target, file, device = c("png", "pdf")) {
  
  device <- match.arg(device)
  
  if (device == "png") {
    png(
      filename = file,
      width = TPS_WIDTH,
      height = TPS_HEIGHT,
      units = "in",
      res = TPS_RES,
      bg = "white"
    )
  }
  
  if (device == "pdf") {
    cairo_pdf(
      filename = file,
      width = TPS_WIDTH,
      height = TPS_HEIGHT,
      bg = "white"
    )
  }
  
  on.exit(dev.off(), add = TRUE)
  
  par(
    mar = c(0, 0, 0, 0),
    oma = c(0, 0, 0, 0),
    xaxs = "i",
    yaxs = "i",
    pty = "s",
    col = TPS_GRID_COL,
    fg  = TPS_GRID_COL,
    lwd = TPS_GRID_LWD,
    cex = TPS_PT_CEX
  )
  
  # Draw TPS grid only: no links, no built-in points
  geomorph::plotRefToTarget(
    ref_raw,
    target,
    method = "TPS",
    mag = TPS_MAG,
    links = NULL,
    main = "",
    gridPars = geomorph::gridPar(
      pt.bg = "transparent",
      pt.size = 0,
      tar.pt.bg = "black",
      tar.pt.size = TPS_PT_CEX,
      
      # TPS grid styling
      grid.lwd = TPS_GRID_LWD,
      n.col.cell = TPS_N_COL_CELL
    )
  )
  
  # Draw ONLY target landmarks manually
 # points(
 #   target[, 1],
 #   target[, 2],
 #   pch = 16,
 #   cex = TPS_PT_CEX,
 #   col = "black"
 # )
}

draw_tps <- function(target, filename_base) {
  
  filename_base <- sub("\\.png$", "", filename_base)
  filename_base <- sub("\\.pdf$", "", filename_base)
  
  draw_tps_one_device(
    target = target,
    file = paste0(filename_base, ".png"),
    device = "png"
  )
  
  draw_tps_one_device(
    target = target,
    file = paste0(filename_base, ".pdf"),
    device = "pdf"
  )
}

# ---------------------------
# Helper: PC synthetic shape
# ---------------------------

make_shape <- function(pc, sign = 1, k = 2) {
  
  X <- geomorph::two.d.array(coords_resid_CT3)
  mean_vec <- colMeans(X)
  
  rot <- pca$rotation
  if (is.null(rot) && !is.null(pca$vectors)) rot <- pca$vectors
  
  sd_pc <- sd(pca$x[, pc])
  
  score_vec <- rep(0, ncol(rot))
  score_vec[pc] <- sign * k * sd_pc
  
  shape_vec <- mean_vec + as.numeric(rot %*% score_vec)
  
  M <- geomorph::arrayspecs(
    matrix(shape_vec, nrow = 1),
    p = dim(coords_resid_CT3)[1],
    k = 2
  )[, , 1]
  
  rownames(M) <- dimnames(coords_resid_CT3)[[1]]
  colnames(M) <- c("x", "y")
  
  M
}

# ---------------------------
# Apply residual deformation to RAW reference
# ---------------------------

apply_resid <- function(shape_resid) {
  ref_raw + (shape_resid - ref_resid)
}

# ============================================================
# TPS: synthetic PC1 / PC2 extremes
# ============================================================

draw_tps(
  apply_resid(make_shape(1, -1)),
  file.path(OUTDIR, "TPS_PC1_minus")
)

draw_tps(
  apply_resid(make_shape(1,  1)),
  file.path(OUTDIR, "TPS_PC1_plus")
)

draw_tps(
  apply_resid(make_shape(2, -1)),
  file.path(OUTDIR, "TPS_PC2_minus")
)

draw_tps(
  apply_resid(make_shape(2,  1)),
  file.path(OUTDIR, "TPS_PC2_plus")
)

# ============================================================
# TPS: mean shape per year
# ============================================================

for (g in unique(gdf_CT3$group)) {
  
  idx_g <- which(gdf_CT3$group == g)
  
  mean_resid <- mshape(coords_resid_CT3[, , idx_g, drop = FALSE])
  mean_resid <- matrix(
    mean_resid,
    ncol = 2,
    dimnames = list(dimnames(coords_resid_CT3)[[1]], c("x", "y"))
  )
  
  target <- apply_resid(mean_resid)
  
  draw_tps(
    target,
    file.path(OUTDIR, paste0("TPS_mean_", g))
  )
}

message("TPS plots complete.")