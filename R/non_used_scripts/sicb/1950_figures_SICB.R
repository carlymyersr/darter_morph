# ============================================================
# 1950_figures_SICB.R
# Clean SICB-style script for 1950 habitat figures:
#   - PCA plot (PC1 vs PC2) for 1950 habitats
#   - TPS plots scaled for figure panels:
#       * synthetic PC1 / PC2 extremes
#       * mean shape per habitat
#
# Built by combining the 1950 PCA/TPS logic with the formatting
# used in CT3_figures_SICB.R.
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
source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")
source("R/methods/01_specimen_sampling_study_design/02_subset_1950.R")

# ---------------------------
# Output directory
# ---------------------------
OUTDIR <- file.path("Figures", "1950_figures")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Preconditions
# ============================================================

if (!exists("coords_1950"))       stop("coords_1950 not found. Ensure R/methods/01_specimen_sampling_study_design/02_subset_1950.R builds it.")
if (!exists("coords_resid_1950")) stop("coords_resid_1950 not found. Ensure R/methods/01_specimen_sampling_study_design/02_subset_1950.R builds it.")
if (!exists("gdf_1950"))          stop("gdf_1950 not found. Ensure R/methods/01_specimen_sampling_study_design/02_subset_1950.R builds it.")
if (!exists("make_hulls"))        stop("make_hulls() not found. Keep it in R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R.")
if (!exists("ref_1950"))          stop("ref_1950 not found. It should be created in R/methods/01_specimen_sampling_study_design/02_subset_1950.R.")
if (!exists("links_1950"))        stop("links_1950 not found. It should be created in R/methods/01_specimen_sampling_study_design/02_subset_1950.R.")

stopifnot(identical(dimnames(coords_resid_1950), dimnames(coords_1950)))
stopifnot(identical(dimnames(coords_1950)[[3]], gdf_1950$specimen))

# If hab_palette is not already defined upstream, define a safe fallback.
if (!exists("hab_palette")) {
  hab_levels <- levels(factor(gdf_1950$habitat))
  hab_palette <- setNames(
    grDevices::hcl.colors(length(hab_levels), palette = "Dark 3"),
    hab_levels
  )
}

# ============================================================
# PCA (RESID) — SICB CLEAN VERSION
# ============================================================

pca <- gm.prcomp(coords_resid_1950)
pct <- 100 * (pca$sdev^2 / sum(pca$sdev^2))

df <- data.frame(
  specimen = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  habitat = factor(gdf_1950$habitat),
  stringsAsFactors = FALSE
)

stopifnot(identical(df$specimen, gdf_1950$specimen))

# ---- Hulls (only groups with >=3 points) ----
df_hull <- df %>%
  dplyr::group_by(habitat) %>%
  dplyr::filter(n() >= 3) %>%
  dplyr::ungroup()

hulls <- make_hulls(
  df = df_hull,
  group_col = "habitat",
  x = "PC1",
  y = "PC2"
)

# ---- PCA plot ----
p <- ggplot(df, aes(PC1, PC2)) +

  geom_polygon(
    data = hulls,
    aes(group = habitat, fill = habitat),
    alpha = 0.15,
    color = NA
  ) +

  geom_polygon(
    data = hulls,
    aes(group = habitat, color = habitat),
    fill = NA,
    linewidth = 0.25
  ) +

  geom_point(
    aes(color = habitat),
    size = 0.65,
    alpha = 0.85
  ) +

  scale_color_manual(values = hab_palette, name = NULL) +
  scale_fill_manual(values = hab_palette, name = NULL) +

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
    x = paste0("PC1 (", round(pct[1], 1), "%)"),
    y = paste0("PC2 (", round(pct[2], 1), "%)")
  )

# ---- Save PCA ----
ggsave(
  file.path(OUTDIR, "PCA_1950_habitats_PC1_PC2_hulls.png"),
  p,
  width = 3,
  height = 2.4,
  units = "in",
  dpi = 1200,
  bg = "white"
)

ggsave(
  file.path(OUTDIR, "PCA_1950_habitats_PC1_PC2_hulls.pdf"),
  p,
  width = 3,
  height = 2.4,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

# ============================================================
# TPS SETTINGS — matched to CT3_figures_SICB.R
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
# Reference shape: RAW 1950 mean + matched RESID 1950 mean
# ---------------------------

ref_raw <- matrix(
  ref_1950,
  ncol = 2,
  dimnames = list(dimnames(coords_1950)[[1]], c("x", "y"))
)

ref_resid <- mshape(coords_resid_1950)
ref_resid <- matrix(
  ref_resid,
  ncol = 2,
  dimnames = list(dimnames(coords_resid_1950)[[1]], c("x", "y"))
)

stopifnot(identical(rownames(ref_raw), rownames(ref_resid)))

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

  # Draw TPS grid only: no links, no built-in points.
  # Target points are currently hidden to match CT3_figures_SICB.R.
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

  # If you decide you want visible target landmarks, uncomment this:
  # points(
  #   target[, 1],
  #   target[, 2],
  #   pch = 16,
  #   cex = TPS_PT_CEX,
  #   col = TPS_PT_COL
  # )
}

draw_tps <- function(target, filename_base) {

  filename_base <- clean_filebase(filename_base)

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

  X <- geomorph::two.d.array(coords_resid_1950)
  mean_vec <- colMeans(X)

  rot <- pca$rotation
  if (is.null(rot) && !is.null(pca$vectors)) rot <- pca$vectors
  if (is.null(rot)) stop("Could not find PCA loadings matrix in pca object (rotation/vectors).")

  sd_pc <- sd(pca$x[, pc])

  score_vec <- rep(0, ncol(rot))
  score_vec[pc] <- sign * k * sd_pc

  shape_vec <- mean_vec + as.numeric(rot %*% score_vec)

  M <- geomorph::arrayspecs(
    matrix(shape_vec, nrow = 1),
    p = dim(coords_resid_1950)[1],
    k = 2
  )[, , 1]

  rownames(M) <- dimnames(coords_resid_1950)[[1]]
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
# TPS: mean shape per habitat
# ============================================================

for (h in levels(factor(gdf_1950$habitat))) {

  idx_h <- which(gdf_1950$habitat == h)
  if (length(idx_h) == 0) next

  mean_resid <- mshape(coords_resid_1950[, , idx_h, drop = FALSE])
  mean_resid <- matrix(
    mean_resid,
    ncol = 2,
    dimnames = list(dimnames(coords_resid_1950)[[1]], c("x", "y"))
  )

  target <- apply_resid(mean_resid)

  h_safe <- gsub("[^A-Za-z0-9]+", "_", h)

  draw_tps(
    target,
    file.path(OUTDIR, paste0("TPS_mean_", h_safe))
  )
}

message("1950 SICB-style PCA and TPS plots complete. Saved to: ", normalizePath(OUTDIR))
