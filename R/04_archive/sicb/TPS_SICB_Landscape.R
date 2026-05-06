# ============================================================
# Scripts/TPS_SICB_Landscape.R
# TPS plots for synthetic PC extremes from the complete landscape
#   Complete landscape = CT time series + 1950 habitats
#   Synthetic extremes for PC1, PC2, PC3, and PC4
#
# Formatting follows the TPS settings/helpers from 1950_figures_SICB.R
# PCA/residual objects follow watershed_landscape_SICB.R
#
# Saves TPS plots to:
#   Figures/TPS_SICB_Landscape/
# ============================================================

suppressPackageStartupMessages({
  library(geomorph)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(showtext)

  # Match the 1950 SICB TPS script formatting as closely as possible.
  # If Arial is unavailable on another machine, the script continues.
  try({
    font_add("Arial", "/System/Library/Fonts/Supplemental/Arial.ttf")
    showtext_auto()
  }, silent = TRUE)
})

# ---------------------------
# Load canonical data
# ---------------------------
source("R/00_setup_morpho.R")
source("R/01_build_metadata.R")
source("R/04_subset_CT_timeseries_plus_1950habitats.R")

# ---------------------------
# Output directory
# ---------------------------
OUTDIR <- file.path("Figures", "TPS_SICB_Landscape")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Preconditions
# ============================================================

if (!exists("coords_resid_Fig6")) {
  stop("coords_resid_Fig6 not found. Ensure R/04_subset_CT_timeseries_plus_1950habitats.R builds it.")
}

if (!exists("gdf_Fig6")) {
  stop("gdf_Fig6 not found. Ensure R/04_subset_CT_timeseries_plus_1950habitats.R builds it.")
}

# The complete-landscape raw coordinates are usually coords_Fig6 in this workflow.
# If your upstream script uses a different name, add it here.
if (exists("coords_Fig6")) {
  coords_raw_landscape <- coords_Fig6
} else if (exists("coords_raw_Fig6")) {
  coords_raw_landscape <- coords_raw_Fig6
} else {
  stop("Raw complete-landscape coordinates not found. Expected coords_Fig6 or coords_raw_Fig6 from R/04_subset_CT_timeseries_plus_1950habitats.R.")
}

stopifnot(identical(dimnames(coords_resid_Fig6), dimnames(coords_raw_landscape)))
stopifnot(identical(dimnames(coords_resid_Fig6)[[3]], gdf_Fig6$specimen))

# ============================================================
# PCA on residual coordinates for complete landscape
# ============================================================

pca_landscape <- geomorph::gm.prcomp(coords_resid_Fig6)
pct_landscape <- 100 * (pca_landscape$sdev^2 / sum(pca_landscape$sdev^2))

message("Percent variance explained:")
for (pc in 1:4) {
  message("  PC", pc, ": ", round(pct_landscape[pc], 2), "%")
}

# ============================================================
# TPS SETTINGS — copied from 1950_figures_SICB.R
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

# ============================================================
# Reference shape: RAW complete-landscape mean + matched RESID mean
# ============================================================

ref_raw <- geomorph::mshape(coords_raw_landscape)
ref_raw <- matrix(
  ref_raw,
  ncol = 2,
  dimnames = list(dimnames(coords_raw_landscape)[[1]], c("x", "y"))
)

ref_resid <- geomorph::mshape(coords_resid_Fig6)
ref_resid <- matrix(
  ref_resid,
  ncol = 2,
  dimnames = list(dimnames(coords_resid_Fig6)[[1]], c("x", "y"))
)

stopifnot(identical(rownames(ref_raw), rownames(ref_resid)))

# ============================================================
# Helpers
# ============================================================

# Prevent .png.png / .png.pdf
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
  # This matches the SICB TPS style in 1950_figures_SICB.R.
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

  X <- geomorph::two.d.array(coords_resid_Fig6)
  mean_vec <- colMeans(X)

  rot <- pca_landscape$rotation
  if (is.null(rot) && !is.null(pca_landscape$vectors)) rot <- pca_landscape$vectors
  if (is.null(rot)) stop("Could not find PCA loadings matrix in pca_landscape object (rotation/vectors).")

  sd_pc <- sd(pca_landscape$x[, pc])

  score_vec <- rep(0, ncol(rot))
  score_vec[pc] <- sign * k * sd_pc

  shape_vec <- mean_vec + as.numeric(rot %*% score_vec)

  M <- geomorph::arrayspecs(
    matrix(shape_vec, nrow = 1),
    p = dim(coords_resid_Fig6)[1],
    k = 2
  )[, , 1]

  rownames(M) <- dimnames(coords_resid_Fig6)[[1]]
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
# TPS: synthetic PC1 / PC2 / PC3 / PC4 extremes
# ============================================================

for (pc in 1:4) {

  draw_tps(
    apply_resid(make_shape(pc, -1)),
    file.path(OUTDIR, paste0("TPS_PC", pc, "_minus"))
  )

  draw_tps(
    apply_resid(make_shape(pc,  1)),
    file.path(OUTDIR, paste0("TPS_PC", pc, "_plus"))
  )
}

message("Complete-landscape SICB TPS plots complete. Saved to: ", normalizePath(OUTDIR, winslash = "/"))
