# ============================================================
# scripts/1950_TPS_residual_PC23_extremes.R
# RESIDUAL TPS (1950): PC2 and PC3 synthetic axis shapes + extreme individuals
#
#   - Uses coords_resid_1950 from R/02_subset_1950.R
#   - Computes gm.prcomp(coords_resid_1950)
#   - Builds synthetic shapes at +/- k SD on PC2 and PC3
#   - Finds extreme individuals on PC2 and PC3
#   - Displays residual deformations on RAW mean shape
#   - Saves PNG (+ optional PDF) to figures/1950_TPS_RESID_PC23/
# ============================================================

# ---------------------------
# Load canonical project objects
# ---------------------------
source("R/00_setup_morpho.R")
source("R/01_build_metadata.R")
source("R/02_subset_1950.R")

suppressPackageStartupMessages({
  library(geomorph)
})

# ---------------------------
# User settings
# ---------------------------
TPS_DIR   <- file.path("figures", "1950_TPS_RESID_PC23")
DIR_OK    <- dir.exists(TPS_DIR) || dir.create(TPS_DIR, recursive = TRUE)
stopifnot(DIR_OK)

PC_K      <- 2
TPS_W     <- 4
TPS_H     <- 4
TPS_RES   <- 600
TPS_MAG   <- 2
TPS_PTCEX <- 0.7
TPS_LWD   <- 0.6

SAVE_PDF  <- FALSE

# ============================================================
# Preconditions / sanity
# ============================================================
if (!exists("coords_resid_1950")) {
  stop("coords_resid_1950 not found. Ensure R/02_subset_1950.R builds it.")
}
stopifnot(identical(dimnames(coords_resid_1950), dimnames(coords_1950)))
stopifnot(is.matrix(ref_1950), ncol(ref_1950) == 2)
stopifnot(identical(dimnames(coords_1950)[[3]], gdf_1950$specimen))

# Residual-space reference (mean residual shape)
ref_resid_1950 <- mshape(coords_resid_1950)
ref_resid_1950 <- matrix(
  ref_resid_1950,
  ncol = 2,
  dimnames = list(dimnames(coords_resid_1950)[[1]], c("x", "y"))
)
stopifnot(identical(rownames(ref_resid_1950), rownames(ref_1950)))

# ============================================================
# Helpers
# ============================================================

save_tps_png <- function(ref, target, links, filename, title,
                         mag = TPS_MAG, width = TPS_W, height = TPS_H, res = TPS_RES,
                         pt_cex = TPS_PTCEX, lwd = TPS_LWD) {
  stopifnot(is.matrix(ref), is.matrix(target))
  stopifnot(all(dim(ref) == dim(target)), ncol(ref) == 2)
  
  png(filename, width = width, height = height, units = "in", res = res)
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  par(mar = c(2, 2, 3, 1), cex = pt_cex, lwd = lwd)
  
  geomorph::plotRefToTarget(
    ref, target,
    method = "TPS",
    mag    = mag,
    links  = links,
    main   = title
  )
}

save_tps_pdf <- function(ref, target, links, filename, title,
                         mag = TPS_MAG, width = TPS_W, height = TPS_H,
                         pt_cex = TPS_PTCEX, lwd = TPS_LWD) {
  stopifnot(is.matrix(ref), is.matrix(target))
  stopifnot(all(dim(ref) == dim(target)), ncol(ref) == 2)
  
  grDevices::cairo_pdf(filename, width = width, height = height)
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  par(mar = c(2, 2, 3, 1), cex = pt_cex, lwd = lwd)
  
  geomorph::plotRefToTarget(
    ref, target,
    method = "TPS",
    mag    = mag,
    links  = links,
    main   = title
  )
}

# Convert a residual-space target into a RAW-looking target by applying
# residual deformation to RAW ref
resid_target_on_raw <- function(target_resid, ref_resid, ref_raw) {
  stopifnot(all(dim(target_resid) == dim(ref_resid)))
  stopifnot(all(dim(ref_raw) == dim(ref_resid)))
  ref_raw + (target_resid - ref_resid)
}

# Build synthetic shape along a PC axis at +/- k SD
make_pc_axis_shape <- function(coords_arr, pca_obj, pc_i = 1, sign = c(-1, 1), k_sd = 2) {
  sign <- as.numeric(sign)
  
  X <- geomorph::two.d.array(coords_arr)
  X_mean <- colMeans(X)
  
  rot <- pca_obj$rotation
  if (is.null(rot) && !is.null(pca_obj$vectors)) rot <- pca_obj$vectors
  if (is.null(rot)) stop("Could not find PCA loadings matrix in pca object (rotation/vectors).")
  
  pc_sd <- apply(pca_obj$x, 2, sd)
  
  score_vec <- rep(0, ncol(rot))
  score_vec[pc_i] <- sign * k_sd * pc_sd[pc_i]
  
  shape_vec <- X_mean + as.numeric(rot %*% score_vec)
  
  p <- dim(coords_arr)[1]
  k <- dim(coords_arr)[2]
  
  arr <- geomorph::arrayspecs(matrix(shape_vec, nrow = 1), p = p, k = k)
  M <- arr[, , 1, drop = TRUE]
  rownames(M) <- dimnames(coords_arr)[[1]]
  colnames(M) <- c("x", "y")
  M
}

# ============================================================
# RESIDUAL PCA (1950)
# ============================================================
pca_resid_1950 <- gm.prcomp(coords_resid_1950)
pct_resid_1950 <- 100 * (pca_resid_1950$sdev^2 / sum(pca_resid_1950$sdev^2))
pc_scores <- pca_resid_1950$x

if (ncol(pc_scores) < 3) {
  stop("PCA has fewer than 3 axes; cannot generate PC3 plots.")
}

# ============================================================
# Synthetic axis shapes on PC2 and PC3
# ============================================================
pc_axis_targets_resid <- list(
  PC2_minus = make_pc_axis_shape(coords_resid_1950, pca_resid_1950, pc_i = 2, sign = -1, k_sd = PC_K),
  PC2_plus  = make_pc_axis_shape(coords_resid_1950, pca_resid_1950, pc_i = 2, sign =  1, k_sd = PC_K),
  PC3_minus = make_pc_axis_shape(coords_resid_1950, pca_resid_1950, pc_i = 3, sign = -1, k_sd = PC_K),
  PC3_plus  = make_pc_axis_shape(coords_resid_1950, pca_resid_1950, pc_i = 3, sign =  1, k_sd = PC_K)
)

for (nm in names(pc_axis_targets_resid)) {
  axis_num <- if (grepl("^PC2", nm)) 2 else 3
  pct      <- pct_resid_1950[axis_num]
  axis_lab <- paste0("PC", axis_num)
  dirn     <- if (grepl("minus$", nm)) "-" else "+"
  
  title <- paste0(
    "SIZE-CORRECTED PCA axis: ", axis_lab, " ", dirn, PC_K, " SD (",
    sprintf("%.1f", pct), "%)\n",
    "Displayed on RAW mean (residual deformation applied)"
  )
  
  target_on_raw <- resid_target_on_raw(
    pc_axis_targets_resid[[nm]],
    ref_resid_1950,
    ref_1950
  )
  
  fn_png <- file.path(TPS_DIR, paste0("TPS_RESID_axisOnRaw_", nm, "_", PC_K, "SD.png"))
  save_tps_png(ref_1950, target_on_raw, links_1950, fn_png, title)
  
  if (isTRUE(SAVE_PDF)) {
    fn_pdf <- file.path(TPS_DIR, paste0("TPS_RESID_axisOnRaw_", nm, "_", PC_K, "SD.pdf"))
    save_tps_pdf(ref_1950, target_on_raw, links_1950, fn_pdf, title)
  }
}

# ============================================================
# Extreme individuals on PC2 and PC3
# ============================================================
pc_extreme_idx <- list(
  PC2_min = which.min(pc_scores[, 2]),
  PC2_max = which.max(pc_scores[, 2]),
  PC3_min = which.min(pc_scores[, 3]),
  PC3_max = which.max(pc_scores[, 3])
)

pc_extreme_shapes_resid <- lapply(
  pc_extreme_idx,
  function(i) coords_resid_1950[, , i, drop = TRUE]
)

for (nm in names(pc_extreme_shapes_resid)) {
  axis_num <- if (grepl("^PC2", nm)) 2 else 3
  pct      <- pct_resid_1950[axis_num]
  axis_lab <- paste0("PC", axis_num)
  ext_lab  <- if (grepl("min$", nm)) "min" else "max"
  
  specimen_id <- dimnames(coords_resid_1950)[[3]][pc_extreme_idx[[nm]]]
  
  title <- paste0(
    "SIZE-CORRECTED PCA: ", axis_lab, " ", ext_lab, " individual (",
    sprintf("%.1f", pct), "%)\n",
    specimen_id, "\n",
    "Displayed on RAW mean (residual deformation applied)"
  )
  
  target_on_raw <- resid_target_on_raw(
    pc_extreme_shapes_resid[[nm]],
    ref_resid_1950,
    ref_1950
  )
  
  fn_png <- file.path(TPS_DIR, paste0("TPS_RESID_extremeIndividualOnRaw_", nm, ".png"))
  save_tps_png(ref_1950, target_on_raw, links_1950, fn_png, title)
  
  if (isTRUE(SAVE_PDF)) {
    fn_pdf <- file.path(TPS_DIR, paste0("TPS_RESID_extremeIndividualOnRaw_", nm, ".pdf"))
    save_tps_pdf(ref_1950, target_on_raw, links_1950, fn_pdf, title)
  }
}

message("Saved residual PC2/PC3 synthetic + extreme TPS outputs to: ", normalizePath(TPS_DIR))