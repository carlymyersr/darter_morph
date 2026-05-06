# ============================================================
# scripts/1950_TPS_residual_synthetic.R
# RESIDUAL TPS (1950):
#   (A) synthetic shapes along PCA axes (±k SD) in residual space, shown on RAW mean
#   (B) optional extreme individuals on RESIDUAL PCA, shown on RAW mean
#   (C) habitat MEAN shapes (RESIDUAL), shown on RAW mean   <-- UPDATED (was median)
#
#   - Uses coords_resid_1950 from R/02_subset_1950.R
#   - Computes gm.prcomp(coords_resid_1950)
#   - Saves PNG (+ optional PDF) to figures/TPS_resid_1950/
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
TPS_DIR   <- file.path("figures", "TPS_resid_1950")
DIR_OK    <- dir.exists(TPS_DIR) || dir.create(TPS_DIR, recursive = TRUE)
stopifnot(DIR_OK)

PC_K      <- 2      # +/- k SD along each PC
TPS_W     <- 4
TPS_H     <- 4
TPS_RES   <- 600
TPS_MAG   <- 2
TPS_PTCEX <- 0.7
TPS_LWD   <- 0.6

SAVE_PDF  <- FALSE     # optional
DO_EXTREME_INDIVIDUALS <- TRUE
DO_HABITAT_MEANS       <- TRUE   # <-- UPDATED name (was DO_HABITAT_MEDIANS)

# ============================================================
# Preconditions / sanity
# ============================================================
if (!exists("coords_resid_1950")) stop("coords_resid_1950 not found. Ensure R/02_subset_1950.R builds it.")
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

# Build synthetic shape along a PC axis at +/- k SD (in whatever space coords_arr lives in)
make_pc_axis_shape <- function(coords_arr, pca_obj, pc_i = 1, sign = c(-1, 1), k_sd = 2) {
  sign <- as.numeric(sign)
  
  X <- geomorph::two.d.array(coords_arr)   # n x (p*k)
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

# Convert a residual-space target into a RAW-looking target by applying residual deformation to RAW ref
resid_target_on_raw <- function(target_resid, ref_resid, ref_raw) {
  stopifnot(all(dim(target_resid) == dim(ref_resid)))
  stopifnot(all(dim(ref_raw) == dim(ref_resid)))
  ref_raw + (target_resid - ref_resid)
}

# Habitat MEANS helper (define only if not already defined somewhere)
if (!exists("habitat_mean_shapes")) {
  habitat_mean_shapes <- function(coords_arr, groups) {
    stopifnot(length(dim(coords_arr)) == 3)
    stopifnot(dim(coords_arr)[2] == 2)
    stopifnot(length(groups) == dim(coords_arr)[3])
    
    groups <- factor(groups)
    levs <- levels(groups)
    
    out <- setNames(vector("list", length(levs)), levs)
    
    for (g in levs) {
      idx <- which(groups == g)
      if (length(idx) == 0) next
      M <- apply(coords_arr[, , idx, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
      colnames(M) <- c("x", "y")
      rownames(M) <- dimnames(coords_arr)[[1]]
      out[[g]] <- M
    }
    out
  }
}

# ============================================================
# RESIDUAL PCA (1950)
# ============================================================

pca_resid_1950 <- gm.prcomp(coords_resid_1950)
pct_resid_1950 <- 100 * (pca_resid_1950$sdev^2 / sum(pca_resid_1950$sdev^2))

# ============================================================
# (A) Synthetic TPS: PC axis shapes at +/- k SD (resid), shown on RAW mean
# ============================================================

pc_axis_targets_resid <- list(
  PC1_minus = make_pc_axis_shape(coords_resid_1950, pca_resid_1950, pc_i = 1, sign = -1, k_sd = PC_K),
  PC1_plus  = make_pc_axis_shape(coords_resid_1950, pca_resid_1950, pc_i = 1, sign =  1, k_sd = PC_K),
  PC2_minus = make_pc_axis_shape(coords_resid_1950, pca_resid_1950, pc_i = 2, sign = -1, k_sd = PC_K),
  PC2_plus  = make_pc_axis_shape(coords_resid_1950, pca_resid_1950, pc_i = 2, sign =  1, k_sd = PC_K)
)

for (nm in names(pc_axis_targets_resid)) {
  pct  <- if (grepl("^PC1", nm)) pct_resid_1950[1] else pct_resid_1950[2]
  axis <- if (grepl("^PC1", nm)) "PC1" else "PC2"
  dirn <- if (grepl("minus$", nm)) "-" else "+"
  
  title <- paste0(
    "SIZE-CORRECTED PCA axis: ", axis, " ", dirn, PC_K, " SD (", sprintf("%.1f", pct), "%)\n",
    "Displayed on RAW mean (residual deformation applied)"
  )
  
  target_on_raw <- resid_target_on_raw(pc_axis_targets_resid[[nm]], ref_resid_1950, ref_1950)
  
  fn_png <- file.path(TPS_DIR, paste0("TPS_RESID_axisOnRaw_", nm, "_", PC_K, "SD.png"))
  save_tps_png(ref_1950, target_on_raw, links_1950, fn_png, title)
  
  if (isTRUE(SAVE_PDF)) {
    fn_pdf <- file.path(TPS_DIR, paste0("TPS_RESID_axisOnRaw_", nm, "_", PC_K, "SD.pdf"))
    save_tps_pdf(ref_1950, target_on_raw, links_1950, fn_pdf, title)
  }
}

# ============================================================
# (B) Optional: extreme INDIVIDUAL TPS on RESIDUAL PCA (shown on RAW mean)
# ============================================================

if (isTRUE(DO_EXTREME_INDIVIDUALS)) {
  pc_scores <- pca_resid_1950$x
  
  pc_extreme_idx <- list(
    PC1_min = which.min(pc_scores[, 1]),
    PC1_max = which.max(pc_scores[, 1]),
    PC2_min = which.min(pc_scores[, 2]),
    PC2_max = which.max(pc_scores[, 2])
  )
  
  # Extreme individuals in residual space
  pc_extreme_shapes_resid <- lapply(pc_extreme_idx, function(i) coords_resid_1950[, , i, drop = TRUE])
  
  for (nm in names(pc_extreme_shapes_resid)) {
    pct  <- if (grepl("^PC1", nm)) pct_resid_1950[1] else pct_resid_1950[2]
    axis <- if (grepl("^PC1", nm)) "PC1" else "PC2"
    ext  <- if (grepl("min$", nm)) "min" else "max"
    
    title <- paste0(
      "SIZE-CORRECTED PCA: ", axis, " ", ext, " individual (", sprintf("%.1f", pct), "%)\n",
      "Displayed on RAW mean (residual deformation applied)"
    )
    
    target_on_raw <- resid_target_on_raw(pc_extreme_shapes_resid[[nm]], ref_resid_1950, ref_1950)
    
    fn_png <- file.path(TPS_DIR, paste0("TPS_RESID_extremeIndividualOnRaw_", nm, ".png"))
    save_tps_png(ref_1950, target_on_raw, links_1950, fn_png, title)
    
    if (isTRUE(SAVE_PDF)) {
      fn_pdf <- file.path(TPS_DIR, paste0("TPS_RESID_extremeIndividualOnRaw_", nm, ".pdf"))
      save_tps_pdf(ref_1950, target_on_raw, links_1950, fn_pdf, title)
    }
  }
}

# ============================================================
# (C) RESIDUAL habitat MEANS (shown on RAW mean)
# ============================================================

if (isTRUE(DO_HABITAT_MEANS)) {
  hab_vec <- factor(gdf_1950$habitat)
  
  hab_mean_resid <- habitat_mean_shapes(coords_resid_1950, hab_vec)
  
  for (h in names(hab_mean_resid)) {
    if (is.null(hab_mean_resid[[h]])) next
    
    h_safe <- gsub("[^A-Za-z0-9]+", "_", h)
    title  <- paste0(
      "SIZE-CORRECTED habitat mean: ", h, " (1950)\n",
      "Displayed on RAW mean (residual deformation applied)"
    )
    
    target_on_raw <- resid_target_on_raw(hab_mean_resid[[h]], ref_resid_1950, ref_1950)
    
    fn_png <- file.path(TPS_DIR, paste0("TPS_RESID_habMeanOnRaw_", h_safe, ".png"))
    save_tps_png(ref_1950, target_on_raw, links_1950, fn_png, title, mag = 2)
    
    if (isTRUE(SAVE_PDF)) {
      fn_pdf <- file.path(TPS_DIR, paste0("TPS_RESID_habMeanOnRaw_", h_safe, ".pdf"))
      save_tps_pdf(ref_1950, target_on_raw, links_1950, fn_pdf, title, mag = 2)
    }
  }
}

message("Saved RESIDUAL TPS outputs to: ", normalizePath(TPS_DIR))