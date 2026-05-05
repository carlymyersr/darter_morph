# ============================================================
# scripts/1950_TPS_raw_synthetic.R
# RAW TPS (1950):
#   (A) synthetic shapes along PCA axes (±k SD)
#   (B) optional extreme individuals on RAW PCA
#   (C) habitat median shapes (RAW) relative to single global ref (ref_1950)
#   - Uses coords_1950 (GPA shapes, not size-corrected residuals)
#   - Saves PNG (+ optional PDF) to figures/TPS_raw_1950/
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
TPS_DIR   <- file.path("figures", "TPS_raw_1950")
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
DO_HABITAT_MEDIANS     <- TRUE

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

# Build synthetic shape along a PC axis at +/- k SD
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

# Habitat medians helper (define only if not already defined somewhere)
if (!exists("habitat_median_shapes")) {
  habitat_median_shapes <- function(coords_arr, groups) {
    stopifnot(length(dim(coords_arr)) == 3)
    stopifnot(dim(coords_arr)[2] == 2)
    stopifnot(length(groups) == dim(coords_arr)[3])
    
    groups <- factor(groups)
    levs <- levels(groups)
    
    out <- setNames(vector("list", length(levs)), levs)
    
    for (g in levs) {
      idx <- which(groups == g)
      if (length(idx) == 0) next
      M <- apply(coords_arr[, , idx, drop = FALSE], c(1, 2), median, na.rm = TRUE)
      colnames(M) <- c("x", "y")
      rownames(M) <- dimnames(coords_arr)[[1]]
      out[[g]] <- M
    }
    out
  }
}

# ============================================================
# RAW PCA (1950)
# ============================================================

pca_raw_1950 <- gm.prcomp(coords_1950)
pct_raw_1950 <- 100 * (pca_raw_1950$sdev^2 / sum(pca_raw_1950$sdev^2))

# Basic sanity
stopifnot(is.matrix(ref_1950), ncol(ref_1950) == 2)
stopifnot(identical(rownames(ref_1950), dimnames(coords_1950)[[1]]))
stopifnot(identical(dimnames(coords_1950)[[3]], gdf_1950$specimen))

# ============================================================
# (A) Synthetic TPS: PC axis shapes at +/- k SD
# ============================================================

pc_axis_targets_raw <- list(
  PC1_minus = make_pc_axis_shape(coords_1950, pca_raw_1950, pc_i = 1, sign = -1, k_sd = PC_K),
  PC1_plus  = make_pc_axis_shape(coords_1950, pca_raw_1950, pc_i = 1, sign =  1, k_sd = PC_K),
  PC2_minus = make_pc_axis_shape(coords_1950, pca_raw_1950, pc_i = 2, sign = -1, k_sd = PC_K),
  PC2_plus  = make_pc_axis_shape(coords_1950, pca_raw_1950, pc_i = 2, sign =  1, k_sd = PC_K)
)

for (nm in names(pc_axis_targets_raw)) {
  pct  <- if (grepl("^PC1", nm)) pct_raw_1950[1] else pct_raw_1950[2]
  axis <- if (grepl("^PC1", nm)) "PC1" else "PC2"
  dirn <- if (grepl("minus$", nm)) "-" else "+"
  
  title <- paste0("RAW PCA axis: ", axis, " ", dirn, PC_K, " SD (", sprintf("%.1f", pct), "%)")
  
  fn_png <- file.path(TPS_DIR, paste0("TPS_RAW_axis_", nm, "_", PC_K, "SD.png"))
  save_tps_png(ref_1950, pc_axis_targets_raw[[nm]], links_1950, fn_png, title)
  
  if (isTRUE(SAVE_PDF)) {
    fn_pdf <- file.path(TPS_DIR, paste0("TPS_RAW_axis_", nm, "_", PC_K, "SD.pdf"))
    save_tps_pdf(ref_1950, pc_axis_targets_raw[[nm]], links_1950, fn_pdf, title)
  }
}

# ============================================================
# (B) Optional: extreme INDIVIDUAL TPS on RAW PCA
# ============================================================

if (isTRUE(DO_EXTREME_INDIVIDUALS)) {
  pc_scores <- pca_raw_1950$x
  pc_extreme_idx <- list(
    PC1_min = which.min(pc_scores[, 1]),
    PC1_max = which.max(pc_scores[, 1]),
    PC2_min = which.min(pc_scores[, 2]),
    PC2_max = which.max(pc_scores[, 2])
  )
  
  pc_extreme_shapes <- lapply(pc_extreme_idx, function(i) coords_1950[, , i, drop = TRUE])
  
  for (nm in names(pc_extreme_shapes)) {
    pct  <- if (grepl("^PC1", nm)) pct_raw_1950[1] else pct_raw_1950[2]
    axis <- if (grepl("^PC1", nm)) "PC1" else "PC2"
    ext  <- if (grepl("min$", nm)) "min" else "max"
    
    title <- paste0("RAW PCA: ", axis, " ", ext, " individual (", sprintf("%.1f", pct), "%)")
    
    fn_png <- file.path(TPS_DIR, paste0("TPS_RAW_extremeIndividual_", nm, ".png"))
    save_tps_png(ref_1950, pc_extreme_shapes[[nm]], links_1950, fn_png, title)
    
    if (isTRUE(SAVE_PDF)) {
      fn_pdf <- file.path(TPS_DIR, paste0("TPS_RAW_extremeIndividual_", nm, ".pdf"))
      save_tps_pdf(ref_1950, pc_extreme_shapes[[nm]], links_1950, fn_pdf, title)
    }
  }
}

# ============================================================
# (C) RAW habitat medians (single global ref = ref_1950)
# ============================================================

if (isTRUE(DO_HABITAT_MEDIANS)) {
  hab_vec <- factor(gdf_1950$habitat)
  
  hab_median_raw <- habitat_median_shapes(coords_1950, hab_vec)
  
  for (h in names(hab_median_raw)) {
    if (is.null(hab_median_raw[[h]])) next
    
    h_safe <- gsub("[^A-Za-z0-9]+", "_", h)
    title  <- paste0("RAW habitat median: ", h, " (1950)")
    
    fn_png <- file.path(TPS_DIR, paste0("TPS_RAW_habMedian_", h_safe, ".png"))
    save_tps_png(ref_1950, hab_median_raw[[h]], links_1950, fn_png, title, mag = 3)
    
    if (isTRUE(SAVE_PDF)) {
      fn_pdf <- file.path(TPS_DIR, paste0("TPS_RAW_habMedian_", h_safe, ".pdf"))
      save_tps_pdf(ref_1950, hab_median_raw[[h]], links_1950, fn_pdf, title, mag = 3)
    }
  }
}

message("Saved RAW TPS outputs to: ", normalizePath(TPS_DIR))