# ============================================================
# Scripts/TPS_CT3_residual_synthetic.R
# RESIDUAL TPS (CT time series: 1950, 1956, 1970):
#   (A) synthetic shapes along PCA axes (±k SD) in residual space, shown on RAW reference
#   (B) optional extreme individuals on RESIDUAL PCA, shown on RAW reference
#   (C) year/group median shapes (RESIDUAL), shown on RAW reference
#
# Requires canonical objects created by:
#   R/00_setup_morpho.R           (coords_gpa, links_all, helpers)
#   R/01_build_metadata.R         (gdf, allometry_residuals, etc.)
#   R/03_subset_CT_timeseries.R   (gdf_CT3, coords_CT3, coords_resid_CT3)
#
# Saves PNG (+ optional PDF) to: figures/TPS_resid_CT3/
# ============================================================

# ---------------------------
# Load canonical project objects
# ---------------------------
source("R/00_setup_morpho.R")
source("R/01_build_metadata.R")
source("R/03_subset_CT_timeseries.R")

suppressPackageStartupMessages({
  library(geomorph)
  library(dplyr)
})

# ---------------------------
# User settings
# ---------------------------
TPS_DIR   <- file.path("figures", "TPS_resid_CT3")
DIR_OK    <- dir.exists(TPS_DIR) || dir.create(TPS_DIR, recursive = TRUE)
stopifnot(DIR_OK)

PC_K      <- 2      # +/- k SD along each PC
TPS_W     <- 4
TPS_H     <- 4
TPS_RES   <- 600
TPS_MAG   <- 2
TPS_PTCEX <- 0.7
TPS_LWD   <- 0.6

# NEW: 80% opacity drawing color (used for TPS grid in most geomorph versions)
TPS_COL   <- grDevices::adjustcolor("black", alpha.f = 0.8)

SAVE_PDF  <- TRUE
DO_EXTREME_INDIVIDUALS <- TRUE
DO_GROUP_MEDIANS       <- TRUE

# Which RAW reference to display on?
#   "CT_1950" = mean RAW shape of CT_1950 only  (recommended for time-series visualization)
#   "CT3_ALL" = mean RAW shape of all CT3 specimens
RAW_REF_MODE <- "CT_1950"

# Optional: bump magnification for medians
MEDIAN_MAG <- 2

# ============================================================
# Preconditions / sanity
# ============================================================
if (!exists("coords_CT3"))       stop("coords_CT3 not found. Ensure R/03_subset_CT_timeseries.R builds it.")
if (!exists("coords_resid_CT3")) stop("coords_resid_CT3 not found. Ensure R/03_subset_CT_timeseries.R builds it.")
if (!exists("gdf_CT3"))          stop("gdf_CT3 not found. Ensure R/03_subset_CT_timeseries.R builds it.")
if (!exists("links_all"))        stop("links_all not found. It should be created in R/00_setup_morpho.R")

stopifnot(identical(dimnames(coords_CT3)[[3]], gdf_CT3$specimen))
stopifnot(identical(dimnames(coords_resid_CT3), dimnames(coords_CT3)))
stopifnot(is.factor(gdf_CT3$group) || !is.null(gdf_CT3$group))

# ============================================================
# Reference shapes: choose RAW ref + matched RESID ref
#   IMPORTANT: ref_resid must correspond to the SAME specimen set as ref_raw
# ============================================================

idx_ref <- switch(
  RAW_REF_MODE,
  "CT_1950" = which(gdf_CT3$group == "CT_1950"),
  "CT3_ALL" = seq_len(nrow(gdf_CT3)),
  stop("RAW_REF_MODE must be one of: 'CT_1950' or 'CT3_ALL'")
)

stopifnot(length(idx_ref) >= 3)

ref_raw_CT3 <- mshape(coords_CT3[, , idx_ref, drop = FALSE])
ref_raw_CT3 <- matrix(
  ref_raw_CT3,
  ncol = 2,
  dimnames = list(dimnames(coords_CT3)[[1]], c("x", "y"))
)

ref_resid_CT3 <- mshape(coords_resid_CT3[, , idx_ref, drop = FALSE])
ref_resid_CT3 <- matrix(
  ref_resid_CT3,
  ncol = 2,
  dimnames = list(dimnames(coords_resid_CT3)[[1]], c("x", "y"))
)

stopifnot(identical(rownames(ref_raw_CT3), rownames(ref_resid_CT3)))

# Use global links (indices align because landmark order is identical)
links_CT3 <- links_all

# ============================================================
# Helpers (mirrors 1950 TPS script)
# ============================================================

save_tps_png <- function(ref, target, links, filename, title,
                         mag = TPS_MAG, width = TPS_W, height = TPS_H, res = TPS_RES,
                         pt_cex = TPS_PTCEX, lwd = TPS_LWD,
                         tps_col = TPS_COL) {
  stopifnot(is.matrix(ref), is.matrix(target))
  stopifnot(all(dim(ref) == dim(target)), ncol(ref) == 2)
  
  png(filename, width = width, height = height, units = "in", res = res)
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  
  # NEW: set default drawing colors so TPS grid inherits opacity
  par(
    mar = c(2, 2, 3, 1),
    cex = pt_cex,
    lwd = lwd,
    col = tps_col,  # default line/text color
    fg  = tps_col   # foreground color used by many base drawing functions
  )
  
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
                         pt_cex = TPS_PTCEX, lwd = TPS_LWD,
                         tps_col = TPS_COL) {
  stopifnot(is.matrix(ref), is.matrix(target))
  stopifnot(all(dim(ref) == dim(target)), ncol(ref) == 2)
  
  grDevices::cairo_pdf(filename, width = width, height = height)
  op <- par(no.readonly = TRUE)
  on.exit({ par(op); dev.off() }, add = TRUE)
  
  # NEW: set default drawing colors so TPS grid inherits opacity
  par(
    mar = c(2, 2, 3, 1),
    cex = pt_cex,
    lwd = lwd,
    col = tps_col,
    fg  = tps_col
  )
  
  geomorph::plotRefToTarget(
    ref, target,
    method = "TPS",
    mag    = mag,
    links  = links,
    main   = title
  )
}

# Synthetic shape along a PC axis at +/- k SD in whatever space coords_arr lives in
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

# Apply residual deformation onto a RAW reference:
# raw_target = ref_raw + (target_resid - ref_resid)
resid_target_on_raw <- function(target_resid, ref_resid, ref_raw) {
  stopifnot(all(dim(target_resid) == dim(ref_resid)))
  stopifnot(all(dim(ref_raw) == dim(ref_resid)))
  ref_raw + (target_resid - ref_resid)
}

# Group medians in a coords array (median at each landmark coordinate)
group_median_shapes <- function(coords_arr, groups) {
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

# ============================================================
# RESIDUAL PCA (CT3)
# ============================================================

pca_resid_CT3 <- gm.prcomp(coords_resid_CT3)
pct_resid_CT3 <- 100 * (pca_resid_CT3$sdev^2 / sum(pca_resid_CT3$sdev^2))

# ============================================================
# (A) Synthetic TPS: PC axis shapes at +/- k SD (resid), shown on RAW reference
# ============================================================

pc_axis_targets_resid <- list(
  PC1_minus = make_pc_axis_shape(coords_resid_CT3, pca_resid_CT3, pc_i = 1, sign = -1, k_sd = PC_K),
  PC1_plus  = make_pc_axis_shape(coords_resid_CT3, pca_resid_CT3, pc_i = 1, sign =  1, k_sd = PC_K),
  PC2_minus = make_pc_axis_shape(coords_resid_CT3, pca_resid_CT3, pc_i = 2, sign = -1, k_sd = PC_K),
  PC2_plus  = make_pc_axis_shape(coords_resid_CT3, pca_resid_CT3, pc_i = 2, sign =  1, k_sd = PC_K)
)

for (nm in names(pc_axis_targets_resid)) {
  pct  <- if (grepl("^PC1", nm)) pct_resid_CT3[1] else pct_resid_CT3[2]
  axis <- if (grepl("^PC1", nm)) "PC1" else "PC2"
  dirn <- if (grepl("minus$", nm)) "-" else "+"
  
  title <- paste0(
    "CT time-series SIZE-CORRECTED PCA axis: ", axis, " ", dirn, PC_K, " SD (", sprintf("%.1f", pct), "%)\n",
    "Displayed on RAW reference (", RAW_REF_MODE, ")"
  )
  
  target_on_raw <- resid_target_on_raw(pc_axis_targets_resid[[nm]], ref_resid_CT3, ref_raw_CT3)
  
  fn_png <- file.path(TPS_DIR, paste0("TPS_RESID_axisOnRaw_", RAW_REF_MODE, "_", nm, "_", PC_K, "SD.png"))
  save_tps_png(ref_raw_CT3, target_on_raw, links_CT3, fn_png, title)
  
  if (isTRUE(SAVE_PDF)) {
    fn_pdf <- file.path(TPS_DIR, paste0("TPS_RESID_axisOnRaw_", RAW_REF_MODE, "_", nm, "_", PC_K, "SD.pdf"))
    save_tps_pdf(ref_raw_CT3, target_on_raw, links_CT3, fn_pdf, title)
  }
}

# ============================================================
# (B) Optional: extreme INDIVIDUAL TPS on RESIDUAL PCA (shown on RAW reference)
# ============================================================

if (isTRUE(DO_EXTREME_INDIVIDUALS)) {
  pc_scores <- pca_resid_CT3$x
  
  pc_extreme_idx <- list(
    PC1_min = which.min(pc_scores[, 1]),
    PC1_max = which.max(pc_scores[, 1]),
    PC2_min = which.min(pc_scores[, 2]),
    PC2_max = which.max(pc_scores[, 2])
  )
  
  pc_extreme_shapes_resid <- lapply(pc_extreme_idx, function(i) coords_resid_CT3[, , i, drop = TRUE])
  
  for (nm in names(pc_extreme_shapes_resid)) {
    pct  <- if (grepl("^PC1", nm)) pct_resid_CT3[1] else pct_resid_CT3[2]
    axis <- if (grepl("^PC1", nm)) "PC1" else "PC2"
    ext  <- if (grepl("min$", nm)) "min" else "max"
    
    sp_id <- dimnames(coords_resid_CT3)[[3]][ pc_extreme_idx[[nm]] ]
    grp   <- as.character(gdf_CT3$group[ pc_extreme_idx[[nm]] ])
    
    title <- paste0(
      "CT time-series SIZE-CORRECTED PCA: ", axis, " ", ext, " individual (", sprintf("%.1f", pct), "%)\n",
      "Specimen: ", sp_id, " | Group: ", grp, " | Displayed on RAW reference (", RAW_REF_MODE, ")"
    )
    
    target_on_raw <- resid_target_on_raw(pc_extreme_shapes_resid[[nm]], ref_resid_CT3, ref_raw_CT3)
    
    fn_png <- file.path(TPS_DIR, paste0("TPS_RESID_extremeIndividualOnRaw_", RAW_REF_MODE, "_", nm, ".png"))
    save_tps_png(ref_raw_CT3, target_on_raw, links_CT3, fn_png, title)
    
    if (isTRUE(SAVE_PDF)) {
      fn_pdf <- file.path(TPS_DIR, paste0("TPS_RESID_extremeIndividualOnRaw_", RAW_REF_MODE, "_", nm, ".pdf"))
      save_tps_pdf(ref_raw_CT3, target_on_raw, links_CT3, fn_pdf, title)
    }
  }
}

# ============================================================
# (C) RESIDUAL group/year medians (shown on RAW reference)
# ============================================================

if (isTRUE(DO_GROUP_MEDIANS)) {
  grp_vec <- gdf_CT3$group
  grp_median_resid <- group_median_shapes(coords_resid_CT3, grp_vec)
  
  for (g in names(grp_median_resid)) {
    if (is.null(grp_median_resid[[g]])) next
    
    title <- paste0(
      "CT time-series SIZE-CORRECTED group median: ", g, "\n",
      "Displayed on RAW reference (", RAW_REF_MODE, ")"
    )
    
    target_on_raw <- resid_target_on_raw(grp_median_resid[[g]], ref_resid_CT3, ref_raw_CT3)
    
    fn_png <- file.path(TPS_DIR, paste0("TPS_RESID_groupMedianOnRaw_", RAW_REF_MODE, "_", g, ".png"))
    save_tps_png(ref_raw_CT3, target_on_raw, links_CT3, fn_png, title, mag = MEDIAN_MAG)
    
    if (isTRUE(SAVE_PDF)) {
      fn_pdf <- file.path(TPS_DIR, paste0("TPS_RESID_groupMedianOnRaw_", RAW_REF_MODE, "_", g, ".pdf"))
      save_tps_pdf(ref_raw_CT3, target_on_raw, links_CT3, fn_pdf, title, mag = MEDIAN_MAG)
    }
  }
}

message("Saved CT3 RESIDUAL TPS outputs to: ", normalizePath(TPS_DIR))