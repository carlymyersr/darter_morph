# ============================================================
# R/02_subset_1950.R
# Canonical 1950 subset objects (geometry + references + allometry)
#
# Requires:
#   from R/00_setup_morpho.R:
#     coords_gpa, curve_endpoints
#     helper fns: make_curve_index, make_links_from_index
#     (optional but recommended): trim_dimnames_landmarks, subset_coords_to_gdf
#
#   from R/01_build_metadata.R:
#     gdf, hab_palette, allometry_residuals
#
# Outputs:
#   gdf_1950, coords_1950, ref_1950, links_1950,
#   fit_allo_1950, coords_resid_1950
#   + SUBSET1950_OBJECTS (character vector of object names created)
# ============================================================

# ---------------------------
# Flags (inherits if already defined)
# ---------------------------
if (!exists("DO_SANITY_PLOTS")) DO_SANITY_PLOTS <- FALSE
if (!exists("VERBOSE"))         VERBOSE         <- TRUE

if (VERBOSE) cat("Running R/02_subset_1950.R...\n")

# ---------------------------
# Preconditions
# ---------------------------
if (!exists("gdf"))         stop("gdf not found. Run source('R/01_build_metadata.R') first.")
if (!exists("coords_gpa"))  stop("coords_gpa not found. Run source('R/00_setup_morpho.R') first.")
if (!exists("hab_palette")) stop("hab_palette not found. Run source('R/01_build_metadata.R') first.")
if (!exists("curve_endpoints")) stop("curve_endpoints not found. Run source('R/00_setup_morpho.R') first.")
if (!exists("allometry_residuals")) stop("allometry_residuals() not found. Run source('R/01_build_metadata.R') first.")

suppressPackageStartupMessages({
  library(dplyr)
  library(geomorph)
})

# ---------------------------
# Minimal fallbacks (only if missing)
# ---------------------------
if (!exists("subset_coords_to_gdf")) {
  subset_coords_to_gdf <- function(coords_array, gdf_meta) {
    ids_all <- dimnames(coords_array)[[3]]
    idx <- match(gdf_meta$specimen, ids_all)
    stopifnot(!any(is.na(idx)))
    
    out <- coords_array[, , idx, drop = FALSE]
    dimnames(out)[[3]] <- gdf_meta$specimen
    stopifnot(identical(dimnames(out)[[3]], gdf_meta$specimen))
    out
  }
}

if (!exists("trim_dimnames_landmarks")) {
  trim_names <- function(x) {
    x <- as.character(x)
    x <- trimws(x)
    gsub("\\s+", " ", x)
  }
  trim_dimnames_landmarks <- function(coords_arr) {
    dn <- dimnames(coords_arr)
    dn[[1]] <- trim_names(dn[[1]])
    dimnames(coords_arr) <- dn
    coords_arr
  }
}

if (!exists("make_curve_index")) {
  stop("make_curve_index not found. Keep it in R/00_setup_morpho.R (Helper functions section).")
}
if (!exists("make_links_from_index")) {
  stop("make_links_from_index not found. Keep it in R/00_setup_morpho.R (Helper functions section).")
}

# ============================================================
# 12) Canonical subset: 1950 only, valid habitats
# ============================================================

gdf_1950 <- gdf %>%
  filter(year == 1950, !is.na(habitat)) %>%
  droplevels()

stopifnot(nrow(gdf_1950) >= 3)

# Subset coords to 1950 specimens (aligned to gdf_1950 order)
coords_1950 <- subset_coords_to_gdf(coords_gpa, gdf_1950)

# Clean up landmark names (trailing spaces etc.)
coords_1950 <- trim_dimnames_landmarks(coords_1950)

# Enforce alignment
stopifnot(identical(dimnames(coords_1950)[[3]], gdf_1950$specimen))

# ============================================================
# 12b) Allometry correction (centroid size default)
#   Model: shape ~ size_for_allometry (typically logCsize)
# Outputs:
#   fit_allo_1950, coords_resid_1950
# ============================================================

if (!("size_for_allometry" %in% names(gdf_1950))) {
  stop("gdf_1950 is missing size_for_allometry. Check R/01_build_metadata.R.")
}

# Use names() so helper can enforce specimen alignment robustly
size_vec_1950 <- gdf_1950$size_for_allometry
names(size_vec_1950) <- gdf_1950$specimen

allo_1950 <- allometry_residuals(coords_1950, size_vec_1950)

coords_resid_1950 <- allo_1950$residuals
fit_allo_1950     <- allo_1950$fit

stopifnot(identical(dimnames(coords_resid_1950), dimnames(coords_1950)))

if (VERBOSE) {
  # Optional label if you stored it in gdf (as in your 01_build_metadata.R)
  if ("size_label" %in% names(gdf) && length(unique(gdf$size_label)) >= 1) {
    cat("\nAllometry model (1950): shape ~ ", unique(gdf$size_label)[1], "\n", sep = "")
  } else {
    cat("\nAllometry model (1950): shape ~ size_for_allometry\n")
  }
}

# ============================================================
# Reference (mean) shape for 1950 (RAW)
# ============================================================

ref_1950 <- mshape(coords_1950)
ref_1950 <- matrix(
  ref_1950,
  ncol = 2,
  dimnames = list(dimnames(coords_1950)[[1]], c("x", "y"))
)

# ============================================================
# Build 1950 links for wireframe / TPS (RAW landmark set)
# ============================================================

pt_names_1950 <- dimnames(coords_1950)[[1]]

idx_hp_1950 <- make_curve_index(
  "hyoid_pelvic",
  curve_endpoints$hyoid_pelvic,
  pt_names_1950
)

idx_co_1950 <- make_curve_index(
  "cranium_orbital",
  curve_endpoints$cranium_orbital,
  pt_names_1950
)

links_1950 <- rbind(
  make_links_from_index(idx_hp_1950),
  make_links_from_index(idx_co_1950)
)

# Ensure palette covers all habitats present in 1950
missing_cols <- setdiff(levels(factor(gdf_1950$habitat)), names(hab_palette))
if (length(missing_cols) > 0) {
  stop("hab_palette missing colors for: ", paste(missing_cols, collapse = ", "))
}

# ============================================================
# Optional sanity plots
# ============================================================

if (DO_SANITY_PLOTS) {
  plot(ref_1950, asp = 1, pch = 16, main = "ref_1950 (mean RAW shape)")
  segments(
    x0 = ref_1950[links_1950[, 1], 1],
    y0 = ref_1950[links_1950[, 1], 2],
    x1 = ref_1950[links_1950[, 2], 1],
    y1 = ref_1950[links_1950[, 2], 2]
  )
}

# ============================================================
# Canonical object list (setup-style)
# ============================================================

SUBSET1950_OBJECTS <- c(
  "gdf_1950",
  "coords_1950",
  "coords_resid_1950",
  "fit_allo_1950",
  "ref_1950",
  "pt_names_1950",
  "idx_hp_1950",
  "idx_co_1950",
  "links_1950"
)

if (VERBOSE) {
  cat("\nObjects created by R/02_subset_1950.R:\n")
  print(SUBSET1950_OBJECTS)
  cat("\n1950 subset build complete.\n")
}