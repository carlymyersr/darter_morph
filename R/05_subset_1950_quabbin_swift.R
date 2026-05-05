# ============================================================
# R/05_subset_1950_quabbin_swift.R
# Canonical Quabbin + Swift River 1950 subset.
#
# Requires:
#   R/02_subset_1950.R
#
# Creates:
#   gdf_1950_qs
#   coords_1950_qs
#   coords_resid_1950_qs
# ============================================================

if (!exists("gdf_1950", inherits = TRUE) ||
    !exists("coords_1950", inherits = TRUE) ||
    !exists("coords_resid_1950", inherits = TRUE)) {
  source("R/00_setup_morpho.R")
  source("R/01_build_metadata.R")
  source("R/02_subset_1950.R")
}

qs_idx <- gdf_1950$habitat %in% c("Quabbin", "Swift River")

gdf_1950_qs <- gdf_1950[qs_idx, , drop = FALSE]
coords_1950_qs <- coords_1950[, , qs_idx, drop = FALSE]
coords_resid_1950_qs <- coords_resid_1950[, , qs_idx, drop = FALSE]

stopifnot(identical(dimnames(coords_1950_qs)[[3]], gdf_1950_qs$specimen))
stopifnot(identical(dimnames(coords_resid_1950_qs)[[3]], gdf_1950_qs$specimen))

QS1950_OBJECTS <- c("gdf_1950_qs", "coords_1950_qs", "coords_resid_1950_qs")

if (!exists("VERBOSE") || isTRUE(VERBOSE)) {
  cat("\nObjects created by R/05_subset_1950_quabbin_swift.R:\n")
  print(QS1950_OBJECTS)
}
