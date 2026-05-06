# ============================================================
# R/methods/01_specimen_sampling_study_design/05_subset_1950_quabbin_swift.R
# Canonical Quabbin + Swift River 1950 subset.
#
# Requires:
#   R/methods/01_specimen_sampling_study_design/02_subset_1950.R
#
# Creates:
#   gdf_1950_qs
#   coords_1950_qs
#   coords_resid_1950_qs
# ============================================================

if (!exists("gdf_1950", inherits = TRUE) ||
    !exists("coords_1950", inherits = TRUE) ||
    !exists("coords_resid_1950", inherits = TRUE)) {
  source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
  source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")
  source("R/methods/01_specimen_sampling_study_design/02_subset_1950.R")
}

qs_idx <- gdf_1950$habitat %in% c("Quabbin", "Swift River")

gdf_1950_qs <- gdf_1950[qs_idx, , drop = FALSE]
coords_1950_qs <- coords_1950[, , qs_idx, drop = FALSE]
coords_resid_1950_qs <- coords_resid_1950[, , qs_idx, drop = FALSE]

stopifnot(identical(dimnames(coords_1950_qs)[[3]], gdf_1950_qs$specimen))
stopifnot(identical(dimnames(coords_resid_1950_qs)[[3]], gdf_1950_qs$specimen))

QS1950_OBJECTS <- c("gdf_1950_qs", "coords_1950_qs", "coords_resid_1950_qs")

if (!exists("VERBOSE") || isTRUE(VERBOSE)) {
  cat("\nObjects created by R/methods/01_specimen_sampling_study_design/05_subset_1950_quabbin_swift.R:\n")
  print(QS1950_OBJECTS)
}
