# ============================================================
# R/methods/01_specimen_sampling_study_design/03_subset_CT_timeseries.R
# Subset: Connecticut River ONLY (1950, 1956, 1970)
#
# Creates:
#   gdf_CT3          (metadata subset; aligned to coords order)
#   coords_CT3       (RAW GPA coords subset)
#   allo_CT3         (allometry fit object from allometry_residuals)
#   coords_resid_CT3 (RESID coords subset)
#
# Requires canonical objects from:
#   R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R  (coords_gpa, helpers)
#   R/methods/01_specimen_sampling_study_design/01_build_metadata.R (gdf, allometry_residuals)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(geomorph)
})

if (!exists("coords_gpa", inherits = TRUE)) source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
if (!exists("gdf",       inherits = TRUE)) source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")

# ---- 0) Define CT groups (no month filtering) ----
gdf_CT3 <- gdf %>%
  dplyr::filter(
    !is.na(habitat),
    habitat == "Connecticut River",
    year %in% c(1950, 1956, 1970)
  ) %>%
  dplyr::mutate(
    group = dplyr::case_when(
      year == 1950 ~ "CT_1950",
      year == 1956 ~ "CT_1956",
      year == 1970 ~ "CT_1970",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = c("CT_1950", "CT_1956", "CT_1970"))
  ) %>%
  dplyr::filter(!is.na(group)) %>%
  droplevels()

stopifnot(nrow(gdf_CT3) >= 3)
cat("\nCT time-series group counts:\n")
print(table(gdf_CT3$group))

# ---- 1) Subset coords to match metadata order ----
if (!exists("subset_coords_to_gdf", inherits = TRUE)) {
  stop("subset_coords_to_gdf() not found. It should be defined in R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
}

coords_CT3 <- subset_coords_to_gdf(coords_gpa, gdf_CT3)
stopifnot(identical(dimnames(coords_CT3)[[3]], gdf_CT3$specimen))

# ---- 2) Allometry residuals on the CT time-series subset ----
if (!exists("allometry_residuals", inherits = TRUE)) {
  stop("allometry_residuals() not found. It should be defined in R/methods/01_specimen_sampling_study_design/01_build_metadata.R")
}

size_vec <- setNames(gdf_CT3$size_for_allometry, gdf_CT3$specimen)
allo_CT3 <- allometry_residuals(coords_CT3, size_vec)

coords_resid_CT3 <- allo_CT3$residuals

cat("\nAllometry model (CT time-series): shape ~ ", unique(gdf_CT3$size_label), "\n", sep = "")
cat("\nObjects created by R/methods/01_specimen_sampling_study_design/03_subset_CT_timeseries.R:\n")
print(c("gdf_CT3", "coords_CT3", "allo_CT3", "coords_resid_CT3"))
cat("\nCT time-series subset build complete.\n")