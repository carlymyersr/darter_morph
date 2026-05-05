# ============================================================
# R/04_subset_CT_timeseries_plus_1950habitats.R
# Subset: CT time series (1950/1956/1970) + ALL 1950 habitats
# Plus: allometry residuals on logCsize for this combined subset
#
# Creates:
#   gdf_Fig6            (metadata subset with group)
#   coords_Fig6         (RAW GPA coords subset aligned to gdf_Fig6)
#   allo_Fig6           (list from allometry_residuals: residuals + fit)
#   coords_resid_Fig6   (RESID coords subset)
#
# Requires canonical objects:
#   coords_gpa, gdf, subset_coords_to_gdf()
#   allometry_residuals() from R/01_build_metadata.R
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(geomorph)
})

if (!exists("coords_gpa", inherits = TRUE)) source("R/00_setup_morpho.R")
if (!exists("gdf",       inherits = TRUE)) source("R/01_build_metadata.R")

# ---- 0) Define combined inclusion rule ----
gdf_Fig6 <- gdf %>%
  dplyr::filter(
    !is.na(habitat),
    (
      (habitat == "Connecticut River" & year %in% c(1950, 1956, 1970)) |
        (year == 1950)
    )
  ) %>%
  dplyr::mutate(
    group = dplyr::case_when(
      habitat == "Connecticut River" & year == 1950 ~ "CT_1950",
      habitat == "Connecticut River" & year == 1956 ~ "CT_1956",
      habitat == "Connecticut River" & year == 1970 ~ "CT_1970",
      year == 1950 ~ paste0(habitat, "_1950"),
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(group)) %>%
  dplyr::mutate(
    group = factor(
      group,
      levels = c(
        "CT_1950", "CT_1956", "CT_1970",
        "Fort River_1950", "Quabbin_1950", "Sawmill River_1950", "Swift River_1950"
      )
    )
  ) %>%
  droplevels()

stopifnot(nrow(gdf_Fig6) >= 3)

cat("\nFigure 6/7 group counts:\n")
print(table(gdf_Fig6$group, useNA = "ifany"))

# ---- 1) Subset coords to match metadata order ----
if (!exists("subset_coords_to_gdf", inherits = TRUE)) {
  stop("subset_coords_to_gdf() not found. It should be defined in R/00_setup_morpho.R")
}
coords_Fig6 <- subset_coords_to_gdf(coords_gpa, gdf_Fig6)
stopifnot(identical(dimnames(coords_Fig6)[[3]], gdf_Fig6$specimen))

# ---- 2) Allometry residuals on logCsize within this combined subset ----
if (!exists("allometry_residuals", inherits = TRUE)) {
  stop("allometry_residuals() not found. It should be defined in R/01_build_metadata.R")
}

size_vec <- setNames(gdf_Fig6$size_for_allometry, gdf_Fig6$specimen)
allo_Fig6 <- allometry_residuals(coords_Fig6, size_vec)

coords_resid_Fig6 <- allo_Fig6$residuals

cat("\nAllometry model (Fig6/7 combined): shape ~ ", unique(gdf_Fig6$size_label), "\n", sep = "")
cat("\nObjects created by R/04_subset_CT_timeseries_plus_1950habitats.R:\n")
print(c("gdf_Fig6", "coords_Fig6", "allo_Fig6", "coords_resid_Fig6"))
cat("\nFigure 6/7 subset build complete.\n")