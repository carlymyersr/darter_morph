# ============================================================
# Combined allometry models for CT time series + 1950 habitats
#
# Methods outline requirement:
#   - shape ~ logCsize
#   - shape ~ logCsize * group
#   - combined CT time series and 1950 waterbodies/habitats
# ============================================================

source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")
source("R/methods/01_specimen_sampling_study_design/04_subset_CT_timeseries_plus_1950habitats.R")

suppressPackageStartupMessages({
  library(geomorph)
})

OUT_DIR <- file.path("Outputs", "combined_allometry_procD")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!exists("coords_Fig6")) stop("coords_Fig6 not found.")
if (!exists("gdf_Fig6")) stop("gdf_Fig6 not found.")

meta <- gdf_Fig6
coords <- coords_Fig6

if (!"group" %in% names(meta)) stop("gdf_Fig6$group not found.")
if (!"logCsize" %in% names(meta)) stop("gdf_Fig6$logCsize not found.")

meta$group <- droplevels(factor(meta$group))
stopifnot(identical(dimnames(coords)[[3]], meta$specimen))

fit_size <- procD.lm(
  coords ~ logCsize,
  data = meta,
  iter = 999,
  RRPP = TRUE
)

fit_size_group <- procD.lm(
  coords ~ logCsize * group,
  data = meta,
  iter = 999,
  RRPP = TRUE
)

capture.output(
  {
    cat("Combined allometry model: shape ~ logCsize\n\n")
    print(summary(fit_size))
  },
  file = file.path(OUT_DIR, "combined_shape_by_logCsize_summary.txt")
)

capture.output(
  {
    cat("Combined allometry interaction model: shape ~ logCsize * group\n\n")
    print(summary(fit_size_group))
  },
  file = file.path(OUT_DIR, "combined_shape_by_logCsize_x_group_summary.txt")
)

write.csv(
  data.frame(
    group = names(table(meta$group)),
    n = as.integer(table(meta$group)),
    stringsAsFactors = FALSE
  ),
  file.path(OUT_DIR, "combined_allometry_group_counts.csv"),
  row.names = FALSE
)

saveRDS(fit_size, file.path(OUT_DIR, "combined_shape_by_logCsize.rds"))
saveRDS(fit_size_group, file.path(OUT_DIR, "combined_shape_by_logCsize_x_group.rds"))

cat("\nCombined allometry procD outputs saved to:\n")
cat("  ", normalizePath(OUT_DIR), "\n")
