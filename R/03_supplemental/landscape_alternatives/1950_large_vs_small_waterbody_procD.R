# ============================================================
# Binary 1950 contrast:
# Connecticut River + Quabbin  vs  all other 1950 habitats
# Size-corrected shape (residuals)
# ============================================================

source("R/00_setup_morpho.R")
source("R/01_build_metadata.R")
source("R/02_subset_1950.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(geomorph)
  library(RRPP)
})

# ---------------------------
# Sanity checks
# ---------------------------
stopifnot(exists("gdf_1950"), exists("coords_resid_1950"))
stopifnot(identical(dimnames(coords_resid_1950)[[3]], gdf_1950$specimen))

# Confirm habitat labels present
print(table(gdf_1950$habitat, useNA = "ifany"))

# ---------------------------
# Define binary grouping
# ---------------------------
gdf_1950 <- gdf_1950 %>%
  mutate(
    habitat_binary = case_when(
      habitat %in% c("Connecticut River", "Quabbin") ~ "CT_Quabbin",
      habitat %in% c("Sawmill River", "Swift River", "Fort River") ~ "Other_1950_habitats",
      TRUE ~ NA_character_
    ),
    habitat_binary = factor(
      habitat_binary,
      levels = c("CT_Quabbin", "Other_1950_habitats")
    )
  )

# Drop any unexpected NA rows just in case
keep <- !is.na(gdf_1950$habitat_binary)
gdf_bin <- gdf_1950[keep, , drop = FALSE]
coords_bin <- coords_resid_1950[, , keep, drop = FALSE]

stopifnot(identical(dimnames(coords_bin)[[3]], gdf_bin$specimen))

# Inspect grouping
print(table(gdf_bin$habitat, gdf_bin$habitat_binary))
print(table(gdf_bin$habitat_binary))

# ---------------------------
# Run Procrustes ANOVA
# ---------------------------
f_bin <- coords_bin ~ habitat_binary
fit_bin <- geomorph::procD.lm(
  f1   = f_bin,
  data = gdf_bin,
  iter = 999,
  RRPP = TRUE
)

# Print results
print(fit_bin)
anova_bin <- anova(fit_bin)
print(anova_bin)

# Optional pairwise object (not necessary for 2 groups, but fine)
pw_bin <- RRPP::pairwise(fit_bin, groups = gdf_bin$habitat_binary)
print(pw_bin)

# ---------------------------
# Save outputs
# ---------------------------
out_dir <- file.path("Outputs", "binary_CT_Quabbin_vs_other_1950")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(fit_bin, file.path(out_dir, "fit_bin_resid_shape_by_binary_habitat.rds"))
saveRDS(anova_bin, file.path(out_dir, "anova_bin_resid_shape_by_binary_habitat.rds"))
saveRDS(pw_bin, file.path(out_dir, "pairwise_bin_resid_shape_by_binary_habitat.rds"))

writeLines(
  c(
    "===== fit print =====",
    capture.output(print(fit_bin)),
    "",
    "===== anova(fit_bin) =====",
    capture.output(print(anova_bin)),
    "",
    "===== habitat by binary group =====",
    capture.output(print(table(gdf_bin$habitat, gdf_bin$habitat_binary))),
    "",
    "===== binary group counts =====",
    capture.output(print(table(gdf_bin$habitat_binary)))
  ),
  file.path(out_dir, "binary_CT_Quabbin_vs_other_1950_summary.txt")
)

message("Done. Outputs written to: ", normalizePath(out_dir))