# ============================================================
# Scripts/1950_procD_variation_source.R
# Procrustes ANOVA on size-corrected 1950 shapes
# Variation source grouping:
#   - Connecticut River
#   - Quabbin
#   - Other rivers
#
# Requires canonical project scripts:
#   R/00_setup_morpho.R
#   R/01_build_metadata.R
#   R/02_subset_1950.R
# ============================================================

# ---------------------------
# Load canonical project objects
# ---------------------------
source("R/00_setup_morpho.R")
source("R/01_build_metadata.R")
source("R/02_subset_1950.R")

suppressPackageStartupMessages({
  library(geomorph)
  library(dplyr)
  library(RRPP)
})

# ---------------------------
# Output directory
# ---------------------------
OUT_DIR <- file.path("Outputs", "1950_variationSource_procD")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ---------------------------
# Preconditions
# ---------------------------
if (!exists("coords_resid_1950")) stop("coords_resid_1950 not found. Did R/02_subset_1950.R run fully?")
if (!exists("gdf_1950"))         stop("gdf_1950 not found. Did R/02_subset_1950.R run fully?")
stopifnot(identical(dimnames(coords_resid_1950)[[3]], gdf_1950$specimen))

# Nice label for output text
size_label <- if ("size_label" %in% names(gdf_1950) && !all(is.na(gdf_1950$size_label))) {
  unique(na.omit(gdf_1950$size_label))[1]
} else if ("size_label" %in% names(gdf) && !all(is.na(gdf$size_label))) {
  unique(na.omit(gdf$size_label))[1]
} else {
  "size_for_allometry"
}

# ---------------------------
# Define variation source grouping
# ---------------------------
gdf_1950 <- gdf_1950 %>%
  mutate(
    variation_source = case_when(
      habitat == "Connecticut River" ~ "Connecticut River",
      habitat == "Quabbin"           ~ "Quabbin",
      habitat %in% c("Sawmill River", "Swift River", "Fort River") ~ "Other rivers",
      TRUE ~ NA_character_
    ),
    variation_source = factor(
      variation_source,
      levels = c("Connecticut River", "Quabbin", "Other rivers")
    )
  )

keep <- !is.na(gdf_1950$variation_source)
gdf_vs <- gdf_1950[keep, , drop = FALSE]
coords_vs <- coords_resid_1950[, , keep, drop = FALSE]

stopifnot(identical(dimnames(coords_vs)[[3]], gdf_vs$specimen))

cat("Variation-source counts:\n")
print(table(gdf_vs$variation_source))
cat("\nHabitat x variation-source table:\n")
print(table(gdf_vs$habitat, gdf_vs$variation_source))

# ---------------------------
# ProcD ANOVA on size-corrected shapes
# Response is already size-corrected residual shape
# ---------------------------
fit_vs <- geomorph::procD.lm(
  coords_vs ~ variation_source,
  data = gdf_vs,
  RRPP = TRUE,
  iter = 999
)

fit_sum <- summary(fit_vs)

# ---------------------------
# Pairwise tests among group means
# ---------------------------
pair_vs <- RRPP::pairwise(fit_vs, groups = gdf_vs$variation_source)
pair_sum <- summary(pair_vs)

# ---------------------------
# Save results
# ---------------------------
saveRDS(fit_vs,  file.path(OUT_DIR, "fit_procD_variationSource_1950.rds"))
saveRDS(pair_vs, file.path(OUT_DIR, "pairwise_procD_variationSource_1950.rds"))
saveRDS(gdf_vs,  file.path(OUT_DIR, "gdf_variationSource_1950.rds"))

sink(file.path(OUT_DIR, "summary_procD_variationSource_1950.txt"))
cat("============================================================\n")
cat("1950 variation-source Procrustes ANOVA (size-corrected)\n")
cat("============================================================\n\n")
cat("Allometry correction applied upstream in R/02_subset_1950.R\n")
cat("Residual model: shape ~ ", size_label, "\n", sep = "")
cat("ProcD model tested here: residual shape ~ variation_source\n\n")
cat("Group counts:\n")
print(table(gdf_vs$variation_source))
cat("\nHabitat x variation-source table:\n")
print(table(gdf_vs$habitat, gdf_vs$variation_source))
cat("\n--- Procrustes ANOVA summary ---\n")
print(fit_sum)
cat("\n--- Pairwise comparisons ---\n")
print(pair_sum)
cat("\nOutput directory:\n", normalizePath(OUT_DIR), "\n", sep = "")
sink()

# Optional CSV of ANOVA table if possible
anova_tab <- tryCatch(as.data.frame(fit_sum$table), error = function(e) NULL)
if (!is.null(anova_tab)) {
  write.csv(anova_tab,
            file = file.path(OUT_DIR, "anova_table_variationSource_1950.csv"),
            row.names = TRUE)
}

# ---------------------------
# Console report
# ---------------------------
cat("\nSaved outputs to:\n  ", normalizePath(OUT_DIR), "\n", sep = "")
cat("Main summary file:\n  ", normalizePath(file.path(OUT_DIR, "summary_procD_variationSource_1950.txt")), "\n", sep = "")
