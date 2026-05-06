# ============================================================
# QS_context_dependence_analysis.R
# Tests whether Quabbin vs Swift structure depends on PCA context
# ============================================================

suppressPackageStartupMessages({
  library(geomorph)
  library(dplyr)
  library(ggplot2)
})

# ============================================================
# LOAD PIPELINE OBJECTS
# ============================================================

source("R/00_setup_morpho.R")
source("R/01_build_metadata.R")
source("R/02_subset_1950.R")
source("R/04_subset_CT_timeseries_plus_1950habitats.R")
source("R/05_subset_1950_quabbin_swift.R")  # :contentReference[oaicite:0]{index=0}

# ============================================================
# 1) QS ONLY (1950)
# ============================================================

cat("\n================ QS ONLY (1950) ================\n")

pca_qs <- gm.prcomp(coords_resid_1950_qs)

fit_qs <- procD.lm(coords_resid_1950_qs ~ habitat,
                   data = gdf_1950_qs,
                   iter = 999)

print(summary(fit_qs))   # <-- IMPORTANT

# ============================================================
# 2) QS within 1950 landscape
# ============================================================

cat("\n================ QS IN 1950 LANDSCAPE ================\n")

pca_1950 <- gm.prcomp(coords_resid_1950)

qs_idx_1950 <- gdf_1950$habitat %in% c("Quabbin", "Swift River")

coords_qs_1950 <- coords_resid_1950[, , qs_idx_1950]
gdf_qs_1950    <- gdf_1950[qs_idx_1950, ]

fit_qs_1950 <- procD.lm(coords_qs_1950 ~ habitat,
                        data = gdf_qs_1950,
                        iter = 999)

print(summary(fit_qs_1950))  # <-- IMPORTANT

# ============================================================
# 3) QS within FULL dataset
# ============================================================

cat("\n================ QS IN FULL DATASET ================\n")

pca_full <- gm.prcomp(coords_resid_Fig6)

qs_idx_full <- gdf_Fig6$habitat %in% c("Quabbin", "Swift River")

coords_qs_full <- coords_resid_Fig6[, , qs_idx_full]
gdf_qs_full    <- gdf_Fig6[qs_idx_full, ]

fit_qs_full <- procD.lm(coords_qs_full ~ habitat,
                        data = gdf_qs_full,
                        iter = 999)

print(summary(fit_qs_full))  # <-- IMPORTANT

# ============================================================
# Pairwise distances (THIS is often what you actually want)
# ============================================================

cat("\n================ PAIRWISE DISTANCES (QS ONLY) ================\n")

pw_qs <- pairwise(fit_qs, groups = gdf_1950_qs$habitat)
print(summary(pw_qs))

# ============================================================
# OPTIONAL: Compare mean shapes directly (robust metric)
# ============================================================

cat("\n================ PAIRWISE DISTANCES ================\n")

pairwise_qs <- pairwise(fit_qs, groups = gdf_1950_qs$habitat)
summary(pairwise_qs)

# ============================================================
# SAVE PCA SCORES FOR PLOTTING
# ============================================================

write.csv(pca_qs$x, "Outputs/PCA_QS_only_scores.csv")
write.csv(pca_1950$x, "Outputs/PCA_1950_scores.csv")
write.csv(pca_full$x, "Outputs/PCA_full_scores.csv")

cat("\nAnalysis complete.\n")