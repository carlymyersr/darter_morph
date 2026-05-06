# ============================================================
# CT reference distribution and rarefaction checks
#
# Implements the updated outline requirement:
#   - pooled CT time points define the reference distribution
#   - PC1-PC4 scores define shared morphospace
#   - Mahalanobis distances relative to pooled CT centroid
#   - empirical 95th percentile threshold
#   - repeated CT subsampling to matched smaller n
# ============================================================

source("R/00_setup_morpho.R")
source("R/01_build_metadata.R")
source("R/04_subset_CT_timeseries_plus_1950habitats.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(geomorph)
})

OUT_DIR <- file.path("Outputs", "ct_reference_distribution")
FIG_DIR <- file.path("figures", "ct_reference_distribution")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

N_ITER <- as.integer(Sys.getenv("DARTER_RAREFACTION_ITER", "999"))
SEED <- as.integer(Sys.getenv("DARTER_RAREFACTION_SEED", "20260506"))
set.seed(SEED)

if (!exists("coords_resid_Fig6")) stop("coords_resid_Fig6 not found.")
if (!exists("gdf_Fig6")) stop("gdf_Fig6 not found.")

meta <- gdf_Fig6 %>%
  mutate(
    reference_group = case_when(
      habitat == "Connecticut River" ~ "Mainstem",
      habitat %in% c("Quabbin", "Swift River") ~ "Reservoir System",
      habitat %in% c("Fort River", "Sawmill River") ~ "Tributaries",
      TRUE ~ NA_character_
    ),
    reference_group = factor(
      reference_group,
      levels = c("Mainstem", "Reservoir System", "Tributaries")
    )
  ) %>%
  filter(!is.na(reference_group)) %>%
  droplevels()

coords <- subset_coords_to_gdf(coords_resid_Fig6, meta)
stopifnot(identical(dimnames(coords)[[3]], meta$specimen))

pca <- gm.prcomp(coords)
pc_scores <- as.data.frame(pca$x[, 1:4, drop = FALSE])
names(pc_scores) <- paste0("PC", 1:4)
pc_scores <- cbind(meta, pc_scores)

ct_scores <- pc_scores %>% filter(reference_group == "Mainstem")
if (nrow(ct_scores) < 5) stop("Need at least 5 CT specimens for PC1-PC4 Mahalanobis reference.")

ct_mat <- as.matrix(ct_scores[, paste0("PC", 1:4)])
ct_center <- colMeans(ct_mat)
ct_cov <- cov(ct_mat)

mahal_d2 <- mahalanobis(as.matrix(pc_scores[, paste0("PC", 1:4)]), center = ct_center, cov = ct_cov)
pc_scores$mahal_d2_mainstem_pc1_pc4 <- mahal_d2
pc_scores$mahal_d_mainstem_pc1_pc4 <- sqrt(mahal_d2)

ct_threshold_d2 <- unname(stats::quantile(
  pc_scores$mahal_d2_mainstem_pc1_pc4[pc_scores$reference_group == "Mainstem"],
  probs = 0.95,
  na.rm = TRUE
))

pc_scores$outside_mainstem_95 <- pc_scores$mahal_d2_mainstem_pc1_pc4 > ct_threshold_d2

mahal_summary <- pc_scores %>%
  group_by(reference_group) %>%
  summarise(
    n = n(),
    n_outside_mainstem_95 = sum(outside_mainstem_95, na.rm = TRUE),
    prop_outside_mainstem_95 = n_outside_mainstem_95 / n,
    mean_mahal_d = mean(mahal_d_mainstem_pc1_pc4, na.rm = TRUE),
    median_mahal_d = median(mahal_d_mainstem_pc1_pc4, na.rm = TRUE),
    max_mahal_d = max(mahal_d_mainstem_pc1_pc4, na.rm = TRUE),
    .groups = "drop"
  )

metric_set <- function(mat, threshold_d2 = NULL) {
  center <- colMeans(mat)
  cov_mat <- cov(mat)
  d2 <- mahalanobis(mat, center = center, cov = cov_mat)
  out <- data.frame(
    n = nrow(mat),
    mean_distance_to_centroid = mean(sqrt(rowSums(scale(mat, center = center, scale = FALSE)^2))),
    max_distance_to_centroid = max(sqrt(rowSums(scale(mat, center = center, scale = FALSE)^2))),
    pc1_pc4_generalized_variance = det(cov_mat),
    mahal_95_d2 = unname(stats::quantile(d2, 0.95, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
  if (!is.null(threshold_d2)) {
    out$prop_outside_threshold <- mean(d2 > threshold_d2)
  }
  out
}

group_counts <- table(pc_scores$reference_group)
target_n <- min(group_counts[group_counts > 0])
if (target_n < 5) stop("Rarefaction target n is too small for PC1-PC4 covariance.")

ct_idx <- which(pc_scores$reference_group == "Mainstem")
rare <- do.call(rbind, lapply(seq_len(N_ITER), function(i) {
  idx <- sample(ct_idx, target_n, replace = FALSE)
  m <- as.matrix(pc_scores[idx, paste0("PC", 1:4)])
  cbind(iteration = i, metric_set(m, threshold_d2 = ct_threshold_d2))
}))

observed_groups <- do.call(rbind, lapply(levels(pc_scores$reference_group), function(g) {
  sub <- pc_scores %>% filter(reference_group == g)
  m <- as.matrix(sub[, paste0("PC", 1:4)])
  cbind(reference_group = g, metric_set(m, threshold_d2 = ct_threshold_d2))
}))

rare_summary <- rare %>%
  summarise(
    target_n = target_n,
    iterations = n(),
    mean_generalized_variance = mean(pc1_pc4_generalized_variance, na.rm = TRUE),
    lo_generalized_variance = quantile(pc1_pc4_generalized_variance, 0.025, na.rm = TRUE),
    hi_generalized_variance = quantile(pc1_pc4_generalized_variance, 0.975, na.rm = TRUE),
    mean_mahal_95_d2 = mean(mahal_95_d2, na.rm = TRUE),
    lo_mahal_95_d2 = quantile(mahal_95_d2, 0.025, na.rm = TRUE),
    hi_mahal_95_d2 = quantile(mahal_95_d2, 0.975, na.rm = TRUE),
    mean_max_distance = mean(max_distance_to_centroid, na.rm = TRUE),
    lo_max_distance = quantile(max_distance_to_centroid, 0.025, na.rm = TRUE),
    hi_max_distance = quantile(max_distance_to_centroid, 0.975, na.rm = TRUE)
  )

write.csv(pc_scores, file.path(OUT_DIR, "ct_reference_pc_scores_mahalanobis.csv"), row.names = FALSE)
write.csv(mahal_summary, file.path(OUT_DIR, "ct_reference_mahalanobis_summary_by_group.csv"), row.names = FALSE)
write.csv(rare, file.path(OUT_DIR, "ct_reference_rarefaction_iterations.csv"), row.names = FALSE)
write.csv(rare_summary, file.path(OUT_DIR, "ct_reference_rarefaction_summary.csv"), row.names = FALSE)
write.csv(observed_groups, file.path(OUT_DIR, "ct_reference_observed_group_occupancy_metrics.csv"), row.names = FALSE)

p_rare <- ggplot(rare, aes(x = pc1_pc4_generalized_variance)) +
  geom_histogram(bins = 40, fill = "grey70", color = "white") +
  geom_vline(
    data = observed_groups,
    aes(xintercept = pc1_pc4_generalized_variance, color = reference_group),
    linewidth = 0.8
  ) +
  labs(
    x = "PC1-PC4 generalized variance",
    y = "Rarefied CT subsamples",
    color = "Observed group"
  ) +
  theme_classic(base_size = 9)

ggsave(file.path(FIG_DIR, "ct_reference_rarefaction_generalized_variance.pdf"), p_rare, width = 6, height = 4)
ggsave(file.path(FIG_DIR, "ct_reference_rarefaction_generalized_variance.png"), p_rare, width = 6, height = 4, dpi = 400)

capture.output(
  {
    cat("CT reference distribution: PC1-PC4 Mahalanobis + rarefaction\n\n")
    cat("Iterations:", N_ITER, "\n")
    cat("Rarefaction target n:", target_n, "\n")
    cat("Mainstem empirical 95th percentile Mahalanobis D2:", ct_threshold_d2, "\n\n")
    cat("Mahalanobis summary by group:\n")
    print(mahal_summary)
    cat("\nRarefaction summary:\n")
    print(rare_summary)
    cat("\nObserved occupancy metrics:\n")
    print(observed_groups)
  },
  file = file.path(OUT_DIR, "ct_reference_mahalanobis_rarefaction_summary.txt")
)

cat("\nCT reference Mahalanobis + rarefaction outputs saved to:\n")
cat("  ", normalizePath(OUT_DIR), "\n")
