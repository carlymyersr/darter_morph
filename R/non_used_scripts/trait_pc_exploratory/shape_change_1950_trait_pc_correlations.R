# ============================================================
# shape_change_along_pcs_1950.R
#
# Goal:
#   Relate PC1-PC4 from the 1950 size-corrected habitat-only PCA
#   to all measured traits currently defined in:
#     - landmark_distance_measurements.R
#     - curve_shape_metrics_raw_then_sizecorrected.R
#
# 1950 habitats included:
#   - Connecticut River
#   - Quabbin
#   - Swift River
#   - Fort River
#   - Sawmill River
#
# Outputs saved to:
#   Outputs/shape_change_along_pcs_1950/
#
# Main outputs:
#   1) PC scores for every specimen (PC1-PC4)
#   2) Specimen-level trait table (all measured traits)
#   3) Trait-by-PC regression summary table
#   4) Trait-by-PC correlation summary table
#   5) Wide tables of R2 and correlation for quick reading
#   6) Heatmaps of signed correlation and R2
#   7) Scatterplots for the top trait associations for each PC
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(geomorph)
  library(tibble)
})

# ------------------------------------------------------------
# 0) Canonical project objects
# ------------------------------------------------------------
source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")
source("R/methods/01_specimen_sampling_study_design/02_subset_1950.R")

# ------------------------------------------------------------
# 1) Output directory
# ------------------------------------------------------------
out_dir <- file.path("Outputs", "shape_change_along_pcs_1950")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

plot_dir <- file.path(out_dir, "plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ------------------------------------------------------------
# 2) Define 1950 habitat groups
# ------------------------------------------------------------
habitat_levels_1950 <- c(
  "Connecticut River",
  "Quabbin",
  "Swift River",
  "Fort River",
  "Sawmill River"
)

gdf_pc <- gdf_1950 %>%
  filter(habitat %in% habitat_levels_1950) %>%
  mutate(group = factor(habitat, levels = habitat_levels_1950)) %>%
  droplevels()

stopifnot(nrow(gdf_pc) > 0)

keep_idx <- match(gdf_pc$specimen, dimnames(coords_resid_1950)[[3]])
stopifnot(!any(is.na(keep_idx)))

coords_resid_pc <- coords_resid_1950[, , keep_idx, drop = FALSE]
stopifnot(identical(dimnames(coords_resid_pc)[[3]], gdf_pc$specimen))

# ------------------------------------------------------------
# 3) PCA on 1950 residual shape coordinates
# ------------------------------------------------------------
pca_1950 <- gm.prcomp(coords_resid_pc)

scores_1950 <- as.data.frame(pca_1950$x[, 1:4, drop = FALSE])
colnames(scores_1950) <- c("PC1", "PC2", "PC3", "PC4")

scores_1950 <- scores_1950 %>%
  rownames_to_column("specimen") %>%
  left_join(
    gdf_pc %>%
      select(specimen, habitat, year, group, size_for_allometry, size_label),
    by = "specimen"
  )

# variance explained table
var_expl <- (pca_1950$d^2) / sum(pca_1950$d^2)
var_table <- tibble(
  PC = paste0("PC", seq_along(var_expl)),
  variance_explained = var_expl,
  percent_explained = 100 * var_expl
)

# ------------------------------------------------------------
# 4) Helper functions for raw-coordinate trait extraction
# ------------------------------------------------------------
subset_coords_to_gdf_raw <- function(coords_array, gdf_meta) {
  ids_all <- dimnames(coords_array)[[3]]
  idx <- match(gdf_meta$specimen, ids_all)
  stopifnot(!any(is.na(idx)))
  
  out <- coords_array[, , idx, drop = FALSE]
  dimnames(out)[[3]] <- gdf_meta$specimen
  stopifnot(identical(dimnames(out)[[3]], gdf_meta$specimen))
  out
}

get_distance_measure <- function(coords_arr, landmark_a, landmark_b) {
  pt_names <- dimnames(coords_arr)[[1]]
  
  if (!(landmark_a %in% pt_names)) stop("Landmark not found in coords_arr: ", landmark_a)
  if (!(landmark_b %in% pt_names)) stop("Landmark not found in coords_arr: ", landmark_b)
  
  A <- coords_arr[landmark_a, , , drop = FALSE]
  B <- coords_arr[landmark_b, , , drop = FALSE]
  
  ax <- as.numeric(A[1, 1, ])
  ay <- as.numeric(A[1, 2, ])
  bx <- as.numeric(B[1, 1, ])
  by <- as.numeric(B[1, 2, ])
  
  sqrt((ax - bx)^2 + (ay - by)^2)
}

extract_distance_measurements <- function(coords_arr, gdf_meta) {
  out <- gdf_meta %>%
    select(specimen, habitat, year, group, size_for_allometry, size_label)
  
  out$Eye_width <- get_distance_measure(coords_arr, "orbit_1", "orbit_2")
  out$Body_depth <- get_distance_measure(coords_arr, "premaxilla", "maxilla")
  out$Operculum_width <- get_distance_measure(coords_arr, "max_curve_preoperculum", "operculum")
  out$Jaw_muscle_length <- get_distance_measure(coords_arr, "max_curve_preoperculum", "preoperculum")
  
  out
}

size_correct_measurements <- function(df_wide, measurement_cols) {
  df_out <- df_wide
  
  for (m in measurement_cols) {
    if (any(df_out[[m]] <= 0, na.rm = TRUE)) {
      stop("Non-positive values found in measurement ", m, "; cannot log-transform.")
    }
    
    log_y <- log(df_out[[m]])
    fit <- lm(log_y ~ size_for_allometry, data = df_out)
    corrected <- exp(residuals(fit) + mean(log_y, na.rm = TRUE))
    df_out[[paste0(m, "_sizecorr")]] <- corrected
  }
  
  df_out
}

curve_defs <- list(
  snout = c(
    "cranium_orbital_start",
    "cranium_orbital_sl1", "cranium_orbital_sl2", "cranium_orbital_sl3",
    "cranium_orbital_sl4", "cranium_orbital_sl5", "cranium_orbital_sl6",
    "cranium_orbital_sl7", "cranium_orbital_sl8",
    "cranium_orbital_end"
  ),
  hyoid = c(
    "hyoid_pelvic_start",
    "hyoid_pelvic_sl1", "hyoid_pelvic_sl2", "hyoid_pelvic_sl3",
    "hyoid_pelvic_sl4", "hyoid_pelvic_sl5", "hyoid_pelvic_sl6",
    "hyoid_pelvic_sl7", "hyoid_pelvic_sl8",
    "hyoid_pelvic_end"
  )
)

curve_length <- function(mat) {
  seg <- mat[-1, , drop = FALSE] - mat[-nrow(mat), , drop = FALSE]
  sum(sqrt(rowSums(seg^2)))
}

chord_length <- function(mat) {
  sqrt(sum((mat[nrow(mat), ] - mat[1, ])^2))
}

curve_tortuosity <- function(mat) {
  cl <- chord_length(mat)
  if (isTRUE(all.equal(cl, 0))) return(NA_real_)
  curve_length(mat) / cl
}

max_deviation_from_chord <- function(mat) {
  p1 <- mat[1, ]
  p2 <- mat[nrow(mat), ]
  v <- p2 - p1
  v_norm <- sqrt(sum(v^2))
  if (isTRUE(all.equal(v_norm, 0))) return(NA_real_)
  
  perp_dist <- apply(mat, 1, function(p) {
    abs(v[1] * (p1[2] - p[2]) - (p1[1] - p[1]) * v[2]) / v_norm
  })
  
  max(perp_dist, na.rm = TRUE)
}

extract_curve_metrics_one <- function(specimen_id, coords_arr) {
  out <- lapply(names(curve_defs), function(curve_nm) {
    pts <- curve_defs[[curve_nm]]
    mat <- coords_arr[pts, , specimen_id, drop = FALSE][, , 1]
    
    tibble(
      specimen = specimen_id,
      curve = curve_nm,
      metric = c("tortuosity", "max_deviation"),
      value_raw = c(
        curve_tortuosity(mat),
        max_deviation_from_chord(mat)
      )
    )
  })
  
  bind_rows(out)
}

# ------------------------------------------------------------
# 5) Extract all measured traits for the 1950 habitat subset
# ------------------------------------------------------------
coords_pc_raw <- subset_coords_to_gdf_raw(coords_all, gdf_pc)
stopifnot(identical(dimnames(coords_pc_raw)[[3]], gdf_pc$specimen))

# Distance traits
measurement_cols <- c("Eye_width", "Body_depth", "Operculum_width", "Jaw_muscle_length")

distance_wide <- extract_distance_measurements(coords_pc_raw, gdf_pc)
distance_wide <- size_correct_measurements(distance_wide, measurement_cols)

distance_traits <- distance_wide %>%
  select(
    specimen, habitat, year, group,
    Eye_width_sizecorr,
    Body_depth_sizecorr,
    Operculum_width_sizecorr,
    Jaw_muscle_length_sizecorr
  ) %>%
  rename(
    Eye_width = Eye_width_sizecorr,
    Body_depth = Body_depth_sizecorr,
    Operculum_width = Operculum_width_sizecorr,
    Jaw_muscle_length = Jaw_muscle_length_sizecorr
  )

# Curve metrics
curve_metrics_raw <- map_dfr(dimnames(coords_pc_raw)[[3]], extract_curve_metrics_one, coords_arr = coords_pc_raw)

curve_metrics_raw <- curve_metrics_raw %>%
  left_join(
    gdf_pc %>%
      select(specimen, habitat, year, group, logCsize),
    by = "specimen"
  )

curve_metrics_sc <- curve_metrics_raw %>%
  group_by(curve, metric) %>%
  group_modify(~{
    dat <- .x
    fit <- lm(value_raw ~ logCsize, data = dat)
    dat %>% mutate(value_size_corrected = resid(fit))
  }) %>%
  ungroup() %>%
  mutate(
    trait = case_when(
      curve == "snout" & metric == "tortuosity"    ~ "Snout_tortuosity",
      curve == "snout" & metric == "max_deviation" ~ "Snout_max_deviation",
      curve == "hyoid" & metric == "tortuosity"    ~ "Hyoid_tortuosity",
      curve == "hyoid" & metric == "max_deviation" ~ "Hyoid_max_deviation",
      TRUE ~ NA_character_
    )
  )

curve_traits <- curve_metrics_sc %>%
  select(specimen, trait, value_size_corrected) %>%
  pivot_wider(names_from = trait, values_from = value_size_corrected)

# Combined trait table
all_traits_wide <- scores_1950 %>%
  select(specimen, habitat, year, group, PC1, PC2, PC3, PC4) %>%
  left_join(distance_traits, by = c("specimen", "habitat", "year", "group")) %>%
  left_join(curve_traits, by = "specimen")

trait_names <- c(
  "Eye_width",
  "Body_depth",
  "Operculum_width",
  "Jaw_muscle_length",
  "Snout_tortuosity",
  "Snout_max_deviation",
  "Hyoid_tortuosity",
  "Hyoid_max_deviation"
)

# ------------------------------------------------------------
# 6) Relate each trait to PC1-PC4
# ------------------------------------------------------------
run_pc_trait_lm <- function(df, pc_name, trait_name) {
  dat <- df %>%
    select(specimen, group, all_of(pc_name), all_of(trait_name)) %>%
    rename(pc = all_of(pc_name), trait = all_of(trait_name)) %>%
    filter(is.finite(pc), is.finite(trait))
  
  if (nrow(dat) < 3) {
    return(tibble(
      PC = pc_name,
      trait = trait_name,
      n = nrow(dat),
      slope = NA_real_,
      intercept = NA_real_,
      r_squared = NA_real_,
      adj_r_squared = NA_real_,
      p_value = NA_real_,
      cor = NA_real_,
      abs_cor = NA_real_
    ))
  }
  
  fit <- lm(pc ~ trait, data = dat)
  sm <- summary(fit)
  coef_tab <- sm$coefficients
  cor_val <- suppressWarnings(cor(dat$pc, dat$trait, use = "complete.obs"))
  
  tibble(
    PC = pc_name,
    trait = trait_name,
    n = nrow(dat),
    slope = unname(coef(fit)["trait"]),
    intercept = unname(coef(fit)["(Intercept)"]),
    r_squared = unname(sm$r.squared),
    adj_r_squared = unname(sm$adj.r.squared),
    p_value = unname(coef_tab["trait", "Pr(>|t|)"]),
    cor = cor_val,
    abs_cor = abs(cor_val)
  )
}

pc_names <- c("PC1", "PC2", "PC3", "PC4")

pc_trait_results <- purrr::map_dfr(pc_names, function(pc) {
  purrr::map_dfr(trait_names, function(tr) run_pc_trait_lm(all_traits_wide, pc, tr))
}) %>%
  arrange(PC, desc(r_squared), p_value)

pc_trait_cor_table <- pc_trait_results %>%
  select(PC, trait, n, cor, abs_cor, r_squared, adj_r_squared, p_value)

pc_trait_r2_wide <- pc_trait_results %>%
  select(PC, trait, r_squared) %>%
  pivot_wider(names_from = PC, values_from = r_squared)

pc_trait_cor_wide <- pc_trait_results %>%
  select(PC, trait, cor) %>%
  pivot_wider(names_from = PC, values_from = cor)

# ------------------------------------------------------------
# 7) Save tables
# ------------------------------------------------------------
write.csv(scores_1950, file.path(out_dir, "pc_scores_1950_habitats.csv"), row.names = FALSE)
write.csv(var_table, file.path(out_dir, "pc_variance_explained_1950.csv"), row.names = FALSE)
write.csv(all_traits_wide, file.path(out_dir, "traits_and_pc_scores_specimen_level_1950.csv"), row.names = FALSE)
write.csv(pc_trait_results, file.path(out_dir, "pc_trait_lm_results_1950.csv"), row.names = FALSE)
write.csv(pc_trait_cor_table, file.path(out_dir, "pc_trait_correlation_summary_1950.csv"), row.names = FALSE)
write.csv(pc_trait_r2_wide, file.path(out_dir, "pc_trait_r2_wide_1950.csv"), row.names = FALSE)
write.csv(pc_trait_cor_wide, file.path(out_dir, "pc_trait_cor_wide_1950.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 8) Heatmaps for quick interpretation
# ------------------------------------------------------------
heat_cor <- pc_trait_results %>%
  mutate(
    trait = factor(trait, levels = rev(trait_names)),
    PC = factor(PC, levels = pc_names)
  ) %>%
  ggplot(aes(x = PC, y = trait, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", cor)), size = 3) +
  scale_fill_gradient2(low = "steelblue4", mid = "white", high = "firebrick3", midpoint = 0) +
  labs(title = "Signed correlation of traits with PC scores (1950 habitats)", x = NULL, y = NULL, fill = "r") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank())

heat_r2 <- pc_trait_results %>%
  mutate(
    trait = factor(trait, levels = rev(trait_names)),
    PC = factor(PC, levels = pc_names)
  ) %>%
  ggplot(aes(x = PC, y = trait, fill = r_squared)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", r_squared)), size = 3) +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  labs(title = expression(R^2 ~ "for trait" ~ rightleftharpoons ~ "PC models (1950 habitats)"), x = NULL, y = NULL, fill = expression(R^2)) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank())

ggsave(file.path(plot_dir, "heatmap_trait_pc_correlations_1950.png"), heat_cor, width = 7.5, height = 5.5, dpi = 400)
ggsave(file.path(plot_dir, "heatmap_trait_pc_correlations_1950.pdf"), heat_cor, width = 7.5, height = 5.5)
ggsave(file.path(plot_dir, "heatmap_trait_pc_r2_1950.png"), heat_r2, width = 7.5, height = 5.5, dpi = 400)
ggsave(file.path(plot_dir, "heatmap_trait_pc_r2_1950.pdf"), heat_r2, width = 7.5, height = 5.5)

# ------------------------------------------------------------
# 9) Scatterplots for strongest trait association per PC
# ------------------------------------------------------------
top_trait_per_pc <- pc_trait_results %>%
  group_by(PC) %>%
  slice_max(order_by = r_squared, n = 1, with_ties = FALSE) %>%
  ungroup()

make_pc_trait_scatter <- function(df, pc_name, trait_name, subtitle_text = NULL) {
  plot_dat <- df %>%
    select(specimen, group, all_of(pc_name), all_of(trait_name)) %>%
    rename(pc = all_of(pc_name), trait = all_of(trait_name)) %>%
    filter(is.finite(pc), is.finite(trait))
  
  ggplot(plot_dat, aes(x = trait, y = pc)) +
    geom_point(aes(shape = group), size = 2, alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(
      title = paste0(pc_name, " vs ", trait_name, " (1950 habitats)"),
      subtitle = subtitle_text,
      x = trait_name,
      y = pc_name,
      shape = "Group"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")
}

for (i in seq_len(nrow(top_trait_per_pc))) {
  pc_i <- top_trait_per_pc$PC[i]
  trait_i <- top_trait_per_pc$trait[i]
  subtitle_i <- paste0(
    "Top single-trait association for ", pc_i,
    "; r = ", sprintf("%.2f", top_trait_per_pc$cor[i]),
    ", R2 = ", sprintf("%.2f", top_trait_per_pc$r_squared[i]),
    ", p = ", formatC(top_trait_per_pc$p_value[i], format = "e", digits = 2)
  )
  
  p_i <- make_pc_trait_scatter(all_traits_wide, pc_i, trait_i, subtitle_i)
  ggsave(file.path(plot_dir, paste0(pc_i, "_top_trait_scatter_1950.png")), p_i, width = 7.5, height = 5.5, dpi = 400)
  ggsave(file.path(plot_dir, paste0(pc_i, "_top_trait_scatter_1950.pdf")), p_i, width = 7.5, height = 5.5)
}

# ------------------------------------------------------------
# 10) Optional: all trait-by-PC scatterplots
# ------------------------------------------------------------
all_scatter_dir <- file.path(plot_dir, "all_pc_trait_scatterplots_1950")
if (!dir.exists(all_scatter_dir)) dir.create(all_scatter_dir, recursive = TRUE)

for (pc_i in pc_names) {
  for (trait_i in trait_names) {
    row_i <- pc_trait_results %>% filter(PC == pc_i, trait == trait_i)
    subtitle_i <- paste0(
      "r = ", sprintf("%.2f", row_i$cor),
      "; R2 = ", sprintf("%.2f", row_i$r_squared),
      "; p = ", formatC(row_i$p_value, format = "e", digits = 2)
    )
    p_i <- make_pc_trait_scatter(all_traits_wide, pc_i, trait_i, subtitle_i)
    file_stub <- paste0(pc_i, "__", trait_i, "_1950")
    ggsave(file.path(all_scatter_dir, paste0(file_stub, ".png")), p_i, width = 7, height = 5.25, dpi = 300)
  }
}

# ------------------------------------------------------------
# 11) Console summary
# ------------------------------------------------------------
cat("\nSaved outputs to:\n", normalizePath(out_dir), "\n")
cat("\nKey files created:\n")
cat("  - pc_scores_1950_habitats.csv\n")
cat("  - pc_variance_explained_1950.csv\n")
cat("  - traits_and_pc_scores_specimen_level_1950.csv\n")
cat("  - pc_trait_lm_results_1950.csv\n")
cat("  - pc_trait_correlation_summary_1950.csv\n")
cat("  - pc_trait_r2_wide_1950.csv\n")
cat("  - pc_trait_cor_wide_1950.csv\n")
cat("  - plots/heatmap_trait_pc_correlations_1950.(png/pdf)\n")
cat("  - plots/heatmap_trait_pc_r2_1950.(png/pdf)\n")
cat("  - plots/PC1_top_trait_scatter_1950.(png/pdf) ... PC4_top_trait_scatter_1950.(png/pdf)\n")
cat("  - plots/all_pc_trait_scatterplots_1950/*.png\n")

cat("\nTop trait association per PC:\n")
print(top_trait_per_pc)