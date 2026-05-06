# ============================================================
# Scripts/SICB_curves_distances_figures.R
#
# Recreate curve metric and landmark-distance figures
# for SICB formatting.
#
# Saves:
#   Figures/SICB_curves_distances/
#     curve_shape_metrics_habitat_1950_SICB.pdf
#     curve_shape_metrics_CT_timeseries_SICB.pdf
#     Fig_distance_measurements_1950_habitats_SICB.pdf
#     Fig_distance_measurements_CT_timeseries_SICB.pdf
#
# All PDFs are 4 inches wide.
# ============================================================

# ---------------------------
# Load canonical project objects
# ---------------------------
source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")
source("R/methods/01_specimen_sampling_study_design/02_subset_1950.R")
source("R/methods/01_specimen_sampling_study_design/03_subset_CT_timeseries.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(purrr)
})

# ---------------------------
# Output directory
# ---------------------------
OUTDIR <- file.path("Figures", "SICB_curves_distances")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# PART 1: CURVE SHAPE METRICS
# Based on: curve_shape_metrics_raw_then_sizecorrected.R
# ============================================================

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

all_needed_pts <- unique(unlist(curve_defs))
missing_pts <- setdiff(all_needed_pts, dimnames(coords_all)[[1]])

if (length(missing_pts) > 0) {
  stop("These curve landmarks are missing from coords_all:\n",
       paste(missing_pts, collapse = ", "))
}

gdf_raw <- gdf %>%
  filter(specimen %in% dimnames(coords_all)[[3]]) %>%
  slice(match(dimnames(coords_all)[[3]], specimen))

stopifnot(identical(gdf_raw$specimen, dimnames(coords_all)[[3]]))

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

extract_curve_metrics_one <- function(specimen_id) {
  out <- lapply(names(curve_defs), function(curve_nm) {
    pts <- curve_defs[[curve_nm]]
    mat <- coords_all[pts, , specimen_id, drop = FALSE][, , 1]
    
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

curve_metrics_raw <- map_dfr(dimnames(coords_all)[[3]], extract_curve_metrics_one)

curve_metrics_raw <- curve_metrics_raw %>%
  left_join(
    gdf_raw %>%
      select(specimen, habitat, year, logCsize),
    by = "specimen"
  )

curve_metrics_sc <- curve_metrics_raw %>%
  group_by(curve, metric) %>%
  group_modify(~{
    dat <- .x
    fit <- lm(value_raw ~ logCsize, data = dat)
    
    dat %>%
      mutate(
        value_size_corrected = resid(fit),
        value_plot = value_size_corrected
      )
  }) %>%
  ungroup() %>%
  mutate(
    curve_label = recode(
      curve,
      "snout" = "Snout",
      "hyoid" = "Hyoid"
    ),
    metric_label = recode(
      metric,
      "tortuosity" = "Tortuosity",
      "max_deviation" = "Max deviation from chord"
    )
  )

habitat_levels <- c(
  "Connecticut River",
  "Quabbin",
  "Swift River",
  "Fort River",
  "Sawmill River"
)

habitat_dat <- curve_metrics_sc %>%
  filter(
    !is.na(habitat),
    year == 1950,
    habitat %in% habitat_levels
  ) %>%
  mutate(group = factor(habitat, levels = habitat_levels))

ct_year_levels <- c(1950, 1956, 1970)

ct_dat <- curve_metrics_sc %>%
  filter(
    habitat == "Connecticut River",
    year %in% ct_year_levels
  ) %>%
  mutate(group = factor(year, levels = ct_year_levels))

make_curve_metric_plot_SICB <- function(dat, xlab) {
  ggplot(dat, aes(x = group, y = value_plot)) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.25,
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_boxplot(
      width = 0.55,
      outlier.shape = NA,
      fill = "white",
      color = "black",
      linewidth = 0.25
    ) +
    geom_jitter(
      width = 0.10,
      height = 0,
      alpha = 0.8,
      size = 0.8
    ) +
    facet_grid(metric_label ~ curve_label, scales = "free_y") +
    labs(
      x = xlab,
      y = "Size-corrected, mean-centered value"
    ) +
    theme_bw(base_size = 6, base_family = "Arial") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.15),
      panel.border = element_rect(linewidth = 0.25),
      strip.background = element_rect(fill = "grey95", linewidth = 0.25),
      strip.text = element_text(size = 6),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 5),
      axis.text.x = element_text(angle = 25, hjust = 1),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      plot.margin = margin(2, 2, 2, 2, unit = "pt")
    )
}

p_curve_habitat <- make_curve_metric_plot_SICB(
  habitat_dat,
  xlab = "Habitat"
)

p_curve_ct <- make_curve_metric_plot_SICB(
  ct_dat,
  xlab = "Year"
)

ggsave(
  filename = file.path(OUTDIR, "curve_shape_metrics_habitat_1950_SICB.pdf"),
  plot = p_curve_habitat,
  width = 4,
  height = 3.2,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(OUTDIR, "curve_shape_metrics_CT_timeseries_SICB.pdf"),
  plot = p_curve_ct,
  width = 4,
  height = 3.2,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

# ============================================================
# PART 2: LANDMARK DISTANCE MEASUREMENTS
# Based on: landmark_distance_measurements.R
# ============================================================

stopifnot(identical(dimnames(coords_all)[[3]], gdf$specimen))

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
  
  if (!(landmark_a %in% pt_names)) {
    stop("Landmark not found in coords_arr: ", landmark_a)
  }
  if (!(landmark_b %in% pt_names)) {
    stop("Landmark not found in coords_arr: ", landmark_b)
  }
  
  A <- coords_arr[landmark_a, , , drop = FALSE]
  B <- coords_arr[landmark_b, , , drop = FALSE]
  
  ax <- as.numeric(A[1, 1, ])
  ay <- as.numeric(A[1, 2, ])
  bx <- as.numeric(B[1, 1, ])
  by <- as.numeric(B[1, 2, ])
  
  sqrt((ax - bx)^2 + (ay - by)^2)
}

extract_measurements <- function(coords_arr, gdf_meta) {
  out <- gdf_meta %>%
    dplyr::select(specimen, habitat, year, size_for_allometry, size_label)
  
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

make_long_measure_df <- function(df_wide, group_var, corrected = TRUE) {
  suffix <- if (corrected) "_sizecorr" else ""
  
  measure_cols <- c(
    paste0("Eye_width", suffix),
    paste0("Body_depth", suffix),
    paste0("Operculum_width", suffix),
    paste0("Jaw_muscle_length", suffix)
  )
  
  df_wide %>%
    dplyr::select(specimen, habitat, year, !!rlang::sym(group_var), all_of(measure_cols)) %>%
    tidyr::pivot_longer(
      cols = all_of(measure_cols),
      names_to = "measurement",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      measurement = stringr::str_remove(measurement, "_sizecorr$"),
      measurement = factor(
        measurement,
        levels = c("Eye_width", "Body_depth", "Operculum_width", "Jaw_muscle_length"),
        labels = c("Eye width", "Body depth", "Operculum width", "Jaw muscle length")
      )
    )
}

plot_measure_boxpoints_SICB <- function(df_long, group_var) {
  ggplot(df_long, aes(x = .data[[group_var]], y = value)) +
    geom_boxplot(
      width = 0.6,
      outlier.shape = NA,
      fill = "grey85",
      color = "black",
      linewidth = 0.25
    ) +
    geom_point(
      position = position_jitter(width = 0.10, height = 0),
      size = 0.8,
      alpha = 0.8
    ) +
    facet_wrap(~ measurement, scales = "free_y", ncol = 2) +
    labs(
      x = NULL,
      y = "Size-corrected distance"
    ) +
    theme_classic(base_family = "Arial", base_size = 6) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 6, face = "bold"),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 5),
      axis.text.x = element_text(angle = 25, hjust = 1),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      plot.margin = margin(2, 2, 2, 2, unit = "pt")
    )
}

gdf_1950_dist <- gdf %>%
  dplyr::filter(year == 1950, !is.na(habitat)) %>%
  dplyr::mutate(
    habitat = factor(habitat, levels = habitat_levels)
  ) %>%
  droplevels()

coords_1950_raw <- subset_coords_to_gdf_raw(coords_all, gdf_1950_dist)

dist_1950_wide <- extract_measurements(coords_1950_raw, gdf_1950_dist)

measurement_cols <- c(
  "Eye_width",
  "Body_depth",
  "Operculum_width",
  "Jaw_muscle_length"
)

dist_1950_wide <- size_correct_measurements(dist_1950_wide, measurement_cols)

dist_1950_long <- dist_1950_wide %>%
  dplyr::mutate(
    group = factor(habitat, levels = habitat_levels)
  ) %>%
  make_long_measure_df(group_var = "group", corrected = TRUE)

p_dist_1950 <- plot_measure_boxpoints_SICB(
  dist_1950_long,
  group_var = "group"
)

gdf_CT_dist <- gdf %>%
  dplyr::filter(
    habitat == "Connecticut River",
    year %in% c(1950, 1956, 1970)
  ) %>%
  dplyr::mutate(
    group = factor(
      paste0("CT_", year),
      levels = c("CT_1950", "CT_1956", "CT_1970")
    )
  ) %>%
  droplevels()

coords_CT_raw <- subset_coords_to_gdf_raw(coords_all, gdf_CT_dist)

dist_CT_wide <- extract_measurements(coords_CT_raw, gdf_CT_dist)

dist_CT_wide <- size_correct_measurements(dist_CT_wide, measurement_cols)

dist_CT_long <- dist_CT_wide %>%
  dplyr::mutate(
    group = factor(
      paste0("CT_", year),
      levels = c("CT_1950", "CT_1956", "CT_1970")
    )
  ) %>%
  make_long_measure_df(group_var = "group", corrected = TRUE)

p_dist_CT <- plot_measure_boxpoints_SICB(
  dist_CT_long,
  group_var = "group"
)

ggsave(
  filename = file.path(OUTDIR, "Fig_distance_measurements_1950_habitats_SICB.pdf"),
  plot = p_dist_1950,
  width = 4,
  height = 3.2,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(OUTDIR, "Fig_distance_measurements_CT_timeseries_SICB.pdf"),
  plot = p_dist_CT,
  width = 4,
  height = 3.2,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)



# ============================================================
# PART 3: COMBINED 1950 HABITATS + CT TIMEPOINTS
# Curves + landmark distances
#
# Order:
#   CT 1970, CT 1956, CT 1950, Quabbin, Swift, Fort, Sawmill
#
# Curves: remove tortuosity; keep max deviation only
# ============================================================

combined_levels <- c(
  "CT 1970",
  "CT 1956",
  "CT 1950",
  "Quabbin",
  "Swift",
  "Fort",
  "Sawmill"
)

# ---------------------------
# Combined curve metrics
# ---------------------------

combined_curve_dat <- curve_metrics_sc %>%
  dplyr::filter(
    metric == "max_deviation",
    (
      habitat == "Connecticut River" & year %in% c(1950, 1956, 1970)
    ) |
      (
        year == 1950 & habitat %in% c(
          "Quabbin",
          "Swift River",
          "Fort River",
          "Sawmill River"
        )
      )
  ) %>%
  dplyr::mutate(
    group = dplyr::case_when(
      habitat == "Connecticut River" & year == 1970 ~ "CT 1970",
      habitat == "Connecticut River" & year == 1956 ~ "CT 1956",
      habitat == "Connecticut River" & year == 1950 ~ "CT 1950",
      habitat == "Quabbin" ~ "Quabbin",
      habitat == "Swift River" ~ "Swift",
      habitat == "Fort River" ~ "Fort",
      habitat == "Sawmill River" ~ "Sawmill",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = combined_levels)
  ) %>%
  dplyr::filter(!is.na(group))

p_curve_combined <- ggplot(combined_curve_dat, aes(x = group, y = value_plot)) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.25,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.25
  ) +
  geom_jitter(
    width = 0.10,
    height = 0,
    alpha = 0.8,
    size = 0.8
  ) +
  facet_wrap(~ curve_label, scales = "free_y", ncol = 2) +
  labs(
    x = NULL,
    y = "Size-corrected, mean-centered max deviation"
  ) +
  theme_bw(base_size = 6, base_family = "Arial") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.15),
    panel.border = element_rect(linewidth = 0.25),
    strip.background = element_rect(fill = "grey95", linewidth = 0.25),
    strip.text = element_text(size = 6),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 5),
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.ticks = element_line(linewidth = 0.25),
    plot.margin = margin(2, 2, 2, 2, unit = "pt")
  )

ggsave(
  filename = file.path(
    OUTDIR,
    "curve_shape_metrics_combined_CT_timepoints_1950_habitats_no_tortuosity_SICB.pdf"
  ),
  plot = p_curve_combined,
  width = 4,
  height = 2.6,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

# ---------------------------
# Combined landmark distances
# ---------------------------

gdf_combined_dist <- gdf %>%
  dplyr::filter(
    (
      habitat == "Connecticut River" & year %in% c(1950, 1956, 1970)
    ) |
      (
        year == 1950 & habitat %in% c(
          "Quabbin",
          "Swift River",
          "Fort River",
          "Sawmill River"
        )
      )
  ) %>%
  dplyr::mutate(
    group = dplyr::case_when(
      habitat == "Connecticut River" & year == 1970 ~ "CT 1970",
      habitat == "Connecticut River" & year == 1956 ~ "CT 1956",
      habitat == "Connecticut River" & year == 1950 ~ "CT 1950",
      habitat == "Quabbin" ~ "Quabbin",
      habitat == "Swift River" ~ "Swift",
      habitat == "Fort River" ~ "Fort",
      habitat == "Sawmill River" ~ "Sawmill",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = combined_levels)
  ) %>%
  dplyr::filter(!is.na(group)) %>%
  droplevels()

coords_combined_raw <- subset_coords_to_gdf_raw(coords_all, gdf_combined_dist)

dist_combined_wide <- extract_measurements(
  coords_combined_raw,
  gdf_combined_dist
)

# 🔑 ADD GROUP BACK IN
dist_combined_wide$group <- gdf_combined_dist$group

dist_combined_wide <- size_correct_measurements(
  dist_combined_wide,
  measurement_cols
)

dist_combined_long <- dist_combined_wide %>%
  make_long_measure_df(group_var = "group", corrected = TRUE)

p_dist_combined <- ggplot(dist_combined_long, aes(x = group, y = value)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    fill = "grey85",
    color = "black",
    linewidth = 0.25
  ) +
  geom_point(
    position = position_jitter(width = 0.10, height = 0),
    size = 0.8,
    alpha = 0.8
  ) +
  facet_wrap(~ measurement, scales = "free_y", ncol = 2) +
  labs(
    x = NULL,
    y = "Size-corrected distance"
  ) +
  theme_classic(base_family = "Arial", base_size = 6) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 6, face = "bold"),
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 5),
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.line = element_line(linewidth = 0.25),
    axis.ticks = element_line(linewidth = 0.25),
    plot.margin = margin(2, 2, 2, 2, unit = "pt")
  )

ggsave(
  filename = file.path(
    OUTDIR,
    "Fig_distance_measurements_combined_CT_timepoints_1950_habitats_SICB.pdf"
  ),
  plot = p_dist_combined,
  width = 4,
  height = 3.2,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)



# ============================================================
# PART 4: SICB STATS FOR COMBINED GROUPS
# Pairwise comparisons:
# CT 1970, CT 1956, CT 1950, Quabbin, Swift, Fort, Sawmill
# ============================================================

STAT_OUTDIR <- file.path("Outputs", "distance_curves_stats_SICB")
dir.create(STAT_OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Generic stats helper
# ---------------------------

run_pairwise_stats_SICB <- function(dat, value_col, facet_cols, dataset_label) {
  
  value_sym <- rlang::sym(value_col)
  
  combos <- dat %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(facet_cols)))
  
  anova_list <- list()
  pairwise_list <- list()
  summary_list <- list()
  
  for (i in seq_len(nrow(combos))) {
    
    sub <- dat
    
    for (fc in facet_cols) {
      sub <- sub %>%
        dplyr::filter(.data[[fc]] == combos[[fc]][i])
    }
    
    sub <- sub %>%
      dplyr::filter(is.finite(!!value_sym), !is.na(group)) %>%
      dplyr::mutate(group = factor(group, levels = combined_levels)) %>%
      droplevels()
    
    if (nrow(sub) < 3 || dplyr::n_distinct(sub$group) < 2) next
    
    # Summary stats
    summary_out <- sub %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        n = sum(is.finite(!!value_sym)),
        mean = mean(!!value_sym, na.rm = TRUE),
        sd = sd(!!value_sym, na.rm = TRUE),
        min = min(!!value_sym, na.rm = TRUE),
        median = median(!!value_sym, na.rm = TRUE),
        max = max(!!value_sym, na.rm = TRUE),
        .groups = "drop"
      )
    
    for (fc in facet_cols) summary_out[[fc]] <- combos[[fc]][i]
    summary_out$dataset <- dataset_label
    summary_list[[length(summary_list) + 1]] <- summary_out
    
    # ANOVA
    fit <- aov(stats::as.formula(paste(value_col, "~ group")), data = sub)
    aov_tab <- as.data.frame(summary(fit)[[1]])
    aov_tab$term <- rownames(aov_tab)
    rownames(aov_tab) <- NULL
    
    for (fc in facet_cols) aov_tab[[fc]] <- combos[[fc]][i]
    aov_tab$dataset <- dataset_label
    
    anova_list[[length(anova_list) + 1]] <- aov_tab
    
    # All pairwise comparisons: pairwise t-test, BH adjusted
    pw <- pairwise.t.test(
      x = sub[[value_col]],
      g = sub$group,
      p.adjust.method = "BH",
      pool.sd = FALSE
    )
    
    pw_df <- as.data.frame(as.table(pw$p.value))
    colnames(pw_df) <- c("group1", "group2", "p_BH")
    
    pw_df <- pw_df %>%
      dplyr::filter(!is.na(p_BH)) %>%
      dplyr::mutate(
        dataset = dataset_label,
        p_adjust_method = "BH",
        test = "pairwise.t.test"
      )
    
    for (fc in facet_cols) pw_df[[fc]] <- combos[[fc]][i]
    
    pairwise_list[[length(pairwise_list) + 1]] <- pw_df
  }
  
  list(
    summary = dplyr::bind_rows(summary_list),
    anova = dplyr::bind_rows(anova_list),
    pairwise = dplyr::bind_rows(pairwise_list)
  )
}

# ---------------------------
# Curve stats: max deviation only, by curve
# ---------------------------

curve_stats_SICB <- run_pairwise_stats_SICB(
  dat = combined_curve_dat,
  value_col = "value_plot",
  facet_cols = c("curve_label"),
  dataset_label = "combined_curve_max_deviation"
)

write.csv(
  curve_stats_SICB$summary,
  file.path(STAT_OUTDIR, "combined_curve_max_deviation_summary_SICB.csv"),
  row.names = FALSE
)

write.csv(
  curve_stats_SICB$anova,
  file.path(STAT_OUTDIR, "combined_curve_max_deviation_anova_SICB.csv"),
  row.names = FALSE
)

write.csv(
  curve_stats_SICB$pairwise,
  file.path(STAT_OUTDIR, "combined_curve_max_deviation_pairwise_BH_SICB.csv"),
  row.names = FALSE
)

# ---------------------------
# Distance stats: by linear measurement
# ---------------------------

distance_stats_SICB <- run_pairwise_stats_SICB(
  dat = dist_combined_long,
  value_col = "value",
  facet_cols = c("measurement"),
  dataset_label = "combined_landmark_distances"
)

write.csv(
  distance_stats_SICB$summary,
  file.path(STAT_OUTDIR, "combined_landmark_distances_summary_SICB.csv"),
  row.names = FALSE
)

write.csv(
  distance_stats_SICB$anova,
  file.path(STAT_OUTDIR, "combined_landmark_distances_anova_SICB.csv"),
  row.names = FALSE
)

write.csv(
  distance_stats_SICB$pairwise,
  file.path(STAT_OUTDIR, "combined_landmark_distances_pairwise_BH_SICB.csv"),
  row.names = FALSE
)

cat("\nSaved SICB stats to:\n")
cat(normalizePath(STAT_OUTDIR), "\n")



cat("\nSaved SICB stats to:\n")
cat(normalizePath(STAT_OUTDIR), "\n")
cat("\nStats files created:\n")
print(list.files(STAT_OUTDIR))


# ============================================================
# PART 5: GRAPH-ONLY NON-SIGNIFICANT PAIRWISE NETWORKS
# Edges = pairwise comparisons with p_BH >= 0.05
# No PCA space
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

NETWORK_OUTDIR <- file.path("Figures", "SICB_curves_distances", "graph_only_networks")
dir.create(NETWORK_OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Fixed node layout
# ---------------------------
# This puts groups in your preferred order, arranged in a circle.
# You can manually tweak x/y later if desired.

node_layout <- data.frame(
  group = factor(combined_levels, levels = combined_levels),
  angle = seq(pi / 2, pi / 2 - 2 * pi + 2 * pi / length(combined_levels),
              length.out = length(combined_levels))
) %>%
  dplyr::mutate(
    x = cos(angle),
    y = sin(angle)
  )

# ---------------------------
# Helper: create graph edges from pairwise table
# ---------------------------

make_graph_edges <- function(pairwise_dat, trait_col, trait_value, alpha = 0.05) {
  
  pairwise_dat %>%
    dplyr::filter(
      .data[[trait_col]] == trait_value,
      p_BH >= alpha
    ) %>%
    dplyr::transmute(
      from = factor(group1, levels = combined_levels),
      to   = factor(group2, levels = combined_levels),
      p_BH = p_BH
    ) %>%
    left_join(node_layout, by = c("from" = "group")) %>%
    rename(x = x, y = y) %>%   # REMOVE THIS LINE COMPLETELY
    left_join(node_layout, by = c("to" = "group"), suffix = c("", "_end")) %>%
    rename(xend = x_end, yend = y_end)
}

# ---------------------------
# Helper: plot one graph-only network
# ---------------------------

plot_graph_only_network <- function(edges_plot, title_text, file_stub) {
  
  p <- ggplot() +
    geom_segment(
      data = edges_plot,
      aes(x = x, y = y, xend = xend, yend = yend),
      color = "grey35",
      linewidth = 0.35,
      alpha = 0.65
    ) +
    geom_point(
      data = node_layout,
      aes(x = x, y = y),
      shape = 21,
      fill = "white",
      color = "black",
      size = 4,
      stroke = 0.4
    ) +
    geom_text(
      data = node_layout,
      aes(x = x, y = y, label = group),
      size = 2.4,
      vjust = -1.15
    ) +
    coord_equal(clip = "off") +
    labs(
      title = title_text,
      subtitle = "Edges connect groups that are not significantly different (BH-adjusted p ≥ 0.05)",
      x = NULL,
      y = NULL
    ) +
    theme_void(base_size = 6, base_family = "Arial") +
    theme(
      plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 5.5, hjust = 0.5),
      plot.margin = margin(8, 8, 8, 8, unit = "pt")
    )
  
  ggsave(
    filename = file.path(NETWORK_OUTDIR, paste0(file_stub, "_graph_network_SICB.pdf")),
    plot = p,
    width = 4,
    height = 3.2,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )
  
  p
}

# ============================================================
# Distance networks: one PDF per measurement
# ============================================================

distance_pairwise_df <- distance_stats_SICB$pairwise

for (m in unique(distance_pairwise_df$measurement)) {
  
  edges_m <- make_graph_edges(
    pairwise_dat = distance_pairwise_df,
    trait_col = "measurement",
    trait_value = m,
    alpha = 0.05
  )
  
  plot_graph_only_network(
    edges_plot = edges_m,
    title_text = paste("Non-significant network:", m),
    file_stub = paste0("distance_", gsub("[^A-Za-z0-9]+", "_", as.character(m)))
  )
}

# ============================================================
# Curve networks: one PDF per curve
# ============================================================

curve_pairwise_df <- curve_stats_SICB$pairwise

for (c in unique(curve_pairwise_df$curve_label)) {
  
  edges_c <- make_graph_edges(
    pairwise_dat = curve_pairwise_df,
    trait_col = "curve_label",
    trait_value = c,
    alpha = 0.05
  )
  
  plot_graph_only_network(
    edges_plot = edges_c,
    title_text = paste("Non-significant network:", c, "max deviation"),
    file_stub = paste0("curve_", gsub("[^A-Za-z0-9]+", "_", as.character(c)))
  )
}

cat("\nSaved graph-only networks to:\n")
cat(normalizePath(NETWORK_OUTDIR), "\n")
print(list.files(NETWORK_OUTDIR))

# ============================================================
# PART 6: Compact letter displays from BH-adjusted pairwise tests
# Uses multcompView on pairwise.t.test BH p-values
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(multcompView)
})

LETTER_OUTDIR <- file.path("Outputs", "distance_curves_stats_SICB")
dir.create(LETTER_OUTDIR, recursive = TRUE, showWarnings = FALSE)

make_letters_from_pairwise <- function(pairwise_df, trait_col, alpha = 0.05) {
  
  traits <- unique(pairwise_df[[trait_col]])
  
  letter_list <- lapply(traits, function(trait_value) {
    
    sub <- pairwise_df %>%
      filter(.data[[trait_col]] == trait_value)
    
    # Create named vector of TRUE/FALSE significant comparisons
    sig_vec <- sub$p_BH < alpha
    
    names(sig_vec) <- paste(sub$group1, sub$group2, sep = "-")
    
    # multcompLetters expects TRUE = significantly different
    letters <- multcompView::multcompLetters(sig_vec)
    
    tibble(
      !!trait_col := trait_value,
      group = names(letters$Letters),
      letter = letters$Letters
    ) %>%
      mutate(
        group = factor(group, levels = combined_levels)
      ) %>%
      arrange(group)
  })
  
  bind_rows(letter_list)
}

# ---------------------------
# Landmark distance letters
# ---------------------------

distance_letters_SICB <- make_letters_from_pairwise(
  pairwise_df = distance_stats_SICB$pairwise,
  trait_col = "measurement",
  alpha = 0.05
)

write.csv(
  distance_letters_SICB,
  file.path(LETTER_OUTDIR, "combined_landmark_distances_letters_BH_SICB.csv"),
  row.names = FALSE
)

print(distance_letters_SICB)

# ---------------------------
# Curve max deviation letters
# ---------------------------

curve_letters_SICB <- make_letters_from_pairwise(
  pairwise_df = curve_stats_SICB$pairwise,
  trait_col = "curve_label",
  alpha = 0.05
)

write.csv(
  curve_letters_SICB,
  file.path(LETTER_OUTDIR, "combined_curve_max_deviation_letters_BH_SICB.csv"),
  row.names = FALSE
)

print(curve_letters_SICB)

# ============================================================
# Final message
# ============================================================

cat("\nSaved SICB PDFs to:\n")
cat(normalizePath(OUTDIR), "\n\n")

cat("Files created:\n")
cat("  - curve_shape_metrics_habitat_1950_SICB.pdf\n")
cat("  - curve_shape_metrics_CT_timeseries_SICB.pdf\n")
cat("  - Fig_distance_measurements_1950_habitats_SICB.pdf\n")
cat("  - Fig_distance_measurements_CT_timeseries_SICB.pdf\n")

cat("\nAdditional combined SICB PDFs created:\n")
cat("  - curve_shape_metrics_combined_CT_timepoints_1950_habitats_no_tortuosity_SICB.pdf\n")
cat("  - Fig_distance_measurements_combined_CT_timepoints_1950_habitats_SICB.pdf\n")
