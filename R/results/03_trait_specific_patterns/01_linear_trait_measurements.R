# ============================================================
# Scripts/landmark_distance_measurements.R
#
# Landmark-to-landmark distance measurements from RAW coordinates
# (pre-GPA, pre-size-standardization), followed by size correction.
#
# Outputs:
#   1) 1950 habitat comparison
#   2) Connecticut River time-series comparison
#
# Figures:
#   - barplots with points, faceted by measurement
#
# Tables:
#   - summary tables (n, mean, sd, min, median, max)
#   - raw individual-level tables
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
})

# ---------------------------
# Output directories
# ---------------------------
FIG_DIR <- "figures"
TAB_DIR <- file.path("Outputs", "distance_measurements")

if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)
if (!dir.exists(TAB_DIR)) dir.create(TAB_DIR, recursive = TRUE)

# ---------------------------
# Preconditions
# ---------------------------
if (!exists("coords_all")) stop("coords_all not found. Run R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R.")
if (!exists("gdf")) stop("gdf not found. Run R/methods/01_specimen_sampling_study_design/01_build_metadata.R.")

stopifnot(identical(dimnames(coords_all)[[3]], gdf$specimen))

# ---------------------------
# Helper: subset coords to metadata order
# ---------------------------
subset_coords_to_gdf_raw <- function(coords_array, gdf_meta) {
  ids_all <- dimnames(coords_array)[[3]]
  idx <- match(gdf_meta$specimen, ids_all)
  stopifnot(!any(is.na(idx)))
  
  out <- coords_array[, , idx, drop = FALSE]
  dimnames(out)[[3]] <- gdf_meta$specimen
  stopifnot(identical(dimnames(out)[[3]], gdf_meta$specimen))
  out
}

# ---------------------------
# Helper: compute Euclidean distance for one landmark pair
# ---------------------------
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

# ---------------------------
# Helper: apply landmark-distance set
# ---------------------------
extract_measurements <- function(coords_arr, gdf_meta) {
  
  out <- gdf_meta %>%
    dplyr::select(specimen, habitat, year, size_for_allometry, size_label)
  
  out$Eye_width <- get_distance_measure(coords_arr, "orbit_1", "orbit_2")
  out$Body_depth <- get_distance_measure(coords_arr, "premaxilla", "maxilla")
  out$Operculum_width <- get_distance_measure(coords_arr, "max_curve_preoperculum", "operculum")
  out$Jaw_muscle_length <- get_distance_measure(coords_arr, "max_curve_preoperculum", "preoperculum")
  
  out
}

# ---------------------------
# Helper: size-correct a set of measurements
# log(measure) ~ size_for_allometry
# then back-transform centered residuals into original units
# ---------------------------
size_correct_measurements <- function(df_wide, measurement_cols) {
  
  df_out <- df_wide
  
  for (m in measurement_cols) {
    if (any(df_out[[m]] <= 0, na.rm = TRUE)) {
      stop("Non-positive values found in measurement ", m, "; cannot log-transform.")
    }
    
    log_y <- log(df_out[[m]])
    fit <- lm(log_y ~ size_for_allometry, data = df_out)
    
    # centered residuals placed back on original-unit scale
    corrected <- exp(residuals(fit) + mean(log_y, na.rm = TRUE))
    
    df_out[[paste0(m, "_sizecorr")]] <- corrected
  }
  
  df_out
}

# ---------------------------
# Helper: wide -> long
# ---------------------------
make_long_measure_df <- function(df_wide, group_var, corrected = TRUE) {
  
  suffix <- if (corrected) "_sizecorr" else ""
  
  measure_cols <- c(
    paste0("Eye_width", suffix),
    paste0("Body_depth", suffix),
    paste0("Operculum_width", suffix),
    paste0("Jaw_muscle_length", suffix)
  )
  
  out <- df_wide %>%
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
  
  out
}

# ---------------------------
# Helper: summary table
# ---------------------------
make_summary_table <- function(df_long, group_var) {
  df_long %>%
    dplyr::group_by(.data[[group_var]], measurement) %>%
    dplyr::summarise(
      n = sum(!is.na(value)),
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      min = min(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      max = max(value, na.rm = TRUE),
      .groups = "drop"
    )
}

# ---------------------------
# Helper: plotting
# ---------------------------
plot_measure_boxpoints <- function(df_long, group_var, title_text, subtitle_text) {
  
  ggplot(df_long, aes(x = .data[[group_var]], y = value)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, fill = "grey85", color = "black") +
    geom_point(
      position = position_jitter(width = 0.12, height = 0),
      size = 1.8,
      alpha = 0.8
    ) +
    facet_wrap(~ measurement, scales = "free_y", ncol = 2) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = NULL,
      y = "Size-corrected distance"
    ) +
    theme_classic(base_family = "Helvetica", base_size = 10) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1)
    )
}

# ============================================================
# 1) 1950 habitat analysis
# ============================================================

habitat_order <- c(
  "Connecticut River",
  "Quabbin",
  "Swift River",
  "Fort River",
  "Sawmill River"
)

gdf_1950_dist <- gdf %>%
  dplyr::filter(year == 1950, !is.na(habitat)) %>%
  dplyr::mutate(
    habitat = factor(habitat, levels = habitat_order)
  ) %>%
  droplevels()

coords_1950_raw <- subset_coords_to_gdf_raw(coords_all, gdf_1950_dist)

dist_1950_wide <- extract_measurements(coords_1950_raw, gdf_1950_dist)

measurement_cols <- c("Eye_width", "Body_depth", "Operculum_width", "Jaw_muscle_length")
dist_1950_wide <- size_correct_measurements(dist_1950_wide, measurement_cols)

dist_1950_long <- dist_1950_wide %>%
  dplyr::mutate(
    group = factor(habitat, levels = habitat_order)
  ) %>%
  make_long_measure_df(group_var = "group", corrected = TRUE)

summary_1950 <- make_summary_table(dist_1950_long, group_var = "group")

p_1950 <- plot_measure_boxpoints(
  dist_1950_long,
  group_var = "group",
  title_text = "Landmark distance measurements by habitat (1950)",
  subtitle_text = paste0("Measured on raw coordinates, then size-corrected using ", unique(gdf_1950_dist$size_label)[1])
)

print(p_1950)

ggsave(
  file.path(FIG_DIR, "Fig_distance_measurements_1950_habitats.pdf"),
  p_1950, width = 8, height = 6, units = "in"
)
ggsave(
  file.path(FIG_DIR, "Fig_distance_measurements_1950_habitats.png"),
  p_1950, width = 8, height = 6, units = "in", dpi = 600
)

write.csv(
  dist_1950_wide,
  file.path(TAB_DIR, "distance_measurements_1950_habitats_individuals.csv"),
  row.names = FALSE
)

write.csv(
  summary_1950,
  file.path(TAB_DIR, "distance_measurements_1950_habitats_summary.csv"),
  row.names = FALSE
)

# ============================================================
# 2) Connecticut River time-series analysis
# ============================================================

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
  dplyr::mutate(group = factor(
    paste0("CT_", year),
    levels = c("CT_1950", "CT_1956", "CT_1970")
  )) %>%
  make_long_measure_df(group_var = "group", corrected = TRUE)

summary_CT <- make_summary_table(dist_CT_long, group_var = "group")

p_CT <- plot_measure_boxpoints(
  dist_CT_long,
  group_var = "group",
  title_text = "Landmark distance measurements within Connecticut River through time",
  subtitle_text = paste0("Measured on raw coordinates, then size-corrected using ", unique(gdf_CT_dist$size_label)[1])
)

print(p_CT)

ggsave(
  file.path(FIG_DIR, "Fig_distance_measurements_CT_timeseries.pdf"),
  p_CT, width = 8, height = 6, units = "in"
)
ggsave(
  file.path(FIG_DIR, "Fig_distance_measurements_CT_timeseries.png"),
  p_CT, width = 8, height = 6, units = "in", dpi = 600
)

write.csv(
  dist_CT_wide,
  file.path(TAB_DIR, "distance_measurements_CT_timeseries_individuals.csv"),
  row.names = FALSE
)

write.csv(
  summary_CT,
  file.path(TAB_DIR, "distance_measurements_CT_timeseries_summary.csv"),
  row.names = FALSE
)

# ============================================================
# 3) Statistical tests (ANOVA + pairwise) and saving outputs
# ============================================================

STAT_DIR <- file.path("Outputs", "distance_measurements", "stats")
if (!dir.exists(STAT_DIR)) dir.create(STAT_DIR, recursive = TRUE)

run_stats <- function(df_long, group_var, dataset_name) {
  
  measurements <- unique(df_long$measurement)
  
  all_results <- list()
  all_pairwise <- list()
  
  for (m in measurements) {
    
    df_sub <- df_long %>%
      dplyr::filter(measurement == m)
    
    # ---------------------------
    # ANOVA
    # ---------------------------
    formula_str <- paste("value ~", group_var)
    fit <- aov(as.formula(formula_str), data = df_sub)
    
    anova_table <- as.data.frame(summary(fit)[[1]])
    anova_table$measurement <- m
    
    # ---------------------------
    # Pairwise t-tests (BH corrected)
    # ---------------------------
    pw <- pairwise.t.test(
      df_sub$value,
      df_sub[[group_var]],
      p.adjust.method = "BH"
    )
    
    pw_df <- as.data.frame(as.table(pw$p.value))
    colnames(pw_df) <- c("group1", "group2", "p_value")
    pw_df <- pw_df %>%
      dplyr::filter(!is.na(p_value)) %>%
      dplyr::mutate(measurement = m)
    
    # store
    all_results[[m]] <- anova_table
    all_pairwise[[m]] <- pw_df
    
    # ---------------------------
    # Save model summary text
    # ---------------------------
    sink(file.path(STAT_DIR, paste0(dataset_name, "_", gsub(" ", "_", m), "_anova.txt")))
    cat("Measurement:", m, "\n")
    print(summary(fit))
    sink()
  }
  
  # Combine + save
  anova_out <- dplyr::bind_rows(all_results)
  pairwise_out <- dplyr::bind_rows(all_pairwise)
  
  write.csv(
    anova_out,
    file.path(STAT_DIR, paste0(dataset_name, "_anova_summary.csv")),
    row.names = FALSE
  )
  
  write.csv(
    pairwise_out,
    file.path(STAT_DIR, paste0(dataset_name, "_pairwise_tests.csv")),
    row.names = FALSE
  )
}

# ---------------------------
# Run for both datasets
# ---------------------------

run_stats(dist_1950_long, "group", "1950_habitats")
run_stats(dist_CT_long, "group", "CT_timeseries")

# ============================================================
# Console output
# ============================================================

cat("\nSaved figures:\n")
cat("  ", normalizePath(file.path(FIG_DIR, "Fig_distance_measurements_1950_habitats.pdf")), "\n")
cat("  ", normalizePath(file.path(FIG_DIR, "Fig_distance_measurements_1950_habitats.png")), "\n")
cat("  ", normalizePath(file.path(FIG_DIR, "Fig_distance_measurements_CT_timeseries.pdf")), "\n")
cat("  ", normalizePath(file.path(FIG_DIR, "Fig_distance_measurements_CT_timeseries.png")), "\n")

cat("\nSaved tables:\n")
cat("  ", normalizePath(file.path(TAB_DIR, "distance_measurements_1950_habitats_individuals.csv")), "\n")
cat("  ", normalizePath(file.path(TAB_DIR, "distance_measurements_1950_habitats_summary.csv")), "\n")
cat("  ", normalizePath(file.path(TAB_DIR, "distance_measurements_CT_timeseries_individuals.csv")), "\n")
cat("  ", normalizePath(file.path(TAB_DIR, "distance_measurements_CT_timeseries_summary.csv")), "\n")