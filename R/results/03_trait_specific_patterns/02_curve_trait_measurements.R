# ============================================================
# Scripts/curve_shape_metrics_raw_then_sizecorrected.R
#
# Goal:
#   Measure shape properties of the two semilandmark curves
#   from RAW landmark coordinates (coords_all), then remove size.
#
# Curves:
#   1) snout = cranium_orbital curve
#   2) hyoid = hyoid_pelvic curve
#
# Metrics:
#   - tortuosity = curve length / chord length
#   - max deviation from chord
#
# Outputs:
#   - PNG boxplots for 1950 habitat comparison
#   - PNG boxplots for Connecticut River time series
#   - CSV summaries by group
#   - CSV ANOVA tables
#   - CSV pairwise Tukey comparisons
#
# Notes:
#   - Metrics are extracted from RAW coords_all (pre-GPA)
#   - Size correction is applied afterward using residuals from:
#         metric ~ logCsize
#   - Residuals are mean-centered around 0, so plots show
#     size-corrected centered values
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(purrr)
})

# ------------------------------------------------------------
# 0) Load required objects if needed
# ------------------------------------------------------------
if (!exists("coords_all", inherits = TRUE)) source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
if (!exists("gdf", inherits = TRUE))        source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")

# ------------------------------------------------------------
# 1) Output directory
# ------------------------------------------------------------
out_dir <- file.path("Figures", "curve_shape_metrics")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("\nSaving outputs to:\n", normalizePath(out_dir), "\n")

# ------------------------------------------------------------
# 2) Define the two curves
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# 3) Sanity checks
# ------------------------------------------------------------
all_needed_pts <- unique(unlist(curve_defs))
missing_pts <- setdiff(all_needed_pts, dimnames(coords_all)[[1]])

if (length(missing_pts) > 0) {
  stop("These curve landmarks are missing from coords_all:\n",
       paste(missing_pts, collapse = ", "))
}

if (!all(dimnames(coords_all)[[3]] %in% gdf$specimen)) {
  stop("Some coords_all specimen IDs are missing from gdf.")
}

# Align metadata to coords_all order
gdf_raw <- gdf %>%
  filter(specimen %in% dimnames(coords_all)[[3]]) %>%
  slice(match(dimnames(coords_all)[[3]], specimen))

stopifnot(identical(gdf_raw$specimen, dimnames(coords_all)[[3]]))

# ------------------------------------------------------------
# 4) Geometry helper functions
# ------------------------------------------------------------

# Total polyline length
curve_length <- function(mat) {
  seg <- mat[-1, , drop = FALSE] - mat[-nrow(mat), , drop = FALSE]
  sum(sqrt(rowSums(seg^2)))
}

# Straight-line distance between endpoints
chord_length <- function(mat) {
  sqrt(sum((mat[nrow(mat), ] - mat[1, ])^2))
}

# Tortuosity = curve length / chord length
curve_tortuosity <- function(mat) {
  cl <- chord_length(mat)
  if (isTRUE(all.equal(cl, 0))) return(NA_real_)
  curve_length(mat) / cl
}

# Maximum perpendicular deviation from the endpoint-to-endpoint chord
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

# Estimated mode for continuous data using kernel density peak
estimate_mode <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  d <- density(x)
  d$x[which.max(d$y)]
}

# ------------------------------------------------------------
# 5) Extract raw curve metrics specimen by specimen
# ------------------------------------------------------------
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

# ------------------------------------------------------------
# 6) Size-correct each metric AFTER extraction
#    value_raw ~ logCsize
# ------------------------------------------------------------
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
  ungroup()

curve_metrics_sc <- curve_metrics_sc %>%
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

# ------------------------------------------------------------
# 7) Ordered grouping variables
# ------------------------------------------------------------

# Habitat plot: 1950 habitats only
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
  mutate(
    group = factor(habitat, levels = habitat_levels)
  )

# CT time series plot: Connecticut River only, 1950/1956/1970
ct_year_levels <- c(1950, 1956, 1970)

ct_dat <- curve_metrics_sc %>%
  filter(
    habitat == "Connecticut River",
    year %in% ct_year_levels
  ) %>%
  mutate(
    group = factor(year, levels = ct_year_levels)
  )

# ------------------------------------------------------------
# 8) Summary tables
# ------------------------------------------------------------
make_summary_table <- function(dat) {
  dat %>%
    group_by(group, curve_label, metric_label) %>%
    summarise(
      n      = sum(is.finite(value_plot)),
      mean   = mean(value_plot, na.rm = TRUE),
      median = median(value_plot, na.rm = TRUE),
      mode   = estimate_mode(value_plot),
      min    = min(value_plot, na.rm = TRUE),
      max    = max(value_plot, na.rm = TRUE),
      .groups = "drop"
    )
}

habitat_summary <- make_summary_table(habitat_dat)
ct_summary <- make_summary_table(ct_dat)

write.csv(
  habitat_summary,
  file.path(out_dir, "curve_shape_metrics_habitat_1950_summary.csv"),
  row.names = FALSE
)

write.csv(
  ct_summary,
  file.path(out_dir, "curve_shape_metrics_CT_timeseries_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9) Plotting function
# ------------------------------------------------------------
make_curve_metric_plot <- function(dat, xlab, title_text) {
  ggplot(dat, aes(x = group, y = value_plot)) +
    geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dashed", color = "grey50") +
    geom_boxplot(
      width = 0.55,
      outlier.shape = NA,
      fill = "white",
      color = "black"
    ) +
    geom_jitter(
      width = 0.12,
      height = 0,
      alpha = 0.8,
      size = 2
    ) +
    facet_grid(metric_label ~ curve_label, scales = "free_y") +
    labs(
      x = xlab,
      y = "Size-corrected, mean-centered value",
      title = title_text
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95"),
      axis.text.x = element_text(angle = 25, hjust = 1)
    )
}

p_habitat <- make_curve_metric_plot(
  habitat_dat,
  xlab = "Habitat",
  title_text = "Curve shape metrics by habitat (1950 only)"
)

p_ct <- make_curve_metric_plot(
  ct_dat,
  xlab = "Year",
  title_text = "Curve shape metrics across Connecticut River time series"
)

ggsave(
  filename = file.path(out_dir, "curve_shape_metrics_habitat_1950.png"),
  plot = p_habitat,
  width = 10,
  height = 8,
  dpi = 400
)

ggsave(
  filename = file.path(out_dir, "curve_shape_metrics_CT_timeseries.png"),
  plot = p_ct,
  width = 8,
  height = 8,
  dpi = 400
)

# ------------------------------------------------------------
# 10) Statistical tests + pairwise comparisons
# ------------------------------------------------------------
run_stats <- function(dat, label) {
  results_list <- list()
  pairwise_list <- list()
  
  combos <- dat %>%
    distinct(curve_label, metric_label)
  
  cat("\n============================\n")
  cat("Running stats for:", label, "\n")
  cat("============================\n")
  
  for (i in seq_len(nrow(combos))) {
    curve_i  <- combos$curve_label[i]
    metric_i <- combos$metric_label[i]
    
    sub <- dat %>%
      filter(
        curve_label == curve_i,
        metric_label == metric_i
      ) %>%
      filter(is.finite(value_plot), !is.na(group)) %>%
      droplevels()
    
    cat("\n---", label, "|", curve_i, "|", metric_i, "---\n")
    print(table(sub$group))
    
    if (nrow(sub) < 3 || n_distinct(sub$group) < 2) {
      cat("Skipping: not enough observations or groups\n")
      next
    }
    
    fit_lm <- try(lm(value_plot ~ group, data = sub), silent = TRUE)
    if (inherits(fit_lm, "try-error")) {
      cat("lm() failed\n")
      print(fit_lm)
      next
    }
    
    aov_tab <- try(anova(fit_lm), silent = TRUE)
    if (inherits(aov_tab, "try-error")) {
      cat("anova() failed\n")
      print(aov_tab)
      next
    }
    
    aov_out <- data.frame(
      dataset = label,
      curve   = curve_i,
      metric  = metric_i,
      term    = rownames(aov_tab),
      aov_tab,
      row.names = NULL,
      check.names = FALSE
    )
    
    results_list[[length(results_list) + 1]] <- aov_out
    cat("ANOVA OK\n")
    
    fit_aov <- try(aov(value_plot ~ group, data = sub), silent = TRUE)
    if (inherits(fit_aov, "try-error")) {
      cat("aov() failed for Tukey\n")
      print(fit_aov)
      next
    }
    
    tuk <- try(TukeyHSD(fit_aov, "group"), silent = TRUE)
    if (inherits(tuk, "try-error")) {
      cat("TukeyHSD() failed\n")
      print(tuk)
      next
    }
    
    if (!("group" %in% names(tuk))) {
      cat("Tukey output had no 'group' component\n")
      next
    }
    
    tuk_df <- as.data.frame(tuk$group)
    tuk_df$comparison <- rownames(tuk_df)
    
    tuk_df <- tuk_df %>%
      mutate(
        dataset = label,
        curve   = curve_i,
        metric  = metric_i
      ) %>%
      select(dataset, curve, metric, comparison, everything())
    
    pairwise_list[[length(pairwise_list) + 1]] <- tuk_df
    cat("Tukey OK\n")
  }
  
  anova_df <- if (length(results_list) > 0) bind_rows(results_list) else data.frame()
  pairwise_df <- if (length(pairwise_list) > 0) bind_rows(pairwise_list) else data.frame()
  
  cat("\nFinished", label, "\n")
  cat("ANOVA rows:", nrow(anova_df), "\n")
  cat("Pairwise rows:", nrow(pairwise_df), "\n")
  
  list(
    anova = anova_df,
    pairwise = pairwise_df
  )
}

habitat_stats <- run_stats(habitat_dat, "habitat_1950")
ct_stats      <- run_stats(ct_dat, "ct_timeseries")

cat("\nSaving stats files...\n")

if (nrow(habitat_stats$anova) > 0) {
  write.csv(
    habitat_stats$anova,
    file.path(out_dir, "curve_metrics_habitat_1950_anova.csv"),
    row.names = FALSE
  )
  cat("Saved habitat ANOVA\n")
} else {
  cat("No habitat ANOVA rows to save\n")
}

if (nrow(habitat_stats$pairwise) > 0) {
  write.csv(
    habitat_stats$pairwise,
    file.path(out_dir, "curve_metrics_habitat_1950_pairwise.csv"),
    row.names = FALSE
  )
  cat("Saved habitat pairwise\n")
} else {
  cat("No habitat pairwise rows to save\n")
}

if (nrow(ct_stats$anova) > 0) {
  write.csv(
    ct_stats$anova,
    file.path(out_dir, "curve_metrics_CT_timeseries_anova.csv"),
    row.names = FALSE
  )
  cat("Saved CT ANOVA\n")
} else {
  cat("No CT ANOVA rows to save\n")
}

if (nrow(ct_stats$pairwise) > 0) {
  write.csv(
    ct_stats$pairwise,
    file.path(out_dir, "curve_metrics_CT_timeseries_pairwise.csv"),
    row.names = FALSE
  )
  cat("Saved CT pairwise\n")
} else {
  cat("No CT pairwise rows to save\n")
}

# ------------------------------------------------------------
# 11) Save full specimen-level table too
# ------------------------------------------------------------
write.csv(
  curve_metrics_sc,
  file.path(out_dir, "curve_shape_metrics_all_specimens_size_corrected.csv"),
  row.names = FALSE
)



# ------------------------------------------------------------
# 12) Final reporting
# ------------------------------------------------------------
cat("\nSaved files to:\n", normalizePath(out_dir), "\n")
cat("\nFiles created:\n")
cat("  - curve_shape_metrics_habitat_1950.png\n")
cat("  - curve_shape_metrics_CT_timeseries.png\n")
cat("  - curve_shape_metrics_habitat_1950_summary.csv\n")
cat("  - curve_shape_metrics_CT_timeseries_summary.csv\n")
cat("  - curve_metrics_habitat_1950_anova.csv\n")
cat("  - curve_metrics_habitat_1950_pairwise.csv\n")
cat("  - curve_metrics_CT_timeseries_anova.csv\n")
cat("  - curve_metrics_CT_timeseries_pairwise.csv\n")
cat("  - curve_shape_metrics_all_specimens_size_corrected.csv\n")




