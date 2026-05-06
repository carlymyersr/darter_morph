# ============================================================
# Scripts/mouth_angle_3_point_analysis.R
#
# Compute and analyze a 3-point mouth-to-body angle using:
#   mouth_tip
#   mouth_base
#   body_axis_anterior
#
# Angle definition:
#   angle at mouth_base between vectors:
#     mouth_base -> mouth_tip
#     mouth_base -> body_axis_anterior
#
# Runs three group comparisons:
#   1) 1950 waterbodies:
#      CT 1950, Quabbin 1950, Swift 1950, Fort 1950, Sawmill 1950
#   2) CT time series:
#      CT 1950, CT 1956, CT 1970
#   3) Full landscape:
#      CT 1950, CT 1956, CT 1970, Quabbin 1950, Swift 1950, Fort 1950, Sawmill 1950
#
# Outputs are saved to:
#   outputs_mouth_to_body_angle/
# ============================================================

# ---------------------------
# Root detection
# ---------------------------
get_script_dir <- function() {
  this <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NULL)
  if (is.null(this) || is.na(this)) return(NULL)
  dirname(this)
}

script_dir <- get_script_dir()
if (!is.null(script_dir)) {
  project_root <- normalizePath(file.path(script_dir, ".."))
  if (dir.exists(file.path(project_root, "trait_measurements"))) {
    if (getwd() != project_root) setwd(project_root)
    cat("Project root set to:", project_root, "\n")
  }
}

# ---------------------------
# Libraries
# ---------------------------
suppressPackageStartupMessages({
  library(StereoMorph)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(multcompView)
  library(igraph)
  library(ggraph)
})

set.seed(123)

# ---------------------------
# User settings
# ---------------------------
ANGLE_SIGN_MULTIPLIER <- 1
# After inspecting extremes, switch to -1 if the biological direction is reversed.

alpha <- 0.05
trait_col <- "mouth_to_body_angle_deg"

# ---------------------------
# Inputs / outputs
# ---------------------------
trait_shapes_dir <- file.path("trait_measurements", "mouth_angle_shapes")
out_root <- "outputs_mouth_to_body_angle"
fig_dir <- file.path(out_root, "figures")
tables_dir <- file.path(out_root, "tables")
text_dir <- file.path(out_root, "text")
rds_dir <- file.path(out_root, "rds")

for (d in c(out_root, fig_dir, tables_dir, text_dir, rds_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

raw_out <- file.path(tables_dir, "mouth_to_body_angle_raw.csv")
merged_out <- file.path(tables_dir, "mouth_to_body_angle_with_metadata.csv")
extremes_out <- file.path(tables_dir, "mouth_to_body_angle_extremes.csv")
missing_out <- file.path(tables_dir, "mouth_to_body_angle_missing_landmark_check.csv")

if (!dir.exists(trait_shapes_dir)) {
  stop("Trait shapes folder not found: ", trait_shapes_dir,
       "\nRun R/10_digitize_trait_landmarks.R first, or update trait_shapes_dir.")
}

# ---------------------------
# Load metadata
# ---------------------------
source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")

# ---------------------------
# Read trait landmarks
# ---------------------------
trait_shapes <- StereoMorph::readShapes(trait_shapes_dir)

if (is.null(trait_shapes$landmarks.scaled)) {
  stop("No scaled landmarks found in: ", trait_shapes_dir)
}

coords <- trait_shapes$landmarks.scaled
ids <- dimnames(coords)[[3]]

required_lms <- c("mouth_tip", "mouth_base", "body_axis_anterior")
missing_lms <- setdiff(required_lms, dimnames(coords)[[1]])

if (length(missing_lms) > 0) {
  stop("Missing required trait landmarks: ", paste(missing_lms, collapse = ", "))
}

# ---------------------------
# Helper functions
# ---------------------------
get_point <- function(coords_arr, landmark, specimen) {
  as.numeric(coords_arr[landmark, , specimen])
}

signed_angle_degrees <- function(reference_vec, target_vec) {
  if (any(is.na(reference_vec)) || any(is.na(target_vec))) return(NA_real_)

  ref_norm <- sqrt(sum(reference_vec^2))
  tar_norm <- sqrt(sum(target_vec^2))

  if (is.na(ref_norm) || is.na(tar_norm) || ref_norm == 0 || tar_norm == 0) {
    return(NA_real_)
  }

  dot <- sum(reference_vec * target_vec)
  cross <- reference_vec[1] * target_vec[2] - reference_vec[2] * target_vec[1]

  atan2(cross, dot) * 180 / pi
}

angle_magnitude_degrees <- function(reference_vec, target_vec) {
  if (any(is.na(reference_vec)) || any(is.na(target_vec))) return(NA_real_)

  denom <- sqrt(sum(reference_vec^2)) * sqrt(sum(target_vec^2))
  if (is.na(denom) || denom == 0) return(NA_real_)

  cos_theta <- sum(reference_vec * target_vec) / denom
  cos_theta <- max(min(cos_theta, 1), -1)

  acos(cos_theta) * 180 / pi
}

# ---------------------------
# Missing landmark diagnostics
# ---------------------------
missing_trait_check <- purrr::map_dfr(ids, function(id) {
  vals <- coords[required_lms, , id, drop = FALSE]
  missing_by_lm <- apply(vals, 1, function(z) any(is.na(z)))

  tibble(
    specimen = id,
    has_missing_trait_landmarks = any(missing_by_lm),
    missing_landmarks = paste(required_lms[missing_by_lm], collapse = "; ")
  )
})

write.csv(missing_trait_check, missing_out, row.names = FALSE)

# ---------------------------
# Compute 3-point mouth-to-body angle
# ---------------------------
angle_df <- purrr::map_dfr(ids, function(id) {
  mouth_tip <- get_point(coords, "mouth_tip", id)
  mouth_base <- get_point(coords, "mouth_base", id)
  body_axis_anterior <- get_point(coords, "body_axis_anterior", id)

  # Vertex is mouth_base.
  # Reference vector points from mouth_base toward body_axis_anterior.
  # Target vector points from mouth_base toward mouth_tip.
  reference_vec <- body_axis_anterior - mouth_base
  target_vec <- mouth_tip - mouth_base

  raw_signed_angle <- signed_angle_degrees(reference_vec, target_vec)
  signed_angle <- ANGLE_SIGN_MULTIPLIER * raw_signed_angle
  unsigned_angle <- angle_magnitude_degrees(reference_vec, target_vec)

  tibble(
    specimen = id,
    raw_mouth_to_body_angle_deg = raw_signed_angle,
    mouth_to_body_angle_deg = signed_angle,
    mouth_to_body_angle_abs_deg = unsigned_angle,
    angle_sign_multiplier = ANGLE_SIGN_MULTIPLIER
  )
})

write.csv(angle_df, raw_out, row.names = FALSE)

mouth_df <- gdf %>%
  left_join(angle_df, by = "specimen") %>%
  mutate(year = as.integer(year))

write.csv(mouth_df, merged_out, row.names = FALSE)

extremes_df <- bind_rows(
  mouth_df %>%
    filter(!is.na(.data[[trait_col]])) %>%
    arrange(.data[[trait_col]]) %>%
    slice_head(n = 10) %>%
    mutate(extreme_type = "most_negative"),
  mouth_df %>%
    filter(!is.na(.data[[trait_col]])) %>%
    arrange(desc(.data[[trait_col]])) %>%
    slice_head(n = 10) %>%
    mutate(extreme_type = "most_positive")
) %>%
  select(
    extreme_type, specimen, habitat, year,
    mouth_to_body_angle_deg,
    raw_mouth_to_body_angle_deg,
    mouth_to_body_angle_abs_deg,
    angle_sign_multiplier
  )

write.csv(extremes_df, extremes_out, row.names = FALSE)

# ---------------------------
# Palettes
# ---------------------------
group_palette_full <- c(
  "CT_1950" = "steelblue",
  "CT_1956" = "deepskyblue3",
  "CT_1970" = "navy",
  "Fort_1950" = "black",
  "Quabbin_1950" = "darkgoldenrod2",
  "Sawmill_1950" = "orchid3",
  "Swift_1950" = "tomato"
)

# ---------------------------
# Analysis helpers
# ---------------------------
summarize_trait <- function(df, group_col, trait_col) {
  group_col <- rlang::ensym(group_col)
  trait_col <- rlang::ensym(trait_col)

  df %>%
    group_by(!!group_col) %>%
    summarize(
      n = sum(!is.na(!!trait_col)),
      mean = mean(!!trait_col, na.rm = TRUE),
      sd = sd(!!trait_col, na.rm = TRUE),
      se = sd / sqrt(n),
      median = median(!!trait_col, na.rm = TRUE),
      min = min(!!trait_col, na.rm = TRUE),
      max = max(!!trait_col, na.rm = TRUE),
      q25 = quantile(!!trait_col, 0.25, na.rm = TRUE),
      q75 = quantile(!!trait_col, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

pairwise_wilcox_bh <- function(df, group_col, trait_col) {
  group_col_name <- rlang::as_name(rlang::ensym(group_col))
  trait_col_name <- rlang::as_name(rlang::ensym(trait_col))

  df2 <- df %>%
    filter(!is.na(.data[[group_col_name]]), !is.na(.data[[trait_col_name]])) %>%
    mutate(.group = droplevels(factor(.data[[group_col_name]])))

  levs <- levels(df2$.group)
  if (length(levs) < 2) {
    return(tibble(group1 = character(), group2 = character(), p_adj = numeric()))
  }

  pw <- pairwise.wilcox.test(
    x = df2[[trait_col_name]],
    g = df2$.group,
    p.adjust.method = "BH",
    exact = FALSE
  )

  pairwise_df <- as.data.frame(as.table(pw$p.value), stringsAsFactors = FALSE) %>%
    filter(!is.na(Freq)) %>%
    transmute(
      group1 = as.character(Var1),
      group2 = as.character(Var2),
      p_adj = as.numeric(Freq)
    ) %>%
    arrange(p_adj)

  summary_df <- df2 %>%
    group_by(.group) %>%
    summarize(
      mean = mean(.data[[trait_col_name]], na.rm = TRUE),
      median = median(.data[[trait_col_name]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(group = .group)

  pairwise_df %>%
    left_join(summary_df %>% select(group1 = group, mean1 = mean, median1 = median), by = "group1") %>%
    left_join(summary_df %>% select(group2 = group, mean2 = mean, median2 = median), by = "group2") %>%
    mutate(
      mean_diff = mean1 - mean2,
      median_diff = median1 - median2,
      significant_BH_0.05 = p_adj < 0.05
    )
}

make_letters_from_pairwise <- function(pairwise_df, group_levels, p_col = "p_adj", alpha = 0.05) {
  group_levels <- as.character(group_levels)

  if (nrow(pairwise_df) == 0) {
    return(tibble(group = factor(group_levels, levels = group_levels), letters = "a"))
  }

  pvals <- pairwise_df[[p_col]]
  names(pvals) <- paste(pairwise_df$group1, pairwise_df$group2, sep = "-")

  letters <- multcompView::multcompLetters(pvals, threshold = alpha)$Letters
  out <- tibble(group = names(letters), letters = as.character(letters))

  missing_groups <- setdiff(group_levels, out$group)
  if (length(missing_groups) > 0) {
    out <- bind_rows(out, tibble(group = missing_groups, letters = "a"))
  }

  out %>%
    mutate(group = factor(group, levels = group_levels)) %>%
    arrange(group)
}

make_letter_positions <- function(df, group_col, trait_col, pad_frac = 0.08) {
  group_col <- rlang::ensym(group_col)
  trait_col <- rlang::ensym(trait_col)

  y_range <- range(dplyr::pull(df, !!trait_col), na.rm = TRUE)
  pad <- diff(y_range) * pad_frac
  if (!is.finite(pad) || pad == 0) pad <- 1

  df %>%
    group_by(!!group_col) %>%
    summarize(y_pos = max(!!trait_col, na.rm = TRUE) + pad, .groups = "drop") %>%
    rename(group = !!group_col)
}

make_boxplot_with_letters <- function(df, group_col, trait_col, letters_df, palette,
                                      title, subtitle, outfile_base) {
  group_col <- rlang::ensym(group_col)
  trait_col <- rlang::ensym(trait_col)
  group_col_name <- rlang::as_name(group_col)
  trait_col_name <- rlang::as_name(trait_col)

  plot_df <- df %>%
    filter(!is.na(!!group_col), !is.na(!!trait_col)) %>%
    mutate(!!group_col_name := droplevels(factor(.data[[group_col_name]])))

  group_levels <- levels(plot_df[[group_col_name]])

  letter_pos <- make_letter_positions(plot_df, !!group_col, !!trait_col) %>%
    left_join(letters_df, by = "group")

  y_range <- range(plot_df[[trait_col_name]], na.rm = TRUE)
  y_pad <- diff(y_range) * 0.22
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 1

  p <- ggplot(plot_df, aes(x = !!group_col, y = !!trait_col, fill = !!group_col)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, alpha = 0.8) +
    geom_boxplot(width = 0.72, outlier.shape = NA, alpha = 0.9, linewidth = 0.35) +
    geom_jitter(width = 0.12, size = 1.7, alpha = 0.75) +
    geom_text(
      data = letter_pos,
      aes(x = group, y = y_pos, label = letters),
      inherit.aes = FALSE,
      size = 4
    ) +
    scale_fill_manual(values = palette[group_levels], drop = FALSE) +
    coord_cartesian(ylim = c(y_range[1] - diff(y_range) * 0.05, y_range[2] + y_pad)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = "Mouth-to-body angle (degrees)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.title = element_text(face = "bold")
    )

  ggsave(file.path(fig_dir, paste0(outfile_base, ".png")), p, width = 7.2, height = 5, dpi = 400, bg = "white")
  ggsave(file.path(fig_dir, paste0(outfile_base, ".pdf")), p, width = 7.2, height = 5, bg = "white")

  p
}

make_nonsig_network_plot <- function(pairwise_df, group_levels, palette, title, outfile_base,
                                     alpha = 0.05, layout_type = "fr") {
  group_levels <- as.character(group_levels)

  edges <- pairwise_df %>%
    filter(!is.na(p_adj), p_adj >= alpha) %>%
    transmute(from = group1, to = group2)

  nodes <- tibble(name = group_levels)

  graph <- igraph::graph_from_data_frame(
    d = edges,
    vertices = nodes,
    directed = FALSE
  )

  p <- ggraph(graph, layout = layout_type) +
    geom_edge_link(linewidth = 0.55, alpha = 0.60, color = "gray50") +
    geom_node_point(aes(fill = name), shape = 21, size = 6.5, color = "black", stroke = 0.7) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3.2, color = "black") +
    scale_fill_manual(values = palette[group_levels], drop = FALSE) +
    theme_void(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 11),
      plot.caption = element_text(hjust = 0.5, size = 8),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    labs(
      title = title,
      caption = "Edges connect groups that are not significantly different (BH-adjusted p >= 0.05)"
    )

  ggsave(file.path(fig_dir, paste0(outfile_base, ".png")), p, width = 4, height = 4, dpi = 400, bg = "white")
  ggsave(file.path(fig_dir, paste0(outfile_base, ".pdf")), p, width = 4, height = 4, bg = "white")

  p
}

write_stats_txt <- function(df, group_col, trait_col, summary_df, pairwise_df, letters_df,
                            analysis_name, outfile) {
  group_col_name <- rlang::as_name(rlang::ensym(group_col))
  trait_col_name <- rlang::as_name(rlang::ensym(trait_col))

  df2 <- df %>%
    filter(!is.na(.data[[group_col_name]]), !is.na(.data[[trait_col_name]])) %>%
    mutate(.group = droplevels(factor(.data[[group_col_name]])))

  dat <- data.frame(y = df2[[trait_col_name]], group = df2$.group)
  fit <- lm(y ~ group, data = dat)

  sink(outfile)
  cat("============================================================\n")
  cat(analysis_name, "\n")
  cat("============================================================\n\n")
  cat("Trait: ", trait_col_name, "\n", sep = "")
  cat("Angle definition: 3-point angle at mouth_base between mouth_base->body_axis_anterior and mouth_base->mouth_tip.\n\n")
  cat("Grouping variable: ", group_col_name, "\n", sep = "")
  cat("N after removing missing trait values: ", nrow(df2), "\n\n", sep = "")

  cat("Group counts:\n")
  print(table(df2$.group))

  cat("\nSummary statistics:\n")
  print(summary_df)

  cat("\nOverall one-way ANOVA:\n")
  print(anova(fit))

  cat("\nLinear model summary:\n")
  print(summary(fit))

  cat("\nKruskal-Wallis test:\n")
  print(kruskal.test(y ~ group, data = dat))

  cat("\nPairwise Wilcoxon rank-sum tests, BH-adjusted:\n")
  print(pairwise_df)

  cat("\nCompact letter display from BH-adjusted pairwise Wilcoxon tests:\n")
  print(letters_df)

  cat("\nInterpretation of letters:\n")
  cat("Groups sharing a letter are not significantly different at alpha = 0.05.\n")

  cat("\nNon-significant network interpretation:\n")
  cat("Edges connect groups that are not significantly different at alpha = 0.05 after BH adjustment.\n")
  sink()
}

run_mouth_to_body_angle_analysis <- function(df, group_col, analysis_label, outfile_prefix, palette) {
  group_col_name <- rlang::as_name(rlang::ensym(group_col))

  df_use <- df %>%
    filter(!is.na(.data[[group_col_name]]), !is.na(.data[[trait_col]])) %>%
    mutate(!!group_col_name := droplevels(factor(.data[[group_col_name]])))

  group_levels <- levels(df_use[[group_col_name]])

  summary_df <- summarize_trait(df_use, !!rlang::sym(group_col_name), !!rlang::sym(trait_col))
  pairwise_df <- pairwise_wilcox_bh(df_use, !!rlang::sym(group_col_name), !!rlang::sym(trait_col))
  letters_df <- make_letters_from_pairwise(pairwise_df, group_levels, "p_adj", alpha)

  write.csv(summary_df, file.path(tables_dir, paste0(outfile_prefix, "_summary.csv")), row.names = FALSE)
  write.csv(pairwise_df, file.path(tables_dir, paste0(outfile_prefix, "_pairwise_wilcox_BH.csv")), row.names = FALSE)
  write.csv(letters_df, file.path(tables_dir, paste0(outfile_prefix, "_letters_BH.csv")), row.names = FALSE)

  write_stats_txt(
    df = df_use,
    group_col = !!rlang::sym(group_col_name),
    trait_col = !!rlang::sym(trait_col),
    summary_df = summary_df,
    pairwise_df = pairwise_df,
    letters_df = letters_df,
    analysis_name = analysis_label,
    outfile = file.path(text_dir, paste0(outfile_prefix, "_statistics.txt"))
  )

  p_box <- make_boxplot_with_letters(
    df = df_use,
    group_col = !!rlang::sym(group_col_name),
    trait_col = !!rlang::sym(trait_col),
    letters_df = letters_df,
    palette = palette,
    title = analysis_label,
    subtitle = "Letters show BH-adjusted pairwise Wilcoxon groupings; dashed line = zero",
    outfile_base = paste0(outfile_prefix, "_boxplot_letters")
  )

  p_net <- make_nonsig_network_plot(
    pairwise_df = pairwise_df,
    group_levels = group_levels,
    palette = palette,
    title = paste0("Non-significant network: ", analysis_label),
    outfile_base = paste0(outfile_prefix, "_network_nonsignificant"),
    alpha = alpha
  )

  list(
    data = df_use,
    summary = summary_df,
    pairwise = pairwise_df,
    letters = letters_df,
    boxplot = p_box,
    network = p_net
  )
}

# ---------------------------
# Build analysis datasets
# ---------------------------
hab_1950_df <- mouth_df %>%
  filter(
    !is.na(habitat),
    year == 1950,
    habitat %in% c("Connecticut River", "Quabbin", "Swift River", "Fort River", "Sawmill River")
  ) %>%
  mutate(
    mouth_group_1950 = case_when(
      habitat == "Connecticut River" ~ "CT_1950",
      habitat == "Quabbin" ~ "Quabbin_1950",
      habitat == "Swift River" ~ "Swift_1950",
      habitat == "Fort River" ~ "Fort_1950",
      habitat == "Sawmill River" ~ "Sawmill_1950",
      TRUE ~ NA_character_
    ),
    mouth_group_1950 = factor(
      mouth_group_1950,
      levels = c("CT_1950", "Quabbin_1950", "Swift_1950", "Fort_1950", "Sawmill_1950")
    )
  )

ct_time_df <- mouth_df %>%
  filter(
    habitat == "Connecticut River",
    year %in% c(1950, 1956, 1970)
  ) %>%
  mutate(
    mouth_group_ct_time = case_when(
      year == 1950 ~ "CT_1950",
      year == 1956 ~ "CT_1956",
      year == 1970 ~ "CT_1970",
      TRUE ~ NA_character_
    ),
    mouth_group_ct_time = factor(
      mouth_group_ct_time,
      levels = c("CT_1950", "CT_1956", "CT_1970")
    )
  )

full_landscape_df <- mouth_df %>%
  filter(
    !is.na(habitat),
    (
      (habitat == "Connecticut River" & year %in% c(1950, 1956, 1970)) |
        (habitat != "Connecticut River" & year == 1950)
    ),
    habitat %in% c("Connecticut River", "Quabbin", "Swift River", "Fort River", "Sawmill River")
  ) %>%
  mutate(
    mouth_group_full = case_when(
      habitat == "Connecticut River" & year == 1950 ~ "CT_1950",
      habitat == "Connecticut River" & year == 1956 ~ "CT_1956",
      habitat == "Connecticut River" & year == 1970 ~ "CT_1970",
      habitat == "Quabbin" & year == 1950 ~ "Quabbin_1950",
      habitat == "Swift River" & year == 1950 ~ "Swift_1950",
      habitat == "Fort River" & year == 1950 ~ "Fort_1950",
      habitat == "Sawmill River" & year == 1950 ~ "Sawmill_1950",
      TRUE ~ NA_character_
    ),
    mouth_group_full = factor(
      mouth_group_full,
      levels = c("CT_1970", "CT_1956", "CT_1950", "Quabbin_1950", "Swift_1950", "Fort_1950", "Sawmill_1950")
    )
  ) %>%
  filter(!is.na(mouth_group_full))

cat("\n1950 habitat counts:\n")
print(table(hab_1950_df$mouth_group_1950, useNA = "ifany"))

cat("\nCT time series counts:\n")
print(table(ct_time_df$mouth_group_ct_time, useNA = "ifany"))

cat("\nFull landscape counts:\n")
print(table(full_landscape_df$mouth_group_full, useNA = "ifany"))

# ---------------------------
# Run analyses
# ---------------------------
res_1950 <- run_mouth_to_body_angle_analysis(
  df = hab_1950_df,
  group_col = mouth_group_1950,
  analysis_label = "Mouth-to-body angle across 1950 waterbodies",
  outfile_prefix = "1950_waterbodies_mouth_to_body_angle",
  palette = group_palette_full
)

res_ct_time <- run_mouth_to_body_angle_analysis(
  df = ct_time_df,
  group_col = mouth_group_ct_time,
  analysis_label = "Mouth-to-body angle across Connecticut River time series",
  outfile_prefix = "CT_timeseries_mouth_to_body_angle",
  palette = group_palette_full
)

res_full <- run_mouth_to_body_angle_analysis(
  df = full_landscape_df,
  group_col = mouth_group_full,
  analysis_label = "Mouth-to-body angle across CT time series and 1950 waterbodies",
  outfile_prefix = "full_landscape_mouth_to_body_angle",
  palette = group_palette_full
)

saveRDS(
  list(
    angle_data = mouth_df,
    res_1950 = res_1950,
    res_ct_time = res_ct_time,
    res_full = res_full
  ),
  file.path(rds_dir, "mouth_to_body_angle_analysis_results.rds")
)


# ============================================================
# SICB-style faceted mouth angle figure
# Add after res_full is created
# ============================================================

facet_fig_dir <- file.path(out_root, "figures", "SICB_facet_style")
dir.create(facet_fig_dir, recursive = TRUE, showWarnings = FALSE)

facet_levels <- c(
  "CT_1970", "CT_1956", "CT_1950",
  "Quabbin_1950", "Swift_1950", "Fort_1950", "Sawmill_1950"
)

facet_labels <- c(
  "CT_1970" = "CT 1970",
  "CT_1956" = "CT 1956",
  "CT_1950" = "CT 1950",
  "Quabbin_1950" = "Quabbin",
  "Swift_1950" = "Swift",
  "Fort_1950" = "Fort",
  "Sawmill_1950" = "Sawmill"
)

mouth_facet_dat <- full_landscape_df %>%
  filter(
    !is.na(mouth_group_full),
    is.finite(.data[[trait_col]])
  ) %>%
  mutate(
    group = factor(as.character(mouth_group_full), levels = facet_levels),
    group_label = factor(facet_labels[as.character(group)],
                         levels = facet_labels[facet_levels]),
    facet_label = factor(
      "Upper mouth angle",
      levels = c("Upper mouth angle", " ")
    )
  )

# Empty second facet to preserve two-column SICB layout
empty_facet_dat <- mouth_facet_dat %>%
  slice(0) %>%
  mutate(
    facet_label = factor(" ", levels = c("Upper mouth angle", " "))
  )

mouth_facet_dat2 <- bind_rows(mouth_facet_dat, empty_facet_dat)

p_mouth_facet <- ggplot(
  mouth_facet_dat2,
  aes(x = group_label, y = .data[[trait_col]])
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
  facet_wrap(
    ~ facet_label,
    scales = "free_y",
    ncol = 2,
    drop = FALSE
  ) +
  labs(
    x = NULL,
    y = "Upper mouth angle (degrees)"
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
  filename = file.path(facet_fig_dir, "upper_mouth_angle_SICB_facet.pdf"),
  plot = p_mouth_facet,
  width = 4,
  height = 2.6,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = file.path(facet_fig_dir, "upper_mouth_angle_SICB_facet.png"),
  plot = p_mouth_facet,
  width = 4,
  height = 2.6,
  units = "in",
  dpi = 400,
  bg = "white"
)

cat("\nSaved SICB-style faceted upper mouth angle figure to:\n")
cat(normalizePath(facet_fig_dir), "\n")





# ---------------------------
# Manifest and run info
# ---------------------------
manifest <- tibble::tribble(
  ~category, ~file,
  "table", raw_out,
  "table", merged_out,
  "table", extremes_out,
  "table", missing_out,
  "text", file.path(text_dir, "1950_waterbodies_mouth_to_body_angle_statistics.txt"),
  "text", file.path(text_dir, "CT_timeseries_mouth_to_body_angle_statistics.txt"),
  "text", file.path(text_dir, "full_landscape_mouth_to_body_angle_statistics.txt"),
  "figure", file.path(fig_dir, "1950_waterbodies_mouth_to_body_angle_boxplot_letters.png"),
  "figure", file.path(fig_dir, "CT_timeseries_mouth_to_body_angle_boxplot_letters.png"),
  "figure", file.path(fig_dir, "full_landscape_mouth_to_body_angle_boxplot_letters.png"),
  "figure", file.path(fig_dir, "1950_waterbodies_mouth_to_body_angle_network_nonsignificant.png"),
  "figure", file.path(fig_dir, "CT_timeseries_mouth_to_body_angle_network_nonsignificant.png"),
  "figure", file.path(fig_dir, "full_landscape_mouth_to_body_angle_network_nonsignificant.png"),
  "rds", file.path(rds_dir, "mouth_to_body_angle_analysis_results.rds")
)

write.csv(manifest, file.path(out_root, "MANIFEST.csv"), row.names = FALSE)

sink(file.path(out_root, "RUN_INFO.txt"))
cat("Working directory:", getwd(), "\n")
cat("Date:", as.character(Sys.time()), "\n")
cat("Trait shapes directory:", trait_shapes_dir, "\n")
cat("Angle sign multiplier:", ANGLE_SIGN_MULTIPLIER, "\n")
cat("Angle definition: 3-point angle at mouth_base between mouth_base->body_axis_anterior and mouth_base->mouth_tip.\n\n")

cat("Missing mouth-to-body angle values:\n")
print(table(is.na(mouth_df[[trait_col]]), useNA = "ifany"))

cat("\n1950 habitat counts:\n")
print(table(hab_1950_df$mouth_group_1950, useNA = "ifany"))

cat("\nCT time series counts:\n")
print(table(ct_time_df$mouth_group_ct_time, useNA = "ifany"))

cat("\nFull landscape counts:\n")
print(table(full_landscape_df$mouth_group_full, useNA = "ifany"))

cat("\nSession info:\n")
print(sessionInfo())
sink()

cat("\n============================================================\n")
cat("Done.\n")
cat("Outputs written to:\n")
cat("  ", out_root, "\n")
cat("============================================================\n")
