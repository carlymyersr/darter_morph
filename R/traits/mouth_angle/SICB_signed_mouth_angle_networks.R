# ============================================================
# Scripts/SICB_signed_mouth_angle_networks.R
#
# SICB-style graph-only NON-SIGNIFICANT pairwise networks
# for signed mouth angle.
#
# Matches the network style used in:
#   SICB_curves_distances_figures.R
#
# Edges connect groups that are NOT significantly different:
#   BH-adjusted p >= 0.05
#
# Input:
#   trait_measurements/mouth_angle_signed_with_metadata.csv
#
# Before running:
#   source("R/build_trait_measurements_signed.R")
#
# Outputs:
#   Figures/SICB_signed_mouth_angle_networks/
#   Outputs/SICB_signed_mouth_angle_stats/
#
# Analyses:
#   1) CT time series: CT 1970, CT 1956, CT 1950
#   2) 1950 habitats: CT 1950, Quabbin, Swift, Fort, Sawmill
#   3) Combined: CT 1970, CT 1956, CT 1950, Quabbin, Swift, Fort, Sawmill
#   4) CT pooled reference: CT time series, Quabbin, Swift, Fort, Sawmill
#
# Interpretation:
#   signed_mouth_angle_deg should be oriented so that:
#     positive = more upturned
#     negative = more downturned
#   If reversed, change ANGLE_SIGN_MULTIPLIER in
#   R/build_trait_measurements_signed.R and rerun first.
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
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(multcompView)
  library(stringr)
  library(tibble)
})

# ---------------------------
# Inputs / outputs
# ---------------------------
INFILE <- file.path("trait_measurements", "mouth_angle_signed_with_metadata.csv")

if (!file.exists(INFILE)) {
  stop(
    "Could not find: ", INFILE,
    "\nRun source('R/build_trait_measurements_signed.R') first."
  )
}

OUTDIR <- file.path("Figures", "SICB_signed_mouth_angle_networks")
STAT_OUTDIR <- file.path("Outputs", "SICB_signed_mouth_angle_stats")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STAT_OUTDIR, recursive = TRUE, showWarnings = FALSE)

alpha <- 0.05
trait_col <- "signed_mouth_angle_deg"

mouth_df <- read.csv(INFILE, stringsAsFactors = FALSE)

required_cols <- c("specimen", "habitat", "year", trait_col)
missing_cols <- setdiff(required_cols, names(mouth_df))
if (length(missing_cols) > 0) {
  stop("Input file missing required columns: ", paste(missing_cols, collapse = ", "))
}

mouth_df <- mouth_df %>%
  mutate(
    year = as.integer(year),
    signed_mouth_angle_deg = as.numeric(signed_mouth_angle_deg)
  )

cat("\nSigned mouth angle missing check:\n")
print(table(is.na(mouth_df[[trait_col]]), useNA = "ifany"))

# ============================================================
# Shared helper functions
# ============================================================

run_pairwise_stats_signed_angle <- function(dat, value_col, dataset_label, group_levels) {

  dat <- dat %>%
    filter(is.finite(.data[[value_col]]), !is.na(group)) %>%
    mutate(group = factor(group, levels = group_levels)) %>%
    droplevels()

  if (nrow(dat) < 3 || dplyr::n_distinct(dat$group) < 2) {
    stop("Not enough data/groups for: ", dataset_label)
  }

  # Summary stats
  summary_out <- dat %>%
    group_by(group) %>%
    summarise(
      n = sum(is.finite(.data[[value_col]])),
      mean = mean(.data[[value_col]], na.rm = TRUE),
      sd = sd(.data[[value_col]], na.rm = TRUE),
      se = sd / sqrt(n),
      min = min(.data[[value_col]], na.rm = TRUE),
      q25 = quantile(.data[[value_col]], 0.25, na.rm = TRUE),
      median = median(.data[[value_col]], na.rm = TRUE),
      q75 = quantile(.data[[value_col]], 0.75, na.rm = TRUE),
      max = max(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(dataset = dataset_label)

  # ANOVA
  fit <- aov(stats::as.formula(paste(value_col, "~ group")), data = dat)
  aov_tab <- as.data.frame(summary(fit)[[1]])
  aov_tab$term <- rownames(aov_tab)
  rownames(aov_tab) <- NULL
  aov_tab$dataset <- dataset_label

  # Kruskal-Wallis
  kw <- kruskal.test(stats::as.formula(paste(value_col, "~ group")), data = dat)
  kw_tab <- tibble(
    dataset = dataset_label,
    statistic = as.numeric(kw$statistic),
    parameter = as.numeric(kw$parameter),
    p_value = as.numeric(kw$p.value),
    method = kw$method
  )

  # Pairwise t-test, BH adjusted, same style as SICB_curves_distances_figures.R
  pw_t <- pairwise.t.test(
    x = dat[[value_col]],
    g = dat$group,
    p.adjust.method = "BH",
    pool.sd = FALSE
  )

  pw_t_df <- as.data.frame(as.table(pw_t$p.value))
  colnames(pw_t_df) <- c("group1", "group2", "p_BH")

  pw_t_df <- pw_t_df %>%
    filter(!is.na(p_BH)) %>%
    mutate(
      dataset = dataset_label,
      p_adjust_method = "BH",
      test = "pairwise.t.test"
    )

  # Pairwise Wilcoxon, BH adjusted, saved as additional robustness check
  pw_w <- pairwise.wilcox.test(
    x = dat[[value_col]],
    g = dat$group,
    p.adjust.method = "BH",
    exact = FALSE
  )

  pw_w_df <- as.data.frame(as.table(pw_w$p.value))
  colnames(pw_w_df) <- c("group1", "group2", "p_BH")

  pw_w_df <- pw_w_df %>%
    filter(!is.na(p_BH)) %>%
    mutate(
      dataset = dataset_label,
      p_adjust_method = "BH",
      test = "pairwise.wilcox.test"
    )

  list(
    data = dat,
    summary = summary_out,
    anova = aov_tab,
    kruskal = kw_tab,
    pairwise_t = pw_t_df,
    pairwise_wilcox = pw_w_df
  )
}

make_letters_from_pairwise <- function(pairwise_df, group_levels, alpha = 0.05) {

  if (nrow(pairwise_df) == 0) {
    return(tibble(
      group = factor(group_levels, levels = group_levels),
      letter = "a"
    ))
  }

  # Create named vector where TRUE = significantly different
  sig_vec <- pairwise_df$p_BH < alpha
  names(sig_vec) <- paste(pairwise_df$group1, pairwise_df$group2, sep = "-")

  letters <- multcompView::multcompLetters(sig_vec)

  out <- tibble(
    group = names(letters$Letters),
    letter = letters$Letters
  )

  missing_groups <- setdiff(group_levels, out$group)
  if (length(missing_groups) > 0) {
    out <- bind_rows(
      out,
      tibble(group = missing_groups, letter = "a")
    )
  }

  out %>%
    mutate(group = factor(group, levels = group_levels)) %>%
    arrange(group)
}

# ============================================================
# SICB-style graph-only network helpers
# ============================================================

make_node_layout <- function(group_levels) {
  data.frame(
    group = factor(group_levels, levels = group_levels),
    angle = seq(
      pi / 2,
      pi / 2 - 2 * pi + 2 * pi / length(group_levels),
      length.out = length(group_levels)
    )
  ) %>%
    mutate(
      x = cos(angle),
      y = sin(angle)
    )
}

make_graph_edges <- function(pairwise_dat, group_levels, alpha = 0.05) {

  node_layout <- make_node_layout(group_levels)

  pairwise_dat %>%
    filter(p_BH >= alpha) %>%
    transmute(
      from = factor(group1, levels = group_levels),
      to   = factor(group2, levels = group_levels),
      p_BH = p_BH
    ) %>%
    left_join(node_layout, by = c("from" = "group")) %>%
    rename(x = x, y = y) %>%
    left_join(node_layout, by = c("to" = "group"), suffix = c("", "_end")) %>%
    rename(xend = x_end, yend = y_end)
}

plot_graph_only_network <- function(edges_plot, group_levels, title_text, file_stub) {

  node_layout <- make_node_layout(group_levels)

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
      subtitle = "Edges connect groups that are not significantly different (BH-adjusted p \u2265 0.05)",
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
    filename = file.path(OUTDIR, paste0(file_stub, "_graph_network_SICB.pdf")),
    plot = p,
    width = 4,
    height = 3.2,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )

  ggsave(
    filename = file.path(OUTDIR, paste0(file_stub, "_graph_network_SICB.png")),
    plot = p,
    width = 4,
    height = 3.2,
    units = "in",
    dpi = 400,
    bg = "white"
  )

  p
}

# ============================================================
# Optional SICB-style boxplot helper
# ============================================================

plot_signed_angle_boxplot_SICB <- function(dat, group_levels, letters_df, title_text, file_stub) {

  dat <- dat %>%
    mutate(group = factor(group, levels = group_levels)) %>%
    filter(!is.na(group), is.finite(.data[[trait_col]])) %>%
    droplevels()

  letter_pos <- dat %>%
    group_by(group) %>%
    summarise(
      y_pos = max(.data[[trait_col]], na.rm = TRUE) +
        0.08 * diff(range(dat[[trait_col]], na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    left_join(letters_df, by = "group")

  y_range <- range(dat[[trait_col]], na.rm = TRUE)
  y_pad <- 0.18 * diff(y_range)
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 1

  p <- ggplot(dat, aes(x = group, y = .data[[trait_col]])) +
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
    geom_text(
      data = letter_pos,
      aes(x = group, y = y_pos, label = letter),
      inherit.aes = FALSE,
      size = 2.5
    ) +
    coord_cartesian(ylim = c(y_range[1] - 0.05 * diff(y_range), y_range[2] + y_pad)) +
    labs(
      title = title_text,
      x = NULL,
      y = "Signed mouth angle (degrees)"
    ) +
    theme_classic(base_family = "Arial", base_size = 6) +
    theme(
      plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 5),
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      plot.margin = margin(2, 2, 2, 2, unit = "pt")
    )

  ggsave(
    filename = file.path(OUTDIR, paste0(file_stub, "_boxplot_letters_SICB.pdf")),
    plot = p,
    width = 4,
    height = 3.2,
    units = "in",
    device = cairo_pdf,
    bg = "white"
  )

  ggsave(
    filename = file.path(OUTDIR, paste0(file_stub, "_boxplot_letters_SICB.png")),
    plot = p,
    width = 4,
    height = 3.2,
    units = "in",
    dpi = 400,
    bg = "white"
  )

  p
}

write_stats_text <- function(stats_obj, letters_df, outfile, title_text) {

  sink(outfile)
  cat("============================================================\n")
  cat(title_text, "\n")
  cat("============================================================\n\n")

  cat("Trait: signed_mouth_angle_deg\n")
  cat("Interpretation: positive values = more upturned; negative values = more downturned.\n")
  cat("Zero line indicates alignment with the cranial anterior-posterior reference axis.\n\n")

  cat("Group counts:\n")
  print(table(stats_obj$data$group))

  cat("\nSummary statistics:\n")
  print(stats_obj$summary)

  cat("\nANOVA:\n")
  print(stats_obj$anova)

  cat("\nKruskal-Wallis test:\n")
  print(stats_obj$kruskal)

  cat("\nPairwise t-tests with BH adjustment:\n")
  print(stats_obj$pairwise_t)

  cat("\nPairwise Wilcoxon tests with BH adjustment:\n")
  print(stats_obj$pairwise_wilcox)

  cat("\nCompact letter display from pairwise t-tests, BH-adjusted:\n")
  print(letters_df)

  cat("\nNetwork interpretation:\n")
  cat("Edges connect groups that are not significantly different based on pairwise t-tests with BH-adjusted p >= 0.05.\n")

  sink()
}

run_one_analysis <- function(dat, group_levels, dataset_label, file_stub, title_text) {

  dat <- dat %>%
    mutate(group = factor(group, levels = group_levels)) %>%
    filter(!is.na(group), is.finite(.data[[trait_col]])) %>%
    droplevels()

  stats_obj <- run_pairwise_stats_signed_angle(
    dat = dat,
    value_col = trait_col,
    dataset_label = dataset_label,
    group_levels = group_levels
  )

  letters_df <- make_letters_from_pairwise(
    pairwise_df = stats_obj$pairwise_t,
    group_levels = group_levels,
    alpha = alpha
  )

  # Save stats
  write.csv(
    stats_obj$summary,
    file.path(STAT_OUTDIR, paste0(file_stub, "_summary_SICB.csv")),
    row.names = FALSE
  )

  write.csv(
    stats_obj$anova,
    file.path(STAT_OUTDIR, paste0(file_stub, "_anova_SICB.csv")),
    row.names = FALSE
  )

  write.csv(
    stats_obj$kruskal,
    file.path(STAT_OUTDIR, paste0(file_stub, "_kruskal_SICB.csv")),
    row.names = FALSE
  )

  write.csv(
    stats_obj$pairwise_t,
    file.path(STAT_OUTDIR, paste0(file_stub, "_pairwise_t_BH_SICB.csv")),
    row.names = FALSE
  )

  write.csv(
    stats_obj$pairwise_wilcox,
    file.path(STAT_OUTDIR, paste0(file_stub, "_pairwise_wilcox_BH_SICB.csv")),
    row.names = FALSE
  )

  write.csv(
    letters_df,
    file.path(STAT_OUTDIR, paste0(file_stub, "_letters_t_BH_SICB.csv")),
    row.names = FALSE
  )

  write_stats_text(
    stats_obj = stats_obj,
    letters_df = letters_df,
    outfile = file.path(STAT_OUTDIR, paste0(file_stub, "_stats_SICB.txt")),
    title_text = title_text
  )

  # Plot network
  edges <- make_graph_edges(
    pairwise_dat = stats_obj$pairwise_t,
    group_levels = group_levels,
    alpha = alpha
  )

  p_net <- plot_graph_only_network(
    edges_plot = edges,
    group_levels = group_levels,
    title_text = paste("Non-significant network:", title_text),
    file_stub = file_stub
  )

  # Plot boxplot with letters
  p_box <- plot_signed_angle_boxplot_SICB(
    dat = dat,
    group_levels = group_levels,
    letters_df = letters_df,
    title_text = title_text,
    file_stub = file_stub
  )

  list(
    data = dat,
    stats = stats_obj,
    letters = letters_df,
    edges = edges,
    network = p_net,
    boxplot = p_box
  )
}

# ============================================================
# Build analysis datasets
# ============================================================

# 1) CT time series
ct_time_levels <- c("CT 1970", "CT 1956", "CT 1950")

ct_time_dat <- mouth_df %>%
  filter(
    habitat == "Connecticut River",
    year %in% c(1950, 1956, 1970),
    is.finite(.data[[trait_col]])
  ) %>%
  mutate(
    group = case_when(
      year == 1970 ~ "CT 1970",
      year == 1956 ~ "CT 1956",
      year == 1950 ~ "CT 1950",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = ct_time_levels)
  ) %>%
  filter(!is.na(group))

# 2) 1950 habitats
habitat_levels <- c(
  "CT 1950",
  "Quabbin",
  "Swift",
  "Fort",
  "Sawmill"
)

habitat_1950_dat <- mouth_df %>%
  filter(
    year == 1950,
    habitat %in% c("Connecticut River", "Quabbin", "Swift River", "Fort River", "Sawmill River"),
    is.finite(.data[[trait_col]])
  ) %>%
  mutate(
    group = case_when(
      habitat == "Connecticut River" ~ "CT 1950",
      habitat == "Quabbin" ~ "Quabbin",
      habitat == "Swift River" ~ "Swift",
      habitat == "Fort River" ~ "Fort",
      habitat == "Sawmill River" ~ "Sawmill",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = habitat_levels)
  ) %>%
  filter(!is.na(group))

# 3) Combined individual CT timepoints + 1950 habitats
combined_levels <- c(
  "CT 1970",
  "CT 1956",
  "CT 1950",
  "Quabbin",
  "Swift",
  "Fort",
  "Sawmill"
)

combined_dat <- mouth_df %>%
  filter(
    (
      habitat == "Connecticut River" & year %in% c(1950, 1956, 1970)
    ) |
      (
        year == 1950 & habitat %in% c("Quabbin", "Swift River", "Fort River", "Sawmill River")
      ),
    is.finite(.data[[trait_col]])
  ) %>%
  mutate(
    group = case_when(
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
  filter(!is.na(group))

# 4) CT pooled reference vs 1950 habitats
reference_levels <- c(
  "CT time series",
  "Quabbin",
  "Swift",
  "Fort",
  "Sawmill"
)

reference_dat <- mouth_df %>%
  filter(
    (
      habitat == "Connecticut River" & year %in% c(1950, 1956, 1970)
    ) |
      (
        year == 1950 & habitat %in% c("Quabbin", "Swift River", "Fort River", "Sawmill River")
      ),
    is.finite(.data[[trait_col]])
  ) %>%
  mutate(
    group = case_when(
      habitat == "Connecticut River" & year %in% c(1950, 1956, 1970) ~ "CT time series",
      habitat == "Quabbin" ~ "Quabbin",
      habitat == "Swift River" ~ "Swift",
      habitat == "Fort River" ~ "Fort",
      habitat == "Sawmill River" ~ "Sawmill",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = reference_levels)
  ) %>%
  filter(!is.na(group))

cat("\nGroup counts: CT time series\n")
print(table(ct_time_dat$group))

cat("\nGroup counts: 1950 habitats\n")
print(table(habitat_1950_dat$group))

cat("\nGroup counts: combined\n")
print(table(combined_dat$group))

cat("\nGroup counts: CT pooled reference\n")
print(table(reference_dat$group))

# ============================================================
# Run analyses
# ============================================================

res_ct_time <- run_one_analysis(
  dat = ct_time_dat,
  group_levels = ct_time_levels,
  dataset_label = "signed_mouth_angle_CT_time_series",
  file_stub = "signed_mouth_angle_CT_time_series",
  title_text = "Signed mouth angle: CT time series"
)

res_habitat_1950 <- run_one_analysis(
  dat = habitat_1950_dat,
  group_levels = habitat_levels,
  dataset_label = "signed_mouth_angle_1950_habitats",
  file_stub = "signed_mouth_angle_1950_habitats",
  title_text = "Signed mouth angle: 1950 habitats"
)

res_combined <- run_one_analysis(
  dat = combined_dat,
  group_levels = combined_levels,
  dataset_label = "signed_mouth_angle_combined_CT_timepoints_1950_habitats",
  file_stub = "signed_mouth_angle_combined_CT_timepoints_1950_habitats",
  title_text = "Signed mouth angle: CT timepoints + 1950 habitats"
)

res_reference <- run_one_analysis(
  dat = reference_dat,
  group_levels = reference_levels,
  dataset_label = "signed_mouth_angle_CT_reference_1950_habitats",
  file_stub = "signed_mouth_angle_CT_reference_1950_habitats",
  title_text = "Signed mouth angle: CT reference + 1950 habitats"
)

# ============================================================
# Save combined RDS and manifest
# ============================================================

saveRDS(
  list(
    res_ct_time = res_ct_time,
    res_habitat_1950 = res_habitat_1950,
    res_combined = res_combined,
    res_reference = res_reference
  ),
  file.path(STAT_OUTDIR, "signed_mouth_angle_SICB_results.rds")
)

manifest <- tibble::tribble(
  ~category, ~file,
  "figure", file.path(OUTDIR, "signed_mouth_angle_CT_time_series_graph_network_SICB.pdf"),
  "figure", file.path(OUTDIR, "signed_mouth_angle_1950_habitats_graph_network_SICB.pdf"),
  "figure", file.path(OUTDIR, "signed_mouth_angle_combined_CT_timepoints_1950_habitats_graph_network_SICB.pdf"),
  "figure", file.path(OUTDIR, "signed_mouth_angle_CT_reference_1950_habitats_graph_network_SICB.pdf"),
  "figure", file.path(OUTDIR, "signed_mouth_angle_CT_time_series_boxplot_letters_SICB.pdf"),
  "figure", file.path(OUTDIR, "signed_mouth_angle_1950_habitats_boxplot_letters_SICB.pdf"),
  "figure", file.path(OUTDIR, "signed_mouth_angle_combined_CT_timepoints_1950_habitats_boxplot_letters_SICB.pdf"),
  "figure", file.path(OUTDIR, "signed_mouth_angle_CT_reference_1950_habitats_boxplot_letters_SICB.pdf"),
  "text", file.path(STAT_OUTDIR, "signed_mouth_angle_CT_time_series_stats_SICB.txt"),
  "text", file.path(STAT_OUTDIR, "signed_mouth_angle_1950_habitats_stats_SICB.txt"),
  "text", file.path(STAT_OUTDIR, "signed_mouth_angle_combined_CT_timepoints_1950_habitats_stats_SICB.txt"),
  "text", file.path(STAT_OUTDIR, "signed_mouth_angle_CT_reference_1950_habitats_stats_SICB.txt")
)

write.csv(
  manifest,
  file.path(STAT_OUTDIR, "MANIFEST_signed_mouth_angle_SICB.csv"),
  row.names = FALSE
)

cat("\nSaved SICB signed mouth angle outputs to:\n")
cat(normalizePath(OUTDIR), "\n\n")

cat("Saved SICB signed mouth angle stats to:\n")
cat(normalizePath(STAT_OUTDIR), "\n\n")

cat("Files created:\n")
print(list.files(OUTDIR))
cat("\nStats files created:\n")
print(list.files(STAT_OUTDIR))

cat("\nDone.\n")
