# ============================================================
# Scripts/mouth_angle_trait_analysis_signed.R
#
# Analyze SIGNED mouth angle across:
#   1) Connecticut River time series: CT_1950, CT_1956, CT_1970
#   2) 1950 waterbodies/habitats: CT, Quabbin, Swift, Fort, Sawmill
#   3) CT time series pooled as reference vs. 1950 non-CT habitats
#
# Input:
#   trait_measurements/mouth_angle_signed_with_metadata.csv
#
# Outputs:
#   Figures/mouth_angle_trait_analysis_signed/
#   Outputs/mouth_angle_trait_analysis_signed_<runid>/
#
# Notes:
#   - Boxplots include compact letters from BH-adjusted pairwise Wilcoxon tests.
#   - A dashed zero line marks the transition between signed directions.
#   - Network plots show NON-significant pairwise connections, matching your SICB style:
#       edges connect groups that are not significantly different (BH-adjusted p >= 0.05).
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
# Inputs / outputs
# ---------------------------
infile <- file.path("trait_measurements", "mouth_angle_signed_with_metadata.csv")

if (!file.exists(infile)) {
  stop(
    "Could not find: ", infile,
    "\nRun source('R/build_trait_measurements_signed.R') first."
  )
}

mouth_df <- read.csv(infile, stringsAsFactors = FALSE)

required_cols <- c("specimen", "habitat", "year", "signed_mouth_angle_deg")
missing_cols <- setdiff(required_cols, names(mouth_df))
if (length(missing_cols) > 0) {
  stop("Input file missing required columns: ", paste(missing_cols, collapse = ", "))
}

mouth_df <- mouth_df %>%
  mutate(
    year = as.integer(year),
    signed_mouth_angle_deg = as.numeric(signed_mouth_angle_deg)
  )

run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")

fig_dir <- file.path("Figures", "mouth_angle_trait_analysis_signed")
out_root <- file.path("Outputs", paste0("mouth_angle_trait_analysis_signed_", run_id))
tables_dir <- file.path(out_root, "tables")
text_dir <- file.path(out_root, "text")
rds_dir <- file.path(out_root, "rds")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(text_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

alpha <- 0.05
trait_col <- "signed_mouth_angle_deg"

cat("\nSigned mouth angle merge check, all metadata:\n")
print(table(is.na(mouth_df[[trait_col]]), useNA = "ifany"))

# ---------------------------
# Palettes
# ---------------------------
group_palette_full <- c(
  "CT_1950" = "steelblue",
  "CT_1956" = "deepskyblue3",
  "CT_1970" = "navy",
  "Fort River_1950" = "black",
  "Quabbin_1950" = "darkgoldenrod2",
  "Sawmill River_1950" = "orchid3",
  "Swift River_1950" = "tomato"
)

habitat_palette <- c(
  "Connecticut River" = "steelblue",
  "Quabbin" = "darkgoldenrod2",
  "Swift River" = "tomato",
  "Fort River" = "black",
  "Sawmill River" = "orchid3"
)

landscape_palette <- c(
  "CT_timeseries" = "steelblue4",
  "Quabbin_1950" = "darkgoldenrod2",
  "Swift River_1950" = "tomato",
  "Fort River_1950" = "black",
  "Sawmill River_1950" = "orchid3"
)

# ============================================================
# Helper functions
# ============================================================

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
    summarize(
      y_pos = max(!!trait_col, na.rm = TRUE) + pad,
      .groups = "drop"
    ) %>%
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
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4, alpha = 0.8) +
    geom_boxplot(width = 0.72, outlier.shape = NA, alpha = 0.9) +
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
      y = "Signed mouth angle relative to cranial axis (degrees)"
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 35, hjust = 1)
    )

  ggsave(file.path(fig_dir, paste0(outfile_base, ".png")), p, width = 7.2, height = 5, dpi = 400, bg = "white")
  ggsave(file.path(fig_dir, paste0(outfile_base, ".pdf")), p, width = 7.2, height = 5, bg = "white")

  p
}

make_nonsig_network_plot <- function(pairwise_df, group_levels, palette, title, outfile_base,
                                     alpha = 0.05, layout_type = "fr") {
  group_levels <- as.character(group_levels)

  # Edges connect groups that are NOT significantly different.
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
    geom_edge_link(
      linewidth = 0.55,
      alpha = 0.60,
      color = "gray50"
    ) +
    geom_node_point(
      aes(fill = name),
      shape = 21,
      size = 6.5,
      color = "black",
      stroke = 0.7
    ) +
    geom_node_text(
      aes(label = name),
      repel = TRUE,
      size = 3.2,
      color = "black"
    ) +
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
      caption = "Edges connect groups that are not significantly different (BH-adjusted p \u2265 0.05)"
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
  cat("Trait interpretation:\n")
  cat("  Positive values = more upturned mouth orientation, after sign check.\n")
  cat("  Negative values = more downturned mouth orientation, after sign check.\n")
  cat("  Values near zero = mouth axis more closely aligned with cranial axis.\n\n")
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

run_signed_mouth_angle_analysis <- function(df, group_col, analysis_label, outfile_prefix, palette) {
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
    title = "Non-significant network: Signed mouth angle",
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

# ============================================================
# Build analysis datasets
# ============================================================

# Analysis 1: Connecticut River time series
ct_time_df <- mouth_df %>%
  filter(
    habitat == "Connecticut River",
    year %in% c(1950, 1956, 1970)
  ) %>%
  mutate(
    ct_time_group = case_when(
      year == 1950 ~ "CT_1950",
      year == 1956 ~ "CT_1956",
      year == 1970 ~ "CT_1970",
      TRUE ~ NA_character_
    ),
    ct_time_group = factor(
      ct_time_group,
      levels = c("CT_1950", "CT_1956", "CT_1970")
    )
  )

# Analysis 2: 1950 habitats/waterbodies
hab_1950_df <- mouth_df %>%
  filter(
    !is.na(habitat),
    year == 1950,
    habitat %in% c("Connecticut River", "Quabbin", "Swift River", "Fort River", "Sawmill River")
  ) %>%
  mutate(
    habitat_1950 = factor(
      habitat,
      levels = c("Connecticut River", "Quabbin", "Swift River", "Fort River", "Sawmill River")
    )
  )

# Analysis 3: Pooled CT time series vs 1950 habitats
landscape_df <- mouth_df %>%
  filter(
    !is.na(habitat),
    (
      (habitat == "Connecticut River" & year %in% c(1950, 1956, 1970)) |
        (habitat != "Connecticut River" & year == 1950)
    )
  ) %>%
  mutate(
    landscape_group = case_when(
      habitat == "Connecticut River" & year %in% c(1950, 1956, 1970) ~ "CT_timeseries",
      habitat == "Quabbin" & year == 1950 ~ "Quabbin_1950",
      habitat == "Swift River" & year == 1950 ~ "Swift River_1950",
      habitat == "Fort River" & year == 1950 ~ "Fort River_1950",
      habitat == "Sawmill River" & year == 1950 ~ "Sawmill River_1950",
      TRUE ~ NA_character_
    ),
    landscape_group = factor(
      landscape_group,
      levels = c(
        "CT_timeseries",
        "Quabbin_1950",
        "Swift River_1950",
        "Fort River_1950",
        "Sawmill River_1950"
      )
    )
  ) %>%
  filter(!is.na(landscape_group))

cat("\nCT time series counts:\n")
print(table(ct_time_df$ct_time_group, useNA = "ifany"))

cat("\n1950 habitat counts:\n")
print(table(hab_1950_df$habitat_1950, useNA = "ifany"))

cat("\nLandscape counts:\n")
print(table(landscape_df$landscape_group, useNA = "ifany"))

# ============================================================
# Run analyses
# ============================================================

res_ct_time <- run_signed_mouth_angle_analysis(
  df = ct_time_df,
  group_col = ct_time_group,
  analysis_label = "Signed mouth angle across Connecticut River time series",
  outfile_prefix = "CT_timeseries_signed_mouth_angle",
  palette = group_palette_full
)

res_hab_1950 <- run_signed_mouth_angle_analysis(
  df = hab_1950_df,
  group_col = habitat_1950,
  analysis_label = "Signed mouth angle across 1950 waterbodies",
  outfile_prefix = "habitats_1950_signed_mouth_angle",
  palette = habitat_palette
)

res_landscape <- run_signed_mouth_angle_analysis(
  df = landscape_df,
  group_col = landscape_group,
  analysis_label = "Signed mouth angle: CT time series reference vs. 1950 habitats",
  outfile_prefix = "CT_timeseries_vs_1950habitats_signed_mouth_angle",
  palette = landscape_palette
)

saveRDS(
  list(
    res_ct_time = res_ct_time,
    res_hab_1950 = res_hab_1950,
    res_landscape = res_landscape
  ),
  file.path(rds_dir, "signed_mouth_angle_analysis_results.rds")
)

# ============================================================
# Manifest / run info
# ============================================================

manifest <- tibble::tribble(
  ~category, ~file,
  "text",  file.path(text_dir, "CT_timeseries_signed_mouth_angle_statistics.txt"),
  "text",  file.path(text_dir, "habitats_1950_signed_mouth_angle_statistics.txt"),
  "text",  file.path(text_dir, "CT_timeseries_vs_1950habitats_signed_mouth_angle_statistics.txt"),
  "table", file.path(tables_dir, "CT_timeseries_signed_mouth_angle_summary.csv"),
  "table", file.path(tables_dir, "CT_timeseries_signed_mouth_angle_pairwise_wilcox_BH.csv"),
  "table", file.path(tables_dir, "CT_timeseries_signed_mouth_angle_letters_BH.csv"),
  "table", file.path(tables_dir, "habitats_1950_signed_mouth_angle_summary.csv"),
  "table", file.path(tables_dir, "habitats_1950_signed_mouth_angle_pairwise_wilcox_BH.csv"),
  "table", file.path(tables_dir, "habitats_1950_signed_mouth_angle_letters_BH.csv"),
  "table", file.path(tables_dir, "CT_timeseries_vs_1950habitats_signed_mouth_angle_summary.csv"),
  "table", file.path(tables_dir, "CT_timeseries_vs_1950habitats_signed_mouth_angle_pairwise_wilcox_BH.csv"),
  "table", file.path(tables_dir, "CT_timeseries_vs_1950habitats_signed_mouth_angle_letters_BH.csv"),
  "figure", file.path(fig_dir, "CT_timeseries_signed_mouth_angle_boxplot_letters.png"),
  "figure", file.path(fig_dir, "habitats_1950_signed_mouth_angle_boxplot_letters.png"),
  "figure", file.path(fig_dir, "CT_timeseries_vs_1950habitats_signed_mouth_angle_boxplot_letters.png"),
  "figure", file.path(fig_dir, "CT_timeseries_signed_mouth_angle_network_nonsignificant.png"),
  "figure", file.path(fig_dir, "habitats_1950_signed_mouth_angle_network_nonsignificant.png"),
  "figure", file.path(fig_dir, "CT_timeseries_vs_1950habitats_signed_mouth_angle_network_nonsignificant.png"),
  "rds", file.path(rds_dir, "signed_mouth_angle_analysis_results.rds")
)

write.csv(manifest, file.path(out_root, "MANIFEST.csv"), row.names = FALSE)

sink(file.path(out_root, "RUN_INFO.txt"))
cat("Run ID:", run_id, "\n")
cat("Working directory:", getwd(), "\n")
cat("Input file:", infile, "\n")
cat("Date:", as.character(Sys.time()), "\n\n")

cat("Signed mouth angle missing check, all specimens:\n")
print(table(is.na(mouth_df[[trait_col]]), useNA = "ifany"))

cat("\nCT time series counts:\n")
print(table(ct_time_df$ct_time_group, useNA = "ifany"))

cat("\n1950 habitat counts:\n")
print(table(hab_1950_df$habitat_1950, useNA = "ifany"))

cat("\nLandscape counts:\n")
print(table(landscape_df$landscape_group, useNA = "ifany"))

cat("\nSession info:\n")
print(sessionInfo())
sink()

cat("\n============================================================\n")
cat("Done.\n")
cat("Figures written to:\n")
cat("  ", fig_dir, "\n")
cat("Outputs written to:\n")
cat("  ", out_root, "\n")
cat("============================================================\n")
