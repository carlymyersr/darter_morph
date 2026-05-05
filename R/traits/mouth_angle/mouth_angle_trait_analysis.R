# ============================================================
# Scripts/mouth_angle_trait_analysis.R
#
# Analyze mouth angle across:
#   1) Connecticut River time series: CT_1950, CT_1956, CT_1970
#   2) 1950 waterbodies/habitats: CT, Quabbin, Swift, Fort, Sawmill
#   3) CT time series pooled as reference vs. 1950 non-CT habitats
#
# Inputs:
#   trait_measurements/mouth_angle_raw.csv
#
# Outputs:
#   Figures/mouth_angle_trait_analysis/
#   Outputs/mouth_angle_trait_analysis_<runid>/
# ============================================================

mouth_df <- read.csv("trait_measurements/mouth_angle_with_metadata.csv")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(multcompView)
})

set.seed(123)
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
fig_dir <- file.path("Figures", "mouth_angle_trait_analysis")
out_root <- file.path("Outputs", paste0("mouth_angle_trait_analysis_", run_id))
tables_dir <- file.path(out_root, "tables")
text_dir <- file.path(out_root, "text")
rds_dir <- file.path(out_root, "rds")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(text_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

alpha <- 0.05

source("R/00_setup_morpho.R")
source("R/01_build_metadata.R")
source("R/04_subset_CT_timeseries_plus_1950habitats.R")

angle_raw_file <- file.path("trait_measurements", "mouth_angle_raw.csv")
angle_merged_file <- file.path("trait_measurements", "mouth_angle_with_metadata.csv")

if (file.exists(angle_raw_file)) {
  angle_df <- read.csv(angle_raw_file, stringsAsFactors = FALSE)
} else if (file.exists(angle_merged_file)) {
  angle_df <- read.csv(angle_merged_file, stringsAsFactors = FALSE) %>%
    dplyr::select(specimen, mouth_angle_deg)
} else {
  stop("No mouth angle file found. Run R/build_trait_measurements.R first.")
}

if (!all(c("specimen", "mouth_angle_deg") %in% names(angle_df))) {
  stop("Mouth angle file must contain columns: specimen, mouth_angle_deg")
}

angle_df <- angle_df %>%
  mutate(specimen = as.character(specimen),
         mouth_angle_deg = as.numeric(mouth_angle_deg))

gdf_angle <- gdf %>% left_join(angle_df, by = "specimen")
gdf_Fig6_angle <- gdf_Fig6 %>% left_join(angle_df, by = "specimen")

cat("\nMouth angle merge check, all metadata:\n")
print(table(is.na(gdf_angle$mouth_angle_deg), useNA = "ifany"))
cat("\nMouth angle merge check, Fig6/Fig7 subset:\n")
print(table(is.na(gdf_Fig6_angle$mouth_angle_deg), useNA = "ifany"))

write.csv(gdf_angle, file.path(tables_dir, "metadata_with_mouth_angle_all_specimens.csv"), row.names = FALSE)
write.csv(gdf_Fig6_angle, file.path(tables_dir, "metadata_with_mouth_angle_CTtimeseries_1950habitats.csv"), row.names = FALSE)

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

  pw <- pairwise.wilcox.test(
    x = df2[[trait_col_name]],
    g = df2$.group,
    p.adjust.method = "BH",
    exact = FALSE
  )

  as.data.frame(as.table(pw$p.value), stringsAsFactors = FALSE) %>%
    filter(!is.na(Freq)) %>%
    transmute(group1 = as.character(Var1),
              group2 = as.character(Var2),
              p_adj = as.numeric(Freq)) %>%
    arrange(p_adj)
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
  out %>% mutate(group = factor(group, levels = group_levels)) %>% arrange(group)
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

make_boxplot_with_letters <- function(df, group_col, trait_col, letters_df, palette, title, subtitle, outfile_base) {
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
  y_pad <- diff(y_range) * 0.18
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 1

  p <- ggplot(plot_df, aes(x = !!group_col, y = !!trait_col, fill = !!group_col)) +
    geom_boxplot(width = 0.72, outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.12, size = 1.7, alpha = 0.75) +
    geom_text(data = letter_pos, aes(x = group, y = y_pos, label = letters),
              inherit.aes = FALSE, size = 4) +
    scale_fill_manual(values = palette[group_levels], drop = FALSE) +
    coord_cartesian(ylim = c(y_range[1], y_range[2] + y_pad)) +
    labs(title = title, subtitle = subtitle, x = NULL,
         y = "Mouth angle relative to cranial axis (degrees)") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 35, hjust = 1))

  ggsave(file.path(fig_dir, paste0(outfile_base, ".png")), p, width = 7.2, height = 5, dpi = 400)
  ggsave(file.path(fig_dir, paste0(outfile_base, ".pdf")), p, width = 7.2, height = 5)
  p
}

make_pairwise_network <- function(pairwise_df, group_levels, palette, title, outfile_base, alpha = 0.05) {
  group_levels <- as.character(group_levels)
  n <- length(group_levels)
  theta <- seq(0, 2*pi, length.out = n + 1)[1:n]
  nodes <- tibble(group = group_levels, x = cos(theta), y = sin(theta))

  edges <- pairwise_df %>%
    filter(!is.na(p_adj), p_adj < alpha) %>%
    left_join(nodes %>% rename(group1 = group, x1 = x, y1 = y), by = "group1") %>%
    left_join(nodes %>% rename(group2 = group, x2 = x, y2 = y), by = "group2") %>%
    mutate(p_label = ifelse(p_adj < 0.001, "p<0.001", paste0("p=", signif(p_adj, 2))))

  p <- ggplot() + coord_equal() + theme_void(base_size = 12) + labs(title = title)

  if (nrow(edges) > 0) {
    p <- p +
      geom_segment(data = edges, aes(x = x1, y = y1, xend = x2, yend = y2),
                   linewidth = 0.8, alpha = 0.75) +
      geom_text(data = edges, aes(x = (x1 + x2)/2, y = (y1 + y2)/2, label = p_label),
                size = 3, vjust = -0.4)
  } else {
    p <- p + annotate("text", x = 0, y = 0,
                      label = paste0("No BH-adjusted pairwise differences at alpha = ", alpha),
                      size = 4)
  }

  p <- p +
    geom_point(data = nodes, aes(x = x, y = y, fill = group),
               shape = 21, size = 8, color = "black", stroke = 0.8) +
    geom_text(data = nodes, aes(x = x * 1.22, y = y * 1.22, label = group), size = 3.4) +
    scale_fill_manual(values = palette[group_levels], drop = FALSE) +
    theme(legend.position = "none")

  ggsave(file.path(fig_dir, paste0(outfile_base, ".png")), p, width = 6.5, height = 6.5, dpi = 400)
  ggsave(file.path(fig_dir, paste0(outfile_base, ".pdf")), p, width = 6.5, height = 6.5)
  p
}

write_stats_txt <- function(df, group_col, trait_col, summary_df, pairwise_df, letters_df, analysis_name, outfile) {
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
  sink()
}

run_mouth_angle_analysis <- function(df, group_col, analysis_label, outfile_prefix, palette) {
  group_col_name <- rlang::as_name(rlang::ensym(group_col))

  df_use <- df %>%
    filter(!is.na(.data[[group_col_name]]), !is.na(mouth_angle_deg)) %>%
    mutate(!!group_col_name := droplevels(factor(.data[[group_col_name]])))

  group_levels <- levels(df_use[[group_col_name]])
  summary_df <- summarize_trait(df_use, !!rlang::sym(group_col_name), mouth_angle_deg)
  pairwise_df <- pairwise_wilcox_bh(df_use, !!rlang::sym(group_col_name), mouth_angle_deg)
  letters_df <- make_letters_from_pairwise(pairwise_df, group_levels, "p_adj", alpha)

  write.csv(summary_df, file.path(tables_dir, paste0(outfile_prefix, "_summary.csv")), row.names = FALSE)
  write.csv(pairwise_df, file.path(tables_dir, paste0(outfile_prefix, "_pairwise_wilcox_BH.csv")), row.names = FALSE)
  write.csv(letters_df, file.path(tables_dir, paste0(outfile_prefix, "_letters_BH.csv")), row.names = FALSE)

  write_stats_txt(df_use, !!rlang::sym(group_col_name), mouth_angle_deg,
                  summary_df, pairwise_df, letters_df, analysis_label,
                  file.path(text_dir, paste0(outfile_prefix, "_statistics.txt")))

  p_box <- make_boxplot_with_letters(
    df_use, !!rlang::sym(group_col_name), mouth_angle_deg, letters_df, palette,
    analysis_label, "Letters show BH-adjusted pairwise Wilcoxon groupings",
    paste0(outfile_prefix, "_boxplot_letters")
  )

  p_net <- make_pairwise_network(
    pairwise_df, group_levels, palette,
    paste0(analysis_label, ": significant pairwise differences"),
    paste0(outfile_prefix, "_network_significant_pairwise"), alpha
  )

  list(data = df_use, summary = summary_df, pairwise = pairwise_df,
       letters = letters_df, boxplot = p_box, network = p_net)
}

# ============================================================
# Analysis 1: Connecticut River time series
# ============================================================

ct_time_df <- gdf_angle %>%
  filter(habitat == "Connecticut River", year %in% c(1950, 1956, 1970)) %>%
  mutate(
    ct_time_group = case_when(
      year == 1950 ~ "CT_1950",
      year == 1956 ~ "CT_1956",
      year == 1970 ~ "CT_1970",
      TRUE ~ NA_character_
    ),
    ct_time_group = factor(ct_time_group, levels = c("CT_1950", "CT_1956", "CT_1970"))
  )

res_ct_time <- run_mouth_angle_analysis(
  ct_time_df, ct_time_group,
  "Mouth angle across Connecticut River time series",
  "CT_timeseries_mouth_angle",
  group_palette_full
)

# ============================================================
# Analysis 2: 1950 habitats / waterbodies
# ============================================================

hab_1950_df <- gdf_angle %>%
  filter(!is.na(habitat), year == 1950,
         habitat %in% c("Connecticut River", "Quabbin", "Swift River", "Fort River", "Sawmill River")) %>%
  mutate(habitat_1950 = factor(
    habitat,
    levels = c("Connecticut River", "Quabbin", "Swift River", "Fort River", "Sawmill River")
  ))

res_hab_1950 <- run_mouth_angle_analysis(
  hab_1950_df, habitat_1950,
  "Mouth angle across 1950 waterbodies",
  "habitats_1950_mouth_angle",
  habitat_palette
)

# ============================================================
# Analysis 3: Pooled CT time series vs. 1950 habitats
# ============================================================

landscape_df <- gdf_angle %>%
  filter(!is.na(habitat),
         ((habitat == "Connecticut River" & year %in% c(1950, 1956, 1970)) |
            (habitat != "Connecticut River" & year == 1950))) %>%
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
      levels = c("CT_timeseries", "Quabbin_1950", "Swift River_1950",
                 "Fort River_1950", "Sawmill River_1950")
    )
  ) %>%
  filter(!is.na(landscape_group))

res_landscape <- run_mouth_angle_analysis(
  landscape_df, landscape_group,
  "Mouth angle: CT time series reference vs. 1950 habitats",
  "CT_timeseries_vs_1950habitats_mouth_angle",
  landscape_palette
)

saveRDS(
  list(res_ct_time = res_ct_time, res_hab_1950 = res_hab_1950, res_landscape = res_landscape),
  file.path(rds_dir, "mouth_angle_analysis_results.rds")
)




manifest <- tibble::tribble(
  ~category, ~file,
  "table", file.path(tables_dir, "metadata_with_mouth_angle_all_specimens.csv"),
  "table", file.path(tables_dir, "metadata_with_mouth_angle_CTtimeseries_1950habitats.csv"),
  "text",  file.path(text_dir, "CT_timeseries_mouth_angle_statistics.txt"),
  "text",  file.path(text_dir, "habitats_1950_mouth_angle_statistics.txt"),
  "text",  file.path(text_dir, "CT_timeseries_vs_1950habitats_mouth_angle_statistics.txt"),
  "figure", file.path(fig_dir, "CT_timeseries_mouth_angle_boxplot_letters.png"),
  "figure", file.path(fig_dir, "habitats_1950_mouth_angle_boxplot_letters.png"),
  "figure", file.path(fig_dir, "CT_timeseries_vs_1950habitats_mouth_angle_boxplot_letters.png"),
  "figure", file.path(fig_dir, "CT_timeseries_mouth_angle_network_significant_pairwise.png"),
  "figure", file.path(fig_dir, "habitats_1950_mouth_angle_network_significant_pairwise.png"),
  "figure", file.path(fig_dir, "CT_timeseries_vs_1950habitats_mouth_angle_network_significant_pairwise.png"),
  "rds", file.path(rds_dir, "mouth_angle_analysis_results.rds")
)

write.csv(manifest, file.path(out_root, "MANIFEST.csv"), row.names = FALSE)

sink(file.path(out_root, "RUN_INFO.txt"))
cat("Run ID:", run_id, "\n")
cat("Working directory:", getwd(), "\n")
cat("Date:", as.character(Sys.time()), "\n\n")
cat("Mouth angle merge check, all specimens:\n")
print(table(is.na(gdf_angle$mouth_angle_deg), useNA = "ifany"))
cat("\nMouth angle merge check, Fig6 subset:\n")
print(table(is.na(gdf_Fig6_angle$mouth_angle_deg), useNA = "ifany"))
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
