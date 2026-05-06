# ============================================================
# size_analysis_SICB.R
# ============================================================
# SIZE ANALYSIS: CENTROID SIZE (Csize)
# 1) 1950 waterbodies
# 2) Connecticut River time series
# 3) Full combined dataset
#
# Includes:
#   - boxplots
#   - Kruskal-Wallis tests
#   - pairwise Wilcoxon tests with BH correction
#   - compact letter displays
#   - min / median / max summaries
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(multcompView)
})

# ---- Load base objects ----
if (!exists("gdf")) source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")

OUTDIR <- file.path("Outputs", "size_analysis")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Helper functions
# ============================================================

get_size_summary <- function(df, group_col, dataset_label) {
  
  df %>%
    filter(!is.na(Csize), !is.na(.data[[group_col]])) %>%
    group_by(across(all_of(group_col))) %>%
    summarise(
      n = n(),
      min_Csize = min(Csize, na.rm = TRUE),
      median_Csize = median(Csize, na.rm = TRUE),
      max_Csize = max(Csize, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(group = all_of(group_col)) %>%
    mutate(dataset = dataset_label, .before = 1)
}

run_size_tests <- function(df, group_col, out_prefix) {
  
  df <- df %>%
    filter(!is.na(Csize), !is.na(.data[[group_col]])) %>%
    mutate(group_tmp = factor(.data[[group_col]]))
  
  all_groups <- levels(df$group_tmp)
  
  # ---- Overall Kruskal-Wallis test ----
  kw <- kruskal.test(Csize ~ group_tmp, data = df)
  
  capture.output(
    kw,
    file = file.path(OUTDIR, paste0(out_prefix, "_kruskal_wallis.txt"))
  )
  
  # ---- Pairwise Wilcoxon tests with BH correction ----
  pw <- pairwise.wilcox.test(
    x = df$Csize,
    g = df$group_tmp,
    p.adjust.method = "BH",
    exact = FALSE
  )
  
  capture.output(
    pw,
    file = file.path(OUTDIR, paste0(out_prefix, "_pairwise_wilcoxon_BH.txt"))
  )
  
  # ---- Compact letter display ----
  pmat <- pw$p.value
  
  # multcompLetters can work from the pairwise p-value matrix,
  # but the returned letters sometimes omit one group because the
  # pairwise matrix is triangular. This section forces all groups back in.
  letters_raw <- multcompView::multcompLetters(
    pmat,
    threshold = 0.05
  )$Letters
  
  letters_df <- data.frame(
    group = names(letters_raw),
    letters = unname(letters_raw),
    stringsAsFactors = FALSE
  )
  
  # Ensure every group receives a row and a letter
  letters_df <- data.frame(
    group = all_groups,
    stringsAsFactors = FALSE
  ) %>%
    left_join(letters_df, by = "group") %>%
    mutate(
      letters = ifelse(is.na(letters), "a", letters)
    )
  
  write.csv(
    letters_df,
    file.path(OUTDIR, paste0(out_prefix, "_compact_letters.csv")),
    row.names = FALSE
  )
  
  # ---- Save diagnostic showing all groups were represented ----
  diagnostics <- list(
    all_groups = all_groups,
    pvalue_matrix_rownames = rownames(pmat),
    pvalue_matrix_colnames = colnames(pmat),
    compact_letters = letters_df
  )
  
  capture.output(
    diagnostics,
    file = file.path(OUTDIR, paste0(out_prefix, "_letter_diagnostics.txt"))
  )
  
  return(list(
    kruskal = kw,
    pairwise = pw,
    letters = letters_df
  ))
}

make_boxplot <- function(df, group_col, title, outfile, letters_df = NULL) {
  
  df <- df %>%
    filter(!is.na(Csize), !is.na(.data[[group_col]])) %>%
    mutate(group_plot = factor(.data[[group_col]]))
  
  p <- ggplot(df, aes(x = group_plot, y = Csize)) +
    geom_boxplot(outlier.size = 1.2, linewidth = 0.35) +
    geom_jitter(width = 0.12, alpha = 0.55, size = 1.4) +
    theme_classic(base_size = 9) +
    labs(
      title = title,
      x = "",
      y = "Centroid size"
    ) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(size = 10, face = "bold")
    )
  
  if (!is.null(letters_df)) {
    
    y_positions <- df %>%
      group_by(group_plot) %>%
      summarise(
        y = max(Csize, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        y = y + 0.05 * diff(range(df$Csize, na.rm = TRUE)),
        group = as.character(group_plot)
      ) %>%
      select(group, y)
    
    letter_plot_df <- letters_df %>%
      mutate(
        group = as.character(group),
        group_plot = factor(group, levels = levels(df$group_plot))
      ) %>%
      left_join(y_positions, by = "group")
    
    p <- p +
      geom_text(
        data = letter_plot_df,
        aes(x = group_plot, y = y, label = letters),
        inherit.aes = FALSE,
        size = 3.2
      )
  }
  
  ggsave(outfile, p, width = 4.2, height = 3.2, dpi = 300)
  return(p)
}

# ============================================================
# 1) 1950 WATERBODIES ONLY
# ============================================================

df_1950 <- gdf %>%
  filter(year == 1950, !is.na(habitat)) %>%
  mutate(
    habitat = factor(
      habitat,
      levels = c(
        "Connecticut River",
        "Fort River",
        "Quabbin",
        "Sawmill River",
        "Swift River"
      )
    )
  )

summary_1950 <- get_size_summary(
  df_1950,
  group_col = "habitat",
  dataset_label = "1950_waterbodies"
)

write.csv(
  summary_1950,
  file.path(OUTDIR, "summary_1950_waterbodies.csv"),
  row.names = FALSE
)

tests_1950 <- run_size_tests(
  df_1950,
  group_col = "habitat",
  out_prefix = "1950_waterbodies"
)

make_boxplot(
  df_1950,
  group_col = "habitat",
  title = "Centroid size across 1950 waterbodies",
  outfile = file.path(OUTDIR, "boxplot_1950_waterbodies_with_letters.png"),
  letters_df = tests_1950$letters
)

# ============================================================
# 2) CONNECTICUT RIVER TIME SERIES ONLY
# ============================================================

df_ct <- gdf %>%
  filter(
    habitat == "Connecticut River",
    year %in% c(1950, 1956, 1970)
  ) %>%
  mutate(
    year = factor(year, levels = c(1950, 1956, 1970))
  )

summary_ct <- get_size_summary(
  df_ct,
  group_col = "year",
  dataset_label = "CT_time_series"
)

write.csv(
  summary_ct,
  file.path(OUTDIR, "summary_CT_time_series.csv"),
  row.names = FALSE
)

tests_ct <- run_size_tests(
  df_ct,
  group_col = "year",
  out_prefix = "CT_time_series"
)

make_boxplot(
  df_ct,
  group_col = "year",
  title = "Centroid size across Connecticut River time series",
  outfile = file.path(OUTDIR, "boxplot_CT_time_series_with_letters.png"),
  letters_df = tests_ct$letters
)

# ============================================================
# 3) FULL COMBINED DATASET
# ============================================================

df_full <- gdf %>%
  filter(!is.na(habitat)) %>%
  mutate(
    group = case_when(
      habitat == "Connecticut River" & year == 1950 ~ "CT_1950",
      habitat == "Connecticut River" & year == 1956 ~ "CT_1956",
      habitat == "Connecticut River" & year == 1970 ~ "CT_1970",
      year == 1950 ~ paste0(habitat, "_1950"),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group)) %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "CT_1950",
        "CT_1956",
        "CT_1970",
        "Fort River_1950",
        "Quabbin_1950",
        "Sawmill River_1950",
        "Swift River_1950"
      )
    )
  )

summary_full <- get_size_summary(
  df_full,
  group_col = "group",
  dataset_label = "full_combined_dataset"
)

write.csv(
  summary_full,
  file.path(OUTDIR, "summary_full_combined_dataset.csv"),
  row.names = FALSE
)

tests_full <- run_size_tests(
  df_full,
  group_col = "group",
  out_prefix = "full_combined_dataset"
)

make_boxplot(
  df_full,
  group_col = "group",
  title = "Centroid size across full combined dataset",
  outfile = file.path(OUTDIR, "boxplot_full_combined_dataset_with_letters.png"),
  letters_df = tests_full$letters
)

# ============================================================
# DONE
# ============================================================

cat("\nCentroid size analysis complete.\n")
cat("Outputs saved to:", OUTDIR, "\n")