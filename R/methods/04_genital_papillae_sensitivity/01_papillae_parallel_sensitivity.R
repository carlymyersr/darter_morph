# ============================================================
# papillae_parallel_analyses.R
#
# Goal:
#   Run two parallel papillae scoring analyses:
#     1) conservative high-confidence-only analysis
#     2) inclusive sensitivity analysis
#
# For each analysis, run:
#   - all data with papillae scores
#   - 1950 habitats
#   - CT time series
#
# For each dataset, run:
#   - raw GPA shape
#   - size-corrected shape residuals
#   - PCA hull plots
#   - procD.lm tests of shape ~ papilla_group
#
# Outputs are saved under:
#   papillae/
#
# Expected papillae CSV columns:
#   name, papilla_presence, papilla_confidence
#
# Example scoring assumptions:
#   Conservative:
#     Present + High  -> papilla_present
#     Absent  + High  -> no_papilla_observed
#     Low confidence rows excluded
#
#   Inclusive sensitivity:
#     Present High/Low -> papilla_present
#     Absent  High/Low -> no_papilla_observed
#
# Notes:
#   - This script assumes specimen names in the papillae CSV match the
#     specimen IDs in coords_gpa/gdf, e.g. F2380_1_side.
#   - If you later use "PH", "PL", "AH", "not observed", etc.,
#     the normalization functions below should still handle them.
# ============================================================

# ---------------------------
# User settings
# ---------------------------

PAPILLAE_CSV <- "papillae.csv"   # change if needed
OUTDIR       <- "papillae"
N_PERM       <- 999
SAVE_WIDTH   <- 3
SAVE_HEIGHT  <- 2.36
SAVE_DPI     <- 600

# Which PC panels to save
PC_PANELS <- list(
  c("PC1", "PC2"),
  c("PC2", "PC3"),
  c("PC3", "PC4")
)

# ---------------------------
# Libraries
# ---------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(geomorph)
  library(RRPP)
  library(tibble)
  library(purrr)
})

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Source canonical morphometric pipeline
# ---------------------------

if (!exists("coords_gpa", inherits = TRUE)) {
  source("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R")
}

if (!exists("gdf", inherits = TRUE)) {
  source("R/methods/01_specimen_sampling_study_design/01_build_metadata.R")
}

# Some older scripts expect this helper, but it may not be defined globally.
subset_coords_to_gdf <- function(coords_arr, gdf_sub) {
  stopifnot("specimen" %in% names(gdf_sub))
  missing_ids <- setdiff(gdf_sub$specimen, dimnames(coords_arr)[[3]])
  if (length(missing_ids) > 0) {
    stop("These specimens are in metadata but missing from coords array: ",
         paste(missing_ids, collapse = ", "))
  }
  coords_arr[, , gdf_sub$specimen, drop = FALSE]
}

# Robust allometry helper, in case the sourced metadata script is unavailable
# or you want this file to be maximally portable.
if (!exists("allometry_residuals", inherits = TRUE)) {
  allometry_residuals <- function(coords_arr, size_vec) {
    stopifnot(length(dim(coords_arr)) == 3)
    stopifnot(length(size_vec) == dim(coords_arr)[3])

    size_vec <- as.numeric(size_vec)
    if (!is.null(names(size_vec))) {
      stopifnot(identical(names(size_vec), dimnames(coords_arr)[[3]]))
    }

    dat <- data.frame(size = size_vec)
    fit <- geomorph::procD.lm(coords_arr ~ size, data = dat, iter = N_PERM)

    coords_resid <- geomorph::arrayspecs(
      residuals(fit),
      p = dim(coords_arr)[1],
      k = dim(coords_arr)[2]
    )

    dimnames(coords_resid) <- dimnames(coords_arr)

    list(residuals = coords_resid, fit = fit)
  }
}

# ---------------------------
# Read and normalize papillae scores
# ---------------------------

if (!file.exists(PAPILLAE_CSV)) {
  stop("Could not find PAPILLAE_CSV: ", PAPILLAE_CSV,
       "\nSet PAPILLAE_CSV to the correct path.")
}

pap <- readr::read_csv(PAPILLAE_CSV, show_col_types = FALSE)
names(pap) <- tolower(gsub("[^A-Za-z0-9]+", "_", names(pap)))
names(pap) <- gsub("_+$", "", names(pap))

required_cols <- c("name", "papilla_presence", "papilla_confidence")
missing_cols <- setdiff(required_cols, names(pap))
if (length(missing_cols) > 0) {
  stop("Papillae CSV is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

normalize_presence <- function(x) {
  x0 <- stringr::str_to_lower(stringr::str_trim(as.character(x)))

  dplyr::case_when(
    x0 %in% c("present", "p", "ph", "pl", "papilla_present", "papilla present",
              "present high", "present low") ~ "present",
    x0 %in% c("absent", "a", "ah", "al", "not observed", "not_observed",
              "no papilla", "no_papilla", "papilla_absent", "papilla absent",
              "absent high", "absent low") ~ "absent",
    TRUE ~ NA_character_
  )
}

normalize_confidence <- function(x) {
  x0 <- stringr::str_to_lower(stringr::str_trim(as.character(x)))

  dplyr::case_when(
    x0 %in% c("high", "h", "high confidence", "high_confidence", "ph", "ah") ~ "high",
    x0 %in% c("low", "l", "low confidence", "low_confidence", "pl", "al",
              "not observed", "not_observed") ~ "low",
    TRUE ~ NA_character_
  )
}

pap <- pap %>%
  dplyr::transmute(
    specimen = as.character(name),
    papilla_presence_raw = papilla_presence,
    papilla_confidence_raw = papilla_confidence,
    papilla_presence = normalize_presence(papilla_presence),
    papilla_confidence = normalize_confidence(papilla_confidence)
  )

if (any(is.na(pap$papilla_presence))) {
  warning("Some papilla_presence values could not be normalized. They will be excluded.")
}

if (any(is.na(pap$papilla_confidence))) {
  warning("Some papilla_confidence values could not be normalized. They will be excluded where confidence is required.")
}

# Avoid duplicate rows silently changing sample sizes.
dupes <- pap %>% dplyr::count(specimen) %>% dplyr::filter(n > 1)
if (nrow(dupes) > 0) {
  stop("Duplicate specimen IDs in papillae CSV: ",
       paste(dupes$specimen, collapse = ", "))
}

# ---------------------------
# Merge metadata + papillae scores
# ---------------------------

gdf_pap <- gdf %>%
  dplyr::left_join(pap, by = "specimen") %>%
  dplyr::filter(!is.na(papilla_presence))

# Save join diagnostics
join_diagnostics <- list(
  n_coords_metadata = nrow(gdf),
  n_papillae_csv = nrow(pap),
  n_joined_with_scores = nrow(gdf_pap),
  papillae_not_in_gdf = list(setdiff(pap$specimen, gdf$specimen)),
  gdf_without_papillae_score = list(setdiff(gdf$specimen, pap$specimen))
)

readr::write_csv(
  tibble::tibble(
    field = names(join_diagnostics),
    value = purrr::map_chr(join_diagnostics, ~ paste(unlist(.x), collapse = "; "))
  ),
  file.path(OUTDIR, "join_diagnostics.csv")
)

readr::write_csv(gdf_pap, file.path(OUTDIR, "metadata_with_papillae_scores.csv"))

# ---------------------------
# Build analysis-specific scoring datasets
# ---------------------------

make_papilla_group <- function(df, analysis = c("conservative", "inclusive")) {
  analysis <- match.arg(analysis)

  out <- df

  if (analysis == "conservative") {
    out <- out %>%
      dplyr::filter(
        papilla_confidence == "high",
        papilla_presence %in% c("present", "absent")
      )
  }

  if (analysis == "inclusive") {
    out <- out %>%
      dplyr::filter(
        papilla_presence %in% c("present", "absent")
      )
  }

  out %>%
    dplyr::mutate(
      papilla_group = dplyr::case_when(
        papilla_presence == "present" ~ "papilla_present",
        papilla_presence == "absent"  ~ "no_papilla_observed",
        TRUE ~ NA_character_
      ),
      papilla_group = factor(
        papilla_group,
        levels = c("no_papilla_observed", "papilla_present")
      )
    ) %>%
    dplyr::filter(!is.na(papilla_group)) %>%
    droplevels()
}

# ---------------------------
# Define biological/figure subsets
# ---------------------------

make_dataset_subset <- function(df, dataset = c("all_data", "habitats_1950", "ct_time_series")) {
  dataset <- match.arg(dataset)

  if (dataset == "all_data") {
    return(df %>% dplyr::filter(!is.na(habitat)))
  }

  if (dataset == "habitats_1950") {
    return(df %>%
             dplyr::filter(!is.na(habitat), year == 1950))
  }

  if (dataset == "ct_time_series") {
    return(df %>%
             dplyr::filter(
               habitat == "Connecticut River",
               year %in% c(1950, 1956, 1970)
             ))
  }
}

# ---------------------------
# Plotting helpers
# ---------------------------

safe_hulls <- function(pc_df, x, y, group_col = "papilla_group") {
  pc_df %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::filter(dplyr::n() >= 3) %>%
    dplyr::slice(chull(.data[[x]], .data[[y]])) %>%
    dplyr::ungroup()
}

plot_pca_hull <- function(pc_df, pct_var, x = "PC1", y = "PC2",
                          title = NULL, subtitle = NULL) {
  hull_df <- safe_hulls(pc_df, x, y)

  x_lab <- paste0(x, " (", round(pct_var[x], 1), "%)")
  y_lab <- paste0(y, " (", round(pct_var[y], 1), "%)")

  ggplot(pc_df, aes(x = .data[[x]], y = .data[[y]], color = papilla_group)) +
    geom_polygon(
      data = hull_df,
      aes(fill = papilla_group, group = papilla_group),
      alpha = 0.12,
      color = NA,
      show.legend = FALSE
    ) +
    geom_path(
      data = hull_df,
      aes(group = papilla_group),
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    geom_point(size = 1.1, alpha = 0.85) +
    coord_equal() +
    labs(
      title = title,
      subtitle = subtitle,
      x = x_lab,
      y = y_lab,
      color = "Papilla score"
    ) +
    theme_classic(base_size = 7) +
    theme(
      plot.title = element_text(size = 7.5, face = "bold"),
      plot.subtitle = element_text(size = 6.5),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 6),
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6),
      legend.key.size = unit(0.35, "lines"),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25)
    )
}

run_shape_analysis <- function(coords_arr, meta, outdir, dataset_name, analysis_name, shape_type) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # Drop incomplete rows for size/metadata.
  meta <- meta %>%
    dplyr::filter(
      specimen %in% dimnames(coords_arr)[[3]],
      !is.na(papilla_group)
    ) %>%
    droplevels()

  coords_arr <- subset_coords_to_gdf(coords_arr, meta)

  # Basic group counts
  counts <- meta %>%
    dplyr::count(papilla_group, name = "n") %>%
    dplyr::mutate(
      analysis = analysis_name,
      dataset = dataset_name,
      shape_type = shape_type
    )

  readr::write_csv(counts, file.path(outdir, "group_counts.csv"))
  readr::write_csv(meta, file.path(outdir, "metadata_used.csv"))

  if (nrow(meta) < 4 || length(unique(meta$papilla_group)) < 2 || any(table(meta$papilla_group) < 2)) {
    msg <- paste0(
      "Skipping procD/PCA for ", analysis_name, " / ", dataset_name, " / ", shape_type,
      " because sample size is too small or one papilla group has fewer than 2 specimens.\n",
      "Group counts:\n",
      paste(capture.output(print(table(meta$papilla_group))), collapse = "\n")
    )
    writeLines(msg, file.path(outdir, "SKIPPED_too_few_specimens.txt"))
    message(msg)
    return(invisible(NULL))
  }

  # procD.lm: shape ~ papilla_group
  procD_fit <- geomorph::procD.lm(
    coords_arr ~ papilla_group,
    data = meta,
    iter = N_PERM
  )

  sink(file.path(outdir, "procD_shape_by_papilla_group.txt"))
  cat("Analysis:", analysis_name, "\n")
  cat("Dataset:", dataset_name, "\n")
  cat("Shape type:", shape_type, "\n\n")
  print(summary(procD_fit))
  sink()

  saveRDS(procD_fit, file.path(outdir, "procD_shape_by_papilla_group.rds"))

  # Pairwise test if both groups are present
  pair_fit <- RRPP::pairwise(
    procD_fit,
    groups = meta$papilla_group
  )

  sink(file.path(outdir, "pairwise_papilla_group.txt"))
  cat("Analysis:", analysis_name, "\n")
  cat("Dataset:", dataset_name, "\n")
  cat("Shape type:", shape_type, "\n\n")
  print(summary(pair_fit))
  sink()

  saveRDS(pair_fit, file.path(outdir, "pairwise_papilla_group.rds"))

  # PCA
  pca <- geomorph::gm.prcomp(coords_arr)
  pc_scores <- as.data.frame(pca$x) %>%
    tibble::rownames_to_column("specimen") %>%
    dplyr::left_join(meta, by = "specimen")

  eig <- pca$sdev^2
  pct_var_vec <- eig / sum(eig) * 100
  pct_var <- setNames(pct_var_vec, paste0("PC", seq_along(pct_var_vec)))

  readr::write_csv(pc_scores, file.path(outdir, "pca_scores.csv"))
  readr::write_csv(
    tibble::tibble(
      PC = names(pct_var),
      percent_variance = as.numeric(pct_var)
    ),
    file.path(outdir, "pca_percent_variance.csv")
  )

  saveRDS(pca, file.path(outdir, "gm_prcomp.rds"))

  # Save PCA hull plots
  for (panel in PC_PANELS) {
    x <- panel[1]
    y <- panel[2]

    if (!all(c(x, y) %in% names(pc_scores))) {
      next
    }

    p <- plot_pca_hull(
      pc_scores,
      pct_var = pct_var,
      x = x,
      y = y,
      title = paste(analysis_name, dataset_name, shape_type, sep = " | "),
      subtitle = paste0(x, " vs ", y)
    )

    fname_base <- paste0("PCA_hull_", x, "_", y)
    ggsave(
      filename = file.path(outdir, paste0(fname_base, ".pdf")),
      plot = p,
      width = SAVE_WIDTH,
      height = SAVE_HEIGHT,
      units = "in"
    )
    ggsave(
      filename = file.path(outdir, paste0(fname_base, ".png")),
      plot = p,
      width = SAVE_WIDTH,
      height = SAVE_HEIGHT,
      units = "in",
      dpi = SAVE_DPI
    )
  }

  invisible(list(
    procD = procD_fit,
    pairwise = pair_fit,
    pca = pca,
    pc_scores = pc_scores,
    counts = counts
  ))
}

run_one_dataset <- function(meta_scored, dataset_name, analysis_name) {
  meta_dataset <- make_dataset_subset(meta_scored, dataset_name)

  dataset_dir <- file.path(OUTDIR, analysis_name, dataset_name)
  dir.create(dataset_dir, recursive = TRUE, showWarnings = FALSE)

  readr::write_csv(meta_dataset, file.path(dataset_dir, "metadata_scored_subset.csv"))

  # Summary table by habitat/year/papilla group
  summary_counts <- meta_dataset %>%
    dplyr::count(habitat, year, papilla_group, papilla_confidence, name = "n") %>%
    dplyr::arrange(habitat, year, papilla_group, papilla_confidence)

  readr::write_csv(summary_counts, file.path(dataset_dir, "summary_counts_by_context.csv"))

  if (nrow(meta_dataset) == 0) {
    writeLines("No specimens after filtering.", file.path(dataset_dir, "SKIPPED_no_specimens.txt"))
    return(invisible(NULL))
  }

  # Raw shape
  coords_raw <- subset_coords_to_gdf(coords_gpa, meta_dataset)

  run_shape_analysis(
    coords_arr = coords_raw,
    meta = meta_dataset,
    outdir = file.path(dataset_dir, "raw_shape"),
    dataset_name = dataset_name,
    analysis_name = analysis_name,
    shape_type = "raw_shape"
  )

  # Size-corrected shape: residualize within the dataset being tested.
  # This keeps the size correction matched to the actual analysis subset.
  size_vec <- setNames(meta_dataset$size_for_allometry, meta_dataset$specimen)

  if (any(is.na(size_vec))) {
    writeLines(
      "Size-corrected analysis skipped because size_for_allometry contains NA.",
      file.path(dataset_dir, "size_corrected_SKIPPED_missing_size.txt")
    )
    return(invisible(NULL))
  }

  allo <- allometry_residuals(coords_raw, size_vec)
  coords_resid <- allo$residuals

  sink(file.path(dataset_dir, "allometry_model_shape_by_size.txt"))
  cat("Analysis:", analysis_name, "\n")
  cat("Dataset:", dataset_name, "\n")
  cat("Allometry model: shape ~ size_for_allometry\n")
  cat("Size variable:", unique(meta_dataset$size_label), "\n\n")
  print(summary(allo$fit))
  sink()

  saveRDS(allo$fit, file.path(dataset_dir, "allometry_model_shape_by_size.rds"))

  run_shape_analysis(
    coords_arr = coords_resid,
    meta = meta_dataset,
    outdir = file.path(dataset_dir, "size_corrected_shape"),
    dataset_name = dataset_name,
    analysis_name = analysis_name,
    shape_type = "size_corrected_shape"
  )

  invisible(NULL)
}

# ---------------------------
# Run all analyses
# ---------------------------

analysis_names <- c("conservative_high_confidence_only", "inclusive_sensitivity")
dataset_names <- c("all_data", "habitats_1950", "ct_time_series")

all_run_counts <- list()

for (analysis_name in analysis_names) {
  message("\n==============================")
  message("Running analysis: ", analysis_name)
  message("==============================")

  scoring_mode <- ifelse(
    analysis_name == "conservative_high_confidence_only",
    "conservative",
    "inclusive"
  )

  meta_scored <- make_papilla_group(gdf_pap, scoring_mode)

  analysis_dir <- file.path(OUTDIR, analysis_name)
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

  readr::write_csv(meta_scored, file.path(analysis_dir, "metadata_scored.csv"))

  analysis_counts <- meta_scored %>%
    dplyr::count(papilla_group, papilla_presence, papilla_confidence, name = "n")

  readr::write_csv(analysis_counts, file.path(analysis_dir, "analysis_group_counts.csv"))

  for (dataset_name in dataset_names) {
    message("  Dataset: ", dataset_name)
    run_one_dataset(meta_scored, dataset_name, analysis_name)

    tmp <- make_dataset_subset(meta_scored, dataset_name) %>%
      dplyr::count(papilla_group, name = "n") %>%
      dplyr::mutate(
        analysis = analysis_name,
        dataset = dataset_name,
        .before = 1
      )

    all_run_counts[[paste(analysis_name, dataset_name, sep = "_")]] <- tmp
  }
}

final_counts <- dplyr::bind_rows(all_run_counts)
readr::write_csv(final_counts, file.path(OUTDIR, "FINAL_group_counts_all_runs.csv"))

cat("\nDone. Outputs saved to: ", normalizePath(OUTDIR), "\n", sep = "")
cat("\nFinal group counts:\n")
print(final_counts)
