# ============================================================
# R/run_manuscript_pipeline.R
#
# Master runner organized by the methods/results outline in:
#   ~/Documents/darter_figures/outline methods and results.txt
#
# Usage:
#   Rscript R/run_manuscript_pipeline.R
#   Rscript R/run_manuscript_pipeline.R --check-only
#   Rscript R/run_manuscript_pipeline.R --stop-on-error
#   Rscript R/run_manuscript_pipeline.R --include-archive
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

CHECK_ONLY <- "--check-only" %in% args
STOP_ON_ERROR <- "--stop-on-error" %in% args
INCLUDE_ARCHIVE <- "--include-archive" %in% args

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x)) y else x

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "R/run_manuscript_pipeline.R")
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "R"))) repo_root <- normalizePath(getwd(), mustWork = TRUE)
setwd(repo_root)

reference_root <- Sys.getenv(
  "DARTER_FIGURES_REFERENCE",
  file.path(path.expand("~"), "Documents", "darter_figures")
)

data_root <- Sys.getenv(
  "DARTER_DATA_ROOT",
  file.path(path.expand("~"), "Documents", "darter_morphometrics")
)

outline_file <- Sys.getenv(
  "DARTER_OUTLINE",
  file.path(path.expand("~"), "Documents", "darter_figures", "outline methods and results.txt")
)

prepare_input_links <- function(data_root) {
  required <- c("darter_curves.txt", "landmarks_ref.txt", "side_shapes")
  missing_required <- file.path(data_root, required)[!file.exists(file.path(data_root, required))]
  if (length(missing_required) > 0) {
    stop("Required input data are missing from DARTER_DATA_ROOT: ",
         paste(missing_required, collapse = ", "))
  }

  link_names <- c("darter_curves.txt", "landmarks_ref.txt", "side_shapes", "photos",
                  "trait_measurements", "papillae.csv")
  for (name in link_names) {
    target <- file.path(data_root, name)
    if (!file.exists(target) || file.exists(name)) next
    ok <- file.symlink(target, name)
    if (!ok && !file.exists(name)) {
      warning("Could not create local input link: ", name, " -> ", target)
    }
  }
}

literal_sources <- function(path) {
  lines <- readLines(path, warn = FALSE)
  source_lines <- lines[grepl("^\\s*(if\\s*\\([^)]*\\)\\s*)?source\\([\"']", lines)]
  if (length(source_lines) == 0) return(character())
  sub(".*source\\([\"']([^\"']+)[\"'].*", "\\1", source_lines)
}

check_source_paths <- function(paths) {
  out <- do.call(rbind, lapply(paths, function(path) {
    srcs <- literal_sources(path)
    if (length(srcs) == 0) return(NULL)
    data.frame(
      script = path,
      source_path = srcs,
      exists = file.exists(srcs),
      stringsAsFactors = FALSE
    )
  }))
  if (is.null(out)) {
    out <- data.frame(script = character(), source_path = character(), exists = logical())
  }
  out
}

list_outputs <- function(root = ".") {
  dirs <- c("Figures", "figures", "Outputs", "trait_measurements", "papillae",
            "outputs_mouth_to_body_angle")
  dirs <- dirs[dir.exists(file.path(root, dirs))]
  if (length(dirs) == 0) return(character())
  normalizePath(
    unlist(lapply(dirs, function(d) {
      list.files(file.path(root, d), recursive = TRUE, full.names = TRUE, all.files = FALSE)
    })),
    winslash = "/",
    mustWork = FALSE
  )
}

reference_outputs <- function(root) {
  if (!dir.exists(root)) return(character())
  exts <- c("png", "pdf", "jpg", "jpeg", "csv", "txt", "rds")
  pattern <- paste0("\\.(", paste(exts, collapse = "|"), ")$")
  list.files(root, recursive = TRUE, full.names = TRUE, pattern = pattern, ignore.case = TRUE)
}

task <- function(script, phase, outline, required = TRUE, archive = FALSE) {
  data.frame(
    phase = phase,
    outline_section = outline,
    script = script,
    required = required,
    archive = archive,
    stringsAsFactors = FALSE
  )
}

tasks <- rbind(
  task("R/00_setup_morpho.R", "00_core", "Methods 2: landmark acquisition and GPA"),
  task("R/01_build_metadata.R", "00_core", "Methods 1: sampling and study design"),
  task("R/02_subset_1950.R", "00_core", "Methods 1: 1950 waterbody subsets"),
  task("R/03_subset_CT_timeseries.R", "00_core", "Methods 1: CT time-series subset"),
  task("R/04_subset_CT_timeseries_plus_1950habitats.R", "00_core", "Methods 1: combined landscape subset"),
  task("R/05_subset_1950_quabbin_swift.R", "00_core", "Methods 1: Quabbin/Swift subset"),

  task("R/01_methods/01_size_allometry/02_size_boxplots_and_tests.R", "01_methods_size", "Methods 3: size correction and allometry", required = FALSE),
  task("R/01_methods/01_size_allometry/03_1950_raw_pca_size_colored.R", "01_methods_size", "Methods 3: 1950 raw size-colored PCA", required = FALSE),
  task("R/01_methods/01_size_allometry/04_ct_timeseries_raw_pca_size_colored.R", "01_methods_size", "Methods 3: CT raw size-colored PCA", required = FALSE),
  task("R/01_methods/01_size_allometry/01_combined_raw_pca_size_colored.R", "01_methods_size", "Methods 3: combined raw size-colored PCA", required = FALSE),
  task("R/01_methods/02_papillae_sensitivity/01_papillae_parallel_sensitivity.R", "01_methods_papillae", "Methods 4: papillae sensitivity", required = FALSE),

  task("R/02_results/01_morphospace_structure/01_1950_waterbody_pca_residual_pc12.R", "02_results_01_morphospace", "Results 1: 1950 PC1-PC2 residual PCA"),
  task("R/02_results/01_morphospace_structure/02_ct_timeseries_pca_residual_pc12.R", "02_results_01_morphospace", "Results 1: CT PC1-PC2 residual PCA"),
  task("R/02_results/01_morphospace_structure/03_1950_tps_residual_pc12_extremes.R", "02_results_01_morphospace", "Results 1: 1950 TPS synthetic extremes"),
  task("R/02_results/01_morphospace_structure/04_ct_timeseries_tps_residual_pc12_extremes.R", "02_results_01_morphospace", "Results 1: CT TPS synthetic extremes"),

  task("R/02_results/02_mean_shape_differentiation/01_mean_positions_and_tps_by_group.R", "02_results_02_mean_shape", "Results 2: mean positions and TPS plots"),
  task("R/02_results/02_mean_shape_differentiation/02_1950_waterbody_procD_pairwise.R", "02_results_02_mean_shape", "Results 2: 1950 procD and pairwise"),
  task("R/02_results/02_mean_shape_differentiation/03_ct_timeseries_procD_pairwise.R", "02_results_02_mean_shape", "Results 2: CT procD and pairwise"),

  task("R/02_results/03_trait_patterns/01_linear_trait_measurements.R", "02_results_03_traits", "Results 3: linear trait measurements", required = FALSE),
  task("R/02_results/03_trait_patterns/02_curve_trait_measurements.R", "02_results_03_traits", "Results 3: hyoid curve trait", required = FALSE),
  task("R/02_results/03_trait_patterns/03_mouth_angle_measurements_and_tests.R", "02_results_03_traits", "Results 3: mouth angle trait", required = FALSE),
  task("R/02_results/03_trait_patterns/04_trait_boxplots_combined_groups.R", "02_results_03_traits", "Results 3: combined trait boxplots", required = FALSE),
  task("R/02_results/03_trait_patterns/05_mouth_angle_network_figure.R", "02_results_03_traits", "Results 3: mouth angle supporting figure", required = FALSE),
  task("R/02_results/03_trait_patterns/06_trait_faceted_summary_figure.R", "02_results_03_traits", "Results 3: trait summary figure", required = FALSE),

  task("R/02_results/04_within_group_variation/01_disparity_and_dispersion_tests.R", "02_results_04_variation", "Results 4: disparity and dispersion", required = FALSE),

  task("R/02_results/05_ct_reference_distribution/01_full_landscape_residual_pc12.R", "02_results_05_ct_reference", "Results 5: full landscape PC1-PC2"),
  task("R/02_results/05_ct_reference_distribution/02_full_landscape_residual_pc23.R", "02_results_05_ct_reference", "Results 5: full landscape PC2-PC3", required = FALSE),

  task("R/02_results/06_persistent_local_divergence/01_quabbin_swift_context_dependence_procD.R", "02_results_06_quabbin_swift", "Results 6: Quabbin/Swift context-dependence", required = FALSE),
  task("R/02_results/07_reservoir_internal_structure/01_quabbin_swift_sampling_location_pca.R", "02_results_07_reservoir", "Results 7: reservoir internal structure", required = FALSE),

  task("R/02_results/08_hydrology_structure/01_hydrology_groups_pca_procD_mahalanobis.R", "02_results_08_hydrology", "Results 8: hydrology groups, procD, Mahalanobis"),

  task("R/02_results/09_modularity_integration/TestA_5modules.R", "02_results_09_modularity", "Results 9: 1950 five-module test", required = FALSE),
  task("R/02_results/09_modularity_integration/TestB_4modules.R", "02_results_09_modularity", "Results 9: 1950 four-module test", required = FALSE),
  task("R/02_results/09_modularity_integration/TestC_3modules.R", "02_results_09_modularity", "Results 9: 1950 three-module test", required = FALSE),
  task("R/02_results/09_modularity_integration/TestD_2modules_AP.R", "02_results_09_modularity", "Results 9: 1950 AP partition", required = FALSE),
  task("R/02_results/09_modularity_integration/TestE_2modules_DV.R", "02_results_09_modularity", "Results 9: 1950 DV partition", required = FALSE),
  task("R/02_results/09_modularity_integration/TestA_5modules_CT3.R", "02_results_09_modularity", "Results 9: CT five-module test", required = FALSE),
  task("R/02_results/09_modularity_integration/TestB_4modules_CT3.R", "02_results_09_modularity", "Results 9: CT four-module test", required = FALSE),
  task("R/02_results/09_modularity_integration/TestC_3modules_CT3.R", "02_results_09_modularity", "Results 9: CT three-module test", required = FALSE),
  task("R/02_results/09_modularity_integration/TestD_2modules_AP_CT3.R", "02_results_09_modularity", "Results 9: CT AP partition", required = FALSE),
  task("R/02_results/09_modularity_integration/TestE_2modules_DV_CT3.R", "02_results_09_modularity", "Results 9: CT DV partition", required = FALSE),
  task("R/02_results/09_modularity_integration/99_multiple_testing_correction_optional.R", "02_results_09_modularity", "Results 9: multiple-testing correction after manual clean-output curation", required = FALSE, archive = TRUE),

  task("R/02_results/01_morphospace_structure/05_1950_waterbody_pca_residual_pc23_supplemental.R", "03_supplemental", "Supplement: 1950 PC2-PC3", required = FALSE),
  task("R/02_results/01_morphospace_structure/06_1950_waterbody_pca_residual_pc34_supplemental.R", "03_supplemental", "Supplement: 1950 PC3-PC4", required = FALSE),
  task("R/02_results/01_morphospace_structure/07_1950_tps_residual_pc23_extremes_supplemental.R", "03_supplemental", "Supplement: 1950 TPS PC2-PC3", required = FALSE),
  task("R/02_results/01_morphospace_structure/08_ct_timeseries_pca_residual_pc23_supplemental.R", "03_supplemental", "Supplement: CT PC2-PC3", required = FALSE),
  task("R/02_results/01_morphospace_structure/09_ct_timeseries_tps_residual_pc23_extremes_supplemental.R", "03_supplemental", "Supplement: CT TPS PC2-PC3", required = FALSE),
  task("R/03_supplemental/quabbin_swift/swift_quabbin_1950_pca_clean.R", "03_supplemental", "Supplement: Quabbin/Swift clean PCA", required = FALSE),
  task("R/03_supplemental/landscape_alternatives/1950_hydrology_groups_pca_procD.R", "03_supplemental", "Supplement: 1950-only hydrology grouping", required = FALSE),
  task("R/03_supplemental/landscape_alternatives/1950_variation_source_pca.R", "03_supplemental", "Supplement: 1950 variation source PCA", required = FALSE),
  task("R/03_supplemental/landscape_alternatives/1950_variation_source_pc23.R", "03_supplemental", "Supplement: 1950 variation source PC2-PC3", required = FALSE),
  task("R/03_supplemental/landscape_alternatives/1950_variation_source_procD.R", "03_supplemental", "Supplement: 1950 variation source procD", required = FALSE),
  task("R/03_supplemental/landscape_alternatives/1950_large_vs_small_waterbody_pca.R", "03_supplemental", "Supplement: large/small waterbody PCA", required = FALSE),
  task("R/03_supplemental/landscape_alternatives/1950_large_vs_small_waterbody_procD.R", "03_supplemental", "Supplement: large/small waterbody procD", required = FALSE),
  task("R/03_supplemental/landscape_alternatives/hydrology_groups_colored_by_mainstem_year.R", "03_supplemental", "Supplement: hydrology groups colored by CT year", required = FALSE),
  task("R/03_supplemental/landscape_alternatives/hydrology_groups_sicb_short.R", "03_supplemental", "Supplement: hydrology groups SICB version", required = FALSE),
  task("R/03_supplemental/size_allometry/1950_raw_tps_synthetic.R", "03_supplemental", "Supplement: raw TPS size check", required = FALSE),
  task("R/03_supplemental/size_allometry/1950_raw_and_residual_stats_legacy.R", "03_supplemental", "Supplement: legacy 1950 stats", required = FALSE),
  task("R/03_supplemental/trait_pc_exploratory/shape_change_1950_trait_pc_correlations.R", "03_supplemental", "Not included: exploratory trait-PC correlation", required = FALSE, archive = TRUE),
  task("R/03_supplemental/trait_pc_exploratory/shape_change_full_landscape_trait_pc_correlations.R", "03_supplemental", "Not included: exploratory trait-PC correlation", required = FALSE, archive = TRUE)
)

prepare_input_links(data_root)

MASTER_RUN_ID <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_dir <- file.path("Outputs", "manuscript_runs", MASTER_RUN_ID)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

message("Repository root: ", normalizePath(repo_root, winslash = "/"))
message("Input data root: ", normalizePath(data_root, winslash = "/", mustWork = FALSE))
message("Reference output root: ", normalizePath(reference_root, winslash = "/", mustWork = FALSE))
message("Outline file: ", normalizePath(outline_file, winslash = "/", mustWork = FALSE))
message("Manuscript run log: ", normalizePath(log_dir, winslash = "/"))

all_repo_scripts <- list.files("R", recursive = TRUE, full.names = TRUE, pattern = "\\.[Rr]$")
source_check <- check_source_paths(all_repo_scripts)
write.csv(source_check, file.path(log_dir, "source_path_check.csv"), row.names = FALSE)
missing_sources <- source_check[!source_check$exists, , drop = FALSE]
if (nrow(missing_sources) > 0) {
  write.csv(missing_sources, file.path(log_dir, "missing_source_paths.csv"), row.names = FALSE)
  stop("Missing source paths detected. See ", file.path(log_dir, "missing_source_paths.csv"))
}

missing_task_files <- tasks$script[!file.exists(tasks$script)]
if (length(missing_task_files) > 0) {
  stop("Master task files missing: ", paste(missing_task_files, collapse = ", "))
}

planned <- normalizePath(tasks$script, winslash = "/", mustWork = FALSE)
allowed_unplanned <- c(
  "R/run_manuscript_pipeline.R",
  "R/master_reproduce_figures.R",
  "R/00_GPA.R",
  "R/01_metadata.R",
  "R/helpers_angles.R",
  "R/01_helpers_angles.R",
  "R/01_build_angle_measurements.R"
)
allowed_unplanned <- normalizePath(allowed_unplanned, winslash = "/", mustWork = FALSE)
unplanned <- setdiff(normalizePath(all_repo_scripts, winslash = "/", mustWork = FALSE),
                     c(planned, allowed_unplanned))
unplanned <- unplanned[!grepl("/R/04_archive/", unplanned, fixed = TRUE)]
unplanned <- unplanned[!grepl("/R/02_results/09_modularity_integration/_run_", unplanned, fixed = TRUE)]
if (length(unplanned) > 0) {
  writeLines(unplanned, file.path(log_dir, "unplanned_scripts.txt"))
  warning("Some scripts are not in the manuscript execution plan. See unplanned_scripts.txt")
}

if (!INCLUDE_ARCHIVE) {
  tasks <- tasks[!tasks$archive, , drop = FALSE]
}

write.csv(tasks, file.path(log_dir, "manuscript_task_plan.csv"), row.names = FALSE)

if (CHECK_ONLY) {
  message("Check-only mode complete. Task plan written.")
  quit(status = 0)
}

results <- data.frame(
  phase = character(),
  outline_section = character(),
  script = character(),
  required = logical(),
  status = character(),
  elapsed_sec = numeric(),
  message = character(),
  stringsAsFactors = FALSE
)

before_outputs <- list_outputs(".")

for (i in seq_len(nrow(tasks))) {
  row <- tasks[i, ]
  message("\n[", i, "/", nrow(tasks), "] ", row$phase, " :: ", row$script)
  start <- proc.time()[["elapsed"]]
  err <- NULL
  status <- "ok"
  tryCatch(
    source(row$script, local = .GlobalEnv),
    error = function(e) {
      err <<- conditionMessage(e)
      status <<- if (isTRUE(row$required)) "failed_required" else "failed_optional"
    }
  )
  elapsed <- proc.time()[["elapsed"]] - start
  results <- rbind(results, data.frame(
    phase = row$phase,
    outline_section = row$outline_section,
    script = row$script,
    required = row$required,
    status = status,
    elapsed_sec = elapsed,
    message = err %||% "",
    stringsAsFactors = FALSE
  ))
  write.csv(results, file.path(log_dir, "manuscript_results.csv"), row.names = FALSE)
  if (status == "failed_required" && STOP_ON_ERROR) {
    stop("Required task failed: ", row$script, "\n", err)
  }
}

after_outputs <- list_outputs(".")
created_outputs <- setdiff(after_outputs, before_outputs)
ref_outputs <- reference_outputs(reference_root)
comparison <- data.frame(
  generated_file = created_outputs,
  generated_basename = basename(created_outputs),
  basename_in_reference = basename(created_outputs) %in% basename(ref_outputs),
  stringsAsFactors = FALSE
)
write.csv(comparison, file.path(log_dir, "generated_vs_reference_basenames.csv"), row.names = FALSE)

message("\nManuscript pipeline complete.")
message("Tasks: ", nrow(results))
message("Successful: ", sum(results$status == "ok"))
message("Required failures: ", sum(results$status == "failed_required"))
message("Optional failures: ", sum(results$status == "failed_optional"))
message("Generated outputs: ", length(created_outputs))
message("Generated output basenames found in reference: ", sum(comparison$basename_in_reference))
message("Log directory: ", normalizePath(log_dir, winslash = "/"))

if (any(results$status == "failed_required")) quit(status = 1)
