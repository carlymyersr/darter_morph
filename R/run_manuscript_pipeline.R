# ============================================================
# R/run_manuscript_pipeline.R
#
# Master runner organized by the current methods/results outline in:
#   ~/Documents/darter_figures/outline methods and results 12.19.16 PM.txt
#
# Usage:
#   Rscript R/run_manuscript_pipeline.R
#   Rscript R/run_manuscript_pipeline.R --check-only
#   Rscript R/run_manuscript_pipeline.R --stop-on-error
#   Rscript R/run_manuscript_pipeline.R --include-non-used
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

CHECK_ONLY <- "--check-only" %in% args
STOP_ON_ERROR <- "--stop-on-error" %in% args
INCLUDE_NON_USED <- "--include-non-used" %in% args

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
  file.path(path.expand("~"), "Documents", "darter_figures", "outline methods and results 12.19.16 PM.txt")
)
if (!file.exists(outline_file)) {
  candidates <- list.files(
    file.path(path.expand("~"), "Documents", "darter_figures"),
    pattern = "^outline methods and results 12[.]19[.]16.*PM[.]txt$",
    full.names = TRUE
  )
  if (length(candidates) > 0) outline_file <- candidates[1]
}

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

output_snapshot <- function(root = ".") {
  paths <- list_outputs(root)
  if (length(paths) == 0) {
    return(data.frame(path = character(), mtime = numeric(), stringsAsFactors = FALSE))
  }
  info <- file.info(paths)
  data.frame(
    path = paths,
    mtime = as.numeric(info$mtime),
    stringsAsFactors = FALSE
  )
}

changed_outputs <- function(before, after) {
  if (nrow(after) == 0) return(character())
  before_mtime <- setNames(before$mtime, before$path)
  is_new <- !after$path %in% before$path
  is_updated <- after$path %in% before$path &
    after$mtime > before_mtime[after$path] + 0.5
  after$path[is_new | is_updated]
}

classify_output <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("png", "pdf", "jpg", "jpeg", "tif", "tiff", "svg")) return("figures")
  if (ext %in% c("csv", "txt", "tsv", "xlsx")) return("tables")
  if (ext %in% c("rds", "rda", "rdata")) return("models")
  "other"
}

copy_section_outputs <- function(paths, section_dir, task_index, task_script, output_root) {
  if (length(paths) == 0) return(character())
  copied <- character()
  script_stub <- tools::file_path_sans_ext(basename(task_script))
  script_stub <- gsub("[^A-Za-z0-9._-]+", "_", script_stub)
  for (path in paths) {
    bucket <- classify_output(path)
    out_dir <- file.path(output_root, section_dir, bucket)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    dest_name <- paste0(sprintf("%02d", task_index), "_", script_stub, "__", basename(path))
    dest <- file.path(out_dir, dest_name)
    if (file.exists(dest)) {
      dest <- file.path(out_dir, paste0(sprintf("%02d", task_index), "_", script_stub, "__",
                                        as.integer(Sys.time()), "_", basename(path)))
    }
    ok <- file.copy(path, dest, overwrite = TRUE)
    if (ok) copied <- c(copied, normalizePath(dest, winslash = "/", mustWork = FALSE))
  }
  copied
}

reference_outputs <- function(root) {
  if (!dir.exists(root)) return(character())
  exts <- c("png", "pdf", "jpg", "jpeg", "csv", "txt", "rds")
  pattern <- paste0("\\.(", paste(exts, collapse = "|"), ")$")
  list.files(root, recursive = TRUE, full.names = TRUE, pattern = pattern, ignore.case = TRUE)
}

task <- function(script, phase, section_dir, outline, required = TRUE, non_used = FALSE) {
  data.frame(
    phase = phase,
    section_dir = section_dir,
    outline_section = outline,
    script = script,
    required = required,
    non_used = non_used,
    stringsAsFactors = FALSE
  )
}

tasks <- rbind(
  task("R/methods/02_landmark_acquisition_geometric_morphometrics/00_setup_morpho.R", "methods_02_gpa", "methods/02_landmark_acquisition_geometric_morphometrics", "Methods 2: landmark acquisition and geometric morphometrics"),
  task("R/methods/01_specimen_sampling_study_design/01_build_metadata.R", "methods_01_sampling", "methods/01_specimen_sampling_study_design", "Methods 1: specimen sampling and study design"),
  task("R/methods/01_specimen_sampling_study_design/02_subset_1950.R", "methods_01_sampling", "methods/01_specimen_sampling_study_design", "Methods 1: 1950 waterbody subsets"),
  task("R/methods/01_specimen_sampling_study_design/03_subset_CT_timeseries.R", "methods_01_sampling", "methods/01_specimen_sampling_study_design", "Methods 1: Connecticut River time-series subset"),
  task("R/methods/01_specimen_sampling_study_design/04_subset_CT_timeseries_plus_1950habitats.R", "methods_01_sampling", "methods/01_specimen_sampling_study_design", "Methods 1: combined watershed subset"),
  task("R/methods/01_specimen_sampling_study_design/05_subset_1950_quabbin_swift.R", "methods_01_sampling", "methods/01_specimen_sampling_study_design", "Methods 1: Quabbin/Swift subset"),

  task("R/methods/03_size_correction_allometry/02_size_boxplots_and_tests.R", "methods_03_size", "methods/03_size_correction_allometry", "Methods 3: size correction and allometry", required = FALSE),
  task("R/methods/03_size_correction_allometry/03_1950_raw_pca_size_colored.R", "methods_03_size", "methods/03_size_correction_allometry", "Methods 3: 1950 raw size-colored PCA", required = FALSE),
  task("R/methods/03_size_correction_allometry/04_ct_timeseries_raw_pca_size_colored.R", "methods_03_size", "methods/03_size_correction_allometry", "Methods 3: CT raw size-colored PCA", required = FALSE),
  task("R/methods/03_size_correction_allometry/01_combined_raw_pca_size_colored.R", "methods_03_size", "methods/03_size_correction_allometry", "Methods 3: combined raw size-colored PCA", required = FALSE),
  task("R/methods/03_size_correction_allometry/05_combined_allometry_procD.R", "methods_03_size", "methods/03_size_correction_allometry", "Methods 3: combined procD allometry and logCsize x group", required = FALSE),
  task("R/methods/04_genital_papillae_sensitivity/01_papillae_parallel_sensitivity.R", "methods_04_papillae", "methods/04_genital_papillae_sensitivity", "Methods 4: genital papillae sensitivity", required = FALSE),

  task("R/results/01_structure_variation_morphospace/01_1950_waterbody_pca_residual_pc12.R", "results_01_morphospace", "results/01_structure_variation_morphospace", "Results 1: 1950 waterbody PC1-PC2 residual PCA"),
  task("R/results/01_structure_variation_morphospace/02_ct_timeseries_pca_residual_pc12.R", "results_01_morphospace", "results/01_structure_variation_morphospace", "Results 1: CT time-series PC1-PC2 residual PCA"),
  task("R/results/01_structure_variation_morphospace/03_1950_tps_residual_pc12_extremes.R", "results_01_morphospace", "results/01_structure_variation_morphospace", "Results 1: 1950 TPS synthetic extremes"),
  task("R/results/01_structure_variation_morphospace/04_ct_timeseries_tps_residual_pc12_extremes.R", "results_01_morphospace", "results/01_structure_variation_morphospace", "Results 1: CT TPS synthetic extremes"),

  task("R/results/02_mean_shape_differentiation/01_mean_positions_and_tps_by_group.R", "results_02_mean_shape", "results/02_mean_shape_differentiation", "Results 2: mean positions and TPS plots"),
  task("R/results/02_mean_shape_differentiation/02_1950_waterbody_procD_pairwise.R", "results_02_mean_shape", "results/02_mean_shape_differentiation", "Results 2: 1950 procD and pairwise"),
  task("R/results/02_mean_shape_differentiation/03_ct_timeseries_procD_pairwise.R", "results_02_mean_shape", "results/02_mean_shape_differentiation", "Results 2: CT procD and pairwise"),

  task("R/results/03_trait_specific_patterns/01_linear_trait_measurements.R", "results_03_traits", "results/03_trait_specific_patterns", "Results 3: linear trait measurements", required = FALSE),
  task("R/results/03_trait_specific_patterns/02_curve_trait_measurements.R", "results_03_traits", "results/03_trait_specific_patterns", "Results 3: hyoid curve trait", required = FALSE),
  task("R/results/03_trait_specific_patterns/03_mouth_angle_measurements_and_tests.R", "results_03_traits", "results/03_trait_specific_patterns", "Results 3: mouth angle trait", required = FALSE),
  task("R/results/03_trait_specific_patterns/04_trait_boxplots_combined_groups.R", "results_03_traits", "results/03_trait_specific_patterns", "Results 3: combined trait boxplots", required = FALSE),
  task("R/results/03_trait_specific_patterns/05_mouth_angle_network_figure.R", "results_03_traits", "results/03_trait_specific_patterns", "Results 3: mouth angle supporting figure", required = FALSE),
  task("R/results/03_trait_specific_patterns/06_trait_faceted_summary_figure.R", "results_03_traits", "results/03_trait_specific_patterns", "Results 3: trait summary figure", required = FALSE),

  task("R/results/04_within_group_variation/01_disparity_and_dispersion_tests.R", "results_04_variation", "results/04_within_group_variation", "Results 4: distribution of within-group variation", required = FALSE),

  task("R/results/05_ct_reference_distribution/01_full_landscape_residual_pc12.R", "results_05_ct_reference", "results/05_ct_reference_distribution", "Results 5: shared empirical landscape PC1-PC2"),
  task("R/results/05_ct_reference_distribution/02_full_landscape_residual_pc23.R", "results_05_ct_reference", "results/05_ct_reference_distribution", "Results 5: shared empirical landscape PC2-PC3", required = FALSE),
  task("R/results/05_ct_reference_distribution/03_ct_reference_mahalanobis_rarefaction.R", "results_05_ct_reference", "results/05_ct_reference_distribution", "Results 5: CT reference Mahalanobis and rarefaction"),

  task("R/results/06_persistent_local_divergence/01_quabbin_swift_context_dependence_procD.R", "results_06_quabbin_swift", "results/06_persistent_local_divergence", "Results 6: persistent Swift/Quabbin local divergence", required = FALSE),
  task("R/results/07_hydrology_based_structure/01_hydrology_groups_pca_procD_mahalanobis.R", "results_07_hydrology", "results/07_hydrology_based_structure", "Results 7: hydrology-based structure"),
  task("R/results/08_reservoir_internal_structure/01_quabbin_swift_sampling_location_pca.R", "results_08_reservoir", "results/08_reservoir_internal_structure", "Results 8: internal reservoir structure", required = FALSE),

  task("R/results/09_modularity_integration/TestA_5modules.R", "results_09_modularity", "results/09_modularity_integration", "Results 9: 1950 five-module test", required = FALSE),
  task("R/results/09_modularity_integration/TestB_4modules.R", "results_09_modularity", "results/09_modularity_integration", "Results 9: 1950 four-module test", required = FALSE),
  task("R/results/09_modularity_integration/TestC_3modules.R", "results_09_modularity", "results/09_modularity_integration", "Results 9: 1950 three-module test", required = FALSE),
  task("R/results/09_modularity_integration/TestD_2modules_AP.R", "results_09_modularity", "results/09_modularity_integration", "Results 9: 1950 AP partition", required = FALSE),
  task("R/results/09_modularity_integration/TestE_2modules_DV.R", "results_09_modularity", "results/09_modularity_integration", "Results 9: 1950 DV partition", required = FALSE),
  task("R/results/09_modularity_integration/TestA_5modules_CT3.R", "results_09_modularity", "results/09_modularity_integration", "Results 9: CT five-module test", required = FALSE),
  task("R/results/09_modularity_integration/TestB_4modules_CT3.R", "results_09_modularity", "results/09_modularity_integration", "Results 9: CT four-module test", required = FALSE),
  task("R/results/09_modularity_integration/TestC_3modules_CT3.R", "results_09_modularity", "results/09_modularity_integration", "Results 9: CT three-module test", required = FALSE),
  task("R/results/09_modularity_integration/TestD_2modules_AP_CT3.R", "results_09_modularity", "results/09_modularity_integration", "Results 9: CT AP partition", required = FALSE),
  task("R/results/09_modularity_integration/TestE_2modules_DV_CT3.R", "results_09_modularity", "results/09_modularity_integration", "Results 9: CT DV partition", required = FALSE),
  task("R/results/09_modularity_integration/99_modularity_multiple_testing_correction.R", "results_09_modularity", "results/09_modularity_integration", "Results 9: Bonferroni/BH corrections", required = FALSE),

  task("R/supplemental/01_structure_variation_extra_axes/05_1950_waterbody_pca_residual_pc23_supplemental.R", "supplemental_01_morphospace", "supplemental/01_structure_variation_extra_axes", "Supplement: 1950 PC2-PC3", required = FALSE),
  task("R/supplemental/01_structure_variation_extra_axes/06_1950_waterbody_pca_residual_pc34_supplemental.R", "supplemental_01_morphospace", "supplemental/01_structure_variation_extra_axes", "Supplement: 1950 PC3-PC4", required = FALSE),
  task("R/supplemental/01_structure_variation_extra_axes/07_1950_tps_residual_pc23_extremes_supplemental.R", "supplemental_01_morphospace", "supplemental/01_structure_variation_extra_axes", "Supplement: 1950 TPS PC2-PC3", required = FALSE),
  task("R/supplemental/01_structure_variation_extra_axes/08_ct_timeseries_pca_residual_pc23_supplemental.R", "supplemental_01_morphospace", "supplemental/01_structure_variation_extra_axes", "Supplement: CT PC2-PC3", required = FALSE),
  task("R/supplemental/01_structure_variation_extra_axes/09_ct_timeseries_tps_residual_pc23_extremes_supplemental.R", "supplemental_01_morphospace", "supplemental/01_structure_variation_extra_axes", "Supplement: CT TPS PC2-PC3", required = FALSE),
  task("R/supplemental/06_persistent_local_divergence/swift_quabbin_1950_pca_clean.R", "supplemental_06_quabbin_swift", "supplemental/06_persistent_local_divergence", "Supplement: Swift/Quabbin clean PCA", required = FALSE),
  task("R/supplemental/07_hydrology_alternatives/1950_hydrology_groups_pca_procD.R", "supplemental_07_hydrology", "supplemental/07_hydrology_alternatives", "Supplement: 1950-only hydrology grouping", required = FALSE),
  task("R/supplemental/07_hydrology_alternatives/1950_variation_source_pca.R", "supplemental_07_hydrology", "supplemental/07_hydrology_alternatives", "Supplement: 1950 variation source PCA", required = FALSE),
  task("R/supplemental/07_hydrology_alternatives/1950_variation_source_pc23.R", "supplemental_07_hydrology", "supplemental/07_hydrology_alternatives", "Supplement: 1950 variation source PC2-PC3", required = FALSE),
  task("R/supplemental/07_hydrology_alternatives/1950_variation_source_procD.R", "supplemental_07_hydrology", "supplemental/07_hydrology_alternatives", "Supplement: 1950 variation source procD", required = FALSE),
  task("R/supplemental/07_hydrology_alternatives/1950_large_vs_small_waterbody_pca.R", "supplemental_07_hydrology", "supplemental/07_hydrology_alternatives", "Supplement: large/small waterbody PCA", required = FALSE),
  task("R/supplemental/07_hydrology_alternatives/1950_large_vs_small_waterbody_procD.R", "supplemental_07_hydrology", "supplemental/07_hydrology_alternatives", "Supplement: large/small waterbody procD", required = FALSE),
  task("R/supplemental/07_hydrology_alternatives/hydrology_groups_colored_by_mainstem_year.R", "supplemental_07_hydrology", "supplemental/07_hydrology_alternatives", "Supplement: hydrology groups colored by CT year", required = FALSE),
  task("R/supplemental/07_hydrology_alternatives/hydrology_groups_sicb_short.R", "supplemental_07_hydrology", "supplemental/07_hydrology_alternatives", "Supplement: hydrology groups SICB version", required = FALSE),
  task("R/supplemental/03_size_correction_allometry/1950_raw_tps_synthetic.R", "supplemental_03_size", "supplemental/03_size_correction_allometry", "Supplement: raw TPS size check", required = FALSE),
  task("R/supplemental/03_size_correction_allometry/1950_raw_and_residual_stats_legacy.R", "supplemental_03_size", "supplemental/03_size_correction_allometry", "Supplement: legacy 1950 stats", required = FALSE)
)

prepare_input_links(data_root)

MASTER_RUN_ID <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_dir <- file.path("Outputs", "manuscript_runs", MASTER_RUN_ID)
section_output_root <- file.path("manuscript_outputs", MASTER_RUN_ID)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
if (!CHECK_ONLY) {
  dir.create(section_output_root, recursive = TRUE, showWarnings = FALSE)
}

message("Repository root: ", normalizePath(repo_root, winslash = "/"))
message("Input data root: ", normalizePath(data_root, winslash = "/", mustWork = FALSE))
message("Reference output root: ", normalizePath(reference_root, winslash = "/", mustWork = FALSE))
message("Outline file: ", normalizePath(outline_file, winslash = "/", mustWork = FALSE))
message("Manuscript run log: ", normalizePath(log_dir, winslash = "/"))
if (!CHECK_ONLY) {
  message("Sectioned run outputs: ", normalizePath(section_output_root, winslash = "/", mustWork = FALSE))
}

all_repo_scripts <- list.files("R", recursive = TRUE, full.names = TRUE, pattern = "\\.[Rr]$")
if (!INCLUDE_NON_USED) {
  all_repo_scripts <- all_repo_scripts[!grepl("^R/non_used_scripts/", all_repo_scripts)]
}
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
  "R/run_manuscript_pipeline.R"
)
allowed_unplanned <- normalizePath(allowed_unplanned, winslash = "/", mustWork = FALSE)
unplanned <- setdiff(normalizePath(all_repo_scripts, winslash = "/", mustWork = FALSE),
                     c(planned, allowed_unplanned))
unplanned <- unplanned[!grepl("/R/non_used_scripts/", unplanned, fixed = TRUE)]
unplanned <- unplanned[!grepl("/R/results/09_modularity_integration/_run_", unplanned, fixed = TRUE)]
if (length(unplanned) > 0) {
  writeLines(unplanned, file.path(log_dir, "unplanned_scripts.txt"))
  warning("Some scripts are not in the manuscript execution plan. See unplanned_scripts.txt")
}

if (!INCLUDE_NON_USED) {
  tasks <- tasks[!tasks$non_used, , drop = FALSE]
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
  section_output_files = integer(),
  stringsAsFactors = FALSE
)

before_outputs <- output_snapshot(".")
all_copied_outputs <- character()

for (i in seq_len(nrow(tasks))) {
  row <- tasks[i, ]
  message("\n[", i, "/", nrow(tasks), "] ", row$phase, " :: ", row$script)
  before_task_outputs <- output_snapshot(".")
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
  after_task_outputs <- output_snapshot(".")
  task_outputs <- changed_outputs(before_task_outputs, after_task_outputs)
  copied_outputs <- copy_section_outputs(
    task_outputs,
    row$section_dir,
    i,
    row$script,
    section_output_root
  )
  all_copied_outputs <- c(all_copied_outputs, copied_outputs)
  results <- rbind(results, data.frame(
    phase = row$phase,
    section_dir = row$section_dir,
    outline_section = row$outline_section,
    script = row$script,
    required = row$required,
    status = status,
    elapsed_sec = elapsed,
    message = err %||% "",
    section_output_files = length(copied_outputs),
    stringsAsFactors = FALSE
  ))
  write.csv(results, file.path(log_dir, "manuscript_results.csv"), row.names = FALSE)
  if (status == "failed_required" && STOP_ON_ERROR) {
    stop("Required task failed: ", row$script, "\n", err)
  }
}

after_outputs <- output_snapshot(".")
created_outputs <- changed_outputs(before_outputs, after_outputs)
ref_outputs <- reference_outputs(reference_root)
comparison <- data.frame(
  generated_file = created_outputs,
  generated_basename = basename(created_outputs),
  basename_in_reference = basename(created_outputs) %in% basename(ref_outputs),
  stringsAsFactors = FALSE
)
write.csv(comparison, file.path(log_dir, "generated_vs_reference_basenames.csv"), row.names = FALSE)
write.csv(
  data.frame(sectioned_output_file = all_copied_outputs, stringsAsFactors = FALSE),
  file.path(log_dir, "sectioned_output_manifest.csv"),
  row.names = FALSE
)

message("\nManuscript pipeline complete.")
message("Tasks: ", nrow(results))
message("Successful: ", sum(results$status == "ok"))
message("Required failures: ", sum(results$status == "failed_required"))
message("Optional failures: ", sum(results$status == "failed_optional"))
message("Generated outputs: ", length(created_outputs))
message("Copied into sectioned run outputs: ", length(all_copied_outputs))
message("Generated output basenames found in reference: ", sum(comparison$basename_in_reference))
message("Log directory: ", normalizePath(log_dir, winslash = "/"))
message("Sectioned output directory: ", normalizePath(section_output_root, winslash = "/"))

if (any(results$status == "failed_required")) quit(status = 1)
