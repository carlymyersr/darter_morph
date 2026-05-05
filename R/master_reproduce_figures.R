# ============================================================
# R/master_reproduce_figures.R
#
# Master runner for the darter morphometrics repository.
# It sources the canonical pipeline, runs figure/stat scripts in sequence,
# and compares generated outputs against the existing Documents/darter_figures
# outputs by filename.
#
# Usage:
#   Rscript R/master_reproduce_figures.R
#   Rscript R/master_reproduce_figures.R --check-only
#   Rscript R/master_reproduce_figures.R --stop-on-error
#   Rscript R/master_reproduce_figures.R --include-interactive
# ============================================================

args <- commandArgs(trailingOnly = TRUE)

CHECK_ONLY <- "--check-only" %in% args
STOP_ON_ERROR <- "--stop-on-error" %in% args
INCLUDE_INTERACTIVE <- "--include-interactive" %in% args

`%||%` <- function(x, y) if (is.null(x) || is.na(x)) y else x

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "R/master_reproduce_figures.R")
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = FALSE)
if (!dir.exists(file.path(repo_root, "R"))) {
  repo_root <- normalizePath(getwd(), mustWork = TRUE)
}
setwd(repo_root)

reference_root <- Sys.getenv(
  "DARTER_FIGURES_REFERENCE",
  file.path(path.expand("~"), "Documents", "darter_figures")
)

data_root <- Sys.getenv(
  "DARTER_DATA_ROOT",
  file.path(path.expand("~"), "Documents", "darter_morphometrics")
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

prepare_input_links(data_root)

MASTER_RUN_ID <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_dir <- file.path("Outputs", "master_runs", MASTER_RUN_ID)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

message("Repository root: ", normalizePath(repo_root, winslash = "/"))
message("Input data root: ", normalizePath(data_root, winslash = "/", mustWork = FALSE))
message("Reference output root: ", normalizePath(reference_root, winslash = "/", mustWork = FALSE))
message("Master run log: ", normalizePath(log_dir, winslash = "/"))

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
  dirs <- c("Figures", "figures", "Outputs", "trait_measurements")
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

task <- function(script, phase, interactive = FALSE, required = TRUE) {
  data.frame(
    phase = phase,
    script = script,
    interactive = interactive,
    required = required,
    stringsAsFactors = FALSE
  )
}

tasks <- rbind(
  task("R/00_setup_morpho.R", "core"),
  task("R/01_build_metadata.R", "core"),
  task("R/02_subset_1950.R", "core"),
  task("R/03_subset_CT_timeseries.R", "core"),
  task("R/04_subset_CT_timeseries_plus_1950habitats.R", "core"),
  task("R/05_subset_1950_quabbin_swift.R", "core"),

  task("R/figures/first_scale_1950/1950_PCA_residual_Fig3.R", "main_figures"),
  task("R/figures/first_scale_1950/1950_Locality_Resid_PC3.R", "main_figures"),
  task("R/figures/first_scale_1950/1950_Locality_Resid_PC4.R", "main_figures"),
  task("R/figures/first_scale_1950/1950_TPS_residual_synthetic.R", "main_figures"),
  task("R/figures/first_scale_1950/1950_TPS_residual_PC23_extremes.R", "main_figures"),
  task("R/figures/first_scale_CT/Figure5_PCA_CT_timeseries_RESID_logCsize.R", "main_figures"),
  task("R/figures/first_scale_CT/Figure5_PCA_CT_timeseries_RESID_logCsize_PC2_PC3.R", "main_figures"),
  task("R/figures/first_scale_CT/TPS_CT3_residual_synthetic.R", "main_figures"),
  task("R/figures/first_scale_CT/CT3_TPS_residual_PC23_extremes.R", "main_figures"),
  task("R/figures/sicb/1950_figures_SICB.R", "main_figures"),
  task("R/figures/sicb/CT3_figures_SICB.R", "main_figures"),
  task("R/figures/mean_shapes/shared_PCA_MEANS_ONLY_SICB_colored.R", "main_figures"),
  task("R/figures/landscape/Figure6_PCA_CTtimeseries_plus_1950habitats_RAW.R", "main_figures"),
  task("R/figures/landscape/Figure7_RESID_CTtimeseries_plus_1950habitats.R", "main_figures"),
  task("R/figures/landscape/Figure7_RESID_CTtimeseries_plus_1950habitats_PC23.R", "main_figures"),
  task("R/figures/landscape/watershed_landscape_SICB.R", "main_figures"),
  task("R/figures/landscape/TPS_SICB_Landscape.R", "main_figures"),

  task("R/traits/mouth_angle/digitize_trait_landmarks.R", "trait_build", interactive = TRUE, required = FALSE),
  task("R/traits/mouth_angle/build_trait_measurements_signed.R", "trait_build", required = FALSE),
  task("R/traits/mouth_angle/mouth_angle_trait_analysis_signed.R", "trait_figures", required = FALSE),
  task("R/traits/mouth_angle/SICB_signed_mouth_angle_networks.R", "trait_figures", required = FALSE),
  task("R/traits/mouth_angle/facet_mouth_angle_SICB_figure.R", "trait_figures", required = FALSE),
  task("R/traits/mouth_angle/mouth_angle_3_point_analysis.R", "trait_figures", required = FALSE),
  task("R/traits/distances_curves/landmark_distance_measurements.R", "trait_figures", required = FALSE),
  task("R/traits/distances_curves/curve_shape_metrics_raw_then_sizecorrected.R", "trait_figures", required = FALSE),
  task("R/traits/distances_curves/SICB_curves_distances_figures.R", "trait_figures", required = FALSE),

  task("R/figures/first_scale_1950/1950_stats_models.R", "stats"),
  task("R/figures/first_scale_CT/CT3_stats_models.R", "stats"),
  task("R/stats/disparity_variance_tests_with_stats_letters.R", "stats", required = FALSE),
  task("R/stats/disparity_variance_tests.R", "stats", required = FALSE),
  task("R/supplemental/size/size_analysis_SICB.R", "stats", required = FALSE),
  task("R/supplemental/sex_effects/papillae_parallel_analyses.R", "stats", required = FALSE),

  task("R/modularity/tests_1950/TestA_5modules.R", "modularity", required = FALSE),
  task("R/modularity/tests_1950/TestB_4modules.R", "modularity", required = FALSE),
  task("R/modularity/tests_1950/TestC_3modules.R", "modularity", required = FALSE),
  task("R/modularity/tests_1950/TestD_2modules_AP.R", "modularity", required = FALSE),
  task("R/modularity/tests_1950/TestE_2modules_DV.R", "modularity", required = FALSE),
  task("R/modularity/tests_CT3/TestA_5modules_CT3.R", "modularity", required = FALSE),
  task("R/modularity/tests_CT3/TestB_4modules_CT3.R", "modularity", required = FALSE),
  task("R/modularity/tests_CT3/TestC_3modules_CT3.R", "modularity", required = FALSE),
  task("R/modularity/tests_CT3/TestD_2modules_AP_CT3.R", "modularity", required = FALSE),
  task("R/modularity/tests_CT3/TestE_2modules_DV_CT3.R", "modularity", required = FALSE),
  task("R/modularity/Bonferroni_modularity_integration.R", "modularity", required = FALSE),

  task("R/supplemental/size/1950_raw__centroid_PCA_CsizeColor_habHull.R", "supplemental", required = FALSE),
  task("R/supplemental/size/1950_TPS_raw_synthetic.R", "supplemental", required = FALSE),
  task("R/supplemental/size/Figure4_PCA_CT_timeseries_RAW.R", "supplemental", required = FALSE),
  task("R/supplemental/landscape/1950_bigsmall_figure.R", "supplemental", required = FALSE),
  task("R/supplemental/landscape/bigwatersmallwaterprocD.R", "supplemental", required = FALSE),
  task("R/supplemental/landscape/1950_procD_variation_source.R", "supplemental", required = FALSE),
  task("R/supplemental/landscape/1950_variation_source_figure.R", "supplemental", required = FALSE),
  task("R/supplemental/landscape/1950_variation_source_pc3_figure.R", "supplemental", required = FALSE),
  task("R/supplemental/landscape/1950_variation_source_alt_mainstem_vs_smalltrib.R", "supplemental", required = FALSE),
  task("R/supplemental/landscape/mainstem_alltime_vs_tribandres_color.R.R", "supplemental", required = FALSE),
  task("R/supplemental/landscape/mainsteam_alltime_vs_tribandres.R", "supplemental", required = FALSE),
  task("R/supplemental/landscape/mainstem_alltime_vs_tribandres_SICB.R", "supplemental", required = FALSE),
  task("R/supplemental/quabbin_swift/1950_Quabbin_vs_Swift_PC12_hullsHabitat_pointsByCollection.R", "supplemental", required = FALSE),
  task("R/supplemental/quabbin_swift/Swift_Quabbin_1950_PCA_clean.R", "supplemental", required = FALSE),
  task("R/supplemental/quabbin_swift/QS_context_dependence_analysis.R", "supplemental", required = FALSE),
  task("R/supplemental/trait_pc/shape_change_1950.R", "supplemental", required = FALSE),
  task("R/supplemental/trait_pc/shape_change_along_pcs.R", "supplemental", required = FALSE)
)

all_repo_scripts <- list.files("R", recursive = TRUE, full.names = TRUE, pattern = "\\.[Rr]$")
all_repo_scripts <- all_repo_scripts[!basename(all_repo_scripts) %in% c(
  "master_reproduce_figures.R",
  "00_GPA.R", "01_metadata.R", "helpers_angles.R", "01_helpers_angles.R",
  "01_build_angle_measurements.R", "build_trait_measurements.R",
  "mouth_angle_trait_analysis.R", "_run_modularity_integration_1950.R",
  "_run_modularity_integration_CT3.R", "1950_stats_models.R"
)]
unplanned <- setdiff(normalizePath(all_repo_scripts, winslash = "/", mustWork = FALSE),
                     normalizePath(tasks$script, winslash = "/", mustWork = FALSE))

source_check <- check_source_paths(list.files("R", recursive = TRUE, full.names = TRUE, pattern = "\\.[Rr]$"))
write.csv(source_check, file.path(log_dir, "source_path_check.csv"), row.names = FALSE)
missing_sources <- source_check[!source_check$exists, , drop = FALSE]

if (length(unplanned) > 0) {
  writeLines(unplanned, file.path(log_dir, "unplanned_scripts.txt"))
  warning("Some imported scripts are not in the master execution plan. See unplanned_scripts.txt")
}

if (nrow(missing_sources) > 0) {
  write.csv(missing_sources, file.path(log_dir, "missing_source_paths.csv"), row.names = FALSE)
  stop("Missing source paths detected. See ", file.path(log_dir, "missing_source_paths.csv"))
}

missing_task_files <- tasks$script[!file.exists(tasks$script)]
if (length(missing_task_files) > 0) {
  stop("Master task files missing: ", paste(missing_task_files, collapse = ", "))
}

if (!INCLUDE_INTERACTIVE) {
  tasks <- tasks[!tasks$interactive, , drop = FALSE]
}

if (CHECK_ONLY) {
  write.csv(tasks, file.path(log_dir, "master_task_plan.csv"), row.names = FALSE)
  message("Check-only mode complete. Task plan written.")
  quit(status = 0)
}

results <- data.frame(
  phase = character(),
  script = character(),
  status = character(),
  seconds = numeric(),
  error = character(),
  outputs_created = integer(),
  stringsAsFactors = FALSE
)
all_created_outputs <- character()

for (i in seq_len(nrow(tasks))) {
  script <- tasks$script[i]
  phase <- tasks$phase[i]
  message("\n[", i, "/", nrow(tasks), "] ", phase, ": ", script)

  before <- list_outputs()
  start <- Sys.time()
  err <- NULL
  status <- "ok"

  tryCatch(
    source(script, local = .GlobalEnv),
    error = function(e) {
      err <<- conditionMessage(e)
      status <<- if (isTRUE(tasks$required[i])) "failed_required" else "failed_optional"
    }
  )

  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  after <- list_outputs()
  created <- setdiff(after, before)

  if (length(created) > 0) {
    all_created_outputs <- unique(c(all_created_outputs, created))
    out_file <- file.path(log_dir, paste0(sprintf("%03d", i), "_outputs.txt"))
    writeLines(created, out_file)
  }

  results <- rbind(results, data.frame(
    phase = phase,
    script = script,
    status = status,
    seconds = elapsed,
    error = err %||% "",
    outputs_created = length(created),
    stringsAsFactors = FALSE
  ))
  write.csv(results, file.path(log_dir, "master_results.csv"), row.names = FALSE)

  if (status == "failed_required" && STOP_ON_ERROR) {
    stop("Required task failed: ", script, "\n", err)
  }
}

generated <- all_created_outputs
expected <- reference_outputs(reference_root)

comparison <- data.frame(
  generated_file = generated,
  basename = basename(generated),
  matched_reference_count = vapply(basename(generated), function(b) sum(basename(expected) == b), integer(1)),
  stringsAsFactors = FALSE
)
write.csv(comparison, file.path(log_dir, "generated_vs_reference_basenames.csv"), row.names = FALSE)

summary_lines <- c(
  paste("run_id:", MASTER_RUN_ID),
  paste("tasks:", nrow(tasks)),
  paste("ok:", sum(results$status == "ok")),
  paste("failed_required:", sum(results$status == "failed_required")),
  paste("failed_optional:", sum(results$status == "failed_optional")),
  paste("generated_outputs:", length(generated)),
  paste("generated_outputs_matching_reference_basenames:", sum(comparison$matched_reference_count > 0))
)
writeLines(summary_lines, file.path(log_dir, "SUMMARY.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

if (any(results$status == "failed_required")) {
  quit(status = 1)
}
